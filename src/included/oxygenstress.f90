! File VersionID:
!   $Id: oxygenstress.f90 378 2018-05-08 13:50:52Z heine003 $
!
! ----------------------------------------------------------------------

! ## MH: December 2017
! Changes to optimize performance of routine OxygenStress.
! In summary:
!  * certain repeatedly calculated constants have been predefined as constant parameters for further use:
!        Fac3230 = 32.0d0/30.0d0
!        d_o2inwater_ref = 1.0d-5*1.2d0/(100.0d0*100.0d0)
!        SecsPerDay = 24.0d0*3600.0d0
!        Tref3 = 1.0d0/(293.0d0**3)
!        pi4 = 4.0d0*pi
!        h1        =   25000.0d0
!        h2        =  762500.0d0
!        h3        = 1500000.0d0
!        log10h1   = dlog10(h1)
!        log10h2   = dlog10(h2)
!        log10h3   = dlog10(h3)
!        log10h3h2 = log10h3 - log10h2
!        FourThird = 4.0d0/3.0d0
!  * Since division is slower than multiplication, where relevant changes were made; e.g. /10.0 = 0.1*
!  * Other change: ()**0.5 replaced by sqrt()
!  * For the calculation of d_soil, constant properties per node were repeatedly calcuated.
!        New: this has been done only once (routine calc_ini_pars) and these properties were stored
!             for all nodes (1..numnod)
!  * From a profiling analysis it became evident that the numerical integration in the routine 
!    waterfilmthickness was a high-demanding CPU unit.
!        The originally chosen method, trapezoidal rule, is straightforward and robust. 
!        If we assume that our function is sufficiently smooth and has no singularities (also not at the end points)
!        then it is better to use the Romberg Integration method (see Press et al., 1987, Numerical Recipes). 
!        This method takes "many, many fewer" function evaluations than the trapezoidal rule alone, 
!        even though the Romberg method makes use of the trapezoidal rule itself.
!        Implementation of this method decreased the CPU time drastically for a loam case study.
!        Currently both methods are still present in the code, and can be switched by a logical parameter flag
!        that is defined in the subroutine waterfilmthickness: 
!           flUseQromb = .true. (use Romberg integration; else use orginal trapezoidal integration)
!  * For convenience, filling the table ResultsOxStr has been disabled, and the code has been moved 
!    to separate routines (FillOxygenStress1, FillOxygenStress2). First reason: the output was never used
!    in the remainder of SWAP, so it was assumed that this table was only used for debugging purposes.
!    Second remark: this table was defined in file Variables.f90, and again locally here. This could be improved
!    if still needed in the future.
!  * Function FUNC makes use of the derivative of the water retention curve. The equation provided in
!    Bartholomeus et al. (2008; J. Hydrol 360:147-165; their appendix A.4) can be further simplified given by:
!        C = (sat_water_cont - res_water_cont) * alpha * ( (ax)**(gen_n - 1.d0) ) * ( ( 1.d0 + axn )**(-gen_m - 1.d0) )* gen_m*gen_n
!           where ax = alpha*x and axn = (alpha*x)**gen_n
!    Therefore, initially three constants per node are calculated ands stored for later use:
!        Capac_term(i) = (WCs-WCr)*Alpha*N*M
!        Nmin1(i)      = N-1
!        Mplus1(i)     = M+1
!  * Two constants initially set in Calc_ini_pars: shape_factor_microbialr, shape_factor_rootr
!  * An alternative solution for function SOLVE has been implemented: using ZBREND from Numerical Recipes (Press et al., 1987).
!        In function SOLVE the programmer can switch between original approach (Newton-Raphson) or 
!        the new approach (ZBREND) via a logical parameter (useZBREND). 
!        New method resulted in negligable differences, but is somewhat faster.
!        For the implementation of the ZBREND method with the required defintion of a function (myfunc)
!        it was helpful to move a large number of variables into a module (O2_pars). In this way it was
!        prevented to carry-over a large number of arguments from SOLVE to ZBREND to myfunc.
!        In principle, O2_pars can also bue used in other routines/functions resulting in less arguments
!        in their calls. This has not yet been implemented.
    
module O2_pars
   use variables, only: c_mroot,f_senes,max_resp_factor,q10_root,q10_microbial,shape_factor_rootr,specific_resp_humus,ztopcp
   real(8), save :: w_root,w_root_z0
   real(8), save :: soil_temp
   real(8), save :: sat_water_cont,gas_filled_porosity
   real(8), save :: d_o2inwater,d_root,perc_org_mat,soil_density
   real(8), save :: d_soil   
   real(8), save :: depth
   real(8), save :: shape_factor_microbialr,root_radius
   real(8), save :: r_microbial_z0
   real(8), save :: waterfilm_thickness,bunsencoeff
   real(8), save :: c_min_micro,c_macro,ctopnode
end module O2_pars

! ## MH      subroutine OxygenStress(node,rwu_factor,ResultsOxStr) 
      subroutine OxygenStress(node,rwu_factor) 
! ----------------------------------------------------------------------
!     Last modified      : January 2014              
!     Purpose            : calculates oxygen stress according to Bartholomeus et al. (2008)
! ----------------------------------------------------------------------
      use variables
      use O2_pars, only: w_root,w_root_z0, soil_temp, sat_water_cont,gas_filled_porosity, d_o2inwater,d_root,        &
                         perc_org_mat,soil_density,depth, shape_factor_microbialr,root_radius, waterfilm_thickness,  &
                         bunsencoeff, c_min_micro, c_macro,ctopnode,r_microbial_z0, d_soil
      
      implicit none

! --- local
      integer glit,lay,node,i,j
      real(8) resp_factor,rwu_factor,xi,accuracy
      real(8) theta0
      real(8) matric_potential
      real(8) air_temp,alpha,d_gassfreeair,gen_n,o2_atmosphere,percentage_sand,surface_tension_water
      real(8) SOLVE,pi
      real(8) soilphystab(7,matab)
      real(8) diff_water_cap_actual
      integer numrec_tab
      real(8) rdepth,rdens,afgen   !RB20140110
      real(8) rdepth_top,rdens_top !RB20140114 

! ## MH : why ResultsOxStr as argument and locally stored, whereas ResultsOxygenStress was already globally defined in varibles ???
! ## MH      real(8) ResultsOxStr(19,macp)    ! result oxygen stress; tabulated for each model compartment
! ## MH : calcluations ResultsOxStr disabled because it is never used anywhere in SWAP
! ## MH : alternative use switch (flag) to choose whether or not to calculate ResultsOxStr (for debugging purposes)
      
      real(8) top1,top2
      
      parameter (pi = 3.1415926535d0)
      real(8), parameter         :: Fac3230 = 32.0d0/30.0d0      ! ## MH
      real(8), dimension(macp)   :: d_soil_term1, d_soil_term2, gfp100
      real(8), dimension(macp)   :: Capac_term, Nmin1, Mplus1
      logical                    :: ini = .true.

!     save values of locals
      save
      !save   waterfilm_thickness,r_microbial_z0,d_soil
      !save   d_soil_term1, d_soil_term2, ini, gfp100

!## MH : some initial calculations      
      if (ini) then
         if (iHWCKmodel(layer(node)) == 3) then
            call fatalerr ('OxygenStress', 'Combination of OxygenStress and bi-modal MvG (iHWCKmodel=3) is not (yet) possible!')
         end if
         call calc_ini_pars (numnod)
         ini = .false.
      end if
!## MH: end

! --- Get max_resp_factor, i.e. the ratio between total respiration and maintenance respiration
      if (node.eq.1) call GET_MAX_RESP_FACTOR(max_resp_factor)

! --- initialize
      c_min_micro = 0.1d0
      c_macro     = 0.2d0
      resp_factor = 1.0d0

! --- soil layer in SWAP
      lay = layer(node)

! --- dry weight of root per unit length of root [kg/m]
      w_root = 1.0d0/SRL
! --- root radius [m]                  ## MH: ()**0.5 replaced by dsqrt()
      if (swrootradius .eq. 1) then
        root_radius = dsqrt((w_root/(pi*dry_mat_cont_roots*             &
     &              (1-air_filled_root_por)*spec_weight_root_tissue))-  &
     &              (var_a))
      endif
      if (swrootradius .eq. 2) then
        root_radius = root_radiusO2
      endif

! --- RB20140117 get wofost parameters
      if ((croptype(icrop) .eq. 2).or.(croptype(icrop) .eq. 3)) then
          q10_root = q10
          c_mroot = rmr*Fac3230 !CH2O --> O2
      endif
      if (croptype(icrop) .eq. 2) then
          f_senes=afgen(rfsetb,30,dvs)    
      endif
      if (croptype(icrop) .eq. 3) then
          f_senes=afgen(rfsetb,30,rid)    
      endif
   
! --- extract a number of variables from Swap for local use in module OxygenStress

! --- set soil density [kg m-3]      
      soil_density = bdens(lay) 
! --- set parameter n of soil hydraulic functions           
      gen_n = cofgen(6,node) 
! --- set saturated water content [-]      
      sat_water_cont = cofgen(2,node) 
! --- set parameter alpha [1/Pa] of soil hydraulic functions, so divide main swap alpha by 100     ## MH: =0.01* 
      alpha = 0.01d0*cofgen(4,node)
! --- set percentage organic matter [%]
      perc_org_mat = orgmat(lay)*100.0d0 
! --- set percentage sand in % of total soil            
      percentage_sand = (psand(lay)*(1.0d0-orgmat(lay)))*100.0d0 
! --- get soil moisture content as defined in further calculations within this routine [-]      
      theta0 = theta(node) 
! --- gas filled porosity      
      gas_filled_porosity = sat_water_cont-theta0    
      if (h(node) >= 0.0d0) gas_filled_porosity = 0.0d0     ! be sure when saturated that gas_filled_porosity = 0

! --- thickness of the soil compartment [m]     ## MH: =0.01* 
      depth = 0.01d0*dz(node)
! --- temperature in the soil compartment [K]
      soil_temp = tsoil(node)+273.d0  
! --- dry weight of roots at nodal depth
!RB20140109 start new calculation of w_root_z0
!previous:  w_root_z0 = w_root_ss * exp(0.01*z(node)/shape_factor_rootr)  
!new:  
! --- static crop. w_root_z0 relative to value of top layer              
      if (croptype(icrop) .eq. 1) then
            rdepth_top = -ztopcp(1)/rd ! (-z(1)-0.5d0*dz(1))/rd
            rdens_top  = afgen(rdctb,22,rdepth_top)
            rdepth     = -ztopcp(node)/rd ! (-z(node)-0.5d0*dz(node))/rd
            rdens      = afgen(rdctb,22,rdepth) 
            w_root_z0  = w_root_ss * rdens/rdens_top !static crop
      endif
! --- calculate wrootz0 [kg/m3] at top of the compartments !adj RB 20171201
! --- dynamic crop. wrt [kg/ha] = 10-4 kg/m2; 
      if ((croptype(icrop) .eq. 2) .or. (croptype(icrop) .eq. 3)) then
        top1 = dabs(ztopcp(node) / rd) ! relative depth top
        top2 = top1 + 1.0d-6 ! define 'infinite' thin layer; fraction
        
        w_root_z0 = 1.0d6*                                   & ! rescale fraction to 1 (top 2)
     &   (afgen(cumdens,202,top2)-afgen(cumdens,202,top1)) * & ! fraction
     &            (wrt*0.0001d0*(1.0d0/(0.01d0*rd)))           ! wrt kg/ha --> kg/m2; rd cm -> m 
      endif

!JKRO20171114_temp output for testing only
!      write(99,'(f15.7,a1,i5,a1,i5,a1,f10.3,a1,e12.3)')                   &
!     &  t1900,",",node,",",nint(dz(node)),",",z(node),",",w_root_z0


! --- Calculate matric potential [Pa]
       matric_potential = -100.0d0 * h(node)          
! --- if gas filled porosity = 0, then root water uptake = 0. Store results and go to end of routine.        
      if (gas_filled_porosity .lt. 1.0d-4) then !RB20131106 .eq. 0
         ! RB 20140106 start if statement added; if max_resp_factor = 1 then no stress so rwufactor = 1        
          if (max_resp_factor .gt. 1.0d0) then
            rwu_factor = 0.0d0
          else
            rwu_factor = 1.0d0
          endif
          c_macro = 0.d0 !RB20140825 added
! --- store value for output results
          call FillOxygenStress1 ()
      
      else   ! (if (gas_filled_porosity .lt. 1.0d-6)) 
        
! --- In case of tabular soil hydraulic functions
        if(swsophy.eq.1) then
! --- Get tabular soil hydraulic function for node
          do i = 1, 7
            do j = 1, numtablay(lay)  !check this
              soilphystab(i,j) = sptab(i,node,j)
            end do
          end do
! --- Get differential water capacity at actual node, (/L --> /Pa)
          diff_water_cap_actual = 0.01d0*dimoca(node)
          numrec_tab = j-1
        endif

! --- atmosphere oxygen concentration [kg/m3] according to general gas law
        if (node.eq.1) then
! ---   atmospheric temperature [K]      
          air_temp = tav + 273.0d0
          o2_atmosphere = (672.0d0) / (8.314472d0 * air_temp)
! ---   PLAATS O2_atmosphere IN DE VECTOR VOOR C_TOP 
          C_top(1) = o2_atmosphere
        endif    

! --- Calculate temperature dependent parameters
        call TEMP_DEPENDENT_PARAMETERS (d_o2inwater,d_root,             &
            d_gassfreeair,surface_tension_water,bunsencoeff,soil_temp)

! --- Calculate diffusivity Dsoil
!## MH: d_soil_term1, d_soil_term2 and gfp100 were initially calculated and stored per node
        d_soil = d_gassfreeair*d_soil_term1(node) *                     &
     &          ((gas_filled_porosity/gfp100(node))**d_soil_term2(node))
!## MH: end
        
! --- Calculate the thickness of the water film that surrounds the roots      
        call waterfilmthickness (waterfilm_thickness,                   &
     &    matric_potential,Capac_term(node),Nmin1(node),Mplus1(node),   &
     &       alpha,gen_n,surface_tension_water,glit,                    &
     &          soilphystab,diff_water_cap_actual,numrec_tab)

! --- Calculate microbial respiration rate
        call microbial_resp (r_microbial_z0,soil_temp,perc_org_mat,     &
     &      soil_density,percentage_sand,matric_potential,              &
     &         specific_resp_humus,q10_microbial)

! --- Calculate sink term variable
! --- Define input for the solving procedure
        xi       = 0.5d0*max_resp_factor
        accuracy = 1.0d-4
        ctopnode = c_top(node)

! --- Calculate actual respiration factor from the solving procedure
        resp_factor =  SOLVE(xi,accuracy)

! --- PLAATS C_MACRO IN DE VECTOR VOOR C_TOP. BEREKENDE WAARDE IS INPUT VOOR VOLGENDE COMPARTIMENT
       C_top(node+1) = C_macro

! --- Calculate the sink term (Root Water Uptake) variable due to oxygen stress.
! --- The decrease in root water uptake is assumed proportional to the decrease
! --- in respiration (given by the maximum and actual respiration factor).

! RB 20140106 start if statement added; if max_resp_factor = 1 then no stress so rwufactor = 1        
          if (max_resp_factor .gt. 1.0d0) then
            rwu_factor = (1.0d0/(max_resp_factor-1.0d0)) * resp_factor  &
     &                 - (1.0d0/(max_resp_factor-1.0d0))
          else
! Ruud Bartholomeus: 19-12-2017
! Wat betreft zuurstofstress als er geen enkele groei van de wortels meer is:
! -  De lineaire afname van RWU met afname van groeirespiratie van de wortels moeten we handhaven. Rwu_factor = 1 als groeirespiratie is optimaal, rwu_factor = 0 als geen groeirespiratie.
! -  Als er geen groeirespiratie van de wortels meer is, maar wel groei van bovengrondse delen is er nog wel transpiratie en dus wateropname van de wortels
! -  Als er geen groeirespiratie van de wortels meer nodig is, maar er wel voldoende zuurstof is voor onderhoudsrespiratie, dan functioneren de wortels dus zoals het moet en nemen ook water op. Als er onvoldoende zuurstof is voor onderhoudsrespiratie van de wortels, dan sterven deze. Ze nemen dan geen water op. 
! -  Voorstel voor verbetering in de code: 
!        o  Als max_resp_factor = 1 en resp_factor < 1 DAN rwu_factor = 0
!        o  Hier zit dus geen geleidelijke afname: onvoldoende zuurstof voor onderhoud = stop 

             if (resp_factor < 1.0d0) then
               rwu_factor = 0.0d0      !## MH: suggested by RB 19-12-2017: rwu_factor = 0  : orignal code before 19-12-2017: rwu_factor = 1.0d0
            else
               rwu_factor = 1.0d0
            end if
          endif
! RB 20140106 end if statement added

          if (rwu_factor .gt. 1.0d0) then
             rwu_factor = 1.0d0
          endif
          if (rwu_factor .lt. 0.d0) then
              rwu_factor = 0.d0
          endif
      
      endif !if (gas_filled_porosity .lt. 1.0d-6) !RB20131216 goto removed
      
! --- store value for output results
      call FillOxygenStress2 ()
      return
      
   contains
   
   subroutine calc_ini_pars (numnod)
!## MH : these constant can be calculated only once at initial call to OxygenStress
!## MH : results are stored and saved per node
!## MH : since it is part of subroutine OxygenStress (via the previous statement "contains")
!## MH : the output variables are known at this level and not need to be part of the list of arguments
   integer :: numnod, i
   real(8), parameter :: h100      = -100.0d0
   real(8), parameter :: h500      = -500.0d0
   real(8), parameter :: log10h100 = dlog10(-h100)
   real(8), parameter :: log10h500 = dlog10(-h500)
   real(8) theta100,theta500,campbell_b,watcon

   do i = 1, numnod
! --- get theta at h=-100 cm and at h=-500 cm; for diffusion coef
      theta100        = watcon(i,h100)
      theta500        = watcon(i,h500)
      sat_water_cont  = cofgen(2,i)
      gfp100(i)       = sat_water_cont - theta100
      campbell_b      = (log10h500-log10h100) / (dlog10(theta100)-dlog10(theta500))
      d_soil_term1(i) = 2.0d0*(gfp100(i)**3)+0.04d0*gfp100(i)
      d_soil_term2(i) = 2.0d0+3.0d0/campbell_b
!     for use in FUNC
!     cofgen(x,i): x = 1...12
!     1 = thetar; 2 = thetas, 3 = Ksatfit; 4 = alpha; 5 = lambda; 6 = n; 7 = m (=1-1/n); 
!     8 = dummy; 9 = h_enpr; 10 = Ksatexm; 11 = relsatthr; 12 = Ksatthr
!     Calculate (sat_water_cont - res_water_cont) * alpha * gen_m * gen_n
      Capac_term(i)   = (cofgen(2,i) - cofgen(1,i)) * 0.01d0*cofgen(4,i) * cofgen(6,i) * cofgen(7,i)
      Nmin1(i)        = cofgen(6,i) - 1.0d0
      Mplus1(i)       = cofgen(7,i) + 1.0d0 
   end do
! --- microbial respiration calculated from organic matter content in actual soil compartment; keep this value fixed
   shape_factor_microbialr = 0.9d0 
! --- microbial respiration calculated from organic matter content in actual soil compartment; keep this value fixed
   shape_factor_rootr = 0.9d0 
   return
   end subroutine calc_ini_pars
   
   subroutine FillOxygenStress1 ()
!!!!!        if (node.eq.1) ResultsOxStr = -999.d0 !RB20140114
!!!!!          ResultsOxStr(1,node) = node
!!!!!          ResultsOxStr(2,node) = z(node)*0.01d0
!!!!!          ResultsOxStr(3,node) = h(node)
!!!!!          ResultsOxStr(4,node) = matric_potential
!!!!!          ResultsOxStr(5,node) = gas_filled_porosity
!!!!!          ResultsOxStr(6,node) = -999.d0
!!!!!          ResultsOxStr(7,node) = -999.d0
!!!!!          ResultsOxStr(8,node) = rwu_factor
!!!!!! --- Rpotz0 kg/m3/d
!!!!!          ResultsOxStr(9,node)= f_senes * (c_mroot * w_root_z0 *        &
!!!!!     &      max_resp_factor) *                                          &
!!!!!     &      ( q10_root ** ( (soil_temp - 298.0d0) / 10.0d0 ) )
!!!!!! --- RpotLay kg/m2/d; integrate over function potresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
!!!!!          ResultsOxStr(10,node)= -1.0d0*((ResultsOxStr(9,node)*         &
!!!!!     &      shape_factor_rootR*dexp(-depth/shape_factor_rootr))-      &
!!!!!     &      (ResultsOxStr(9,node)*shape_factor_rootr*                   &
!!!!!     &      dexp((0.0d0)/shape_factor_rootr)))
!!!!!! --- Ractz0 kg/m3/d
!!!!!          if (resp_factor.gt.max_resp_factor) then
!!!!!            resp_factor = max_resp_factor
!!!!!          endif
!!!!!          if (resp_factor.lt.1.0d0) then
!!!!!            resp_factor = 0.0d0
!!!!!          endif      
!!!!!          ResultsOxStr(11,node)= f_senes * (c_mroot * w_root_z0 *       &
!!!!!     &      resp_factor) *                                              &
!!!!!     &      ( q10_root ** ( (soil_temp - 298.0d0 ) / 10.d0 ) ) 
!!!!!! --- RactLay kg/m2/d; integrate over function actresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
!!!!!        ResultsOxStr(12,node)= -1.d0*((ResultsOxStr(11,node)*           &
!!!!!     &      shape_factor_rootR*dexp(-depth/shape_factor_rootr))-        &
!!!!!     &      (ResultsOxStr(11,node)*shape_factor_rootr*                  &
!!!!!     &    dexp((0.d0)/shape_factor_rootr)))
!!!!!! --- RredLay kg/m2/d
!!!!!          ResultsOxStr(13,node)= ResultsOxStr(10,node)-                 &
!!!!!     &      ResultsOxStr(12,node)    
!!!!!! --- RsoilActLay kg/m2/d     
!!!!!          ResultsOxStr(14,node)= -999.d0
!!!!!! --- dsoil
!!!!!          ResultsOxStr(15,node) = -999.d0     
!!!!!! --- c_macro
!!!!!          ResultsOxStr(16,node)=-999.d0 !kg/m3
!!!!!          ResultsOxStr(17,node)=-999.d0       
!!!!!! --- Max resp factor RB20140115
!!!!!          ResultsOxStr(18,node)=max_resp_factor 
!!!!!! --- wrt from wofost
!!!!!          ResultsOxStr(19,node)=-999.d0 !kg/ha          
   return
   end subroutine FillOxygenStress1
   subroutine FillOxygenStress2 ()
!!!!!      if (node.eq.1) ResultsOxStr = -999.0d0 !RB20140114
!!!!!      ResultsOxStr(1,node) = node
!!!!!      ResultsOxStr(2,node) = z(node)*0.01d0
!!!!!      ResultsOxStr(3,node) = h(node)
!!!!!      ResultsOxStr(4,node) = matric_potential
!!!!!      ResultsOxStr(5,node) = gas_filled_porosity
!!!!!      ResultsOxStr(6,node) = w_root_z0
!!!!!      ResultsOxStr(7,node) = waterfilm_thickness
!!!!!      ResultsOxStr(8,node) = rwu_factor
!!!!!! --- Rpotz0 kg/m3/d
!!!!!      ResultsOxStr(9,node)= f_senes * (c_mroot * w_root_z0 *            &
!!!!!     &    max_resp_factor) *                                            &
!!!!!     &    ( q10_root ** ( (soil_temp - 298.0d0 ) / 10.0d0 ) )
!!!!!! --- RpotLay kg/m2/d; integrate over function potresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
!!!!!      ResultsOxStr(10,node)= -1.0d0*((ResultsOxStr(9,node)*             &
!!!!!     &  shape_factor_rootR*dexp(-depth/shape_factor_rootr))-            &
!!!!!     &  (ResultsOxStr(9,node)*shape_factor_rootr*                       &
!!!!!     &  dexp((0.0d0)/shape_factor_rootr)))
!!!!!! --- Ractz0 kg/m3/d
!!!!!      if (resp_factor.gt.max_resp_factor) then
!!!!!        resp_factor = max_resp_factor
!!!!!      endif
!!!!!      if (resp_factor.lt.1.0d0) then
!!!!!        resp_factor = 0.0d0
!!!!!      endif      
!!!!!      ResultsOxStr(11,node)= f_senes * (c_mroot * w_root_z0 *           &
!!!!!     &  resp_factor) *                                                  &
!!!!!     &  ( q10_root ** ( (soil_temp - 298.0d0 ) / 10.0d0 ) ) 
!!!!!! --- RactLay kg/m2/d; integrate over function actresp_top_of_compartment*exp(x/shape_factor_rootR); analytical solution
!!!!!      ResultsOxStr(12,node)= -1.d0*((ResultsOxStr(11,node)*             &
!!!!!     &  shape_factor_rootR*dexp(-depth/shape_factor_rootr))-            &
!!!!!     &  (ResultsOxStr(11,node)*shape_factor_rootr*                      &
!!!!!     &  dexp((0.d0)/shape_factor_rootr)))
!!!!!! --- RredLay kg/m2/d
!!!!!      ResultsOxStr(13,node)= ResultsOxStr(10,node)-                     &
!!!!!     & ResultsOxStr(12,node)    
!!!!!! --- RsoilActLay kg/m2/d     
!!!!!      ResultsOxStr(14,node)= -1.d0*((r_microbial_z0*                    &
!!!!!     &   shape_factor_microbialr*dexp(-depth/shape_factor_microbialr))- &
!!!!!     &     (r_microbial_z0*shape_factor_microbialr*                     &
!!!!!     &      dexp((0.d0)/shape_factor_microbialr)))
!!!!!! --- dsoil
!!!!!      ResultsOxStr(15,node) = d_soil     
!!!!!! --- c_macro
!!!!!      ResultsOxStr(16,node)=c_macro !kg/m3
!!!!!      ResultsOxStr(17,node)=100.d0*                                     &
!!!!!     & c_macro/((0.032*1d5)/(8.314472d0*soil_temp)) !% or kPa       
!!!!!! --- Max resp factor RB20140115
!!!!!      ResultsOxStr(18,node)=max_resp_factor    
!!!!!! --- wrt from wofost RBf20140120
!!!!!      if ((croptype(icrop) .eq. 2).or.(croptype(icrop) .eq. 3)) then
!!!!!        ResultsOxStr(19,node)=wrt    
!!!!!      endif
   return
   end subroutine FillOxygenStress2
      
      end subroutine OxygenStress

! --- End of main module OxygenStress ---------------------------------------------------------------------------------------
      
      subroutine GET_MAX_RESP_FACTOR (max_resp_factor_gmrf)
      use Variables
      implicit none
! --- Procedure to derive max_resp_factor, 
! --- i.e. the ratio between total respiration and maintenance respiration [-]
! --- This ratio is either given in the input file (for a static crop) 
! --- or calculated from a series of equations taken from WOFOST (for a dynamic crop)

      real(8) afgen
      real(8) rmres_gmrf,teff_gmrf,mres_gmrf,asrc_gmrf
      real(8) fr_gmrf,fl_gmrf,fs_gmrf,fo_gmrf   
      real(8) cvf_gmrf
      real(8) Froots, Rg_roots,Rm_roots,Max_resp_factor_gmrf        
! --- static crop        
! --- static crop: max_resp_factor is given in the input file
      if (croptype(icrop) .eq. 1) then
        max_resp_factor_gmrf = max_resp_factor
      endif !if (croptype(icrop) .eq. 1)

! --- dynamic crop: max_resp_factor is calculated following the procedure
! --- for the calculation of root maintenance respiration and root growth respiration as 
! --- used in WOFOST.
! --- dynamic crop, not grass 
      if (croptype(icrop) .eq. 2) then

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration, based on actual plant state variables
        rmres_gmrf = (rmr*wrt+rml*wlv+rms*wst+rmo*wso)*                 &
     &            afgen(rfsetb,30,dvs)
        !teff_gmrf = q10**((tsoil(10)-25.0d0)/10.0d0) !TEMPORARY!!!! ONLY TO CHECK EFFECT OF USING TSOIL INSTEAD OF TAV; ## MH: /10 = 0.1*
        teff_gmrf = q10**(0.1d0*(tav-25.0d0))      

        mres_gmrf = min(pgass,rmres_gmrf*teff_gmrf)  ! ## MM 2018-05-07
        asrc_gmrf = pgass - mres_gmrf                ! ## MM 2018-05-07
! --- partitioning factors
        fr_gmrf = afgen(frtb,30,dvs) !rid for grass, dvs for wofost
        fl_gmrf = afgen(fltb,30,dvs)
        fs_gmrf = afgen(fstb,30,dvs)
        fo_gmrf = afgen(fotb,30,dvs)
! --- dry matter increase, only part in which cvf is calculated
        cvf_gmrf = 1.0d0/((fl_gmrf /cvl+fs_gmrf /cvs+fo_gmrf/cvo)*      &
     &  (1.0d0-fr_gmrf)+fr_gmrf/cvr)
        
! --- cvf: factor used in wofost to calculate the increase in biomass (dmi) from the
! ---  net assimilation of the whole plant (asrc); dmi = cvf*asrc. 
! ---  What is left is the growth respiration (i.e. asrc = dmi + growth respiration). 
! ---  Therefore, growth respiration of the whole plant = asrc*(1-cvf)
! --- Froots: contribution of the roots to cvf
        Froots = (fr_gmrf/cvr)*cvf_gmrf
! --- Rg_roots: growth respiration roots        
        Rg_roots = Froots*(1.0d0-cvf_gmrf)*asrc_gmrf
! --- Rm_roots: maintenance respiration roots        
        Rm_roots = min(Froots*(1.0d0-cvf_gmrf)*pgass,                   &
     &      rmr*wrt*afgen(rfsetb,30,dvs)*teff_gmrf)
! --- Max_resp_factor: ratio total respiration / maintenance respiration        
        if (Rm_roots.gt.0.0d0) then
            Max_resp_factor_gmrf = (Rg_roots+Rm_roots)/Rm_roots
        else 
            Max_resp_factor_gmrf = 1.0d0  
        endif         
      endif !if (croptype(icrop) .eq. 2)
        
! --- dynamic crop, grass 
      if (croptype(icrop) .eq. 3) then        
        Max_resp_factor_gmrf = 1.0d0  !RB20140317
! --- skip in case of regrowth, equal to wofost detailed grass
! --- note: daycrop.ge.idregrpot (wofost) --> daycrop.gt.idregrpot, because idregrpot is result of wofost of previous day
        if (daycrop.eq.0 .or.daycrop.gt.idregr) then           

! --- respiration and partitioning of carbohydrates between growth and
! --- maintenance respiration, based on actual plant state variables
          rmres_gmrf = (rmr*wrt+rml*wlv+rms*wst)*afgen(rfsetb,30,rid) 
!        teff_gmrf = q10**((tsoil(10)-25.0d0)/10.0d0) !TEMPORARY!!!! ONLY TO CHECK EFFECT OF USING TSOIL INSTEAD OF TAV; ## MH: /10=*0.1
          teff_gmrf = q10**(0.1d0*(tav-25.0d0))

          mres_gmrf = min(pgass,rmres_gmrf*teff_gmrf)  ! ## MM 2018-05-07
          asrc_gmrf = pgass - mres_gmrf                ! ## MM 2018-05-07
! --- partitioning factors
          fr_gmrf = afgen(frtb,30,rid) !rid for grass, dvs for wofost
          fl_gmrf = afgen(fltb,30,rid)
          fs_gmrf = afgen(fstb,30,rid)
! --- dry matter increase, only part in which cvf is calculated
          cvf_gmrf = 1.0d0/((fl_gmrf /cvl+fs_gmrf /cvs)*                &
     &        (1.0d0-fr_gmrf)+fr_gmrf/cvr)
! --- cvf: factor used in wofost to calculate the increase in biomass (dmi) from the
! ---  net assimilation of the whole plant (asrc); dmi = cvf*asrc. 
! ---  What is left is the growth respiration (i.e. asrc = dmi + growth respiration). 
! ---  Therefore, growth respiration of the whole plant = asrc*(1-cvf)
! --- Froots: contribution of the roots to cvf
          Froots = (fr_gmrf/cvr)*cvf_gmrf
! --- Rg_roots: growth respiration roots            
          Rg_roots = Froots*(1.0d0-cvf_gmrf)*asrc_gmrf
! --- Rm_roots: maintenance respiration roots     
          Rm_roots = min(Froots*(1.0d0-cvf_gmrf)*pgass,                &
     &        rmr*wrt*afgen(rfsetb,30,rid)*teff_gmrf)
! --- Max_resp_factor: ratio total respiration / maintenance respiration        
          if (Rm_roots.gt.0.d0) then
              Max_resp_factor_gmrf = (Rg_roots+Rm_roots)/Rm_roots
          else 
              Max_resp_factor_gmrf = 1.0d0  
          endif                  
        endif !RB20140317 #skip in case of regrowth                   
      endif !if (croptype(icrop) .eq. 3)

      return
      end
      
      subroutine TEMP_DEPENDENT_PARAMETERS (d_o2inwater,d_root,         &
     &   d_gassfreeair,surface_tension_water,bunsencoeff,soil_temp)

      implicit none
      real(8) d_o2inwater,d_root,d_gassfreeair
      real(8) surface_tension_water,bunsencoeff,soil_temp
!     ## MH: d_o2inwater_ref = 24.0d0*3600.0d0*((1.0d-5*1.2d0/(100.0d0*100.0d0))
      real(8), parameter :: d_o2inwater_ref = 1.0d-5*1.2d0/(100.0d0*100.0d0)
      real(8), parameter :: SecsPerDay = 24.0d0*3600.0d0
!     ## MH: Tref3 = 1/(293.0d0**3)
      real(8), parameter :: Tref3 = 1.0d0/(293.0d0**3)
      
! --- Calculate parameters that depend on temperature
! --- diffusion coefficient for oxygen in water [m2/d] Lango et al 1996
!     ## MH:
      d_o2inwater = SecsPerDay*(d_o2inwater_ref * dexp(0.026d0*(soil_temp-273.0d0)))
! --- diffusion coefficient for oxygen in root tissue [m2/d]
! --- scaled to d_o2inwater, according the value given by
! --- van noordwijk & de willigen 1987 at 293 k
      d_root = 0.4d0*d_o2inwater
! --- diffusion coefficient for oxygen in free air [m2/d]
! --- hirschfelder et al 1964 molecular theory of gases and liquids
!     ## MH: Tref3 = 1/(293.0d0**3)
      d_gassfreeair = 1.74528d0 * (soil_temp**3)*Tref3
! --- surface tension of water [n/m] eotvos rule
      surface_tension_water = 0.07275d0 *                               &
     &                      ( 1.0d0 - 0.002d0 * (soil_temp - 291.0d0 ) )
! --- bunsen solubility coeff for oxygen Lango et al 1996
      bunsencoeff =                                                     &
!## MH     &         (1413.d0*(exp(-144.397d0+7775.18d0*(soil_temp**(-1.d0))  &
     &         (1413.d0*(dexp(-144.397d0+7775.18d0/soil_temp            &
     &              +18.3977d0*dlog(soil_temp)+0.0094437d0*soil_temp))) &
     &              *273.15d0/soil_temp 
      return
      end

      real(8) function FUNC(x,Capac_term,Nmin1,Mplus1,alpha,gen_n,      &
     &                      surface_tension_water)

      implicit none
      real(8) x
      real(8) Capac_term,Nmin1,Mplus1
      real(8) alpha,gen_n,surface_tension_water
      real(8) pi
      parameter (pi = 3.1415926535897932d0)
      real(8), parameter :: pi4 = 4.0d0*pi         ! ## MH
      real(8)  :: ax, axn                          ! ## MH
! --- function needed for calculation of length density of gas filled pores in
! ---   subroutine 'water film thickness'                                          
      ax  = alpha * x
      axn = ax**gen_n
      func= Capac_term * ax**Nmin1 * ( ( 1.d0 + axn )**(-Mplus1) ) /    &
     &      (pi4 * surface_tension_water**2/x**2)
     !!!!! func= -(( -(sat_water_cont - res_water_cont) * alpha *            &    
     !!!!!&      ( (ax)**(gen_n - 1.d0) ) *                                  &
     !!!!!&      ( ( 1.d0 + axn )**(gen_m - 1.d0) )*                         &
     !!!!!&      gen_m*gen_n ) /                                             &
     !!!!!&      (((( axn ) + 1.d0 )**gen_m)**2) ) /                         &
     !!!!!&      (pi4*(((surface_tension_water**2)/(x**2))))
     !!! func= -( -(sat_water_cont - res_water_cont) * alpha *            &    
     !!!&      ( (ax)**(gen_n - 1.d0) ) *                                  &
     !!!&      ( ( 1.d0 + axn )**(-gen_m - 1.d0) )*                        &
     !!!&      gen_m*gen_n ) /                                             &
     !!!&      (pi4*(((surface_tension_water**2)/(x**2))))
      return
      end

      real(8) function FUNCtab(x,diff_water_cap,surface_tension_water)

      implicit none
      real(8) x
      real(8) diff_water_cap,surface_tension_water
      real(8) pi
      parameter (pi = 3.1415926535897932d0)
      real(8), parameter :: pi4 = 4.0d0*pi         ! ## MH
! --- function needed for calculation of length density of gas filled pores in
! ---   subroutine 'water film thickness'                                          
      functab= diff_water_cap / (pi4*(surface_tension_water**2)/(x**2))
      return
   end     

!## MH: Romberg could be a better (=faster) method; ask Ruud if he considered this
      subroutine TRAPZD(a,b,s,n,Capac_term,Nmin1,Mplus1,alpha,gen_n,    &
     &                  surface_tension_water,glit)
     
      implicit none
      real(8) FUNC
      real(8) a,b,s
      integer n
! --- procedure from numerical recipes. calculate integral numerically
! --- Needed for procedure 'water film thickness'                                
! --- programs calling trapzd must provide a function
! --- func(x:real(8)):real(8)which is to be integrated. they must
! --- also define the variable
! --- var glit: integer;
! --- in the main routine. *)
!RB20140312 sum replaced by mysum
      integer j,glit
      real(8) x,tnm,mysum,del
      real(8) Capac_term,Nmin1,Mplus1
      real(8) alpha,gen_n,surface_tension_water
     
! ---  initialize j
      j=0
      if (n .eq. 1) then
         s = 0.5d0*(b-a)*(func(a,Capac_term,Nmin1,Mplus1,alpha,gen_n,   &
     &                         surface_tension_water)                   &
     &                  + func(b,Capac_term,Nmin1,Mplus1,alpha,gen_n,   &
     &                         surface_tension_water))
         glit = 1
      else
         tnm = dble(glit)
         del = (b-a)/tnm
         x = a+0.5d0*del
         mysum = 0.0d0
         do while (j .lt. glit) 
            mysum = mysum+func(x,Capac_term,Nmin1,Mplus1,alpha,gen_n,   &
     &                         surface_tension_water)
            x = x+del
            j = j + 1
         enddo
         s = 0.5d0*(s+(b-a)*mysum/tnm)
         glit = 2*glit
      endif
      return
      end

      subroutine TRAPZDtab(a,b,s,n,diff_water_cap,                      &
     &                     surface_tension_water,glit)
     
      implicit none
      real(8) functab
      real(8) a,b,s
      integer n
! --- procedure from numerical recipes. calculate integral numerically
! --- Needed for procedure 'water film thickness'                                
! --- programs calling trapzd must provide a function
! --- func(x:real(8)):real(8)which is to be integrated. they must
! --- also define the variable
! --- var glit: integer;
! --- in the main routine. *)
!RB20140312 sum replaced by mysum
      integer j,glit
      real(8) x,tnm,mysum,del
      real(8) diff_water_cap,surface_tension_water
! ---  initialize j
      j=0
      if (n .eq. 1) then
         s = 0.5d0*(b-a)*(functab(a,diff_water_cap,                     &
     &                            surface_tension_water)                &
     &                   +functab(b,diff_water_cap,                     &
     &                            surface_tension_water))
         glit = 1
      else
         tnm = dble(glit)
         del = (b-a)/tnm
         x = a+0.5d0*del
         mysum = 0.0d0
         do while (j .lt. glit) 
            mysum = mysum+functab(x,diff_water_cap,                     &
     &                            surface_tension_water)
            x = x+del
            j = j + 1
         enddo
         s = 0.5d0*(s+(b-a)*mysum/tnm)
         glit = 2*glit
      endif
      return
      end

      subroutine waterfilmthickness (waterfilm_thickness,               &
     &   matric_potential,Capac_term, Nmin1, Mplus1,                    &
     &      alpha,gen_n,surface_tension_water,glit,                     &
     &         soilphystab,diff_water_cap_actual,numrec_tab)
! --- calculate water film thickness. method according to simojoki 2000
      use Variables
      use doln
      implicit none
      
      real(8) waterfilm_thickness, matric_potential
      real(8) alpha,gen_n,surface_tension_water,Capac_term, Nmin1, Mplus1
      real(8) soilphystab(7,matab)
      real(8) diff_water_cap_actual
      real(8) lowlim,upplim,diff_water_cap,htab
      real(8) length_density_gas_pores_sub
      integer glit    
      integer i,ii,numrec_tab
      real(8) new,old,s
      real(8) length_density_gas_pores
      logical, parameter :: flUseQromb = .true. !## MH: (use Romberg integration; else use orginal trapezoidal integration)
      
      save s

      real(8) pi
      parameter (pi = 3.1415926535897932d0)
      
!     Extremely dry; to prevent problems in integration below (observed for coarse sand) we set waterfilm_thickness to som low value
!        This will ensure no oxygen related stress (?)
      if (matric_potential > 1.0d7) then
         waterfilm_thickness = 1.0d-8
         return
      end if
      
      if (swsophy.eq.0) then
! --- calculate length density air filled (gas) pores. [number per m2]
! --- subroutine trapdz is used to solve the integral defined in
! --- function 'func'.
        i = 1
        if (flUseQromb) then
           call QROMBD(1.d-10,matric_potential,s,Capac_term,            &
     &         Nmin1,Mplus1,alpha,gen_n,surface_tension_water,glit)
        else
           call TRAPZD(1.d-10,matric_potential,s,i,Capac_term,          &
     &         Nmin1,Mplus1,alpha,gen_n,surface_tension_water,glit)
        end if
        new = s
        old = new + new

        do while (dabs(new-old) .gt. 0.00001d0*new) 
          i = i + 1
          if (flUseQromb) then
            call QROMBD(1.d-10,matric_potential,s,Capac_term,           &
     &         Nmin1,Mplus1,alpha,gen_n,surface_tension_water,glit)
          else
            call TRAPZD(1.d-10,matric_potential,s,i,Capac_term,         &
     &         Nmin1,Mplus1,alpha,gen_n,surface_tension_water,glit)
          end if
          old = new
          new = s
        enddo
! --- result from trapzd and func:      
        length_density_gas_pores = new
! --- calculated water film thickness [m]
        waterfilm_thickness = 2.d0 *                                    &
     &   ( ( dsqrt( 1.0d0 / ( pi * length_density_gas_pores ) ) )       &  ! ## MH: replaced ()**0.5 by dsqrt()
     &   - 2.0d0 * surface_tension_water / matric_potential )
      endif
      !!!!sptab(1,node,ii) = soilphystab(1,ii)
      if(swsophy.eq.1) then
          ii = numrec_tab !34  !run over points in soil hydraulic table
! --- lower value of matric potential for integration interval
         lowlim = 0.000001d0 !-100*soilphystab(1,ii) !initial value should be zero !RB20140725, very close to zero, otherwise/0 in functab
! --- upper value of matric potential for integration interval
         if (do_ln_trans) then
            upplim = -100.d0*0.5d0*(-(dexp(-soilphystab(1,ii))-1.0d0) - (dexp(-soilphystab(1,ii-1))-1.0d0)) !middle of points 1 and 2
         else
            upplim = -100.d0*0.5d0*(soilphystab(1,ii)+soilphystab(1,ii-1)) !middle of points 1 and 2
         end if
         upplim = MIN(upplim, matric_potential) !upplim can be higher than matpot          
! --- initialize length density gas filled pores
         length_density_gas_pores = 0.0d0 !initial value
! --- start loop, i.e. run over full integration interval;
! --- in each loop, the length density of gas filled pores of that interval is calculated;
! --- The sum of all loops gives the total length density.
! --- The integration interval runs from 0 to the matric potential at the actual node
         do while ((upplim.lt.matric_potential) .and. (ii.gt.2)) !adjusted RB20150714 ii.gt.2 added
            if (do_ln_trans) then
               htab = -(dexp(-soilphystab(1,ii))-1.0d0)
               diff_water_cap = 0.01d0*soilphystab(4,ii)/(-htab+1.0d0) !at matric potential centre of lowlim and upplim, except for first point
            else
               diff_water_cap = 0.01d0*soilphystab(4,ii) !at matric potential centre of lowlim and upplim, except for first point
            end if
            if (flUseQromb) then
               call QROMBDtab(lowlim,upplim,s,diff_water_cap,           & !n=1 as integral is over a linear function, i.e. convergence is reached after one step
     &                        surface_tension_water,glit)
            else
               call TRAPZDtab(lowlim,upplim,s,1,diff_water_cap,         & !n=1 as integral is over a linear function, i.e. convergence is reached after one step
     &                        surface_tension_water,glit)
            end if
            length_density_gas_pores_sub = s
            length_density_gas_pores = length_density_gas_pores +       &
     &                                 length_density_gas_pores_sub
            ii = ii - 1
            lowlim = upplim
!##MH            upplim=-100.d0*0.5d0*(soilphystab(1,ii)+soilphystab(1,ii-1))
            if (do_ln_trans) then
               upplim = -100.d0*0.5d0*(-(dexp(-soilphystab(1,ii))-1.0d0) - (dexp(-soilphystab(1,ii-1))-1.0d0)) !middle of points 1 and 2
            else
               upplim = -100.d0*0.5d0*(soilphystab(1,ii)+soilphystab(1,ii-1)) !middle of points 1 and 2
            end if
          enddo
! --- in the last step (upplim >= matric_potential), always use the differential water capacity corresponding to the actual matric_potential of the node
          diff_water_cap = diff_water_cap_actual 
          upplim = MIN(upplim, matric_potential)
          if (flUseQromb) then
            call QROMBDtab(lowlim,upplim,s,diff_water_cap,              & !n=1 as integral is over a linear function, i.e. convergence is reached after one step
     &                     surface_tension_water,glit)
          else
            call TRAPZDtab(lowlim,upplim,s,1,diff_water_cap,            &  ! n=1 as integral is over a linear function, i.e. convergence is reached after one step
     &                     surface_tension_water,glit)
          end if
          length_density_gas_pores_sub = s
          length_density_gas_pores = length_density_gas_pores +         &
     &                               length_density_gas_pores_sub
! --- calculated water film thickness [m]
          waterfilm_thickness = 2.0d0 *                                 &
     &     ( ( dsqrt( 1.d0 / ( pi * length_density_gas_pores ) ) )      &  ! ## MH: replaced ()**0.5 by dsqrt()
     &       - 2.d0 * surface_tension_water / matric_potential )
      endif
      !ResultsOxStr(7,node) = waterfilm_thickness !RB20140115
      return
      end

      subroutine microbial_resp (r_microbial_z0,soil_temp,perc_org_mat, &
     &   soil_density,percentage_sand,matric_potential,                 &
     &      specific_resp_humus,q10_microbial)
! --- Calculate microbial respiration rate, based on available amount of
! --- organic matter, soil moisture and temperature

      implicit none
      real(8) r_microbial_z0,soil_temp,perc_org_mat
      real(8) soil_density,percentage_sand,matric_potential
      real(8) specific_resp_humus,q10_microbial
      real(8) f_moisture_humus,t_humus,carbon_humuspools
      real(8) saturated_matric_potential

!     ## MH : define constants
      real(8), parameter :: h1        =   25000.0d0
      real(8), parameter :: h2        =  762500.0d0
      real(8), parameter :: h3        = 1500000.0d0
      real(8), parameter :: log10h1   = dlog10(h1)
      real(8), parameter :: log10h2   = dlog10(h2)
      real(8), parameter :: log10h3   = dlog10(h3)
      real(8), parameter :: log10h3h2 = log10h3 - log10h2
!     ## MH: end

! --- humus temperature [K]
      t_humus = soil_temp
! --- available amount of organic carbon [kg/m3]            ## MH: /100 = *0.01
      carbon_humuspools = 0.48d0 * ( 0.01d0*perc_org_mat ) *soil_density
! --- saturated matric potential [pa] calculated according to cosby et al. 1984
! --- * 100: cm -> pa
      saturated_matric_potential =(10.d0**(-0.0131d0*percentage_sand+   &
     &                            1.88d0))*100.0d0
! --- reduction function for soil moisture:
      if ( matric_potential .lt. saturated_matric_potential) then
         f_moisture_humus = 0.5d0
      endif
      if (( matric_potential .ge. saturated_matric_potential)           &
     &                      .and. (matric_potential .lt. h1)) then
         f_moisture_humus = 1.d0 - 0.5d0 *                              &
     &                      ((log10h1-dlog10(matric_potential))/        &
     &                      (log10h1-                                   &
     &                      dlog10(saturated_matric_potential) ) )
      endif
      if ((matric_potential .ge. h1)                                    &
     &                    .and. (matric_potential .le. h2)) then
         f_moisture_humus = 1.0d0
      endif
      if ((matric_potential .gt. h2)                                    &
     &                    .and. (matric_potential .le. h3)) then
         f_moisture_humus = 1.d0 -                                      &
     &                     ((dlog10(matric_potential)-log10h2)/         &
     &                     ( log10h3h2 ) )
      endif
      if (matric_potential .gt. h3) then
         f_moisture_humus = 0.0d0
      endif
! --- microbial respiration rate at soil surface [kg/m3/d]
      r_microbial_z0 = specific_resp_humus * carbon_humuspools *        &
     &             ( q10_microbial**( 0.1d0*(t_humus-298.d0) ) )        &
     &             * f_moisture_humus
      return
      end

      subroutine MICRO (c_mroot,w_root,f_senes,q10_root,                &
     &                 soil_temp,sat_water_cont,gas_filled_porosity,    &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, resp_factor)
! --- calculate minimum oxygen concentration in the gas phase of the soil
! --- that is needed to provide all cells within a plant root with a
! --- sufficient amount of oxygen.                                                                     *)

      implicit none
      real(8) c_min_micro,c_mroot,w_root,f_senes,q10_root,soil_temp
      real(8) sat_water_cont,gas_filled_porosity
      real(8) d_o2inwater,d_root,perc_org_mat,soil_density
      real(8) specific_resp_humus,q10_microbial,depth
      real(8) shape_factor_microbialr,root_radius
      real(8) waterfilm_thickness,bunsencoeff
      real(8) resp_factor

      real(8) r_mref,r_mroot,waterfilm_porosity,d_waterfilm,lambda
      real(8) r_waterfilm_lengthroot,alpha_alpha
      real(8) f_moisture_humus,t_humus,carbon_humuspools
      real(8) c_min_micro_interphase

      real(8) r_microbial_volumetric_wf ,r_microbial_z0_wf
      real(8) pi
      parameter (pi = 3.1415926535897932d0)
      real(8), parameter :: FourThird = 4.0d0/3.0d0      ! ## MH

! --- calc respiration per unit length of root [kg/m/d]
      r_mref = c_mroot * w_root * resp_factor

      r_mroot = f_senes * r_mref *                                      &
     &          ( q10_root ** ( 0.1d0*(soil_temp - 298.d0 ) ) )
! --- calc porosity of water film [-]
      waterfilm_porosity = ( sat_water_cont - gas_filled_porosity ) /   &
     &                ( 1.0d0 - gas_filled_porosity )
      waterfilm_porosity = max(0.075d0,waterfilm_porosity) !RB20180507, set minimum value for very dry conditions
! --- calc diffusion coeff of water film [m2/d]
      d_waterfilm = d_o2inwater * ( waterfilm_porosity**FourThird )
! --- calc ratio of d_root and d_waterfilm, input for c_min_micro_interphase
      lambda = d_root / d_waterfilm
! --- start microbial respiration rate in water film
! --- temperature humus [k]
      t_humus = soil_temp
! --- amount of carbon in humus [kg/m3]
      carbon_humuspools = 0.48d0 * (0.01d0*perc_org_mat) * soil_density
! --- reduction factor for moisture [-] (saturation)
      f_moisture_humus = 0.5d0
! --- microbial respiration rate at soil surface [kg/m3/d]
      r_microbial_z0_wf = specific_resp_humus * carbon_humuspools *     &
     &                 ( q10_microbial**( 0.1d0*(t_humus-298.d0) ) ) *  &
     &                 f_moisture_humus
! --- volumetric microbial resp rate at depth z [kg/m3/d]
      r_microbial_volumetric_wf = r_microbial_z0_wf*                    &
     &                          dexp(-depth/shape_factor_microbialr)
! --- respiration rate in water film per unit length of root [kg/m/d]
      r_waterfilm_lengthroot = pi *                                     &
     &                     ( (root_radius + waterfilm_thickness)**2 -   &
     &                     root_radius**2 )*r_microbial_volumetric_wf
! --- end microbial respiration rate in water film*)

! --- calculation procedure following de willigen & van noordwijk 1983
! --- p.220-221
! --- ratio of rhizosphere respiration to total respiration [-]
      alpha_alpha = r_waterfilm_lengthroot /                            &
     &              (r_waterfilm_lengthroot+r_mroot)
! --- volumetric respiration rate of the root + rhizosphere, but attributed
! ---   to the root (following de willigen & van noordwijk 1983) [kg/m3/d]
      c_min_micro_interphase = ((r_waterfilm_lengthroot + r_mroot)/     &
     &          (2.d0*pi*d_root)) *                                     &
     &          ( 0.5d0 + ( ( lambda - 1.d0 ) * alpha_alpha / 2.d0 ) +  &
     &          lambda * dlog(1.d0 + waterfilm_thickness / root_radius)-&
     &          ( lambda * alpha_alpha *                                &
     &          ( 1.d0 + waterfilm_thickness / root_radius )**2 ) *     &
     &          dlog(1.d0 + waterfilm_thickness / root_radius ) /       &
     &          ( ( waterfilm_thickness / root_radius ) *               &
     &          (2.d0 + waterfilm_thickness / root_radius) )  )

      c_min_micro = c_min_micro_interphase/bunsencoeff
! --- end calculation procedure following de willigen & van noordwijk 1983*)
      return
      end

      subroutine MACRO (c_macro,depth,resp_factor,                      &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
     
! --- calculate oxygen concentration in the soil at a specific depth with respect
! --- to the soil surface. diffusion from atmosphere into soil is considered,
! --- including microbial and root respiration.

      implicit none
      real(8) c_macro, depth, resp_factor
      real(8) c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctop
      real(8) shape_factor_microbialr,shape_factor_rootr,r_microbial_z0

      real(8) dum,r_mref_z0,r_mroot_z0,d_soil,l,lnew,lnew_initial,fi,   &
     &        fi_a
      integer counterMacro, counterMacroSub
      character(len=200) error_messag

! --- reference respiration at soil surface per unit volume of roots [kg/m3] *)
      r_mref_z0 = c_mroot * w_root_z0 * resp_factor

! --- respiration at soil surface per unit volume of roots [kg/m3]
      r_mroot_z0 = f_senes * r_mref_z0 *                                &
     &             ( q10_root ** ( 0.1d0*(soil_temp - 298.d0 ) ) )

! --- calculate oxygen concentration at certain depth (c_macro [kg/m3]
! --- in the soil. method according to cook 1995
! --- calculate criterium (dum) to switch to specific equation (if then else)
      dum = (( shape_factor_microbialr**2) * r_microbial_z0 /d_soil)+   &
     &      ( ( shape_factor_rootr**2 ) * r_mroot_z0 / d_soil )
! --- as z goes to infinity, the oxygen concentration asymptotically approaches
! ---    a constant non-zero value
      if (dum .lt. ctop) then
         c_macro = ctop -                                               &
     &     ( (shape_factor_microbialr**2) * r_microbial_z0 / d_soil)*   &
     &     ( 1.d0- dexp( - depth / shape_factor_microbialr ) ) -        &
     &     ( (shape_factor_rootr**2) * r_mroot_z0 / d_soil ) *          &
     &     ( 1.d0- dexp( - depth / shape_factor_rootr ) )
! --- at z = l the oxygen concentration goes to zero. find l through
! ---    newton-raphson method
      else
         counterMacro = 0
         counterMacroSub = 0
         lnew_initial = 0.1d0 !initialize
         fi = 1.d0 !initialize
             fi_a = 1.d0 !initialize
         lnew = lnew_initial
         do while ((dabs(fi).gt.1.d-8) .AND.                            &
     &                 (dabs(fi_a).gt.0.d0) .AND.                       &
     &             (lnew.gt.1.d-4))
            counterMacro = counterMacro + 1
            counterMacroSub = counterMacroSub + 1
            l = lnew
            fi = ctop -                                                 &
     &         ((shape_factor_microbialr**2)*r_microbial_z0/d_soil)*    &
     &         ( 1.d0- (l/shape_factor_microbialr)*                     &
     &         dexp( - l / shape_factor_microbialr ) -                  &
     &         dexp( - l / shape_factor_microbialr ) ) -                &
     &         ( (shape_factor_rootr**2) * r_mroot_z0 / d_soil ) *      &
     &         (1.d0-(l/shape_factor_rootr)*                            &
     &         dexp(-l/shape_factor_rootr )-                            &
     &         dexp( - l / shape_factor_rootr ) )
! --- derivative of fi to l
            fi_a = -r_microbial_z0/d_soil*l*                            &
     &             dexp(-l/shape_factor_microbialr)-                    &
     &             r_mroot_z0/d_soil*l*dexp(-l/shape_factor_rootr)
               
            if (dabs(fi_a) .gt. 0.d0) then
                lnew = dabs(l - (fi / fi_a))
            endif
            
! --- for convergence --> depends on initial value of lnew
            if (lnew.gt.1.d3                                            &
     &           .or.                                                   &
     &           counterMacroSub.gt.100) then
               !restart do while loop with new lnew_initial value
               lnew_initial = lnew_initial + 0.1d0
               lnew = lnew_initial
               counterMacroSub = 0
            endif
! ---       fatal error if too many iterations --> 
            if (counterMacro .gt. 1.d6) then
              error_messag = '1 Too much iterations for macroscopic '   &
     &               //' oxygen diffusion.'
!D              call warn ('rootextraction',error_messag,logf,swscre)
              call fatalerr ('rootextraction',error_messag)
            endif            
         enddo
        
         if (depth .lt. l) then
            c_macro = ctop -                                            &
     &         ((shape_factor_microbialr**2)*r_microbial_z0/d_soil )*   &
     &         ( 1.d0- (depth/shape_factor_microbialr)*                 &
     &         dexp( - l / shape_factor_microbialr ) -                  &
     &         dexp( - depth / shape_factor_microbialr ) ) -            &
     &         ( (shape_factor_rootr**2) * r_mroot_z0 / d_soil ) *      &
     &         ( 1.d0- (depth/shape_factor_rootr)*                      &
     &         dexp( - l / shape_factor_rootr ) -                       &
     &         dexp( - depth / shape_factor_rootr ) )
         else
            c_macro = 0.0d0
         endif
      endif
      return
      end

      real(8) function SOLVE (xi,accuracy)
     
      use O2_pars
      implicit none
      
      real(8) xiplus1,ximin1
      real(8) delta
      real(8) xi,accuracy
      real(8) xx,fxi,xiplusdelta,ximindelta,fxiplusdelta,fximindelta,   &
     &        dfxi
      integer counterSolve
      character(len=200) error_messag

!      logical, parameter :: UseZBREND = .false.
      logical, parameter :: UseZBREND = .true.
      real(8)  :: a, b, Dif_a, Dif_b, ZBREND, myfunc
      external :: myfunc
      
      if (UseZBREND) then
! ## MH: start 
         a     = 0.0d0 ! 1.0d-6
         Dif_a = myfunc(a)
         b     = max_resp_factor
         Dif_b = myfunc(b)
      
         if (Dif_a*Dif_b < 0.0d0) then
            SOLVE = ZBREND (myFUNC,a,b,Dif_a,Dif_b,accuracy)
            return
         else
            ! same sign
            if (Dif_a > 0.0d0) then
               if (Dif_b < Dif_a) then
                  xx = max_resp_factor
               else
                  xx = 0.0d0
               end if
            else
               if (Dif_b < Dif_a) then
                  xx = 0.0d0
               else
                  xx = max_resp_factor
               end if
            end if

            call MACRO (c_macro,depth,xx,                               &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode, &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
            SOLVE = xx
            return
         end if
         return
! ## MH: end

      else !## MH: (not use ZBREND; i.e. use original Newton-Raphson approach)
      
! --- iterative procedure to find the RESPIRATION FACTOR for which holds !RB20131106
! --- that c_macro = c_micro. reference value = c_micro
! --- Newton-Raphson method

! --- speed up simulations: cut of if ctop = 0 and if waterfilm_thickness is extremely high (due to very low gas filled porosity)     !RB20140826
      if (dabs(ctopnode) .lt. 1.0d-6 .OR.                               &
     &    waterfilm_thickness .gt. 1.d0) then 
         xx = 0.d0
         C_macro = 0.0d0
         SOLVE = xx
         return
      endif

      delta= 1.d-8
      fxi = 100.d0
      counterSolve = 1
      do while (fxi .gt. accuracy)
         call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,         &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, xi)
         call MACRO (c_macro,depth,xi,                                  &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
         fxi = dabs(C_min_micro-C_macro)
         xiplusdelta = xi + delta
         ximindelta = xi - delta
         
         if (dabs(xiplusdelta-ximindelta).lt.1.d-20) then !if (dabs(xiplusdelta-ximindelta).lt.1.d-12) then !if (xiplusdelta .eq. ximindelta) then 
            xiplusdelta = xi + 1.d-6
            ximindelta = xi - 1.d-6
         endif
         
         call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,         &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, xiplusdelta)
         call MACRO (c_macro,depth,xiplusdelta,                         &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
     
         fxiplusdelta = dabs(C_min_micro-C_macro)
         call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,         &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro, ximindelta)
         call MACRO (c_macro,depth,ximindelta,                          &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)    
     
         fximindelta = dabs(C_min_micro-C_macro)  

         dfxi = (fxiplusdelta - fximindelta) / (xiplusdelta-ximindelta)
         
        if (dabs(dfxi) .lt. 1.d-20) then !if (dabs(dfxi) .lt. 1.d-15) then !RB20140312
          xiplus1 = xi
        else
          xiplus1 = max(1.0d-8,(xi - fxi/dfxi)) !xiplus1 = max(1.0d-6,(xi - fxi/dfxi)) !RB20131106, never lt 0
        endif

      !prevent endless sign change without conversion !RB20140827
      if (mod(counterSolve,2)>1.d-12) then !if even than store xi
        ximin1 = xi
      endif
!D      write(*,*) counterSolve,xi,ximin1,xiplus1 !DEBUG
      if (dabs(xiplus1-ximin1) .lt. 1.d-6) then  !no change in xi value, than take average of both values between which is iterated
        xx = (xi+xiplus1)*0.5d0
        call MACRO (c_macro,depth,xx,                                   &
     &       c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode,         &
     &       shape_factor_microbialr,shape_factor_rootr,                &
     &       r_microbial_z0,d_soil)     
        SOLVE = xx
        return
      endif

! --- speed up simulations: cut of if RWU_factor will be >1 or <1e-6 (=0)       
         if (xi .lt. 1.0d-6 .AND. xi .gt. 0.0d0) then !RB20131106 greater than 0 is required
            xx = 0.d0
            C_macro = 0.0d0
            SOLVE = xx
            return
         endif
         if (xi .gt. max_resp_factor) then
            xx = max_resp_factor
            call MACRO (c_macro,depth,xx,                               &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode,     &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)     
            SOLVE = xx
            return
         endif
      
         xi = xiplus1
         
! ---    error if max iterations reached--> 
         counterSolve = counterSolve + 1
         if (counterSolve .gt. 100) then
           error_messag = '1 Max iterations for solver oxygen'          &
     &            //' stress reached.'
            xx = xi
           SOLVE = xx
           call fatalerr ('rootextraction',error_messag)
           return          
        endif            
      enddo
      xx = xi
      SOLVE = xx
      return

      end if ! (UseZBREND)

      end function SOLVE

! ----------------------------------------------------------------------
      subroutine oxygen_dat (SwTopSub,NrStaring,OxygenSlope,            &
     &                       OxygenIntercept)
! ----------------------------------------------------------------------
!     date               : december 2009
!     purpose            : set parameter values of metafunction for oxygenstress
! ----------------------------------------------------------------------

! --- global
      integer SwTopSub,NrStaring
      real(8) OxygenSlope(6),OxygenIntercept(6)

! --- local
      integer i
      real(8) xtop(18,6),xsub(18,6),ytop(18,6),ysub(18,6)

      data (xtop(1,i),i=1,6)                                            &
     & /5.07d-03,2.40d+02,-4.39d+00,-6.31d+02,1.67d+00,9.08d+02/
      data (xtop(2,i),i=1,6)                                            &
     & /1.20d-02,3.26d+02,-8.91d+00,-9.29d+02,2.49d+00,1.64d+03/
      data (xtop(3,i),i=1,6)                                            &
     & /1.21d-02,3.64d+02,-9.09d+00,-9.90d+02,2.60d+00,1.70d+03/
      data (xtop(4,i),i=1,6)                                            &
     & /1.67d-02,4.42d+02,-1.21d+01,-1.26d+03,3.34d+00,2.20d+03/
      data (xtop(5,i),i=1,6)                                            &
     & /4.17d-03,1.93d+02,-3.72d+00,-5.28d+02,1.44d+00,7.77d+02/
      data (xtop(6,i),i=1,6)                                            &
     & /2.11d-02,5.02d+02,-1.51d+01,-1.40d+03,3.69d+00,2.71d+03/
      data (xtop(7,i),i=1,6)                                            &
     & /2.86d-02,5.57d+02,-1.97d+01,-1.67d+03,4.50d+00,3.41d+03/
      data (xtop(8,i),i=1,6)                                            &
     & /1.75d-02,5.55d+02,-1.32d+01,-1.34d+03,3.31d+00,2.48d+03/
      data (xtop(9,i),i=1,6)                                            &
     & /2.34d-02,6.00d+02,-1.70d+01,-1.40d+03,3.42d+00,3.09d+03/
      data (xtop(10,i),i=1,6)                                           &
     & /3.12d-02,6.62d+02,-2.20d+01,-1.53d+03,3.65d+00,3.91d+03/
      data (xtop(11,i),i=1,6)                                           &
     & /2.58d-02,6.42d+02,-1.83d+01,-1.61d+03,4.00d+00,3.26d+03/
      data (xtop(12,i),i=1,6)                                           &
     & /2.50d-02,6.53d+02,-1.79d+01,-1.45d+03,3.40d+00,3.21d+03/
      data (xtop(13,i),i=1,6)                                           &
     & /2.53d-02,5.97d+02,-1.79d+01,-1.72d+03,4.58d+00,3.18d+03/
      data (xtop(14,i),i=1,6)                                           &
     & /2.82d-02,7.10d+02,-2.04d+01,-1.54d+03,3.56d+00,3.71d+03/
      data (xtop(15,i),i=1,6)                                           &
     & /2.08d-02,4.72d+02,-1.46d+01,-1.37d+03,3.65d+00,2.57d+03/
      data (xtop(16,i),i=1,6)                                           &
     & /1.99d-02,4.80d+02,-1.39d+01,-1.34d+03,3.47d+00,2.45d+03/
      data (xtop(17,i),i=1,6)                                           &
     & /2.27d-02,6.31d+02,-1.62d+01,-1.58d+03,3.91d+00,2.91d+03/
      data (xtop(18,i),i=1,6)                                           &
     & /2.23d-02,6.50d+02,-1.60d+01,-1.74d+03,4.45d+00,2.89d+03/

      data (xsub(1,i),i=1,6)                                            &
     & /7.21d-04,1.76d+02,-1.59d+00,-4.33d+02,1.16d+00,4.52d+02/
      data (xsub(2,i),i=1,6)                                            &
     & /3.75d-03,2.25d+02,-3.61d+00,-5.96d+02,1.59d+00,7.91d+02/
      data (xsub(3,i),i=1,6)                                            &
     & /7.06d-03,2.86d+02,-5.91d+00,-7.52d+02,2.00d+00,1.19d+03/
      data (xsub(4,i),i=1,6)                                            &
     & /1.29d-02,3.74d+02,-9.74d+00,-1.03d+03,2.72d+00,1.82d+03/
      data (xsub(5,i),i=1,6)                                            &
     & /2.92d-03,1.86d+02,-3.06d+00,-5.07d+02,1.39d+00,6.85d+02/
      data (xsub(6,i),i=1,6)                                            &
     & /3.00d-02,5.41d+02,-2.06d+01,-1.69d+03,4.61d+00,3.55d+03/
      data (xsub(7,i),i=1,6)                                            &
     & /2.76d-02,6.63d+02,-1.95d+01,-1.67d+03,4.15d+00,3.48d+03/
      data (xsub(8,i),i=1,6)                                            &
     & /1.85d-02,4.88d+02,-1.34d+01,-1.34d+03,3.50d+00,2.43d+03/
      data (xsub(9,i),i=1,6)                                            &
     & /2.08d-02,5.46d+02,-1.50d+01,-1.50d+03,3.92d+00,2.72d+03/
      data (xsub(10,i),i=1,6)                                           &
     & /1.85d-02,5.94d+02,-1.39d+01,-1.45d+03,3.60d+00,2.60d+03/
      data (xsub(11,i),i=1,6)                                           &
     & /2.29d-02,6.33d+02,-1.67d+01,-1.65d+03,4.21d+00,3.05d+03/
      data (xsub(12,i),i=1,6)                                           &
     & /2.60d-02,5.89d+02,-1.83d+01,-1.33d+03,3.12d+00,3.26d+03/
      data (xsub(13,i),i=1,6)                                           &
     & /3.36d-02,6.27d+02,-2.29d+01,-1.40d+03,3.26d+00,3.95d+03/
      data (xsub(14,i),i=1,6)                                           &
     & /3.84d-02,6.74d+02,-2.71d+01,-1.40d+03,3.21d+00,4.83d+03/
      data (xsub(15,i),i=1,6)                                           &
     & /2.25d-02,6.16d+02,-1.66d+01,-1.37d+03,3.22d+00,3.05d+03/
      data (xsub(16,i),i=1,6)                                           &
     & /1.88d-02,4.75d+02,-1.32d+01,-1.29d+03,3.31d+00,2.33d+03/
      data (xsub(17,i),i=1,6)                                           &
     & /2.19d-02,5.40d+02,-1.53d+01,-1.50d+03,3.87d+00,2.70d+03/
      data (xsub(18,i),i=1,6)                                           &
     & /1.98d-02,5.00d+02,-1.41d+01,-1.40d+03,3.67d+00,2.52d+03/

      data (ytop(1,i),i=1,6)                                            &
     & /1.89d-04,-3.91d-01,-8.65d-02,2.11d+01,-8.61d-02,7.12d+00/
      data (ytop(2,i),i=1,6)                                            &
     & /5.02d-05,-7.56d-02,-1.01d-02,1.80d+01,-7.27d-02,-2.81d+00/
      data (ytop(3,i),i=1,6)                                            &
     & /5.36d-06,3.54d-01,1.20d-02,1.73d+01,-7.09d-02,-5.51d+00/
      data (ytop(4,i),i=1,6)                                            &
     & /-2.67d-05,4.87d-02,3.03d-02,1.76d+01,-7.02d-02,-7.90d+00/
      data (ytop(5,i),i=1,6)                                            &
     & /2.60d-04,-1.74d+00,-1.16d-01,2.28d+01,-8.99d-02,9.67d+00/
      data (ytop(6,i),i=1,6)                                            &
     & /-1.23d-04,-6.11d-02,8.42d-02,1.66d+01,-6.63d-02,-1.51d+01/
      data (ytop(7,i),i=1,6)                                            &
     & /-1.02d-04,8.54d-03,7.12d-02,1.76d+01,-7.06d-02,-1.33d+01/
      data (ytop(8,i),i=1,6)                                            &
     & /-8.09d-05,2.23d-01,5.65d-02,1.62d+01,-6.61d-02,-1.08d+01/
      data (ytop(9,i),i=1,6)                                            &
     & /-9.68d-05,1.56d-01,6.56d-02,1.67d+01,-6.76d-02,-1.21d+01/
      data (ytop(10,i),i=1,6)                                           &
     & /-1.37d-04,-4.78d-02,8.94d-02,1.81d+01,-7.09d-02,-1.55d+01/
      data (ytop(11,i),i=1,6)                                           &
     & /-1.71d-04,-2.25d-01,1.08d-01,1.68d+01,-6.46d-02,-1.78d+01/
      data (ytop(12,i),i=1,6)                                           &
     & /-1.36d-04,-5.72d-01,8.92d-02,1.02d+01,-3.91d-02,-1.54d+01/
      data (ytop(13,i),i=1,6)                                           &
     & /-1.19d-04,-8.68d-03,8.28d-02,1.77d+01,-6.97d-02,-1.54d+01/
      data (ytop(14,i),i=1,6)                                           &
     & /-9.91d-05,-9.32d-01,7.13d-02,1.35d+01,-5.10d-02,-1.35d+01/
      data (ytop(15,i),i=1,6)                                           &
     & /-7.90d-05,2.75d-02,5.85d-02,1.55d+01,-6.20d-02,-1.16d+01/
      data (ytop(16,i),i=1,6)                                           &
     & /-7.91d-05,2.32d-01,5.71d-02,1.10d+01,-4.40d-02,-1.12d+01/
      data (ytop(17,i),i=1,6)                                           &
     & /-2.21d-04,-5.23d-01,1.39d-01,1.44d+01,-5.42d-02,-2.27d+01/
      data (ytop(18,i),i=1,6)                                           &
     & /-1.16d-04,-3.60d-01,7.85d-02,1.05d+01,-3.96d-02,-1.41d+01/

      data (ysub(1,i),i=1,6)                                            &
     & /4.17d-04,-1.41d+00,-2.06d-01,2.42d+01,-9.84d-02,2.20d+01/
      data (ysub(2,i),i=1,6)                                            &
     & /2.56d-04,-4.40d-01,-1.22d-01,2.19d+01,-8.92d-02,1.16d+01/
      data (ysub(3,i),i=1,6)                                            &
     & /1.93d-04,-2.84d-01,-8.88d-02,1.96d+01,-8.08d-02,7.68d+00/
      data (ysub(4,i),i=1,6)                                            &
     & /5.69d-05,-7.44d-02,-1.48d-02,1.80d+01,-7.36d-02,-2.01d+00/
      data (ysub(5,i),i=1,6)                                            &
     & /5.21d-04,-1.81d+00,-2.61d-01,1.96d+01,-7.85d-02,3.01d+01/
      data (ysub(6,i),i=1,6)                                            &
     & /-1.10d-04,6.88d-02,7.90d-02,1.77d+01,-7.11d-02,-1.48d+01/
      data (ysub(7,i),i=1,6)                                            &
     & /-1.76d-04,-3.81d-01,1.11d-01,1.79d+01,-6.90d-02,-1.85d+01/
      data (ysub(8,i),i=1,6)                                            &
     & /-2.38d-05,5.72d-01,2.47d-02,1.64d+01,-6.82d-02,-6.50d+00/
      data (ysub(9,i),i=1,6)                                            &
     & /-4.68d-05,4.19d-01,3.84d-02,1.74d+01,-7.16d-02,-8.59d+00/
      data (ysub(10,i),i=1,6)                                           &
     & /-1.10d-04,1.03d-01,7.32d-02,1.69d+01,-6.77d-02,-1.32d+01/
      data (ysub(11,i),i=1,6)                                           &
     & /-1.48d-04,-3.36d-01,9.58d-02,1.78d+01,-6.94d-02,-1.64d+01/
      data (ysub(12,i),i=1,6)                                           &
     & /-1.58d-04,1.36d-01,9.92d-02,1.56d+01,-6.18d-02,-1.65d+01/
      data (ysub(13,i),i=1,6)                                           &
     & /-2.69d-04,-6.12d-01,1.66d-01,1.27d+01,-4.80d-02,-2.65d+01/
      data (ysub(14,i),i=1,6)                                           &
     & /-6.34d-05,-4.27d-01,5.05d-02,1.73d+01,-6.82d-02,-1.06d+01/
      data (ysub(15,i),i=1,6)                                           &
     & /3.08d-05,5.25d-02,-6.04d-03,8.44d+00,-3.54d-02,-2.00d+00/
      data (ysub(16,i),i=1,6)                                           &
     & /-6.20d-05,3.40d-01,4.76d-02,9.54d+00,-3.85d-02,-9.98d+00/
      data (ysub(17,i),i=1,6)                                           &
     & /-2.90d-05,3.99d-01,2.75d-02,8.07d+00,-3.30d-02,-6.73d+00/
      data (ysub(18,i),i=1,6)                                           &
     & /-6.76d-05,4.41d-01,5.01d-02,1.50d+01,-6.10d-02,-1.01d+01/

      if (SwTopSub .eq. 1) then
        do i = 1,6
          OxygenSlope(i)     = xtop(NrStaring,i)
          OxygenIntercept(i) = ytop(NrStaring,i)
        enddo
      else
        do i = 1,6
          OxygenSlope(i)     = xsub(NrStaring,i)
          OxygenIntercept(i) = ysub(NrStaring,i)
        enddo
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine OxygenReproFunction (OxygenSlope,OxygenIntercept,      &
     &   theta,thetas,tsoil,node,z,dz,rwu_factor)
! ----------------------------------------------------------------------
!     date               : January 2010
!     purpose            : Calculate oxygen stress according to reproduction function
! ----------------------------------------------------------------------
      use variables, only: zbotcp
      implicit none
      include 'arrays.fi'

! --- global
      integer node,i
      real(8) OxygenSlope(6),OxygenIntercept(6),theta(macp),thetas(macp)
      real(8) tsoil(macp),z(macp),dz(macp)

! --- local
      real(8) intercept,slope,sum_porosity
      real(8) gas_filled_porosity
      real(8) soil_temp,depth_ss,mean_gas_filled_porosity
      real(8) rwu_factor

      gas_filled_porosity = thetas(node) - theta(node)
      soil_temp = tsoil(node) + 273.d0
      depth_ss = -z(node) * 0.01d0

      if (gas_filled_porosity .lt. 1.d-10) then
         rwu_factor = 0.d0 
         return
      endif

! --- mean gas filled porosity
      sum_porosity = 0.0d0
      do i = 1,node
        sum_porosity = sum_porosity +                                   &
     &                 (thetas(i) - theta(i)) * dz(i)
      enddo
      mean_gas_filled_porosity = sum_porosity /(-zbotcp(node))

      intercept = OxygenIntercept(1)*soil_temp**2 +                     &
     &            OxygenIntercept(2)*depth_ss**2 +                      &
     &            OxygenIntercept(3)*soil_temp +                        &
     &            OxygenIntercept(4)*depth_ss +                         &
     &            OxygenIntercept(5)*soil_temp*depth_ss +               &
     &            OxygenIntercept(6)

      slope = OxygenSlope(1)*soil_temp**2 +                             &
     &        OxygenSlope(2)*depth_ss**2 +                              &
     &        OxygenSlope(3)*soil_temp +                                &
     &        OxygenSlope(4)*depth_ss +                                 &
     &        OxygenSlope(5)*soil_temp*depth_ss +                       &
     &        OxygenSlope(6)

! --- Calculate the sink term (Root Water Uptake) variable due to oxygen stress.
      rwu_factor = intercept + slope*mean_gas_filled_porosity
      if (rwu_factor .gt. 1.d0) then
         rwu_factor = 1.d0
      endif
      if (rwu_factor .lt. 0.d0) then
         rwu_factor = 0.d0
      endif
      return

   end

   
!-----------------------------------------------------------------------*
! From W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vettering,    *
! 1986. Numerical recipes. The art of scientific computing. Cambridge   *
! University Press, Cambridge, 818 pp.                                  *
! Double Precision version                                              *
!-----------------------------------------------------------------------*
!      SUBROUTINE QROMBD (FUNC,A,B,SS)
      SUBROUTINE QROMBD (a,b,ss,Capac_term,Nmin1,Mplus1,                &
     &                   alpha,gen_n,surface_tension_water,glit)

      IMPLICIT NONE
      real(8) a,b,ss
      integer glit
      real(8) Capac_term,Nmin1,Mplus1
      real(8) alpha,gen_n,surface_tension_water
      INTEGER          JMAX,JMAXP,K,KM,J,L
      REAL(8)          EPS,S,H,DSS
!      REAL(8)          FUNC

      PARAMETER (EPS=1.D-5,JMAX=20,JMAXP=JMAX+1,K=5,KM=K-1)
!      PARAMETER (EPS=1.D-10,JMAX=100,JMAXP=JMAX+1,K=5,KM=K-1)
      DIMENSION S(JMAXP),H(JMAXP)
!      EXTERNAL  FUNC

      H(1) = 1.D0; s(1) = 0.0d0
      DO 11 J = 1, JMAX
!        CALL TRAPZDD (FUNC,A,B,S(J),J)
        CALL TRAPZD (a,b,s(j),j,Capac_term,Nmin1,Mplus1,                &
     &               alpha,gen_n,surface_tension_water,glit)
        IF (J .GE. K) THEN
          L = J - KM
!          CALL POLINTD (H(L),S(L),K,0.D0,SS,DSS)
          CALL POLINTD (H(L:J),S(L:J),K,0.D0,SS,DSS)
          IF (DABS(DSS) .LT. EPS*DABS(SS)) RETURN
        END IF
        S(J+1) = S(J)
        H(J+1) = 0.25D0*H(J)
11    CONTINUE
!      PAUSE 'Too many steps.'
      write (*,*) 'QROMBD Too many steps.'
      !read (*,*)
   END

!-----------------------------------------------------------------------*
! From W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vettering,    *
! 1986. Numerical recipes. The art of scientific computing. Cambridge   *
! University Press, Cambridge, 818 pp.                                  *
! Double Precision version                                              *
!-----------------------------------------------------------------------*
      SUBROUTINE QROMBDtab (a,b,ss,dif_water_cap,surface_tension_water,glit)

      IMPLICIT NONE
      real(8) a,b,ss
      integer glit
      real(8) dif_water_cap, surface_tension_water
      INTEGER          JMAX,JMAXP,K,KM,J,L
      REAL(8)          EPS,S,H,DSS
!      REAL(8)          FUNC

      PARAMETER (EPS=1.D-5,JMAX=100,JMAXP=JMAX+1,K=5,KM=K-1)
!      PARAMETER (EPS=1.D-10,JMAX=100,JMAXP=JMAX+1,K=5,KM=K-1)
      DIMENSION S(JMAXP),H(JMAXP)
!      EXTERNAL  FUNC

      H(1) = 1.D0; s(1) = 0.0d0
      DO 11 J = 1, JMAX
!        CALL TRAPZDD (FUNC,A,B,S(J),J)
        CALL TRAPZDtab (a,b,s(j),j,dif_water_cap,surface_tension_water,glit)
        IF (J .GE. K) THEN
          L = J - KM
!          CALL POLINTD (H(L),S(L),K,0.D0,SS,DSS)
          CALL POLINTD (H(L:J),S(L:J),K,0.D0,SS,DSS)
          IF (DABS(DSS) .LT. EPS*DABS(SS)) RETURN
        END IF
        S(J+1) = S(J)
        H(J+1) = 0.25D0*H(J)
11    CONTINUE
!      PAUSE 'Too many steps.'
      write (*,*) 'QROMBDtab Too many steps.'
      read (*,*)
   END

!-----------------------------------------------------------------------*
! From W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vettering,    *
! 1986. Numerical recipes. The art of scientific computing. Cambridge   *
! University Press, Cambridge, 818 pp.                                  *
! Double Precision version                                              *
!-----------------------------------------------------------------------*
      SUBROUTINE POLINTD (XA,YA,N,X,Y,DY)

      IMPLICIT NONE
      INTEGER          N,NS,I,M,NMAX
      REAL(8)         XA,YA,C,D,X,Y,DY,DIF,DIFT,HO,HP,DEN,W

      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

      NS=1
      DIF=DABS(X-XA(1))
      DO 11 I=1,N
        DIFT=DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        END IF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.D0) then
!             PAUSE 'NR_POLINTD: DEN = 0'
             write (*,*) 'NR_POLINTD: DEN = 0'
             read (*,*)
          end if
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        END IF
        Y=Y+DY
13    CONTINUE
      RETURN
   END

   real(8) function myfunc(x)
   use O2_pars
   real(8) :: x
      call MICRO (c_mroot,w_root,f_senes,q10_root,soil_temp,            &
     &                 sat_water_cont,gas_filled_porosity,              &
     &                 d_o2inwater,d_root,perc_org_mat,soil_density,    &
     &                 specific_resp_humus,q10_microbial,depth,         &
     &                 shape_factor_microbialr,root_radius,             &
     &                 waterfilm_thickness,bunsencoeff,                 &
     &                 c_min_micro,x)
      call MACRO (c_macro,depth,x,                                      &
     &           c_mroot,w_root_z0,f_senes,q10_root,soil_temp,ctopnode, &
     &           shape_factor_microbialr,shape_factor_rootr,            &
     &           r_microbial_z0,d_soil)
      myfunc = c_macro - c_min_micro
   end function myfunc
   
!     Adapted: initial FA and FB are input
      DOUBLE PRECISION FUNCTION ZBREND (FUNC,X1,X2,FA,FB,TOL)
      IMPLICIT NONE
      INTEGER          ITMAX,ITER
      DOUBLE PRECISION EPS,X1,X2,TOL,A,B,FA,FB,FC,C,D,E,XM,TOL1,P,Q,R,S
      DOUBLE PRECISION FUNC
!      PARAMETER (ITMAX=100,EPS=3.D-8)
      PARAMETER (ITMAX=100,EPS=3.D-16)
      A=X1
      B=X2
!      FA=FUNC(A)
!      FB=FUNC(B)
!      IF (FB*FA.GT.0.D0) PAUSE 'Root must be bracketed for ZBREND.'
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.D0) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(DABS(FC).LT.DABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.D0*EPS*DABS(B)+0.5D0*TOL
        XM=.5D0*(C-B)
        IF(DABS(XM).LE.TOL1 .OR. FB.EQ.0.D0)THEN
          ZBREND=B
          RETURN
        ENDIF
        IF(DABS(E).GE.TOL1 .AND. DABS(FA).GT.DABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.D0*XM*S
            Q=1.D0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.D0*XM*Q*(Q-R)-(B-A)*(R-1.D0))
            Q=(Q-1.D0)*(R-1.D0)*(S-1.D0)
          ENDIF
          IF(P.GT.0.D0) Q=-Q
          P=DABS(P)
          IF(2.D0*P .LT. DMIN1(3.D0*XM*Q-DABS(TOL1*Q),DABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(DABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
!      PAUSE 'ZBREND exceeding maximum iterations.'
      write (*,*) 'ZBREND exceeding maximum iterations.'
      read (*,*)
      ZBREND=B
      RETURN
      END
