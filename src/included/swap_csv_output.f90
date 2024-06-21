module SWAP_csv_output

   use variables, only: igsnow,igird,iintc,irunon,iruno,ipeva,ievap,iqdra,iqbot,                                                 &
                        gwl,pond,iptra,iqrot,iqreddry,iqredwet,iqredsol,iqredfrs,tsum,dvs,pgasspot,pgass,                        &
                        cwdmpot,cwdm,wsopot,wso,wlvpot,wlv,wstpot,wst,wrtpot,wrt,dwso,dwlv,dwlvpot,dwst,dwstpot,dwrt,dwrtpot,    &
                        ch,cf,laipot,lai,rdpot,rd,tagppot,tagp,tagptpot,tagpt,cuptgrazpot,cuptgraz,plossdm,lossdm,               &
                        imsqprec,imsqirrig,imsqbot,imsqdra,imdectot,imrottot,sampro,solbal,                                      &
                        wc10,Runoff_CN,iqtdo,iqtup,iqinfmax,TeTop,TeBot,ies0,iet0,iew0,inrai,inird,volact,ssnow,iqssdi,          &
                        flprintshort, date, t1900, dz, numnod, zbotcp, ztopcp, H, theta, K, Tsoil, cml, cmsy, inq, inqrot,       &
                        igrai, isnrai, iqmpoutdrrap, irunocn, isubl, c_top, heacap, heacon, inqssdi, inqdra, nrlevs, frarmtrx,   &
                        iqdo, iqup, inqdra_in, inqdra_out, pathwork, outfil, project, InList_csv, macp, madr

   implicit none

   private
   public :: csv_out

   ! local list
   integer,              parameter           :: M      = 103            ! total number of variables present in predefined list
   integer,              parameter           :: Mnodes =  10            ! Maximum number of nodes (Mnodes) or sub-regions (Mnodes/2; since lower and upper node are needed); MUST BE EVEN !!!

   type output
      character(len=12), dimension(M)        :: name                    ! base name of variable
      character(len=12), dimension(M)        :: unit                    ! units of variable
      character(len=24), dimension(Mnodes,M) :: head                    ! specific header name of variable (when necessary: including depth/node number between [])
      integer,           dimension(M)        :: iyes                    ! yes (1) or no (1) variable is requested for output
      real(8),           dimension(Mnodes,M) :: value                   ! value of variable; up to maximum of nodes [Mnodes] or in case of sub-regions [Mnodes/2]
      integer,           dimension(M)        :: Nnodes                  ! actual number of nodes/subregions; if non-nodal information then Nnodes = 1
      integer,           dimension(Mnodes,M) :: nodes                   ! node numbers for [] output (nodes or sub-region)
   end type output

   type local_pointer
      logical                        :: fldo
      integer                        :: jpos
      real(8), dimension(Mnodes)     :: vals
   end type local_pointer

   integer,              parameter           :: Msn    =   7            ! Maximum number of shortnames
   type shortnames
      integer                                :: nsn
      integer,           dimension(Msn)      :: n
      character(len=20), dimension(Msn)      :: alias
      character(len=20), dimension(Msn,M)    :: names
   end type shortnames

   ! all variables related information
   type(output),                save :: vars

   ! short names for input
   type(shortnames),            save :: shorts

   ! local pointers: per compartment
   type(local_pointer),         save :: lp_H, lp_WC, lp_K
   type(local_pointer),         save :: lp_C, lp_CA, lp_O2
   type(local_pointer),         save :: lp_T, lp_HCA, lp_HCO
   type(local_pointer),         save :: lp_DRA, lp_RWU, lp_Q, lp_SI
   ! local pointers: sub-regions
   type(local_pointer),         save :: lp_WTOT, lp_QTRA, lp_QTOP, lp_QBOT, lp_QDRA
   type(local_pointer),         save :: lp_QTIN, lp_QTOU, lp_QBIN, lp_QBOU
   type(local_pointer),         save :: lp_QDIN, lp_QDOU

   ! local help pointers
   real(8)                           :: dstor, baldev
   ! other local help
   real(8),                     save :: VolOld, PondOld, SnowOld

   ! user input names
   character(len=2096)               :: InList
   character(len=512), dimension(M)  :: userlist
   integer                           :: nlist

   ! dummy
   integer                           :: i

   ! initial; allowed list of input and units that apply; this sequence determines sequence in output file
   ! MUST be uppercase
   ! Sequence here determines sequence in output file
   data (vars%name(i),  vars%unit(i), i = 1, M) /&
         'RAIN',        '(cm)',        &
         'RAIN_NET',    '(cm)',        &
         'SNOW',        '(cm)',        &
         'IRRIG',       '(cm)',        &
         'IRRIG_NET',   '(cm)',        &
         'INTERC',      '(cm)',        &
         'RUNON',       '(cm)',        &
         'RUNOFF',      '(cm)',        &
         'EPOT',        '(cm)',        &
         'EACT',        '(cm)',        &
         'SUBLIM',      '(cm)',        &
         'DRAINAGE',    '(cm)',        &
         'QBOTTOM',     '(cm)',        &
         'GWL',         '(cm)',        &
         'POND',        '(cm)',        &
         'SSNOW',       '(cm)',        &
         'TPOT',        '(cm)',        &
         'TACT',        '(cm)',        &
         'TREDDRY',     '(cm)',        &
         'TREDWET',     '(cm)',        &
         'TREDSOL',     '(cm)',        &
         'TREDFRS',     '(cm)',        &
         'ES0',         '(cm)',        &
         'ET0',         '(cm)',        &
         'EW0',         '(cm)',        &
         'DSTOR',       '(cm)',        &
         'BALDEV',      '(cm)',        &
         'VOLACT',      '(cm)',        &
         'QSSDI',       '(cm)',        &
         'TSUM',        '(deg C)',     &
         'DVS',         '(-)',         &
         'PGASSPOT',    '(kgch/ha)',   &
         'PGASS',       '(kgch/ha)',   &
         'CPWDM',       '(kg/ha)',     &
         'CWDM',        '(kg/ha)',     &
         'CPWSO',       '(kg/ha)',     &
         'CWSO',        '(kg/ha)',     &
         'PWLV',        '(kg/ha)',     &
         'WLV',         '(kg/ha)',     &
         'PWST',        '(kg/ha)',     &
         'WST',         '(kg/ha)',     &
         'PWRT',        '(kg/ha)',     &
         'WRT',         '(kg/ha)',     &
         'DWSO',        '(kg/ha)',     &
         'DWLV',        '(kg/ha)',     &
         'DWLVPOT',     '(kg/ha)',     &
         'DWST',        '(kg/ha)',     &
         'DWSTPOT',     '(kg/ha)',     &
         'DWRT',        '(kg/ha)',     &
         'DWRTPOT',     '(kg/ha)',     &
         'HEIGHT',      '(cm)',        &
         'CRPFAC',      '(-)',         &
         'LAIPOT',      '(m2/m2)',     &
         'LAI',         '(m2/m2)',     &
         'RDPOT',       '(cm)',        &
         'RD',          '(cm)',        &
         'PGRASSDM',    '(kg/ha)',     &
         'GRASSDM',     '(kg/ha)',     &
         'PMOWDM',      '(kg/ha)',     &
         'MOWDM',       '(kg/ha)',     &
         'PGRAZDM',     '(kg/ha)',     &
         'GRAZDM',      '(kg/ha)',     &
         'PLOSSDM',     '(kg/ha)',     &
         'LOSSDM',      '(kg/ha)',     &
         'SQPREC',      '(g/cm2)',     &
         'SQIRRIG',     '(g/cm2)',     &
         'SQBOT',       '(g/cm2)',     &
         'SQDRA',       '(g/cm2)',     &
         'DECTOT',      '(g/cm2)',     &
         'ROTTOT',      '(g/cm2)',     &
         'SAMPRO',      '(g/cm2)',     &
         'SOLBAL',      '(g/cm2)',     &
         'WC10',        '(cm3/cm3)',   &
         'RUNOFFCN',    '(cm) ',       &
         'QTOPIN',      '(cm) ',       &
         'QTOPOUT',     '(cm) ',       &
         'QINFMAX',     '(cm)',        &
         'H[',          '(cm)',        &
         'WC[',         '(cm3/cm3)',   &
         'TEMP[',       '(deg C)',     &
         'K[',          '(cm/d)',      &
         'CONC[',       '(g/cm3 w)',   &
         'CONCADS[',    '(g/cm3)',     &
         'O2TOP[',      '(kg/m3)',     &
         'HEACAP[',     '(J/cm3/K)',   &
         'HEACON[',     '(J/cm/K/d)',  &
         'TETOP',       '(deg C)',     &
         'TEBOT',       '(deg C)',     &
         'DRAIN[',      '(cm)',        &
         'RWU[',        '(cm)',        &
         'FLUX[',       '(cm)',        &
         'SSDI[',       '(cm)',        &
         'WTOT[',       '(cm)',        &
         'QTRANS[',     '(cm)',        &
         'QTOP[',       '(cm)',        &
         'QBOT[',       '(cm)',        &
         'QDRA[',       '(cm)',        &
         'QTOPIN[',     '(cm)',        &
         'QTOPOUT[',    '(cm)',        &
         'QBOTIN[',     '(cm)',        &
         'QBOTOUT[',    '(cm)',        &
         'QDRAININ[',   '(cm)',        &
         'QDRAINOUT[',  '(cm)'         &
         /

   ! must take care on number of items provided
   data shorts%nsn           /7/
   data shorts%n(1:7)        /17, 8, 8, 33, 21, 5, 11/                                    ! shorts%nsn values
   data shorts%alias(1:7)    /"WATBAL", "SOLBAL", "ETTERMS", "CROP", "GRASS",       &     ! shorts%nsn values
                              "SUBREG_MIN", "SUBREG_ALL" /
   data shorts%names(1,1:17) /"RAIN", "SNOW", "INTERC", "RUNON", "RUNOFF",          &     ! shorts%n(1) values
                              "IRRIG", "TPOT", "TACT", "EPOT", "EACT",              &
                              "SUBLIM", "QSSDI", "QBOTTOM", "DRAINAGE", "VOLACT",   &
                              "DSTOR", "BALDEV"/
   data shorts%names(2,1:8)  /"SQPREC", "SQIRRIG", "SQBOT", "DECTOT", "ROTTOT",     &     ! shorts%n(2) values
                              "SQDRA", "SAMPRO", "SOLBAL"/
   data shorts%names(3,1:8)  /"TPOT", "TACT", "TREDDRY", "TREDWET", "TREDSOL",      &     ! shorts%n(3) values
                              "TREDFRS", "EPOT", "EACT"/
   data shorts%names(4,1:33) /"TPOT", "TACT", "TREDDRY", "TREDWET", "TREDSOL",      &     ! shorts%n(4) values
                              "TREDFRS", "TSUM", "DVS", "PGASSPOT", "PGASS",        &
                              "CPWDM", "CWDM", "CPWSO", "CWSO", "PWLV",             &
                              "WLV", "PWST", "WST", "PWRT", "WRT",                  &
                              "DWSO", "DWLV", "DWLVPOT", "DWST", "DWSTPOT",         &
                              "DWRT", "DWRTPOT", "HEIGHT", "CRPFAC", "LAIPOT",      &
                              "LAI", "RDPOT", "RD"/
   data shorts%names(5,1:21) /"TPOT", "TACT", "TREDDRY", "TREDWET", "TREDSOL",      &     ! shorts%n(5) values
                              "TREDFRS", "TSUM", "HEIGHT", "CRPFAC", "LAIPOT",      &
                              "LAI", "RDPOT", "RD", "PGRASSDM", "GRASSDM",          &
                              "PMOWDM", "MOWDM", "PGRAZDM", "GRAZDM", "PLOSSDM",    &
                              "LOSSDM"/
   data shorts%names(6,1:5)  /"WTOT", "QTRANS", "QTOP", "QBOT", "QDRA"/                   ! shorts%n(6) values
   data shorts%names(7,1:11) /"WTOT", "QTRANS", "QTOP", "QBOT", "QDRA",             &     ! shorts%n(7) values
                              "QTOPIN", "QTOPOUT", "QBOTIN", "QBOTOUT", "QDRAININ", &
                              "QDRAINOUT"/

   contains

   subroutine set_values
   implicit none
   integer :: i
   ! this must be hard-coded
   ! this provides the actual link between user-supplied name and internal variable name (MUST be uppercase)
   do concurrent (i = 1:M, vars%iyes(i) == 1)
      ! single valued vars
      if (vars%name(i) == 'RAIN')        vars%value(1,i) = igrai + isnrai
      if (vars%name(i) == 'RAIN_NET')    vars%value(1,i) = inrai
      if (vars%name(i) == 'SNOW')        vars%value(1,i) = igsnow
      if (vars%name(i) == 'IRRIG')       vars%value(1,i) = igird
      if (vars%name(i) == 'IRRIG_NET')   vars%value(1,i) = inird
      if (vars%name(i) == 'INTERC')      vars%value(1,i) = iintc
      if (vars%name(i) == 'RUNON')       vars%value(1,i) = irunon
      if (vars%name(i) == 'RUNOFF')      vars%value(1,i) = iruno
      if (vars%name(i) == 'EPOT')        vars%value(1,i) = ipeva
      if (vars%name(i) == 'EACT')        vars%value(1,i) = ievap
      if (vars%name(i) == 'SUBLIM')      vars%value(1,i) = isubl
      if (vars%name(i) == 'DRAINAGE')    vars%value(1,i) = iqdra + iQMpOutDrRap
      if (vars%name(i) == 'QBOTTOM')     vars%value(1,i) = iqbot
      if (vars%name(i) == 'GWL')         vars%value(1,i) = gwl
      if (vars%name(i) == 'POND')        vars%value(1,i) = pond
      if (vars%name(i) == 'TPOT')        vars%value(1,i) = iptra
      if (vars%name(i) == 'TACT')        vars%value(1,i) = iqrot
      if (vars%name(i) == 'TREDDRY')     vars%value(1,i) = iqreddry
      if (vars%name(i) == 'TREDWET')     vars%value(1,i) = iqredwet
      if (vars%name(i) == 'TREDSOL')     vars%value(1,i) = iqredsol
      if (vars%name(i) == 'TREDFRS')     vars%value(1,i) = iqredfrs
      if (vars%name(i) == 'ES0')         vars%value(1,i) = ies0
      if (vars%name(i) == 'ET0')         vars%value(1,i) = iet0
      if (vars%name(i) == 'EW0')         vars%value(1,i) = iew0
      if (vars%name(i) == 'DSTOR')       vars%value(1,i) = dstor
      if (vars%name(i) == 'BALDEV')      vars%value(1,i) = baldev
      if (vars%name(i) == 'VOLACT')      vars%value(1,i) = volact
      if (vars%name(i) == 'SSNOW')       vars%value(1,i) = ssnow
      if (vars%name(i) == 'QSSDI')       vars%value(1,i) = iqssdi
      if (vars%name(i) == 'TSUM')        vars%value(1,i) = tsum
      if (vars%name(i) == 'DVS')         vars%value(1,i) = dvs
      if (vars%name(i) == 'PGASSPOT')    vars%value(1,i) = pgasspot
      if (vars%name(i) == 'PGASS')       vars%value(1,i) = pgass
      if (vars%name(i) == 'CPWDM')       vars%value(1,i) = cwdmpot
      if (vars%name(i) == 'CWDM')        vars%value(1,i) = cwdm
      if (vars%name(i) == 'CPWSO')       vars%value(1,i) = wsopot
      if (vars%name(i) == 'CWSO')        vars%value(1,i) = wso
      if (vars%name(i) == 'PWLV')        vars%value(1,i) = wlvpot
      if (vars%name(i) == 'WLV')         vars%value(1,i) = wlv
      if (vars%name(i) == 'PWST')        vars%value(1,i) = wstpot
      if (vars%name(i) == 'WST')         vars%value(1,i) = wst
      if (vars%name(i) == 'PWRT')        vars%value(1,i) = wrtpot
      if (vars%name(i) == 'WRT')         vars%value(1,i) = wrt
      if (vars%name(i) == 'DWSO')        vars%value(1,i) = dwso
      if (vars%name(i) == 'DWLV')        vars%value(1,i) = dwlv
      if (vars%name(i) == 'DWLVPOT')     vars%value(1,i) = dwlvpot
      if (vars%name(i) == 'DWST')        vars%value(1,i) = dwst
      if (vars%name(i) == 'DWSTPOT')     vars%value(1,i) = dwstpot
      if (vars%name(i) == 'DWRT')        vars%value(1,i) = dwrt
      if (vars%name(i) == 'DWRTPOT')     vars%value(1,i) = dwrtpot
      if (vars%name(i) == 'HEIGHT')      vars%value(1,i) = ch
      if (vars%name(i) == 'CRPFAC')      vars%value(1,i) = cf
      if (vars%name(i) == 'LAIPOT')      vars%value(1,i) = laipot
      if (vars%name(i) == 'LAI')         vars%value(1,i) = lai
      if (vars%name(i) == 'RDPOT')       vars%value(1,i) = rdpot
      if (vars%name(i) == 'RD')          vars%value(1,i) = rd
      if (vars%name(i) == 'PGRASSDM')    vars%value(1,i) = tagppot
      if (vars%name(i) == 'GRASSDM')     vars%value(1,i) = tagp
      if (vars%name(i) == 'PMOWDM')      vars%value(1,i) = tagptpot
      if (vars%name(i) == 'MOWDM')       vars%value(1,i) = tagpt
      if (vars%name(i) == 'PGRAZDM')     vars%value(1,i) = cuptgrazpot
      if (vars%name(i) == 'GRAZDM')      vars%value(1,i) = cuptgraz
      if (vars%name(i) == 'PLOSSDM')     vars%value(1,i) = plossdm
      if (vars%name(i) == 'LOSSDM')      vars%value(1,i) = lossdm
      if (vars%name(i) == 'SQPREC')      vars%value(1,i) = imsqprec
      if (vars%name(i) == 'SQIRRIG')     vars%value(1,i) = imsqirrig
      if (vars%name(i) == 'SQBOT')       vars%value(1,i) = imsqbot
      if (vars%name(i) == 'SQDRA')       vars%value(1,i) = imsqdra
      if (vars%name(i) == 'DECTOT')      vars%value(1,i) = imdectot
      if (vars%name(i) == 'ROTTOT')      vars%value(1,i) = imrottot
      if (vars%name(i) == 'SAMPRO')      vars%value(1,i) = sampro
      if (vars%name(i) == 'SOLBAL')      vars%value(1,i) = solbal
      if (vars%name(i) == 'WC10')        vars%value(1,i) = wc10
      if (vars%name(i) == 'RUNOFFCN')    vars%value(1,i) = Runoff_CN
      if (vars%name(i) == 'QTOPIN')      vars%value(1,i) = iqtdo
      if (vars%name(i) == 'QTOPOUT')     vars%value(1,i) = iqtup
      if (vars%name(i) == 'QINFMAX')     vars%value(1,i) = iqinfmax
      if (vars%name(i) == 'TETOP')       vars%value(1,i) = TeTop
      if (vars%name(i) == 'TEBOT')       vars%value(1,i) = TeBot

      ! vars per layer
      if (vars%name(i) == 'H[')          vars%value(1:Mnodes,i) = lp_H%vals(1:Mnodes)
      if (vars%name(i) == 'WC[')         vars%value(1:Mnodes,i) = lp_WC%vals(1:Mnodes)
      if (vars%name(i) == 'TEMP[')       vars%value(1:Mnodes,i) = lp_T%vals(1:Mnodes)
      if (vars%name(i) == 'K[')          vars%value(1:Mnodes,i) = lp_K%vals(1:Mnodes)
      if (vars%name(i) == 'CONC[')       vars%value(1:Mnodes,i) = lp_C%vals(1:Mnodes)
      if (vars%name(i) == 'CONCADS[')    vars%value(1:Mnodes,i) = lp_CA%vals(1:Mnodes)
      if (vars%name(i) == 'O2TOP[')      vars%value(1:Mnodes,i) = lp_O2%vals(1:Mnodes)
      if (vars%name(i) == 'HEACAP[')     vars%value(1:Mnodes,i) = lp_HCA%vals(1:Mnodes)
      if (vars%name(i) == 'HEACON[')     vars%value(1:Mnodes,i) = lp_HCO%vals(1:Mnodes)
      if (vars%name(i) == 'DRAIN[')      vars%value(1:Mnodes,i) = lp_DRA%vals(1:Mnodes)
      if (vars%name(i) == 'RWU[')        vars%value(1:Mnodes,i) = lp_RWU%vals(1:Mnodes)
      if (vars%name(i) == 'FLUX[')       vars%value(1:Mnodes,i) = lp_Q%vals(1:Mnodes)
      if (vars%name(i) == 'SSDI[')       vars%value(1:Mnodes,i) = lp_SI%vals(1:Mnodes)

      ! values for subregions
      if (vars%name(i) == 'WTOT[')       vars%value(1:Mnodes,i) = lp_WTOT%vals(1:Mnodes)
      if (vars%name(i) == 'QTRANS[')     vars%value(1:Mnodes,i) = lp_QTRA%vals(1:Mnodes)
      if (vars%name(i) == 'QTOP[')       vars%value(1:Mnodes,i) = lp_QTOP%vals(1:Mnodes)
      if (vars%name(i) == 'QBOT[')       vars%value(1:Mnodes,i) = lp_QBOT%vals(1:Mnodes)
      if (vars%name(i) == 'QDRA[')       vars%value(1:Mnodes,i) = lp_QDRA%vals(1:Mnodes)
      if (vars%name(i) == 'QTOPIN[')     vars%value(1:Mnodes,i) = lp_QTIN%vals(1:Mnodes)
      if (vars%name(i) == 'QTOPOUT[')    vars%value(1:Mnodes,i) = lp_QTOU%vals(1:Mnodes)
      if (vars%name(i) == 'QBOTIN[')     vars%value(1:Mnodes,i) = lp_QBIN%vals(1:Mnodes)
      if (vars%name(i) == 'QBOTOUT[')    vars%value(1:Mnodes,i) = lp_QBOU%vals(1:Mnodes)
      if (vars%name(i) == 'QDRAININ[')   vars%value(1:Mnodes,i) = lp_QDIN%vals(1:Mnodes)
      if (vars%name(i) == 'QDRAINOUT[')  vars%value(1:Mnodes,i) = lp_QDOU%vals(1:Mnodes)

   end do
   end subroutine set_values

   subroutine csv_out(iTask)
   implicit none

   ! global
   integer, intent(in)  :: iTask

   ! local
   integer              :: j, il, iuncsv
   character(len=19)    :: datexti
   character(len=1024)  :: line
   character(len=300)   :: filcsv
   character(len=20)    :: form
   character(len=30)    :: myval
   ! functions
   integer              :: getun
   save                 :: iuncsv

   select case (iTask)
   case (1)
      ! initial values in vars
      vars%iyes(1:M)           =    0
      vars%value(1:Mnodes,1:M) = -999.9d0
      vars%Nnodes(1:M)         =    1
      vars%nodes(1:Mnodes,1:M) =    0
      
      ! Set format for output
      call what_form(1)

      ! to be based on user input information
!      InList = 'H[-10.0,4,5], RAIN, WC[2,3,6], H[4], WTOT[0:-15.0,4:6], Tact, WC[1], WC[-2], WC[1], rain,GWL, RWU[1,2,3,4,5,6], QTRANS[1:3,4:6], qtopin[0:-30]'
      Inlist = InList_csv

      ! extract all individual vars from InList; remove double entries; merge same vars with [] in single var with []; determine which vars are require
      call make_userlist(InList)
      call merge
      call det_which_vars

      ! output file; write header
      iuncsv = getun (500, 900)
      filcsv = trim(pathwork)//trim(outfil)//'_output.csv'
      open(unit=iuncsv, file=filcsv, status="unknown")
      call makeheader(iuncsv, filcsv)

      ! store inital values
      VolOld  = volact
      PondOld = pond
      SnowOld = ssnow

   case (2)
      ! dynamic
      ! some local pointers must be filled here first; then values are transferred to vars%values
      call fill_values
      call set_values

      ! line contains results in comma-separated format; il is its length
      ! time is first value
      line = ""; il = 0
      if (.not. flprintshort) then
         call addstr(line, il, trim(date)); call addstr(line, il, ",")
      else
         ! determine date-time
         call dtdpst ('year-month-day hour:minute:seconds', t1900, datexti)
         call addstr(line, il, trim(datexti)); call addstr(line, il, ",")
      end if

      ! next: add all desired values
      do j = 1, M
         if (vars%iyes(j) == 1) then
            do i = 1, vars%Nnodes(j)
               call what_form(2, vars%value(i,j), form)
               write (myval, form) vars%value(i,j)
               line = line(1:il) // trim(adjustl(myval)); il = len_trim(line)
            end do
         end if
      end do

      ! write result (skip last character which is a comma)
      write(iuncsv,'(A)') line(1:il-1)

   case (3)
      close (unit=iuncsv)

   case default
      call fatalerr("csv_out","Illegal iTask")

   end select

   end subroutine csv_out

   subroutine what_form(iTask, var, form)
   implicit none
   ! global
   integer,          intent(in)            :: iTask
   real(8),          intent(in),  optional :: var
   character(len=*), intent(out), optional :: form
   ! local
   ! when to automatically swith from F to E formatting
   integer,           parameter           :: num_d = 5         ! # of decimals; later: user input?
   integer,           parameter           :: num_w = num_d+7   ! total width of format, for E-formatting: 7 positions are needed for "-x."at start and "E+00" at end
   integer,           parameter           :: expo = 4
   real(8),           parameter           :: t1 = 1.0d0/(10.0d0**expo)
   real(8),           parameter           :: t2 = 10.0d0**expo
   character(len=2)                       :: cval1, cval2
   character(len=20),           save      :: form_rea_F1, form_rea_E1,form_rea_F0

   if (iTask == 1) then
      ! some basic info
      !  Note: field width of zero in I and F edit descriptors is allowed as of Fortran95 to ensure as little space usage
      !        in output file as possible (Metcalf et al., 2004, Fortran 95/2003 explained, Oxford Univ. Prrss, p. 199)
      write (cval1,'(I0)') num_d
      write (cval2,'(I0)') num_w
      form_rea_F0 = '(F0.0,",")'
      form_rea_F1 = '(F0.' // trim(cval1) // ',",")'
      form_rea_E1 = '(1P,E' // trim(cval2) // '.' // trim(cval1) // ',",")'
   else
      if (abs(var) < 1.0D-10) then
         form = form_rea_F0
      else if (abs(var) < t1 .OR. abs(var) > t2) then
         form = form_rea_E1
      else
         form = form_rea_F1
      end if
   end if
   end subroutine what_form

   subroutine fill_values
   implicit none
   integer :: i

   ! change in storage and deviation in mass balance
   dstor        = (volact + pond + ssnow) - (VolOld + PondOld + SnowOld) 
   baldev       = (igrai+isnrai+igsnow+igird+irunon + iqssdi) - dstor -             &
                  (iintc+iruno+irunoCN+iqrot+ievap+isubl+iQMpOutDrRap+iqdra+(-1.0d0*iqbot))

   ! replace Old values by current values (needed for next output moment)
   VolOld  = volact
   PondOld = pond
   SnowOld = ssnow

   ! H, WC, TEMP, K, CONC, CONCADS, O2TOP
   if (lp_H%fldo)   lp_H%vals(1:vars%Nnodes(lp_H%jpos))     =       H(vars%nodes(1:vars%Nnodes(lp_H%jpos), lp_H%jpos))
   if (lp_WC%fldo)  lp_WC%vals(1:vars%Nnodes(lp_WC%jpos))   =   theta(vars%nodes(1:vars%Nnodes(lp_WC%jpos),lp_WC%jpos))
   if (lp_K%fldo)   lp_K%vals(1:vars%Nnodes(lp_K%jpos))     =       K(vars%nodes(1:vars%Nnodes(lp_K%jpos), lp_K%jpos))
   if (lp_T%fldo)   lp_T%vals(1:vars%Nnodes(lp_T%jpos))     =   Tsoil(vars%nodes(1:vars%Nnodes(lp_T%jpos), lp_T%jpos))
   if (lp_C%fldo)   lp_C%vals(1:vars%Nnodes(lp_C%jpos))     =     cml(vars%nodes(1:vars%Nnodes(lp_C%jpos), lp_C%jpos))
   if (lp_CA%fldo)  lp_CA%vals(1:vars%Nnodes(lp_CA%jpos))   =    cmsy(vars%nodes(1:vars%Nnodes(lp_CA%jpos),lp_CA%jpos))
   if (lp_O2%fldo)  lp_O2%vals(1:vars%Nnodes(lp_O2%jpos))   =   c_top(vars%nodes(1:vars%Nnodes(lp_O2%jpos),lp_O2%jpos))

   ! HEACAP, HEACON, DRAIN, RWU, FLUX, SSDI
   if (lp_HCA%fldo) lp_HCA%vals(1:vars%Nnodes(lp_HCA%jpos)) =  HEACAP(vars%nodes(1:vars%Nnodes(lp_HCA%jpos),lp_HCA%jpos))
   if (lp_HCO%fldo) lp_HCO%vals(1:vars%Nnodes(lp_HCO%jpos)) =  HEACON(vars%nodes(1:vars%Nnodes(lp_HCO%jpos),lp_HCO%jpos))
   if (lp_RWU%fldo) lp_RWU%vals(1:vars%Nnodes(lp_RWU%jpos)) =  inqrot(vars%nodes(1:vars%Nnodes(lp_RWU%jpos),lp_RWU%jpos))
   if (lp_Q%fldo)   lp_Q%vals(1:vars%Nnodes(lp_Q%jpos))     =     inq(vars%nodes(1:vars%Nnodes(lp_Q%jpos),  lp_Q%jpos))
   if (lp_SI%fldo)  lp_SI%vals(1:vars%Nnodes(lp_SI%jpos))   = inqssdi(vars%nodes(1:vars%Nnodes(lp_SI%jpos), lp_SI%jpos))

   ! summed for all drainage levels
   if (lp_DRA%fldo) then
      do i = 1, vars%Nnodes(lp_DRA%jpos); lp_DRA%vals(i) = sum(inqdra(1:nrlevs,i)); end do
   end if

   ! handle subregions: WTOT, QTRANS, QTOP, QBOT, QDRA, lp_QTIN, lp_QTOU, lp_QBIN, lp_QBOU, lp_QDIN, lp_QDOU
   if (lp_WTOT%fldo) call fill_1(vars%Nnodes(lp_WTOT%jpos), lp_WTOT%jpos, vars%nodes, lp_WTOT%vals, 1, theta, dz, FrArMtrx)

   ! summed sink-terms inside subregion
   if (lp_QTRA%fldo) call fill_1(vars%Nnodes(lp_QTRA%jpos), lp_QTRA%jpos, vars%nodes, lp_QTRA%vals, 0, inqrot)
   if (lp_QDRA%fldo) call fill_2(vars%Nnodes(lp_QDRA%jpos), lp_QDRA%jpos, vars%nodes, lp_QDRA%vals, nrlevs, inqdra)
   if (lp_QDIN%fldo) call fill_2(vars%Nnodes(lp_QDIN%jpos), lp_QDIN%jpos, vars%nodes, lp_QDIN%vals, nrlevs, inqdra_in)
   if (lp_QDOU%fldo) call fill_2(vars%Nnodes(lp_QDOU%jpos), lp_QDOU%jpos, vars%nodes, lp_QDOU%vals, nrlevs, inqdra_out)

   ! fluxes at top and bottom boudanries of subregions
   if (lp_QTOP%fldo) call fill_3(vars%Nnodes(lp_QTOP%jpos), lp_QTOP%jpos, 1, 0, vars%nodes, lp_QTOP%vals, -inq)
   if (lp_QBOT%fldo) call fill_3(vars%Nnodes(lp_QBOT%jpos), lp_QBOT%jpos, 0, 1, vars%nodes, lp_QBOT%vals, -inq)
   if (lp_QTIN%fldo) call fill_3(vars%Nnodes(lp_QTIN%jpos), lp_QTIN%jpos, 1, 0, vars%nodes, lp_QTIN%vals, iqdo)
   if (lp_QTOU%fldo) call fill_3(vars%Nnodes(lp_QTOU%jpos), lp_QTOU%jpos, 1, 0, vars%nodes, lp_QTOU%vals, iqup)
   if (lp_QBIN%fldo) call fill_3(vars%Nnodes(lp_QBIN%jpos), lp_QBIN%jpos, 0, 1, vars%nodes, lp_QBIN%vals, iqup)
   if (lp_QBOU%fldo) call fill_3(vars%Nnodes(lp_QBOU%jpos), lp_QBOU%jpos, 0, 1, vars%nodes, lp_QBOU%vals, iqdo)

   end subroutine fill_values

   subroutine fill_1(Nnodes, jpos, nodes, vals, iwtot, swap_vals, dz, fr)
   ! global
   integer,                        intent(in)            :: Nnodes, jpos, iwtot
   integer, dimension(Mnodes,M),   intent(in)            :: nodes
   real(8), dimension(Nnodes),     intent(out)           :: vals
   real(8), dimension(macp),       intent(in)            :: swap_vals
   real(8), dimension(macp),       intent(in), optional  :: dz, fr
   ! local
   integer                                               :: i, j
   do i = 1, Nnodes
      vals(i) = 0.0d0
      do j = nodes(2*i-1,jpos), nodes(2*i,jpos)
         if (iwtot /= 0) then
            vals(i) = vals(i) + swap_vals(j)*dz(j)*fr(j)
         else
            vals(i) = vals(i) + swap_vals(j)
         end if
      end do
   end do
   end subroutine fill_1

   subroutine fill_2(Nnodes, jpos, nodes, vals, nrlevs, swap_vals)
   ! global
   integer,                        intent(in)   :: Nnodes, jpos, nrlevs
   integer, dimension(Mnodes,M),   intent(in)   :: nodes
   real(8), dimension(Nnodes),     intent(out)  :: vals
   real(8), dimension(madr,macp),  intent(in)   :: swap_vals
   ! local
   integer                                      :: i, j
   do i = 1, Nnodes
      vals(i) = 0.0d0
      do j = nodes(2*i-1,jpos), nodes(2*i,jpos)
         vals(i) = vals(i) + sum(swap_vals(1:nrlevs,j))
      end do
   end do
   end subroutine fill_2

   subroutine fill_3(Nnodes, jpos, low, raise, nodes, vals, swap_vals)
   ! global
   integer,                        intent(in)   :: Nnodes, jpos, low, raise
   integer, dimension(Mnodes,M),   intent(in)   :: nodes
   real(8), dimension(Nnodes),     intent(out)  :: vals
   real(8), dimension(macp),       intent(in)   :: swap_vals
   ! local
   integer                                      :: i, j
   do i = 1, Nnodes
      j = nodes(2*i-low,jpos)
      vals(i) = swap_vals(j+raise)
   end do
   end subroutine fill_3

   subroutine makeheader(iuncsv, filcsv)
   implicit none
   integer, intent(in) :: iuncsv
   character(len=*), intent(in) :: filcsv
   character(len=1024) :: header_units, header_names
   integer :: il1, il2, j, k

   ! write universal SWAP header
   call writehead (iuncsv, 1, filcsv, 'specified output data of SWAP', project)

   ! typical csv header lines: (units, names)
   header_units = "* (d)"; il1 = len_trim(header_units)
   header_names = "DATETIME";  il2 = len_trim(header_names)
   do j = 1, M
      if (vars%iyes(j) == 1) then
         do k = 1, vars%Nnodes(j)
            call addstr(header_units, il1, "," // vars%unit(j))
            call addstr(header_names, il2, "," // vars%head(k,j))
         end do
      end if
   end do
   write(iuncsv,'(A)') trim(header_units)
   write(iuncsv,'(A)') trim(header_names)
   end subroutine makeheader

   subroutine det_which_vars()
   implicit none
   ! global
   ! local
   integer :: i, j, ipos1, ipos2
   vars%iyes = 0
   do i = 1, nlist
      do j = 1, M
         ipos1 = index(userlist(i), "[")
         if (ipos1 == 0) then
            if (userlist(i) == vars%name(j)) then
               vars%iyes(j)   = 1
               vars%head(1,j) = vars%name(j)
            end if
         else
            if (userlist(i)(1:ipos1) == vars%name(j)) then
               vars%iyes(j) = 1
               ipos2 = index(userlist(i), "]")
               call handle_sqr_brackets(i, j, ipos1, ipos2)
            end if
         end if
      end do
   end do
   end subroutine det_which_vars

   subroutine handle_sqr_brackets(i, j, ipos1, ipos2)
   implicit none
   ! global
   integer, intent(in)           :: i, j, ipos1, ipos2
   ! words
   integer, parameter            :: ilw = Mnodes
   integer, dimension(ilw)       :: iwbeg, iwend
   integer                       :: ifnd, iwar, k, ktel, ipos3
   character(len=*), parameter   :: colon = ":"
   character(len=*), parameter   :: semicolon = ";"
   character(len=256)            :: str, substr, name
   real                          :: rv
   logical                       :: isSUB

   isSUB = .false.
   name  = userlist(i)(1:ipos1-1)

   ! H, WC, TEMP, K, CONC, CONCADS, O2TOP
   if (name == "H")         then; lp_H%fldo    = .true.; lp_H%jpos    = j; end if
   if (name == "WC")        then; lp_WC%fldo   = .true.; lp_WC%jpos   = j; end if
   if (name == "TEMP")      then; lp_T%fldo    = .true.; lp_T%jpos    = j; end if
   if (name == "K")         then; lp_K%fldo    = .true.; lp_K%jpos    = j; end if
   if (name == "CONC")      then; lp_C%fldo    = .true.; lp_C%jpos    = j; end if
   if (name == "CONCADS")   then; lp_CA%fldo   = .true.; lp_CA%jpos   = j; end if
   if (name == "O2TOP")     then; lp_O2%fldo   = .true.; lp_O2%jpos   = j; end if

   ! HEACAP, HEACON, DRAIN, RWU, FLUX, SSDI
   if (name == "HEACAP")    then; lp_HCA%fldo  = .true.; lp_HCA%jpos  = j; end if
   if (name == "HEACON")    then; lp_HCO%fldo  = .true.; lp_HCO%jpos  = j; end if
   if (name == "DRAIN")     then; lp_DRA%fldo  = .true.; lp_DRA%jpos  = j; end if
   if (name == "RWU")       then; lp_RWU%fldo  = .true.; lp_RWU%jpos  = j; end if
   if (name == "FLUX")      then; lp_Q%fldo    = .true.; lp_Q%jpos    = j; end if
   if (name == "SSDI")      then; lp_SI%fldo   = .true.; lp_SI%jpos   = j; end if

   ! Subregion: WTOT, QTRANS, QTOP, QBOT, QDRA, QTOPIN, QTOPOUT, QBOTIN, QBOTOUT, QDRAININ, QDRAINOUT
   if (name == "WTOT")      then; lp_WTOT%fldo = .true.; lp_WTOT%jpos = j; isSUB = .true.; end if
   if (name == "QTRANS")    then; lp_QTRA%fldo = .true.; lp_QTRA%jpos = j; isSUB = .true.; end if
   if (name == "QTOP")      then; lp_QTOP%fldo = .true.; lp_QTOP%jpos = j; isSUB = .true.; end if
   if (name == "QBOT")      then; lp_QBOT%fldo = .true.; lp_QBOT%jpos = j; isSUB = .true.; end if
   if (name == "QDRA")      then; lp_QDRA%fldo = .true.; lp_QDRA%jpos = j; isSUB = .true.; end if
   if (name == "QTOPIN")    then; lp_QTIN%fldo = .true.; lp_QTIN%jpos = j; isSUB = .true.; end if
   if (name == "QTOPOUT")   then; lp_QTOU%fldo = .true.; lp_QTOU%jpos = j; isSUB = .true.; end if
   if (name == "QBOTIN")    then; lp_QBIN%fldo = .true.; lp_QBIN%jpos = j; isSUB = .true.; end if
   if (name == "QBOTOUT")   then; lp_QBOU%fldo = .true.; lp_QBOU%jpos = j; isSUB = .true.; end if
   if (name == "QDRAININ")  then; lp_QDIN%fldo = .true.; lp_QDIN%jpos = j; isSUB = .true.; end if
   if (name == "QDRAINOUT") then; lp_QDOU%fldo = .true.; lp_QDOU%jpos = j; isSUB = .true.; end if

   str = userlist(i)(ipos1+1:ipos2-1)
   call words (str, ilw, semicolon, iwbeg, iwend, ifnd)
   if (.not. isSUB) then
      ! refers to specific nodes
      if (ifnd > Mnodes) call fatalerr("handle_sqr_brackets", "Too many nodes between [] input")
      vars%Nnodes(j) = ifnd
      do k = 1, ifnd
         call decrea (iwar, str(iwbeg(k):iwend(k)), rv)
         if (rv <= 0.0d0) rv = real(nodenumber(rv, 1))
         vars%nodes(k,j) = nint(rv)
         vars%head(k,j) = trim(name) // "[" // trim(str(iwbeg(k):iwend(k))) // "]"
      end do
      call my_piksrt(ifnd, vars%nodes(1:ifnd,j), vars%head(1:ifnd,j))
   else
      ! refers to sub parts of soil column
      if (ifnd > Mnodes/2) call fatalerr("handle_sqr_brackets", "Too many subregions between [] input")
      vars%Nnodes(j) = ifnd
      ktel = 0
      do k = 1, ifnd
         substr = trim(str(iwbeg(k):iwend(k)))
         ipos3 = index(substr, colon); if (ipos3 == 0) call fatalerr("handle_sqr_brackets", "Expecting : in subregions []")
         call decrea (iwar, substr(1:(ipos3-1)), rv)
         ktel = ktel + 1
         ! to do: if rv < 0; then determine nodenumber at given depth
         if (rv <= 0.0d0) rv = real(nodenumber(rv, 2))
         vars%nodes(ktel,j) = nint(rv)
         call decrea (iwar, substr((ipos3+1):), rv)
         ktel = ktel + 1
         if (rv <= 0.0d0) rv = real(nodenumber(rv, 1))
         vars%nodes(ktel,j) = nint(rv)
         vars%head(k,j) = trim(name) // "[" // trim(substr) // "]"
      end do
   end if
   end subroutine handle_sqr_brackets

   subroutine make_userlist(InList)
   implicit none
   ! global
   character(len=*),               intent(inout) :: InList
   ! local
   integer                   :: i, j, k, ipos1, ipos2, il
   integer, parameter        :: ilw = M
   integer, dimension(ilw)   :: iwbeg, iwend
   character(len=2096)       :: InListShort, InListTmp
   character(len=1024)       :: str_insert
   character(len=256)        :: str_sqr, str_item
   logical                   :: bracket, first
   ! function
   logical :: STRinARSTR

   ! replace , inside [] by ;
   ipos2 = 1
1  continue
   ipos1 = index(trim(InList(ipos2:)), "[") + (ipos2-1)
   if (ipos1 > ipos2) then
      ipos2 = index(trim(InList(ipos1:)),"]") + (ipos1-1)
      if (ipos2 == 0) call fatalerr("make_userlist","Unmatched []")
      do i = ipos1+1,ipos2-1
         if (InList(i:i) == ",") InList(i:i) = ";"
      end do
      go to 1
   end if
   call upperc(Inlist)

   ! replace shortnames (alias) by all accompanying names
   first = .true.
   InListShort = ""
   call words(trim(InList), ilw, ",", iwbeg, iwend, nlist)
   do i = 1, nlist
     
     str_item = Inlist(iwbeg(i):iwend(i))
     
     ! check for brackets
     bracket = .false.
     ipos1 = index(str_item, "[")
     if (ipos1 > 0) then
       bracket = .true.
       str_sqr = str_item(ipos1:len_trim(str_item))
       str_item = str_item(1:(ipos1 - 1))
     end if
     
     ! check all shortnames
     do j = 1, shorts%nsn
       
       ipos1 = index(trim(str_item), trim(shorts%alias(j)))
       if (ipos1 == 1 .and. len_trim(str_item) == len_trim(shorts%alias(j))) then
         
         iwbeg(i) = 0
           
         ! replacement string
         str_insert = ""; il = len_trim(str_insert)
         do k = 1, shorts%n(j)
           call addstr(str_insert, il, shorts%names(j,k))
           if (bracket) call addstr(str_insert, il, str_sqr)
           if (k < shorts%n(j)) call addstr(str_insert, il, ",")
         end do
         
         ! add items to InlistShort
         if (first) then
           first = .false.
           InListShort = trim(str_insert)
         else
           InListShort = trim(InListShort) // "," // trim(str_insert)
         end if
       
       end if
     end do
   end do
   
   ! combine Inlist and InlistShort (if needed)
   if (.not. first) then
     first = .true.
     InListTmp = ""
     do i = 1, nlist
       if (iwbeg(i) > 0) then
         str_item = Inlist(iwbeg(i):iwend(i))
         if (first) then
           first = .false.
           InListTmp = trim(str_item)
         else
           InListTmp = trim(InListTmp) // "," // trim(str_item)
         end if
       end if
     end do
     InList = trim(InListTmp) // "," // trim(InListShort)
   end if
   
   
   
   
!   do j = 1, shorts%nsn
!      ipos1 = index(trim(InList), trim(shorts%alias(j)))
!      if (ipos1 > 0) then
!         ! alias is there
!         ipos2 = index(trim(InList(ipos1:)), ",") + (ipos1-1); if (ipos2 < ipos1) ipos2 = len_trim(InList)
!         str_help = trim(InList(ipos1:ipos2))
!         ipos3 = index(str_help, "[")
!         if (ipos3 > 0) then
!            ! square barckets are there
!            ipos4 = index(str_help, "]")
!            str_sqr = str_help(ipos3:ipos4)
!         end if
!         ! replacement string
!         str_insert = ""; il = len_trim(str_insert)
!         do k = 1, shorts%n(j)
!            call addstr(str_insert, il, shorts%names(j,k))
!            if (ipos3 > 0) call addstr(str_insert, il, str_sqr)
!            if (k < shorts%n(j)) call addstr(str_insert, il, ",")
!         end do
!         ! replace
!         if(ipos2 == len_trim(InList)) then
!            Inlist = Inlist(1:ipos1-1) // trim(str_insert)
!         else
!            Inlist = Inlist(1:ipos1-1) // trim(str_insert) // Inlist(ipos2:)
!         end if
!      end if
!   end do

   ! replace input in single string by list of individual items
   call words(trim(InList), ilw, ",", iwbeg, iwend, nlist)
   do i = 1, nlist
      userlist(i) = trim(adjustl(InList(iwbeg(i):iwend(i))))
      if (.not. STRinARSTR(userlist(i), vars%name, M)) call fatalerr ("csv_write", "Item in userlist not known in vars%name")
   end do

   end subroutine make_userlist

   subroutine merge()
   implicit none
   ! global
   ! local
   integer                            :: i, j, k, m, ipos1, ipos2, ipos3, ifnd
   character(len=256)                 :: help, ctmp
   integer,            parameter      :: ilw = 20
   integer,            dimension(ilw) :: iwbeg, iwend

   ! remove double entries
1  continue
   do i = 1, nlist-1
      do j = i+1, nlist
         if (userlist(i) == userlist(j)) then
            do k = j, nlist-1
               userlist(k) = userlist(k+1)
            end do
            nlist = nlist - 1
            userlist(nlist+1) = ""
            go to 1
         end if
      end do
   end do

   ! merge within []
2  continue
   do i = 1, nlist-1
      ipos1 = index(userlist(i), "[")
      if (ipos1 > 0) then
         do j = i+1, nlist
            if (userlist(i)(1:ipos1) == userlist(j)(1:ipos1)) then
               ipos2 = index(userlist(i), "]")
               ipos3 = index(userlist(j), "]")
               userlist(i) = userlist(i)(1:ipos2-1) // ";" // userlist(j)(ipos1+1:ipos3)
               do k = j, nlist-1
                  userlist(k) = userlist(k+1)
               end do
               nlist = nlist - 1
               userlist(nlist+1) = ""
               go to 2
            end if
         end do
      end if
   end do

   ! check doubles inside []
   do i = 1, nlist-1
      ipos1 = index(userlist(i), "[")
      if (ipos1 > 0) then
         ipos2 = index(userlist(i), "]")
         help = userlist(i)(ipos1+1:ipos2-1)
         call words(trim(help), ilw, ";", iwbeg, iwend, ifnd)
3        continue
         do j = 1, ifnd-1
            do k = j+1, ifnd
               if (help(iwbeg(j):iwend(j)) == help(iwbeg(k):iwend(k))) then
                  do m = k, ifnd-1
                     iwbeg(m) = iwbeg(m+1)
                     iwend(m) = iwend(m+1)
                  end do
                  ifnd = ifnd - 1
                  goto 3
               end if
            end do
         end do
         ctmp = userlist(i)(1:ipos1)
         do k = 1, ifnd
            ctmp = trim(ctmp) // help(iwbeg(k):iwend(k)) // ";"
         end do
         ctmp = trim(ctmp(1:len_trim(ctmp)-1)) // "]"
         userlist(i) = trim(ctmp)
      end if
   end do

   end subroutine merge

   function nodenumber(depth, iway)
!   use variables, only : zbotcp, numnod
   implicit none
   ! global
   real,    intent(in)  :: depth
   integer, intent(in)  :: iway
   integer              :: nodenumber
   ! local
   integer              :: i

   if (iway == 1) then
      i = 1
      do while (real(zbotcp(i)) .gt. (depth + 1.0d-5))
         i = i + 1
         if (i > numnod) call fatalerr("nodenumber", "Depth > zbotcp(numnod)")
      end do
      nodenumber = i
   else
      i = 1
      do while (real(ztopcp(i)) .gt. (depth + 1.0d-5))
         i = i + 1
         if (i > numnod) call fatalerr("nodenumber", "Depth > ztopcp(numnod)")
      end do
      nodenumber = i
   end if

   end function nodenumber

end module SWAP_csv_output

module SWAP_csv_output_tz

private
public :: csv_out_tz

contains

subroutine csv_out_tz (iTask)
! Routine designed for CSV output of user-selected vaiables (provided matching defined variables in this routine).
! Specifically for selected time-depth variables
! Contains help routines: do_write_csv; check_list; remove_sqbr; det_node; Make_Header

! import global variables contianing possible output
use variables, only: pathwork, outfil, project, InList_csv_tz, tz_z1_z2, numnod, z, zbotcp, flprintshort, date, t1900, &
                     h, theta, tsoil, K, cml, cmsy, c_top, HEACAP, HEACON, inqrot

implicit none
! global
integer,          intent(in)        :: iTask

! local (some need to be saved)

! allowed variable names: number of items in Allowed must be exactly equal to Mlist
! Since Forcheck reports an error when elements have different number of characters, all now have same length
! Must be UPPERCASE
integer,                               parameter   :: Mlist   = 10
character(len=10),   dimension(Mlist)              :: Allowed
character(len=10),   dimension(Mlist)              :: Units

integer,             dimension(Mlist), save        :: iCSV, iPOS

integer,             parameter                     :: ilw = Mlist
integer,             dimension(ilw)                :: iWbeg, iWend
integer,                               save        :: iuncsv, Nvars, nod_1, nod_2
integer                                            :: i, j
character(len=*),    parameter                     :: comma   = ','
character(len=300)                                 :: filnam
character(len=160)                                 :: filtext
character(len=2)                                   :: cval1, cval2
character(len=20),                     save        :: formZ, form_rea_E
character(len=1024)                                :: Header, HeaderUnits
character(len=20),   dimension(ilw)                :: listVars
character(len=19)                                  :: datexti

! when to automatically swith from F to E formatting
integer,             parameter                     :: num_d = 5         ! # of decimals; later: user input?
integer,             parameter                     :: num_w = num_d+7   ! total width of format, for E-formatting: 7 positions are needed for "-x."at start and "E+00" at end

! functions
integer                                            :: getun

data (Allowed(i), Units(i), i = 1, Mlist) / &
      'H',        '(cm)',        &
      'WC',       '(cm3/cm3)',   &
      'TEMP',     '(deg C)',     &
      'K',        '(cm/d)',      &
      'CONC',     '(g/cm3 w)',   &
      'CONCADS',  '(g/cm3)',     &
      'O2TOP',    '(kg/m3)',     &
      'HEACAP',   '(J/m3/K)',    &
      'HEACON',   '(W/m/K)',     &
      'RWU',      '(cm/d)'       &
   /

select case (iTask)
case (1)

   ! some basic info
   !  Note: field width of zero in I and F edit descriptors is allowed as of Fortran95 to ensure as little space usage
   !        in output file as possible (Metcalf et al., 2004, Fortran 95/2003 explained, Oxford Univ. Press, p. 199)
   write (cval1,'(I0)') num_d
   write (cval2,'(I0)') num_w
   form_rea_E = '(A1,1P,E' // trim(cval2) // '.' // trim(cval1) // ')'
   formZ = '(A1,F0.3)'

   ! counter for number of records
   !Nlines = 0

   ! make user-supplied list UPPERCASE
   call upperc (InList_csv_tz)

   ! Nvars en ListVars
   call words (InList_csv_tz, ilw, ', ', iWbeg, iWend, Nvars)
   do i = 1, Nvars
      ListVars(i) = InList_csv_tz(iWbeg(i):iWend(i))
   end do

   ! check if ListVars contains valid data; determine iCSV and iPOS
   call check_list_tz ()

   ! sort position in inList_csv in same way as in Allowed; ListVars is changed accordingly
   ! this needed since sequence in output columns are fixed by appearance in Allowed
   call sort_list_tz (InList_csv_tz, ListVars, Nvars)

   ! depth interval to consider
   if (tz_z1_z2(1) > 0.0d0) then
      nod_1 = 1
      nod_2 = numnod
   else
      ! determine layer number of first z
      i = 1
      do while (zbotcp(i) .gt. (tz_z1_z2(1) + 1.0d-5))
         i = i + 1
      end do
      nod_1 = i
      ! determine layer number of second z
      i = 1
      do while (zbotcp(i) .gt. (tz_z1_z2(2) + 1.0d-5))
         i = i + 1
      end do
      nod_2 = i
   end if
   
   ! open file for output; existing file will be overwritten; formatted output
   !    to do: write some basic info at the top of the output file?
   iuncsv = getun (500, 900)
   filnam = trim(pathwork)//trim(outfil)//'_output_tz.csv'
   open (unit=iuncsv, file=filnam, status='replace', form='formatted')

   ! column header
   Header = 'DATE,DEPTH'
   HeaderUnits = '*,(cm)'
   call make_header_tz (ListVars, Nvars, Header)
   call make_headerunits_tz (HeaderUnits)
   filtext = 'specified output data of SWAP'
   call writehead (iuncsv,1,filnam,filtext,project)
   write (iuncsv,'(A)') trim(HeaderUnits)
   write (iuncsv,'(A)') trim(Header)

case (2)

  do j = nod_1, nod_2

    ! date and time
    if (.not. flprintshort) then
       write (iuncsv,'(2A)',advance='no') trim(date)
    else
       ! determine date-time
       call dtdpst ('year-month-day hour:minute:seconds',t1900,datexti)
       write (iuncsv,'(2A)',advance='no') trim(datexti)
    end if

    ! depth
    write (iuncsv,formZ,advance='no') comma, z(j)

    ! all other variables
    ! programmer is responsible for correct correspondence between Names (and their position) in Allowed and
    ! actual SWAP variables as used below

    if (iCSV(1)  == 1) call do_write_csv_tz (h(j))
    if (iCSV(2)  == 1) call do_write_csv_tz (theta(j))
    if (iCSV(3)  == 1) call do_write_csv_tz (tsoil(j))
    if (iCSV(4)  == 1) call do_write_csv_tz (k(j))
    if (iCSV(5)  == 1) call do_write_csv_tz (cml(j))
    if (iCSV(6)  == 1) call do_write_csv_tz (cmsy(j))
    if (iCSV(7)  == 1) call do_write_csv_tz (c_top(j))
    if (iCSV(8)  == 1) call do_write_csv_tz (HEACAP(j)/1.0d-6)    ! from J/cm3/K  to J/m3/K
    if (iCSV(9)  == 1) call do_write_csv_tz (HEACON(j)/864.0d0)   ! from J/cm/K/d to W/m/K
    if (iCSV(10) == 1) call do_write_csv_tz (inqrot(j))

    ! finalize record (advance to next line)
    write (iuncsv,*)

  end do

case (3)

   close (unit=iuncsv)

case default

    call fatalerr ('csv_write_tz','Illegal iTask value; range allowed: [1-3]')

end select

return

contains

subroutine do_write_csv_tz (var)
implicit none
real(8),          intent(in)  :: var

write (iuncsv,form_rea_E,advance='no') comma, real(var)

end subroutine do_write_csv_tz

subroutine check_list_tz ()
implicit none

integer :: i, j
logical :: isthere

iCSV = 0
iPOS = 0
do i = 1, Nvars
   isthere = .false.
   do j = 1, Mlist
      isthere = trim(ListVars(i)) == trim(Allowed(j))
      if (isthere) then
         iCSV(j) = 1    ! indicator for Yes/No for writing
         if (iPOS(j) > 0) call fatalerr ('check_list_tz', 'Double entries in inList_csv_tz are not allowed')
         iPOS(j) = i    ! help vector to be used later for sorting
         exit
      end if
   end do
   ! error if not there; or should we write a warning?
   if (.not. isthere) call fatalerr ('check_list_tz','Illegal variable name in inList_csv_tz for write_csv output')
end do
end subroutine check_list_tz

subroutine sort_list_tz (InList_csv, ListVars, Nvars)

implicit none
! global; NB: on output both variables are arranged in same sorted order as Allowed
integer,                            intent(in)     :: Nvars
character(len=*),                   intent(out)    :: InList_csv
character(len=*), dimension(Nvars), intent(inout)  :: ListVars
! local
integer              :: il, i, j
character(len=1024)  :: temp

il = 0
! iPOS contains for each of the Allowed names the position in the user-supplied list ListVars,
! and adds this to a temporary string
do j = 1, Mlist
   if (iPOS(j) > 0) then
      if (il == 0) then
         temp = trim(listVars(iPOS(j))); il = len_trim(temp)
      else
         temp = temp(1:il) // comma // trim(listVars(iPOS(j))); il = len_trim(temp)
      end if
   end if
end do

! temp now contains the user-supplied names in the same order of appearance is in Allowed
! decompose temp into ListVars (overwrite)
call words (temp, ilw, comma, iWbeg, iWend, Nvars)
do i = 1, Nvars
   ListVars(i) = temp(iWbeg(i):iWend(i))
end do

! on return: InList_csv becomes equal to sorted string temp
InList_csv = trim(temp)

end subroutine sort_list_tz

subroutine Make_Header_tz (ListVars, Nvars, Header)
implicit none
! global
integer,                            intent(in)     :: Nvars
character(len=*), dimension(Nvars), intent(in)     :: ListVars
character(len=*),                   intent(inout)  :: Header
! local
integer :: i, il

il = len_trim(Header)
do i = 1, Nvars
   Header = Header(1:il) // "," // trim(ListVars(i)); il = len_trim(Header)
end do

end subroutine Make_Header_tz

subroutine make_headerunits_tz (HeaderUnits)
implicit none
! global
character(len=*),                   intent(inout)  :: HeaderUnits
! local
integer :: i, il

il = len_trim(HeaderUnits)
do i = 1, Mlist
   if (iCSV(i) == 1) then
      HeaderUnits = HeaderUnits(1:il) // "," // trim(Units(i)); il = len_trim(HeaderUnits)
   end if
end do

end subroutine make_headerunits_tz

end subroutine csv_out_tz

end module SWAP_csv_output_tz

function STRinARSTR(string, arraystring, N)
implicit none
! global
integer,          intent(in)               :: N
character(len=*), intent(in)               :: string
character(len=*), dimension(N), intent(in) :: arraystring
logical                                    :: STRinARSTR
! local
integer                                    :: i, il, ip

STRinARSTR = .false.
do i = 1, N
   il = len_trim(trim(string))
   ip = index(string, "[")
   if (ip > 0) il = ip
   if (trim(string(1:il)) == trim(arraystring(i))) then
      STRinARSTR = .true.
      exit
   end if
end do
end function STRinARSTR

subroutine writeheader(iun)
implicit none
integer, intent(in) :: iun
write(iun,'(A)') '* Project:       '
write(iun,'(A)') '* File content:  '
write(iun,'(A)') '* File name:     '
write(iun,'(A)') '* Model version: '
write(iun,'(A)') '* Generated at:  '
end subroutine writeheader

SUBROUTINE my_piksrt(n,arr,brr)
! adapted after Numerical Recipes routine piksr2 (inputs now 1 integer and 1 char arrays)
   implicit none
   ! global
   INTEGER,                        intent(in)      :: n
   integer,          dimension(n), intent(inout)   :: arr
   character(len=*), dimension(n), intent(inout)   :: brr
   ! local
   INTEGER                                         :: i, j, a
   character(len=80)                               :: b

   do j = 2, n
      a = arr(j)
      b = brr(j)
      do i = j-1, 1, -1
         if (arr(i) <= a) goto 10
         arr(i+1) = arr(i)
         brr(i+1) = brr(i)
      end do
      i = 0
10    arr(i+1) = a
      brr(i+1) = b
   end do
   return
END
