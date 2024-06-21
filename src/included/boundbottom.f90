! File VersionID:
!   $Id: boundbottom.f90 362 2018-01-08 13:08:33Z kroes006 $
! ----------------------------------------------------------------------
!> \brief Determines the soil profile bottom boundary conditions.
!!
!! This subroutine calculates the boundary conditions at the bottom of a soil profile.
!! It handles various cases such as given groundwater level, regional bottom flux, 
!! seepage or infiltration from/to deep groundwater, zero flux, free drainage, and 
!! lysimeter with free drainage.
!!
!! \date August 2004 / Sept 2005
! ----------------------------------------------------------------------
      subroutine BoundBottom
      use variables
      implicit none

! --- local variables
      integer node, nodnumgwl 

      real(8) cvalprof,gwlmean, thetabot, twopi, freq
      real(8) satnodgwl
      real(8) watcon,hconduc,afgen
      character(len=300) messag

! --- Constants
      twopi = 8.0d0*datan(1.0d0)
      freq  = twopi/365.0d0
! ----------------------------------------------------------------------
! --- interpolation between daily values of given groundwaterlevel
      if (swbotb.eq.1) then
          gwlinp = afgen(gwltab,mabbc*2,t1900+dt)
      endif

! --- regional bottom flux is given
      if (abs(swbotb) .eq. 2) then
 
! Comment PietG (8-1-08):
! ---   if the moisture content in the soil profile is depleted by 
!       a combination of inconsistent boundary conditions, and 
!       the pressure head at the bottom tends to very low values,
!       then the choice for swbotb=2 is not appropriate. 

        if(h(numnod).lt.-1.0d+7) then ! oven dry conditions at bottom
           if(swbotb.eq.2) then
              write(messag,'(a)')'Oven dry conditions in lowest'
              write(messag,'(a)')'compartment therefore switched to'
              write(messag,'(a)')'free drainage at date '
              write(messag,*)
              write(messag,'(a11)') date
              call warn ('BoundBottom',messag,logf,0)
           end if
           swbotb = -2
        else
           swbotb = 2
        end if

        if(swbotb .eq. 2)then
           if (sw2 .eq. 1) then
! ---     sine function is used
             qbot = sinave + sinamp * dcos( freq * (t-sinmax))        
           else
! ---     table is used
             qbot = afgen (qbotab,mabbc*2,t1900+dt)
           endif
        end if

! ---   free drainage assumed in case of h(numnod) < -1.0E7
        if (swbotb .eq. -2) qbot = -1.0d0 * kmean(numnod+1)

      endif

! --- seepage or infiltration from/to deep groundwater
      if (swbotb .eq. 3) then
        gwlmean = hdrain + shape*(gwl - hdrain)
! ---   determine hydraulic head of deep aquifer
        if (sw3 .eq. 1) then
          deepgw = aqave + aqamp * dcos( twopi/aqper*(t - aqtmax))
        else
          deepgw = afgen (haqtab,mabbc*2,t1900+dt)
        endif

! ---   determine C-value (vertical resistance) in saturated part of modelled profile 
        if (SwBotb3ResVert.eq.0) then
!     -   find number node with groundwater level
          node = numnod
          do while (gwlmean.gt.ztopcp(node) .and. node.gt.1)
            node = node - 1
          enddo
          nodnumgwl = node
          satnodgwl = gwlmean - zbotcp(nodnumgwl)
          cvalprof = satnodgwl/cofgen(3,nodnumgwl)
          do node = nodnumgwl+1, numnod
            cvalprof = cvalprof + dz(node)/cofgen(3,node)
          enddo
        elseif (SwBotb3ResVert.eq.1) then
          cvalprof = 0.0d0
        endif
!
        qbot   = (deepgw - gwlmean)/(rimlay+cvalprof)

! ---   extra groundwater flux might be added
        if (sw4 .eq. 1) qbot = qbot + afgen (qbotab,mabbc*2,t1900+dt)
      endif

! --- flux calculated as function of h
      if (swbotb .eq. 4 ) then
        if (swqhbot.eq.1) then
          qbot = cofqha * dexp(cofqhb * dabs(gwl))
          if (swcofqhc.eq.1) qbot = qbot + cofqhc
        else if (swqhbot.eq.2) then
          qbot = afgen(qbotab,mabbc*2,dabs(gwl))
        endif
      endif

! --- interpolation between daily values of given pressurehead
      if (swbotb .eq. 5) then
        hbot           = afgen (hbotab,mabbc*2,t1900+dt)
        thetabot = watcon(numnod,hbot)
                   
        kmean(numnod+1)= hconduc(numnod,hbot,thetabot,rfcp(numnod))
        if(flMacroPore) then
          kmean(numnod+1) = FrArMtrx(numnod) * kmean(numnod+1)
        endif
      endif

! --- zero flux at the bottom
      if (swbotb .eq. 6) qbot = 0.0d0

! --- free drainage
      if (swbotb .eq. 7) qbot = -1.0d0 * kmean(numnod+1)

! --- lysimeter with free drainage
      if (swbotb .eq. 8) qbot = 0.0d0

      qbot_nonfrozen = qbot

      return
      end
