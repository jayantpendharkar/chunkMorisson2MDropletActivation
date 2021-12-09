subroutine aerosolnact(ncols,kMax,ntot_amode,ncnst_tot)

implicit none


integer, parameter :: r8  = SELECTED_REAL_KIND(P=13,R=300)

integer :: ncols,kMax,ntot_amode,ncnst_tot

!
REAL(r8), PARAMETER :: r_d=287.04
REAL(r8), PARAMETER :: r_v=461.6
REAL(r8), PARAMETER ::  SVP1=0.6112
REAL(r8), PARAMETER ::  SVP2=17.67
REAL(r8), PARAMETER ::  SVP3=29.65
REAL(r8), PARAMETER ::  SVPT0=273.15
REAL(r8), PARAMETER ::  EP_1=r_v/r_d-1.
REAL(r8), PARAMETER ::  EP_2=r_d/r_v

real(r8) :: nspec_amode(ntot_amode)
real(r8) :: sigmag_amode(ntot_amode)
real(r8) :: dgnumlo_amode(ntot_amode)
real(r8) :: dgnumhi_amode(ntot_amode)
real(r8) :: alogsig(ntot_amode)
real(r8) :: exp45logsig(ntot_amode)
real(r8) :: f1(ntot_amode)
real(r8) :: f2(ntot_amode)
real(r8) :: voltonumblo_amode(ntot_amode)
real(r8) :: voltonumbhi_amode(ntot_amode)

real(r8) :: t0            ! reference temperature
real(r8) :: aten
real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
real(r8) :: alog2, alog3, alogaten
real(r8) :: third, twothird, sixth, zero
real(r8) :: sq2, sqpi

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8 !- default Jayant

! CCN diagnostic fields
integer, parameter :: psat=6    ! number of supersaturations to calc ccn concentration
real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn conc.
                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)
character(len=8) :: ccn_name(psat)= &
                    (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)

! indices in state and pbuf structures
integer :: numliq_idx = -1
integer :: kvh_idx    = -1










end subroutine aerosolnact
