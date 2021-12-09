program main

use mo_fio_utils,      ONLY: getunit_filerw,ReadNameList,invars_fname
use mo_fio_utils,      only: onetond,kMax,nvars,nlin,ncols,deltat
use mo_fio_utils,      only: TRC,LV,EXTLON,EXTLAT
use mo_fio_utils,      only: jact,jnum,jnuc,jqsmall
use Micro_HugMorr,     only: Init_Micro_HugMorr,MP_MORR_TWO_MOMENT

implicit none

integer, parameter :: PBL=1
integer, parameter :: r8  = SELECTED_REAL_KIND(P=13,R=300)
integer, parameter :: sizeofreal=4
!
logical :: lex
integer :: runit,io
integer :: nx,ny
integer :: i,nt,k
integer :: nsteps
real(r8), allocatable :: linp(:)
real(r8), allocatable :: tinp(:,:)
!
character(len=48) :: ifname,ofname
!
real(r8), allocatable, dimension(:,:) :: pc,tc,dzc,wc,wvarc
real(r8), allocatable, dimension(:,:) :: qv,npccn
real(r8), allocatable, dimension(:,:) :: qc,qi,qs,qr,qg
real(r8), allocatable, dimension(:,:) :: nc,ni,ns,nr,ng
real(r8), allocatable, dimension(:,:) :: qin1,qin2,qin3,qin4
real(r8), allocatable, dimension(:,:) :: qcw1,qcw2,qcw3,qcw4
real(r8), allocatable, dimension(:,:) :: cldo,cldn
!
real(r8), allocatable, dimension(:,:) :: qout1,qout2,qout3,qout4
real(r8), allocatable, dimension(:,:) :: qcwout1,qcwout2,qcwout3,qcwout4
!
real(r8), allocatable, dimension(:,:) :: dTcdt
real(r8), allocatable, dimension(:,:) :: dqvdt
real(r8), allocatable, dimension(:,:) :: dqcdt
real(r8), allocatable, dimension(:,:) :: dqrdt
real(r8), allocatable, dimension(:,:) :: dqidt
real(r8), allocatable, dimension(:,:) :: dqsdt
real(r8), allocatable, dimension(:,:) :: dqgdt
real(r8), allocatable, dimension(:,:) :: dnidt
real(r8), allocatable, dimension(:,:) :: dnsdt
real(r8), allocatable, dimension(:,:) :: dnrdt
real(r8), allocatable, dimension(:,:) :: dNGdt
real(r8), allocatable, dimension(:,:) :: dNCdt
real(r8), allocatable, dimension(:,:,:) :: EFFCS
real(r8), allocatable, dimension(:,:,:) :: EFFIS
!
real(r8), allocatable, dimension(:,:) :: qrcuten
real(r8), allocatable, dimension(:,:) :: qscuten
real(r8), allocatable, dimension(:,:) :: qicuten
real(r8), allocatable, dimension(:,:,:) :: refl_10cm,qndrop,tke
real(r8), allocatable, dimension(:,:) :: P,dz,kzh,w
real(r8), allocatable, dimension(:,:,:) :: Tc_mic,qv_mic
real(r8), allocatable, dimension(:,:,:) :: qc_mic,qr_mic,qi_mic,qs_mic,qg_mic
real(r8), allocatable, dimension(:,:,:) :: NC_mic,nr_mic,ni_mic,ns_mic,NG_mic
real(r8), allocatable, dimension(:,:) :: SR,RAINNC,RAINNCV,SNOW
!
real(r8) :: LSRAIN(ncols),LSSNOW(ncols)
real(r8) :: mu(ncols)
LOGICAL :: diagflag
INTEGER :: do_radar_ref  ! GT added for reflectivity calcs
LOGICAL :: F_QNDROP      ! wrf-chem
real(r8) :: dt_in
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ReadNameList

call Init_Micro_HugMorr(jact,jnum,jnuc,jqsmall)

nx=nvars
ny=nlin
nsteps=ny/kMax
dt_in=deltat
!
inquire(file=trim(invars_fname),exist=lex)
if (lex) then
  runit=getunit_filerw(30,100)
else
  write(*,*)'File not found. Exiting!'
  stop
end if

allocate(linp(nx))
allocate(tinp(ny,nx))
open(unit=runit,file=invars_fname,status='old',access='sequential',recl=ny*nx*sizeofreal,form='formatted')
do i=1,ny
  read(runit,iostat=io,fmt='(1e10.2)')linp
  tinp(i,1)=int(linp(1))
  tinp(i,2:nx)=linp(2:nx)
end do
close(runit)
!
! create file for output
write(ifname,'(a7,a4,a1,a4,a1,a9,a1,a6,a4)')'invars_',TRC,'_',LV,'_',EXTLAT,'_',EXTLON,'.dat'
write(ofname,'(a12,i1,a1,a4,a1,a4,a1,a9,a1,a6,a4)')'outvars_iact',jact,'_',TRC,'_',LV,'_',EXTLAT,'_',EXTLON,'.dat'
!
allocate(dTcdt(ncols,kMax))
allocate(dqvdt(ncols,kMax))
allocate(dqcdt(ncols,kMax))
allocate(dqrdt(ncols,kMax))
allocate(dqidt(ncols,kMax))
allocate(dqsdt(ncols,kMax))
allocate(dqgdt(ncols,kMax))
allocate(dnidt(ncols,kMax))
allocate(dnsdt(ncols,kMax))
allocate(dnrdt(ncols,kMax))
allocate(dNGdt(ncols,kMax))
allocate(dNCdt(ncols,kMax))
allocate(EFFCS(ncols,kMax,nsteps))
allocate(EFFIS(ncols,kMax,nsteps))
allocate(qrcuten(ncols,kMax))
allocate(qscuten(ncols,kMax))
allocate(qicuten(ncols,kMax))
allocate(refl_10cm(ncols,kMax,nsteps))
allocate(qndrop(ncols,kMax,nsteps))
allocate(P(ncols,kMax))
allocate(dz(ncols,kMax))
allocate(tke(ncols,kMax,nsteps))
allocate(kzh(ncols,kMax))
allocate(w(ncols,kMax))
allocate(Tc_mic(ncols,kMax,nsteps))
allocate(qv_mic(ncols,kMax,nsteps))
allocate(qc_mic(ncols,kMax,nsteps))
allocate(qr_mic(ncols,kMax,nsteps))
allocate(qi_mic(ncols,kMax,nsteps))
allocate(qs_mic(ncols,kMax,nsteps))
allocate(qg_mic(ncols,kMax,nsteps))
allocate(ni_mic(ncols,kMax,nsteps))
allocate(ns_mic(ncols,kMax,nsteps))
allocate(nr_mic(ncols,kMax,nsteps))
allocate(NG_mic(ncols,kMax,nsteps))
allocate(NC_mic(ncols,kMax,nsteps))
allocate(RAINNC(ncols,nsteps))
allocate(RAINNCV(ncols,nsteps))
allocate(SNOW(ncols,nsteps))
allocate(SR(ncols,nsteps))
!
P=0.0
dz=0.0
tke=1.0
kzh=0.0
w=0.0
dTcdt=0.0
dqvdt=0.0
dqcdt=0.0
dqrdt=0.0
dqidt=0.0
dqsdt=0.0
dqgdt=0.0
dnidt=0.0
dnsdt=0.0
dnrdt=0.0
dNGdt=0.0
dNCdt=0.0
EFFCS=0.0
EFFIS=0.0
qrcuten=0.0
qscuten=0.0
qicuten=0.0
RAINNC =0.0_r8
refl_10cm=0.0_r8
RAINNCV=0.0_r8
SNOW   =0.0_r8
mu     =1.0_r8
F_QNDROP=.FALSE.
qndrop=0.0_r8
diagflag=.TRUE.
!
allocate(pc(kMax,nsteps))
allocate(tc(kMax,nsteps))
allocate(dzc(kMax,nsteps))
allocate(wc(kMax,nsteps))
allocate(wvarc(kMax,nsteps))
allocate(qv(kMax,nsteps))
allocate(npccn(kMax,nsteps))
allocate(qc(kMax,nsteps))
allocate(qi(kMax,nsteps))
allocate(qs(kMax,nsteps))
allocate(qr(kMax,nsteps))
allocate(qg(kMax,nsteps))
allocate(nc(kMax,nsteps))
allocate(ni(kMax,nsteps))
allocate(ns(kMax,nsteps))
allocate(nr(kMax,nsteps))
allocate(ng(kMax,nsteps))
allocate(qin1(kMax,nsteps))
allocate(qin2(kMax,nsteps))
allocate(qin3(kMax,nsteps))
allocate(qin4(kMax,nsteps))
allocate(qcw1(kMax,nsteps))
allocate(qcw2(kMax,nsteps))
allocate(qcw3(kMax,nsteps))
allocate(qcw4(kMax,nsteps))
allocate(qout1(kMax,nsteps))
allocate(qout2(kMax,nsteps))
allocate(qout3(kMax,nsteps))
allocate(qout4(kMax,nsteps))
allocate(qcwout1(kMax,nsteps))
allocate(qcwout2(kMax,nsteps))
allocate(qcwout3(kMax,nsteps))
allocate(qcwout4(kMax,nsteps))
allocate(cldo(kMax,nsteps))
allocate(cldn(kMax,nsteps))

pc=onetond(tinp(:,2),ny,kMax,nsteps)
tc=onetond(tinp(:,3),ny,kMax,nsteps)
dzc=onetond(tinp(:,4),ny,kMax,nsteps)
wc=onetond(tinp(:,5),ny,kMax,nsteps)
wvarc=onetond(tinp(:,6),ny,kMax,nsteps)
qv=onetond(tinp(:,7),ny,kMax,nsteps)
npccn=onetond(tinp(:,18),ny,kMax,nsteps)
qc=onetond(tinp(:,8),ny,kMax,nsteps)
qi=onetond(tinp(:,9),ny,kMax,nsteps)
qs=onetond(tinp(:,10),ny,kMax,nsteps)
qr=onetond(tinp(:,11),ny,kMax,nsteps)
qg=onetond(tinp(:,12),ny,kMax,nsteps)
nc=onetond(tinp(:,13),ny,kMax,nsteps)
ni=onetond(tinp(:,14),ny,kMax,nsteps)
ns=onetond(tinp(:,15),ny,kMax,nsteps)
nr=onetond(tinp(:,16),ny,kMax,nsteps)
ng=onetond(tinp(:,17),ny,kMax,nsteps)
qin1=onetond(tinp(:,21),ny,kMax,nsteps)
qin2=onetond(tinp(:,25),ny,kMax,nsteps)
qin3=onetond(tinp(:,29),ny,kMax,nsteps)
qin4=onetond(tinp(:,33),ny,kMax,nsteps)
qcw1=onetond(tinp(:,22),ny,kMax,nsteps)
qcw2=onetond(tinp(:,26),ny,kMax,nsteps)
qcw3=onetond(tinp(:,30),ny,kMax,nsteps)
qcw4=onetond(tinp(:,34),ny,kMax,nsteps)
!
cldo=onetond(tinp(:,35),ny,kMax,nsteps)
cldn=onetond(tinp(:,36),ny,kMax,nsteps)
!
qout1=onetond(tinp(:,19),ny,kMax,nsteps)
qout2=onetond(tinp(:,23),ny,kMax,nsteps)
qout3=onetond(tinp(:,27),ny,kMax,nsteps)
qout4=onetond(tinp(:,31),ny,kMax,nsteps)
qcwout1=onetond(tinp(:,20),ny,kMax,nsteps)
qcwout2=onetond(tinp(:,24),ny,kMax,nsteps)
qcwout3=onetond(tinp(:,28),ny,kMax,nsteps)
qcwout4=onetond(tinp(:,32),ny,kMax,nsteps)
!
open(23,file=trim(ifname),access='sequential',recl=ny*nx*sizeofreal,form='formatted')
open(29,file=trim(ofname),access='sequential',recl=ny*nx*sizeofreal,form='formatted')
! --- time loop
do nt=1,nsteps
  !
  do k=1,kMax
    do i=1,ncols
      P(i,k) = pc(k,nt)
      dz(i,k)= dzc(k,nt)
      w(i,k) = wc(k,nt)
      kzh(i,k) = wvarc(k,nt)*30.0
      Tc_mic(i,k,nt) = tc(k,nt)
      qv_mic(i,k,nt) = qv(k,nt)
      qc_mic(i,k,nt) = qc(k,nt)
      qr_mic(i,k,nt) = qr(k,nt)
      qi_mic(i,k,nt) = qi(k,nt)
      qs_mic(i,k,nt) = qs(k,nt)
      qg_mic(i,k,nt) = qg(k,nt)
      ni_mic(i,k,nt) = ni(k,nt)
      ns_mic(i,k,nt) = ns(k,nt)
      nr_mic(i,k,nt) = nr(k,nt)
      NG_mic(i,k,nt) = ng(k,nt)
      NC_mic(i,k,nt) = nc(k,nt)
    end do
  end do
  
  CALL MP_MORR_TWO_MOMENT(   &
       nCols               , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                , &!INTEGER      , INTENT(IN   ) :: kMax 
       tc_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: Tc_mic (1:nCols, 1:kMax)
       QV_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: qv_mic (1:nCols, 1:kMax)
       QC_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: qc_mic (1:nCols, 1:kMax)
       QR_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: qr_mic (1:nCols, 1:kMax)
       QI_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: qi_mic (1:nCols, 1:kMax)
       QS_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: qs_mic (1:nCols, 1:kMax)
       QG_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: qg_mic (1:nCols, 1:kMax)
       NI_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: ni_mic (1:nCols, 1:kMax)
       NS_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: ns_mic (1:nCols, 1:kMax)
       NR_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: nr_mic (1:nCols, 1:kMax)
       NG_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: NG_mic (1:nCols, 1:kMax)   
       NC_mic      (1:ncols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: NC_mic (1:nCols, 1:kMax)  
       TKE         (1:nCols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: TKE (1:nCols, 1:kMax)  
       KZH         (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(INOUT) :: KZH (1:nCols, 1:kMax)  
       PBL                          , &!INTEGER(KIND=r8), INTENT(IN   ) :: PBL
       P           (1:nCols, 1:kMax), &!AIR PRESSURE (PA)
       DT_IN                        , &!REAL(KIND=r8), INTENT(IN   ) :: dt_in
       DZ          (1:nCols, 1:kMax), &!hm
       W           (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: w  (1:nCols, 1:kMax) !, tke, nctend, nnColsnd,kzh
       RAINNC      (1:nCols,nt)        , &!REAL(KIND=r8), INTENT(INOUT) :: RAINNC (1:nCols)
       RAINNCV     (1:nCols,nt)        , &!REAL(KIND=r8), INTENT(INOUT) :: RAINNCV(1:nCols)
       SNOW        (1:nCols,nt)        , &!REAL(KIND=r8), INTENT(INOUT) :: SNOW(1:nCols)
       SR          (1:nCols,nt)        , &!REAL(KIND=r8), INTENT(INOUT) :: SR     (1:nCols)
       refl_10cm   (1:nCols, 1:kMax,nt), &!REAL(KIND=r8), INTENT(INOUT) :: refl_10cm(1:nCols, 1:kMax)
       qrcuten     (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: qrcuten(1:nCols, 1:kMax)
       qscuten     (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: qscuten(1:nCols, 1:kMax)
       qicuten     (1:nCols, 1:kMax), &!REAL(KIND=r8), INTENT(IN   ) :: qicuten(1:nCols, 1:kMax)
       mu          (1:nCols)        , &!REAL(KIND=r8), INTENT(IN   ) :: mu     (1:nCols)  ! hm added
       diagflag                     , &!LOGICAL      , OPTIONAL, INTENT(IN) :: diagflag
       do_radar_ref                 , &!INTEGER      , OPTIONAL, INTENT(IN) :: do_radar_ref ! GT added for reflectivity calcs
       F_QNDROP                     , &!LOGICAL      , OPTIONAL, INTENT(IN) :: F_QNDROP  ! wrf-chem
       qndrop      (1:nCols, 1:kMax,nt), &!REAL(KIND=r8), OPTIONAL, INTENT(INOUT):: qndrop(1:nCols, 1:kMax) ! hm added, wrf-chem 
       EFFCS       (1:nCols, 1:kMax,nt), &
       EFFIS       (1:nCols, 1:kMax,nt)  &
       )

  ! saving the variables
  do k=1,kMax
    do i=1,ncols
      write(23,'(i4,1p,15e10.2)')k,pc(k,nt),tc(k,nt),wc(k,nt), &
        qv(k,nt),cldn(k,nt), &
        qc(k,nt),qi(k,nt),qs(k,nt),qr(k,nt),qg(k,nt),   &
        nc(k,nt),ni(k,nt),ns(k,nt),nr(k,nt),ng(k,nt)
      !
      write(29,'(i4,1p,17e10.2)')k,P(i,k),tc_mic(i,k,nt),W(i,k), &
        qv_mic(i,k,nt),cldn(k,nt),      &
        qc_mic(i,k,nt),qi_mic(i,k,nt),qs_mic(i,k,nt),qr_mic(i,k,nt),qg_mic(i,k,nt),  &
        nc_mic(i,k,nt),ni_mic(i,k,nt),ns_mic(i,k,nt),nr_mic(i,k,nt),ng_mic(i,k,nt),  &
        effcs(i,k,nt),effis(i,k,nt)
    end do
  end do

end do
close(23)
close(29)
!
! deallocate
deallocate(pc,tc,dzc,wc,wvarc)
deallocate(qv,npccn)
deallocate(qc,qi,qs,qr,qg)
deallocate(nc,ni,ns,nr,ng)
deallocate(qin1,qin2,qin3,qin4)
deallocate(qcw1,qcw2,qcw3,qcw4)
deallocate(qout1,qout2,qout3,qout4)
deallocate(qcwout1,qcwout2,qcwout3,qcwout4)
deallocate(cldo,cldn)
!
deallocate(dTcdt)
deallocate(dqvdt)
deallocate(dqcdt)
deallocate(dqrdt)
deallocate(dqidt)
deallocate(dqsdt)
deallocate(dqgdt)
deallocate(dnidt)
deallocate(dnsdt)
deallocate(dnrdt)
deallocate(dNGdt)
deallocate(dNCdt)
deallocate(EFFCS)
deallocate(EFFIS)
deallocate(qrcuten)
deallocate(qscuten)
deallocate(qicuten)
!
end program main
