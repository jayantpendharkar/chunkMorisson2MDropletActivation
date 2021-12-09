MODULE mo_fio_utils

IMPLICIT NONE

PRIVATE

integer, public, parameter :: ncols=1
!
integer, public, parameter :: nvars=36

integer, parameter :: r8  = SELECTED_REAL_KIND(P=13,R=300)

integer,  public :: trunc
integer,  public :: kMax
real(r8), public :: latc
real(r8), public :: lonc
integer,  public :: nlin
real(r8), public :: deltat
integer,  public :: jact,jnum,jnuc
real(r8), public :: jqsmall
!
character(len=6), public :: TRC
character(len=4), public :: LV
character(len=10),public :: EXTLAT
character(len=6), public :: EXTLON
character(len=24),public :: invars_fname

public :: ReadNameList
public :: getunit_filerw
public :: onetond

CONTAINS

subroutine ReadNameList

  character(len=64) :: fnamelist 
  integer :: nlunit,ierr
  logical :: lex

  integer :: io
  real(r8)    :: tvar(nvars)
  
  NAMELIST /MODEL_CONF/trunc,kMax,latc,lonc,deltat,jact,jnum,jnuc,jqsmall
  
  fnamelist = 'param.inp'
  
  inquire(file=trim(fnamelist),exist=lex)
  if (lex) then
    nlunit = getunit_filerw(30,100)
  else
    write(*,*)'Namelist file not found, exiting!'
    stop
  end if
  
  open(unit=nlunit,file=fnamelist,action="read",status="old",iostat=ierr)
  read(nlunit,MODEL_CONF,iostat=ierr)
  close(nlunit)
  
  write(TRC,'(a1,i3.3)')'T',trunc
  write(LV,'(a1,i3.3)')'L',kMax
  write(EXTLON,'(a3,i3)')'lon',int(lonc)
  if (latc .lt. 0.) then
    write(EXTLAT,'(a3,f5.2,a1)')'lat',abs(latc),'S'
  else
    write(EXTLAT,'(a3,f5.2,a1)')'lat',abs(latc),'N'
  end if
  
  write(invars_fname,'(a8,a9,a1,a6)')'metvars_',EXTLAT,'_',EXTLON

  nlin=0
  open(unit=11,file=invars_fname,status='old',access='sequential',recl=4,form='formatted')
  do 
    read(11,iostat=io,fmt='(1e10.2)')tvar
    if (io .eq. 0) nlin=nlin+1
    if (io .lt. 0) exit
    end do
  close (11)

end subroutine

FUNCTION getunit_filerw(i1,i2) RESULT(unit)

  integer :: i1,i2,unit,i
  logical :: i_opened

  DO i=i1,i2
    INQUIRE(i,opened=i_opened)
    IF (.NOT.i_opened) THEN
      unit=i
      EXIT
    END IF
  END DO

END FUNCTION getunit_filerw

function onetond(arrin,nlin,x1,x2) result(arrout)

  integer :: nlin,x1,x2
  real(r8)    :: arrin(nlin)
  real(r8)    :: arrout(x1,x2)

  arrout = reshape(arrin,(/x1,x2/))

end function onetond

END MODULE mo_fio_utils
