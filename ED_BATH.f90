!########################################################################
!PROGRAM  : ED_BATH
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_BATH
  USE ED_VARS_GLOBAL
  implicit none

  private

  public :: allocate_bath
  public :: deallocate_bath
  public :: check_bath_dimension
  public :: init_bath_ed
  public :: get_bath_size
  public :: dump_bath
  public :: write_bath
  public :: set_bath
  public :: copy_bath

  public :: delta_and
  public :: delta_bath


  !Bath parameters (to be used in H)
  !=========================================================
  real(8),allocatable,dimension(:,:),public :: ebath
  real(8),allocatable,dimension(:,:),public :: vbath
  logical :: bath_status=.false.

contains

  subroutine allocate_bath()
    allocate(ebath(Nspin,Nbath),vbath(Nspin,Nbath))
    bath_status=.true.
  end subroutine allocate_bath


  subroutine deallocate_bath()
    deallocate(ebath,vbath)
    bath_status=.false.
  end subroutine deallocate_bath


  subroutine check_bath_dimension(bath)
    real(8),dimension(:) :: bath
    integer              :: N_
    N_=size(bath)
    if(N_< 2*Nspin*Nbath)&
         call error("CHECK_BATH_DIMENSION: error")
  end subroutine check_bath_dimension


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function get_bath_size() result(N)
    integer :: N
    N=2*Nspin*Nbath
  end function get_bath_size



  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed
    integer :: i,ispin
    logical :: IOfile
    inquire(file=trim(Hfile),exist=IOfile)
    if(.NOT.IOfile)then
       write(LOGfile,"(A)")bg_red('Generating bath from scratch')
       call guess_bath_params
    else
       write(LOGfile,"(A)")bg_green('Reading bath/xmu from file')
       open(51,file=trim(Hfile))
       do i=1,Nbath
          read(51,"(6(F13.9,1X))")(ebath(ispin,i),vbath(ispin,i),ispin=1,Nspin)
       enddo
       close(51)
    endif
  end subroutine init_bath_ed




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine dump_bath(bath_file)
    character(len=*) :: bath_file
    integer :: i,ispin
    if(.not.bath_status)call error("DUMP_BATH: bath not allocated")
    open(51,file=trim(bath_file))
    call write_bath(unit=51)
    close(51)
  end subroutine dump_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine write_bath(unit)
    integer :: i,unit,ispin 
    if(.not.bath_status)call error("DUMP_BATH: bath not allocated")
    do i=1,Nbath
       write(unit,"(6(F13.9,1X))")(ebath(ispin,i),vbath(ispin,i),ispin=1,Nspin)
    enddo
  end subroutine write_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath)
    real(8),dimension(:)   :: bath
    integer                :: ispin,a1,a2,b1,b2
    do ispin=1,Nspin
       a1=1+(ispin-1)*2*Nbath
       b1=Nbath+(ispin-1)*2*Nbath
       a2=1+Nbath+(ispin-1)*2*Nbath
       b2=2*Nbath+(ispin-1)*2*Ntot
       ebath(ispin,:) = bath(a1:b1)  
       vbath(ispin,:) = bath(a2:b2)
    enddo
  end subroutine set_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function copy_bath() result(bath)
    real(8),dimension(2*Nspin*Nbath) :: bath
    integer                          :: ispin,a1,a2,b1,b2
    if(.not.bath_status)call error("DUMP_BATH: bath not allocated")
    do ispin=1,Nspin
       a1=1+(ispin-1)*2*Nbath
       b1=Nbath+(ispin-1)*2*Nbath
       a2=1+Nbath+(ispin-1)*2*Nbath
       b2=2*Nbath+(ispin-1)*2*Ntot
       bath(a1:b1) = ebath(ispin,1:Nbath)
       bath(a2:b2) = vbath(ispin,1:Nbath)
    enddo
  end function copy_bath




  !+------------------------------------------------------------------+
  !PURPOSE  : Build the parameters for the Hamiltonian
  !+------------------------------------------------------------------+
  subroutine guess_bath_params
    integer :: i,ispin,n2
    real(8) :: ran(Nbath)
    if(.not.bath_status)call error("DUMP_BATH: bath not allocated")
    !n2=Nbath/2;if(n2==0)n2=1
    call random_number(ran(:))
    do ispin=1,Nspin
       do i=1,Nbath
          !i=0,Nbath-1
          !ebath(ispin,i+1)=2.d0*dfloat(i-1-n2)/dfloat(n2)
          !vbath(ispin,i+1)=dsqrt(1.d0/dfloat(Nbath))
          ebath(ispin,i)=(2.d0*ran(i)-1.d0)*real(Nbath,8)/2.d0
          vbath(ispin,i)=1.d0/sqrt(real(Nbath,8))
       enddo
    enddo
  end subroutine guess_bath_params



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  pure function delta_and(x,bath,ichan) result(fg)
    real(8),dimension(2*Nspin*Nbath),intent(in) :: bath
    complex(8),intent(in)                       :: x
    integer,intent(in)                          :: ichan
    complex(8)                                  :: fg
    integer                                     :: i,a1,a2,b1,b2
    real(8),dimension(Nbath)                    :: ee,vv
    a1=1+(ichan-1)*2*Nbath
    b1=Nbath+(ichan-1)*2*Nbath
    a2=1+Nbath+(ichan-1)*2*Nbath
    b2=2*Nbath+(ichan-1)*2*Ntot
    ee(1:Nbath) = bath(a1:b1)  
    vv(1:Nbath) = bath(a2:b2)
    fg=zero
    do i=1,Nbath
       fg=fg+vv(i)**2/(x-ee(i))
    enddo
  end function delta_and



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function delta_bath(x,ichan) result(fg)
    complex(8),intent(in)                       :: x
    integer,intent(in)                          :: ichan
    complex(8)                                  :: fg
    integer                                     :: i
    if(.not.bath_status)call error("DUMP_BATH: bath not allocated")
    fg=zero
    do i=1,Nbath
       fg=fg+vbath(ichan,i)**2/(x-ebath(ichan,i))
    enddo
  end function delta_bath


END MODULE ED_BATH
