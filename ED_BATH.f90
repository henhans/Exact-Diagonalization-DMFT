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
  public :: write_bath
  public :: set_bath
  public :: copy_bath

  public :: delta_and
  public :: delta_bath


  !Bath parameters (to be used in H)
  !The total bath size should be:
  !N=(Norb+1)[energies&hybridizations]*Nspin[# of spins]*Nbath[# of bath sites]
  !=========================================================
  real(8),allocatable,dimension(:,:),public   :: ebath
  real(8),allocatable,dimension(:,:,:),public :: vbath
  logical                                     :: bath_status=.false.

contains

  subroutine allocate_bath()
    allocate(ebath(Nspin,Nbath),vbath(Norb,Nspin,Nbath))
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
    if(N_ < (Norb+1)*Nspin*Nbath)&
         call error("CHECK_BATH_DIMENSION: error")
  end subroutine check_bath_dimension


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function get_bath_size() result(N)
    integer :: N
    N=(Norb+1)*Nspin*Nbath
  end function get_bath_size



  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed
    integer :: i,iorb,ispin,unit
    logical :: IOfile
    real(8) :: ran(Nbath)
    character(len=100) :: foo
    if(bath_status)call deallocate_bath
    call allocate_bath
    inquire(file=trim(Hfile),exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")bg_green('Reading bath/xmu from file')
       unit = free_unit()
       open(unit,file=trim(Hfile))
       read(unit,*)
       do i=1,Nbath
          read(unit,"(90(F13.9,1X))")(ebath(ispin,i),(vbath(iorb,ispin,i),iorb=1,Norb),ispin=1,Nspin)
       enddo
       close(unit)
    else
       write(LOGfile,"(A)")bg_red('Generating bath from scratch')
       call random_number(ran(:))
       do ispin=1,Nspin
          do i=1,Nbath
             ebath(ispin,i)=(2.d0*ran(i)-1.d0)*real(Nbath,8)/2.d0
             vbath(ispin,1:Norb,i)=1.d0/sqrt(real(Nbath,8))
          enddo
       enddo


    endif
  end subroutine init_bath_ed



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine write_bath(unit)
    integer :: i,unit,ispin,iorb
    if(.not.bath_status)call error("WRITE_BATH: bath not allocated")
    write(unit,"(90(A13,1X))")("# ek_s"//trim(adjustl(trim(txtfy(ispin)))),&
         ("Vk^orb"//trim(adjustl(trim(txtfy(iorb))))//"_"//&
         trim(adjustl(trim(txtfy(ispin)))),iorb=1,Norb),ispin=1,Nspin)
    do i=1,Nbath
       write(unit,"(90(F13.9,1X))")(ebath(ispin,i),&
            (vbath(iorb,ispin,i),iorb=1,Norb),ispin=1,Nspin)
    enddo
  end subroutine write_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath)
    real(8),dimension(:) :: bath
    integer              :: iorb,ispin,stride
    integer              :: ae,be
    integer              :: av(Norb),bv(Norb)
    if(.not.bath_status)call error("SET_BATH: bath not allocated")
    do ispin=1,Nspin
       stride=(ispin-1)*(Norb+1)*Nbath
       ae=stride + 1
       be=stride + Nbath
       ebath(ispin,:) = bath(ae:be)
       do iorb=1,Norb
          av(iorb)=stride + iorb*Nbath + 1
          bv(iorb)=stride +(iorb+1)*Nbath
          vbath(iorb,ispin,:) = bath(av(iorb):bv(iorb))
       enddo
    enddo
  end subroutine set_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function copy_bath() result(bath)
    real(8),dimension((Norb+1)*Nspin*Nbath) :: bath
    integer                                 :: iorb,ispin,stride
    integer                                 :: ae,be
    integer                                 :: av(Norb),bv(Norb)
    if(.not.bath_status)call error("COPY_BATH: bath not allocated")
    do ispin=1,Nspin
       stride=(ispin-1)*(Norb+1)*Nbath
       ae=stride + 1
       be=stride + Nbath
       bath(ae:be) = ebath(ispin,:)
       do iorb=1,Norb
          av(iorb)=stride + iorb*Nbath + 1
          bv(iorb)=stride +(iorb+1)*Nbath
          bath(av(iorb):bv(iorb)) = vbath(iorb,ispin,:)
       enddo
    enddo
    ! do ispin=1,Nspin
    !    a1=1+(ispin-1)*2*Nbath
    !    b1=Nbath+(ispin-1)*2*Nbath
    !    a2=1+Nbath+(ispin-1)*2*Nbath
    !    b2=2*Nbath+(ispin-1)*2*Ntot
    !    bath(a1:b1) = ebath(ispin,1:Nbath)
    !    bath(a2:b2) = vbath(ispin,1:Nbath)
    ! enddo
  end function copy_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  pure function delta_and(x,bath,orb1,orb2,spin) result(fg)
    real(8),dimension((Norb+1)*Nspin*Nbath),intent(in) :: bath
    complex(8),intent(in)                              :: x
    integer,intent(in)                                 :: orb1,orb2,spin
    complex(8)                                         :: fg
    integer                                            :: i,iorb,ispin,stride
    integer                                            :: ae,be
    integer                                            :: av(Norb),bv(Norb)
    real(8),dimension(Nspin,Nbath)                     :: epsk
    real(8),dimension(Norb,Nspin,Nbath)                :: vpsk
    do ispin=1,Nspin
       stride=(ispin-1)*(Norb+1)*Nbath
       ae=stride + 1
       be=stride + Nbath
       epsk(ispin,:) = bath(ae:be)
       do iorb=1,Norb
          av(iorb)=stride + iorb*Nbath + 1
          bv(iorb)=stride +(iorb+1)*Nbath
          vpsk(iorb,ispin,:) = bath(av(iorb):bv(iorb))
       enddo
    enddo
    fg=zero
    do i=1,Nbath
       fg=fg + vpsk(orb1,spin,i)*vpsk(orb2,spin,i)/(x-epsk(spin,i))
    enddo
  end function delta_and



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  function delta_bath(x,iorb,jorb,ispin) result(fg)
    complex(8),intent(in)                       :: x
    integer,intent(in)                          :: iorb,jorb,ispin
    complex(8)                                  :: fg
    integer                                     :: i
    if(.not.bath_status)call error("DELTA_BATH: bath not allocated")
    fg=zero
    do i=1,Nbath
       fg=fg + vbath(iorb,ispin,i)*vbath(jorb,ispin,i)/(x-ebath(ispin,i))
    enddo
  end function delta_bath


END MODULE ED_BATH
