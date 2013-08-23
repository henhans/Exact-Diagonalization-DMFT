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
  !The dimensions of the bath components are:
  !ebath =Nspin[# of spins] * Nbath[# of bath sites per orbital] * Norb[# of orbitals]
  !vbath =Nspin[# of spins] * Nbath[# of bath sites per orbital] * Norb[# of orbitals]
  !N=(Norb+1)[energies&hybridizations]*Nspin[# of spins]*Nbath[# of bath sites]
  !=========================================================
  real(8),allocatable,dimension(:,:,:),public :: ebath
  real(8),allocatable,dimension(:,:,:),public :: vbath
  logical                                     :: bath_status=.false.

contains

  subroutine allocate_bath()
    allocate(ebath(Nspin,Norb,Nbath),vbath(Nspin,Norb,Nbath))
    bath_status=.true.
  end subroutine allocate_bath


  subroutine deallocate_bath()
    deallocate(ebath,vbath)
    bath_status=.false.
  end subroutine deallocate_bath


  subroutine check_bath_dimension(bath)
    real(8),dimension(:) :: bath
    integer              :: N_,Ntrue
    N_=size(bath)
    Ntrue = 2*Nspin*Norb*Nbath
    if(N_ /= Ntrue)&
         call error("CHECK_BATH_DIMENSION: wrong dimensions!")
  end subroutine check_bath_dimension


  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_size() result(N)
    integer :: N
    N=2*Nspin*Norb*Nbath
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
          read(unit,"(90(F22.15,1X))")((ebath(ispin,iorb,i),vbath(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
       enddo
       close(unit)
    else
       write(LOGfile,"(A)")bg_red('Generating bath from scratch')
       call random_number(ran(:))
       do ispin=1,Nspin
          do i=1,Nbath
             ebath(ispin,1:Norb,i)=(2.d0*ran(i)-1.d0)*real(Nbath,8)/2.d0
             vbath(ispin,1:Norb,i)=1.d0/sqrt(real(Nbath,8))
          enddo
       enddo
    endif
  end subroutine init_bath_ed



  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_bath(unit)
    integer :: i,unit,ispin,iorb
    if(.not.bath_status)call error("WRITE_BATH: bath not allocated")
    write(unit,"(90(A22,1X))")&
         (("# Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
         "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    do i=1,Nbath
       write(unit,"(90(F22.15,1X))")((ebath(ispin,iorb,i),vbath(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
    enddo
  end subroutine write_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath)
    real(8),dimension(:) :: bath
    integer              :: iorb,ispin,stride_spin,stride_orb
    integer              :: ae,be,av,bv
    if(.not.bath_status)call error("SET_BATH: bath not allocated")
    call check_bath_dimension(bath)
    do ispin=1,Nspin
       stride_spin=(ispin-1)*Norb*Nbath
       do iorb=1,Norb
          stride_orb=(iorb-1)*2*Nbath
          ae=stride_spin + stride_orb + 1
          be=stride_spin + stride_orb + Nbath
          av=stride_spin + stride_orb + Nbath + 1
          bv=stride_spin + stride_orb + Nbath + Nbath
          ebath(ispin,iorb,1:Nbath) = bath(ae:be)
          vbath(ispin,iorb,1:Nbath) = bath(av:bv)
       enddo
    enddo
  end subroutine set_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  function copy_bath() result(bath)
    real(8),dimension(2*Nspin*Norb*Nbath) :: bath
    integer                               :: iorb,ispin,stride_spin,stride_orb
    integer                               :: ae,be,av,bv
    if(.not.bath_status)call error("COPY_BATH: bath not allocated")
    do ispin=1,Nspin
       stride_spin=(ispin-1)*Norb*Nbath
       do iorb=1,Norb
          stride_orb=(iorb-1)*2*Nbath
          ae=stride_spin + stride_orb + 1
          be=stride_spin + stride_orb + Nbath
          av=stride_spin + stride_orb + Nbath + 1
          bv=stride_spin + stride_orb + Nbath + Nbath
          bath(ae:be) = ebath(ispin,iorb,1:Nbath)
          bath(av:bv) = vbath(ispin,iorb,1:Nbath)
       enddo
    enddo
  end function copy_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : given the bath array, compute the hybridization function
  ! for a given spin and orbital indices ispin,iorb at a given point x
  !+-------------------------------------------------------------------+
  pure function delta_and(ispin,iorb,x,bath) result(fg)
    integer,intent(in)                               :: ispin,iorb
    complex(8),intent(in)                            :: x
    real(8),dimension(2*Nspin*Norb*Nbath),intent(in) :: bath
    complex(8)                                       :: fg
    integer                                          :: i,stride_spin,stride_orb
    integer                                          :: ae,be,av,bv
    real(8),dimension(Nbath)                         :: epsk
    real(8),dimension(Nbath)                         :: vpsk
    stride_spin=(ispin-1)*Norb*Nbath
    stride_orb=(iorb-1)*2*Nbath
    ae=stride_spin + stride_orb + 1
    be=stride_spin + stride_orb + Nbath
    av=stride_spin + stride_orb + Nbath + 1
    bv=stride_spin + stride_orb + Nbath + Nbath
    epsk(1:Nbath) = bath(ae:be)
    vpsk(1:Nbath) = bath(av:bv)
    fg=zero
    do i=1,Nbath
       fg=fg + vpsk(i)**2/(x-epsk(i))
    enddo
  end function delta_and




  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function for a given spin and 
  ! orbital indices ispin and iorb at point x, from determined bath 
  ! components ebath,vbath
  !+-------------------------------------------------------------------+
  function delta_bath(ispin,iorb,x) result(fg)
    complex(8),intent(in)                       :: x
    integer,intent(in)                          :: iorb,ispin
    complex(8)                                  :: fg
    integer                                     :: i
    if(.not.bath_status)call error("DELTA_BATH: bath not allocated")
    fg=zero
    do i=1,Nbath
       fg=fg + vbath(ispin,iorb,i)**2/(x-ebath(ispin,iorb,i))
    enddo
  end function delta_bath


END MODULE ED_BATH
