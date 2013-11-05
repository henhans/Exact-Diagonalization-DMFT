!########################################################################
!PROGRAM  : ED_BATH
!AUTHORS  : Adriano Amaricci
! The dimensions of the bath components are:
! e = Nspin[# of spins] * Norb[# of orbitals] * Nbath[# of bath sites per orbital]
! v = Nspin[# of spins] * Norb[# of orbitals] * Nbath[# of bath sites per orbital]
! N = Nspin*(2*Norb)*Nbath
!
! e = Nspin[# of spins] * 1 * Nbath[# of bath sites]
! v = Nspin[# of spins] * Norb[# of orbitals] * Nbath[# of bath sites]
! N = Nspin*(Norb+1)*Nbath
!########################################################################
MODULE ED_BATH
  USE ED_VARS_GLOBAL
  implicit none
  private

  !type(effective_bath) :: dmft_bath

  interface delta_and
     module procedure delta_and_irr,delta_and_red
  end interface delta_and

  interface delta_bath
     module procedure delta_bath_irr,delta_bath_red
  end interface delta_bath

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
  !
  public :: dmft_bath



contains

  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_bath(dmft_bath)
    type(effective_bath) :: dmft_bath
    if(dmft_bath%status)call deallocate_bath(dmft_bath)
    select case(bath_type)
    case('hybrid')
       allocate(dmft_bath%e(Nspin,1,Nbath),dmft_bath%v(Nspin,Norb,Nbath))
    case default
       allocate(dmft_bath%e(Nspin,Norb,Nbath),dmft_bath%v(Nspin,Norb,Nbath))
    end select
    dmft_bath%status=.true.
  end subroutine allocate_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_bath(dmft_bath)
    type(effective_bath) :: dmft_bath
    if(allocated(dmft_bath%e))deallocate(dmft_bath%e)
    if(allocated(dmft_bath%v))deallocate(dmft_bath%v)
    dmft_bath%status=.false.
  end subroutine deallocate_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  subroutine check_bath_dimension(bath_)
    real(8),dimension(:) :: bath_
    integer              :: N_,Ntrue
    N_=size(bath_)
    select case(bath_type)
    case('hybrid')
       Ntrue=Nspin*(Norb+1)*Nbath
    case default
       Ntrue = Nspin*(2*Norb)*Nbath
    end select
    if(N_ /= Ntrue)stop "CHECK_BATH_DIMENSION: wrong dimensions!"
  end subroutine check_bath_dimension


  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_size() result(N)
    integer :: N
    select case(bath_type)
    case('hybrid')
       N = Nspin*(Norb+1)*Nbath
    case default
       N = Nspin*(2*Norb)*Nbath
    end select
  end function get_bath_size



  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed(dmft_bath)
    type(effective_bath) :: dmft_bath
    integer              :: i,iorb,ispin,unit,n2
    logical              :: IOfile
    real(8)              :: ran(Nbath)
    if(.not.dmft_bath%status)stop "init_bath: bath not allocated"
    if(mpiID==0)then
       inquire(file=trim(Hfile),exist=IOfile)
       if(IOfile)then
          write(LOGfile,"(A)")'Reading bath from file'
          write(LOGfile,"(A)")'- - - - - - - - - - - -'
          unit = free_unit()
          open(unit,file=trim(Hfile))
          read(unit,*)
          select case(bath_type)
          case default
             do i=1,Nbath
                read(unit,"(90(F22.15,1X))")((dmft_bath%e(ispin,iorb,i),dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          case ('hybrid')
             do i=1,Nbath
                read(unit,"(90(F13.9,1X))")( dmft_bath%e(ispin,1,i),    (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo

          end select
          close(unit)
       else
          write(LOGfile,"(A)")"Generating bath from scratch"
          write(LOGfile,"(A)")'- - - - - - - - - - - - - - -'
          !call random_number(ran(:))
          N2=Nbath/2
          do i=1,Nbath
             dmft_bath%e(:,:,i)=2.d0*real(i-1-n2,8)/real(n2,8)
             dmft_bath%v(:,:,i)=1.d0/sqrt(real(Nbath,8))
          enddo
       endif
    endif
#ifdef _MPI
    call MPI_BCAST(dmft_bath%e,size(dmft_bath%e),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIerr)
    call MPI_BCAST(dmft_bath%v,size(dmft_bath%v),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIerr)
#endif
  end subroutine init_bath_ed




  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_bath(dmft_bath,unit)
    type(effective_bath) :: dmft_bath
    integer              :: i,unit,ispin,iorb
    if(.not.dmft_bath%status)stop "WRITE_BATH: bath not allocated"
    if(mpiID==0)then
       select case(bath_type)
       case default
          write(unit,"(90(A22,1X))")&
               (("# Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
               "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(F22.15,1X))")((dmft_bath%e(ispin,iorb,i),dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
          enddo

       case('hybrid')
          write(unit,"(90(A22,1X))")("# Ek_s"//reg(txtfy(ispin)),&
               ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(F22.15,1X))")(dmft_bath%e(ispin,1,i),&
                  (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
          enddo
       end select
    endif
  end subroutine write_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath_,dmft_bath)
    real(8),dimension(:) :: bath_
    type(effective_bath) :: dmft_bath
    integer              :: iorb,ispin,stride_spin,stride_orb
    integer              :: ei(2),vi(2)
    if(.not.dmft_bath%status)stop "SET_BATH: bath not allocated"
    call check_bath_dimension(bath_)
    select case(bath_type)
    case default
       do ispin=1,Nspin
          stride_spin=(ispin-1)*Norb*Nbath
          do iorb=1,Norb
             stride_orb =(iorb-1)*2*Nbath 
             ei(1)=stride_spin + stride_orb + 1
             ei(2)=stride_spin + stride_orb + Nbath
             vi(1)=stride_spin + stride_orb + Nbath + 1
             vi(2)=stride_spin + stride_orb + Nbath + Nbath
             dmft_bath%e(ispin,iorb,1:Nbath) = bath_(ei(1):ei(2)) 
             dmft_bath%v(ispin,iorb,1:Nbath) = bath_(vi(1):vi(2)) 
          enddo
       enddo

    case ('hybrid')
       do ispin=1,Nspin
          stride_spin=(ispin-1)*(Norb+1)*Nbath
          ei(1)=stride_spin + 1
          ei(2)=stride_spin + Nbath
          dmft_bath%e(ispin,1,1:Nbath) = bath_(ei(1):ei(2))
          do iorb=1,Norb             
             stride_orb = iorb*Nbath
             vi(1)=stride_spin + stride_orb + 1
             vi(2)=stride_spin + stride_orb + Nbath             
             dmft_bath%v(ispin,iorb,1:Nbath) = bath_(vi(1):vi(2))
          enddo
       enddo

    end select
  end subroutine set_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine copy_bath(dmft_bath,bath_)
    type(effective_bath) :: dmft_bath
    real(8),dimension(:) :: bath_
    integer              :: iorb,ispin
    integer              :: stride_spin,stride_orb
    integer              :: ei(2),vi(2)!ae,be,av,bv

    if(.not.dmft_bath%status)stop "COPY_BATH: bath not allocated"
    call check_bath_dimension(bath_)
    select case(bath_type)
    case default
       do ispin=1,Nspin
          stride_spin=(ispin-1)*Norb*Nbath
          do iorb=1,Norb
             stride_orb =(iorb-1)*2*Nbath 
             ei(1)=stride_spin + stride_orb + 1
             ei(2)=stride_spin + stride_orb + Nbath
             vi(1)=stride_spin + stride_orb + Nbath + 1
             vi(2)=stride_spin + stride_orb + Nbath + Nbath
             bath_(ei(1):ei(2)) = dmft_bath%e(ispin,iorb,1:Nbath)
             bath_(vi(1):vi(2)) = dmft_bath%v(ispin,iorb,1:Nbath)
          enddo
       enddo

    case ('hybrid')
       do ispin=1,Nspin
          stride_spin=(ispin-1)*(Norb+1)*Nbath
          ei(1)=stride_spin + 1
          ei(2)=stride_spin + Nbath
          bath_(ei(1):ei(2)) = dmft_bath%e(ispin,1,1:Nbath)
          do iorb=1,Norb             
             stride_orb = iorb*Nbath
             vi(1)=stride_spin + stride_orb + 1
             vi(2)=stride_spin + stride_orb + Nbath             
             bath_(vi(1):vi(2)) = dmft_bath%v(ispin,iorb,1:Nbath)
          enddo
       enddo

    end select

  end subroutine copy_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : given the bath array, compute the hybridization function
  ! for a given spin and orbital indices ispin,iorb at a given point x
  !+-------------------------------------------------------------------+
  function delta_and_irr(ispin,iorb,x,bath_) result(fg)
    integer,intent(in)                               :: ispin,iorb
    complex(8),intent(in)                            :: x
    real(8),dimension(2*Nspin*Norb*Nbath),intent(in) :: bath_
    complex(8)                                       :: fg
    integer                                          :: i
    type(effective_bath)                             :: dmft_bath_
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    fg=zero
    do i=1,Nbath
       fg=fg + dmft_bath_%v(ispin,iorb,i)**2/(x-dmft_bath_%e(ispin,iorb,i))
    enddo
    call deallocate_bath(dmft_bath_)
  end function delta_and_irr

  function delta_and_red(ispin,iorb,jorb,x,bath_) result(fg)
    integer,intent(in)                               :: ispin,iorb,jorb
    complex(8),intent(in)                            :: x
    real(8),dimension(Nspin*(Norb+1)*Nbath),intent(in) :: bath_
    complex(8)                                       :: fg
    integer                                          :: i
    type(effective_bath)                             :: dmft_bath_
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    fg=zero
    do i=1,Nbath
       fg=fg + dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)/(x-dmft_bath_%e(ispin,1,i))
    enddo
    call deallocate_bath(dmft_bath_)
  end function delta_and_red



  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function for a given spin and 
  ! orbital indices ispin and iorb at point x, from determined bath 
  ! components ebath,vbath
  !+-------------------------------------------------------------------+
  function delta_bath_irr(ispin,iorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,ispin
    complex(8)            :: fg
    fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2/(x-dmft_bath_%e(ispin,iorb,1:Nbath)))
  end function delta_bath_irr

  function delta_bath_red(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg
    fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)&
         /(x-dmft_bath_%e(ispin,1,1:Nbath)))
  end function delta_bath_red

END MODULE ED_BATH
