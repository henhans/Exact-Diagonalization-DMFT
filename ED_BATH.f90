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

  ! interface delta_and
  !    module procedure delta_and_irred,delta_and_irred_,delta_and_hybrd,delta_and_hybrd_
  ! end interface delta_and

  ! interface fdelta_and
  !    module procedure fdelta_and_irred,fdelta_and_irred_,fdelta_and_hybrd,fdelta_and_hybrd_
  ! end interface fdelta_and

  interface delta_bath
     module procedure delta_bath_irred,delta_bath_irred_,delta_bath_hybrd,delta_bath_hybrd_
  end interface delta_bath

  interface fdelta_bath
     module procedure fdelta_bath_irred,fdelta_bath_irred_,fdelta_bath_hybrd,fdelta_bath_hybrd_
  end interface fdelta_bath

  public :: allocate_bath
  public :: deallocate_bath
  public :: check_bath_dimension
  public :: init_bath_ed
  public :: get_bath_size
  public :: write_bath
  public :: set_bath
  public :: copy_bath
  public :: spin_symmetrize_bath
  !  public :: delta_and,fdelta_and
  public :: delta_bath,fdelta_bath
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
    case default
       if(.not.ed_supercond)then
          allocate(dmft_bath%e(Nspin,Norb,Nbath),dmft_bath%v(Nspin,Norb,Nbath))
       else
          allocate(&
               dmft_bath%e(Nspin,Norb,Nbath),&
               dmft_bath%v(Nspin,Norb,Nbath),&
               dmft_bath%d(Nspin,Norb,Nbath))
       endif
    case('hybrid')
       if(.not.ed_supercond)then
          allocate(dmft_bath%e(Nspin,1,Nbath),dmft_bath%v(Nspin,Norb,Nbath))
       else
          allocate(&
               dmft_bath%e(Nspin,1,Nbath),&
               dmft_bath%v(Nspin,Norb,Nbath),&
               dmft_bath%d(Nspin,1,Nbath))
       endif
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
    if(allocated(dmft_bath%d))deallocate(dmft_bath%d)
    dmft_bath%status=.false.
  end subroutine deallocate_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_size() result(dims)
    integer :: dims(2)
    dims(1)=Nspin
    select case(bath_type)
    case default
       if(.not.ed_supercond)then
          dims(2) = (2*Norb)*Nbath
       else
          dims(2) = (3*Norb)*Nbath
       endif
    case('hybrid')
       if(.not.ed_supercond)then
          dims(2) = (Norb+1)*Nbath
       else
          dims(2) = (Norb+2)*Nbath
       endif
    end select
  end function get_bath_size



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_) result(bool)
    real(8),dimension(:,:) :: bath_
    integer                :: N1_,N2_,Ntrue
    logical :: bool
    N1_=size(bath_,1)
    N2_=size(bath_,2)
    select case(bath_type)
    case default
       if(.not.ed_supercond)then
          Ntrue = (2*Norb)*Nbath
       else
          Ntrue = (3*Norb)*Nbath
       endif
    case('hybrid')
       if(.not.ed_supercond)then
          Ntrue = (Norb+1)*Nbath
       else
          Ntrue = (Norb+2)*Nbath
       endif
    end select
    bool = (N1_ == Nspin).AND.(N2_ == Ntrue)
  end function check_bath_dimension





  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed(dmft_bath,hwband_)
    type(effective_bath) :: dmft_bath
    real(8)              :: hwband_,wband_
    integer              :: i,iorb,ispin,unit,flen
    logical              :: IOfile
    real(8)              :: Nh,de !N2,di, 
    if(.not.dmft_bath%status)stop "init_bath: bath not allocated"
    if(mpiID==0)then
       !Generating the bath anyway, then you may want to read it to 
       !update some entries. This way you can restart even 
       !from different Ns calculation, this is better than 
       !start from a complete guess
       write(LOGfile,"(A)")"Generating bath:"
       dmft_bath%e(:,:,1)    =-hwband_ 
       dmft_bath%e(:,:,Nbath)= hwband_ 
       Nh=Nbath/2
       if(mod(Nbath,2)==0)then
          de=hwband_/dble(Nh-1)
          dmft_bath%e(:,:,Nh)  = -1.d-4
          dmft_bath%e(:,:,Nh+1)=  1.d-4
          do i=2,Nh-1
             dmft_bath%e(:,:,i)   =-hwband_ + (i-1)*de
             dmft_bath%e(:,:,Nh+i)= hwband_ - (i-1)*de
          enddo
       else
          de=hwband_/dble(Nh)
          dmft_bath%e(:,:,Nh+1)= 0.d0
          do i=2,Nh
             dmft_bath%e(:,:,i)   =-hwband_ + (i-1)*de
             dmft_bath%e(:,:,Nh+i)= hwband_ - (i-1)*de
          enddo
       endif
       do i=1,Nbath
          dmft_bath%v(:,:,i)=1.d-1/sqrt(dble(Nbath))
       enddo
       if(ed_supercond)then
          dmft_bath%d(:,:,:)=deltasc
       endif


       !Read from file if exist:

       !
       inquire(file=trim(Hfile),exist=IOfile)
       if(IOfile)then
          write(LOGfile,"(A)")'Reading bath from file'
          unit = free_unit()
          flen = file_length(trim(Hfile))
          open(unit,file=trim(Hfile))
          read(unit,*)
          select case(bath_type)
          case default
             if(.not.ed_supercond)then
                do i=1,min(flen,Nbath)
                   read(unit,*)((dmft_bath%e(ispin,iorb,i),&
                        dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
                enddo
             else
                do i=1,min(flen,Nbath)
                   read(unit,*)((dmft_bath%e(ispin,iorb,i),dmft_bath%d(ispin,iorb,i),&
                        dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
                enddo
             endif
          case ('hybrid')
             if(.not.ed_supercond)then
                do i=1,min(flen,Nbath)
                   read(unit,*)(dmft_bath%e(ispin,1,i),&
                        (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
                enddo
             else
                do i=1,min(flen,Nbath)
                   read(unit,*)(dmft_bath%e(ispin,1,i),dmft_bath%d(ispin,1,i),&
                        (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
                enddo
             endif
          end select
          close(unit)
       endif
    endif                       !mpiID==0
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
          if(.not.ed_supercond)then
             write(unit,"(90(A21,1X))")&
                  (("#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")((dmft_bath%e(ispin,iorb,i),&
                     dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          else
             write(unit,"(90(A21,1X))")&
                  (("#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),"#Dk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
                  "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")((dmft_bath%e(ispin,iorb,i),dmft_bath%d(ispin,iorb,i),&
                     dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          endif

       case('hybrid')
          if(.not.ed_supercond)then
             write(unit,"(90(A21,1X))")("#Ek_s"//reg(txtfy(ispin)),&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")( dmft_bath%e(ispin,1,i),&
                     (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          else
             write(unit,"(90(A21,1X))")("#Ek_s"//reg(txtfy(ispin)),"#Dk_s"//reg(txtfy(ispin)),&
                  ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
             do i=1,Nbath
                write(unit,"(90(F21.12,1X))")( dmft_bath%e(ispin,1,i),dmft_bath%d(ispin,1,i),&
                     (dmft_bath%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          endif
       end select
    endif
  end subroutine write_bath





  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath_,dmft_bath)
    real(8),dimension(:,:) :: bath_
    type(effective_bath)   :: dmft_bath
    integer                :: iorb,ispin,stride_spin,stride_orb
    integer                :: ei(2),vi(2),di(2)
    logical                :: check
    if(.not.dmft_bath%status)stop "SET_BATH: bath not allocated"
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_bath: wrong bath dimensions"
    select case(bath_type)
    case default
       if(.not.ed_supercond)then
          do ispin=1,Nspin
             do iorb=1,Norb
                stride_orb =(iorb-1)*2*Nbath
                ei(1)=stride_orb + 1
                ei(2)=stride_orb + Nbath
                vi(1)=stride_orb + Nbath + 1
                vi(2)=stride_orb + Nbath + Nbath
                dmft_bath%e(ispin,iorb,1:Nbath) = bath_(ispin,ei(1):ei(2)) 
                dmft_bath%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2)) 
             enddo
          enddo
       else
          do ispin=1,Nspin
             do iorb=1,Norb
                stride_orb =(iorb-1)*3*Nbath
                ei(1)=stride_orb + 1
                ei(2)=stride_orb + Nbath
                di(1)=stride_orb + Nbath + 1
                di(2)=stride_orb + Nbath + Nbath
                vi(1)=stride_orb + Nbath + Nbath + 1
                vi(2)=stride_orb + Nbath + Nbath + Nbath
                dmft_bath%e(ispin,iorb,1:Nbath) = bath_(ispin,ei(1):ei(2))
                dmft_bath%d(ispin,iorb,1:Nbath) = bath_(ispin,di(1):di(2))
                dmft_bath%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2))
             enddo
          enddo
       endif
    case ('hybrid')
       if(.not.ed_supercond)then
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             dmft_bath%e(ispin,1,1:Nbath) = bath_(ispin,ei(1):ei(2))
             do iorb=1,Norb             
                stride_orb = iorb*Nbath
                vi(1)=stride_orb + 1
                vi(2)=stride_orb + Nbath             
                dmft_bath%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2))
             enddo
          enddo
       else
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             di(1)=Nbath+1
             di(2)=Nbath+Nbath
             dmft_bath%e(ispin,1,1:Nbath) = bath_(ispin,ei(1):ei(2))
             dmft_bath%d(ispin,1,1:Nbath) = bath_(ispin,di(1):di(2))
             do iorb=1,Norb             
                stride_orb = (iorb+1)*Nbath
                vi(1)=stride_orb + 1
                vi(2)=stride_orb + Nbath             
                dmft_bath%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2))
             enddo
          enddo
       endif
    end select
  end subroutine set_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine copy_bath(dmft_bath,bath_)
    type(effective_bath)   :: dmft_bath
    real(8),dimension(:,:) :: bath_
    integer                :: iorb,ispin
    integer                :: stride_spin,stride_orb
    integer                :: ei(2),vi(2),di(2)
    logical                :: check
    if(.not.dmft_bath%status)stop "COPY_BATH: bath not allocated"
    check=check_bath_dimension(bath_)
    if(.not.check)stop "copy_bath: wrong bath dimensions"
    select case(bath_type)
    case default
       if(.not.ed_supercond)then
          do ispin=1,Nspin
             do iorb=1,Norb
                stride_orb =(iorb-1)*2*Nbath 
                ei(1)=stride_orb + 1
                ei(2)=stride_orb + Nbath
                vi(1)=stride_orb + Nbath + 1
                vi(2)=stride_orb + Nbath + Nbath
                bath_(ispin,ei(1):ei(2)) = dmft_bath%e(ispin,iorb,1:Nbath)
                bath_(ispin,vi(1):vi(2)) = dmft_bath%v(ispin,iorb,1:Nbath)
             enddo
          enddo
       else
          do ispin=1,Nspin
             do iorb=1,Norb
                stride_orb =(iorb-1)*3*Nbath
                ei(1)=stride_orb + 1
                ei(2)=stride_orb + Nbath
                di(1)=stride_orb + Nbath + 1
                di(2)=stride_orb + Nbath + Nbath
                vi(1)=stride_orb + Nbath + Nbath + 1
                vi(2)=stride_orb + Nbath + Nbath + Nbath
                bath_(ispin,ei(1):ei(2)) = dmft_bath%e(ispin,iorb,1:Nbath)  
                bath_(ispin,di(1):di(2)) = dmft_bath%d(ispin,iorb,1:Nbath)  
                bath_(ispin,vi(1):vi(2)) = dmft_bath%v(ispin,iorb,1:Nbath)  
             enddo
          enddo
       endif
    case ('hybrid')
       if(.not.ed_supercond)then
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             bath_(ispin,ei(1):ei(2)) = dmft_bath%e(ispin,1,1:Nbath)
             do iorb=1,Norb             
                stride_orb = iorb*Nbath
                vi(1) = stride_orb + 1
                vi(2) = stride_orb + Nbath             
                bath_(ispin,vi(1):vi(2)) = dmft_bath%v(ispin,iorb,1:Nbath)
             enddo
          enddo
       else
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             di(1)=Nbath+1
             di(2)=Nbath+Nbath
             bath_(ispin,ei(1):ei(2)) = dmft_bath%e(ispin,1,1:Nbath)  
             bath_(ispin,di(1):di(2)) = dmft_bath%d(ispin,1,1:Nbath)  
             do iorb=1,Norb             
                stride_orb = (iorb+1)*Nbath
                vi(1)=stride_orb + 1
                vi(2)=stride_orb + Nbath             
                bath_(ispin,vi(1):vi(2)) = dmft_bath%v(ispin,iorb,1:Nbath)
             enddo
          enddo
       endif
    end select
  end subroutine copy_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array set both spin components to have 
  !the same bath, i.e. impose non-magnetic solution
  !+-------------------------------------------------------------------+
  subroutine spin_symmetrize_bath(bath_)
    real(8),dimension(:,:) :: bath_
    type(effective_bath)   :: dmft_bath_
    if(Nspin==1)then
       write(LOGfile,"(A)")"spin_symmetrize_bath: Nspin=1 nothing to symmetrize"
       return
    endif
    if(Nspin>2)stop "spin_symmetrize_bath: Nspin>2..."
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(1,:,:)
    dmft_bath_%v(Nspin,:,:)=dmft_bath_%v(1,:,:)
    if(ed_supercond)dmft_bath_%d(Nspin,:,:)=dmft_bath_%d(1,:,:)
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine spin_symmetrize_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : given the bath array, compute the hybridization function
  ! for a given spin and orbital indices ispin,iorb at a given point x
  !+-------------------------------------------------------------------+
  ! function delta_and_irred(ispin,iorb,x,bath_) result(fg)
  !   integer,intent(in)                               :: ispin,iorb
  !   complex(8),intent(in)                            :: x
  !   real(8),dimension(Nspin,2*Norb*Nbath),intent(in) :: bath_
  !   complex(8)                                       :: fg
  !   integer                                          :: i
  !   type(effective_bath)                             :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   if(.not.ed_supercond)then
  !      do i=1,Nbath
  !         fg=fg + dmft_bath_%v(ispin,iorb,i)**2/(x-dmft_bath_%e(ispin,iorb,i))
  !      enddo
  !   else
  !      do i=1,Nbath
  !         fg=fg - dmft_bath_%v(ispin,iorb,i)**2*(x+dmft_bath_%e(ispin,iorb,i))/&
  !              (dimag(x)**2+dmft_bath_%e(ispin,iorb,i)**2+dmft_bath_%d(ispin,iorb,i)**2)
  !      enddo
  !   endif
  !   call deallocate_bath(dmft_bath_)
  ! end function delta_and_irred

  ! function delta_and_irred_(ispin,iorb,w,eps,bath_) result(fg)
  !   integer,intent(in)                               :: ispin,iorb
  !   real(8),intent(in)                               :: w,eps
  !   real(8),dimension(Nspin,2*Norb*Nbath),intent(in) :: bath_
  !   complex(8)                                       :: fg,x
  !   integer                                          :: i
  !   type(effective_bath)                             :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   x=dcmplx(w,eps)
  !   if(.not.ed_supercond)then
  !      do i=1,Nbath
  !         fg=fg + dmft_bath_%v(ispin,iorb,i)**2/(x-dmft_bath_%e(ispin,iorb,i))
  !      enddo
  !   else
  !      do i=1,Nbath
  !         fg=fg + dmft_bath_%v(ispin,iorb,i)**2*(x+dmft_bath_%e(ispin,iorb,i))/&
  !              (w**2+xi*eps-dmft_bath_%e(ispin,iorb,i)**2-dmft_bath_%d(ispin,iorb,i)**2)
  !      enddo
  !   endif
  !   call deallocate_bath(dmft_bath_)
  ! end function delta_and_irred_

  ! function fdelta_and_irred(ispin,iorb,x,bath_) result(fg)
  !   integer,intent(in)                               :: ispin,iorb
  !   complex(8),intent(in)                            :: x
  !   real(8),dimension(Nspin,2*Norb*Nbath),intent(in) :: bath_
  !   complex(8)                                       :: fg
  !   integer                                          :: i
  !   type(effective_bath)                             :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   do i=1,Nbath
  !      fg=fg + dmft_bath_%d(ispin,iorb,i)*dmft_bath_%v(ispin,iorb,i)**2/&
  !           (dimag(x)**2+dmft_bath_%e(ispin,iorb,i)**2+dmft_bath_%d(ispin,iorb,i)**2)
  !   enddo
  !   call deallocate_bath(dmft_bath_)
  ! end function fdelta_and_irred

  ! function fdelta_and_irred_(ispin,iorb,w,eps,bath_) result(fg)
  !   integer,intent(in)                               :: ispin,iorb
  !   real(8),intent(in)                               :: w,eps
  !   real(8),dimension(Nspin,2*Norb*Nbath),intent(in) :: bath_
  !   complex(8)                                       :: fg
  !   integer                                          :: i
  !   type(effective_bath)                             :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   do i=1,Nbath
  !      fg=fg + dmft_bath_%d(ispin,iorb,i)*dmft_bath_%v(ispin,iorb,i)**2/&
  !           (w**2+xi*eps+dmft_bath_%e(ispin,iorb,i)**2+dmft_bath_%d(ispin,iorb,i)**2)
  !   enddo
  !   call deallocate_bath(dmft_bath_)
  ! end function fdelta_and_irred_







  ! function delta_and_hybrd(ispin,iorb,jorb,x,bath_) result(fg)
  !   integer,intent(in)                                 :: ispin,iorb,jorb
  !   complex(8),intent(in)                              :: x
  !   real(8),dimension(Nspin,(Norb+1)*Nbath),intent(in) :: bath_
  !   complex(8)                                         :: fg
  !   integer                                            :: i
  !   type(effective_bath)                               :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   if(.not.ed_supercond)then
  !      do i=1,Nbath
  !         fg=fg + dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)/(x-dmft_bath_%e(ispin,1,i))
  !      enddo
  !   else
  !      do i=1,Nbath
  !         fg=fg - dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)*(x+dmft_bath_%e(ispin,1,i))/&
  !              (dimag(x)**2+dmft_bath_%e(ispin,1,i)**2+dmft_bath_%d(ispin,1,i)**2)
  !      enddo
  !   endif
  !   call deallocate_bath(dmft_bath_)
  ! end function delta_and_hybrd

  ! function delta_and_hybrd_(ispin,iorb,jorb,w,eps,bath_) result(fg)
  !   integer,intent(in)                                 :: ispin,iorb,jorb
  !   real(8),intent(in)                                 :: w,eps
  !   real(8),dimension(Nspin,(Norb+1)*Nbath),intent(in) :: bath_
  !   complex(8)                                         :: fg,x
  !   integer                                            :: i
  !   type(effective_bath)                               :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   x=dcmplx(w,eps)
  !   if(.not.ed_supercond)then
  !      do i=1,Nbath
  !         fg=fg + dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)/(x-dmft_bath_%e(ispin,1,i))
  !      enddo
  !   else
  !      do i=1,Nbath
  !         fg=fg + dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)*(x+dmft_bath_%e(ispin,1,i))/&
  !              (w**2+xi*eps-dmft_bath_%e(ispin,1,i)**2-dmft_bath_%d(ispin,1,i)**2)
  !      enddo
  !   endif
  !   call deallocate_bath(dmft_bath_)
  ! end function delta_and_hybrd_

  ! function fdelta_and_hybrd(ispin,iorb,jorb,x,bath_) result(fg)
  !   integer,intent(in)                                 :: ispin,iorb,jorb
  !   complex(8),intent(in)                              :: x
  !   real(8),dimension(Nspin,(Norb+1)*Nbath),intent(in) :: bath_
  !   complex(8)                                         :: fg
  !   integer                                            :: i
  !   type(effective_bath)                               :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   do i=1,Nbath
  !      fg=fg + dmft_bath_%d(ispin,1,i)*dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)/&
  !           (dimag(x)**2+dmft_bath_%e(ispin,1,i)**2+dmft_bath_%d(ispin,1,i)**2)
  !   enddo
  !   call deallocate_bath(dmft_bath_)
  ! end function fdelta_and_hybrd

  ! function fdelta_and_hybrd_(ispin,iorb,jorb,w,eps,bath_) result(fg)
  !   integer,intent(in)                                 :: ispin,iorb,jorb
  !   real(8),intent(in)                                 :: w,eps
  !   real(8),dimension(Nspin,(Norb+1)*Nbath),intent(in) :: bath_
  !   complex(8)                                         :: fg
  !   integer                                            :: i
  !   type(effective_bath)                               :: dmft_bath_
  !   call allocate_bath(dmft_bath_)
  !   call set_bath(bath_,dmft_bath_)
  !   fg=zero
  !   do i=1,Nbath
  !      fg=fg + dmft_bath_%d(ispin,1,i)*dmft_bath_%v(ispin,iorb,i)*dmft_bath_%v(ispin,jorb,i)/&
  !           (w**2+xi*eps+dmft_bath_%e(ispin,1,i)**2+dmft_bath_%d(ispin,1,i)**2)
  !   enddo
  !   call deallocate_bath(dmft_bath_)
  ! end function fdelta_and_hybrd_











  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function for a given spin and 
  ! orbital indices ispin and iorb at point x, from determined bath 
  ! components ebath,vbath
  !+-------------------------------------------------------------------+
  function delta_bath_irred(ispin,iorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,ispin
    complex(8)            :: fg
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2/(x-dmft_bath_%e(ispin,iorb,1:Nbath)))
    else
       fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2*(x+dmft_bath_%e(ispin,iorb,1:Nbath))/&
            (dimag(x)**2+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
    endif
  end function delta_bath_irred

  function delta_bath_irred_(ispin,iorb,w,eps,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    real(8),intent(in) :: w,eps
    integer,intent(in)    :: iorb,ispin
    complex(8)            :: fg,x
    x=dcmplx(w,eps)
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2/(x-dmft_bath_%e(ispin,iorb,1:Nbath)))
    else
       ! fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2*(x+dmft_bath_%e(ispin,iorb,1:Nbath))/&
       !      (w**2+xi*eps+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2*(x+dmft_bath_%e(ispin,iorb,1:Nbath))/&
            ( x*(-x)+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
    endif
  end function delta_bath_irred_

  function fdelta_bath_irred(ispin,iorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,ispin
    complex(8)            :: fg
    fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
         (dimag(x)**2+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  end function fdelta_bath_irred

  function fdelta_bath_irred_(ispin,iorb,w,eps,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    real(8),intent(in) :: w,eps
    integer,intent(in)    :: iorb,ispin
    complex(8)            :: fg,x
    x=dcmplx(w,eps)
    ! fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
    !      (w**2+xi*eps+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
    fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
         ( x*(-x) + dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  end function fdelta_bath_irred_









  function delta_bath_hybrd(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)&
            /(x-dmft_bath_%e(ispin,1,1:Nbath)))
    else
       fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)*(x+dmft_bath_%e(ispin,1,1:Nbath))/&
            (dimag(x)**2+dmft_bath_%e(ispin,1,1:Nbath)**2+dmft_bath_%d(ispin,1,1:Nbath)**2))
    endif
  end function delta_bath_hybrd

  function delta_bath_hybrd_(ispin,iorb,jorb,w,eps,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    real(8),intent(in)    :: w,eps
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg,x
    x=dcmplx(w,eps)
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)&
            /(x-dmft_bath_%e(ispin,1,1:Nbath)))
    else
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)*(x+dmft_bath_%e(ispin,1,1:Nbath))/&
            (w**2+xi*eps+dmft_bath_%e(ispin,1,1:Nbath)**2+dmft_bath_%d(ispin,1,1:Nbath)**2))
    endif
  end function delta_bath_hybrd_

  function fdelta_bath_hybrd(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg
    fg =sum(dmft_bath_%d(ispin,1,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/&
         (dimag(x)**2+dmft_bath_%e(ispin,1,1:Nbath)**2+dmft_bath_%d(ispin,1,1:Nbath)**2))
  end function fdelta_bath_hybrd

  function fdelta_bath_hybrd_(ispin,iorb,jorb,w,eps,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    real(8),intent(in)    :: w,eps
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg,x
    x=dcmplx(w,eps)
    fg =sum(dmft_bath_%d(ispin,1,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/&
         (w**2+xi*eps+dmft_bath_%e(ispin,1,1:Nbath)**2+dmft_bath_%d(ispin,1,1:Nbath)**2))
  end function fdelta_bath_hybrd_

END MODULE ED_BATH
