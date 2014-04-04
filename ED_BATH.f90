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
  USE COMMON_VARS
  USE IOTOOLS, only:free_unit,reg,file_length
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  implicit none

  private

  interface delta_bath_mats
     module procedure &
          delta_bath_irred_mats, & 
          delta_bath_hybrd_mats
  end interface delta_bath_mats

  interface delta_bath_real
     module procedure &
          delta_bath_irred_real, &
          delta_bath_hybrd_real
  end interface delta_bath_real

  interface fdelta_bath_mats
     module procedure &
          fdelta_bath_irred_mats, &
          fdelta_bath_hybrd_mats
  end interface fdelta_bath_mats

  interface fdelta_bath_real
     module procedure &
          fdelta_bath_irred_real, &
          fdelta_bath_hybrd_real
  end interface fdelta_bath_real


  public :: allocate_bath
  public :: deallocate_bath
  public :: check_bath_dimension
  public :: init_bath_ed
  public :: get_bath_size
  public :: write_bath
  public :: set_bath
  public :: copy_bath
  public :: spin_symmetrize_bath
  public :: break_symmetry_bath
  public :: delta_bath_mats
  public :: delta_bath_real
  public :: fdelta_bath_mats
  public :: fdelta_bath_real

contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(dmft_bath_%status)call deallocate_bath(dmft_bath_)
    select case(bath_type)
    case default
       if(.not.ed_supercond)then
          allocate(dmft_bath_%e(Nspin,Norb,Nbath),dmft_bath_%v(Nspin,Norb,Nbath))
       else
          allocate(&
               dmft_bath_%e(Nspin,Norb,Nbath),&
               dmft_bath_%v(Nspin,Norb,Nbath),&
               dmft_bath_%d(Nspin,Norb,Nbath))
       endif
    case('hybrid')
       if(.not.ed_supercond)then
          allocate(dmft_bath_%e(Nspin,1,Nbath),dmft_bath_%v(Nspin,Norb,Nbath))
       else
          allocate(&
               dmft_bath_%e(Nspin,1,Nbath),&
               dmft_bath_%v(Nspin,Norb,Nbath),&
               dmft_bath_%d(Nspin,1,Nbath))
       endif
    end select
    dmft_bath_%status=.true.
  end subroutine allocate_bath


  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_bath(dmft_bath_)
    type(effective_bath) :: dmft_bath_
    if(allocated(dmft_bath_%e))deallocate(dmft_bath_%e)
    if(allocated(dmft_bath_%v))deallocate(dmft_bath_%v)
    if(allocated(dmft_bath_%d))deallocate(dmft_bath_%d)
    dmft_bath_%status=.false.
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
  subroutine init_bath_ed(dmft_bath_,hwband_)
    type(effective_bath) :: dmft_bath_
    real(8)              :: hwband_,wband_
    integer              :: i,iorb,ispin,unit,flen
    logical              :: IOfile
    real(8)              :: Nh,de
    if(.not.dmft_bath_%status)stop "init_bath: bath not allocated"
    !Generating the bath anyway, then you may want to read it to 
    !update some entries. This way you can restart even 
    !from different Ns calculation, this is better than 
    !start from a complete guess
    dmft_bath_%e(:,:,1)    =-hwband_ 
    dmft_bath_%e(:,:,Nbath)= hwband_ 
    Nh=Nbath/2
    if(mod(Nbath,2)==0)then
       de=hwband_/dble(Nh-1)
       dmft_bath_%e(:,:,Nh)  = -1.d-4
       dmft_bath_%e(:,:,Nh+1)=  1.d-4
       do i=2,Nh-1
          dmft_bath_%e(:,:,i)   =-hwband_ + (i-1)*de
          dmft_bath_%e(:,:,Nh+i)= hwband_ - (i-1)*de
       enddo
    else
       de=hwband_/dble(Nh)
       dmft_bath_%e(:,:,Nh+1)= 0.d0
       do i=2,Nh
          dmft_bath_%e(:,:,i)   =-hwband_ + (i-1)*de
          dmft_bath_%e(:,:,Nh+i)= hwband_ - (i-1)*de
       enddo
    endif
    do i=1,Nbath
       dmft_bath_%v(:,:,i)=1.d-1/sqrt(dble(Nbath))
    enddo
    if(ed_supercond)then
       dmft_bath_%d(:,:,:)=deltasc
    endif
    !
    !Read from file if exist:
    !
    inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//".restart"
       unit = free_unit()
       flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
       open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
       read(unit,*)
       select case(bath_type)
       case default
          if(.not.ed_supercond)then
             do i=1,min(flen,Nbath)
                read(unit,*)((dmft_bath_%e(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          else
             do i=1,min(flen,Nbath)
                read(unit,*)((dmft_bath_%e(ispin,iorb,i),dmft_bath_%d(ispin,iorb,i),&
                     dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          endif
       case ('hybrid')
          if(.not.ed_supercond)then
             do i=1,min(flen,Nbath)
                read(unit,*)(dmft_bath_%e(ispin,1,i),&
                     (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          else
             do i=1,min(flen,Nbath)
                read(unit,*)(dmft_bath_%e(ispin,1,i),dmft_bath_%d(ispin,1,i),&
                     (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
             enddo
          endif
       end select
       close(unit)
    endif
  end subroutine init_bath_ed



  !+------------------------------------------------------------------+
  !PURPOSE  : breaking symmetry in a given direction using given field
  !+------------------------------------------------------------------+
  subroutine break_symmetry_bath(bath_,field,sign)
    real(8),dimension(:,:) :: bath_
    type(effective_bath)   :: dmft_bath_
    real(8)                :: field
    real(8)                :: sign
    call allocate_bath(dmft_bath_)
    call set_bath(bath_,dmft_bath_)
    dmft_bath_%e(1,:,:)    =dmft_bath_%e(1,:,:)      + sign*field
    dmft_bath_%e(Nspin,:,:)=dmft_bath_%e(Nspin,:,:)  - sign*field
    call copy_bath(dmft_bath_,bath_)
    call deallocate_bath(dmft_bath_)
  end subroutine break_symmetry_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_bath(dmft_bath_,unit)
    type(effective_bath) :: dmft_bath_
    integer              :: i,unit,ispin,iorb
    if(.not.dmft_bath_%status)stop "WRITE_BATH: bath not allocated"
    select case(bath_type)
    case default
       if(.not.ed_supercond)then
          write(unit,"(90(A21,1X))")&
               (("#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
               "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(F21.12,1X))")((dmft_bath_%e(ispin,iorb,i),&
                  dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
          enddo
       else
          write(unit,"(90(A21,1X))")&
               (("#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),"#Dk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
               "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(F21.12,1X))")((dmft_bath_%e(ispin,iorb,i),dmft_bath_%d(ispin,iorb,i),&
                  dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
          enddo
       endif

    case('hybrid')
       if(.not.ed_supercond)then
          write(unit,"(90(A21,1X))")("#Ek_s"//reg(txtfy(ispin)),&
               ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(F21.12,1X))")( dmft_bath_%e(ispin,1,i),&
                  (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
          enddo
       else
          write(unit,"(90(A21,1X))")("#Ek_s"//reg(txtfy(ispin)),"#Dk_s"//reg(txtfy(ispin)),&
               ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
          do i=1,Nbath
             write(unit,"(90(F21.12,1X))")( dmft_bath_%e(ispin,1,i),dmft_bath_%d(ispin,1,i),&
                  (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),ispin=1,Nspin)
          enddo
       endif
    end select
  end subroutine write_bath





  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_bath(bath_,dmft_bath_)
    real(8),dimension(:,:) :: bath_
    type(effective_bath)   :: dmft_bath_
    integer                :: iorb,ispin,stride_spin,stride_orb
    integer                :: ei(2),vi(2),di(2)
    logical                :: check
    if(.not.dmft_bath_%status)stop "SET_BATH: bath not allocated"
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
                dmft_bath_%e(ispin,iorb,1:Nbath) = bath_(ispin,ei(1):ei(2)) 
                dmft_bath_%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2)) 
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
                dmft_bath_%e(ispin,iorb,1:Nbath) = bath_(ispin,ei(1):ei(2))
                dmft_bath_%d(ispin,iorb,1:Nbath) = bath_(ispin,di(1):di(2))
                dmft_bath_%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2))
             enddo
          enddo
       endif
    case ('hybrid')
       if(.not.ed_supercond)then
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             dmft_bath_%e(ispin,1,1:Nbath) = bath_(ispin,ei(1):ei(2))
             do iorb=1,Norb             
                stride_orb = iorb*Nbath
                vi(1)=stride_orb + 1
                vi(2)=stride_orb + Nbath             
                dmft_bath_%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2))
             enddo
          enddo
       else
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             di(1)=Nbath+1
             di(2)=Nbath+Nbath
             dmft_bath_%e(ispin,1,1:Nbath) = bath_(ispin,ei(1):ei(2))
             dmft_bath_%d(ispin,1,1:Nbath) = bath_(ispin,di(1):di(2))
             do iorb=1,Norb             
                stride_orb = (iorb+1)*Nbath
                vi(1)=stride_orb + 1
                vi(2)=stride_orb + Nbath             
                dmft_bath_%v(ispin,iorb,1:Nbath) = bath_(ispin,vi(1):vi(2))
             enddo
          enddo
       endif
    end select
  end subroutine set_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine copy_bath(dmft_bath_,bath_)
    type(effective_bath)   :: dmft_bath_
    real(8),dimension(:,:) :: bath_
    integer                :: iorb,ispin
    integer                :: stride_spin,stride_orb
    integer                :: ei(2),vi(2),di(2)
    logical                :: check
    if(.not.dmft_bath_%status)stop "COPY_BATH: bath not allocated"
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
                bath_(ispin,ei(1):ei(2)) = dmft_bath_%e(ispin,iorb,1:Nbath)
                bath_(ispin,vi(1):vi(2)) = dmft_bath_%v(ispin,iorb,1:Nbath)
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
                bath_(ispin,ei(1):ei(2)) = dmft_bath_%e(ispin,iorb,1:Nbath)  
                bath_(ispin,di(1):di(2)) = dmft_bath_%d(ispin,iorb,1:Nbath)  
                bath_(ispin,vi(1):vi(2)) = dmft_bath_%v(ispin,iorb,1:Nbath)  
             enddo
          enddo
       endif
    case ('hybrid')
       if(.not.ed_supercond)then
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             bath_(ispin,ei(1):ei(2)) = dmft_bath_%e(ispin,1,1:Nbath)
             do iorb=1,Norb             
                stride_orb = iorb*Nbath
                vi(1) = stride_orb + 1
                vi(2) = stride_orb + Nbath             
                bath_(ispin,vi(1):vi(2)) = dmft_bath_%v(ispin,iorb,1:Nbath)
             enddo
          enddo
       else
          do ispin=1,Nspin
             ei(1)=1
             ei(2)=Nbath
             di(1)=Nbath+1
             di(2)=Nbath+Nbath
             bath_(ispin,ei(1):ei(2)) = dmft_bath_%e(ispin,1,1:Nbath)  
             bath_(ispin,di(1):di(2)) = dmft_bath_%d(ispin,1,1:Nbath)  
             do iorb=1,Norb             
                stride_orb = (iorb+1)*Nbath
                vi(1)=stride_orb + 1
                vi(2)=stride_orb + Nbath             
                bath_(ispin,vi(1):vi(2)) = dmft_bath_%v(ispin,iorb,1:Nbath)
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
  !PURPOSE  : compute the hybridization function for a given spin and 
  ! orbital indices ispin and iorb at point x, from determined bath 
  ! components ebath,vbath
  !+-------------------------------------------------------------------+


  !NORMAL/IRREDUCIBLE BATH:
  !Matsubara:
  function delta_bath_irred_mats(ispin,iorb,x,dmft_bath_) result(fg)
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
  end function delta_bath_irred_mats
  !
  function fdelta_bath_irred_mats(ispin,iorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,ispin
    complex(8)            :: fg
    fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
         (dimag(x)**2+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  end function fdelta_bath_irred_mats



  !NORMAL/IRREDUCIBLE BATH:
  !Real axis:
  function delta_bath_irred_real(ispin,iorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    integer,intent(in)    :: iorb,ispin
    complex(8),intent(in) :: x
    complex(8)            :: fg
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2/(x-dmft_bath_%e(ispin,iorb,1:Nbath)))
    else
       fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)**2*(x+dmft_bath_%e(ispin,iorb,1:Nbath))/&
            ( x*(-x)+dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
    endif
  end function delta_bath_irred_real
  !
  function fdelta_bath_irred_real(ispin,iorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    integer,intent(in)    :: iorb,ispin
    complex(8),intent(in) :: x
    complex(8)            :: fg
    fg = sum(dmft_bath_%d(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)**2/&
         ( x*(-x) + dmft_bath_%e(ispin,iorb,1:Nbath)**2+dmft_bath_%d(ispin,iorb,1:Nbath)**2))
  end function fdelta_bath_irred_real




  !<ACTHUNG/TODO: extend the irred expressions for the Delta to the hybrid case.
  !HYBRIDIZED BATH:
  !Matsubara:
  function delta_bath_hybrd_mats(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/(x - dmft_bath_%e(ispin,1,1:Nbath)))
    else
       fg = -sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)*(x + dmft_bath_%e(ispin,1,1:Nbath))/&
            (dimag(x)**2 + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
    endif
  end function delta_bath_hybrd_mats
  !
  function fdelta_bath_hybrd_mats(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    complex(8),intent(in) :: x
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg
    fg =sum(dmft_bath_%d(ispin,1,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/&
         ( dimag(x)**2 + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
  end function fdelta_bath_hybrd_mats



  !HYBRIDIZED BATH:
  !Real-axis:
  function delta_bath_hybrd_real(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg,x
    if(.not.ed_supercond)then
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)&
            /(x-dmft_bath_%e(ispin,1,1:Nbath)))
    else
       fg = sum(dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)*(x+dmft_bath_%e(ispin,1,1:Nbath))/&
            ( -x**2 + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
    endif
  end function delta_bath_hybrd_real
  !
  function fdelta_bath_hybrd_real(ispin,iorb,jorb,x,dmft_bath_) result(fg)
    type(effective_bath)  :: dmft_bath_
    integer,intent(in)    :: iorb,jorb,ispin
    complex(8)            :: fg,x
    fg =sum(dmft_bath_%d(ispin,1,1:Nbath)*dmft_bath_%v(ispin,iorb,1:Nbath)*dmft_bath_%v(ispin,jorb,1:Nbath)/&
         ( -x**2  + dmft_bath_%e(ispin,1,1:Nbath)**2 + dmft_bath_%d(ispin,1,1:Nbath)**2))
  end function fdelta_bath_hybrd_real










END MODULE ED_BATH
