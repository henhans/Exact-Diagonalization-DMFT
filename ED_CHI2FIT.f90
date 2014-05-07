!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!STRUCTURE: 
!
! chi2fit_interface
!
! fit_irred
!  + chi2_anderson/functions
!  + chi2_weiss/functions
!
! fit_irred_SC
!  + chi2_anderson_sc/functions
!  + chi2_weiss_sc/functinos
!
! fit_hybrd
!  + chi2_anderson/functions
!  + chi2_weiss/functions
!
!########################################################################
MODULE ED_CHI2FIT
  USE CONSTANTS
  USE OPTIMIZE, only:fmin_cg
  USE MATRIX,   only:matrix_inverse
  USE IOTOOLS,  only:reg,free_unit,txtfy
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH

  implicit none
  private

  interface chi2_fitgf
     module procedure chi2_fitgf_,chi2_fitgf__
  end interface chi2_fitgf

  public :: chi2_fitgf

  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb
  integer,dimension(:),allocatable      :: getIorb,getJorb
  integer                               :: Orb_indx,Spin_indx


contains



  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta general interface
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_(fg,bath,ispin)
    complex(8),dimension(:,:,:)          :: fg
    real(8),dimension(:,:),intent(inout) :: bath
    integer                              :: ispin
    select case(bath_type)
    case default
       call chi2_fitgf_irred(fg,bath,ispin)
    case ('hybrid')
       call chi2_fitgf_hybrd(fg,bath,ispin)
    end select
  end subroutine chi2_fitgf_

  subroutine chi2_fitgf__(fg,bath,ispin)
    complex(8),dimension(:,:,:,:)          :: fg
    real(8),dimension(:,:),intent(inout) :: bath
    integer                              :: ispin
    select case(bath_type)
    case default
       call chi2_fitgf_irred_sc(fg,bath,ispin)
    case ('hybrid')
       stop 'Error: Hybrid bath + SC is not implemented yet: ask the developer...'
    end select
  end subroutine chi2_fitgf__






  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************








  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Irreducible bath normal phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_irred(fg,bath_,ispin)
    complex(8),dimension(:,:,:)          :: fg
    real(8),dimension(:,:),intent(inout) :: bath_
    integer                              :: ispin
    real(8),dimension(2*Nbath)           :: a
    integer                              :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
    real(8)                              :: chi
    logical                              :: check
    type(effective_bath)                 :: dmft_bath
    complex(8)                           :: fgand
    real(8)                              :: w
    character(len=20)                  :: suffix
    integer                            :: unit

    if(size(fg,1)/=Norb)stop"CHI2_FITGF: wrong dimension 1 in chi2_input"
    if(size(fg,2)/=Norb)stop"CHI2_FITGF: wrong dimension 2 in chi2_input"
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
    Ldelta = Lfit
    if(Ldelta>size(fg,3))Ldelta=size(fg,3)
    !
    allocate(Fdelta(1,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
    !
    select case(Cg_weight)
    case default
       Wdelta=dble(Ldelta)
    case(1)
       Wdelta=1.d0
    case(2)
       Wdelta=(/(real(i,8),i=1,Ldelta)/)
    case(3)
       Wdelta=Xdelta
    end select
    !
    call allocate_bath(dmft_bath)
    call set_bath(bath_,dmft_bath)
    do iorb=1,Norb
       Orb_indx=iorb
       Spin_indx=ispin
       Fdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
       a(1:Nbath)         = dmft_bath%e(ispin,iorb,1:Nbath)
       a(Nbath+1:2*Nbath) = dmft_bath%v(ispin,iorb,1:Nbath)
       if(cg_scheme=='weiss')then
          call fmin_cg(a,chi2_weiss_irred,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.)
       else
          call fmin_cg(a,chi2_delta_irred,dchi2_delta_irred,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.)
       endif
       if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)") 'chi^2|iter'//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
       if(ed_verbose<2)then
          suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
          unit=free_unit()
          open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
          write(unit,"(ES18.9,1x,I5)") chi,iter
          close(unit)
       endif
       dmft_bath%e(ispin,iorb,1:Nbath) = a(1:Nbath)
       dmft_bath%v(ispin,iorb,1:Nbath) = a(Nbath+1:2*Nbath)
    enddo
    if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
    unit=free_unit()
    open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
    call write_bath(dmft_bath,unit)
    close(unit)
    if(ed_verbose<3)call write_fit_result(ispin)
    call copy_bath(dmft_bath,bath_)
    call deallocate_bath(dmft_bath)
    deallocate(Fdelta,Xdelta,Wdelta)
    !
  contains
    !
    subroutine write_fit_result(ispin)
      complex(8)        :: fgand
      integer           :: i,j,iorb,ispin
      real(8)           :: w
      do iorb=1,Norb
         suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
         Fdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
         fgand=zero
         unit=free_unit()
         open(unit,file="fit_delta"//reg(suffix)//".ed")
         do i=1,Ldelta
            w = Xdelta(i)
            if(cg_scheme=='weiss')then
               fgand = xi*w + xmu - hloc(ispin,ispin,iorb,iorb) - delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
               fgand = one/fgand
            else
               fgand = delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
            endif
            write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand),dreal(Fdelta(1,i)),dreal(fgand)
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_irred



  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  ! The gradient is evaluated analytically d\chi^2.
  ! \Delta_Anderson  and its gradient d\Delta_Anderson are evaluated 
  ! as separate functions.  
  ! NORMAL bath.
  ! SPIN & ORBITAL DIAGONAL
  !+-------------------------------------------------------------+
  function chi2_delta_irred(a) result(chi2)
    real(8),dimension(:)         ::  a
    complex(8),dimension(Ldelta) ::  g0
    real(8)                      ::  chi2
    integer                      ::  i,iorb
    chi2 = 0.d0 
    iorb=Orb_indx
    do i=1,Ldelta
       g0(i)   = fg_delta_irred(xdelta(i),iorb,a)
    enddo
    chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
  end function chi2_delta_irred
  ! the analytic GRADIENT of \chi^2
  !+-------------------------------------------------------------+
  function dchi2_delta_irred(a) result(dchi2)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2
    real(8),dimension(size(a))           :: df
    complex(8),dimension(Ldelta)         :: g0
    complex(8),dimension(Ldelta,size(a)) :: dg0
    integer                              :: i,j,iorb
    df=0.d0
    iorb=Orb_indx
    do i=1,Ldelta
       g0(i)    = fg_delta_irred(xdelta(i),iorb,a)
       dg0(i,:) = grad_fg_delta_irred(xdelta(i),iorb,a)
    enddo
    do j=1,size(a)
       df(j)=sum( dreal(Fdelta(1,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
            sum(  dimag(Fdelta(1,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
    enddo
    dchi2 = -2.d0*df
  end function dchi2_delta_irred
  ! the \Delta_Anderson function used in \chi^2 and d\chi^2
  ! \Delta = \sum_l V_l^2/(iw-e_l)
  !+-------------------------------------------------------------+
  function fg_delta_irred(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)                   :: gg,delta
    integer                      :: i,iorb,Nb
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    delta=sum(vps(:)**2/(xi*w-eps(:)))
    gg=delta
  end function fg_delta_irred
  ! the gradient d\Delta_Anderson function used in d\chi^2
  ! d\Delta = \grad_{V_k,e_k}\sum_l V_l^2/(iw-e_l)
  !+-------------------------------------------------------------+
  function grad_fg_delta_irred(w,iorb,a) result(dgz)
    real(8)                         :: w,sgn
    real(8),dimension(:)            :: a
    real(8),dimension(size(a)/2)    :: eps,vps
    complex(8),dimension(size(a))   :: dgz
    complex(8)                      :: gg,delta
    integer                         :: i,iorb,Nb
    dgz=zero
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    do i=1,Nb
       dgz(i)    = vps(i)*vps(i)/(xi*w-eps(i))**2
       dgz(i+Nb) = 2.d0*vps(i)/(xi*w-eps(i))
    enddo
  end function grad_fg_delta_irred


  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distanec of G_0  function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  ! G_0 is evaluated in a function below
  ! NORMAL bath.
  ! SPIN & ORBITAL DIAGONAL
  !+-------------------------------------------------------------+
  function chi2_weiss_irred(a) result(chi2)
    real(8),dimension(:)         ::  a
    complex(8),dimension(Ldelta) ::  g0
    real(8)                      ::  chi2
    integer                      ::  i,iorb
    chi2 = 0.d0 
    iorb=Orb_indx
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = fg_weiss_irred(xdelta(i),iorb,a)
    enddo
    chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
  end function chi2_weiss_irred
  ! the inverse non-interacting GF (~inverse Weiss Field) 
  ! used in \chi^2(\caG_0 - G_0) 
  ! \G_0 = [iw_n + xmu - H_loc - \Delta]^-1
  !+-------------------------------------------------------------+
  function fg_weiss_irred(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    complex(8)                   :: gg,delta
    integer                      :: i,iorb,ispin
    ispin=Spin_indx
    delta=fg_delta_irred(w,iorb,a)
    gg = one/(xi*w+xmu-hloc(ispin,ispin,iorb,iorb)-delta)
  end function fg_weiss_irred









  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************









  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Irreducible bath Superconducting phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_irred_sc(fg,bath_,ispin)
    complex(8),dimension(:,:,:,:)        :: fg
    real(8),dimension(:,:),intent(inout) :: bath_
    integer                              :: ispin
    real(8),dimension(3*Nbath)           :: a
    integer                              :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
    real(8)                              :: chi
    logical                              :: check
    type(effective_bath)                 :: dmft_bath
    complex(8)                           :: fgand
    real(8)                              :: w
    character(len=20)                    :: suffix
    integer                              :: unit
    if(size(fg,1)/=2)stop"CHI2_FITGF: wrong dimension 1 in chi2_input"
    if(size(fg,2)/=Norb)stop"CHI2_FITGF: wrong dimension 2 in chi2_input"
    if(size(fg,3)/=Norb)stop"CHI2_FITGF: wrong dimension 3 in chi2_input"
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
    Ldelta = Lfit
    if(Ldelta>size(fg,4))Ldelta=size(fg,4)
    !
    allocate(Fdelta(2,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
    !
    select case(Cg_weight)
    case default
       Wdelta=dble(Ldelta)
    case(1)
       Wdelta=1.d0
    case(2)
       Wdelta=(/(real(i,8),i=1,Ldelta)/)
    case(3)
       Wdelta=Xdelta
    end select
    !
    call allocate_bath(dmft_bath)
    call set_bath(bath_,dmft_bath)
    do iorb=1,Norb
       Orb_indx=iorb
       Spin_indx=ispin
       Fdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
       Fdelta(2,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
       a(1:Nbath)         = dmft_bath%e(ispin,iorb,1:Nbath) 
       a(Nbath+1:2*Nbath) = dmft_bath%d(ispin,iorb,1:Nbath)
       a(2*Nbath+1:3*Nbath) = dmft_bath%v(ispin,iorb,1:Nbath)
       if(cg_scheme=='weiss')then
          call fmin_cg(a,chi2_weiss_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.)
       else
          call fmin_cg(a,chi2_delta_irred_sc,dchi2_delta_irred_sc,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.)
       endif
       if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5,A)") 'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,"  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
       if(ed_verbose<2)then
          suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
          unit=free_unit()
          open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
          write(unit,"(ES18.9,1x,I5)") chi,iter
          close(unit)
       endif
       dmft_bath%e(ispin,iorb,1:Nbath) = a(1:Nbath)
       dmft_bath%d(ispin,iorb,1:Nbath) = a(Nbath+1:2*Nbath)
       dmft_bath%v(ispin,iorb,1:Nbath) = a(2*Nbath+1:3*Nbath)
    enddo
    if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
    unit=free_unit()
    open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
    call write_bath(dmft_bath,unit)
    close(unit)
    if(ed_verbose<3)call write_fit_result(ispin)
    call copy_bath(dmft_bath,bath_)
    call deallocate_bath(dmft_bath)
    deallocate(Fdelta,Xdelta,Wdelta)
    !
  contains
    !
    subroutine write_fit_result(ispin)
      complex(8)        :: fgand(2),det
      integer           :: i,j,iorb,ispin
      real(8)           :: w
      do iorb=1,Norb
         suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
         Fdelta(1,1:Ldelta) = fg(1,iorb,iorb,1:Ldelta)
         Fdelta(2,1:Ldelta) = fg(2,iorb,iorb,1:Ldelta)
         fgand=zero
         unit=free_unit()
         open(unit,file="fit_delta"//reg(suffix)//".ed")
         do i=1,Ldelta
            w = Xdelta(i)
            if(cg_scheme=='weiss')then
               fgand(1) = xi*w + xmu - hloc(ispin,ispin,iorb,iorb) - delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
               fgand(2) = -fdelta_bath_mats(ispin,iorb,xi*w,dmft_bath)
               det     =  abs(fgand(1))**2 + (fgand(2))**2
               fgand(1) = conjg(fgand(1))/det
               fgand(2) = fgand(2)/det
            else
               fgand(1) = delta_bath_mats(ispin,iorb,xi*w,dmft_bath)
               fgand(2) = fdelta_bath_mats(ispin,iorb,xi*w,dmft_bath)
            endif
            write(unit,"(10F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand(1)),dreal(Fdelta(1,i)),dreal(fgand(1)), &
                 dimag(Fdelta(2,i)),dimag(fgand(2)),dreal(Fdelta(2,i)),dreal(fgand(2))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_irred_sc










  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function 
  !         in the SUPERCONDUCTING case.
  ! The gradient is not evaluated.
  ! NORMAL bath.
  ! SPIN & ORBITAL DIAGONAL
  !+-------------------------------------------------------------+
  function chi2_delta_irred_sc(a) result(chi2_sc)
    real(8),dimension(:)           ::  a
    complex(8),dimension(2,Ldelta) ::  g0
    real(8)                        ::  chi2_sc
    integer                        ::  i,iorb
    chi2_sc = 0.d0 
    iorb=Orb_indx
    do i=1,Ldelta
       g0(:,i)   = fg_delta_irred_sc(xdelta(i),iorb,a)
    enddo
    chi2_sc=sum(abs(Fdelta(1,:)-g0(1,:))**2/Wdelta(:)) + sum(abs(Fdelta(2,:)-g0(2,:))**2/Wdelta(:)) 
  end function chi2_delta_irred_sc
  ! the analytic GRADIENT of \chi^2
  !+-------------------------------------------------------------+
  function dchi2_delta_irred_sc(a) result(dchi2_sc)
    real(8),dimension(:)                   ::  a
    complex(8),dimension(2,Ldelta)         ::  g0
    real(8),dimension(size(a))             ::  dchi2_sc,df
    complex(8),dimension(2,Ldelta,size(a)) ::  dg0
    integer                                ::  i,j,iorb
    dchi2_sc = 0.d0 
    df = 0.d0
    iorb=Orb_indx
    do i=1,Ldelta
       g0(:,i)   = fg_delta_irred_sc(xdelta(i),iorb,a)
       dg0(:,i,:)= grad_fg_delta_irred_sc(xdelta(i),iorb,a)
    enddo
    do j=1,size(a)
       df(j) = &
            sum( dreal(Fdelta(1,:)-g0(1,:))*dreal(dg0(1,:,j))/Wdelta(:) ) + sum(  dimag(Fdelta(1,:)-g0(1,:))*dimag(dg0(1,:,j))/Wdelta(:) ) +&
            sum( dreal(Fdelta(2,:)-g0(2,:))*dreal(dg0(2,:,j))/Wdelta(:) ) + sum(  dimag(Fdelta(2,:)-g0(2,:))*dimag(dg0(2,:,j))/Wdelta(:) )
    enddo
    dchi2_sc=-2.d0*df
  end function dchi2_delta_irred_sc

  ! the \Delta_Anderson function used in \chi^2 and d\chi^2
  ! \Delta = \sum_l V_l^2/(iw-e_l)
  !+-------------------------------------------------------------+
  function fg_delta_irred_sc(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    real(8),dimension(size(a)/3) :: eps,vps,dps
    complex(8)                   :: gg(2),delta(2),x
    integer                      :: i,iorb,Nb
    Nb=size(a)/3
    eps=a(1:Nb)
    dps=a(Nb+1:2*Nb)
    vps=a(2*Nb+1:3*Nb)
    gg = zero
    x = xi*w
    gg(1) = -sum(vps(:)**2*(x+eps(:))/(dimag(x)**2+eps(:)**2+dps(:)**2))
    gg(2) =  sum(dps(:)*vps(:)**2/(dimag(x)**2+eps(:)**2+dps(:)**2))
  end function fg_delta_irred_sc
  ! the gradient d\Delta_Anderson function used in d\chi^2
  ! d\Delta = \grad_{e_k,d_k,V_k}Delta(e_k,d_k,V_k)
  !+-------------------------------------------------------------+
  function grad_fg_delta_irred_sc(w,iorb,a) result(dgz)
    real(8)                         :: w
    real(8),dimension(:)            :: a
    real(8),dimension(size(a)/3)    :: eps,vps,dps
    complex(8),dimension(2,size(a)) :: dgz
    complex(8)                      :: x,den
    integer                         :: i,iorb,Nb
    Nb=size(a)/3
    eps=a(1:Nb)
    dps=a(Nb+1:2*Nb)
    vps=a(2*Nb+1:3*Nb)
    x = xi*w
    do i=1,Nb
       den = dimag(x)**2+eps(i)**2+dps(i)**2
       dgz(1,i) = -vps(i)*vps(i)*(1.d0/den - 2.d0*eps(i)*(x+eps(i))/den**2)
       dgz(1,i+Nb) = 2.d0*vps(i)*vps(i)*dps(i)*(x+eps(i))/den**2
       dgz(1,i+2*Nb) = -2.d0*vps(i)*(x+eps(i))/den
       !
       dgz(2,i) = -2.d0*vps(i)*vps(i)*dps(i)*eps(i)/den**2
       dgz(2,i+Nb) = vps(i)*vps(i)*(1.d0/den - 2.d0*dps(i)*dps(i)/den**2)
       dgz(2,i+2*Nb) = 2.d0*vps(i)*dps(i)/den
    enddo
  end function grad_fg_delta_irred_sc


  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of G_0 function 
  !         in the SUPERCONDUCTING case.
  ! The gradient is not evaluated.
  ! NORMAL bath.
  ! SPIN & ORBITAL DIAGONAL
  !+-------------------------------------------------------------+
  function chi2_weiss_irred_sc(a) result(chi2_sc)
    real(8),dimension(:)           ::  a
    complex(8),dimension(2,Ldelta) ::  g0
    real(8)                        ::  chi2_sc
    integer                        ::  i,iorb
    chi2_sc = 0.d0 
    iorb=Orb_indx
    do i=1,Ldelta
       g0(:,i)   = fg_weiss_irred_sc(xdelta(i),iorb,a)
    enddo
    chi2_sc = sum(abs(Fdelta(1,:)-g0(1,:))**2/Wdelta(:)) 
    chi2_sc = chi2_sc + sum(abs(Fdelta(2,:)-g0(2,:))**2/Wdelta(:)) 
  end function chi2_weiss_irred_sc
  ! the non interacting GF (~ inverse weiss field) 
  !+-------------------------------------------------------------+
  function fg_weiss_irred_sc(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    complex(8)                   :: gg(2),g0(2),delta(2),det
    integer                      :: i,iorb,ispin
    ispin = Spin_indx
    delta = fg_delta_irred_sc(w,iorb,a)
    g0(1) = xi*w + xmu - hloc(ispin,ispin,iorb,iorb) - delta(1)
    g0(2) = -delta(2)
    det   = abs(g0(1))**2 + (g0(2))**2
    gg(1) = conjg(g0(1))/det
    gg(2) = g0(2)/det
  end function fg_weiss_irred_sc






  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************






  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for Hybrid bath normal phase
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_hybrd(fg,bath_,ispin)
    complex(8),dimension(:,:,:)          :: fg
    real(8),dimension(:,:),intent(inout) :: bath_
    integer                              :: ispin
    real(8),dimension((1+Norb)*Nbath)    :: a
    integer                              :: iter,stride,ifirst,ilast,i,j,corb,l
    integer                              :: iorb,jorb
    real(8)                              :: chi
    complex(8),dimension(size(fg,3))     ::  g0
    logical                              :: check
    type(effective_bath)                 :: dmft_bath
    complex(8)                           :: fgand
    real(8)                              :: w
    character(len=20)                    :: suffix
    integer                              :: unit
    if(size(fg,1)/=Norb)stop "CHI2FIT: wrong dimension 1 in chi2_input"
    if(size(fg,2)/=Norb)stop "CHI2FIT: wrong dimension 2 in chi2_input"
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
    allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
    corb=0
    do iorb=1,Norb
       do jorb=iorb,Norb
          corb=corb+1
          getIorb(corb)=iorb
          getJorb(corb)=jorb
       enddo
    enddo
    totNorb=corb
    if(totNorb/=(Norb*(Norb+1)/2))stop "CHI2FIT: Error counting the orbitals"
    !
    Ldelta = Lfit
    if(Ldelta>size(fg,3))Ldelta=size(fg,3)
    allocate(Fdelta(totNorb,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    do i=1,totNorb
       Fdelta(i,1:Ldelta) = fg(getIorb(i),getJorb(i),1:Ldelta)
    enddo
    forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
    select case(Cg_weight)
    case default
       Wdelta=dble(Ldelta)
    case(1)
       Wdelta=1.d0
    case(2)
       Wdelta=(/(real(i,8),i=1,Ldelta)/)
    case(3)
       Wdelta=Xdelta
    end select
    !
    call allocate_bath(dmft_bath)
    Spin_indx=ispin
    a(:) = bath_(ispin,:)
    if(cg_scheme=='weiss')then
       call fmin_cg(a,chi2_weiss_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.)
    else
       call fmin_cg(a,chi2_delta_hybrd,dchi2_delta_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol,iverbose=.false.)
    endif
    bath_(ispin,:) = a(:)
    call set_bath(bath_,dmft_bath)
    if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5)") 'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,"  <--  all Orbs, Spin"//reg(txtfy(ispin))
    if(ed_verbose<2)then
       suffix="_ALLorb_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       unit=free_unit()
       open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
       write(unit,"(ES18.9,1x,I5)") chi,iter
       close(unit)
    endif
    if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
    unit=free_unit()
    open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
    call write_bath(dmft_bath,unit)
    close(unit)
    if(ed_verbose<3)call write_fit_result(ispin)
    call deallocate_bath(dmft_bath)
    deallocate(Fdelta,Xdelta,Wdelta)
    deallocate(getIorb,getJorb)
  contains
    subroutine write_fit_result(ispin)
      integer                              :: i,j,l,m,iorb,jorb,ispin,jspin
      real(8)                              :: w
      character(len=20)                    :: suffix
      integer                              :: unit
      complex(8),dimension(Norb,Norb)      :: gwf
      complex(8),dimension(totNorb,Ldelta) :: fgand
      do i=1,Ldelta
         w=Xdelta(i)
         if(cg_scheme=='weiss')then
            do l=1,Norb
               gwf(l,l)=xi*w + xmu - hloc(ispin,ispin,l,l) - delta_bath_mats(ispin,l,l,xi*w,dmft_bath)
               do m=l+1,Norb
                  gwf(l,m) = - hloc(ispin,ispin,l,m) - delta_bath_mats(ispin,l,m,xi*w,dmft_bath)
                  gwf(m,l) = - hloc(ispin,ispin,m,l) - delta_bath_mats(ispin,m,l,xi*w,dmft_bath)
               enddo
            enddo
            call matrix_inverse(gwf)
         else
            do l=1,Norb
               gwf(l,l)= delta_bath_mats(ispin,l,l,xi*w,dmft_bath)
               do m=l+1,Norb
                  gwf(l,m) = delta_bath_mats(ispin,l,m,xi*w,dmft_bath)
                  gwf(m,l) = delta_bath_mats(ispin,m,l,xi*w,dmft_bath)
               enddo
            enddo
         endif
         do l=1,totNorb
            iorb=getIorb(l)
            jorb=getJorb(l)
            fgand(l,i)=gwf(iorb,jorb)
         enddo
      enddo
      !
      do l=1,totNorb
         iorb=getIorb(l)
         jorb=getJorb(l)
         suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
         unit=free_unit()
         open(unit,file="fit_delta"//reg(suffix)//".ed")
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(l,i)),dimag(fgand(l,i)),&
                 dreal(Fdelta(l,i)),dreal(fgand(l,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
  end subroutine chi2_fitgf_hybrd




  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
  ! The gradient is evaluated analytically d\chi^2.
  ! \Delta_Anderson  and its gradient d\Delta_Anderson are evaluated 
  ! as separate functions.  
  ! HYBRID bath. 
  ! SPIN  DIAGONAL ; ORBITAL NON-DIAGONAL
  !+-------------------------------------------------------------+
  function chi2_delta_hybrd(a) result(chi2)
    real(8),dimension(:)                 :: a
    real(8),dimension(totNorb)           :: chi_orb
    complex(8),dimension(Ldelta)         :: Delta_orb
    real(8)                              :: chi2
    integer                              :: i,l,iorb,jorb
    chi_orb = 0.d0 
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       do i=1,Ldelta
          Delta_orb(i)= fg_delta_hybrd(xdelta(i),iorb,jorb,a)
       enddo
       chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
    enddo
    chi2=sum(chi_orb)
  end function chi2_delta_hybrd
  ! the analytic GRADIENT of \chi^2
  !+-------------------------------------------------------------+
  function dchi2_delta_hybrd(a) result(dchi2)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2
    real(8),dimension(totNorb,size(a))   :: df
    complex(8),dimension(Ldelta)         :: g0
    complex(8),dimension(Ldelta,size(a)) :: dg0
    integer                              :: i,j,l,iorb,jorb
    df=0.d0
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       do i=1,Ldelta
          g0(i)    = fg_delta_hybrd(xdelta(i),iorb,jorb,a)
          dg0(i,:) = grad_fg_delta_hybrd(xdelta(i),iorb,jorb,a)
       enddo
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Fdelta(l,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
               sum( dimag(Fdelta(l,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
       enddo
    enddo
    dchi2 = -2.d0*sum(df,1)     !sum over all orbital indices
  end function dchi2_delta_hybrd
  ! the \Delta_Anderson function used in \chi^2 and d\chi^2
  ! \Delta_ab = \sum_l V^a_l*V^b_l/(iw-e_l)
  !+-------------------------------------------------------------+
  function fg_delta_hybrd(w,orb1,orb2,a) result(gg)
    real(8)                      :: w
    integer                      :: orb1,orb2
    real(8),dimension(:)         :: a
    real(8),dimension(Nbath)     :: eps
    real(8),dimension(Norb,Nbath):: vps
    complex(8)                   :: gg,gwf(Norb,Norb)
    integer                      :: i,l,m
    eps=a(1:Nbath)
    do l=1,Norb
       vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
    enddo
    gg=zero
    do i=1,Nbath
       gg=gg + vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))
    enddo
  end function fg_delta_hybrd
  ! the gradient d\Delta_Anderson function used in d\chi^2
  ! d\Delta_ab = \grad_{V^c_k,e_k}\sum_l V^a_l*V^c_l/(iw-e_l)
  !+-------------------------------------------------------------+
  function grad_fg_delta_hybrd(w,orb1,orb2,a) result(dgz)
    real(8)                         :: w,sgn
    integer                         :: orb1,orb2
    real(8),dimension(:)            :: a
    real(8),dimension(Nbath)        :: eps
    real(8),dimension(Norb,Nbath)   :: vps
    complex(8),dimension(size(a))   :: dgz
    integer                         :: i,l
    dgz=zero
    eps=a(1:Nbath)
    do l=1,Norb
       vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
    enddo
    !
    do i=1,Nbath
       dgz(i)    = vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))**2
       if(orb1==orb2)then
          dgz(orb1*Nbath+i) = 2.d0*Vps(orb1,i)/(xi*w-eps(i))
       else
          dgz(orb1*Nbath+i) = Vps(orb2,i)/(xi*w-eps(i))
          dgz(orb2*Nbath+i) = Vps(orb1,i)/(xi*w-eps(i))
       endif
    enddo
  end function grad_fg_delta_hybrd




  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 distanec of G_0  function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient. 
  ! G_0 is evaluated in a function below
  ! HYBRID bath.
  ! SPIN DIAGONAL & ORBITAL OFF-DIAGONAL
  !+-------------------------------------------------------------+
  function chi2_weiss_hybrd(a) result(chi2)
    real(8),dimension(:)                 :: a
    real(8),dimension(totNorb)           :: chi_orb
    complex(8),dimension(Ldelta)         :: Delta_orb
    real(8)                              :: chi2
    integer                              :: i,l,iorb,jorb
    chi_orb = 0.d0 
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       do i=1,Ldelta
          Delta_orb(i)= fg_weiss_hybrd(xdelta(i),iorb,jorb,a)
       enddo
       chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
    enddo
    chi2=sum(chi_orb)
  end function chi2_weiss_hybrd
  !
  ! the inverse non-interacting GF (~inverse Weiss Field) 
  ! used in \chi^2(\caG_0 - G_0) 
  ! {\GG_0}_ab = {[iw_n + xmu - H_loc - \Delta]^-1}_ab
  !+-------------------------------------------------------------+
  function fg_weiss_hybrd(w,orb1,orb2,a) result(gg)
    real(8)                      :: w
    integer                      :: orb1,orb2
    real(8),dimension(:)         :: a
    complex(8)                   :: gg,gwf(Norb,Norb)
    integer                      :: i,l,m,ispin
    ispin=Spin_indx
    do l=1,Norb
       gwf(l,l) = xi*w + xmu - hloc(ispin,ispin,l,l) - fg_delta_hybrd(w,l,l,a) !sum(vps(l,:)**2/(xi*w-eps(:)))
       do m=l+1,Norb
          gwf(l,m) = -hloc(ispin,ispin,l,m) - fg_delta_hybrd(w,l,m,a) !sum(vps(l,:)*vps(m,:)/(xi*w-eps(:)))
          gwf(m,l) = -hloc(ispin,ispin,m,l) - fg_delta_hybrd(w,m,l,a) !sum(vps(m,:)*vps(l,:)/(xi*w-eps(:)))
       enddo
    enddo
    call matrix_inverse(gwf)
    gg = gwf(orb1,orb2)
  end function fg_weiss_hybrd















   end MODULE ED_CHI2FIT
