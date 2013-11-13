!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE ED_VARS_GLOBAL
  USE ED_BATH
  !
  USE OPTIMIZE, only:fmin_cg
  USE MATRIX,   only:matrix_inverse

  implicit none
  private


  interface chi2_fitgf
     module procedure chi2_fitgf_irred,chi2_fitgf_hybrd
  end interface chi2_fitgf

  public :: chi2_fitgf

  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb
  integer,dimension(:),allocatable      :: getIorb,getJorb
  integer                               :: Orb_indx

contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_irred(fg,bath_,ispin,iverbose)
    complex(8),dimension(:,:)            :: fg
    real(8),dimension(:,:),intent(inout) :: bath_
    integer                              :: ispin
    real(8),dimension(2*Nbath)           :: a
    integer                              :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
    real(8)                              :: chi
    logical                              :: check
    type(effective_bath)                 :: dmft_bath
    logical,optional :: iverbose
    logical          :: iverbose_
    complex(8) :: fgand
    real(8) :: w
    if(mpiID==0)then
       iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
       if(cg_scheme=='weiss')then
          write(LOGfile,"(A)")"CHI2FIT: Weiss field function:"
       else
          write(LOGfile,"(A)")"CHI2FIT: Delta function:"
       endif
       if(size(fg,1)/=Norb)stop"CHI2_FITGF: wrong dimension 1 in chi2_input"
       check= check_bath_dimension(bath_)
       if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
       Ldelta = Nfit
       if(Ldelta>size(fg,2))Ldelta=size(fg,2)
       !
       allocate(Fdelta(1,Ldelta))
       allocate(Xdelta(Ldelta))
       allocate(Wdelta(Ldelta))
       !
       forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
       !
       select case(Cg_weight)
       case default
          Wdelta=(/(1.d0,i=1,Ldelta)/)
       case(1)
          Wdelta=(/(real(i,8),i=1,Ldelta)/)
       case(2)
          Wdelta=Xdelta
       end select
       !
       do iorb=1,Norb
          Orb_indx=iorb
          Fdelta(1,1:Ldelta) = fg(iorb,1:Ldelta)
          call allocate_bath(dmft_bath)
          call set_bath(bath_,dmft_bath)
          a(1:Nbath)         = dmft_bath%e(ispin,iorb,1:Nbath) 
          a(Nbath+1:2*Nbath) = dmft_bath%v(ispin,iorb,1:Nbath)
          if(cg_scheme=='weiss')then
             call fmin_cg(a,chi2_weiss,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
          else
             call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
          endif
          write(LOGfile,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
          dmft_bath%e(ispin,iorb,1:Nbath) = a(1:Nbath)
          dmft_bath%v(ispin,iorb,1:Nbath) = a(Nbath+1:2*Nbath)
          if(iverbose_)call write_bath(dmft_bath,LOGfile)
          call copy_bath(dmft_bath,bath_)
          call deallocate_bath(dmft_bath)
       enddo
       call dump_fit_result(bath_,ispin)
       !
       print*," "
       deallocate(Fdelta,Xdelta,Wdelta)
    endif
#ifdef _MPI
    call MPI_BCAST(bath_,size(bath_),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif
  contains
    subroutine dump_fit_result(bath_,ispin)
      real(8),dimension(:,:)       :: bath_
      complex(8)                   :: fgand
      integer                      :: i,j,iorb,ispin
      real(8)                      :: w
      character(len=20)            :: suffix
      integer :: unit
      do iorb=1,Norb
         suffix="_orb"//reg(txtfy(iorb))//"_"//reg(txtfy(ispin))//".ed"
         fgand=zero
         unit=free_unit()
         open(unit,file="fit_delta"//reg(suffix))
         if(cg_scheme=='weiss')then
            do i=1,Ldelta
               w = Xdelta(i)
               fgand = xi*w + xmu - eloc(iorb) - delta_and(ispin,iorb,xi*w,bath_)
               fgand = one/fgand
               write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand),&
                    dreal(Fdelta(1,i)),dreal(fgand)
            enddo
         else
            do i=1,Ldelta
               w = Xdelta(i)
               fgand = delta_and(ispin,iorb,xi*w,bath_)
               write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand),&
                    dreal(Fdelta(1,i)),dreal(fgand)
            enddo
         endif
         close(unit)
      enddo
    end subroutine dump_fit_result
  end subroutine chi2_fitgf_irred








  !+-------------------------------------------------------------+
  !PURPOSE  : 
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
    !
    if(mpiID==0)then
       if(cg_scheme=='weiss')then
          write(LOGfile,"(A)")"CHI2FIT: Weiss field function:"
       else
          write(LOGfile,"(A)")"CHI2FIT: Delta function:"
       endif
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
       Ldelta = Nfit
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
          Wdelta=(/(1.d0,i=1,Ldelta)/)
       case(1)
          Wdelta=(/(real(i,8),i=1,Ldelta)/)
       case(2)
          Wdelta=Xdelta
       end select
       !
       a(:) = bath_(ispin,:)
       if(cg_scheme=='weiss')then
          call fmin_cg(a,chi2_hybrd_weiss,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
       else
          call fmin_cg(a,chi2_hybrd,dchi2_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
       endif
       bath_(ispin,:) = a(:)
       call dump_fit_result(bath_,ispin)
       write(*,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
       print*," "
       deallocate(Fdelta,Xdelta,Wdelta)
       deallocate(getIorb,getJorb)
    endif
#ifdef _MPI
    call MPI_BCAST(bath_,size(bath_),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif
  contains
    subroutine dump_fit_result(bath_,ispin)
      real(8),dimension(:,:) :: bath_
      complex(8)             :: fgand
      integer                :: i,j,l,iorb,jorb,ispin
      real(8)                :: w
      character(len=20)      :: suffix
      integer                :: unit
      do l=1,totNorb
         iorb=getIorb(l)
         jorb=getJorb(l)
         suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//".ed"
         fgand=zero
         unit=free_unit()
         open(unit,file="fit_delta"//reg(suffix))
         do i=1,Ldelta
            w=Xdelta(i)
            fgand = delta_and(ispin,iorb,jorb,xi*w,bath_)
            write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(l,i)),dimag(fgand),&
                 dreal(Fdelta(l,i)),dreal(fgand)
         enddo
         close(unit)
      enddo
    end subroutine dump_fit_result
  end subroutine chi2_fitgf_hybrd



  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************



  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 and the \Delta_Anderson/Hybridization 
  ! function w/  gradient d\Delta_Anderson to be minimized.  
  ! These functions are for normal (i.e. irreducible) bath.
  !+-------------------------------------------------------------+
  function chi2(a)
    real(8),dimension(:)         ::  a
    complex(8),dimension(Ldelta) ::  g0
    real(8)                      ::  chi2
    integer                      ::  i,iorb
    chi2 = 0.d0 
    iorb=Orb_indx
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = fg_anderson(xdelta(i),iorb,a)
    enddo
    chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
  end function chi2

  ! the analytic GRADIENT of \chi^2
  !+-------------------------------------------------------------+
  function dchi2(a)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2
    real(8),dimension(size(a))           :: df
    complex(8),dimension(Ldelta)         :: g0
    complex(8),dimension(Ldelta,size(a)) :: dg0
    integer                              :: i,j,iorb
    df=0.d0
    iorb=Orb_indx
    do i=1,Ldelta
       g0(i)    = fg_anderson(xdelta(i),iorb,a)
       dg0(i,:) = grad_fg_anderson(xdelta(i),iorb,a)
    enddo
    do j=1,size(a)
       df(j)=sum( dreal(Fdelta(1,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
            sum(  dimag(Fdelta(1,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
    enddo
    dchi2 = -2.d0*df
  end function dchi2

  ! the \Delta_Anderson function used in \chi^2 and d\chi^2
  ! \Delta = \sum_l V_l^2/(iw-e_l)
  !+-------------------------------------------------------------+
  function fg_anderson(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)                   :: gg,delta
    integer                      :: i,iorb,Nb
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    delta=zero
    do i=1,Nb
       delta=delta+vps(i)**2/(xi*w-eps(i))
    enddo
    gg=delta
  end function fg_anderson

  ! the gradient d\Delta_Anderson function used in d\chi^2
  ! d\Delta = \grad_{V_k,e_k}\sum_l V_l^2/(iw-e_l)
  !+-------------------------------------------------------------+
  function grad_fg_anderson(w,iorb,a) result(dgz)
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
    ! select case(cg_scheme)
    ! case default
    do i=1,Nb
       dgz(i)    = vps(i)*vps(i)/(xi*w-eps(i))**2
       dgz(i+Nb) = 2.d0*vps(i)/(xi*w-eps(i))
    enddo
    ! case ('weiss')
    !    !get the WF gg = one/(iw+mu-eloc-\Delta)
    !    delta=zero
    !    do i=1,Nb
    !       delta=delta + vps(i)**2/(xi*w-eps(i))
    !    enddo
    !    gg = (xi*w+xmu-eloc(iorb)-delta)
    !    !get the gradient \partial_a \Delta
    !    do i=1,Nb
    !       dgz(i)    = vps(i)*vps(i)/(xi*w-eps(i))**2
    !       dgz(i+Nb) = 2.d0*vps(i)/(xi*w-eps(i))
    !    enddo
    !    !build the total gradient
    !    !\partial_a G_0 = \partial_a (iw+mu-eloc-\Delta(a))^{-1}
    !    ! = \partial_a\Delta(a)*(iw+mu-eloc-\Delta(a))^{-2}
    !    dgz=dgz/gg/gg
    ! end select
  end function grad_fg_anderson



  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************




  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 and \Delta_Anderson/Hybridization  
  ! function w/ gradient d\Delta_Anderson
  ! These functions are for HYBRID bath.
  !+-------------------------------------------------------------+
  function chi2_hybrd(a) result(chi2)
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
          Delta_orb(i)= fg_anderson_hybrd(xdelta(i),iorb,jorb,a)
       enddo
       chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
    enddo
    chi2=sum(chi_orb)
  end function chi2_hybrd

  ! the analytic GRADIENT of \chi^2
  !+-------------------------------------------------------------+
  function dchi2_hybrd(a) result(dchi2)
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
          g0(i)    = fg_anderson_hybrd(xdelta(i),iorb,jorb,a)
          dg0(i,:) = grad_fg_anderson_hybrd(xdelta(i),iorb,jorb,a)
       enddo
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Fdelta(l,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
               sum( dimag(Fdelta(l,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
       enddo
    enddo
    dchi2 = -2.d0*sum(df,1)     !sum over all orbital indices
  end function dchi2_hybrd

  ! the \Delta_Anderson function used in \chi^2 and d\chi^2
  ! \Delta_ab = \sum_l V^a_l*V^b_l/(iw-e_l)
  !+-------------------------------------------------------------+
  function fg_anderson_hybrd(w,orb1,orb2,a) result(gg)
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
  end function fg_anderson_hybrd

  ! the gradient d\Delta_Anderson function used in d\chi^2
  ! d\Delta_ab = \grad_{V^c_k,e_k}\sum_l V^a_l*V^c_l/(iw-e_l)
  !+-------------------------------------------------------------+
  function grad_fg_anderson_hybrd(w,orb1,orb2,a) result(dgz)
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
  end function grad_fg_anderson_hybrd



  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************




  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 for the \calG_0 Weiss function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient.
  ! NORMAL bath.
  !+-------------------------------------------------------------+
  function chi2_weiss(a) result(chi2)
    real(8),dimension(:)         ::  a
    complex(8),dimension(Ldelta) ::  g0
    real(8)                      ::  chi2
    integer                      ::  i,iorb
    chi2 = 0.d0 
    iorb=Orb_indx
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = fg_weiss(xdelta(i),iorb,a)
    enddo
    chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
  end function chi2_weiss

  ! the inverse non-interacting GF (~inverse Weiss Field) 
  ! used in \chi^2(\caG_0 - G_0) 
  ! \G_0 = [iw_n + xmu - H_loc - \Delta]^-1
  !+-------------------------------------------------------------+
  function fg_weiss(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)                   :: gg,delta
    integer                      :: i,iorb,Nb
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    delta=zero
    do i=1,Nb
       delta=delta + vps(i)**2/(xi*w-eps(i))
    enddo
    gg = one/(xi*w+xmu-eloc(iorb)-delta)
  end function fg_weiss


  !+-------------------------------------------------------------+
  !PURPOSE: Evaluate the \chi^2 for the \calG_0 Weiss function 
  ! The Gradient is not evaluated, so the minimization requires 
  ! a numerical estimate of the gradient.
  !  HYBRID bath.
  !+-------------------------------------------------------------+
  function chi2_hybrd_weiss(a) result(chi2)
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
          Delta_orb(i)= fg_anderson_hybrd(xdelta(i),iorb,jorb,a)
       enddo
       chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
    enddo
    chi2=sum(chi_orb)
  end function chi2_hybrd_weiss

  ! the inverse non-interacting GF (~inverse Weiss Field) 
  ! used in \chi^2(\caG_0 - G_0) 
  ! {\GG_0}_ab = {[iw_n + xmu - H_loc - \Delta]^-1}_ab
  !+-------------------------------------------------------------+
  function fg_hybrd_weiss(w,orb1,orb2,a) result(gg)
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
    gwf=zero
    do l=1,Norb
       gwf(l,l)=xi*w + xmu - eloc(l)
       do m=1,Norb
          do i=1,Nbath
             gwf(l,m) = gwf(l,m) - vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))
          enddo
       enddo
    enddo
    call matrix_inverse(gwf)
    gg = gwf(orb1,orb2)
  end function fg_hybrd_weiss



end MODULE ED_CHI2FIT
