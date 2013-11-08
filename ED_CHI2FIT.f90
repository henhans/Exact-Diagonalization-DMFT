!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE OPTIMIZE
  USE ED_VARS_GLOBAL
  USE ED_BATH
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
  subroutine chi2_fitgf_irred(fg,bath_,ispin)
    complex(8),dimension(:,:)            :: fg
    real(8),dimension(:,:),intent(inout) :: bath_
    integer                              :: ispin
    real(8),dimension(2*Nbath)           :: a
    integer                              :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
    real(8)                              :: chi
    logical                              :: check
    type(effective_bath)                 :: dmft_bath
    if(mpiID==0)then
       write(LOGfile,"(A)")"FIT Delta function:"
       if(size(fg,1)/=Norb)stop"CHI2_FITGF: wrong dimension 1 in Delta_input"
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
          call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
          dmft_bath%e(ispin,iorb,1:Nbath) = a(1:Nbath)
          dmft_bath%v(ispin,iorb,1:Nbath) = a(Nbath+1:2*Nbath)
          call copy_bath(dmft_bath,bath_)
          call deallocate_bath(dmft_bath)
          write(LOGfile,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
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
         do i=1,Ldelta
            w = Xdelta(i)
            fgand = delta_and(ispin,iorb,xi*w,bath_)
            write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand),&
                 dreal(Fdelta(1,i)),dreal(fgand)
         enddo
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
       write(LOGfile,"(A)")"CHI2FIT: Delta function:"
       if(size(fg,1)/=Norb)stop "CHI2FIT: wrong dimension 1 in Delta_input"
       if(size(fg,2)/=Norb)stop "CHI2FIT: wrong dimension 2 in Delta_input"
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
       call fmin_cg(a,chi2_hybrd,dchi2_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
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





  !+-------------------------------------------------------------+
  !PURPOSE  : the CHI^2 to be minimized
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



  !+-------------------------------------------------------------+
  !PURPOSE: the GRADIENT OF CHI^2 used in the CG minimization
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
       ! df(j) = sum(dreal( (Fdelta(1,:)-g0(:))*dg0(:,j)/Wdelta(:) ))
    enddo
    dchi2 = -2.d0*df
  end function dchi2

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





  !+-------------------------------------------------------------+
  !PURPOSE: evaluate the Delta function appearing in the CHI2
  ! formulae. This routine is similar to that used in ED_BATH
  ! but 
  !+-------------------------------------------------------------+
  function fg_anderson(w,iorb,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)                   :: gg
    integer                      :: i,iorb,Nb
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    gg=zero
    do i=1,Nb
       gg=gg+vps(i)**2/(xi*w-eps(i))
    enddo
    if(cg_scheme=='weiss')gg = xi*w + xmu - eloc(iorb) - gg
  end function fg_anderson
  !
  function fg_anderson_hybrd(w,orb1,orb2,a) result(gg)
    real(8)                      :: w
    integer                      :: orb1,orb2
    real(8),dimension(:)         :: a
    real(8),dimension(Nbath)     :: eps
    real(8),dimension(Norb,Nbath):: vps
    complex(8)                   :: gg
    integer                      :: i,l
    eps=a(1:Nbath)
    do l=1,Norb
       vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
    enddo
    gg=zero
    do i=1,Nbath
       gg=gg + vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))
    enddo
    if(cg_scheme=='weiss')then
       if(orb1==orb2)then
          gg = xi*w + xmu - eloc(Orb_indx) - gg
       else
          gg = - gg
       endif
    endif
  end function fg_anderson_hybrd




  !+-------------------------------------------------------------+
  !PURPOSE : Evaluate the gradient of the Delta function 
  ! as appearing in the dCHI2 gradient function.
  !+-------------------------------------------------------------+
  function grad_fg_anderson(w,iorb,a) result(dgz)
    real(8)                         :: w,sgn
    real(8),dimension(:)            :: a
    real(8),dimension(size(a)/2)    :: eps,vps
    complex(8),dimension(size(a))   :: dgz
    integer                         :: i,iorb,Nb
    dgz=zero
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    sgn=1.d0
    if(cg_scheme=='weiss')sgn=-1.d0
    do i=1,Nb
       dgz(i)    = sgn*vps(i)*vps(i)/(xi*w-eps(i))**2
       dgz(i+Nb) = sgn*2.d0*vps(i)/(xi*w-eps(i))
    enddo
  end function grad_fg_anderson
  !
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
    sgn=1.d0
    if(cg_scheme=='weiss')sgn=-1.d0
    do l=1,Norb
       vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
    enddo
    !
    do i=1,Nbath
       dgz(i)    = sgn*vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))**2
       if(orb1==orb2)then
          dgz(orb1*Nbath+i) = sgn*2.d0*Vps(orb1,i)/(xi*w-eps(i))
       else
          dgz(orb1*Nbath+i) = sgn*Vps(orb2,i)/(xi*w-eps(i))
          dgz(orb2*Nbath+i) = sgn*Vps(orb1,i)/(xi*w-eps(i))
       endif
    enddo
  end function grad_fg_anderson_hybrd



end MODULE ED_CHI2FIT
