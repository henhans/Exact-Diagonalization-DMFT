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
     module procedure chi2_fitgf_irr,chi2_fitgf_red
  end interface chi2_fitgf

  public :: chi2_fitgf

  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb
  integer,dimension(:),allocatable      :: getIorb,getJorb


contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_irr(fg,bath_,ispin)
    complex(8),dimension(:,:)          :: fg
    real(8),dimension(:),intent(inout) :: bath_
    integer                            :: ispin
    real(8),dimension(2*Nbath)         :: a
    integer                            :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
    real(8)                            :: chi
    if(mpiID==0)then
       write(LOGfile,"(A)")"FIT Delta function:"
       if(size(fg,1)/=Norb)stop"CHI2_FITGF: wrong dimension 1 in Delta_input"
       call check_bath_dimension(bath_)
       Ldelta = size(fg,2)
       allocate(Fdelta(1,Ldelta))
       allocate(Xdelta(Ldelta))
       allocate(Wdelta(Ldelta))
       forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
       select case(CG_type)
       case(0)
          Wdelta=(/(1.d0,i=1,Ldelta)/)
       case(1)
          Wdelta=(/(real(i,8),i=1,Ldelta)/)
       case(2)
          Wdelta=Xdelta
       end select
       do iorb=1,Norb
          Fdelta(1,:) = fg(iorb,:)
          stride_spin = (ispin-1)*Norb*Nbath
          stride_orb  = (iorb-1)*2*Nbath
          ifirst = stride_spin + stride_orb + 1
          ilast  = stride_spin + stride_orb + Nbath + Nbath
          a(:) = bath_(ifirst:ilast)
          call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
          bath_(ifirst:ilast) = a(:)
          call dump_fit_result(a,ispin,iorb)
          write(LOGfile,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
       enddo
       !
       print*," "
       deallocate(Fdelta,Xdelta,Wdelta)
    endif
#ifdef _MPI
    call MPI_BCAST(bath_,size(bath_),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif
  contains
    subroutine dump_fit_result(bath_,ispin,iorb)
      real(8),dimension(:)         :: bath_
      real(8),dimension(size(bath_)/2)    :: eps,vps
      complex(8),dimension(Ldelta) :: fgand
      integer                      :: i,j,iorb,ispin,Nb
      real(8)                      :: w
      character(len=20)            :: suffix
      integer :: unit
      suffix="_orb"//reg(txtfy(iorb))//"_"//reg(txtfy(ispin))//".ed"
      Nb=size(bath_)/2
      eps=bath_(1:Nb)
      vps=bath_(Nb+1:2*Nb)
      fgand=zero
      unit=free_unit()
      open(unit,file="fit_delta"//reg(suffix))
      do i=1,Ldelta
         w=Xdelta(i)
         do j=1,Nb
            fgand(i)=fgand(i)+vps(j)**2/(xi*w-eps(j))
         enddo
         write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(1,i)),dimag(fgand(i)),&
              dreal(Fdelta(1,i)),dreal(fgand(i))
      enddo
      close(unit)
    end subroutine dump_fit_result
  end subroutine chi2_fitgf_irr





  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_red(fg,bath_,ispin)
    complex(8),dimension(:,:,:)        :: fg
    real(8),dimension(:),intent(inout) :: bath_
    integer                            :: ispin
    real(8),dimension((1+Norb)*Nbath)  :: a
    integer                            :: iter,stride,ifirst,ilast,i,j,corb
    integer :: iorb,jorb
    real(8)                            :: chi
    !
    if(mpiID==0)then
       write(LOGfile,"(A)")"CHI2FIT: Delta function:"
       if(size(fg,1)/=Norb)stop "CHI2FIT: wrong dimension 1 in Delta_input"
       if(size(fg,2)/=Norb)stop "CHI2FIT: wrong dimension 2 in Delta_input"
       call check_bath_dimension(bath_)
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
       Ldelta = size(fg,3)    
       allocate(Fdelta(totNorb,Ldelta))
       allocate(Xdelta(Ldelta))
       allocate(Wdelta(Ldelta))
       do i=1,totNorb
          Fdelta(i,:) = fg(getIorb(i),getJorb(i),:)
       enddo
       forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
       select case(CG_type)
       case(0)
          Wdelta=(/(1.d0,i=1,Ldelta)/)
       case(1)
          Wdelta=(/(real(i,8),i=1,Ldelta)/)
       case(2)
          Wdelta=Xdelta
       end select
       !
       stride=(ispin-1)*(Norb+1)*Nbath
       ifirst=stride+1 ; ilast =stride+(Norb+1)*Nbath
       a(:) = bath_(ifirst:ilast)
       print*,bath_
       stop "ED_CHI2FIT"
       call fmin_cg(a,chi2_red,dchi2_red,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
       bath_(ifirst:ilast) = a(:)
       !
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
      real(8),dimension(:)             :: bath_
      complex(8),dimension(Ldelta)     :: fgand
      integer                          :: i,j,l,iorb,jorb,ispin
      real(8)                          :: w
      character(len=20)                :: suffix
      integer                          :: unit
      do l=1,totNorb
         iorb=getIorb(l)
         jorb=getJorb(l)
         suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//".ed"
         fgand=zero
         unit=free_unit()
         open(unit,file="fit_delta"//reg(suffix))
         do i=1,Ldelta
            w=Xdelta(i)
            fgand(i) = delta_and(ispin,iorb,jorb,xi*w,bath_)
            write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(l,i)),dimag(fgand(i)),&
                 dreal(Fdelta(l,i)),dreal(fgand(i))
         enddo
         close(unit)
      enddo
    end subroutine dump_fit_result
  end subroutine chi2_fitgf_red





  !+-------------------------------------------------------------+
  !PURPOSE  : the CHI^2 to be minimized
  !+-------------------------------------------------------------+
  function chi2(a)
    real(8),dimension(:)         ::  a
    complex(8),dimension(Ldelta) ::  g0
    real(8)                      ::  chi2
    integer                      ::  i
    chi2 = 0.d0 
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = gand(xdelta(i),a)
    enddo
    chi2=sum(abs(Fdelta(1,:)-g0(:))**2/Wdelta(:))
  end function chi2
  !
  function chi2_red(a) result(chi2)
    real(8),dimension(:)                 ::  a
    complex(8),dimension(totNorb,Ldelta) ::  g0
    real(8),dimension(totNorb)           ::  chi_sum
    real(8)                              ::  chi2
    integer                              ::  i,l,iorb,jorb
    chi_sum = 0.d0 
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       print*,l,iorb,jorb
       do i=1,Ldelta
          g0(l,i) = gand_red(xdelta(i),iorb,jorb,a)
       enddo
       chi_sum(l) = sum(abs(Fdelta(l,:)-g0(l,:))**2/Wdelta(:))
       print*,l,chi_sum(l)
    enddo
    chi2=sum(chi_sum)
  end function chi2_red



  !+-------------------------------------------------------------+
  !PURPOSE: the GRADIENT OF CHI^2 used in the CG minimization
  !+-------------------------------------------------------------+
  function dchi2(a)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2
    real(8),dimension(size(a))           :: df
    complex(8),dimension(Ldelta)         :: g0
    complex(8),dimension(Ldelta,size(a)) :: dg0
    integer                              :: i,j
    df=0.d0
    do i=1,Ldelta
       g0(i)    = gand(xdelta(i),a)
       dg0(i,:) = grad_gand(xdelta(i),a)
    enddo
    do j=1,size(a)
       df(j)=sum( dreal(Fdelta(1,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
            sum(  dimag(Fdelta(1,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
    enddo
    dchi2 = -2.d0*df
  end function dchi2
  !
  function dchi2_red(a) result(dchi2)
    real(8),dimension(:)                         :: a
    real(8),dimension(size(a))                   :: dchi2
    real(8),dimension(totNorb,size(a))           :: df
    complex(8),dimension(totNorb,Ldelta)         :: g0
    complex(8),dimension(totNorb,Ldelta,size(a)) :: dg0
    integer                                      :: i,j,l,iorb,jorb
    df=0.d0
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       do i=1,Ldelta
          g0(l,i)    = gand_red(xdelta(i),iorb,jorb,a)
          dg0(l,i,:) = grad_gand_red(xdelta(i),iorb,jorb,a)
       enddo
       !
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Fdelta(l,:)-g0(l,:))*dreal(dg0(l,:,j))/Wdelta(:) ) + &
               sum( dimag(Fdelta(l,:)-g0(l,:))*dimag(dg0(l,:,j))/Wdelta(:) )
       enddo
    enddo
    dchi2 = -2.d0*sum(df,1)     !sum over all orbital indices
  end function dchi2_red






  !+-------------------------------------------------------------+
  !PURPOSE: evaluate the Delta function appearing in the CHI2
  ! formulae. This routine is similar to that used in ED_BATH
  ! but 
  !+-------------------------------------------------------------+
  function gand(w,a) result(gg)
    real(8)                      :: w
    real(8),dimension(:)         :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)                   :: gg
    integer                      :: i,Nb
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    gg=zero
    do i=1,Nb
       gg=gg+vps(i)**2/(xi*w-eps(i))
    enddo
  end function gand
  !
  function gand_red(w,orb1,orb2,a) result(gg)
    real(8)                      :: w
    integer                      :: orb1,orb2
    real(8),dimension(:)         :: a
    real(8),dimension(Nbath)     :: eps
    real(8),dimension(Norb,Nbath):: vps
    complex(8)                   :: gg
    integer                      :: i,l
    eps=a(1:Nbath)
    print*,"size(a)",size(a)
    print*,1,Nbath
    print*,"eps=",eps
    do l=1,Norb
       vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
       print*,l*Nbath+1,(l+1)*Nbath
       print*,size(vps(l,:)),size(a((l*Nbath+1):((l+1)*Nbath)))
       print*,"vps=",vps(l,:)
    enddo
    gg=zero
    do i=1,Nbath
       gg=gg + vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))
    enddo
  end function gand_red




  !+-------------------------------------------------------------+
  !PURPOSE : Evaluate the gradient of the Delta function 
  ! as appearing in the dCHI2 gradient function.
  !+-------------------------------------------------------------+
  function grad_gand(w,a) result(dgz)
    real(8)                         :: w
    real(8),dimension(:)            :: a
    real(8),dimension(size(a)/2)    :: eps,vps
    complex(8),dimension(size(a))   :: dgz
    integer                         :: i,Nb
    dgz=zero
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    do i=1,Nb
       dgz(i)    = vps(i)*vps(i)/(xi*w-eps(i))**2
       dgz(i+Nb) = 2.d0*vps(i)/(xi*w-eps(i))
    enddo
  end function grad_gand
  !
  function grad_gand_red(w,orb1,orb2,a) result(dgz)
    real(8)                         :: w
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
       do l=1,Norb
          dgz(l*Nbath+i) = (&
               dKronecker(l,orb1)*Vps(orb2,i) + &
               dKronecker(l,orb2)*Vps(orb1,i)   &
               )/(xi*w-eps(i))
       enddo
    enddo
  contains
    function dKronecker(i,j) result(delta)
      real(8) :: delta
      integer :: i,j
      delta=0.d0
      if(i==j)delta=1.d0
    end function dKronecker
  end function grad_gand_red



end MODULE ED_CHI2FIT
