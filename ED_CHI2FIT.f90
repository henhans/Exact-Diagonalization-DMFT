!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE OPTIMIZE
  USE ED_VARS_GLOBAL
  USE ED_BATH
  implicit none
  private

  public :: chi2_fitgf

  integer                             :: Ldelta
  complex(8),dimension(:),allocatable :: Fdelta
  real(8),dimension(:),allocatable    :: Xdelta,Wdelta
contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf(fg,bath,ispin)
    complex(8),dimension(:,:)          :: fg
    real(8),dimension(:),intent(inout) :: bath
    integer                            :: ispin
    real(8),dimension(2*Nbath)         :: a
    integer                            :: iter,stride_spin,stride_orb,ifirst,ilast,i,j,iorb
    real(8)                            :: chi
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"FIT Delta function:"
       if(size(fg,1)/=Norb)stop"CHI2_FITGF: wrong dimension 1 in Delta_input"
       call check_bath_dimension(bath)

       Ldelta = size(fg,2)
       allocate(Fdelta(Ldelta))
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
          Fdelta(:) = fg(iorb,:)

          stride_spin = (ispin-1)*Norb*Nbath
          stride_orb  = (iorb-1)*2*Nbath

          ifirst = stride_spin + stride_orb + 1
          ilast  = stride_spin + stride_orb + Nbath + Nbath
          a(:) = bath(ifirst:ilast)
          call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
          bath(ifirst:ilast) = a(:)
          call dump_fit_result(a,ispin,iorb)
          write(LOGfile,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
       enddo
       !
       print*," "
       deallocate(Fdelta,Xdelta,Wdelta)
#ifdef _MPI
    endif
    call MPI_BCAST(bath,size(bath),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif
  end subroutine chi2_fitgf




  !FUNCTION CHI^2 to be minimized:
  function chi2(a)
    real(8),dimension(:) ::  a
    real(8)              ::  chi2
    integer              ::  i
    complex(8)           ::  g0(Ldelta)    
    chi2 = 0.d0 
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = gand(xdelta(i),a)
    enddo
    chi2=sum(abs(fdelta(:)-g0(:))**2/Wdelta(:))
  end function chi2

  !GRADIENT OF CHI^2 used in the CG minimization:
  function dchi2(a)
    real(8),dimension(:)                 :: a
    real(8),dimension(size(a))           :: dchi2,df
    integer                              :: i,j
    complex(8)                           :: g0(Ldelta)
    complex(8),dimension(Ldelta,size(a)) :: dg0
    df=0.d0
    do i=1,Ldelta
       g0(i)    = gand(xdelta(i),a)
       dg0(i,:) = grad_gand(xdelta(i),a)
    enddo
    do j=1,size(a)
       df(j)=sum( dreal(fdelta(:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
            sum(  dimag(fdelta(:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
    enddo
    dchi2 = -2.d0*df
  end function dchi2

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



  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine dump_fit_result(bath,ispin,iorb)
    real(8),dimension(:)         :: bath
    real(8),dimension(size(bath)/2)    :: eps,vps
    complex(8),dimension(Ldelta) :: fgand
    integer                      :: i,j,iorb,ispin,Nb
    real(8)                      :: w
    character(len=20)            :: suffix
    integer :: unit
    suffix="_orb"//reg(txtfy(iorb))//"_"//reg(txtfy(ispin))//".ed"
    Nb=size(bath)/2
    eps=bath(1:Nb)
    vps=bath(Nb+1:2*Nb)
    fgand=zero
    unit=free_unit()
    open(unit,file="fit_delta"//reg(suffix))
    do i=1,Ldelta
       w=Xdelta(i)
       do j=1,Nb
          fgand(i)=fgand(i)+vps(j)**2/(xi*w-eps(j))
       enddo
       write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(i)),dimag(fgand(i)),dreal(Fdelta(i)),dreal(fgand(i))
    enddo
    close(unit)
  end subroutine dump_fit_result

end MODULE ED_CHI2FIT
