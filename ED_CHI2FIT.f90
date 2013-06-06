!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
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
  subroutine chi2_fitgf(fg,bath,ichan)
    complex(8),dimension(:)            :: fg
    real(8),dimension(:),intent(inout) :: bath
    integer                            :: i,ichan
    real(8),dimension(2*Nbath)         :: a
    integer                            :: iter
    real(8)                            :: chi
    !
    call msg("FIT Delta function:",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)

    Ldelta = size(fg)
    allocate(Fdelta(Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    Fdelta = fg
    Xdelta = pi/beta*dble(2*arange(1,Ldelta)-1)
    select case(CGtype)
    case(0)
       Wdelta=(/(1.d0,i=1,Ldelta)/)
    case(1)
       Wdelta=(/(dble(i),i=1,Ldelta)/)
    case(2)
       Wdelta=xdelta
    end select

    a(1:Nbath)           = ebath(ichan,:)
    a(Nbath+1:2*Nbath)   = vbath(ichan,:)
    !
    call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cgNitmax,ftol=cgFtol)
    !
    ebath(ichan,:) = a(1:Nbath)
    vbath(ichan,:) = a(Nbath+1:2*Nbath)

    write(*,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
    print*," "
    call dump_fit_result(ichan)
    bath = copy_bath()
    call deallocate_bath
    deallocate(Fdelta,Xdelta,Wdelta)
  end subroutine chi2_fitgf


  !********************************************************************
  !********************************************************************
  !********************************************************************




  function chi2(x)
    real(8),dimension(:) ::  x
    real(8)              ::  chi2
    integer              ::  i
    complex(8)           ::  g0(Ldelta)    
    chi2 = 0.d0 
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = gand(xdelta(i),x)
    enddo
    chi2=sum(abs(fdelta(:)-g0(:))**2/Wdelta(:))
  end function chi2

  function dchi2(x)
    real(8),dimension(:)                 :: x
    real(8),dimension(size(x))           :: dchi2,df
    integer                              :: i,j
    complex(8)                           :: g0(Ldelta)
    complex(8),dimension(Ldelta,size(x)) :: dg0
    df=0.d0
    do i=1,Ldelta
       g0(i)    = gand(xdelta(i),x)
       dg0(i,:) = grad_gand(xdelta(i),x)
    enddo
    do j=1,size(x)
       df(j)=sum( dreal(fdelta(:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
            sum(  dimag(fdelta(:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
    enddo
    dchi2 = -2.d0*df
  end function dchi2

  function gand(w,a) result(gg)
    real(8)                :: w
    real(8),dimension(:)   :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)             :: gg
    integer                :: i,Nb
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
  subroutine dump_fit_result(ichan)
    complex(8),dimension(Ldelta) :: fgand
    integer                      :: ichan,i,j,Lw
    real(8)                      :: w
    fgand=zero
    do i=1,Ldelta
       w=Xdelta(i)
       fgand(i) = delta_bath(xi*w,ichan)
    enddo
    call splot("fit_delta.ed",Xdelta,dimag(Fdelta),dimag(fgand),dreal(Fdelta),dreal(fgand))
  end subroutine dump_fit_result



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_CHI2FIT
