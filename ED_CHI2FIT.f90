!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, ONLY:delta_and
  implicit none
  private

  public :: chi2_fitgf

  integer                             :: Ldelta
  complex(8),dimension(:),allocatable :: Fdelta
  real(8),dimension(:),allocatable    :: Xdelta
contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf(fg,ichan)
    complex(8),dimension(:)    :: fg
    integer                    :: ichan
    real(8),dimension(2*Nbath) :: a
    integer                    :: iter
    real(8)                    :: chi
    !
    call msg("FIT Delta function:")

    Ldelta = size(fg)
    allocate(fdelta(Ldelta))
    allocate(xdelta(Ldelta))
    Fdelta = fg
    Xdelta = pi/beta*dble(2*arange(1,Ldelta)-1)

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
    deallocate(Fdelta,Xdelta)
  end subroutine chi2_fitgf


  !********************************************************************
  !********************************************************************
  !********************************************************************




  function chi2(x)
    real(8),dimension(:)  ::  x
    real(8)               ::  chi2
    integer               ::  i
    complex(8)            ::  g0(Ldelta)
    chi2 = 0.d0 
    do i=1,Ldelta   !Number of freq. in common to the module
       g0(i)   = gand(xdelta(i),x)
    enddo
    select case(CGtype)
    case(0)
       do i=1,Ldelta   !Number of freq. in common to the module
          chi2 = chi2 + abs(fdelta(i)-g0(i))**2
       end do
    case(1)
       do i=1,Ldelta   !Number of freq. in common to the module
          chi2 = chi2 + abs(fdelta(i)-g0(i))**2/dble(i)
       end do
    case(2)
       do i=1,Ldelta   !Number of freq. in common to the module
          chi2 = chi2 + abs(fdelta(i)-g0(i))**2/dble(xdelta(i))
       end do

    end select
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
    select case(CGtype)
    case(0)
       do j=1,size(x)
          do i=1,Ldelta
             df(j) = df(j)+(dreal(fdelta(i))-dreal(g0(i)))*dreal(dg0(i,j)) + &
                  (dimag(fdelta(i))-dimag(g0(i)))*dimag(dg0(i,j))
          enddo
       enddo
    case(1)
       do j=1,size(x)
          do i=1,Ldelta
             df(j) = df(j)+(dreal(fdelta(i))-dreal(g0(i)))*dreal(dg0(i,j)) + &
                  (dimag(fdelta(i))-dimag(g0(i)))*dimag(dg0(i,j))/dble(i)
          enddo
       enddo
    case(2)
       do j=1,size(x)
          do i=1,Ldelta
             df(j) = df(j)+(dreal(fdelta(i))-dreal(g0(i)))*dreal(dg0(i,j)) + &
                  (dimag(fdelta(i))-dimag(g0(i)))*dimag(dg0(i,j))/dble(xdelta(i))
          enddo
       enddo
    end select
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
       fgand(i) = delta_and(xi*w,ichan)
    enddo
    call splot("fit_delta.ed",Xdelta,dimag(Fdelta),dimag(fgand),dreal(Fdelta),dreal(fgand))
  end subroutine dump_fit_result



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_CHI2FIT
