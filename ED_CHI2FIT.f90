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
  complex(8),dimension(:),allocatable :: fdelta
  complex(8),dimension(:),allocatable :: xdelta
contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf(fg,ee,vv)
    complex(8),dimension(:)         :: fg
    real(8),dimension(size(fg))     :: wm
    real(8),dimension(:)            :: ee,vv
    real(8),dimension(2*size(ee))   :: a
    integer                         :: Nb
    integer                         :: iter
    real(8),parameter               ::   ftol=1.d-12
    real(8)                         :: chi
    !
    wm = pi/beta*dble(2*arange(1,size(fg))-1)
    Nb = size(ee)
    if(size(vv)/=Nb)then
       print*,"chi2_fitgf: size Bath error."
       stop
    endif
    a(1:nb)      = ee
    a(nb+1:2*nb) = vv

    Ldelta = size(fg)
    if(allocated(fdelta))deallocate(fdelta)
    if(allocated(xdelta))deallocate(xdelta)
    allocate(fdelta(Ldelta))
    allocate(xdelta(Ldelta))
    fdelta = fg
    xdelta = xi*wm

    call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cgNitmax,ftol=cgFtol)
    write(*,"(A,E17.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
    ee = a(1:nb)
    vv = a(nb+1:2*nb)

    call dump_fit_result(wm,fg,ee,vv)
    deallocate(fdelta,xdelta)
  end subroutine chi2_fitgf


  !********************************************************************
  !********************************************************************
  !********************************************************************




  function chi2(x)
    real(8),dimension(:)  ::  x
    real(8)               ::  chi2
    integer               ::  i
    complex(8)            ::  g0
    chi2 = 0.d0 
    do i=1,Ldelta   !Number of freq. in common to the module
       g0   = gand(xdelta(i),x)
       chi2 = chi2 + abs(fdelta(i)-g0)**2
    end do
  end function chi2

  function dchi2(x)
    real(8),dimension(:)              :: x
    real(8),dimension(size(x))       :: dchi2,df
    integer                          :: i,j
    complex(8)                       :: g0
    complex(8),dimension(size(x))    :: dg0
    df=0.d0
    do i=1,Ldelta
       g0  = gand(xdelta(i),x)
       dg0 = grad_gand(xdelta(i),x)
       do j=1,size(x)
          df(j) = df(j)+(dreal(fdelta(i))-dreal(g0))*dreal(dg0(j)) + &
               (dimag(fdelta(i))-dimag(g0))*dimag(dg0(j))
       enddo
    enddo
    dchi2 = -2.d0*df
  end function dchi2

  function gand(z,a) result(gg)
    complex(8)             :: z
    real(8),dimension(:)   :: a
    real(8),dimension(size(a)/2) :: eps,vps
    complex(8)             :: gg
    integer                :: i,Nb
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    gg=zero
    do i=1,Nb
       gg=gg+vps(i)**2/(z-eps(i))
    enddo
  end function gand

  function grad_gand(z,a) result(dgz)
    complex(8)                      :: z
    real(8),dimension(:)            :: a
    real(8),dimension(size(a)/2)    :: eps,vps
    complex(8),dimension(size(a))   :: dgz
    integer                         :: i,Nb
    dgz=zero
    Nb=size(a)/2
    eps=a(1:Nb)
    vps=a(Nb+1:2*Nb)
    do i=1,Nb
       dgz(i)    = vps(i)*vps(i)/(z-eps(i))**2
       dgz(i+Nb) = 2.d0*vps(i)/(z-eps(i))
    enddo
  end function grad_gand



  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine dump_fit_result(wm,fg,epsi,vi)
    complex(8),dimension(:)        :: fg
    real(8),dimension(size(fg))    :: wm
    complex(8),dimension(size(fg)) :: fgand
    real(8),dimension(:)           :: epsi,vi
    integer                        :: i,j,Lw
    Lw=size(fg)
    fgand=zero
    do i=1,Lw
       fgand(i) = delta_and(xi*wm(i),epsi,vi)
    enddo
    call splot("fit_delta.ed",wm,dimag(fg),dimag(fgand),dreal(fg),dreal(fgand))
  end subroutine dump_fit_result



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_CHI2FIT
