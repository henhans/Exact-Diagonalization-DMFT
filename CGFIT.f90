  include "cgfit_common.f90"
  include "cgfit_function.f90"
  include "cgfit_minimize.f90"
  module CGFIT
    USE CGFIT_COMMON
    USE CGFIT_FUNCTION
    implicit none
    private
    public  :: fitgreen

  contains

    !+-------------------------------------------------------------------+
    !PURPOSE : fit given Green function with using CG minimization.
    ! This subroutine fits a given GF "g" at frequencies z=iw_n
    ! against the functional form Sum_i=1^nsite b(i)^2/(z-a(i))
    ! using the conjugate gradient method.
    !+-------------------------------------------------------------------+
    subroutine fitgreen(z,g,aa,bb)
      complex(8),dimension(:)       :: g
      real(8),dimension(size(g))    :: z
      real(8),dimension(:)          :: aa,bb
      integer                       :: i,iter,nsite
      real(8),dimension(2*size(aa)) :: a
      real(8),parameter             :: ftol=1.d-12
      real(8)                       :: chi
      write(*,*)"Fit Delta(iw):"
      Nff = size(g)
      allocate(ome(Nff),gf(Nff))
      ome(1:Nff) = z(1:Nff)
      gf(1:Nff)  = g(1:Nff)
      nsite=size(aa)
      forall(i=1:nsite)
         a(i)       = aa(i)
         a(i+nsite) = bb(i)
      end forall
      call conjugate(a,ftol,iter,chi)
      forall(i=1:nsite)
         aa(i) = a(i)
         bb(i) = a(i+nsite)
      end forall
      write(*,"(A,E17.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
      deallocate(ome,gf)
    end subroutine fitgreen



    !********************************************************************
    !********************************************************************
    !********************************************************************






    !+-------------------------------------------------------------------+
    !     PROGRAM  : CONJUGATE
    !     TYPE     : subroutine
    !     PURPOSE  : Minimize the Chi^2 distance using conjugate gradient
    !     Adapted by FRPRM subroutine from NumRec (10.6)
    !     Given a starting point P that is a vector of length N, 
    !     the Fletcher-Reeves-Polak-Ribiere minimisation is performed 
    !     n a functin FUNC,using its gradient as calculated by a 
    !     routine DFUNC. The convergence tolerance on the function 
    !     value is input as FTOL.  
    !     Returned quantities are: 
    !     - P (the location of the minimum), 
    !     - ITER (the number of iterations that were performed), 
    !     - FRET (the minimum value of the function). 
    !     The routine LINMIN is called to perform line minimisations.
    !     Minimisation routines: DFPMIN, LINMIN, MNBRAK, BRENT and F1DIM
    !     come from Numerical Recipes.
    !+-------------------------------------------------------------------+
    subroutine conjugate(p,ftol,iter,fret)
      IMPLICIT NONE
      REAL(8), DIMENSION(:), INTENT(INOUT) :: p
      INTEGER, INTENT(OUT)                 :: iter
      REAL(8), INTENT(IN)                  :: ftol
      REAL(8), INTENT(OUT)                 :: fret
      INTEGER, PARAMETER                   :: ITMAX=500
      REAL(8), PARAMETER                   :: EPS=1.D-9
      INTEGER                              :: its
      REAL(8)                              :: dgg,fp,gam,gg
      REAL(8), DIMENSION(size(p))          :: g,h,xi
      fp=func(p)
      xi=dfunc(p) !  call dfunc(p, xi)
      g=-xi
      h=g
      xi=h
      do its=1,ITMAX
         iter=its
         call linmin(p,xi,fret)
         if (2.0*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
         !fp=fret
         fp = func(p) !========MODIFICATION=======
         xi = dfunc(p)        
         gg=dot_product(g,g)
         !dgg=dot_product(xi,xi)   !Fletcher-Reeves.
         dgg=dot_product(xi+g,xi)  !Polak-Ribiere
         if (gg == 0.0) RETURN
         gam=dgg/gg
         g=-xi
         h=g+gam*h
         xi=h
      end do
      print*, 'Too many iteration in CG'
      return
    END SUBROUTINE conjugate
    !********************************************************************
    !********************************************************************
    !********************************************************************







    !+-------------------------------------------------------------------+
    !     PROGRAM  : LINMIN
    !     TYPE     : subroutine
    !     PURPOSE  : Minimization routine
    !     COMMENTS : 
    !     Given an N dimensional point P and an N dimensional direction 
    !     XI, LINMIN
    !     moves and resets P to where the function FUNC(P) takes on a minimum
    !     along the direction XI from P, and replaces XI by the actual vector
    !     displacement that P was moved.  Also returns FRET the value of 
    !     FUNC at the returned location P.  
    !     This is actually all accomplished by calling the routines 
    !     MNBRAK and BRENT.
    !+-------------------------------------------------------------------+
    SUBROUTINE linmin(p,xi,fret)
      USE F1DIM_MOD
      USE LOCAL_UTILS
      REAL(8), INTENT(OUT) :: fret
      REAL(8), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
      REAL(8), PARAMETER :: TOL=1.0e-4
      REAL(8) :: ax,bx,fa,fb,fx,xmin,xx
      ncom=size(p) ; if(ncom /= size(xi))stop "Error in LinMin"
      pcom=>p
      xicom=>xi
      ax=0.0
      xx=1.0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      !...construct the vector results to return
      xi=xmin*xi
      p=p+xi
      return
    end subroutine linmin
    !********************************************************************
    !********************************************************************
    !********************************************************************
  end module CGFIT











