! File containing the routines to build the functions
! to be fitted:
! - - - - - - - - - - - - - - - -
! func  = the function of minimize (Chi^2)
! G0and = Anderson hybridization function (\Delta)
! - - - - - - - - - - - - - - - -
! dfunc     = the gradient of the function of minimize (\Grad\Chi^2)
! gradG0and = Gradient of the Anderson Hybridization function (\Grad\Delta)
! - - - - - - - - - - - - - - - -
!
! These functions have to be provided by the user
! Here the functions corresponds to the Chi^2 fit 
! for the CG minimization of the Anderson model
! non-interacting Greens function
! \Delta(w) = \sum_k=1,Ns V^2_k/(w-\epsk_k)
module CGFIT_FUNCTION
  USE CGFIT_COMMON
  implicit none
  private 
  public :: func, dfunc
contains
  !+-------------------------------------------------------------------+
  !     PROGRAM  : FUNC
  !     TYPE     : subroutine
  !     PURPOSE  : Build the Chi^2 distance function
  !     COMMENTS : 
  !+-------------------------------------------------------------------+
  function func(a)
    real(8),dimension(:),intent(in)  ::  a
    real(8)                          ::  func
    integer                          :: i
    complex(8)                       :: Gand
    func = 0.d0 
    do i=1,Nff   !Number of freq. in common to the module
       if (abs(Ome(i)) > abs(WMIN)) then
          Gand = G0and(Ome(i), a)
          func = func + abs(Gf(i)-Gand)**2 !Gf in common
       endif
    enddo
    return
  end function func
  !********************************************************************
  !********************************************************************
  !********************************************************************


  !+-------------------------------------------------------------------+
  !     PROGRAM  : DFUNC
  !     TYPE     : subroutine
  !     PURPOSE  : obtain the gradient of Chi^2 function
  !     COMMENTS : 
  !+-------------------------------------------------------------------+
  function dfunc(a)
    real(8),dimension(:),intent(in)  :: a
    real(8),dimension(size(a))       :: dfunc,df
    integer                          :: i,j
    complex(8)                       :: Gand
    complex(8),dimension(size(a))    :: dGand
    df=0.d0
    do i=1,Nff
       if (abs(Ome(i)) > abs(WMIN)) then
          Gand  = G0and(Ome(i), a)
          dGand = GradG0and(Ome(i), a)
          do j=1,size(a)
             df(j) = df(j) + (real(Gf(i))-real(Gand))*real(dGand(j)) +&
                  (aimag(Gf(i))-aimag(Gand))*aimag(dGand(j))
          enddo
       endif
    enddo
    dfunc = -2.d0*df
    return
  end function dfunc
  !********************************************************************
  !********************************************************************
  !********************************************************************




  !+-------------------------------------------------------------------+
  !     PROGRAM  : G0and
  !     TYPE     : subroutine
  !     PURPOSE  : Construct the Anderson function to fit against given GF
  !     COMMENTS : 
  !+-------------------------------------------------------------------+
  function G0and(x,a)
    real(8),dimension(:),intent(in) :: a
    real(8)                         :: x
    complex(8)                      :: G0and,Gz
    integer                         :: i,Ns
    complex(8)                      :: z
    Ns=size(a)/2
    z = cmplx(0.d0,x)
    Gz= (0.d0,0.d0)
    do i=1,Ns
       Gz = Gz + a(i+ns)**2/(z-a(i))
    enddo
    G0and=Gz
    return
  end function G0and
  !********************************************************************
  !********************************************************************
  !********************************************************************



  !+-------------------------------------------------------------------+
  !     PROGRAM  : GRADG0AND
  !     TYPE     : subroutine
  !     PURPOSE  : get the gradient of the Anderson gf
  !     COMMENTS : 
  !+-------------------------------------------------------------------+
  function gradG0and(x, a)
    real(8),dimension(:),intent(in) :: a
    real(8)                         :: x
    complex(8),dimension(size(a))   :: gradG0and,dGz
    integer                         :: i,Ns
    complex(8)                      :: z  
    z = cmplx(0.d0,x)
    Ns=size(a)/2
    do i=1,ns
       dGz(i)    = a(i+ns)*a(i+ns)/(z-a(i))**2
       dGz(i+ns) =      2.d0*a(i+ns)/(z-a(i))
    enddo
    gradG0and = dGz
    return
  end function gradG0and
  !********************************************************************
  !********************************************************************
  !********************************************************************
end module CGFIT_FUNCTION
