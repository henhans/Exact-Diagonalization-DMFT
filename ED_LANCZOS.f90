  include "LANCZOS_ARPACK.f90"
  include "LANCZOS_PLAIN.f90" 
  module ED_LANCZOS
    USE ED_VARS_GLOBAL
    use PLAIN_LANCZOS
    implicit none
    integer :: numzero
    integer,dimension(:),allocatable :: iszero
    type gstate
       real(8)                      :: egs
       real(8),dimension(:),pointer :: vec
    end type gstate
    type(gstate),dimension(:),allocatable :: groundstate
    real(8),dimension(:,:),allocatable :: H0

  contains
    subroutine HtimesV(N,v,Hv)
      real(8),dimension(N) :: v
      real(8),dimension(N) :: Hv
      integer              :: N
      call dgemv('N',n,n,1.d0,H0,n,v,1,0.d0,Hv,1)
    end subroutine HtimesV

  end module ED_LANCZOS
