module ARPACK_LANCZOS
  USE CONSTANTS, only:zero
  USE ED_VARS_GLOBAL
  implicit none
  private 

  interface lanczos_arpack
     module procedure lanczos_arpack_d,lanczos_arpack_c
  end interface lanczos_arpack
  public :: lanczos_arpack

#ifdef _MPI
  interface lanczos_parpack
     module procedure lanczos_parpack_d,lanczos_parpack_c
  end interface lanczos_parpack
  public :: lanczos_parpack
#endif

contains

  !+-------------------------------------------------------------------+
  ! This routine use PARALLEL_/ARPACK to find a few eigenvalues
  ! LAMBDA and corresponding eigenvectors X for the standard
  ! eigenvalue problem:
  !      A * X = LAMBDA * X
  ! where A is an N by N 
  ! real-symmetric (_d.f90) or 
  ! complex-hermitian (_c.f90) cmatrix
  !+-------------------------------------------------------------------+
  include "lanczos_arpack.f90"
#ifdef _MPI
  include "lanczos_parpack.f90"
#endif

end module ARPACK_LANCZOS
