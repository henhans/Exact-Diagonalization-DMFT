module ARPACK_LANCZOS
  USE CONSTANTS, only:zero
  USE ED_VARS_GLOBAL
  implicit none
  private 

  interface lanczos_arpack
     module procedure lanczos_arpack_d,lanczos_arpack_c
  end interface lanczos_arpack
  public :: lanczos_arpack

contains

  !+-------------------------------------------------------------------+
  ! This routine use ARPACK to find a few eigenvalues
  ! LAMBDA and corresponding eigenvectors X for the standard
  ! eigenvalue problem:
  !      A * X = LAMBDA * X
  ! where A is an N by N 
  ! real-symmetric or complex-hermitian cmatrix
  !+-------------------------------------------------------------------+
  include "lanczos_arpack.f90"

end module ARPACK_LANCZOS
