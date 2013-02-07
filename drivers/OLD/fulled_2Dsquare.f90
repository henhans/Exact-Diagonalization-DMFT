program CuO2fulled
  !########################################################################
  !PROGRAM  : HMED
  !TYPE     : main code
  !PURPOSE  : Complete ED with STAR GEOMETRY for CuO2 model (Emery-OKA)
  !AUTHORS  : Adriano Amaricci, L. de'Medici,G.Sordi, M.Rozenberg
  !########################################################################
  !LOCAL:
  USE VARS_GLOBAL
  USE DIAG
  USE TOFITGF
  implicit none
  integer :: success

  call initialize("inputED.in")

  call Build_2DSquareLattice(Nx,Ny,Nk=Lk)
  allocate(epsik(Lk)) ; call get_epsik(epsik,ts,tsp)

  !Store initial parameters:
  xmu0=xmu+ed0;xmu=xmu0


  !Starts DMFT iteration
  do iloop=1,nloop
     call ed_solver(iloop,.true.)
     call get_delta()
     success = scc_fit(iloop)
  enddo
  call finalize(success)
contains
  !+----------------------------------------+
  subroutine get_delta
    call getdelta_2dsquare
  end subroutine get_delta
  !+----------------------------------------+
  include "get_delta_samples.f90"
  !+----------------------------------------+
end program



