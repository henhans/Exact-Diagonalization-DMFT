program fulled_rerun
  USE DMFT_FULLED
  implicit none
  integer                :: Nb
  !Bath:
  real(8),allocatable    :: Bath(:)
  call read_input("used.inputED.in")
  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_ed_solver(bath)
  !Solve the EFFECTIVE IMPURITY PROBLEM only
  call ed_solver(bath) 
end program fulled_rerun



