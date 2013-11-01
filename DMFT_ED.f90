module DMFT_ED
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_CHI2FIT
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  implicit none


contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine init_ed_solver(bath)
    real(8),dimension(:),intent(inout) :: bath
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"INIT SOLVER, SETUP EIGENSPACE"
#ifdef _MPI
    endif
#endif
    call init_ed_structure
    call check_bath_dimension(bath)
    call init_bath_ed
    call write_bath(LOGfile)
    call setup_pointers
    if(ed_type=='full')call setup_eigenspace
    bath = copy_bath()
    call deallocate_bath
  end subroutine init_ed_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_solver(bath)
    real(8),dimension(:),intent(in) :: bath
    integer                         :: unit
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"ED SOLUTION"
#ifdef _MPI
    endif
#endif
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)
    select case(ed_type)
    case default
       call lanc_ed_diag
       call lanc_ed_getgf
       if(chiflag)call lanc_ed_getchi
    case ('full')
       call reset_eigenspace()
       call full_ed_diag
       call full_ed_getgf
       if(chiflag)call full_ed_getchi
    end select
    call ed_getobs
#ifdef _MPI
    if(mpiID==0)then
#endif
       unit=free_unit()
       open(unit,file=trim(Hfile))
       call write_bath(unit)
       close(unit)
#ifdef _MPI
    endif
#endif
    call deallocate_bath
  end subroutine ed_solver


end module DMFT_ED
