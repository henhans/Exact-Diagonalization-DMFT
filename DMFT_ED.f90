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
  subroutine init_ed_solver(bath_)
    real(8),dimension(:),intent(inout) :: bath_
    if(mpiID==0)write(LOGfile,"(A)")"INIT SOLVER, SETUP EIGENSPACE"
    bath_=0.d0
    call init_ed_structure
    call check_bath_dimension(bath_)
    call allocate_bath(dmft_bath)
    call init_bath_ed(dmft_bath)
    call copy_bath(dmft_bath,bath_)
    call write_bath(dmft_bath,LOGfile)
    call setup_pointers
    if(ed_method=='full')call setup_eigenspace
    call deallocate_bath(dmft_bath)
  end subroutine init_ed_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_solver(bath_)
    real(8),dimension(:),intent(in) :: bath_
    integer                         :: unit
    if(mpiID==0)write(LOGfile,"(A)")"ED SOLUTION"
    call check_bath_dimension(bath_)
    call allocate_bath(dmft_bath)
    call set_bath(bath_,dmft_bath)
    select case(ed_method)
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
    if(mpiID==0)then
       unit=free_unit()
       open(unit,file=trim(Hfile))
       call write_bath(dmft_bath,unit)
       close(unit)
    endif
    call deallocate_bath(dmft_bath)
  end subroutine ed_solver


end module DMFT_ED
