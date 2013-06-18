!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine init_full_ed_solver(bath)
  real(8),dimension(:),intent(inout) :: bath
  integer                            :: i   
  call msg("INIT SOLVER, SETUP EIGENSPACE",unit=LOGfile)
  call check_bath_dimension(bath)
  call allocate_bath
  call init_bath_ed
  if(Nspin==2)then
     heff=abs(heff)
     if(mpiID==0)write(LOGfile,"(A,F12.9)")"Symmetry Breaking field = ",heff
     ebath(1,:) = ebath(1,:) + heff
     ebath(2,:) = ebath(2,:) - heff
     heff=0.d0
  endif
  call setup_pointers
  call setup_eigenspace
  call write_bath(LOGfile)
  bath = copy_bath()
  call deallocate_bath
  call msg("SET STATUS TO 0 in ED_SOLVER",unit=LOGfile)
end subroutine init_full_ed_solver

!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_solver(bath)
  real(8),dimension(:),intent(in) :: bath
  integer                         :: unit
  call msg("ED SOLUTION",unit=LOGfile)
  call check_bath_dimension(bath)
  call allocate_bath
  call set_bath(bath)
  call reset_eigenspace()
  call full_ed_diag
  call full_ed_getgf
  if(chiflag)call full_ed_getchi
  call full_ed_getobs
  unit=free_unit()
  open(unit,file=trim(Hfile))
  call write_bath(unit)
  close(unit)
  call deallocate_bath
end subroutine full_ed_solver

!+-------------------------------------------------------------------+
!PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
! GS, build the Green's functions calling all the necessary routines
!+------------------------------------------------------------------+
subroutine full_ed_diag
  integer :: in,is,isector,dim
  real(8),dimension(Nsect) :: e0 
  integer                  :: info,i,j
  integer                  :: lwork
  e0=0.d0
  call msg("Get Hamiltonian:",unit=LOGfile)
  call start_timer
  do isector=startloop,lastloop
     call eta(isector,lastloop,file="ETA_diag.ed")
     dim=getdim(isector)
     call full_ed_geth(isector,espace(isector)%M(:,:))
     call matrix_diagonalize(espace(isector)%M,espace(isector)%e,'V','U')
     if(isector >=startloop)e0(isector)=minval(espace(isector)%e)
  enddo
  call stop_timer
  call findgs(e0)
  return
end subroutine full_ed_diag

!+-------------------------------------------------------------------+
!PURPOSE  : 
!+-------------------------------------------------------------------+
subroutine findgs(e0)
  integer :: i,isector,dim
  real(8) :: egs
  real(8),dimension(Nsect) :: e0 
  egs=minval(e0)
  forall(isector=startloop:lastloop)espace(isector)%e = espace(isector)%e - egs
  !Get the partition function Z and rescale energies
  zeta_function=0.d0;zeta_function=0.d0
  do isector=startloop,lastloop
     dim=getdim(isector)
     do i=1,dim
        zeta_function=zeta_function+exp(-beta*espace(isector)%e(i))
     enddo
  enddo
  call msg("DIAG resume:",unit=LOGfile)
  if(mpiID==0)then
     write(LOGfile,"(A,f18.12)")'egs  =',egs
     write(LOGfile,"(A,f18.12)")'Z    =',zeta_function    
     write(LOGfile,*)""
     open(3,file='egs.ed',access='append')
     write(3,*)egs
     close(3)
  endif
end subroutine findgs

!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine setup_eigenspace
  integer :: isector,dim,jsector
  if(allocated(espace)) deallocate(espace)
  allocate(espace(startloop:lastloop))
  do isector=startloop,lastloop
     dim=getdim(isector)
     allocate(espace(isector)%e(dim),espace(isector)%M(dim,dim))
  enddo
end subroutine setup_eigenspace

!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine reset_eigenspace
  integer :: isector
  forall(isector=startloop:lastloop)
     espace(isector)%e=0.d0
     espace(isector)%M=0.d0
  end forall
end subroutine reset_eigenspace
