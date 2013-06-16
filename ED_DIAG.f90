!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!New ordering of the sites:|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
!########################################################################
module ED_DIAG
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_GETH
  USE ED_GETGF
  USE ED_GETOBS
  implicit none
  private

  public :: init_full_ed_solver
  public :: full_ed_solver
  public :: init_lanc_ed_solver
  public :: lanc_ed_solver

contains

  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  include 'fulled_diag.f90'



  !####################################################################
  !                    LANCZOS DIAGONALIZATION (T=0, GS only)
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine init_lanc_ed_solver(bath)
    real(8),dimension(:),intent(inout) :: bath
    integer                            :: i   
    call msg("INIT SOLVER, SETUP EIGENSPACE",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call init_bath_ed
    if(Nspin==2)then
       heff=abs(heff)
       write(LOGfile,"(A,F12.9)")"Symmetry Breaking field = ",heff
       ebath(1,:) = ebath(1,:) + heff
       ebath(2,:) = ebath(2,:) - heff
       heff=0.d0
    endif
    call setup_pointers
    call write_bath(LOGfile)
    bath = copy_bath()
    call deallocate_bath
    call msg("SET STATUS TO 0 in ED_SOLVER",unit=LOGfile)
  end subroutine init_lanc_ed_solver

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_solver(bath)
    real(8),dimension(:),intent(in) :: bath
    call msg("ED SOLUTION",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)
    call lanc_ed_diag
    call lanc_ed_getgf
    call lanc_ed_getobs
    call dump_bath(Hfile)
    call deallocate_bath
  end subroutine lanc_ed_solver



!+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag
    integer             :: nup,ndw,isector,dim
    integer             :: nup0,ndw0,isect0,dim0,izero
    integer             :: info,i,j
    integer             :: Nitermax,Neigen
    real(8)             :: oldzero,enemin,egs
    real(8),allocatable :: eig_values(:)
    real(8),allocatable :: eig_basis(:,:)
    if(.not.groundstate%status)groundstate=es_init_espace()
    call es_free_espace(groundstate)
    oldzero=1000.d0
    numzero=0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_timer
    do isector=startloop,lastloop
       call eta(isector,lastloop,file="ETA_diag.ed")
       !Get Hamiltonian (to be changed)
       dim=getdim(isector)
       ! !##IF SPARSE_MATRIX:
       call sp_init_matrix(spH0,dim)
       call lanc_ed_geth(isector)
       ! !##ELSE DIRECT H*V PRODUCT:
       ! call set_Hsector(isector)
       select case(dim)
       case default
          Neigen=1
          Nitermax=min(dim,nLancitermax)
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          call lanczos_arpack(dim,Neigen,Nitermax,eig_values,eig_basis,spHtimesV,.false.)
       case (1)
          allocate(eig_values(dim),eig_basis(dim,dim))
          ! !##IF SPARSE_MATRIX:
          eig_values(dim) =sp_get_element(spH0,dim,dim)
          ! !##ELSE DIRECT H*V PRODUCT:
          ! call full_ed_geth(isector,eig_basis)
          ! eig_values(dim)=eig_basis(dim,dim)
          eig_basis(1,1)=1.d0
       end select
       enemin=eig_values(1)  
       if (enemin < oldzero-10.d-9) then
          numzero=1
          oldzero=enemin
          call es_free_espace(groundstate)
          call es_insert_state(groundstate,enemin,eig_basis(1:dim,1),isector)
       elseif(abs(enemin-oldzero) <= 1.d-9)then
          numzero=numzero+1
          if (numzero > Nsect)call error('too many gs')
          oldzero=min(oldzero,enemin)
          call es_insert_state(groundstate,enemin,eig_basis(1:dim,1),isector)
       endif
       deallocate(eig_values,eig_basis)
       !Delete Hamiltonian matrix:
       ! !##IF SPARSE_MATRIX:
       call sp_delete_matrix(spH0)
    enddo
    !
    write(LOGfile,"(A)")"groundstate sector(s):"
    do izero=1,numzero
       isect0= es_get_sector(groundstate,izero)
       egs   = es_get_energy(groundstate,izero)
       nup0    = getnup(isect0)
       ndw0    = getndw(isect0)
       dim0  = getdim(isect0)
       write(LOGfile,"(A,f18.12,2I4)")'egs =',egs,nup0,ndw0
    enddo
    write(LOGfile,"(A,f18.12)")'Z   =',dble(numzero)
    open(3,file='egs.ed',access='append')
    write(3,*)egs
    close(3)
    call stop_timer
  end subroutine lanc_ed_diag






end MODULE ED_DIAG
