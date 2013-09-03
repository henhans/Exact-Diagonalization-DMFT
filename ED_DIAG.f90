!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
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
    call init_bath_ed
    if(Nspin==2)then
       heff=abs(heff)
       write(LOGfile,"(A,F12.9)")"Symmetry Breaking field = ",heff
       ebath(1,:,:)     = ebath(1,:,:)     + heff
       ebath(Nspin,:,:) = ebath(Nspin,:,:) - heff
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
    integer :: unit
    call msg("ED SOLUTION",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)
    call lanc_ed_diag
    call lanc_ed_getgf
    call lanc_ed_getobs
    unit=free_unit()
    open(unit,file=trim(Hfile))
    call write_bath(unit)
    close(unit)
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
    integer             :: Nitermax,Neigen,Nblock
    real(8)             :: oldzero,enemin,egs
    real(8),allocatable :: eig_values(:)
    real(8),allocatable :: eig_basis(:,:)
    logical             :: isolve

    if(.not.groundstate%status)groundstate=es_init_espace()
    call es_free_espace(groundstate)
    oldzero=1000.d0
    numzero=0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_progress(LOGfile)
    do isector=1,Nsect
       call progress(isector,Nsect)
       dim=getdim(isector)
       dim     = getdim(isector)
       Neigen  = min(dim,nLanceigen)
       Nitermax= min(dim,nLancitermax)
       Nblock  = max(nLancblock,5*Neigen+10)
       Nblock  = min(dim,Nblock)
       isolve  = .true.
       If((Neigen==Nitermax).AND.&
            (Neigen==Nblock).AND.&
            (Neigen==dim))isolve=.false.
       allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
       eig_values=0.d0 ; eig_basis=0.d0
       select case(isolve)
       case (.true.)
          ! !##IF SPARSE_MATRIX:
          call sp_init_matrix(spH0,dim)
          call lanc_ed_geth(isector)
          ! !##ELSE DIRECT H*V PRODUCT:
          ! call set_Hsector(isector)
          call lanczos_arpack(dim,Neigen,Nitermax,eig_values,eig_basis,spHtimesV_d,.false.)
       case (.false.)
          call full_ed_geth(isector,eig_basis)
          call matrix_diagonalize(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=1.d0
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
       nup0  = getnup(isect0)
       ndw0  = getndw(isect0)
       dim0  = getdim(isect0)
       write(LOGfile,"(1A6,f20.12,2I4)")'egs =',egs,nup0,ndw0
    enddo
    write(LOGfile,"(1A6,I4)")'Z   =',numzero
    open(3,file='egs.ed',access='append')
    write(3,*)egs
    close(3)
    call stop_progress
  end subroutine lanc_ed_diag



  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
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
       ebath(1,:,:)     = ebath(1,:,:)     + heff
       ebath(Nspin,:,:) = ebath(Nspin,:,:) - heff
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
    call start_progress(LOGfile)
    do isector=1,Nsect
       call progress(isector,Nsect)
       dim=getdim(isector)
       call full_ed_geth(isector,espace(isector)%M(:,:))
       call matrix_diagonalize(espace(isector)%M,espace(isector)%e,'V','U')
       e0(isector)=minval(espace(isector)%e)
    enddo
    call stop_progress
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
    forall(isector=1:Nsect)espace(isector)%e = espace(isector)%e - egs
    !Get the partition function Z and rescale energies
    zeta_function=0.d0;zeta_function=0.d0
    do isector=1,Nsect
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
    allocate(espace(1:Nsect))
    do isector=1,Nsect
       dim=getdim(isector)
       allocate(espace(isector)%e(dim),espace(isector)%M(dim,dim))
    enddo
  end subroutine setup_eigenspace


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine reset_eigenspace
    integer :: isector
    forall(isector=1:Nsect)
       espace(isector)%e=0.d0
       espace(isector)%M=0.d0
    end forall
  end subroutine reset_eigenspace

end MODULE ED_DIAG
