!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  implicit none
  private

  public :: lanc_ed_diag
  public :: full_ed_diag
  public :: setup_eigenspace
  public :: reset_eigenspace
contains


  !+-------------------------------------------------------------------+
  !                    LANCZOS DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag
    integer             :: nup,ndw,isector,dim
    integer             :: nup0,ndw0,isect0,dim0,izero
    integer             :: info,i,j
    integer             :: Nitermax,Neigen,Nblock
    real(8)             :: oldzero,enemin,egs,emax
    real(8),allocatable :: eig_values(:)
    real(8),allocatable :: eig_basis(:,:)
    logical             :: lanc_solve
    !
    real(8) :: e
    !
    if(.not.groundstate%status)groundstate=es_init_espace()
    call es_free_espace(groundstate)
    !<finiteT
    if(.not.state_list%status)state_list=es_init_espace()
    call es_free_espace(state_list)
    !>finiteT
    oldzero=1000.d0
    emax=1000.d0
    numzero=0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_progress(LOGfile)
    sector: do isector=1,Nsect
       !call progress(isector,Nsect)
       dim     = getdim(isector)
       Neigen  = min(dim,neigen_sector(isector))!lanc_neigen)
       print*,isector,Neigen,dim
       Nitermax= min(dim,lanc_niter)
       Nblock  = min(dim,5*Neigen+2)
       !
       allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
       lanc_solve  = .true. ; if(Neigen==dim)lanc_solve=.false.
       !
       if(lanc_solve)then
          eig_values=0.d0 ; eig_basis=0.d0 
          call sp_init_matrix(spH0,dim)
          call ed_geth(isector)
          call lanczos_arpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_d,.false.)
       else
          eig_values=0.d0 ; eig_basis=0.d0 
          call ed_geth(isector,eig_basis)
          call matrix_diagonalize(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=1.d0
       endif

       !try to add each obtained state to the list if its energy is smaller than emax
       do i=1,Neigen
          call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,size=lanc_nstates,verbose=.true.)
       enddo

       !this find the ground state
       enemin=eig_values(1)
       if (enemin < oldzero-10.d-9) then
          numzero=1
          oldzero=enemin
          call es_free_espace(groundstate)
          call es_insert_state(groundstate,enemin,eig_basis(1:dim,1),isector)
       elseif(abs(enemin-oldzero) <= 1.d-9)then
          numzero=numzero+1
          if (numzero > Nsect)stop "ed_diag: too many gs"
          oldzero=min(oldzero,enemin)
          call es_insert_state(groundstate,enemin,eig_basis(1:dim,1),isector)
       endif
       !
       deallocate(eig_values,eig_basis)
       if(spH0%status)call sp_delete_matrix(spH0)
       !
    enddo sector
    call stop_progress
    do i=1,state_list%size
       e      = es_return_energy(state_list,i)
       isect0 = es_return_sector(state_list,i)
       nup0  = getnup(isect0)
       ndw0  = getndw(isect0)
       write(*,"(i3,f25.18,2i3)"),i,e,nup0,ndw0
    enddo
    !
    write(LOGfile,"(A)")"groundstate sector(s):"
    do izero=1,numzero
       isect0= es_return_sector(groundstate,izero)
       egs   = es_return_energy(groundstate,izero)
       nup0  = getnup(isect0)
       ndw0  = getndw(isect0)
       dim0  = getdim(isect0)
       write(LOGfile,"(1A6,f20.12,2I4)")'egs =',egs,nup0,ndw0
    enddo
    write(LOGfile,"(1A6,I4)")'Z   =',numzero
    open(3,file='egs.ed',access='append')
    write(3,*)egs
    close(3)
  end subroutine lanc_ed_diag



  !+-------------------------------------------------------------------+
  !                    FULL DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+-------------------------------------------------------------------+
  subroutine full_ed_diag
    integer                  :: info,i,j,in,is,isector,dim
    real(8),dimension(Nsect) :: e0 
    real(8)                  :: egs
    e0=0.d0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_progress(LOGfile)
    do isector=1,Nsect
       call progress(isector,Nsect)
       dim=getdim(isector)
       call ed_geth(isector,espace(isector)%M(:,:))
       call matrix_diagonalize(espace(isector)%M,espace(isector)%e,'V','U')
       e0(isector)=minval(espace(isector)%e)
    enddo
    call stop_progress
    !
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
    !!<MPI
    ! if(mpiID==0)then
    !!>MPI
    write(LOGfile,"(A,f18.12)")'egs  =',egs
    write(LOGfile,"(A,f18.12)")'Z    =',zeta_function    
    write(LOGfile,*)""
    open(3,file='egs.ed',access='append')
    write(3,*)egs
    close(3)
    !!<MPI
    ! endif
    !!>MPI
    return
  end subroutine full_ed_diag



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
