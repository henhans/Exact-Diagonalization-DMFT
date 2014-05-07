!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!|{ImpUP1,...,ImpUPN},BathUP>|{ImpDW1,...,ImpDWN},BathDW>
!########################################################################
module ED_DIAG
  USE CONSTANTS
  USE MATRIX, only: matrix_diagonalize
  USE TIMER
  USE IOTOOLS, only:reg,free_unit
  USE STATISTICS
  !
  USE ARPACK_LANCZOS
  USE ED_INPUT_VARS
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

  !                    LANCZOS DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag
    logical :: iverbose_
    select case(ed_type)
    case default
       call lanc_ed_diag_d
    case('c')
       call lanc_ed_diag_c
    end select
  end subroutine lanc_ed_diag




  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag_d
    integer             :: nup,ndw,isector,dim
    integer             :: nup0,ndw0,isect0,dim0,izero,sz0
    integer             :: i,j,unit
    integer             :: Nitermax,Neigen,Nblock
    real(8)             :: oldzero,enemin,Egs,Ei,Ec
    real(8),allocatable :: eig_values(:)
    real(8),allocatable :: eig_basis(:,:)
    logical             :: lanc_solve
    type(histogram)     :: hist
    real(8)             :: hist_a,hist_b,hist_w
    integer             :: hist_n
    integer,allocatable :: list_sector(:),count_sector(:)
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    !call es_free_espace(state_list)
    oldzero=1000.d0
    numgs=0
    if(ed_verbose<2)call start_progress(LOGfile)
    sector: do isector=1,Nsect
       if(ed_verbose<1)call progress(isector,Nsect)
       dim     = getdim(isector)
       Neigen  = min(dim,neigen_sector(isector))
       Nitermax= min(dim,lanc_niter)
       Nblock  = min(dim,5*Neigen+10)
       !
       lanc_solve  = .true. ; if(Neigen==dim)lanc_solve=.false.
       if(dim<=128)lanc_solve=.false.
       !
       if(lanc_solve)then
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          eig_values=0.d0 ; eig_basis=0.d0
          call ed_buildH_d(isector)
          call lanczos_arpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_dd,.false.)
       else
          allocate(eig_values(Dim),eig_basis(Dim,dim))
          eig_values=0.d0 ; eig_basis=0.d0 
          call ed_buildH_d(isector,eig_basis)
          call matrix_diagonalize(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=1.d0
       endif
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,size=lanc_nstates_total)
          enddo
       else
          enemin = eig_values(1)
          if (enemin < oldzero-10.d-9) then
             numgs=1
             oldzero=enemin
             call es_free_espace(state_list)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector)
          elseif(abs(enemin-oldzero) <= 1.d-9)then
             numgs=numgs+1
             if (numgs > Nsect)stop "ed_diag: too many gs"
             oldzero=min(oldzero,enemin)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector)
          endif
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       if(spH0%status)call sp_delete_matrix(spH0)
       !
    enddo sector
    if(ed_verbose<2)call stop_progress

    !POST PROCESSING:
    if(ed_verbose<2)then
       unit=free_unit()
       open(unit,file="state_list"//reg(ed_file_suffix)//".ed")
       if(.not.ed_supercond)then
          write(unit,"(A)")"#i       E_i             DE_i             nup ndw Sect  Dim"
       else
          write(unit,"(A)")"#i       E_i             DE_i             Sz    Sect    Dim"
       endif
       do i=1,state_list%size
          Ei     = es_return_energy(state_list,i)
          isect0 = es_return_sector(state_list,i)
          if(.not.ed_supercond)then
             nup0   = getnup(isect0)
             ndw0   = getndw(isect0)
             write(unit,"(i3,f18.12,E18.9,1x,2i3,3x,i3,i10)"),i,Ei,exp(-beta*(Ei-state_list%emin)),nup0,ndw0,isect0,getdim(isect0)
          else
             sz0   = getsz(isect0)
             write(unit,"(i3,f18.12,E18.9,1x,i3,3x,i3,i10)"),i,Ei,exp(-beta*(Ei-state_list%emin)),sz0,isect0,getdim(isect0)
          endif
       enddo
       close(unit)
    endif
    zeta_function=0.d0
    Egs = state_list%emin
    if(finiteT)then
       do i=1,state_list%size
          ei            = es_return_energy(state_list,i)
          zeta_function = zeta_function + exp(-beta*(Ei-Egs))
       enddo
    else
       if(numgs/=state_list%size)stop "ED_DIAG: error in evaluating Numgs!"
       zeta_function=real(numgs,8)
    end if
    !
    if(finiteT)then
       numgs=es_return_groundstates(state_list)
       if(numgs>Nsect)stop "ed_diag: too many gs"
    endif
    do izero=1,numgs
       isect0= es_return_sector(state_list,izero)
       Egs   = es_return_energy(state_list,izero)
       dim0  = getdim(isect0)
       if(.not.ed_supercond)then
          nup0  = getnup(isect0)
          ndw0  = getndw(isect0)
          if(ed_verbose<3)write(LOGfile,"(A,F20.12,2I4)")'egs =',egs,nup0,ndw0
       else
          sz0  = getsz(isect0)
          if(ed_verbose<3)write(LOGfile,"(A,F20.12,I4)")'egs =',egs,sz0
       endif
    enddo
    if(ed_verbose<3)write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    ! open(3,file='egs'//reg(ed_file_suffix)//".ed",access='append')
    ! write(3,*)egs
    ! close(3)


    !Get histogram distribution of the sector contributing to the evaluated spectrum:
    !Go thru states list and update the neigen_sector(isector) sector-by-sector
    if(ed_verbose<2.AND.finiteT)then
       unit=free_unit()
       open(unit,file="histogram_states"//reg(ed_file_suffix)//".ed",access='append')
       hist_n = Nsect
       hist_a = 1.d0
       hist_b = real(Nsect,8)
       hist_w = 1.d0
       hist = histogram_allocate(hist_n)
       call histogram_set_range_uniform(hist,hist_a,hist_b)
       do i=1,state_list%size
          isect0 = es_return_sector(state_list,i)
          call histogram_accumulate(hist,dble(isect0),hist_w)
       enddo
       call histogram_print(hist,unit)
       write(unit,*)""
       close(unit)
       !
       allocate(list_sector(state_list%size),count_sector(Nsect))
       !get the list of actual sectors contributing to the list
       do i=1,state_list%size
          list_sector(i) = es_return_sector(state_list,i)
       enddo
       !count how many times a sector appears in the list
       do i=1,Nsect
          count_sector(i) = count(list_sector==i)
       enddo
       !adapt the number of required Neig for each sector based on how many
       !appeared in the list.
       do i=1,Nsect
          if(any(list_sector==i))then !if(count_sector(i)>1)then
             neigen_sector(i)=neigen_sector(i)+1
          else
             neigen_sector(i)=neigen_sector(i)-1
          endif
          !prevent Neig(i) from growing unbounded but 
          !try to put another state in the list from sector i
          if(neigen_sector(i) > count_sector(i))neigen_sector(i)=count_sector(i)+1
          if(neigen_sector(i) <= 0)neigen_sector(i)=1
       enddo
       !check if the number of states is enough to reach the required accuracy:
       !the condition to fullfill is:
       ! exp(-beta(Ec-Egs)) < \epsilon_c
       Egs = state_list%emin
       Ec  = state_list%emax
       if(exp(-beta*(Ec-Egs)) > cutoff)then
          lanc_nstates_total=lanc_nstates_total + 2!*lanc_nincrement
          write(*,"(A,I4)")"Increasing lanc_nstates_total+2:",lanc_nstates_total
       endif
    endif
  end subroutine lanc_ed_diag_d



  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine lanc_ed_diag_c
    integer                :: nup,ndw,isector,dim
    integer                :: nup0,ndw0,isect0,dim0,izero,sz0
    integer                :: i,j,unit
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin,Egs,Ei,Ec
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:)
    logical                :: lanc_solve
    type(histogram)        :: hist
    real(8)                :: hist_a,hist_b,hist_w
    integer                :: hist_n
    integer,allocatable    :: list_sector(:),count_sector(:)
    if(.not.state_list%status)state_list=es_init_espace()
    call es_free_espace(state_list)
    oldzero=1000.d0
    numgs=0
    write(LOGfile,"(A)")"Get Hamiltonian:"
    if(ed_verbose<2)call start_progress(LOGfile)
    sector: do isector=1,Nsect
       if(ed_verbose<1)call progress(isector,Nsect)
       dim     = getdim(isector)
       Neigen  = min(dim,neigen_sector(isector))
       Nitermax= min(dim,lanc_niter)
       Nblock  = min(dim,5*Neigen+10)
       !
       lanc_solve  = .true. ; if(Neigen==dim)lanc_solve=.false.
       if(dim<=128)lanc_solve=.false.
       !
       if(lanc_solve)then
          allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
          eig_values=0.d0 ; eig_basis=zero
          call ed_buildH_c(isector)
          call lanczos_arpack(dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,spHtimesV_cc,.false.)
       else
          allocate(eig_values(Dim),eig_basis(Dim,dim))
          eig_values=0.d0 ; eig_basis=zero
          call ed_buildH_c(isector,eig_basis)
          call matrix_diagonalize(eig_basis,eig_values,'V','U')
          if(dim==1)eig_basis(dim,dim)=one
       endif
       !
       if(finiteT)then
          do i=1,Neigen
             call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,size=lanc_nstates_total)
          enddo
       else
          enemin = eig_values(1)
          if (enemin < oldzero-10.d-9) then
             numgs=1
             oldzero=enemin
             call es_free_espace(state_list)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector)
          elseif(abs(enemin-oldzero) <= 1.d-9)then
             numgs=numgs+1
             if (numgs > Nsect)stop "ed_diag: too many gs"
             oldzero=min(oldzero,enemin)
             call es_insert_state(state_list,enemin,eig_basis(1:dim,1),isector)
          endif
>>>>>>> devel
       endif
       !
       if(allocated(eig_values))deallocate(eig_values)
       if(allocated(eig_basis))deallocate(eig_basis)
       if(spH0%status)call sp_delete_matrix(spH0)
       !
    enddo sector
    if(ed_verbose<2)call stop_progress
    !POST PROCESSING:
    if(ed_verbose<2)then
       unit=free_unit()
       open(unit,file="state_list"//reg(ed_file_suffix)//".ed")
       write(unit,"(A)")"#i       E_i                nup ndw"
       do i=1,state_list%size
          Ei     = es_return_energy(state_list,i)
          isect0 = es_return_sector(state_list,i)
          if(.not.ed_supercond)then
             nup0   = getnup(isect0)
             ndw0   = getndw(isect0)
             write(unit,"(i3,f25.18,2i3,3x,2i3)"),i,Ei,nup0,ndw0,isect0,getdim(isect0)
          else
             sz0   = getsz(isect0)
             write(unit,"(i3,f25.18,i3,3x,2i3)"),i,Ei,sz0,isect0,getdim(isect0)
          endif
       enddo
       close(unit)
    endif
    zeta_function=0.d0
    Egs = state_list%emin
    if(finiteT)then
       do i=1,state_list%size
          ei            = es_return_energy(state_list,i)
          zeta_function = zeta_function + exp(-beta*(Ei-Egs))
       enddo
    else
       if(numgs/=state_list%size)stop "ED_DIAG: error in evaluating Numgs!"
       zeta_function=real(numgs,8)
    end if
    !
    if(finiteT)then
       numgs=es_return_groundstates(state_list)
       if(numgs>Nsect)stop "ed_diag: too many gs"
    endif
    do izero=1,numgs
       isect0= es_return_sector(state_list,izero)
       Egs   = es_return_energy(state_list,izero)
       dim0  = getdim(isect0)
       if(.not.ed_supercond)then
          nup0  = getnup(isect0)
          ndw0  = getndw(isect0)
          if(ed_verbose<3)write(LOGfile,"(A,F20.12,2I4)")'egs =',egs,nup0,ndw0
       else
          sz0  = getsz(isect0)
          if(ed_verbose<3)write(LOGfile,"(A,F20.12,I4)")'egs =',egs,sz0
       endif
    enddo
    if(ed_verbose<3)write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    ! open(3,file='egs'//reg(ed_file_suffix)//".ed",access='append')
    ! write(3,*)egs
    ! close(3)
    !Get histogram distribution of the sector contributing to the evaluated spectrum:
    !Go thru states list and update the neigen_sector(isector) sector-by-sector
    if(ed_verbose<2.AND.finiteT)then
       unit=free_unit()
       open(unit,file="histogram_states"//reg(ed_file_suffix)//".ed",access='append')
       hist_n = Nsect
       hist_a = 1.d0
       hist_b = real(Nsect,8)
       hist_w = 1.d0!/real(state_list%size,8)
       hist = histogram_allocate(hist_n)
       call histogram_set_range_uniform(hist,hist_a,hist_b)
       do i=1,state_list%size
          isect0 = es_return_sector(state_list,i)
          call histogram_accumulate(hist,dble(isect0),hist_w)
       enddo
       call histogram_print(hist,unit)
       write(unit,*)""
       close(unit)
       allocate(list_sector(state_list%size),count_sector(Nsect))
       !get the list of actual sectors contributing to the list
       do i=1,state_list%size
          list_sector(i) = es_return_sector(state_list,i)
       enddo
       !count how many times a sector appears in the list
       do i=1,Nsect
          count_sector(i) = count(list_sector==i)
       enddo
       !adapt the number of required Neig for each sector based on how many
       !appeared in the list.
       do i=1,Nsect
          if(any(list_sector==i))then !if(count_sector(i)>1)then
             neigen_sector(i)=neigen_sector(i)+1
          else
             neigen_sector(i)=neigen_sector(i)-1
          endif
          !prevent Neig(i) from growing unbounded but 
          !try to put another state in the list from sector i
          if(neigen_sector(i) > count_sector(i))neigen_sector(i)=count_sector(i)+1
          if(neigen_sector(i) <= 0)neigen_sector(i)=1
       enddo
       !check if the number of states is enough to reach the required accuracy:
       !the condition to fullfill is:
       ! exp(-beta(Ec-Egs)) < \epsilon_c
       Egs = state_list%emin
       Ec  = state_list%emax
       if(exp(-beta*(Ec-Egs)) > cutoff)then
          lanc_nstates_total=lanc_nstates_total + 2!*lanc_nincrement
          write(*,"(A,I4)")"Increasing lanc_nstates_total+2:",lanc_nstates_total
       endif
    endif
  end subroutine lanc_ed_diag_c




  !+-------------------------------------------------------------------+
  !                    FULL DIAGONALIZATION
  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+-------------------------------------------------------------------+
  subroutine full_ed_diag
    integer                  :: i,j,in,is,isector,dim
    real(8),dimension(Nsect) :: e0 
    real(8)                  :: egs
    e0=0.d0
    write(LOGfile,"(A)")"Get Hamiltonian:"
    call start_progress(LOGfile)
    do isector=1,Nsect
       call progress(isector,Nsect)
       dim=getdim(isector)
       !call setup_Hv_sector(isector)
       call ed_buildH_d(isector,espace(isector)%M(:,:))
       !call delete_Hv_sector()
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
    write(LOGfile,"(A)")"DIAG resume:"
    write(LOGfile,"(A,f20.12)")'egs  =',egs
    write(LOGfile,"(A,f20.12)")'Z    =',zeta_function    
    open(3,file='egs'//reg(ed_file_suffix)//".ed",access='append')
    write(3,*)egs
    close(3)
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
