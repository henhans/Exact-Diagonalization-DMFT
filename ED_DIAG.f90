!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!New ordering of the sites:|ImpUP,BathUP;,ImpDW,BathDW >
!########################################################################
module ED_DIAG
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_GETH
  USE ED_GETGF
  USE ED_GETOBS
  USE ED_LANCZOS
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
       write(LOGfile,"(A,F12.9)")"Symmetry Breaking field = ",heff
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
    call msg("ED SOLUTION",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)
    call reset_eigenspace()
    call full_ed_diag
    call full_ed_getgf
    if(chiflag)call full_ed_getchi
    call full_ed_getobs
    call dump_bath(Hfile)
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
    integer                  :: in,is,isector,dim
    integer                  :: info,i,j
    real(8),allocatable      :: eval(:),evec(:,:)
    integer                  :: Nitermax,Neigen,n0,s0,isect0,dim0,izero
    real(8)                  :: oldzero,enemin
    if(.not.allocated(iszero))allocate(iszero(Nsect))
    if(.not.allocated(groundstate))allocate(groundstate(Nsect))
    oldzero=1000.d0
    numzero=0
    iszero=0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_timer
    do isector=startloop,lastloop
       call eta(isector,lastloop,file="ETA_diag.ed")
       dim=getdim(isector)
       allocate(H0(dim,dim))
       call full_ed_geth(isector,H0)
       if(dim>1)then
          Neigen=1
          Nitermax=min(dim,512)
          allocate(Eval(Neigen),Evec(Dim,Neigen))
          call lanczos_arpack(Dim,Neigen,Nitermax,eval,evec,HtimesV,.false.)
          enemin=eval(1)         
          if (enemin < oldzero-10.d-9) then
             numzero=1
             iszero(numzero)=isector
             oldzero=enemin
             allocate(groundstate(numzero)%vec(dim))
             groundstate(numzero)%vec(1:dim)=evec(1:dim,1)
             groundstate(numzero)%egs=enemin
          elseif(abs(enemin-oldzero) <= 1.d-9)then
             numzero=numzero+1
             if (numzero > Nsect) stop 'too many gs'
             iszero(numzero)=isector
             oldzero=min(oldzero,enemin)
             allocate(groundstate(numzero)%vec(dim))
             groundstate(numzero)%vec(1:dim)=evec(1:dim,1)
             groundstate(numzero)%egs=enemin
          endif
          deallocate(Eval,Evec)
       endif
       deallocate(H0)
    enddo
    print*,"numzero=",numzero
    print*,"groundstate sector:"
    do izero=1,numzero
       isect0 = iszero(izero)
       n0 = getin(isect0)
       s0 = getis(isect0)
       dim0 = getdim(isect0)
       print*,n0,s0,dim0,groundstate(izero)%egs
    enddo
    call stop_timer
  end subroutine lanc_ed_diag







  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: in,is,dim,isector,jn,js,jsector
    integer                          :: ism
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    allocate(imap(NP),invmap(NN))
    isector=0
    do in=1,N
       ism=in
       if(in.gt.Ns)ism=N-in
       do is=-ism,ism,2
          isector=isector+1
          call build_sector_ns(in,is,dim,imap,invmap)
          getsector(in,is)=isector
          getin(isector)=in
          getis(isector)=is
          getdim(isector)=dim
          Hmap(isector,:)=imap
          invHmap(isector,:)=invmap
       enddo
    enddo
    deallocate(imap,invmap)

    getCUPsector=0
    do isector=1,Nsect
       if(isector < getsector(2,0))cycle
       in=getin(isector);is=getis(isector)
       jn=in-1;js=is-1;if(abs(js) > jn)cycle
       jsector=getsector(jn,js)
       getCUPsector(isector)=jsector
    enddo

    getCDGUPsector=0
    do isector=1,Nsect
       if(isector > getsector(N-1,-1))cycle
       in=getin(isector);is=getis(isector)
       jn=in+1;js=is+1;if(abs(js) > jn)cycle
       jsector=getsector(jn,js)
       getCDGUPsector(isector)=jsector
    enddo

    getCDWsector=0
    do isector=1,Nsect
       if(isector < getsector(2,-2))cycle
       in=getin(isector);is=getis(isector)
       jn=in-1;js=is+1;if(abs(js) > jn)cycle
       jsector=getsector(jn,js)
       getCDWsector(isector)=jsector
    enddo

    getCDGDWsector=0
    do isector=1,Nsect
       if(isector > getsector(N-1,1))cycle
       in=getin(isector);is=getis(isector)
       jn=in+1;js=is-1;if(abs(js) > jn)cycle
       jsector=getsector(jn,js)
       getCDGDWsector(isector)=jsector
    enddo

    startloop=1;lastloop=Nsect
    do isector=startloop,lastloop
       jsector=getCUPsector(isector)
       if(jsector==0)cycle
       if(startloop > jsector)startloop=jsector     
    enddo
    do isector=startloop,lastloop
       jsector=getCDWsector(isector)
       if(jsector==0)cycle
       if(startloop > jsector)startloop=jsector     
    enddo
  end subroutine setup_pointers


  subroutine setup_eigenspace
    integer :: isector,dim,jsector
    if(allocated(espace)) deallocate(espace)
    allocate(espace(startloop:lastloop))
    do isector=startloop,lastloop
       dim=getdim(isector)
       allocate(espace(isector)%e(dim),espace(isector)%M(dim,dim))
    enddo
  end subroutine setup_eigenspace


  subroutine reset_eigenspace
    integer :: isector
    forall(isector=startloop:lastloop)
       espace(isector)%e=0.d0
       espace(isector)%M=0.d0
    end forall
  end subroutine reset_eigenspace
















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
    write(LOGfile,"(A,f18.12)")'egs  =',egs
    write(LOGfile,"(A,f18.12)")'Z    =',zeta_function    
    write(LOGfile,*)""

    open(3,file='egs.ed',access='append')
    write(3,*)egs
    close(3)
  end subroutine findgs



end MODULE ED_DIAG
