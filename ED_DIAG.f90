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

  public :: init_ed_solver
  public :: ed_solver
  public :: init_lanc_solver
  public :: lanc_solver

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine init_ed_solver(bath)
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
    call setup_impurity
    call setup_eigenspace
    call write_bath(LOGfile)
    bath = copy_bath()
    call deallocate_bath
    call msg("SET STATUS TO 0 in ED_SOLVER",unit=LOGfile)
  end subroutine init_ed_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine init_lanc_solver(bath)
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
    call setup_impurity
    call write_bath(LOGfile)
    bath = copy_bath()
    call deallocate_bath
    call msg("SET STATUS TO 0 in ED_SOLVER",unit=LOGfile)
  end subroutine init_lanc_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_solver(bath)
    real(8),dimension(:),intent(in) :: bath
    call msg("ED SOLUTION",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)
    call reset_eigenspace()
    call imp_diag
    call imp_getfunx
    if(chiflag)call imp_getchi
    call imp_getobs
    call dump_bath(Hfile)
    call deallocate_bath
  end subroutine ed_solver


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_solver(bath)
    real(8),dimension(:),intent(in) :: bath
    call msg("ED SOLUTION",unit=LOGfile)
    call check_bath_dimension(bath)
    call allocate_bath
    call set_bath(bath)
    call lanc_diag
    call lanc_getgf
    call lanc_getobs
    call dump_bath(Hfile)
    call deallocate_bath
  end subroutine lanc_solver




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_impurity
    integer                          :: in,is,idg,isloop,jn,js,jsloop
    integer                          :: ism
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    allocate(imap(NP),invmap(NN))
    isloop=0
    do in=1,N
       ism=in
       if(in.gt.Ns)ism=N-in
       do is=-ism,ism,2
          isloop=isloop+1
          call imp_sectorns(in,is,idg,imap,invmap)
          getloop(in,is)=isloop
          getin(isloop)=in
          getis(isloop)=is
          deg(isloop)=idg
          nmap(isloop,:)=imap
          invnmap(isloop,:)=invmap
       enddo
    enddo
    deallocate(imap,invmap)

    getCUPloop=0
    do isloop=1,Nsect
       if(isloop < getloop(2,0))cycle
       in=getin(isloop);is=getis(isloop)
       jn=in-1;js=is-1;if(abs(js) > jn)cycle
       jsloop=getloop(jn,js)
       getCUPloop(isloop)=jsloop
    enddo

    getCDGUPloop=0
    do isloop=1,Nsect
       if(isloop > getloop(N-1,-1))cycle
       in=getin(isloop);is=getis(isloop)
       jn=in+1;js=is+1;if(abs(js) > jn)cycle
       jsloop=getloop(jn,js)
       getCDGUPloop(isloop)=jsloop
    enddo

    getCDWloop=0
    do isloop=1,Nsect
       if(isloop < getloop(2,-2))cycle
       in=getin(isloop);is=getis(isloop)
       jn=in-1;js=is+1;if(abs(js) > jn)cycle
       jsloop=getloop(jn,js)
       getCDWloop(isloop)=jsloop
    enddo

    getCDGDWloop=0
    do isloop=1,Nsect
       if(isloop > getloop(N-1,1))cycle
       in=getin(isloop);is=getis(isloop)
       jn=in+1;js=is-1;if(abs(js) > jn)cycle
       jsloop=getloop(jn,js)
       getCDGDWloop(isloop)=jsloop
    enddo

    startloop=1;lastloop=Nsect
    do isloop=startloop,lastloop
       jsloop=getCUPloop(isloop)
       if(jsloop==0)cycle
       if(startloop > jsloop)startloop=jsloop     
    enddo
    do isloop=startloop,lastloop
       jsloop=getCDWloop(isloop)
       if(jsloop==0)cycle
       if(startloop > jsloop)startloop=jsloop     
    enddo
  end subroutine setup_impurity


  subroutine setup_eigenspace
    integer :: isloop,idg,jsloop
    if(allocated(espace)) deallocate(espace)
    allocate(espace(startloop:lastloop))
    do isloop=startloop,lastloop
       idg=deg(isloop)
       allocate(espace(isloop)%e(idg),espace(isloop)%M(idg,idg))
    enddo
  end subroutine setup_eigenspace


  subroutine reset_eigenspace
    integer :: isloop
    forall(isloop=startloop:lastloop)
       espace(isloop)%e=0.d0
       espace(isloop)%M=0.d0
    end forall
  end subroutine reset_eigenspace






  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine imp_diag
    integer :: in,is,isloop,idg
    real(8),dimension(Nsect) :: e0 
    integer                  :: info,i,j
    integer                  :: lwork
    e0=0.d0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_timer
    do isloop=startloop,lastloop
       call eta(isloop,lastloop,file="ETA_diag.ed")
       idg=deg(isloop)
       call imp_geth(isloop,espace(isloop)%M(:,:))
       call matrix_diagonalize(espace(isloop)%M,espace(isloop)%e,'V','U')
       if(isloop >=startloop)e0(isloop)=minval(espace(isloop)%e)
    enddo
    call stop_timer
    call findgs(e0)
    return
  end subroutine imp_diag




  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine lanc_diag
    integer                  :: in,is,isloop,idg
    integer                  :: info,i,j
    real(8),allocatable      :: eval(:),evec(:,:)
    integer                  :: Nitermax,Neigen,n0,s0,isect0,idg0,izero
    real(8)                  :: oldzero,enemin
    if(.not.allocated(iszero))allocate(iszero(Nsect))
    if(.not.allocated(groundstate))allocate(groundstate(Nsect))
    oldzero=1000.d0
    numzero=0
    iszero=0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_timer
    do isloop=startloop,lastloop
       call eta(isloop,lastloop,file="ETA_diag.ed")
       idg=deg(isloop)
       allocate(H0(idg,idg))
       call imp_geth(isloop,H0)
       if(idg>1)then
          Neigen=1
          Nitermax=min(idg,512)
          allocate(Eval(Neigen),Evec(Idg,Neigen))
          call lanczos_arpack(Idg,Neigen,Nitermax,eval,evec,HtimesV,.false.)
          enemin=eval(1)         
          if (enemin < oldzero-10.d-9) then
             numzero=1
             iszero(numzero)=isloop
             oldzero=enemin
             allocate(groundstate(numzero)%vec(idg))
             groundstate(numzero)%vec(1:idg)=evec(1:idg,1)
             groundstate(numzero)%egs=enemin
          elseif(abs(enemin-oldzero) <= 1.d-9)then
             numzero=numzero+1
             if (numzero > Nsect) stop 'too many gs'
             iszero(numzero)=isloop
             oldzero=min(oldzero,enemin)
             allocate(groundstate(numzero)%vec(idg))
             groundstate(numzero)%vec(1:idg)=evec(1:idg,1)
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
       idg0 = deg(isect0)
       print*,n0,s0,idg0,groundstate(izero)%egs
    enddo
    call stop_timer
  end subroutine lanc_diag






  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine findgs(e0)
    integer :: i,isloop,idg
    real(8) :: egs
    real(8),dimension(Nsect) :: e0 
    egs=minval(e0)
    forall(isloop=startloop:lastloop)espace(isloop)%e = espace(isloop)%e - egs

    !Get the partition function Z and rescale energies
    zeta_function=0.d0;zeta_function=0.d0
    do isloop=startloop,lastloop
       idg=deg(isloop)
       do i=1,idg
          zeta_function=zeta_function+exp(-beta*espace(isloop)%e(i))
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
