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
  implicit none
  private

  public :: init_ed_solver
  public :: ed_solver

contains


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
    call setup_eigenspace
    call write_bath(LOGfile)
    bath = copy_bath()
    call deallocate_bath
    call msg("SET STATUS TO 0 in ED_SOLVER",unit=LOGfile)
  end subroutine init_ed_solver





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
  subroutine setup_eigenspace
    integer :: isloop,idg,jsloop
    call imp_setup
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
    integer                  :: info
    integer                  :: lwork
    e0=0.d0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_timer
    do isloop=startloop,lastloop
       call eta(isloop,lastloop,file="ETA_diag.ed")
       idg=deg(isloop)
       call imp_geth(isloop)
       call matrix_diagonalize(espace(isloop)%M,espace(isloop)%e,'V','U')
       if(isloop >=startloop)e0(isloop)=minval(espace(isloop)%e)
    enddo
    call stop_timer
    call findgs(e0)
    return
  end subroutine imp_diag





  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_setup
    integer :: in,is,idg,isloop,jn,js,jsloop
    integer :: ism
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

    getCDWloop=0
    do isloop=1,Nsect
       if(isloop < getloop(2,-2))cycle
       in=getin(isloop);is=getis(isloop)
       jn=in-1;js=is+1;if(abs(js) > jn)cycle
       jsloop=getloop(jn,js)
       getCDWloop(isloop)=jsloop
    enddo
    return
  end subroutine imp_setup





  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
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
