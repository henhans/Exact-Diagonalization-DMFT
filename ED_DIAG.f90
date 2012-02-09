!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!AUTHORS  : Adriano Amaricci
!New ordering of the sites:|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
!########################################################################
module ED_DIAG
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, ONLY:imp_sectorns,bdecomp,init_bath_ed
  USE ED_GETH
  USE ED_GETGF
  USE ED_GETOBS
  implicit none
  private
  public :: ed_solver

contains

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_solver(status)
    integer,optional :: status
    integer,save :: quo=-2

    if(present(status))quo=status

    select case(quo)
    case(-2)
       call msg(bold_green("status=-2 | INIT SOLVER, SETUP EIGENSPACE"))
       call init_bath_ed
       if(heff/=0.d0)then
          heff=-abs(heff)
          write(*,"(A,F12.9)")"Symmetry Breaking field = ",heff
          epsiup = epsiup - heff
          epsidw = epsidw + heff
          heff=0.d0
       endif
       call setup_eigenspace
       quo=0

    case(-1)
       call msg(bold_green("status=-1 | FINALIZE SOLVER"))
       call flush_eigenspace()
       call imp_diag
       call imp_getfunx
       if(chiflag)call imp_getchi
       call imp_getobs
       call dump_bath(Hfile)

    case default
       call msg(bold_green("status=0 | NORMAL"))
       call flush_eigenspace()
       call imp_diag
       call imp_getfunx
       call imp_getobs
       call dump_bath(Hfile)
    end select

  end subroutine ed_solver


  !*********************************************************************
  !*********************************************************************
  !*********************************************************************



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_eigenspace
    integer :: isloop,idg,jsloop
    call imp_setup
    !call getloop_range(1) 
    startloop=1;lastloop=Nsect
    startloop_mod=startloop;lastloop_mod=lastloop
    do isloop=startloop,lastloop
       jsloop=getCUPloop(isloop)
       if(jsloop==0)cycle
       if(startloop_mod > jsloop)startloop_mod=jsloop     
    enddo
    do isloop=startloop,lastloop
       jsloop=getCDWloop(isloop)
       if(jsloop==0)cycle
       if(startloop_mod > jsloop)startloop_mod=jsloop     
    enddo
    !Each CPU allocate the necessary and sufficient number of eigenspaces
    if(allocated(espace)) deallocate(espace)
    allocate(espace(startloop_mod:lastloop_mod))
    do isloop=startloop_mod,lastloop_mod
       idg=deg(isloop)
       allocate(espace(isloop)%e(idg),espace(isloop)%M(idg,idg))
    enddo
  end subroutine setup_eigenspace


  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


  subroutine flush_eigenspace
    integer :: isloop
    forall(isloop=startloop_mod:lastloop_mod)
       espace(isloop)%e=0.d0
       espace(isloop)%M=0.d0
    end forall
  end subroutine flush_eigenspace






  !+-------------------------------------------------------------------+
  !PURPOSE  : Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+------------------------------------------------------------------+
  subroutine imp_diag
    integer :: in,is,isloop,idg
    real(8),dimension(Nsect) :: e0 
    integer                  :: info
    integer                  :: lwork
    real(8),dimension(3*NP)  :: work !biggest dimension

    lwork=3*NP
    e0=0.d0
    call msg("Get Hamiltonian:")
    call start_timer
    !_mod are the "enlarged" sector range to take care of GF calculation
    do isloop=startloop_mod,lastloop_mod
       call eta(isloop,lastloop_mod,file="diag.eta")
       idg=deg(isloop)
       call imp_geth(isloop)
       call dsyev('V','U',idg,espace(isloop)%M,idg,espace(isloop)%e,work,lwork,info)
       if(info /= 0)print*,info
       !Avoid superpositions of sectors: startloop_mod < startloop
       if(isloop >=startloop)e0(isloop)=minval(espace(isloop)%e)
    enddo
    call stop_timer
    call findgs(e0)
    return
  end subroutine imp_diag



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************





  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine imp_searchmu(loop)
  !   integer          :: loop
  !   if(trim(order)=="co")then 
  !      call imp_searchmu_CO(loop) !model and searchmode NOT implemented yet
  !      return
  !   else
  !      call imp_searchmu_std(loop)
  !      return
  !   endif
  ! end subroutine imp_searchmu
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************



  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : SEARCHMU
  ! !TYPE     : subroutine
  ! !PURPOSE  : Solve the impurity problem by ED fix the density
  ! !+-------------------------------------------------------------------+
  ! !FIXDENSITY routines
  ! include "get_nobj.f90"
  ! include "imp_searchmu.f90"
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************






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
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************





  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine findgs(e0)
    integer :: i,isloop,idg
    real(8) :: egs
    real(8),dimension(Nsect) :: e0 

    !Get the GS energy and BCAST it
    egs=minval(e0)
    forall(isloop=startloop_mod:lastloop_mod)espace(isloop)%e = espace(isloop)%e - egs

    !Get the partition function Z and rescale energies
    zeta_function=0.d0;zeta_function=0.d0
    do isloop=startloop,lastloop
       idg=deg(isloop)
       do i=1,idg
          zeta_function=zeta_function+exp(-beta*espace(isloop)%e(i))
       enddo
    enddo

    write(*,"(A,f14.9)")'egs  =',egs
    write(*,"(A,f14.9)")'Z    =',zeta_function    
    open(3,file='egs.ed',access='append')
    write(3,*)egs
    close(3)
  end subroutine findgs
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine dump_bath(bath_file)
    character(len=*) :: bath_file
    integer :: i
    open(51,file=trim(bath_file))
    write(51,*)xmu
    do i=1,Nbath
       write(51,"(4(F13.9,1X))")epsiup(i),epsidw(i),vup(i),vdw(i)
    enddo
    close(51)
    return
  end subroutine dump_bath


end MODULE ED_DIAG
