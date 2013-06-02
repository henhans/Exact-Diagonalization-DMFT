!########################################################################
!PURPOSE  : Diagonalize the Effective Impurity Problem
!New ordering of the sites:|ImpUP,BathUP;,ImpDW,BathDW >
!########################################################################
include "arpack_lanczos.f90"
include "plain_lanczos.f90" 
module ED_GSVEC
  USE ED_VARS_GLOBAL
  implicit none
  integer :: numzero
  integer,dimension(:),allocatable :: iszero
  type gstate
     real(8)                      :: egs
     real(8),dimension(:),pointer :: vec
  end type gstate
  type(gstate),dimension(:),allocatable :: groundstate
end module ED_GSVEC


module ED_DIAG
  USE ED_VARS_GLOBAL
  USE ED_GSVEC
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_GETH
  USE ED_GETGF
  USE ED_GETOBS
  use lanczos_simple
  implicit none
  private

  public :: init_ed_solver
  public :: ed_solver

  real(8),dimension(:,:),allocatable :: H0
  real(8),dimension(:),allocatable :: wm,tau,wr
  complex(8),dimension(:,:),allocatable :: lanc_Giw,lanc_Gwr

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
    integer                  :: info,i,j
    integer                  :: lwork
    real(8),allocatable :: eval(:),evec(:,:)
    integer :: Nitermax,Neigen,n0,s0,isect0,idg0,izero
    real(8) :: oldzero,enemin
    !
    !
    !
    !#LANCZOS
    if(.not.allocated(iszero))allocate(iszero(Nsect))
    if(.not.allocated(groundstate))allocate(groundstate(Nsect))
    oldzero=1000.d0
    numzero=0
    iszero=0
    !#LANCZOS
    !
    !
    !
    e0=0.d0
    call msg("Get Hamiltonian:",unit=LOGfile)
    call start_timer
    do isloop=startloop,lastloop
       call eta(isloop,lastloop,file="ETA_diag.ed")
       idg=deg(isloop)
       call imp_geth(isloop,espace(isloop)%M(:,:))
       call matrix_diagonalize(espace(isloop)%M,espace(isloop)%e,'V','U')
       if(isloop >=startloop)e0(isloop)=minval(espace(isloop)%e)
       !
       !
       !
       !#LANCZOS
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
          write(*,*)isloop,eval(1),getin(isloop),getis(isloop)
          deallocate(Eval,Evec)
       endif
       deallocate(H0)
       !#LANCZOS
       !
       !
       !
    enddo

    !
    !
    !
    !#LANCZOS    
    print*,"numzero=",numzero
    print*,"groundstate sector:"
    do izero=1,numzero
       isect0 = iszero(izero)
       n0 = getin(isect0)
       s0 = getis(isect0)
       idg0 = deg(isect0)
       print*,n0,s0,idg0
       do j=1,idg0
          print*,groundstate(izero)%vec(j),espace(isect0)%M(j,1)
       enddo
       print*,groundstate(izero)%egs
    enddo
    !#LANCZOS
    !
    !
    !
    call stop_timer
    call findgs(e0)

    !
    !
    !
    !#LANCZOS   
    call lanc_getgf()
    !#LANCZOS
    !
    !
    !
    return
  end subroutine imp_diag

  !
  !
  !
  !#LANCZOS
  subroutine HtimesV(N,v,Hv)
    real(8),dimension(N) :: v
    real(8),dimension(N) :: Hv
    integer              :: N
    call dgemv('N',n,n,1.d0,H0,n,v,1,0.d0,Hv,1)
  end subroutine HtimesV

  subroutine lanc_getgf()
    integer :: i,izero,isect0,jsect0,m,j
    integer :: in0,is0,idg0
    integer :: jn0,js0,jdg0
    real(8) :: norm0,sgn,gs,nup,ndw
    integer :: ib(N),k,r,Nlanc,Nitermax
    real(8) :: factor
    real(8),allocatable :: vvinit(:),alfa_(:),beta_(:),vout(:)

    call plain_lanczos_set_htimesv(HtimesV)

    allocate(lanc_giw(Nspin,NL),lanc_gwr(Nspin,Nw))
    lanc_giw=zero
    lanc_gwr=zero

    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)

    Nitermax=250!min(jdg0,512)
    allocate(alfa_(Nitermax),beta_(Nitermax))

    factor=real(numzero,8)

    nsimp  = 0.d0
    nupimp = 0.d0
    ndwimp = 0.d0
    dimp   = 0.d0
    magimp = 0.d0
    m2imp  = 0.d0

    do izero=1,numzero   
       !GET THE GROUNDSTATE (make some checks)
       norm0=sqrt(dot_product(&
            groundstate(izero)%vec,groundstate(izero)%vec))
       if(norm0-1.d0>1.d-9)print*,"GS",izero,"is not normalized:",norm0
       isect0 = iszero(izero)
       in0    = getin(isect0)
       is0    = getis(isect0)
       idg0   = deg(isect0)


       do i=1,idg0
          m=nmap(isect0,i)
          call bdecomp(m,ib)
          nup=real(ib(1),8)
          ndw=real(ib(1+Ns),8)
          gs=groundstate(izero)%vec(i)
          nsimp  = nsimp  +  (nup+ndw)*gs**2
          nupimp = nupimp +  (nup)*gs**2
          ndwimp = ndwimp +  (ndw)*gs**2
          dimp   = dimp   +  (nup*ndw)*gs**2
          magimp = magimp +  (nup-ndw)*gs**2
          m2imp  = m2imp  +  gs**2*(nup-ndw)**2
       enddo


       !ADD ONE PARTICLE UP:
       !get cdg_up sector informations:
       jsect0 = getCDGUPloop(isect0);if(jsect0==0)cycle
       jdg0   = deg(jsect0)
       jn0    = getin(jsect0)
       js0    = getis(jsect0)
       print*,'sector C^+_up|gs>',jn0,js0,jdg0
       !allocate cdg_ip|gs> vector:
       allocate(vvinit(jdg0));vvinit=0.d0
       !build cdg_up|gs> vector:
       do m=1,idg0              !loop over |gs> components m
          i=nmap(isect0,m)      !map m to full-Hilbert space state i
          call bdecomp(i,ib)    !decompose i into number representation ib=|1/0,1/0,1/0...>
          if(ib(1)==0)then      !if impurity is empty: proceed
             call cdg(1,i,r);sgn=dfloat(r)/dfloat(abs(r));r=abs(r) !apply cdg_up (1), bring from i to r
             j=invnmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
             vvinit(j) = sgn*groundstate(izero)%vec(m)             !build the cdg_up|gs> state
          endif
       enddo
       !normalize the cdg_up|gs> state
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !get cdg_up-sector Hamiltonian
       allocate(H0(jdg0,jdg0))
       call imp_geth(jsect0,H0)
       !Tri-diagonalize w/ Lanczos the resolvant:
       allocate(vout(jdg0))
       vout= 0.d0 ; alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
       call plain_lanczos_step(vvinit,vout,alfa_,beta_,nitermax,nlanc)
       call add_to_lanczos_gf(norm0,groundstate(izero)%egs,nlanc,alfa_,beta_,1,1)
       deallocate(H0,vout,vvinit)


       !REMOVE ONE PARTICLE UP:
       !get c_up sector informations:
       jsect0 = getCUPloop(isect0);if(jsect0==0)cycle
       jdg0   = deg(jsect0)
       jn0    = getin(jsect0)
       js0    = getis(jsect0)
       print*,'sector C_up|gs>',jn0,js0,jdg0
       !allocate c_up|gs> vector:
       allocate(vvinit(jdg0));vvinit=0.d0
       !build c_up|gs> vector:
       do m=1,idg0              !loop over |gs> components m
          i=nmap(isect0,m)      !map m to full-Hilbert space state i
          call bdecomp(i,ib)    !decompose i into number representation ib=|1/0,1/0,1/0...>
          if(ib(1)==1)then      !if impurity is empty: proceed
             call c(1,i,r);sgn=dfloat(r)/dfloat(abs(r));r=abs(r) !apply c_up (1), bring from i to r
             j=invnmap(jsect0,r)                                   !map r back to c_up sector jsect0
             vvinit(j) = sgn*groundstate(izero)%vec(m)             !build the c_up|gs> state
          endif
       enddo
       !normalize the c_up|gs> state
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !get cdg_up-sector Hamiltonian
       allocate(H0(jdg0,jdg0))
       call imp_geth(jsect0,H0)
       !Tri-diagonalize w/ Lanczos the resolvant:
       allocate(vout(jdg0))
       vout= 0.d0 ; alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
       call plain_lanczos_step(vvinit,vout,alfa_,beta_,nitermax,nlanc)
       call add_to_lanczos_gf(norm0,groundstate(izero)%egs,nlanc,alfa_,beta_,-1,1)
       deallocate(H0,vout,vvinit)

    enddo

    nsimp  = nsimp/factor
    nupimp = nupimp/factor
    ndwimp = ndwimp/factor
    dimp   = dimp/factor
    magimp = magimp/factor
    m2imp  = m2imp/factor
    lanc_giw=lanc_giw/factor
    lanc_gwr=lanc_gwr/factor

    print*,nsimp,dimp

    do i=1,Nw
       write(200,*)wr(i),-dimag(lanc_gwr(1,i))/pi
    enddo
    rewind(200)
    do i=1,NL
       write(300,*)wm(i),dimag(lanc_giw(1,i))
    enddo
    rewind(300)
    deallocate(wm,wr)
    deallocate(lanc_giw,lanc_gwr)
  end subroutine lanc_getgf



  subroutine add_to_lanczos_gf(vnorm,emin,nlanc,alanc,blanc,isign,ispin)
    real(8),dimension(:)                         :: alanc,blanc 
    real(8),dimension(size(alanc),size(alanc))   :: Z
    real(8),dimension(size(alanc))               :: diag,subdiag
    real(8) :: vnorm,emin
    integer :: i,j,isign,ispin,ierr,Nlanc
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    do i=1,NL
       do j=1,nlanc
          lanc_giw(ispin,i)=lanc_giw(ispin,i) + vnorm**2*Z(1,j)**2/(xi*wm(i) - isign*(diag(j)-emin))
       enddo
    enddo
    do i=1,Nw
       do j=1,nlanc
          lanc_gwr(ispin,i)=lanc_gwr(ispin,i) + vnorm**2*Z(1,j)**2/(dcmplx(wr(i),eps)-isign*(diag(j)-emin))
       enddo
    enddo
  end subroutine add_to_lanczos_gf
  !#LANCZOS
  !



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
