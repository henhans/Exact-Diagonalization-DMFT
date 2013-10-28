!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  !
  implicit none
  private 

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer :: state_vec
  real(8)                      :: state_e

  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)    :: Chitau
  real(8),allocatable,dimension(:)      :: Chitautot
  complex(8),allocatable,dimension(:,:) :: Chiw,Chiiw
  complex(8),allocatable,dimension(:)   :: Chiwtot,Chiiwtot

  public :: lanc_ed_getgf
  public :: full_ed_getgf
  public :: lanc_ed_getchi
  public :: full_ed_getchi

contains


  !####################################################################
  !                    lanczos DIAGONALIZATION 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions using Lanczos algorithm
  !+------------------------------------------------------------------+
  subroutine lanc_ed_getgf()
    integer :: izero,iorb,jorb,ispin,i
    integer :: isect0,nup0,ndw0,numstates
    real(8) :: norm0
    call allocate_grids
    impGmats=zero
    impGreal=zero
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"Evaluating Green's functions"
#ifdef _MPI
    endif
#endif

    numstates=numgs
    if(finiteT)numstates=state_list%size
#ifdef _MPI
    call lanczos_plain_set_htimesv_d(HtimesV)
#else
    call lanczos_plain_set_htimesv_d(HtimesV)
#endif
    do ispin=1,Nspin
       do iorb=1,Norb
#ifdef _MPI
          if(mpiID==0)then
#endif
             write(LOGfile,"(A)")"Evaluating G_imp_Orb"//&
                  reg(txtfy(iorb))//"_Spin"//reg(txtfy(ispin))
             call start_progress
#ifdef _MPI
          endif
#endif

          do izero=1,numstates
#ifdef _MPI
             if(mpiID==0)then
#endif
                call progress(izero,numstates)
#ifdef _MPI
             endif
#endif
             isect0 =  es_return_sector(state_list,izero)
             state_e    =  es_return_energy(state_list,izero)
             state_vec  => es_return_vector(state_list,izero)
             nup0  = getnup(isect0)
             ndw0  = getndw(isect0)
             norm0=sqrt(dot_product(state_vec,state_vec))
             if(abs(norm0-1.d0)>1.d-9)then
                write(LOGfile,*) "GS"//reg(txtfy(izero))//&
                     "is not normalized:"//reg(txtfy(norm0))
                stop
             endif
             call lanc_ed_buildgf(isect0,iorb,ispin)
             nullify(state_vec)
          enddo
#ifdef _MPI
          if(mpiID==0)then
#endif
             call stop_progress
#ifdef _MPI
          endif
#endif
       enddo
    enddo
    call lanczos_plain_delete_htimesv_d
    impGmats=impGmats/zeta_function
    impGreal=impGreal/zeta_function
    !Print impurity functions:
#ifdef _MPI
    if(mpiID==0)then
#endif
       call print_imp_gf
#ifdef _MPI
    endif
#endif
    deallocate(wm,wr,tau,vm)
  end subroutine lanc_ed_getgf



  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions using Lanczos algorithm
  !+------------------------------------------------------------------+
  subroutine lanc_ed_getchi()
    integer :: izero,iorb,jorb,ispin,i
    integer :: isect0,nup0,ndw0,numstates
    real(8) :: norm0
    allocate(wm(NL))
    wm     = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(vm(0:NL))          !bosonic frequencies
    do i=0,NL
       vm(i) = pi/beta*2.d0*real(i,8)
    enddo
    allocate(wr(Nw))
    wr     = linspace(wini,wfin,Nw)
    allocate(tau(0:Ltau))
    tau(0:)= linspace(0.d0,beta,Ltau+1)
    allocate(Chitau(Norb,0:Ltau),Chiw(Norb,Nw),Chiiw(Norb,0:NL))
    Chitau=0.d0
    Chiw=zero
    Chiiw=zero
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"Evaluating Susceptibility:"
#ifdef _MPI
    endif
#endif

    numstates=numgs
    if(finiteT)numstates=state_list%size
#ifdef _MPI
    call lanczos_plain_set_htimesv_d(HtimesV)
#else
    call lanczos_plain_set_htimesv_d(HtimesV)
#endif
    do iorb=1,Norb
#ifdef _MPI
       if(mpiID==0)then
#endif
          write(LOGfile,"(A)")"Evaluating Chi_Orb"//reg(txtfy(iorb))
          call start_progress
#ifdef _MPI
       endif
#endif
       do izero=1,numstates
          call progress(izero,numstates)
          isect0 =  es_return_sector(state_list,izero)
          state_e    =  es_return_energy(state_list,izero)
          state_vec  => es_return_vector(state_list,izero)
          nup0  = getnup(isect0)
          ndw0  = getndw(isect0)
          norm0=sqrt(dot_product(state_vec,state_vec))
          if(abs(norm0-1.d0)>1.d-9)then
             write(LOGfile,*) "GS"//reg(txtfy(izero))//&
                  "is not normalized:"//reg(txtfy(norm0))
             stop
          endif
          call lanc_ed_buildchi(isect0,iorb)
          nullify(state_vec)
       enddo
#ifdef _MPI
       if(mpiID==0)then
#endif
          call stop_progress
#ifdef _MPI
       endif
#endif
    enddo
    call lanczos_plain_delete_htimesv_d
    Chitau=Chitau/zeta_function
    Chiw=Chiw/zeta_function
    Chiiw=Chiiw/zeta_function
#ifdef _MPI
    if(mpiID==0)then
#endif
       call print_imp_chi()
#ifdef _MPI
    endif
#endif
    deallocate(Chitau,Chiw,Chiiw)
    deallocate(wm,wr,tau,vm)
  end subroutine lanc_ed_getchi




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_buildgf(isect0,iorb,ispin,iverbose)
    integer             :: iorb,ispin,isite,isect0
    integer             :: nlanc,idim0,jsect0
    integer             :: jup0,jdw0,jdim0
    integer             :: ib(Ntot)
    integer             :: m,i,j,r
    real(8)             :: norm0,sgn
    real(8),allocatable :: vvinit(:),alfa_(:),beta_(:)
    integer             :: Nitermax
    logical,optional    :: iverbose
    logical             :: iverbose_
    integer,allocatable,dimension(:)     :: HImap,HJmap    !map of the Sector S to Hilbert space H
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    !
    Nitermax=lanc_nGFiter
    allocate(alfa_(Nitermax),beta_(Nitermax))
    !Get site index of the iorb-impurity:
    isite=impIndex(iorb,ispin)
    idim0  = getdim(isect0)
    !allocate map from isect0 to HS
    allocate(HImap(idim0))
    call build_sector(isect0,HImap)
    !
    !ADD ONE PARTICLE:
    jsect0 = getCDGsector(ispin,isect0)
    if(jsect0/=0)then 
       jdim0  = getdim(jsect0)
       jup0   = getnup(jsect0)
       jdw0   = getndw(jsect0)
#ifdef _MPI
       if(mpiID==0)then
#endif
          if(iverbose_)write(LOGfile,"(A,2I3,I15)")'add particle:',&
               jup0,jdw0,jdim0
#ifdef _MPI
       endif
#endif
       allocate(HJmap(jdim0))
       call build_sector(jsect0,HJmap)
       !
       allocate(vvinit(jdim0))
       vvinit=0.d0
       do m=1,idim0                     !loop over |gs> components m
          i=HImap(m)                    !map m to Hilbert space state i
          call bdecomp(i,ib)            !i into binary representation
          if(ib(isite)==0)then          !if impurity is empty: proceed
             call cdg(isite,i,r,sgn)
             !sgn=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=binary_search(HJmap,r)      !map r back to  jsect0
             vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
          endif
       enddo
       deallocate(HJmap)
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
       call setup_Hv_sector(jsect0)
       call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax)
       call delete_Hv_sector()
       call add_to_lanczos_gf(norm0,state_e,nitermax,alfa_,beta_,1,iorb,ispin)
       deallocate(vvinit)
       if(spH0%status)call sp_delete_matrix(spH0)
    endif
    !
    !
    !REMOVE ONE PARTICLE:
    jsect0 = getCsector(ispin,isect0)
    if(jsect0/=0)then
       jdim0  = getdim(jsect0)
       jup0    = getnup(jsect0)
       jdw0    = getndw(jsect0)
#ifdef _MPI
       if(mpiID==0)then
#endif      
          if(iverbose_)write(LOGfile,"(A,2I3,I15)")'del particle:',&
               jup0,jdw0,jdim0
#ifdef _MPI
       endif
#endif
       allocate(HJmap(jdim0))
       call build_sector(jsect0,HJmap)
       !
       allocate(vvinit(jdim0)) 
       vvinit=0.d0
       do m=1,idim0
          i=HImap(m)
          call bdecomp(i,ib)
          if(ib(isite)==1)then
             call c(isite,i,r,sgn)
             !sgn=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=binary_search(HJmap,r)
             vvinit(j) = sgn*state_vec(m)
          endif
       enddo
       deallocate(HJmap)
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       alfa_=0.d0 ; beta_=0.d0
       call setup_Hv_sector(jsect0)
       call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax)
       call delete_Hv_sector()
       call add_to_lanczos_gf(norm0,state_e,nitermax,alfa_,beta_,-1,iorb,ispin)
       deallocate(vvinit)
       if(spH0%status)call sp_delete_matrix(spH0)
    endif
    deallocate(alfa_,beta_,HImap)
  end subroutine lanc_ed_buildgf



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_buildchi(isect0,iorb,iverbose)
    integer             :: iorb,isite,isect0
    integer             :: nlanc,idim0
    integer             :: iup0,idw0
    integer             :: ib(Ntot)
    integer             :: m,i,j,r
    real(8)             :: norm0,sgn
    real(8),allocatable :: vvinit(:),alfa_(:),beta_(:)
    integer             :: Nitermax
    logical,optional    :: iverbose
    logical             :: iverbose_
    integer,allocatable,dimension(:)   :: HImap    !map of the Sector S to Hilbert space H
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    Nitermax=lanc_nGFiter
    allocate(alfa_(Nitermax),beta_(Nitermax))
    idim0  = getdim(isect0)
    if(isect0/=0)then 
       iup0   = getnup(isect0)
       idw0   = getndw(isect0)
       !allocate map from isect0 to HS
       allocate(HImap(idim0))
       call build_sector(isect0,HImap)
       allocate(vvinit(idim0))
       vvinit=0.d0
       do m=1,idim0                     !loop over |gs> components m
          i=HImap(m)
          call bdecomp(i,ib)
          sgn = real(ib(iorb),8)-real(ib(iorb+Ns),8)
          vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
       enddo
       deallocate(HImap)
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
       call setup_Hv_sector(isect0)
       call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax)
       call delete_Hv_sector()
       call add_to_lanczos_chi(norm0,state_e,nitermax,alfa_,beta_,iorb)
       deallocate(vvinit)
       if(spH0%status)call sp_delete_matrix(spH0)
    endif
    deallocate(alfa_,beta_)
  end subroutine lanc_ed_buildchi






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_gf(vnorm,Ei,nlanc,alanc,blanc,isign,iorb,ispin)
    real(8)                                    :: vnorm,Ei,Egs,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(nlanc)                   :: alanc,blanc 
    integer                                    :: isign,iorb,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    pesoBZ = 1.d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*vnorm**2*Z(1,j)**2
       do i=1,NL
          iw=xi*wm(i)
          impGmats(ispin,iorb,i)=impGmats(ispin,iorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Nw
          iw=cmplx(wr(i),eps,8)
          impGreal(ispin,iorb,i)=impGreal(ispin,iorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_chi(vnorm,Ei,nlanc,alanc,blanc,iorb)
    real(8)                                    :: vnorm,Ei,Egs,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(nlanc)                   :: alanc,blanc 
    integer                                    :: isign,iorb
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    Egs = state_list%emin       !get the gs energy
    pesoBZ = 1.d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       if(de>0.d0)then
          peso = pesoBZ*vnorm**2*Z(1,j)**2
          if(de>cutoff)chiiw(iorb,0)=chiiw(iorb,0) - peso*(exp(-beta*de)-1.d0)/de
          do i=1,NL
             iw=xi*vm(i)
             chiiw(iorb,i)=chiiw(iorb,i) + peso*(exp(-beta*de)-1.d0)/(iw-de)
          enddo
          !
          do i=1,Nw
             iw=cmplx(wr(i),eps,8)
             chiw(iorb,i)=chiw(iorb,i) + peso*(exp(-beta*de)-1.d0)/(iw-de)
          enddo
          !
          do i=0,Ltau 
             chitau(iorb,i)=chitau(iorb,i) + peso*exp(-tau(i)*de)
          enddo
       endif
    enddo
  end subroutine add_to_lanczos_chi







  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : Evalute green's function and suscpetibility
  ! using complete spectral decomposition 
  !+------------------------------------------------------------------+
  subroutine full_ed_getgf()
    integer :: iorb,jorb,ispin
    write(LOGfile,"(A)")"Evaluating Green's functions"
    call allocate_grids
    impGmats   =zero
    impGreal   =zero
#ifdef _MPI
    if(mpiID==0)then
#endif
       call start_timer
#ifdef _MPI
    endif
#endif
    do ispin=1,Nspin
       do iorb=1,Norb
          call full_ed_buildgf(iorb,ispin)
       enddo
    enddo
#ifdef _MPI
    if(mpiID==0)then
#endif
       call print_imp_gf
       call stop_timer
#ifdef _MPI
    endif
#endif
    deallocate(wm,tau,wr,vm)
  end subroutine full_ed_getgf



  !+------------------------------------------------------------------+
  !PURPOSE  :  Evalue the spectral sum for impurity Green's function. 
  !<i|C^+|j>=<in,is,idim|C^+|jn,js,jdim>=C^+_{ij} |
  !+------------------------------------------------------------------+
  subroutine full_ed_buildgf(iorb,ispin)
    integer                 :: iorb,ispin,isite,jsite,nsite
    real(8)                 :: cdgmat,cc
    integer,dimension(Ntot) :: ib
    integer                 :: i,j,k,r,ll,m,in,is
    integer                 :: idim,jdim,isector,jsector,ia
    real(8)                 :: Ei,Ej,matcdg
    real(8)                 :: expterm,peso,de,w0,it,chij1
    complex(8)              :: iw
    integer,allocatable,dimension(:)     :: HJmap,HImap    !map of the Sector S to Hilbert space H


    nsite=1
    isite=impIndex(iorb,ispin)
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"Evaluating G_imp_Orb"//reg(txtfy(iorb))//"_Spin"//reg(txtfy(ispin))
       call start_progress(LOGfile)
#ifdef _MPI
    endif
#endif

    do isector=1,Nsect
       jsector=getCsector(1,isector);if(jsector==0)cycle
#ifdef _MPI
       if(mpiID==0)then
#endif
          call progress(isector,Nsect)
#ifdef _MPI
       endif
#endif
       idim=getdim(isector)     !i-th sector dimension
       jdim=getdim(jsector)     !j-th sector dimension
       allocate(HImap(idim))
       allocate(HJmap(jdim))
       call build_sector(isector,HImap)
       call build_sector(jsector,HJmap)
       do i=1,idim          !loop over the states in the i-th sect.
          do j=1,jdim       !loop over the states in the j-th sect.
             cdgmat=0.d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             do ll=1,jdim              !loop over the component of |j> (IN state!)
                !m=Hmap(jsector)%map(ll)!map from IN state (j) to full Hilbert space
                m = HJmap(ll)
                call bdecomp(m,ib)
                if(ib(isite) == 0)then
                   call cdg(isite,m,k,cc)!;cc=dble(k)/dble(abs(k));k=abs(k)
                   r = binary_search(HImap,k)
                   cdgmat=cdgmat+espace(isector)%M(r,i)*cc*espace(jsector)%M(ll,j)
                endif
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             peso=expterm/zeta_function
             matcdg=peso*cdgmat**2
             !build Matsubara GF
             do m=1,NL
                iw=xi*wm(m)
                impGmats(ispin,iorb,m)=impGmats(ispin,iorb,m)+matcdg/(iw+de)
             enddo
             !build Real-freq. GF
             do m=1,Nw 
                w0=wr(m);iw=cmplx(w0,eps)
                impGreal(ispin,iorb,m)=impGreal(ispin,iorb,m)+matcdg/(iw+de)
             enddo
          enddo
       enddo
       deallocate(HImap,HJmap)
    enddo
#ifdef _MPI
    if(mpiID==0)then
#endif
       call stop_progress
#ifdef _MPI
    endif
#endif
  end subroutine full_ed_buildgf


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine full_ed_getchi()
    real(8)                               :: cdgmat(Norb),chij(Norb),chitot,spin,spintot
    integer,dimension(N)                  :: ib(Ntot)
    integer                               :: i,j,k,r,ll,m,in,is,ispin,iorb,jorb,isector
    integer                               :: idim,ia,unit(3)
    real(8)                               :: Ei,Ej,cc,peso(Norb),pesotot
    real(8)                               :: expterm,de,w0,it
    complex(8)                            :: iw
    integer,allocatable,dimension(:)      :: HImap    !map of the Sector S to Hilbert space H
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"Evaluating Suceptibility:"
#ifdef _MPI
    endif
#endif

    call allocate_grids
    allocate(Chitau(Norb,0:Ltau),Chiw(Norb,Nw),Chiiw(Norb,0:NL))
    allocate(Chitautot(0:Ltau),Chiwtot(Nw),Chiiwtot(0:NL))
    Chitau=0.d0
    Chiw=zero
    Chiiw=zero
    Chitautot=0.d0
    Chiwtot=zero
    Chiiwtot=zero
    !Spin susceptibility \X(tau). |<i|S_z|j>|^2
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"Evaluating Chi_Sz"
       call start_progress(LOGfile)
#ifdef _MPI
    endif
#endif

    do isector=1,Nsect !loop over <i| total particle number
#ifdef _MPI
       if(mpiID==0)then
#endif
          call progress(isector,Nsect)
#ifdef _MPI
       endif
#endif
       idim=getdim(isector)
       allocate(HImap(idim))
       call build_sector(isector,HImap)
       do i=1,idim 
          do j=1,idim
             chij=0.d0;chitot=0.d0
             expterm=exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim 
                ia=HImap(ll)
                call bdecomp(ia,ib)
                spintot=0.d0
                do iorb=1,Norb
                   spin=real(ib(iorb),8)-real(ib(iorb+Ns),8) !nup - ndw
                   chij(iorb)=chij(iorb)+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
                enddo
                spintot=spintot+real(ib(iorb),8)-real(ib(iorb+Ns),8)
                chitot=chitot+espace(isector)%M(ll,i)*spintot*espace(isector)%M(ll,j)
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso=chij**2/zeta_function
             pesotot=chitot**2/zeta_function
             do iorb=1,Norb
                !Matsubara (bosonic) frequency
                if(de>cutoff)chiiw(iorb,0)=chiiw(iorb,0)-peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
                do m=1,NL
                   iw=xi*vm(m)
                   chiiw(iorb,m)=chiiw(iorb,m)+peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
                enddo
                !Imaginary time:
                do m=0,Ltau 
                   it=tau(m)
                   chitau(iorb,m)=chitau(iorb,m) + exp(-it*espace(isector)%e(i))*&
                        exp(-(beta-it)*espace(isector)%e(j))*peso(iorb)
                   ! chitau(iorb,m)=chitau(iorb,m) + peso(iorb)*exp(-beta*Ej)*exp(-it*de)
                enddo
                !Real-frequency
                do m=1,Nw
                   w0=wr(m);iw=cmplx(w0,eps,8)
                   !Time-ordered
                   ! chiw(iorb,m)=chiw(iorb,m)-exp(-beta*espace(isector)%e(j))*&
                   !      (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso(iorb)
                   !Retarded = Commutator = response function
                   chiw(iorb,m)=chiw(iorb,m)+&
                        exp(-beta*Ej)*peso(iorb)*(exp(-beta*de)-1.d0)/(iw-de)
                enddo
             enddo
             !Matsubara (bosonic) frequency
             if(de>cutoff)chiiwtot(0)=chiiwtot(0)-pesotot*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
             do m=1,NL
                iw=xi*vm(m)
                chiiwtot(m)=chiiwtot(m)+pesotot*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
             enddo
             do m=0,Ltau 
                it=tau(m)
                chitautot(m)=chitautot(m) + exp(-it*espace(isector)%e(i))*&
                     exp(-(beta-it)*espace(isector)%e(j))*pesotot
                !chitautot(m)=chitautot(m) + pesotot*exp(-beta*Ej)*exp(-it*de)
             enddo
             do m=1,Nw
                w0=wr(m);iw=cmplx(w0,eps,8)
                ! chiwtot(m)=chiwtot(m)-exp(-beta*espace(isector)%e(j))*&
                !      (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*pesotot
                chiwtot(m)=chiwtot(m)+&
                     exp(-beta*Ej)*pesotot*(exp(-beta*de)-1.d0)/(iw-de)
             enddo
          enddo
       enddo
    enddo
#ifdef _MPI
    if(mpiID==0)then
#endif
       call stop_progress
       call print_imp_chi()
#ifdef _MPI
    endif
#endif
    deallocate(Chitau,Chiw,Chiiw,Chitautot,Chiwtot,Chiiwtot)
    deallocate(wm,tau,wr,vm)
  end subroutine full_ed_getchi







  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    allocate(wm(NL))
    wm     = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(vm(0:NL))          !bosonic frequencies
    do i=0,NL
       vm(i) = pi/beta*2.d0*real(i,8)
    enddo
    allocate(wr(Nw))
    wr     = linspace(wini,wfin,Nw)
    allocate(tau(0:Ltau))
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_gf
    integer                             :: i,j,ispin,unit(6),iorb,jorb
    complex(8)                          :: iw
    complex(8),dimension(Nspin,Norb,NL) :: G0iw
    complex(8),dimension(Nspin,Norb,Nw) :: G0wr
    character(len=20)                   :: prefix
    write(LOGfile,"(A)")"Printing the impurity GF"
    do ispin=1,Nspin
       do iorb=1,Norb
          !Get Weiss Fields (from Bath):
          do i=1,NL
             iw=xi*wm(i)
             G0iw(ispin,iorb,i)= (iw+xmu)-delta_bath(ispin,iorb,iw)
             impSmats(ispin,iorb,i) = G0iw(ispin,iorb,i) - one/impGmats(ispin,iorb,i)
             G0iw(ispin,iorb,i)= delta_bath(ispin,iorb,iw)
          enddo
          do i=1,Nw
             iw=cmplx(wr(i),eps)
             G0wr(ispin,iorb,i)= (iw+xmu)-delta_bath(ispin,iorb,iw)
             impSreal(ispin,iorb,i) = G0wr(ispin,iorb,i) - one/impGreal(ispin,iorb,i)
             G0wr(ispin,iorb,i)= delta_bath(ispin,iorb,iw)
          enddo
       enddo
    enddo
    !Print the impurity functions:
    do iorb=1,Norb
       prefix="_orb"//reg(txtfy(iorb))
       unit(1)=free_unit()
       open(unit(1),file=reg(GFfile)//reg(prefix)//"_iw.ed")
       unit(2)=free_unit()
       open(unit(2),file=reg(GFfile)//reg(prefix)//"_realw.ed")
       unit(3)=free_unit()
       open(unit(3),file="impSigma"//reg(prefix)//"_iw.ed")
       unit(4)=free_unit()
       open(unit(4),file="impSigma"//reg(prefix)//"_realw.ed")
       unit(5)=free_unit()
       open(unit(5),file="impDelta"//reg(prefix)//"_iw.ed")
       unit(6)=free_unit()
       open(unit(6),file="impDelta"//reg(prefix)//"_realw.ed")
       do i=1,NL
          write(unit(1),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impGmats(ispin,iorb,i)),dreal(impGmats(ispin,iorb,i)),ispin=1,Nspin)
          write(unit(3),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impSmats(ispin,iorb,i)),dreal(impSmats(ispin,iorb,i)),ispin=1,Nspin)
          write(unit(5),"(F26.15,6(F26.15))")wm(i),&
               (dimag(G0iw(ispin,iorb,i)),dreal(G0iw(ispin,iorb,i)),ispin=1,Nspin)
       enddo
       do i=1,Nw
          write(unit(2),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impGreal(ispin,iorb,i)),dreal(impGreal(ispin,iorb,i)),ispin=1,Nspin)
          write(unit(4),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impSreal(ispin,iorb,i)),dreal(impSreal(ispin,iorb,i)),ispin=1,Nspin)
          write(unit(6),"(F26.15,6(F26.15))")wr(i),&
               (dimag(G0wr(ispin,iorb,i)),dreal(G0wr(ispin,iorb,i)),ispin=1,Nspin)
       enddo
       do i=1,6
          close(unit(i))
       enddo
    enddo
  end subroutine print_imp_gf



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_chi
    integer                               :: i,j,iorb
    integer                               :: unit(3)
    do iorb=1,Norb
       unit(1)=free_unit()
       open(unit(1),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_tau.ed")
       unit(2)=free_unit()
       open(unit(2),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_realw.ed")
       unit(3)=free_unit()
       open(unit(3),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_iw.ed")
       do i=0,Ltau/2
          write(unit(1),*)tau(i),chitau(iorb,i)
       enddo
       do i=1,Nw
          if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(iorb,i)),dreal(chiw(iorb,i))
       enddo
       do i=0,NL
          write(unit(3),*)vm(i),dimag(chiiw(iorb,i)),dreal(chiiw(iorb,i))
       enddo
       close(unit(1))
       close(unit(2))
       close(unit(3))
    enddo
    ! unit(1)=free_unit()
    ! open(unit(1),file=trim(CHIfile)//"_tot_tau.ed")
    ! unit(2)=free_unit()
    ! open(unit(2),file=trim(CHIfile)//"_tot_realw.ed")
    ! unit(3)=free_unit()
    ! open(unit(3),file=trim(CHIfile)//"_tot_iw.ed")
    ! do i=0,Ltau
    !    write(unit(1),*)tau(i),chitautot(i)
    ! enddo
    ! do i=1,Nw
    !    if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiwtot(i)),dreal(chiwtot(i))
    ! enddo
    ! do i=0,NL
    !    write(unit(3),*)vm(i),dimag(chiiwtot(i)),dreal(chiiwtot(i))
    ! enddo
    ! close(unit(1))
    ! close(unit(2))
    ! close(unit(3))
  end subroutine print_imp_chi

end MODULE ED_GREENS_FUNCTIONS
