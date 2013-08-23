!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GETGF
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_GETH
  !

  implicit none
  private 

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer :: gsvec
  real(8)                      :: egs



  public :: full_ed_getgf
  public :: lanc_ed_getgf
  public :: full_ed_getchi

contains

  subroutine allocate_grids
    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)
    allocate(tau(0:Ltau))
    tau   = linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  include 'include_fulled_getgf.f90'



  !####################################################################
  !                    LANCZOS DIAGONALIZATION (T=0, GS only)
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_getgf()
    integer :: izero,iorb,jorb,ispin
    integer :: isect0
    real(8) :: norm0
    !SET THE LANCZOS H*v method:
    call plain_lanczos_set_htimesv_d(spHtimesV)
    !set grids
    call allocate_grids
    !Set Max GF iterations
    impGmats=zero
    impGreal=zero

    !Zeta:
    zeta_function=real(numzero,8)
    call start_timer
    do izero=1,numzero 
       !get gs-sector information
       isect0 =  es_get_sector(groundstate,izero)
       egs    =  es_get_energy(groundstate,izero)
       gsvec  => es_get_vector(groundstate,izero)
       norm0=sqrt(dot_product(gsvec,gsvec))
       if(abs(norm0-1.d0)>1.d-9)call warning("GS"//reg(txtfy(izero))//"is not normalized:"//txtfy(norm0))
       !
       do ispin=1,Nspin
          do iorb=1,Norb             
             call lanc_ed_buildgf(isect0,iorb,ispin)
          enddo          
       enddo
       !
       nullify(gsvec)
    enddo
    impGmats=impGmats/zeta_function
    impGreal=impGreal/zeta_function
    !Print convenience impurity functions:
    call print_imp_gf
    call stop_timer
    deallocate(wm,wr,tau)
  end subroutine lanc_ed_getgf



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_buildgf(isect0,iorb,ispin)
    integer             :: iorb,ispin,isite,isect0
    integer             :: nlanc,idim0,jsect0
    integer             :: jup0,jdw0,jdim0
    integer             :: ib(Ntot)
    integer             :: m,i,j,r
    real(8)             :: norm0,sgn
    real(8),allocatable :: vvinit(:),alfa_(:),beta_(:)
    integer             :: Nitermax
    call msg("Evaluating diagonal G_imp_Orb"//reg(txtfy(iorb))//reg(txtfy(iorb))//&
         "_Spin"//reg(txtfy(ispin)),unit=LOGfile)
    Nitermax=nGFitermax
    allocate(alfa_(Nitermax),beta_(Nitermax))
    !Get site index of the iorb-impurity:
    isite=impIndex(iorb,ispin)
    !Get dimension of the gs-sector isect0
    idim0  = getdim(isect0)
    !Get the +up particle sector information:
    jsect0 = getCDGsector(ispin,isect0)
    if(jsect0/=0)then 
       jdim0  = getdim(jsect0)
       jup0   = getnup(jsect0)
       jdw0   = getndw(jsect0)
       write(*,"(A,2I3,I15)")'GetGF sector:',jup0,jdw0,jdim0
       allocate(vvinit(jdim0));vvinit=0.d0
       do m=1,idim0                                                !loop over |gs> components m
          i=Hmap(isect0)%map(m)                                    !map m to full-Hilbert space state i
          call bdecomp(i,ib)                                       !decompose i into binary representation
          if(ib(isite)==0)then                                     !if impurity is empty: proceed
             call cdg(isite,i,r)
             sgn=dfloat(r)/dfloat(abs(r));r=abs(r)                 !apply cdg_up (1), bring from i to r
             j=invHmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
             vvinit(j) = sgn*gsvec(m)                                !build the cdg_up|gs> state
          endif
       enddo
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       ! !##IF SPARSE_MATRIX:
       call sp_init_matrix(spH0,jdim0)
       call lanc_ed_geth(jsect0)
       ! !##ELSE DIRECT H*V PRODUCT:
       ! call set_Hsector(jsect0)
       alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
       call plain_lanczos_tridiag(vvinit,alfa_,beta_,nitermax)
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,1,iorb,ispin)
       deallocate(vvinit)
       ! !##IF SPARSE_MATRIX:
       call sp_delete_matrix(spH0)
    endif
    !
    !REMOVE ONE PARTICLE UP:
    jsect0 = getCsector(ispin,isect0)
    if(jsect0/=0)then
       jdim0  = getdim(jsect0)
       jup0    = getnup(jsect0)
       jdw0    = getndw(jsect0)
       write(*,"(A,2I3,I15)")'GetGF: sector:',jup0,jdw0,jdim0
       allocate(vvinit(jdim0)) ; vvinit=0.d0
       do m=1,idim0                                                !loop over |gs> components m
          i=Hmap(isect0)%map(m)                                    !map m to full-Hilbert space state i
          call bdecomp(i,ib)                                       !decompose i into binary representation
          if(ib(isite)==1)then                                     !if impurity is empty: proceed
             call c(isite,i,r)
             sgn=dfloat(r)/dfloat(abs(r));r=abs(r)                 !apply cdg_up (1), bring from i to r
             j=invHmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
             vvinit(j) = sgn*gsvec(m)                                !build the cdg_up|gs> state
          endif
       enddo
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       call sp_init_matrix(spH0,jdim0)
       call lanc_ed_geth(jsect0)
       ! !##ELSE DIRECT H*V PRODUCT:
       ! call set_Hsector(jsect0)
       alfa_=0.d0 ; beta_=0.d0
       call plain_lanczos_tridiag(vvinit,alfa_,beta_,nitermax)
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,-1,iorb,ispin)
       deallocate(vvinit)
       ! !##IF SPARSE_MATRIX:
       call sp_delete_matrix(spH0)
    endif
    deallocate(alfa_,beta_)
  end subroutine lanc_ed_buildgf








  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_gf(vnorm,emin,nlanc,alanc,blanc,isign,iorb,ispin)
    real(8)                                    :: vnorm,emin
    integer                                    :: nlanc
    real(8),dimension(nlanc)                   :: alanc,blanc 
    integer                                    :: isign,iorb,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: cdummy
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    do i=1,NL
       do j=1,nlanc
          impGmats(ispin,i)=impGmats(ispin,i) + vnorm**2*Z(1,j)**2/(xi*wm(i) - isign*(diag(j)-emin))
       enddo
    enddo
    do i=1,Nw
       do j=1,nlanc
          impGreal(ispin,i)=impGreal(ispin,i) + vnorm**2*Z(1,j)**2/(cmplx(wr(i),eps,8)-isign*(diag(j)-emin))
       enddo
    enddo
  end subroutine add_to_lanczos_gf




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_gf
    integer                      :: i,j,ispin,unit(6)
    complex(8)                   :: iw
    complex(8),dimension(1:2,NL) :: G0iw
    complex(8),dimension(1:2,Nw) :: G0wr
    unit(1)=free_unit()
    open(unit(1),file=trim(GMfile)//".ed")
    unit(2)=free_unit()
    open(unit(2),file="impG0_iw.ed")
    unit(3)=free_unit()
    open(unit(3),file="impSigma_iw.ed")
    unit(4)=free_unit()
    open(unit(4),file=trim(GRfile)//".ed")
    unit(5)=free_unit()
    open(unit(5),file="impG0_realw.ed")
    unit(6)=free_unit()
    open(unit(6),file="impSigma_realw.ed")
    do ispin=1,Nspin
       !Get Weiss Fields (from Bath):
       do i=1,NL
          iw=xi*wm(i)
          G0iw(ispin,i)= iw+xmu-delta_bath(iw,ispin)
       enddo
       do i=1,Nw
          iw=cmplx(wr(i),eps)
          G0wr(ispin,i)= iw+xmu-delta_bath(iw,ispin)
       enddo
       impSmats(ispin,:) = G0iw(ispin,:) - one/impGmats(ispin,:)
       impSreal(ispin,:) = G0wr(ispin,:) - one/impGreal(ispin,:)
       !Print GFs
    enddo
    do i=1,NL
       write(unit(1),"(F18.10,6(F18.10))")wm(i),(dimag(impGmats(ispin,i)),dreal(impGmats(ispin,i)),ispin=1,Nspin)
       write(unit(2),"(F18.10,6(F18.10))")wm(i),(dimag(one/G0iw(ispin,i)),dreal(one/G0iw(ispin,i)),ispin=1,Nspin)
       write(unit(3),"(F18.10,6(F18.10))")wm(i),(dimag(impSmats(ispin,i)),dreal(impSmats(ispin,i)),ispin=1,Nspin)
    enddo
    do i=1,Nw
       write(unit(4),"(F18.10,6(F18.10))")wr(i),(dimag(impGreal(ispin,i)),dreal(impGreal(ispin,i)),ispin=1,Nspin)
       write(unit(5),"(F18.10,6(F18.10))")wr(i),(dimag(one/G0wr(ispin,i)),dreal(one/G0wr(ispin,i)),ispin=1,Nspin)
       write(unit(6),"(F18.10,6(F18.10))")wr(i),(dimag(impSreal(ispin,i)),dreal(impSreal(ispin,i)),ispin=1,Nspin)
    enddo
    do i=1,6
       close(unit(i))
    enddo
  end subroutine print_imp_gf





end MODULE ED_GETGF
