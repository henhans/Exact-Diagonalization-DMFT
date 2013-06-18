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

    !Zeta:
    zeta_function=real(numzero,8)
    call start_timer
    do izero=1,numzero 
       !get gs-sector information
       isect0 =  es_get_sector(groundstate,izero)
       gsvec  => es_get_vector(groundstate,izero)
       egs    =  es_get_energy(groundstate,izero)
       norm0=sqrt(dot_product(gsvec,gsvec))
       if(abs(norm0-1.d0)>1.d-9)call warning("GS"//txtfy(izero)//"is not normalized:"//txtfy(norm0))
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             call msg("Evaluating diagonal G_imp_Orb"//trim(txtfy(iorb))//trim(txtfy(iorb))//&
                  "_Spin"//trim(txtfy(ispin))//"_Sect0"//trim(txtfy(izero)),unit=LOGfile)
             call lanc_ed_buildgf(isect0,iorb,ispin)
          enddo
       enddo
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                call lanc_ed_buildgf_mix(isect0,iorb,jorb,ispin)
             enddo
          enddo
       enddo
       !
    enddo
    call stop_timer
    Giw=Giw/zeta_function
    Gwr=Gwr/zeta_function
    !Print convenience impurity functions:
    call print_imp_gf
    deallocate(wm,wr)
    deallocate(gsvec)
  end subroutine lanc_ed_getgf




  subroutine lanc_ed_buildgf(isect0,iorb,ispin)
    integer             :: iorb,ispin,isite,isect0
    integer             :: nlanc,idim0,jsect0
    integer             :: jup0,jdw0,jdim0
    integer             :: ib(Ntot)
    integer             :: m,i,j,r
    real(8)             :: norm0,sgn
    real(8),allocatable :: vvinit(:),alfa_(:),beta_(:)
    integer             :: Nitermax
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
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,1,iorb,iorb,ispin)
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
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,-1,iorb,iorb,ispin)
       deallocate(vvinit)
       ! !##IF SPARSE_MATRIX:
       call sp_delete_matrix(spH0)
    endif
    deallocate(alfa_,beta_)
  end subroutine lanc_ed_buildgf


  subroutine lanc_ed_buildgf_mix(isect0,iorb,jorb,ispin)
    integer             :: iorb,jorb,ispin,isite,jsite,isect0
    integer             :: nlanc,idim0,jsect0
    integer             :: jup0,jdw0,jdim0
    integer             :: ib(Ntot)
    integer             :: m,i,j,r
    real(8)             :: norm0,sgn
    real(8),allocatable :: v1(:),v2(:),vvinit(:),alfa_(:),beta_(:)
    integer             :: Nitermax
    Nitermax=nGFitermax
    allocate(alfa_(Nitermax),beta_(Nitermax))
    !Get site index of the iorb/jorb-impurity:
    isite=impIndex(iorb,ispin)
    jsite=impIndex(jorb,ispin)
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
          if(ib(jsite)==0)then                                     !if impurity is empty: proceed
             call cdg(jsite,i,r)
             sgn=dfloat(r)/dfloat(abs(r));r=abs(r)                 !apply cdg_up (1), bring from i to r
             j=invHmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
             vvinit(j) = vvinit(j) + sgn*gsvec(m)                                !build the cdg_up|gs> state
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
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,1,iorb,jorb,ispin)
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
          if(ib(jsite)==1)then                                     !if impurity is empty: proceed
             call c(jsite,i,r)
             sgn=dfloat(r)/dfloat(abs(r));r=abs(r)                 !apply cdg_up (1), bring from i to r
             j=invHmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
             vvinit(j) = vvinit(j) + sgn*gsvec(m)                                !build the cdg_up|gs> state
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
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,-1,iorb,jorb,ispin)
       deallocate(vvinit)
       ! !##IF SPARSE_MATRIX:
       call sp_delete_matrix(spH0)
    endif
    deallocate(alfa_,beta_)
  end subroutine lanc_ed_buildgf_mix




  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_gf(vnorm,emin,nlanc,alanc,blanc,isign,iorb,jorb,ispin)
    real(8)                                    :: vnorm,emin
    integer                                    :: nlanc
    real(8),dimension(nlanc)                   :: alanc,blanc 
    integer                                    :: isign,iorb,jorb,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8) :: cdummy
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    if(iorb==jorb)then
       do i=1,NL
          cdummy = vnorm**2*sum(Z(1,:)**2/(xi*wm(i)-isign*(diag(:)-emin)))
          Giw(iorb,iorb,ispin,i)=Giw(iorb,iorb,ispin,i)+cdummy          
       enddo
       do i=1,Nw
          cdummy = vnorm**2*sum(Z(1,:)**2/(dcmplx(wr(i),eps)-isign*(diag(:)-emin)))
          Gwr(iorb,iorb,ispin,i)=Gwr(iorb,iorb,ispin,i) + cdummy          
       enddo
    else
       do i=1,NL
          cdummy = vnorm**2*sum(Z(1,:)**2/(xi*wm(i)-isign*(diag(:)-emin)))
          Giw(iorb,jorb,ispin,i) = 0.50d0*(cdummy - Giw(iorb,iorb,ispin,i) - Giw(jorb,jorb,ispin,i))
       enddo
       do i=1,Nw
          cdummy = vnorm**2*sum(Z(1,j)**2/(dcmplx(wr(:),eps)-isign*(diag(:)-emin)))
          Gwr(iorb,jorb,ispin,i) = 0.5d0*(cdummy - Gwr(iorb,iorb,ispin,i) - Gwr(jorb,jorb,ispin,i))
       enddo
    endif
  end subroutine add_to_lanczos_gf




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_gf
    integer                                  :: i,j,ispin,unit(6),iorb,jorb
    complex(8)                               :: iw
    complex(8),dimension(Norb,Norb,Nspin,NL) :: G0iw
    complex(8),dimension(Norb,Norb,Nspin,Nw) :: G0wr
    complex(8),dimension(Norb,Norb)          :: Gfoo
    complex(8)                               :: G0inv(Nspin,NL),G0invr(Nspin,Nw)
    real(8) :: kdelta  
    character(len=20) :: suffix



    !Build the impurity Self-energies:
    do ispin=1,Nspin
       do iorb=1,Norb
          do jorb=1,Norb
             !Get Weiss Fields (from Bath):
             kdelta=0.d0;if(iorb==jorb)kdelta=1.d0
             do i=1,NL
                iw=xi*wm(i)
                G0iw(iorb,jorb,ispin,i)= kdelta*(iw+xmu)-delta_bath(iw,iorb,jorb,ispin)
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                G0wr(iorb,jorb,ispin,i)= kdelta*(iw+xmu)-delta_bath(iw,iorb,jorb,ispin)
             enddo
          enddo
       enddo
       do i=1,NL
          Gfoo = Giw(:,:,ispin,i)
          call matrix_inverse(Gfoo)
          Siw(:,:,ispin,i) = G0iw(:,:,ispin,i) - Gfoo(:,:)           
       enddo
       do i=1,Nw
          Gfoo = Gwr(:,:,ispin,i)
          call matrix_inverse(Gfoo)
          Swr(:,:,ispin,i) = G0wr(:,:,ispin,i) - Gfoo(:,:)
       enddo
    enddo


    !Get the Weiss-Field of the Anderson problem: by inverting calG0^-1 matrix
    do ispin=1,Nspin
       do i=1,NL
          call matrix_inverse(G0iw(:,:,ispin,i))
       enddo
       do i=1,Nw
          call matrix_inverse(G0wr(:,:,ispin,i))
       enddo
    enddo


    !Print the impurity functions:
    do iorb=1,Norb
       do jorb=1,Norb
          suffix="_orb"//reg(txtfy(iorb))//reg(txtfy(jorb))//".ed"

          unit(1)=free_unit()
          open(unit(1),file=trim(GMfile)//reg(suffix))

          unit(2)=free_unit()
          open(unit(2),file="impG0_iw"//reg(suffix))

          unit(3)=free_unit()
          open(unit(3),file="impSigma_iw"//reg(suffix))

          unit(4)=free_unit()
          open(unit(4),file=trim(GRfile)//reg(suffix))

          unit(5)=free_unit()
          open(unit(5),file="impG0_realw"//reg(suffix))

          unit(6)=free_unit()
          open(unit(6),file="impSigma_realw"//reg(suffix))

          do i=1,NL
             write(unit(1),"(F18.10,6(F18.10))")wm(i),&
                  (dimag(Giw(iorb,jorb,ispin,i)),dreal(Giw(iorb,jorb,ispin,i)),ispin=1,Nspin)
             write(unit(2),"(F18.10,6(F18.10))")wm(i),&
                  (dimag(G0iw(iorb,jorb,ispin,i)),dreal(G0iw(iorb,jorb,ispin,i)),ispin=1,Nspin)
             write(unit(3),"(F18.10,6(F18.10))")wm(i),&
                  (dimag(Siw(iorb,jorb,ispin,i)),dreal(Siw(iorb,jorb,ispin,i)),ispin=1,Nspin)
          enddo
          do i=1,Nw
             write(unit(4),"(F18.10,6(F18.10))")wr(i),&
                  (dimag(Gwr(iorb,jorb,ispin,i)),dreal(Gwr(iorb,jorb,ispin,i)),ispin=1,Nspin)
             write(unit(5),"(F18.10,6(F18.10))")wr(i),&
                  (dimag(G0wr(iorb,jorb,ispin,i)),dreal(G0wr(iorb,jorb,ispin,i)),ispin=1,Nspin)
             write(unit(6),"(F18.10,6(F18.10))")wr(i),&
                  (dimag(Swr(iorb,jorb,ispin,i)),dreal(Swr(iorb,jorb,ispin,i)),ispin=1,Nspin)
          enddo
          do i=1,6
             close(unit(i))
          enddo
       enddo
    enddo
  end subroutine print_imp_gf





end MODULE ED_GETGF
