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

  public :: full_ed_getgf
  public :: lanc_ed_getgf
  public :: full_ed_getchi

contains

  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine full_ed_getgf()
    real(8)                      :: cdgmat,matcdg
    integer,dimension(Ntot)      :: ib
    integer                      :: i,j,k,r,ll,m,in,is,ispin,iorb
    integer                      :: idim,jdim,isector,jsector,ia
    real(8)                      :: Ei,Ej
    real(8)                      :: cc,spin1,peso1
    real(8)                      :: expterm,peso,de,w0,it,chij1
    complex(8)                   :: iw
    !----------------------------------------------
    !<i|C^+|j>=<in,is,idim|C^+|jn,js,jdim>=C^+_{ij} |
    !----------------------------------------------
    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(tau(0:Ltau))
    tau   = linspace(0.d0,beta,Ltau+1)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)

    !Initialize some functions
    Giw   =zero
    Gwr   =zero

    do iorb=1,Norb
       do ispin=1,Nspin
          call full_ed_buildgf(iorb,ispin)
       enddo
    enddo

    call print_imp_gf
    deallocate(wm,tau,wr)
  end subroutine full_ed_getgf


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine full_ed_getchi()
    ! USE STATISTICS
    real(8)                  :: cdgmat,matcdg
    integer,dimension(Ntot)  :: ib,ibi,ibj
    integer                  :: i,j,k,r,ll,m,iup,idw,ispin,kk
    integer                  :: idim,jdim,isector,jsector,ia,unit(2)
    real(8)                  :: cc,spin,peso,chij,weigth
    real(8)                  :: de,w0,it,Ei,Ej,expterm,Pchi(Nsect),totPchi
    complex(8)               :: iw
    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(tau(0:Ltau))
    tau   = linspace(0.d0,beta,Ltau+1)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)
    chitau=0.d0
    chiw=zero
    !Spin susceptibility \X(tau). |<i|S_z|j>|^2
    call msg("Evaluating Chi_Sz",unit=LOGfile)
    call start_timer
    do isector=1,Nsect !loop over <i| total particle number
       call eta(isector,lastloop,file="ETA_chi.ed")
       idim=getdim(isector)
       Pchi(isector)=0.d0
       do i=1,idim 
          do j=1,idim
             chij=0.d0
             expterm=exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim 
                ia=Hmap(isector)%map(ll)
                call bdecomp(ia,ib)
                spin=real(ib(1),8)-real(ib(1+Ns),8) !nup - ndw
                chij=chij+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
             enddo
             Pchi(isector)=Pchi(isector)+chij**2
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ej-Ei
             ! if(chij**2 > 1.d-1)then
             !    ia=Hmap(isector)%map(i)
             !    call bdecomp(ia,ibi)
             !    ia=Hmap(isector)%map(j)!Hmap(isector,j)
             !    call bdecomp(ia,ibj)
             !    write(300,"(4I5,2F15.9,2x,14I1,2x,14I1)")in,is,i,j,chij**2,de,&
             !         (ibi(kk),kk=1,Ntot),(ibj(kk),kk=1,Ntot)
             ! endif
             peso=chij**2/zeta_function
             !Imaginary time:
             do m=0,Ltau 
                it=tau(m)
                chitau(m)=chitau(m) + exp(-it*espace(isector)%e(i))*&
                     exp(-(beta-it)*espace(isector)%e(j))*peso
             enddo
             !Real-frequency
             do m=1,Nw
                w0=wr(m);iw=cmplx(w0,eps,8)
                chiw(m)=chiw(m)-exp(-beta*espace(isector)%e(j))*&
                     (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso
             enddo
          enddo
       enddo
    enddo
    call stop_timer
    !rewind(300)

    unit(1)=free_unit()
    open(unit(1),file=trim(CTfile))
    unit(2)=free_unit()
    open(unit(2),file=trim(CWfile))
    do i=0,Ltau
       write(unit(1),*)tau(i),chitau(i)
    enddo
    do i=1,Nw
       if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(i)),dreal(chiw(i))
    enddo
    close(unit(1))
    close(unit(2))
    deallocate(wm,tau,wr)

    totPchi=sum(Pchi)
    Pchi=Pchi/totPchi
    do i=1,Nsect
       write(100,*)getnup(i),getndw(i),Pchi(i)
    enddo
    rewind(100)
  end subroutine full_ed_getchi











  !####################################################################
  !                    LANCZOS DIAGONALIZATION (T=0, GS only)
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_getgf()
    integer                      :: i,izero,isect0,jsect0,m,j
    integer                      :: iorb,ispin
    integer                      :: iup0,idw0,idim0
    integer                      :: jup0,jdw0,jdim0
    real(8)                      :: norm0,sgn,gs,nup,ndw
    integer                      :: ib(Ntot),k,r,Nlanc,Nitermax
    real(8),dimension(:),pointer :: vec
    real(8)                      :: egs
    real(8),allocatable          :: vvinit(:),alfa_(:),beta_(:)
    !
    !SET THE LANCZOS H*v method:
    call plain_lanczos_set_htimesv(spHtimesV)
    !Initialize some functions
    Giw   =zero
    Gwr   =zero
    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)
    !Set Max GF iterations
    Nitermax=nGFitermax
    allocate(alfa_(Nitermax),beta_(Nitermax))
    !Zeta:
    zeta_function=real(numzero,8)

    call start_timer
    do izero=1,numzero 
       !get gs-sector information
       isect0 = es_get_sector(groundstate,izero)
       iup0    = getnup(isect0)
       idw0    = getndw(isect0)
       idim0  = getdim(isect0)
       vec => es_get_vector(groundstate,izero)
       egs =  es_get_energy(groundstate,izero)
       norm0=sqrt(dot_product(vec,vec))
       if(abs(norm0-1.d0)>1.d-9)call warning("GS"//txtfy(izero)//"is not normalized:"//txtfy(norm0))
       !
       do iorb=1,Norb
          do ispin=1,Nspin
             call msg("Evaluating G_imp_Orb"//trim(txtfy(iorb))//&
                  "_Spin"//trim(txtfy(ispin))//&
                  "_Sect0"//trim(txtfy(izero)),unit=LOGfile)
             call lanc_ed_buildgf(iorb,ispin)
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

  contains

    subroutine lanc_ed_buildgf(iorb,ispin)
      integer :: iorb,ispin,isite
      isite=impIndex(iorb,ispin)
      !get sector informations:
      jsect0 = getCDGsector(ispin,isect0)
      if(jsect0/=0)then 
         jdim0  = getdim(jsect0)
         jup0    = getnup(jsect0)
         jdw0    = getndw(jsect0)
         write(*,"(A,2I3,I15)")'GetGF sector:',jup0,jdw0,jdim0
         allocate(vvinit(jdim0));vvinit=0.d0
         do m=1,idim0                                                !loop over |gs> components m
            i=Hmap(isect0)%map(m)                                    !map m to full-Hilbert space state i
            call bdecomp(i,ib)                                       !decompose i into binary representation
            if(ib(isite)==0)then                                     !if impurity is empty: proceed
               call cdg(isite,i,r)
               sgn=dfloat(r)/dfloat(abs(r));r=abs(r)                 !apply cdg_up (1), bring from i to r
               j=invHmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
               vvinit(j) = sgn*vec(m)                                !build the cdg_up|gs> state
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
         call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,1,ispin)
         deallocate(vvinit)
         ! !##IF SPARSE_MATRIX:
         call sp_delete_matrix(spH0)
      endif
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
               vvinit(j) = sgn*vec(m)                                !build the cdg_up|gs> state
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
         call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,-1,ispin)
         deallocate(vvinit)
         ! !##IF SPARSE_MATRIX:
         call sp_delete_matrix(spH0)
      endif
    end subroutine lanc_ed_buildgf


  end subroutine lanc_ed_getgf











  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################

  subroutine full_ed_buildgf(iorb,ispin)
    integer                 :: iorb,ispin,isite
    real(8)                 :: cdgmat,matcdg
    integer,dimension(Ntot) :: ib
    integer                 :: i,j,k,r,ll,m,in,is
    integer                 :: idim,jdim,isector,jsector,ia
    real(8)                 :: Ei,Ej
    real(8)                 :: cc,spin1,peso1
    real(8)                 :: expterm,peso,de,w0,it,chij1
    complex(8)              :: iw
    isite=impIndex(iorb,ispin)
    call msg("Evaluating G_imp_Orb"//trim(txtfy(iorb))//"_Spin"//trim(txtfy(ispin)),unit=LOGfile)
    call start_timer
    do isector=startloop,lastloop
       !if(isector < minCsector(ispin))cycle
       jsector=getCsector(1,isector);if(jsector==0)cycle
       call eta(isector,lastloop,file="ETA_GF_Orb"//trim(txtfy(iorb))//"_Spin"//trim(txtfy(ispin))//".ed")
       idim=getdim(isector)     !i-th sector dimension
       jdim=getdim(jsector)     !j-th sector dimension
       do i=1,idim          !loop over the states in the i-th sect.
          do j=1,jdim       !loop over the states in the j-th sect.
             cdgmat=0.d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm > cutoff)then
                do ll=1,jdim              !loop over the component of |j> (IN state!)
                   m=Hmap(jsector)%map(ll)!map from IN state (j) 2 full Hilbert space
                   call bdecomp(m,ib)
                   if(ib(isite) == 0)then
                      call cdg(isite,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invHmap(isector,k)
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
                   Giw(ispin,m)=Giw(ispin,m)+matcdg/(iw+de)
                enddo
                !build Real-freq. GF
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps)
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg/(iw+de)
                enddo
             endif
          enddo
       enddo
    enddo
    call stop_timer
  end subroutine full_ed_buildgf


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_gf(vnorm,emin,nlanc,alanc,blanc,isign,ispin)
    real(8),dimension(nlanc)                     :: alanc,blanc 
    real(8),dimension(size(alanc),size(alanc))   :: Z
    real(8),dimension(size(alanc))               :: diag,subdiag
    real(8)                                      :: vnorm,emin
    integer                                      :: i,j,isign,ispin,ierr,Nlanc
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    do i=1,NL
       do j=1,nlanc
          Giw(ispin,i)=Giw(ispin,i) + vnorm**2*Z(1,j)**2/(xi*wm(i) - isign*(diag(j)-emin))
       enddo
    enddo
    do i=1,Nw
       do j=1,nlanc
          Gwr(ispin,i)=Gwr(ispin,i) + vnorm**2*Z(1,j)**2/(dcmplx(wr(i),eps)-isign*(diag(j)-emin))
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
    open(unit(1),file=trim(GMfile))
    unit(2)=free_unit()
    open(unit(2),file="impG0_iw.ed")
    unit(3)=free_unit()
    open(unit(3),file="impSigma_iw.ed")
    unit(4)=free_unit()
    open(unit(4),file=trim(GRfile))
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
       Siw(ispin,:) = G0iw(ispin,:) - one/Giw(ispin,:)
       Swr(ispin,:) = G0wr(ispin,:) - one/Gwr(ispin,:)
       !Print GFs
    enddo
    do i=1,NL
       write(unit(1),"(F18.10,6(F18.10))")wm(i),(dimag(Giw(ispin,i)),dreal(Giw(ispin,i)),ispin=1,Nspin)
       write(unit(2),"(F18.10,6(F18.10))")wm(i),(dimag(one/G0iw(ispin,i)),dreal(one/G0iw(ispin,i)),ispin=1,Nspin)
       write(unit(3),"(F18.10,6(F18.10))")wm(i),(dimag(Siw(ispin,i)),dreal(Siw(ispin,i)),ispin=1,Nspin)
    enddo
    do i=1,Nw
       write(unit(4),"(F18.10,6(F18.10))")wr(i),(dimag(Gwr(ispin,i)),dreal(Gwr(ispin,i)),ispin=1,Nspin)
       write(unit(5),"(F18.10,6(F18.10))")wr(i),(dimag(one/G0wr(ispin,i)),dreal(one/G0wr(ispin,i)),ispin=1,Nspin)
       write(unit(6),"(F18.10,6(F18.10))")wr(i),(dimag(Swr(ispin,i)),dreal(Swr(ispin,i)),ispin=1,Nspin)
    enddo
    do i=1,6
       close(unit(i))
    enddo
  end subroutine print_imp_gf

end MODULE ED_GETGF
