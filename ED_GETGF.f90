!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GETGF
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_GETH
  !USE ED_LANCZOS
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
    integer,dimension(N)         :: ib(N)
    integer                      :: i,j,k,r,ll,m,in,is,ispin
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

    !Paramagnetic .OR. spin=UP :
    !==========================================================================
    call msg("Evaluating G_imp_s1",unit=LOGfile)
    call start_timer
    ispin=1
    do isector=startloop,lastloop
       if(isector < getsector(2,0))cycle
       jsector=getCUPsector(isector);if(jsector==0)cycle
       call eta(isector,lastloop,file="ETA_Gimp_s1.ed")
       idim=getdim(isector)     !i-th sector dimension
       jdim=getdim(jsector)     !j-th sector dimension
       do i=1,idim          !loop over the states in the i-th sect.
          do j=1,jdim       !loop over the states in the j-th sect.
             cdgmat=0.d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm > cutoff)then
                do ll=1,jdim   !loop over the component of |j> (IN state!)
                   m=Hmap(jsector,ll) !map from IN state (j) 2 full Hilbert space
                   call bdecomp(m,ib)
                   if(ib(1) == 0)then
                      call cdg(1,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invHmap(isector,k) !map back to OUT sector (i)
                      cdgmat=cdgmat+&
                           espace(isector)%M(r,i)*cc*espace(jsector)%M(ll,j)
                   endif
                enddo
                Ei=espace(isector)%e(i)
                Ej=espace(jsector)%e(j)
                de=Ej-Ei
                peso=expterm/zeta_function
                matcdg=peso*cdgmat**2
                !
                do m=1,NL
                   iw=xi*wm(m)
                   Giw(ispin,m)=Giw(ispin,m)+matcdg/(iw+de)
                enddo
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps)
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg/(iw+de)
                enddo
                !
             endif
          enddo
       enddo
    enddo
    call stop_timer


    !spin=DW :
    !==========================================================================
    if(Nspin==2)then
       call msg("Evaluating Gimp_s2",unit=LOGfile)
       call start_timer
       ispin=2
       do isector=startloop,lastloop
          call eta(isector,lastloop,file="ETA_Gimp_s2.ed")
          if(isector < getsector(2,-2))cycle
          jsector=getCDWsector(isector);if(jsector==0)cycle
          idim=getdim(isector)     !i-th sector dimension
          jdim=getdim(jsector)     !j-th sector dimension
          do i=1,idim          !loop over the states in the i-th sect.
             do j=1,jdim       !loop over the states in the j-th sect.
                cdgmat=0.d0
                expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))     
                if(expterm < cutoff)cycle
                !
                do ll=1,jdim   !loop over the component of |j> (IN state!)
                   m=Hmap(jsector,ll) !map from IN state (j) 2 full Hilbert space
                   call bdecomp(m,ib)
                   if(ib(1+Ns) == 0)then
                      call cdg(1+Ns,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invHmap(isector,k) !map back to OUT sector (i)
                      cdgmat=cdgmat+&
                           espace(isector)%M(r,i)*cc*espace(jsector)%M(ll,j)
                   endif
                enddo
                Ei=espace(isector)%e(i)
                Ej=espace(jsector)%e(j)
                de=Ej-Ei
                peso=expterm/zeta_function
                matcdg=peso*cdgmat**2
                !
                do m=1,NL
                   iw=xi*wm(m)
                   Giw(ispin,m)=Giw(ispin,m)+matcdg/(iw+de)
                enddo
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps)
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg/(iw+de)
                enddo
                !
             enddo
          enddo
       enddo
       call stop_timer
    endif
    !###########################PRINTING######################################
    !Print convenience impurity functions:
    call print_imp_gf
    deallocate(wm,tau,wr)
  end subroutine full_ed_getgf




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine full_ed_getchi()
    USE STATISTICS
    real(8)                  :: cdgmat,matcdg
    integer,dimension(N)     :: ib(N),ibi(N),ibj(N)
    integer                  :: i,j,k,r,ll,m,in,is,ispin,kk
    integer                  :: idim,jdim,isector,jsector,ia,unit(2)
    real(8)                  :: cc,spin,peso,chij,weigth
    real(8)                  :: de,w0,it,Ei,Ej,expterm,Pchi(Nsect),totPchi
    complex(8)               :: iw
    type(histogram)      :: hist
    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(tau(0:Ltau))
    tau   = linspace(0.d0,beta,Ltau+1)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)
    chitau=0.d0
    chiw=zero

    hist = histogram_allocate(50)



    !Imaginary time susceptibility \X(tau). |<i|S_z|j>|^2
    call msg("Evaluating Chi_Sz",unit=LOGfile)
    call start_timer
    do isector=1,Nsect !loop over <i| total particle number
       call eta(isector,lastloop,file="ETA_chi.ed")
       in=getin(isector)
       is=getis(isector)
       idim=getdim(isector)
       Pchi(isector)=0.d0
       ! call histogram_reset(hist)
       ! call histogram_set_range_uniform(hist, 0.d0, 10.d0)
       do i=1,idim 
          do j=1,idim
             chij=0.d0
             expterm=exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim 
                ia=Hmap(isector,ll)
                call bdecomp(ia,ib)
                spin=real(ib(1),8)-real(ib(1+Ns),8) !nup - ndw
                chij=chij+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
             enddo



             Pchi(isector)=Pchi(isector)+chij**2
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ej-Ei
             if(chij**2 > 1.d-1)then
                ia=Hmap(isector,i)
                call bdecomp(ia,ibi)
                ia=Hmap(isector,j)
                call bdecomp(ia,ibj)
                write(300,"(4I5,2F15.9,2x,14I1,2x,14I1)")in,is,i,j,chij**2,de,&
                     (ibi(kk),kk=1,N),(ibj(kk),kk=1,N)

             endif
             peso=chij**2/zeta_function
             !
             ! call histogram_accumulate(hist,abs(de),chij**2)
             !
             do m=0,Ltau 
                it=tau(m)
                chitau(m)=chitau(m) + exp(-it*espace(isector)%e(i))*&
                     exp(-(beta-it)*espace(isector)%e(j))*peso
             enddo
             !
             do m=1,Nw
                w0=wr(m);iw=cmplx(w0,eps,8)
                chiw(m)=chiw(m)-exp(-beta*espace(isector)%e(j))*&
                     (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso
             enddo
          enddo
       enddo
       ! call histogram_print(hist,200)
       ! write(200,*)""
       ! write(201,*)isector,in,is
    enddo
    call stop_timer
    ! rewind(200)
    ! rewind(201)
    rewind(300)
    !###########################PRINTING######################################
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
       write(100,*)getin(i),getis(i),Pchi(i)
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
    integer                      :: in0,is0,idim0
    integer                      :: jn0,js0,jdim0
    real(8)                      :: norm0,sgn,gs,nup,ndw
    integer                      :: ib(N),k,r,Nlanc,Nitermax
    real(8)                      :: factor
    real(8),dimension(:),pointer :: vec
    real(8)                      :: egs
    real(8),allocatable          :: vvinit(:),alfa_(:),beta_(:),vout(:)

    call plain_lanczos_set_htimesv(HtimesV)
    !Initialize some functions
    Giw   =zero
    Gwr   =zero

    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)

    Nitermax=nGFitermax
    allocate(alfa_(Nitermax),beta_(Nitermax))

    factor=real(numzero,8)

    do izero=1,numzero   
       !GET THE GROUNDSTATE (make some checks)
       isect0=es_get_sector(groundstate,izero)
       isect0 = iszero(izero)
       in0    = getin(isect0)
       is0    = getis(isect0)
       idim0  = getdim(isect0)
       vec => es_get_vector(groundstate,izero)
       egs =  es_get_energy(groundstate,izero)
       norm0=sqrt(dot_product(vec,vec))
       if(norm0-1.d0>1.d-9)print*,"GS",izero,"is not normalized:",norm0

       !ADD ONE PARTICLE UP:
       !get cdg_up sector informations:
       jsect0 = getCDGUPsector(isect0);if(jsect0==0)cycle
       jdim0   = getdim(jsect0)
       jn0    = getin(jsect0)
       js0    = getis(jsect0)
       print*,'GetGF: sector C^+_up|gs>',jn0,js0,jdim0
       !allocate cdg_ip|gs> vector:
       allocate(vvinit(jdim0));vvinit=0.d0
       !build cdg_up|gs> vector:
       do m=1,idim0              !loop over |gs> components m
          i=Hmap(isect0,m)      !map m to full-Hilbert space state i
          call bdecomp(i,ib)    !decompose i into number representation ib=|1/0,1/0,1/0...>
          if(ib(1)==0)then      !if impurity is empty: proceed
             call cdg(1,i,r);sgn=dfloat(r)/dfloat(abs(r));r=abs(r) !apply cdg_up (1), bring from i to r
             j=invHmap(jsect0,r)                                   !map r back to cdg_up sector jsect0
             vvinit(j) = sgn*vec(m)                                !build the cdg_up|gs> state
          endif
       enddo
       !normalize the cdg_up|gs> state
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       !get cdg_up-sector Hamiltonian
       allocate(H0(jdim0,jdim0))
       call full_ed_geth(jsect0,H0)
       !Tri-diagonalize w/ Lanczos the resolvant:
       allocate(vout(jdim0))
       vout= 0.d0 ; alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
       call plain_lanczos_tridiag(vvinit,alfa_,beta_,nitermax)
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,1,1)
       deallocate(H0,vout,vvinit)


       !REMOVE ONE PARTICLE UP:
       jsect0 = getCUPsector(isect0);if(jsect0==0)cycle
       jdim0   = getdim(jsect0)
       jn0    = getin(jsect0)
       js0    = getis(jsect0)
       print*,'GetGF: sector C_up|gs>',jn0,js0,jdim0
       allocate(vvinit(jdim0)) ; vvinit=0.d0
       do m=1,idim0              !loop over |gs> components m
          i=Hmap(isect0,m)      !map m to full-Hilbert space state i
          call bdecomp(i,ib)    !decompose i into number representation ib=|1/0,1/0,1/0...>
          if(ib(1)==1)then      !if impurity is empty: proceed
             call c(1,i,r);sgn=dfloat(r)/dfloat(abs(r));r=abs(r) !apply c_up (1), bring from i to r
             j=invHmap(jsect0,r)                                 !map r back to c_up sector jsect0
             vvinit(j) = sgn*vec(m)                              !build the c_up|gs> state
          endif
       enddo
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       allocate(H0(jdim0,jdim0))
       call full_ed_geth(jsect0,H0)
       allocate(vout(jdim0))
       vout= 0.d0 ; alfa_=0.d0 ; beta_=0.d0
       call plain_lanczos_tridiag(vvinit,alfa_,beta_,nitermax)
       call add_to_lanczos_gf(norm0,egs,nitermax,alfa_,beta_,-1,1)
       deallocate(H0,vout,vvinit)


    enddo

    Giw=Giw/factor
    Gwr=Gwr/factor

    !###########################PRINTING######################################
    !Print convenience impurity functions:
    call print_imp_gf
    deallocate(wm,wr)

  end subroutine lanc_ed_getgf

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_gf(vnorm,emin,nlanc,alanc,blanc,isign,ispin)
    real(8),dimension(nlanc)                         :: alanc,blanc 
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
          Giw(ispin,i)=Giw(ispin,i) + vnorm**2*Z(1,j)**2/(xi*wm(i) - isign*(diag(j)-emin))
       enddo
    enddo
    do i=1,Nw
       do j=1,nlanc
          Gwr(ispin,i)=Gwr(ispin,i) + vnorm**2*Z(1,j)**2/(dcmplx(wr(i),eps)-isign*(diag(j)-emin))
       enddo
    enddo
  end subroutine add_to_lanczos_gf





  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
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
