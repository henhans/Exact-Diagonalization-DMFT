!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GETGF
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_GETH
  USE ED_LANCZOS
  !
  implicit none
  private 

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr

  public :: imp_getfunx
  public :: lanc_getgf
  public :: imp_getchi

contains


  subroutine lanc_getgf()
    integer                      :: i,izero,isect0,jsect0,m,j
    integer                      :: in0,is0,idg0
    integer                      :: jn0,js0,jdg0
    real(8)                      :: norm0,sgn,gs,nup,ndw
    integer                      :: ib(N),k,r,Nlanc,Nitermax
    real(8)                      :: factor
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
       norm0=sqrt(dot_product(groundstate(izero)%vec,groundstate(izero)%vec))
       if(norm0-1.d0>1.d-9)print*,"GS",izero,"is not normalized:",norm0
       isect0 = iszero(izero)
       in0    = getin(isect0)
       is0    = getis(isect0)
       idg0   = deg(isect0)
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
       call plain_lanczos_tridiag(vvinit,alfa_,beta_,nitermax)
       call add_to_lanczos_gf(norm0,groundstate(izero)%egs,nitermax,alfa_,beta_,1,1)
       deallocate(H0,vout,vvinit)


       !REMOVE ONE PARTICLE UP:
       jsect0 = getCUPloop(isect0);if(jsect0==0)cycle
       jdg0   = deg(jsect0)
       jn0    = getin(jsect0)
       js0    = getis(jsect0)
       print*,'sector C_up|gs>',jn0,js0,jdg0
       allocate(vvinit(jdg0));vvinit=0.d0
       do m=1,idg0              !loop over |gs> components m
          i=nmap(isect0,m)      !map m to full-Hilbert space state i
          call bdecomp(i,ib)    !decompose i into number representation ib=|1/0,1/0,1/0...>
          if(ib(1)==1)then      !if impurity is empty: proceed
             call c(1,i,r);sgn=dfloat(r)/dfloat(abs(r));r=abs(r) !apply c_up (1), bring from i to r
             j=invnmap(jsect0,r)                                   !map r back to c_up sector jsect0
             vvinit(j) = sgn*groundstate(izero)%vec(m)             !build the c_up|gs> state
          endif
       enddo
       norm0=sqrt(dot_product(vvinit,vvinit))
       vvinit=vvinit/norm0
       allocate(H0(jdg0,jdg0))
       call imp_geth(jsect0,H0)
       allocate(vout(jdg0))
       vout= 0.d0 ; alfa_=0.d0 ; beta_=0.d0
       call plain_lanczos_tridiag(vvinit,alfa_,beta_,nitermax)
       call add_to_lanczos_gf(norm0,groundstate(izero)%egs,nitermax,alfa_,beta_,-1,1)
       deallocate(H0,vout,vvinit)
    enddo

    Giw=Giw/factor
    Gwr=Gwr/factor

    !###########################PRINTING######################################
    !Print convenience impurity functions:
    call print_imp_gf
    deallocate(wm,wr)

  end subroutine lanc_getgf






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_getfunx()
    real(8)                      :: cdgmat(2),matcdg(2)
    integer,dimension(N)         :: ib(N)
    integer                      :: i,j,k,r,ll,m,in,is,ispin
    integer                      :: idg,jdg,isloop,jsloop,ia
    real(8)                      :: Ei,Ej,cc,spin1,spin2,spin12,peso1,peso2,peso12
    real(8)                      :: expterm,peso,de,w0,it,chij1,chij2,chij12
    complex(8)                   :: iw
    !----------------------------------------------
    !<i|C^+|j>=<in,is,idg|C^+|jn,js,jdg>=C^+_{ij} |
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
    if(Nimp==2)then
       G2iw = zero  
       G2wr=zero
    endif

    !Paramagnetic .OR. spin=UP :
    !==========================================================================
    call msg("Evaluating G_imp_s1",unit=LOGfile)
    call start_timer
    ispin=1
    do isloop=startloop,lastloop
       if(isloop < getloop(2,0))cycle
       jsloop=getCUPloop(isloop);if(jsloop==0)cycle
       call eta(isloop,lastloop,file="ETA_Gimp_s1.ed")
       idg=deg(isloop)     !i-th sector dimension
       jdg=deg(jsloop)     !j-th sector dimension
       do i=1,idg          !loop over the states in the i-th sect.
          do j=1,jdg       !loop over the states in the j-th sect.
             cdgmat=0.d0
             expterm=exp(-beta*espace(isloop)%e(i))+exp(-beta*espace(jsloop)%e(j))
             if(expterm > cutoff)then
                do ll=1,jdg   !loop over the component of |j> (IN state!)
                   m=nmap(jsloop,ll) !map from IN state (j) 2 full Hilbert space
                   call bdecomp(m,ib)
                   if(ib(1) == 0)then
                      call cdg(1,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invnmap(isloop,k) !map back to OUT sector (i)
                      cdgmat(1)=cdgmat(1)+&
                           espace(isloop)%M(r,i)*cc*espace(jsloop)%M(ll,j)
                   endif
                   if(Nimp==2)then
                      if(ib(2) == 0)then
                         call cdg(2,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                         r=invnmap(isloop,k) !map back to OUT sector (i)
                         cdgmat(2)=cdgmat(2)+&
                              espace(isloop)%M(r,i)*cc*espace(jsloop)%M(ll,j)
                      endif
                   endif
                enddo
                Ei=espace(isloop)%e(i)
                Ej=espace(jsloop)%e(j)
                de=Ej-Ei
                peso=expterm/zeta_function
                matcdg=peso*cdgmat**2
                !
                do m=1,NL
                   iw=xi*wm(m)
                   Giw(ispin,m)=Giw(ispin,m)+matcdg(1)/(iw+de)
                enddo
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps)
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg(1)/(iw+de)
                enddo
                !
                if(Nimp==2)then
                   do m=1,NL
                      iw=xi*wm(m)
                      G2iw(ispin,m)=G2iw(ispin,m)+matcdg(2)/(de+iw)
                   enddo
                   do m=1,Nw 
                      w0=wr(m);iw=cmplx(w0,eps)
                      G2wr(ispin,m)=G2wr(ispin,m)+matcdg(2)/(de+iw)
                   enddo
                endif
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
       do isloop=startloop,lastloop
          call eta(isloop,lastloop,file="ETA_Gimp_s2.ed")
          if(isloop < getloop(2,-2))cycle
          jsloop=getCDWloop(isloop);if(jsloop==0)cycle
          idg=deg(isloop)     !i-th sector dimension
          jdg=deg(jsloop)     !j-th sector dimension
          do i=1,idg          !loop over the states in the i-th sect.
             do j=1,jdg       !loop over the states in the j-th sect.
                cdgmat=0.d0
                expterm=exp(-beta*espace(isloop)%e(i))+exp(-beta*espace(jsloop)%e(j))     
                if(expterm < cutoff)cycle
                !
                do ll=1,jdg   !loop over the component of |j> (IN state!)
                   m=nmap(jsloop,ll) !map from IN state (j) 2 full Hilbert space
                   call bdecomp(m,ib)
                   if(ib(1+Ns) == 0)then
                      call cdg(1+Ns,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invnmap(isloop,k) !map back to OUT sector (i)
                      cdgmat(1)=cdgmat(1)+&
                           espace(isloop)%M(r,i)*cc*espace(jsloop)%M(ll,j)
                   endif
                   if(Nimp==2 .AND. ib(2+Ns) == 0)then
                      call cdg(2+Ns,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invnmap(isloop,k) !map back to OUT sector (i)
                      cdgmat(2)=cdgmat(2)+&
                           espace(isloop)%M(r,i)*cc*espace(jsloop)%M(ll,j)
                   endif
                enddo
                Ei=espace(isloop)%e(i)
                Ej=espace(jsloop)%e(j)
                de=Ej-Ei
                peso=expterm/zeta_function
                matcdg=peso*cdgmat**2
                !
                do m=1,NL
                   iw=xi*wm(m)
                   Giw(ispin,m)=Giw(ispin,m)+matcdg(1)/(iw+de)
                enddo
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps)
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg(1)*(one/(iw+de))
                enddo
                !
                if(Nimp==2)then
                   do m=1,NL
                      iw=xi*wm(m)
                      G2iw(ispin,m)=G2iw(ispin,m)+matcdg(2)/(de+iw)
                   enddo
                   do m=1,Nw 
                      w0=wr(m);iw=cmplx(w0,eps)
                      G2wr(ispin,m)=G2wr(ispin,m)+matcdg(2)*(one/(de+iw))
                   enddo
                endif
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

  end subroutine imp_getfunx



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


    if(Nimp==2)then
       unit(1)=free_unit()
       open(unit(1),file=trim(GM2file))
       unit(2)=free_unit()
       open(unit(2),file=trim(GR2file))
       do i=1,NL
          write(unit(1),"(F18.10,6(F18.10))")wm(i),(dimag(G2iw(ispin,i)),dreal(G2iw(ispin,i)),ispin=1,Nspin)
       enddo
       do i=1,Nw
          write(unit(2),"(F18.10,6(F18.10))")wr(i),(dimag(G2wr(ispin,i)),dreal(G2wr(ispin,i)),ispin=1,Nspin)
       enddo
       close(unit(1))
       close(unit(2))
    end if



  end subroutine print_imp_gf

  !*********************************************************************
  !*********************************************************************
  !*********************************************************************





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_getchi()
    real(8)                  :: cdgmat(2),matcdg(2)
    integer,dimension(N)     :: ib(N)
    integer                  :: i,j,k,r,ll,m,in,is,ispin
    integer                  :: idg,jdg,isloop,jsloop,ia,unit(6)
    real(8)                  :: Ei,Ej,cc,spin1,spin2,spin12,peso1,peso2,peso12
    real(8)                  :: expterm,peso,de,w0,it,chij1,chij2,chij12
    complex(8)               :: iw

    !Freq. arrays
    allocate(wm(NL))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(tau(0:Ltau))
    tau   = linspace(0.d0,beta,Ltau+1)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw)

    Chitau=0.d0
    Chiw=zero
    if(Nimp==2)then
       Chi2tau=0.d0
       Chi2w=zero
       Chi12tau=0.d0
       Chi12w=zero
    endif
    !Imaginary time susceptibility \X(tau). |<i|S_z|j>|^2
    call msg("Evaluating Chi(tau)",unit=LOGfile)
    call start_timer
    do isloop=1,Nsect !loop over <i| total particle number
       call eta(isloop,lastloop,file="ETA_chi.ed")
       in=getin(isloop)
       is=getis(isloop)
       idg=deg(isloop)
       do i=1,idg 
          do j=1,idg
             chij1=0.d0
             chij2=0.d0
             chij12=0.d0
             expterm=exp(-beta*espace(isloop)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idg 
                ia=nmap(isloop,ll)
                call bdecomp(ia,ib)
                spin1=dble(ib(1))-dble(ib(1+Ns))                
                chij1=chij1+espace(isloop)%M(ll,i)*spin1*espace(isloop)%M(ll,j)
                if(Nimp==2)then
                   spin2=dble(ib(2))-dble(ib(2+Ns))
                   spin12=spin1+spin2
                   chij2=chij2+espace(isloop)%M(ll,i)*spin2*espace(isloop)%M(ll,j)
                   chij12=chij12+espace(isloop)%M(ll,i)*spin12*espace(isloop)%M(ll,j)
                endif
             enddo
             Ei=espace(isloop)%e(i)
             Ej=espace(isloop)%e(j)
             de=Ej-Ei
             peso1=chij1**2/zeta_function
             do m=0,Ltau 
                it=tau(m)
                chitau(m)=chitau(m) + exp(-it*espace(isloop)%e(i))*&
                     exp(-(beta-it)*espace(isloop)%e(j))*peso1
             enddo
             !
             do m=1,Nw
                w0=wr(m);iw=cmplx(w0,eps,8)
                chiw(m)=chiw(m)-exp(-beta*espace(isloop)%e(j))*&
                     (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso1
             enddo
             if(Nimp==2)then
                peso2=chij2**2/zeta_function
                peso12=chij12**2/zeta_function
                do m=0,Ltau 
                   it=tau(m)
                   chi2tau(m)=chi2tau(m)+exp(-it*espace(isloop)%e(i))*&
                        exp(-(beta-it)*espace(isloop)%e(j))*peso2
                   chi12tau(m)=chi12tau(m)+exp(-it*espace(isloop)%e(i))*&
                        exp(-(beta-it)*espace(isloop)%e(j))*peso12
                enddo
                do m=1,Nw
                   w0=wr(m);iw=cmplx(w0,eps,8)
                   chi2w(m)=chi2w(m)-exp(-beta*espace(isloop)%e(j))*&
                        (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso2
                   chi12w(m)=chi12w(m)-exp(-beta*espace(isloop)%e(j))*&
                        (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso12
                enddo
             endif
          enddo
       enddo
    enddo
    call stop_timer

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

    if(Nimp==2)then
       unit(3)=free_unit()
       open(unit(3),file=trim(CT2file))
       unit(4)=free_unit()
       open(unit(4),file=trim(CW2file))
       unit(5)=free_unit()
       open(unit(5),file=trim(CT12file))
       unit(6)=free_unit()
       open(unit(6),file=trim(CW12file))
       do i=0,Ltau
          write(unit(3),*)tau(i),chi2tau(i)
          write(unit(5),*)tau(i),chi12tau(i)
       enddo
       do i=1,Nw
          if(wr(i)>=0.d0)write(unit(4),*)wr(i),dimag(chi2w(i)),dreal(chi2w(i))
          if(wr(i)>=0.d0)write(unit(6),*)wr(i),dimag(chi12w(i)),dreal(chi12w(i))
       enddo
       close(unit(3))
       close(unit(4))
       close(unit(5))
       close(unit(6))
    endif

    deallocate(wm,tau,wr)
  end subroutine imp_getchi




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




end MODULE ED_GETGF
