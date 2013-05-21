!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GETGF
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, only:c,cdg,bdecomp,delta_and
  implicit none
  private 

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr

  public :: imp_getfunx,imp_getchi


contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_getfunx()
    real(8)                  :: cdgmat,matcdg
    integer,dimension(N)     :: ib(N)
    integer                  :: i,j,k,r,ll,m,in,is,ispin,unit(6)
    integer                  :: idg,jdg,isloop,jsloop,ia
    real(8)                  :: Ei,Ej
    real(8)                  :: cc,spin1,peso1
    real(8)                  :: expterm,peso,de,w0,it,chij1
    complex(8)               :: iw
    complex(8),dimension(1:2,NL) :: G0iw
    complex(8),dimension(1:2,Nw) :: G0wr

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

    !Paramagnetic .OR. spin=UP :
    !==========================================================================
    call msg("Evaluating G_imp_s1")
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
                      cdgmat=cdgmat+&
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
       call msg("Evaluating Gimp_s2")
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
                      cdgmat=cdgmat+&
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
          G0iw(ispin,i)= iw+xmu-delta_and(iw,ispin)
       enddo
       do i=1,Nw
          iw=cmplx(wr(i),eps)
          G0wr(ispin,i)= iw+xmu-delta_and(iw,ispin)
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
    deallocate(wm,tau,wr)
  end subroutine imp_getfunx





  !*********************************************************************
  !*********************************************************************
  !*********************************************************************





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_getchi()
    real(8)                  :: cdgmat,matcdg
    integer,dimension(N)     :: ib(N)
    integer                  :: i,j,k,r,ll,m,in,is,ispin
    integer                  :: idg,jdg,isloop,jsloop,ia,unit(2)
    real(8)                  :: cc,spin,peso,chij
    real(8)                  :: de,w0,it,Ei,Ej,expterm
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
    !Imaginary time susceptibility \X(tau). |<i|S_z|j>|^2
    call msg("Evaluating Chi(tau)")
    call start_timer
    do isloop=1,Nsect !loop over <i| total particle number
       call eta(isloop,lastloop,file="ETA_chi.ed")
       in=getin(isloop)
       is=getis(isloop)
       idg=deg(isloop)
       do i=1,idg 
          do j=1,idg
             chij=0.d0
             expterm=exp(-beta*espace(isloop)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idg 
                ia=nmap(isloop,ll)
                call bdecomp(ia,ib)
                spin=dble(ib(1))-dble(ib(1+Ns))                
                chij=chij+espace(isloop)%M(ll,i)*spin*espace(isloop)%M(ll,j)
             enddo
             Ei=espace(isloop)%e(i)
             Ej=espace(isloop)%e(j)
             de=Ej-Ei
             peso=chij**2/zeta_function
             do m=0,Ltau 
                it=tau(m)
                chitau(m)=chitau(m) + exp(-it*espace(isloop)%e(i))*&
                     exp(-(beta-it)*espace(isloop)%e(j))*peso
             enddo
             !
             do m=1,Nw
                w0=wr(m);iw=cmplx(w0,eps,8)
                chiw(m)=chiw(m)-exp(-beta*espace(isloop)%e(j))*&
                     (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso
             enddo
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
    deallocate(wm,tau,wr)
  end subroutine imp_getchi


  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_GETGF
