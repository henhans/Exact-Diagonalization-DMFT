!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GETGF
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, only:c,cdg,bdecomp,delta_and
  implicit none
  private 
  public :: imp_getfunx,imp_getchi

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_getfunx()
    real(8)                  :: cdgmat(2),matcdg(2)
    integer,dimension(N)     :: ib(N)
    integer                  :: i,j,k,r,ll,m,in,is,ispin
    integer                  :: idg,jdg,isloop,jsloop,ia
    real(8)                  :: alfa=0.001d0,cc,spin1,spin2,spin12,peso1,peso2,peso12
    real(8)                  :: expterm,peso,de,w0,it,chij1,chij2,chij12
    complex(8)               :: iw
    complex(8),dimension(1:2,NL) :: G0iw
    complex(8),dimension(1:2,Nw) :: G0wr

    !----------------------------------------------
    !<i|C^+|j>=<in,is,idg|C^+|jn,js,jdg>=C^+_{ij} |
    !----------------------------------------------

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
       call eta(isloop,lastloop,file="Gimp_s1.eta")
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
                enddo
                peso=expterm/zeta_function
                de=espace(jsloop)%e(j)-espace(isloop)%e(i)
                matcdg=peso*cdgmat**2
                !
                do m=1,NL
                   iw=xi*wm(m)
                   Giw(ispin,m)=Giw(ispin,m)+matcdg(1)/(de+iw)
                enddo
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps+alfa*abs(w0))
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg(1)/(de+iw)
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
          call eta(isloop,lastloop,file="Gimp_s2.eta")
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
                enddo
                peso=expterm/zeta_function
                de=espace(jsloop)%e(j)-espace(isloop)%e(i)
                matcdg=peso*cdgmat**2
                !
                do m=1,NL
                   iw=xi*wm(m)
                   Giw(ispin,m)=Giw(ispin,m)+matcdg(1)/(de+iw)
                enddo
                do m=1,Nw 
                   w0=wr(m);iw=cmplx(w0,eps+alfa*abs(w0))
                   Gwr(ispin,m)=Gwr(ispin,m)+matcdg(1)*(one/(de+iw))
                enddo
                !
             enddo
          enddo
       enddo
       call stop_timer
    endif

    !###########################PRINTING######################################
    !Print convenience impurity functions:
    call system("rm -fv impG0_*.ed impSigma_*.ed "//trim(GMimp1file)//" "//trim(GRimp1file))
    do ispin=1,Nspin
       !Get Weiss Fields (from Bath):
       do i=1,NL
          iw=xi*wm(i)
          G0iw(ispin,i)= iw+xmu-ed0-delta_and(iw,ebath(ispin,:),vbath(ispin,:))
       enddo
       do i=1,Nw
          iw=cmplx(wr(i),eps)
          G0wr(ispin,i)= iw+xmu-ed0-delta_and(iw,ebath(ispin,:),vbath(ispin,:))
       enddo
       Siw(ispin,:) = G0iw(ispin,:) - one/Giw(ispin,:)
       Swr(ispin,:) = G0wr(ispin,:) - one/Gwr(ispin,:)
       !Print GFs
       call splot(trim(GMimp1file),wm,Giw(ispin,:),append=.true.)
       call splot(trim(GRimp1file),wr,Gwr(ispin,:),append=.true.)
       call splot("impG0_iw.ed",wm,one/G0iw(ispin,:),append=.true.)
       call splot("impG0_realw.ed",wr,one/G0wr(ispin,:),append=.true.)
       call splot("impSigma_iw.ed",wm,Siw(ispin,:),append=.true.)
       call splot("impSigma_realw.ed",wr,Swr(ispin,:),append=.true.)
    enddo

  end subroutine imp_getfunx





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
    integer                  :: idg,jdg,isloop,jsloop,ia
    real(8)                  :: cc,spin1,spin2,spin12,peso1,peso2,peso12
    real(8)                  :: expterm,peso,de,w0,it,chij1,chij2,chij12
    complex(8)               :: iw

    chitau=0.d0
    !Imaginary time susceptibility \X(tau). |<i|S_z|j>|^2
    call msg("Evaluating Chi(tau)")
    call start_timer
    do isloop=1,Nsect !loop over <i| total particle number
       call eta(isloop,lastloop,file="chi.eta")
       in=getin(isloop)
       is=getis(isloop)
       idg=deg(isloop)
       do i=1,idg 
          do j=1,idg
             chij1=0.d0
             chij2=0.d0
             chij12=0.d0
             do ll=1,idg 
                ia=nmap(isloop,ll)
                call bdecomp(ia,ib)
                spin1=dble(ib(1))-dble(ib(1+Ns))                
                chij1=chij1+espace(isloop)%M(ll,i)*spin1*espace(isloop)%M(ll,j)
             enddo
             peso1=chij1**2/zeta_function
             do m=0,Ltau 
                it=tau(m)
                chitau(m)=chitau(m)+&
                     exp(-it*espace(isloop)%e(i))*&
                     exp(-(beta-it)*espace(isloop)%e(j))*peso1/zeta_function
             enddo
             !
             ! do m=1,Nw 
             !    w0=wr(m);iw=cmplx(w0,eps+alfa*abs(w0))
             !    !Gwr(ispin,m)=Gwr(ispin,m)+matcdg(1)*(one/(de+iw))
             ! enddo
          enddo
       enddo
    enddo
    call stop_timer


    !###########################PRINTING######################################
    call splot(trim(CTimp1file),tau,chitau)
  end subroutine imp_getchi


  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_GETGF
