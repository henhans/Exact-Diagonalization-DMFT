!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GETGF
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, only:c,cdg,bdecomp,delta_and
  implicit none
  private 
  public :: imp_getfunx
contains

  subroutine imp_getfunx()
    real(8)                  :: cdgmat(2),matcdg(2)
    integer,dimension(N)     :: ib(N)
    integer                  :: i,j,k,r,ll,m,in,is,ispin
    integer                  :: idg,jdg,isloop,jsloop,ia
    real(8)                  :: alfa,cc,spin1,spin2,spin12,peso1,peso2,peso12
    real(8)                  :: expterm,peso,de,w0,it,chij1,chij2,chij12
    complex(8)               :: iw
    complex(8),dimension(1:2,NL) :: G0iw,Siw
    complex(8),dimension(1:2,Nw) :: G0wr,Swr

    !----------------------------------------------
    !<i|C^+|j>=<in,is,idg|C^+|jn,js,jdg>=C^+_{ij} |
    !----------------------------------------------

    alfa=0.001d0

    !Initialize some functions
    Giw   =zero
    Gwr   =zero
    chitau=0.d0

    if(Nimp==2)then
       G2iw = zero  
       G2wr=zero
    endif


    !Paramagnetic .OR. spin=UP :
    !==========================================================================
    call msg("Evaluating G_imp_1")
    call start_timer
    ispin=1
    do isloop=startloop,lastloop
       if(isloop < getloop(2,0))cycle
       jsloop=getCUPloop(isloop);if(jsloop==0)cycle
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
                   if(Nimp==2 .AND. ib(2) == 0)then
                      call cdg(2,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                      r=invnmap(isloop,k) !map back to OUT sector (i)
                      cdgmat(2)=cdgmat(2)+&
                           espace(isloop)%M(r,i)*cc*espace(jsloop)%M(ll,j)
                   endif
                enddo
                peso=expterm/zeta
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
                if(Nimp==2)then
                   do m=1,NL
                      iw=xi*wm(m)
                      G2iw(ispin,m)=G2iw(ispin,m)+matcdg(2)/(de+iw)
                   enddo
                   do m=1,Nw 
                      w0=wr(m);iw=cmplx(w0,eps+alfa*abs(w0))
                      G2wr(ispin,m)=G2wr(ispin,m)+matcdg(2)*(one/(de+iw))
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
       call msg("Evaluating Gimp_2")
       call start_timer
       ispin=2
       do isloop=startloop,lastloop
          call eta(isloop,lastloop,file="ETA.Gimp_2")
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
                peso=expterm/zeta
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
                if(Nimp==2)then
                   do m=1,NL
                      iw=xi*wm(m)
                      G2iw(ispin,m)=G2iw(ispin,m)+matcdg(2)/(de+iw)
                   enddo
                   do m=1,Nw 
                      w0=wr(m);iw=cmplx(w0,eps+alfa*abs(w0))
                      G2wr(ispin,m)=G2wr(ispin,m)+matcdg(2)*(one/(de+iw))
                   enddo
                endif
                !
             enddo
          enddo
       enddo
       call stop_timer
    endif


    !Imaginary time susceptibility \X(tau). |<i|S_z|j>|^2
    if(chiflag)then
       call msg("Evaluating \Chi(tau)")
       call start_timer
       do isloop=1,Nsect !loop over <i| total particle number
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
                   if(Nimp==2)then
                      spin2=dble(ib(2))-dble(ib(2+Ns))
                      spin12=spin1+spin2
                      chij2=chij2+espace(isloop)%M(ll,i)*spin2*espace(isloop)%M(ll,j)
                      chij12=chij12+espace(isloop)%M(ll,i)*spin12*espace(isloop)%M(ll,j)
                   endif
                enddo
                peso1=chij1**2/zeta
                do m=0,Ltau 
                   it=tau(m)
                   chitau(m)=chitau(m)+&
                        exp(-it*espace(isloop)%e(i))*&
                        exp(-(beta-it)*espace(isloop)%e(j))*peso/zeta
                enddo
                !
                if(Nimp==2)then
                   peso2=chij2**2/zeta
                   peso12=chij12**2/zeta
                   do m=1,Ltau 
                      it=tau(m)
                      chi2tau(m)=chi2tau(m)+&
                           exp(-it*espace(isloop)%e(i))*&
                           exp(-(beta-it)*espace(isloop)%e(j))*peso2/zeta
                      chi12tau(m)=chi12tau(m)+&
                           exp(-it*espace(isloop)%e(i))*&
                           exp(-(beta-it)*espace(isloop)%e(j))*peso12/zeta
                   enddo
                endif
             enddo
          enddo
       enddo
       call stop_timer
    endif


    !Print convenience impurity functions:
    select case(Nspin)
    case default
       call splot(trim(GMimp1file),wm,Giw(1,:))
       call splot(trim(GRimp1file),wr,Gwr(1,:))
       if(Nimp==2)then
          call splot(trim(GMimp2file),wm,G2iw(1,:))
          call splot(trim(GRimp2file),wr,G2wr(1,:))
       endif

       do i=1,NL          
          iw=xi*wm(i)
          G0iw(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
          Siw(1,i) = G0iw(1,i) - one/Giw(1,i)
       enddo
       do i=1,Nw
          iw=cmplx(wr(i),eps)
          G0wr(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
          Swr(1,i) = G0wr(1,i) - one/Gwr(1,i)
       enddo
       call splot("impG0_iw.ed",wm(1:Nfit),one/G0iw(1,1:Nfit))
       call splot("impSigma_iw.ed",wm(1:Nfit),Siw(1,1:Nfit))
       call splot("impG0_realw.ed",wr,one/G0wr(1,:))
       call splot("impSigma_realw.ed",wr,Swr(1,:))


    case (2)
       call splot(trim(GMimp1file),wm,Giw(1,:))
       call splot(trim(GMimp1file),wm,Giw(2,:),append=TT)
       call splot(trim(GRimp1file),wr,Gwr(1,:))
       call splot(trim(GRimp1file),wr,Gwr(2,:),append=TT)
       if(Nimp==2)then
          call splot(trim(GMimp2file),wm,G2iw(1,:))
          call splot(trim(GMimp2file),wm,G2iw(2,:),append=TT)
          call splot(trim(GRimp2file),wr,G2wr(1,:))
          call splot(trim(GRimp2file),wr,G2wr(2,:),append=TT)
       endif

       do i=1,NL
          iw=xi*wm(i)
          G0iw(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
          G0iw(2,i)= iw + xmu - ed0 - delta_and(iw,epsidw,vdw)
       enddo
       Siw = G0iw - one/Giw

       do i=1,Nw
          iw=cmplx(wr(i),eps)
          G0wr(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
          G0wr(2,i)= iw + xmu - ed0 - delta_and(iw,epsidw,vdw)
       enddo
       Swr = G0wr - one/Gwr
       call splot("impG0_iw.ed",wm(1:Nfit),one/G0iw(1,1:Nfit))
       call splot("impG0_iw.ed",wm(1:Nfit),one/G0iw(2,1:Nfit),append=TT)
       call splot("impSigma_iw.ed",wm(1:Nfit),Siw(1,1:Nfit))
       call splot("impSigma_iw.ed",wm(1:Nfit),Siw(2,1:Nfit),append=TT)
       call splot("impG0_realw.ed",wr,one/G0wr(1,:))
       call splot("impG0_realw.ed",wr,one/G0wr(2,:),append=TT)
       call splot("impSigma_realw.ed",wr,Swr(1,:))
       call splot("impSigma_realw.ed",wr,Swr(2,:),append=TT)
    end select
    !
    if(chiflag)then
       call splot(trim(CTimp1file),tau,chitau)
       if(Nimp==2)then
          call splot(trim(CTimp2file),tau,chi2tau)
          call splot(trim(CTimpAfile),tau,chi12tau)
       endif
    endif
  end subroutine imp_getfunx



end MODULE ED_GETGF
