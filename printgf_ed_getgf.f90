  
  !Print convenience impurity functions:
  select case(Nspin)
  case default
     !Print GFs
     call splot(trim(GMimp1file),wm,Giw(1,:))
     call splot(trim(GRimp1file),wr,Gwr(1,:))
     if(Nimp==2)then
        call splot(trim(GMimp2file),wm,G2iw(1,:))
        call splot(trim(GRimp2file),wr,G2wr(1,:))
     endif


     !Get Weiss Fields (from Bath):
     do i=1,NL
        iw=xi*wm(i)
        G0iw(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
     enddo
     do i=1,Nw
        iw=cmplx(wr(i),eps)
        G0wr(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
     enddo
     call splot("impG0_iw.ed",wm(1:Nfit),one/G0iw(1,1:Nfit))
     call splot("impG0_realw.ed",wr,one/G0wr(1,:))

     !If 1imp print Sigma (PAM is special and this Sigma is not correct)
     if(Nimp==1)then
        Siw(1,:) = G0iw(1,:) - one/Giw(1,:)
        Swr(1,:) = G0wr(1,:) - one/Gwr(1,:)
        call splot("impSigma_iw.ed",wm(1:Nfit),Siw(1,1:Nfit))
        call splot("impSigma_realw.ed",wr,Swr(1,:))
     endif

  case (2)
     !Print GFs
     call splot(trim(GMimp1file),wm,Giw(1,:))
     call splot(trim(GMimp1file),wm,Giw(2,:),append=.true.)
     call splot(trim(GRimp1file),wr,Gwr(1,:))
     call splot(trim(GRimp1file),wr,Gwr(2,:),append=.true.)
     if(Nimp==2)then
        call splot(trim(GMimp2file),wm,G2iw(1,:))
        call splot(trim(GMimp2file),wm,G2iw(2,:),append=.true.)
        call splot(trim(GRimp2file),wr,G2wr(1,:))
        call splot(trim(GRimp2file),wr,G2wr(2,:),append=.true.)
     endif

     !Get Weiss Fields (from Bath):
     do i=1,NL
        iw=xi*wm(i)
        G0iw(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
        G0iw(2,i)= iw + xmu - ed0 - delta_and(iw,epsidw,vdw)
     enddo
     do i=1,Nw
        iw=cmplx(wr(i),eps)
        G0wr(1,i)= iw + xmu - ed0 - delta_and(iw,epsiup,vup)
        G0wr(2,i)= iw + xmu - ed0 - delta_and(iw,epsidw,vdw)
     enddo
     call splot("impG0_iw.ed",wm(1:Nfit),one/G0iw(1,1:Nfit))
     call splot("impG0_iw.ed",wm(1:Nfit),one/G0iw(2,1:Nfit),append=.true.)
     call splot("impG0_realw.ed",wr,one/G0wr(1,:))
     call splot("impG0_realw.ed",wr,one/G0wr(2,:),append=.true.)
     if(Nimp==1)then
        Siw = G0iw - one/Giw
        Swr = G0wr - one/Gwr
        call splot("impSigma_iw.ed",wm(1:Nfit),Siw(1,1:Nfit))
        call splot("impSigma_iw.ed",wm(1:Nfit),Siw(2,1:Nfit),append=.true.)
        call splot("impSigma_realw.ed",wr,Swr(1,:))
        call splot("impSigma_realw.ed",wr,Swr(2,:),append=.true.)
     endif
  end select
  !

  !
  if(chiflag)then
     call splot(trim(CTimp1file),tau,chitau)
     if(Nimp==2)then
        call splot(trim(CTimp2file),tau,chi2tau)
        call splot(trim(CTimpAfile),tau,chi12tau)
     endif
  endif
  !
