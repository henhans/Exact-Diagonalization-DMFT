!###################################################################
!PROGRAM  : FULLED
!TYPE     : main code
!PURPOSE  : Complete ED solution of DMFT equations (STD interface HM/PAM)
!AUTHORS  : A. Amaricci
!COMMENTS : USE THIS INTERFACE TO BUILD YOUR OWN CODE 
!###################################################################
program FULLED_AF_PAM_BETHE
  USE DMFT_FULLED
  USE FFTGF
  implicit none
  logical             :: converged
  !Human readable chemical potential (rescaled)
  real(8)             :: gzero,gmu,xmu0,gmu0,ntot
  real(8),allocatable :: wt(:),epsik(:),aa(:),bb(:)
  integer             :: Lk

  call read_input("inputED.in")

  !this shift contain |ep0-ed0| already so you don't need to add/remove it!
  gmu=xmu
  if(ed0 > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*tpd**2))
  if(ed0 < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*tpd**2))
  if(ed0 /=0.d0)xmu=gmu+gzero !true ED chemical potential
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  D=2.d0*ts; Lk=Nx**2 ; allocate(wt(Lk),epsik(Lk))
  call bethe_lattice(wt,epsik,Lk,D)

  !setup solver
  call ed_solver(status=-2)
  allocate(aa(size(epsiup)),bb(size(vup)))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver() 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe_af_pam

     !Fit the new bath, starting from the old bath + the supplied delta
     !aa=epsiup;bb=vup
     call chi2_fitgf(delta(1,:),epsiup,vup)
     !epsiup= weight*epsiup+(1.d0-weight)*aa
     !vup   = weight*vup   +(1.d0-weight)*bb
     !
     !aa=epsidw;bb=vdw
     call chi2_fitgf(delta(2,:),epsidw,vdw)
     !epsidw= weight*epsidw+(1.d0-weight)*aa
     !vdw   = weight*vdw   +(1.d0-weight)*bb
     !
     converged = check_convergence(delta(:,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo

  !finalize calculation
  call ed_solver(status=-1)


contains


  !+----------------------------------------+
  subroutine get_delta_bethe_af_pam
    integer                         :: i,j,ik
    real(8)                         :: w,nimp2UP,nimp2DW
    complex(8)                      :: iw,zita(Nspin),alpha(Nspin),gamma(Nspin)
    complex(8),dimension(Nspin,NL)  :: gd,gp,sp,sd,g0d
    complex(8),dimension(Nspin,Nw)  :: gdr,gpr,spr,sdr,g0dr,deltar
    real(8),dimension(Nspin,0:Ltau) :: gptau,gdtau
    real(8)                         :: mu(Nspin)
    delta=zero
    mu(1)=xmu-heff
    mu(2)=xmu+heff
    gp =zero ; gd=zero
    gpr=zero ; gdr=zero
    do i=1,NL
       iw=xi*wm(i)
       g0d(1,i) = iw + mu(1) -ed0 - delta_and(iw,epsiup,vup)
       g0d(2,i) = iw + mu(2) -ed0 - delta_and(iw,epsidw,vdw)
       sd(:,i)  = g0d(:,i) - one/Giw(:,i)
       alpha(:) = iw + mu(:) - ep0
       gamma(:) = iw + mu(:) - ed0 - sd(:,i)
       sp(:,i)  = tpd**2/gamma(:)
       zita     = alpha - sp(:,i)
       !
       do ik=1,Lk
          gp(:,i)=gp(:,i)+ wt(ik)/(zita(1)*zita(2)-epsik(ik)**2)
          gd(1,i)=gd(1,i)+ wt(ik)*(zita(2)-epsik(ik)**2/alpha(1))/&
               (zita(1)*zita(2)-epsik(ik)**2)
          gd(2,i)=gd(2,i)+ wt(ik)*(zita(1)-epsik(ik)**2/alpha(2))/&
               (zita(1)*zita(2)-epsik(ik)**2)
       enddo
       gp(1,i)=zita(2)*gp(1,i)
       gp(2,i)=zita(1)*gp(2,i)
       gd(:,i) = alpha/gamma*gd(:,i)
       !
       delta(1,i) =  tpd**2/(iw+mu(1)-ep0-gp(2,i)*D**2/4.d0)
       delta(2,i) =  tpd**2/(iw+mu(2)-ep0-gp(1,i)*D**2/4.d0)
       !
    enddo

    call fftgf_iw2tau(gp(1,:),gptau(1,0:),beta)
    call fftgf_iw2tau(gp(2,:),gptau(2,0:),beta)
    call fftgf_iw2tau(gd(1,:),gdtau(1,0:),beta)
    call fftgf_iw2tau(gd(2,:),gdtau(2,0:),beta)
    nimp2UP=-gptau(1,Ltau)
    nimp2DW=-gptau(2,Ltau)
    nimp2=nimp2UP+nimp2DW
    ntot=nimp1+nimp2

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0dr(1,i) = wr(i) + mu(1) - ed0 - delta_and(wr(i)+zero,epsiup,vup)
       g0dr(2,i) = wr(i) + mu(2) - ed0 - delta_and(wr(i)+zero,epsidw,vdw)
       sdr(:,i)  = g0dr(:,i) - one/Gwr(:,i)
       alpha(:)  = iw + mu(:) - ep0
       gamma(:)  = iw + mu(:) - ed0 - sdr(:,i)
       spr(:,i)  = tpd**2/gamma(:) !(iw + mu(:) - ed0 - sdr(:,i))
       zita      = alpha(:) - spr(:,i)
       !
       do ik=1,Lk
          gpr(:,i)=gpr(:,i)+ wt(ik)/(zita(1)*zita(2)-epsik(ik)**2)
          gdr(1,i)=gdr(1,i)+ wt(ik)*(zita(2)-epsik(ik)**2/alpha(1))/&
               (zita(1)*zita(2)-epsik(ik)**2)
          gdr(2,i)=gdr(2,i)+ wt(ik)*(zita(1)-epsik(ik)**2/alpha(2))/&
               (zita(1)*zita(2)-epsik(ik)**2)
       enddo
       gpr(1,i)=zita(2)*gpr(1,i)
       gpr(2,i)=zita(1)*gpr(2,i)
       gdr(:,i) = alpha/gamma*gdr(:,i)
       deltar(:,i) =  tpd**2/(iw+mu(:)-ep0-gpr(:,i)*D**2/4.d0)
    enddo
    call splot("Delta_iw.ed",wm,delta(1,:))
    call splot("Delta_iw.ed",wm,delta(2,:),append=.true.)
    call splot("Gpp_iw.ed",wm,gp(1,:))
    call splot("Gpp_iw.ed",wm,gp(2,:),append=.true.)
    call splot("Gdd_iw.ed",wm,gd(1,:))
    call splot("Gdd_iw.ed",wm,gd(2,:),append=.true.)
    call splot("Sigmapp_iw.ed",wm,sp(1,:))
    call splot("Sigmapp_iw.ed",wm,sp(2,:),append=.true.)
    call splot("Sigmadd_iw.ed",wm,sd(1,:))
    call splot("Sigmadd_iw.ed",wm,sd(2,:),append=.true.)
    !
    call splot("Gpp_tau.ed",tau,gptau(1,0:))
    call splot("Gpp_tau.ed",tau,gptau(2,0:),append=.true.)
    !
    call splot("Delta_realw.ed",wr,deltar(1,:))
    call splot("Delta_realw.ed",wr,deltar(2,:),append=.true.)
    call splot("Gpp_realw.ed",wr,gpr(1,:))
    call splot("Gpp_realw.ed",wr,gpr(2,:),append=.true.)
    call splot("Gdd_realw.ed",wr,gdr(1,:))
    call splot("Gdd_realw.ed",wr,gdr(2,:),append=.true.)
    call splot("Sigmapp_realw.ed",wr,spr(1,:))
    call splot("Sigmapp_realw.ed",wr,spr(2,:),append=.true.)
    call splot("Sigmadd_realw.ed",wr,sdr(1,:))
    call splot("Sigmadd_realw.ed",wr,sdr(2,:),append=.true.)
    !
    call splot("npup.npdw.np.ntot.ed",iloop,nimp2UP,nimp2DW,nimp2,ntot,append=TT)

  end subroutine get_delta_bethe_af_pam
  !+----------------------------------------+

end program FULLED_AF_PAM_BETHE



