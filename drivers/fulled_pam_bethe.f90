!###################################################################
!PROGRAM  : FULLED
!TYPE     : main code
!PURPOSE  : Complete ED solution of DMFT equations (STD interface HM/PAM)
!AUTHORS  : A. Amaricci
!COMMENTS : USE THIS INTERFACE TO BUILD YOUR OWN CODE 
!###################################################################
program fullED
  USE DMFT_FULLED
  USE FFTGF
  implicit none
  logical :: converged

  !Human readable chemical potential (rescaled)
  real(8) :: gzero,gmu,xmu0,gmu0,ntot

  call read_input("inputED.in")

  !this shift contain |ep0-ed0| already so you don't need to add/remove it!
  gmu=xmu
  if(ed0 > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*tpd**2))
  if(ed0 < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*tpd**2))
  if(ed0 /=0.d0)xmu=gmu+gzero !true ED chemical potential
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero


  !setup solver
  call ed_solver(status=-2)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver() 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe_pam

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(1,:),epsiup,vup)
     epsidw=epsiup;vdw=vup
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo

  !finalize calculation
  call ed_solver(status=-1)


contains


  !+----------------------------------------+
  subroutine get_delta_bethe_pam
    integer    :: i,j
    real(8)    :: w
    complex(8) :: iw,zita
    complex(8) :: gp(NL),gd(NL),sp(NL),sd(NL),g0p(NL),g0d(NL)
    complex(8) :: gpr(Nw),gdr(Nw),spr(Nw),sdr(Nw),g0pw(Nw),g0dr(Nw)
    real(8)    :: gptau(0:Ltau),gdtau(0:Ltau)

    delta=zero
    do i=1,NL
       iw=xi*wm(i)
       g0d(i) = iw + xmu -ed0 - delta_and(iw,epsiup,vup)
       sd(i)  = g0d(i) - one/Giw(1,i)
       sp(i)  = tpd**2/(iw + xmu - ed0 - sd(i))
       zita   = iw + xmu -ep0 - sp(i)
       gp(i)  = gfbethe(wm(i),zita,d)
       gd(i)  = one/(iw+xmu-ed0-sd(i)) + tpd**2/(iw+xmu-ed0-sd(i))**2*gp(i)
       !
       delta(1,i) =  tpd**2/(iw+xmu-d**2*gp(i)/4.d0)
       !
    enddo

    call fftgf_iw2tau(gp,gptau,beta)
    call fftgf_iw2tau(gd,gdtau,beta)
    nimp2=-2.d0*real(gptau(Ltau),8)
    ntot=nimp1+nimp2

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0dr(i) = wr(i) + xmu - ed0 - delta_and(wr(i)+zero,epsiup,vup)
       sdr(i)  = g0dr(i) - one/Gwr(1,i)
       spr(i)  = tpd**2/(iw + xmu - ed0 - sdr(i))
       zita    = iw + xmu -ep0 - spr(i)
       gpr(i)  = gfbether(wr(i),zita,D)
       gdr(i)  = one/(iw+xmu-ed0-sdr(i)) + tpd**2/(iw+xmu-ed0-sdr(i))**2*gpr(i)
    enddo
    call splot("Delta_iw.ed",wm,delta(1,:))
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("G_tau_ddpp.ed",tau,gdtau,gptau)
    call splot("G0dd_realw.ed",wr,one/g0dr)
    call splot("Gdd_realw.ed",wr,gdr)
    call splot("Gpp_realw.ed",wr,gpr)
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)
    call splot("Sigmapp_realw.ed",wr,spr)
    call splot("Sigmadd_realw.ed",wr,sdr)
    call splot("np.ntot.ed",nimp2,ntot,append=TT)
  end subroutine get_delta_bethe_pam
  !+----------------------------------------+

end program fullED



