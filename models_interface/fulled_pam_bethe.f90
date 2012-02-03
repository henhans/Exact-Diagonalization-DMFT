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


  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     call ed_solver() !Solve the EFFECTIVE IMPURITY PROBLEM
     call get_delta_bethe_pam
     call chi2_fitgf(delta(1,:),epsiup,vup)
     epsidw=epsiup;vdw=vup
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     call end_loop
  enddo
  call ed_solver()


contains


  !+----------------------------------------+
  subroutine get_delta_bethe_pam
    integer    :: i,j
    real(8)    :: w
    complex(8) :: iw,g0d,sd,sp,zita,gp
    complex(8) :: gploc(NL),gdloc(NL)
    complex(8) :: gprloc(Nw),gdrloc(Nw)
    real(8)    :: gptau(0:Ltau)
    delta=zero
    do i=1,NL
       iw=xi*wm(i)
       g0d= iw + xmu -ed0 -delta_and(iw,epsiup,vup)
       sd  = g0d - one/Giw(1,i)
       sp  = tpd**2/(iw + xmu - ed0 - sd)
       zita= iw + xmu -ep0 - sp
       gploc(i)  = gfbethe(wm(i),zita,d)
       gdloc(i)  = one/(iw+xmu-ed0-sd) + tpd**2/(iw+xmu-ed0-sd)**2*gploc(i)
       Delta(1,i) =  tpd**2/(iw+xmu-d**2*gploc(i)/4.d0)
    enddo

    call fftgf_iw2tau(gploc,gptau,beta)
    nimp2=-2.d0*real(gptau(Ltau),8)
    ntot=nimp1+nimp2

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0d = iw + xmu - ed0 -delta_and(iw,epsiup,vup)
       sd  = g0d - one/Gwr(1,i)
       sp  = tpd**2/(iw + xmu - ed0 - sd)
       zita= iw + xmu -ep0 - sp
       gprloc(i)  = gfbether(wr(i),zita,D)
       gdrloc(i)  = one/(iw+xmu-ed0-sd) + tpd**2/(iw+xmu-ed0-sd)**2*gprloc(i)
    enddo
    call splot("Delta_iw.data",wm,delta(1,:),append=TT)
    call splot("locGdd_iw.data",wm,gdloc)
    call splot("locGpp_iw.data",wm,gploc)
    call splot("locGdd_realw.data",wr,gdrloc)
    call splot("locGpp_realw.data",wr,gprloc)
  end subroutine get_delta_bethe_pam
  !+----------------------------------------+

end program fullED



