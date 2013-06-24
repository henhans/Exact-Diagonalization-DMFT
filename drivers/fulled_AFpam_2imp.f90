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
  complex(8),allocatable,dimension(:),save :: dold

  call read_input("inputED.in")

  !this shift contain |ep0-ed0| already so you don't need to add/remove it!
  gmu=xmu ;gzero=0.d0
  if(ed0 > 0.d0)gzero=0.5*(ep0+ed0+sqrt((ep0-ed0)**2 + 4*tpd**2))
  if(ed0 < 0.d0)gzero=0.5*(ep0+ed0-sqrt((ep0-ed0)**2 + 4*tpd**2))
  if(ed0 /=0.d0)xmu=gmu+gzero !true ED chemical potential
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero


  call ed_solver(status=-2)

  iloop=0;converged=.false.
  allocate(dold(NL));dold=zero
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(status=0)

     !Get the Weiss field/Delta function to be fitted (user defined)

     call get_delta_pam_2imp() 

     !Fit the new bath, starting from the old bath + the supplied delta
     if(iloop>1)delta(1,:) = weigth*delta(1,:) + (1.d0-weigth)*dold
     dold=delta(1,:)

     call chi2_fitgf(delta(1,:),epsiup,vup)
     epsidw=epsiup;vdw=vup
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     call end_loop
  enddo

  call ed_solver(status=-1)


contains

  subroutine get_delta_pam_2imp
    integer :: i
    complex(8) :: iw,zeta
    complex(8) :: g0p(NL),g0d(NL),sd(NL),sp(NL),gp(NL),gd(NL)
    complex(8) :: g0pr(Nw),g0dr(Nw),sdr(Nw),spr(Nw),gpr(Nw),gdr(Nw)
    real(8)    :: gptau(0:Ltau),gdtau(0:Ltau),zdd,zpp

    do i=1,NL
       iw = xi*wm(i)
       g0p(i) = iw + xmu - ep0 - delta_and(iw,epsiup,vup)
       g0d(i) = iw + xmu - ed0 - (tpd**2)/(iw + xmu - ep0 - 0.25*d**2*G2iw(1,i))
       sd(i)  = g0d(i) - one/Giw(1,i)
       sp(i)  = (tpd**2)/(iw + xmu - ed0 - sd(i))
       zeta   = iw+xmu-ep0-sp(i)
       gp(i)  = gfbethe(wm(i),zeta,d)
       gd(i)  = one/(iw+xmu-ed0-sd(i)) + gp(i)*(tpd/(iw+xmu-ed0-sd(i)))**2
    enddo
    call splot("Delta_iw.ed",wm,delta(1,:),append=TT)
    call splot("G0dd_iw.ed",wm,one/g0d)
    call splot("G0pp_iw.ed",wm,one/g0p)
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)
    !
    delta(1,:) = 0.25*d**2*G2iw(1,:)
    !

    ntot=nimp1+nimp2
    zdd=1.d0+abs(sd(1))/wm(1);zdd=1.d0/zdd
    zpp=1.d0+abs(sp(1))/wm(1);zpp=1.d0/zpp
    call splot("ntot.zdd.zpp.ed",ntot,zdd,zpp,append=TT)
  end subroutine get_delta_pam_2imp

end program fullED



