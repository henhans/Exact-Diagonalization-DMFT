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
  gmu=xmu ;gzero=0.d0
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
     delta(1,:)=(d**2/4.d0)*G2iw(1,:) !the p=2 function!!
     call chi2_fitgf(delta(1,:),epsiup,vup)
     epsidw=epsiup;vdw=vup
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     call end_loop
  enddo
  call ed_solver()

end program fullED



