program fullED
  USE DMFT_ED
  USE FFTGF
  implicit none
  integer                :: iloop,Nb
  logical                :: converged
  real(8)                :: gzero,gzerop,gzerom,gmu,ntot,npimp
  real(8)                :: wband
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  !The local hybridization function:
  !Bath:
  real(8),allocatable    :: Bath(:)
  complex(8),allocatable :: Delta(:)
  integer                :: ntype
  real(8)                :: nobj

  call read_input("inputED.in")
  call parse_cmd_variable(wband,"WBAND",default=1.d0)
  call parse_cmd_variable(ntype,"NTYPE",default=0)


  allocate(wm(NL),wr(Nw),tau(0:Ltau))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)
  tau = linspace(0.d0,beta,Ltau+1)

  !this shift contain |ep0-ed0|
  gzero=0.d0
  gmu=xmu
  gzerop=0.5d0*(ep0 + sqrt(ep0**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0 - sqrt(ep0**2 + 4.d0*tpd**2))
  if(ep0 < 0.d0)gzero=gzerop
  if(ep0 > 0.d0)gzero=gzerom
  if(ep0/= 0.d0)xmu=gmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_lanc_solver(bath)

  !allocate delta function
  allocate(delta(NL))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call lanc_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_pam_2imp() 

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(:),bath,ichan=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nobj,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains

  subroutine get_delta_pam_2imp
    integer                   :: i,j
    complex(8)                :: iw,zita
    complex(8),dimension(NL)  :: g0p,g0d,sd,sp,gp,gd
    complex(8),dimension(Nw)  :: g0pr,g0dr,sdr,spr,gpr,gdr
    real(8),dimension(0:Ltau) :: gptau,gdtau
    real(8)                   :: zdd,zpp

    do i=1,NL
       iw     = xi*wm(i)
       g0p(i) = iw + xmu - ep0 - delta_and(iw,bath,1)
       g0d(i) = iw + xmu - (tpd**2)/(iw + xmu - ep0 - 0.25d0*wband**2*G2iw(1,i))
       sd(i)  = g0d(i) - one/Giw(1,i)
       sp(i)  = tpd**2/(iw + xmu - sd(i))
       zita   = iw + xmu - ep0 - sp(i)
       gp(i)  = gfbethe(wm(i),zita,Wband)
       gd(i)  = one/(iw+xmu-sd(i)) + tpd**2/(iw+xmu-sd(i))**2*gp(i)
       !
       delta(i) = 0.25*Wband**2*G2iw(1,i)
       !
    enddo
    call splot("Delta_iw.ed",wm,delta)
    call splot("G0dd_iw.ed",wm,one/g0d)
    call splot("G0pp_iw.ed",wm,one/g0p)
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)

    ntot=nsimp+nimp2
    zdd=1.d0+abs(sd(1))/wm(1);zdd=1.d0/zdd
    zpp=1.d0+abs(sp(1))/wm(1);zpp=1.d0/zpp
    call splot("ntot.zdd.zpp.ed",iloop,ntot,zdd,zpp,gmu,append=.true.)

    if(ntype==1)then
       nobj=nsimp
    elseif(ntype==2)then
       nobj=npimp
    else
       nobj=ntot
    endif
  end subroutine get_delta_pam_2imp

end program fullED



