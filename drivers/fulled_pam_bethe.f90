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
  call init_ed_solver(bath)

  !allocate delta function
  allocate(delta(NL))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe_pam

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(:),bath,ichan=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nobj,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains


  !+----------------------------------------+
  subroutine get_delta_bethe_pam
    integer                   :: i,j
    complex(8)                :: iw,zita
    complex(8),dimension(NL)  :: gp,gd,sp,sd,g0d
    complex(8),dimension(Nw)  :: gpr,gdr,spr,sdr,g0dr
    real(8),dimension(0:Ltau) :: gptau,gdtau

    do i=1,NL
       iw     = xi*wm(i)
       g0d(i) = iw + xmu - delta_and(iw,bath,1)
       sd(i)  = g0d(i) - one/Giw(1,i)
       sp(i)  = tpd**2/(iw + xmu - sd(i))
       zita   = iw + xmu - ep0 - sp(i)
       gp(i)  = gfbethe(wm(i),zita,Wband)
       gd(i)  = one/(iw+xmu-sd(i)) + tpd**2/(iw+xmu-sd(i))**2*gp(i)
       !
       delta(i) =  iw+xmu-one/gd(i)-sd(i)
       !delta(i) =  tpd**2/(iw+xmu-ep0-d**2*gp(i)/4.d0)
       !
    enddo
    call fftgf_iw2tau(gp,gptau,beta)
    call fftgf_iw2tau(gd,gdtau,beta)
    npimp=-2.d0*real(gptau(Ltau),8)
    ntot=nsimp+npimp

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0dr(i) = wr(i) + xmu - delta_and(iw,bath,1)!delta_and(wr(i)+zero,bath,1)
       sdr(i)  = g0dr(i) - one/Gwr(1,i)
       spr(i)  = tpd**2/(iw + xmu - sdr(i))
       zita    = iw + xmu -ep0 - spr(i)
       gpr(i)  = gfbether(wr(i),zita,wband)
       gdr(i)  = one/(iw+xmu-sdr(i)) + tpd**2/(iw+xmu-sdr(i))**2*gpr(i)
    enddo
    !Print:
    call splot("Delta_iw.ed",wm,delta(:))
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    !
    call splot("G0dd_realw.ed",wr,one/g0dr)
    call splot("Gdd_realw.ed",wr,gdr)
    call splot("Gpp_realw.ed",wr,gpr)
    call splot("DOSdd.ed",wr,-dimag(gdr)/pi)
    call splot("DOSpp.ed",wr,-dimag(gpr)/pi)
    !
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)
    call splot("Sigmapp_realw.ed",wr,spr)
    call splot("Sigmadd_realw.ed",wr,sdr)
    !
    call splot("G_tau_ddpp.ed",tau,gdtau,gptau)
    call splot("np.ntot.ed",iloop,npimp,ntot,gmu,append=.true.)

    if(ntype==1)then
       nobj=nsimp
    elseif(ntype==2)then
       nobj=npimp
    else
       nobj=ntot
    endif
  end subroutine get_delta_bethe_pam
  !+----------------------------------------+

end program fullED



