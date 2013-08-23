program fullED
  USE DMFT_ED
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: iloop,Nb
  logical                :: converged
  real(8)                :: gzero,gzerop,gzerom,gmu,ntotal,npimp
  real(8)                :: wband,ep0,tpd
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  !The local hybridization function:
  real(8),allocatable    :: bath(:)
  complex(8),allocatable :: Delta(:)
  !
  namelist/vars/wband,ep0,tpd

  call read_input("inputED.in")
  call parse_cmd_variable(wband,"wband","D",default=1.d0)
  call parse_cmd_variable(tpd,"TPD",default=0.4d0)
  call parse_cmd_variable(ep0,"EP0",default=0.d0)


  allocate(wm(NL),wr(Nw),tau(0:Ltau))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)
  tau = linspace(0.d0,beta,Ltau+1)

  !this shift contain |ep0-ed0|
  gmu=xmu
  gzerop=0.5d0*(ep0 + sqrt(ep0**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0 - sqrt(ep0**2 + 4.d0*tpd**2))
  if(ep0 < 0.d0)gzero=gzerop
  if(ep0 > 0.d0)gzero=gzerom
  if(ep0/= 0.d0)xmu=gmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  allocate(delta(NL))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_lanc_ed_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call lanc_ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe_pam

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(:),bath(:),1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntotal,converged)
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
       sd(i)  = g0d(i) - one/impGmats(1,i)
       sp(i)  = tpd**2/(iw + xmu - sd(i))
       zita   = iw + xmu - ep0 - sp(i)
       gp(i)  = gfbethe(wm(i),zita,Wband)
       gd(i)  = one/(iw+xmu-sd(i)) + tpd**2/(iw+xmu-sd(i))**2*gp(i)
       !
       delta(i) =  iw+xmu-one/gd(i)-sd(i)
       !delta(i) =  tpd**2/(iw+xmu-ep0-wband**2*gp(i)/4.d0)
       !
    enddo
    call fftgf_iw2tau(gp,gptau,beta)
    call fftgf_iw2tau(gd,gdtau,beta)
    npimp=-2.d0*real(gptau(Ltau),8)
    ntotal=nimp+npimp

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0dr(i) = wr(i)+xmu-delta_and(wr(i)+zero,bath,1)
       sdr(i)  = g0dr(i) - one/impGreal(1,i)
       spr(i)  = tpd**2/(iw + xmu - sdr(i))
       zita    = iw + xmu -ep0 - spr(i)
       gpr(i)  = gfbether(wr(i),zita,wband)
       gdr(i)  = one/(iw+xmu-sdr(i)) + tpd**2/(iw+xmu-sdr(i))**2*gpr(i)
    enddo
    !Print:
    call splot("Delta_iw.ed",wm,delta(:),append=.true.)
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
    call splot("np.ntot.ed",iloop,npimp,ntotal,gmu,append=.true.)
  end subroutine get_delta_bethe_pam
  !+----------------------------------------+

end program fullED



