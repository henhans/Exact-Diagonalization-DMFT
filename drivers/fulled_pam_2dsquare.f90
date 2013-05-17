program fullED
  USE DMFT_FULLED
  USE SQUARE_LATTICE
  USE FFTGF
  implicit none
  integer                :: i,iloop
  logical                :: converged
  real(8)                :: x,kint,eint,dw
  real(8)                :: gzero,gzerop,gzerom,gmu,ntot,npimp
  real(8),allocatable    :: wm(:),wr(:),tau(:),dos2d(:)
  !variables for the model:
  real(8)                :: ts,wband,ep0,tpd
  !The local hybridization function:
  complex(8),allocatable :: Delta(:)
  !
  namelist/vars/ts,ep0,tpd

  !Read inputs:
  call read_input("inputED.in")
  ts=1.d0
  ep0=0.d0
  tpd=0.4d0
  open(10,file="inputED.in")
  read(10,nml=vars)
  close(10) 
  call parse_cmd_variable(ts,"TS")
  call parse_cmd_variable(tpd,"TPD")
  call parse_cmd_variable(ep0,"EP0")
  write(*,nml=vars)

  !this shift contain |ep0-ed0|
  gmu=xmu
  gzerop=0.5d0*(ep0 + sqrt(ep0**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0 - sqrt(ep0**2 + 4.d0*tpd**2))
  if(ep0 < 0.d0)gzero=gzerop
  if(ep0 > 0.d0)gzero=gzerom
  if(ep0/= 0.d0)xmu=gmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  !Allocate:
  allocate(delta(NL))
  allocate(wm(NL),wr(Nw),tau(0:Ltau))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw,mesh=dw)
  tau = linspace(0.d0,beta,Ltau+1)

  !Analytic 2Dsquare DOS:
  wband=4.d0*ts
  allocate(dos2d(Nw))
  do i=1,Nw
     x=(wr(i)/abs(wfin-wini))**2-1.d0
     call comelp(x,kint,eint)
     dos2d(i)=2.d0/wband/pi**2*kint*heaviside(wband-abs(wr(i)))
  enddo
  dos2d=dos2d/sum(dos2d)/dw
  call splot("dos2d.ipt",wr,dos2d)
  dos2d=dos2d*dw

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
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(:),ebath(1,:),vbath(1,:))

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo

  !finalize calculation
  call ed_solver(status=-1)


contains


  !+----------------------------------------+
  subroutine get_delta
    integer                   :: i,j
    complex(8)                :: iw,zita
    complex(8),dimension(NL)  :: gp,gd,sp,sd,g0d
    complex(8),dimension(Nw)  :: gpr,gdr,spr,sdr,g0dr
    real(8),dimension(0:Ltau) :: gptau,gdtau

    delta=zero
    do i=1,NL
       iw     = xi*wm(i)
       g0d(i) = iw + xmu - delta_and(iw,ebath(1,:),vbath(1,:))
       sd(i)  = g0d(i) - one/Giw(1,i)
       sp(i)  = tpd**2/(iw + xmu - sd(i))
       zita   = iw + xmu -ep0 - sp(i)
       gp(i)  = sum_overk_zeta(zita,wr,dos2d)
       gd(i)  = one/(iw+xmu-sd(i)) + tpd**2/(iw+xmu-sd(i))**2*gp(i)
       !
       delta(i) = iw+xmu-one/gd(i)-sd(i)
       !
    enddo
    !
    call fftgf_iw2tau(gp,gptau,beta)
    call fftgf_iw2tau(gd,gdtau,beta)
    npimp=-2.d0*gptau(Ltau)
    ntot=nsimp+npimp

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0dr(i) = wr(i) + xmu - delta_and(wr(i)+zero,ebath(1,:),vbath(1,:))
       sdr(i)  = g0dr(i) - one/Gwr(1,i)
       spr(i)  = tpd**2/(iw + xmu - sdr(i))
       zita    = iw + xmu -ep0 - spr(i)
       gpr(i)  =  sum_overk_zeta(zita,wr,dos2d)
       gdr(i)  = one/(iw+xmu-sdr(i)) + tpd**2/(iw+xmu-sdr(i))**2*gpr(i)
    enddo
    call splot("Delta_iw.ed",wm,delta)
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("G_tau_ddpp.ed",tau,gdtau,gptau)
    call splot("G0dd_realw.ed",wr,one/g0dr)
    call splot("Gdd_realw.ed",wr,gdr)
    call splot("Gpp_realw.ed",wr,gpr)
    call splot("DOSdd.ed",wr,-dimag(gdr)/pi)
    call splot("DOSpp.ed",wr,-dimag(gpr)/pi)
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)
    call splot("Sigmapp_realw.ed",wr,spr)
    call splot("Sigmadd_realw.ed",wr,sdr)
    call splot("np.ntot.ed",npimp,ntot,append=TT)
  end subroutine get_delta
  !+----------------------------------------+

end program fullED



