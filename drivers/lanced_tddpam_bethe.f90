program lanctddpam
  USE DMFT_ED
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: iloop,Nb
  logical                :: converged
  real(8)                :: gzero,gzerop,gzerom,gmu,ntotal,npimp
  real(8)                :: wband,alpha
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  !The local hybridization function:
  real(8),allocatable    :: bath(:)
  complex(8),allocatable :: Delta(:,:,:)
  integer                :: ntype
  real(8)                :: nobj
  namelist/xtravars/wband,alpha,tpd,ep0

  call read_input("inputED.in")
  call parse_cmd_variable(wband,"WBAND","D",default=1.d0)
  call parse_cmd_variable(alpha,"ALPHA","A",default=0.25d0)
  open(100,file="used.xtravars.in")
  write(100,nml=xtravars)
  close(100)


  allocate(wm(NL),wr(Nw),tau(0:Ltau))
  wm  = pi/beta*real(2*arange(1,NL)-1,8)
  wr  = linspace(wini,wfin,Nw)
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

  allocate(delta(Norb,Norb,NL))

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
     call get_delta_bethe_tddpam

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta,bath(:),1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(Norb,Norb,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntotal,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains


  !+----------------------------------------+
  subroutine get_delta_bethe_tddpam
    integer                   :: i,j,ie
    complex(8)                :: iw,zita,gamma,foo(3)
    complex(8),dimension(NL)  :: g0d,sd,gp,gd,gpd
    complex(8),dimension(Nw)  :: g0dr,sdr,gpr,gdr,gpdr
    real(8),dimension(0:Ltau) :: gptau,gdtau
    real(8)                   :: dosb(Nw),eb(Nw),de

    de=2.d0*wband/dfloat(Nw)
    do i=1,Nw
       eb(i)  = -Wband + real(i-1,8)*de
       dosb(i)= dens_bethe(eb(i),Wband)*de
    enddo

    do i=1,NL
       iw     = xi*wm(i)
       g0d(i) = iw + xmu - delta_and(iw,bath,1,1,1)
       sd(i)  = g0d(i) - one/impGmats(1,1,1,i)
       foo=zero
       do ie=1,Nw
          gamma = iw+xmu-alpha*eb(ie)-sd(i)
          zita  = iw+xmu-ep0-eb(ie)
          foo(1) = foo(1) + dosb(ie)/(gamma - tpd**2/zita)
          foo(2) = foo(2) + dosb(ie)/(zita  - tpd**2/gamma)
          foo(3) = foo(3) + dosb(ie)/(gamma*zita - tpd**2)
       enddo
       gd(i)  = foo(1)
       gp(i)  = foo(2)
       gpd(i) = tpd*foo(3)
       !
       delta(Norb,Norb,i) =  iw+xmu-one/gd(i)-sd(i)
       !
    enddo

    npimp=get_density_fromFFT(gp,beta)
    ntotal=nimp(1)+npimp
    write(*,*)"np  =",npimp
    write(*,*)"ntot=",ntotal

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0dr(i) = wr(i)+xmu-delta_and(wr(i)+zero,bath,1,1,1)
       sdr(i)  = g0dr(i) - one/impGreal(1,1,1,i)
       foo=zero
       do ie=1,Nw
          gamma = iw+xmu-alpha*eb(ie)-sdr(i)
          zita  = iw+xmu-ep0-eb(ie)
          foo(1) = foo(1) + dosb(ie)/(gamma - tpd**2/zita)
          foo(2) = foo(2) + dosb(ie)/(zita  - tpd**2/gamma)
          foo(3) = foo(3) + dosb(ie)/(gamma*zita - tpd**2)
       enddo
       gdr(i)  = foo(1)
       gpr(i)  = foo(2)
       gpdr(i) = tpd*foo(3)
    enddo

    !Print:
    call splot("Delta_iw.ed",wm,delta(Norb,Norb,:),append=.true.)
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("Gpd_iw.ed",wm,gpd)
    !
    call splot("G0dd_realw.ed",wr,one/g0dr)
    call splot("Gdd_realw.ed",wr,gdr)
    call splot("Gpp_realw.ed",wr,gpr)
    call splot("Gpd_realw.ed",wr,gpdr)
    call splot("DOSdd.ed",wr,-dimag(gdr)/pi)
    call splot("DOSpp.ed",wr,-dimag(gpr)/pi)
    call splot("DOSpd.ed",wr,-dimag(gpdr)/pi)
    !
    call splot("Sigmadd_iw.ed",wm,sd)
    call splot("Sigmadd_realw.ed",wr,sdr)
    !
    call splot("np.ntot.ed",iloop,npimp,ntotal,gmu,append=.true.)
  end subroutine get_delta_bethe_tddpam
  !+----------------------------------------+


  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT

end program lanctddpam



