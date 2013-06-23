program fullED
  USE DMFT_ED
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: iloop,Nb
  logical                :: converged
  real(8)                :: gzero,gzerop,gzerom,gmu,npimp
  real(8)                :: wband
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  !The local hybridization function:
  !Bath:
  real(8),allocatable    :: Bath(:)
  complex(8),allocatable :: Delta(:,:,:)
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
  call init_lanc_ed_solver(bath)

  !allocate delta function
  allocate(delta(Norb,Norb,NL))

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call lanc_ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_pam_2imp() 

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(2,2,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nobj,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains

  subroutine get_delta_pam_2imp
    integer                             :: i,j
    complex(8)                          :: iw,zita   
    complex(8),dimension(:),allocatable :: g0p,g0d,sd,sp,gp,gd
    real(8)                             :: zdd,zpp

    delta=zero
    allocate(g0p(NL),g0d(NL),sd(NL),sp(NL),gp(NL),gd(NL))
    do i=1,NL
       iw     = xi*wm(i)
       g0p(i) = iw + xmu - ep0 - delta_and(iw,bath,2,2,1)
       g0d(i) = iw + xmu - (tpd**2)/(iw + xmu - ep0 - 0.25d0*wband**2*impGmats(2,2,1,i))
       sd(i)  = g0d(i) - one/impGmats(1,1,1,i)
       sp(i)  = tpd**2/(iw + xmu - sd(i))
       zita   = iw + xmu - ep0 - sp(i)
       gp(i)  = gfbethe(wm(i),zita,Wband)
       gd(i)  = one/(iw+xmu-sd(i)) + tpd**2/(iw+xmu-sd(i))**2*gp(i)
       !
       delta(2,2,i) = 0.25d0*Wband**2*impGmats(2,2,1,i)
       !
    enddo
    call splot("Delta_iw.ed",wm,delta)
    call splot("G0dd_iw.ed",wm,one/g0d)
    call splot("G0pp_iw.ed",wm,one/g0p)
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)

    zdd=1.d0+abs(sd(1))/wm(1);zdd=1.d0/zdd
    zpp=1.d0+abs(sp(1))/wm(1);zpp=1.d0/zpp
    call splot("ntot.zdd.zpp.ed",iloop,nimp(1)+nimp(2),zdd,zpp,gmu,append=.true.)
    deallocate(g0p,g0d,sd,sp,gp,gd)

    ! allocate(g0p(Nw),g0d(Nw),sd(Nw),sp(Nw),gp(Nw),gd(Nw))
    ! do i=1,Nw
    !    iw     = wr(i)+xi*eps
    !    g0p(i) = iw + xmu - ep0 - delta_and(iw,bath,2,2,1)
    !    g0d(i) = iw + xmu - (tpd**2)/(iw + xmu - ep0 - 0.25d0*wband**2*impGreal(2,2,1,i))
    !    sd(i)  = g0d(i) - one/impGreal(1,1,1,i)
    !    sp(i)  = tpd**2/(iw + xmu - sd(i))
    !    zita   = iw + xmu - ep0 - sp(i)
    !    gp(i)  = gfbethe(wr(i),zita,Wband)
    !    gd(i)  = one/(iw+xmu-sd(i)) + tpd**2/(iw+xmu-sd(i))**2*gp(i)
    ! enddo
    ! call splot("G0dd_realw.ed",wr,one/g0d)
    ! call splot("G0pp_realw.ed",wr,one/g0p)
    ! call splot("Gdd_realw.ed",wr,gd)
    ! call splot("Gpp_realw.ed",wr,gp)
    ! call splot("DOSdd.ed",wr,-dimag(gd)/pi)
    ! call splot("DOSpp.ed",wr,-dimag(gp)/pi)
    ! call splot("Sigmapp_realw.ed",wr,sp)
    ! call splot("Sigmadd_realw.ed",wr,sd)
    ! deallocate(g0p,g0d,sd,sp,gp,gd)

    if(ntype==1)then
       nobj=nimp(1)
    elseif(ntype==2)then
       nobj=nimp(2)
    else
       nobj=nimp(1)+nimp(2)
    endif
  end subroutine get_delta_pam_2imp

end program fullED



