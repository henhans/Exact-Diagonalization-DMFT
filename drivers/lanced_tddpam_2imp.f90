program fullED
  USE DMFT_ED
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: iloop,Nb
  logical                :: converged
  real(8)                :: gzero,gzerop,gzerom,gmu,npimp
  real(8)                :: wband,alpha
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  !The local hybridization function:
  !Bath:
  real(8),allocatable    :: Bath(:)
  complex(8),allocatable :: Delta(:,:,:)
  integer                :: ntype
  real(8)                :: nobj
  namelist/xtravars/wband,alpha,tpd,ep0

  call read_input("inputED.in")
  call parse_cmd_variable(ntype,"NTYPE",default=0)
  call parse_cmd_variable(wband,"WBAND","D",default=1.d0)
  call parse_cmd_variable(alpha,"ALPHA","A",default=0.25d0)
  open(100,file="used.xtravars.in")
  write(100,nml=xtravars)
  close(100)

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

  !allocate delta function
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
     call get_delta_tddpam_2imp() 

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:)+delta(2,2,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nobj,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains

  subroutine get_delta_tddpam_2imp
    integer                             :: i,j,iorb,jorb,ie
    complex(8)                          :: iw,zeta(Norb)   
    complex(8)                          :: calG0(Norb,Norb),Gloc(Norb,Norb),Sloc(Norb,Norb),Foo(Norb,Norb)
    complex(8),dimension(:),allocatable :: sd,sp,gp,gd,gpd
    real(8)                             :: zdd,zpp,Hk(Norb,Norb)
    real(8)                             :: dosb(Nw),eb(Nw),de

    de=2.d0*wband/dfloat(Nw)
    do i=1,Nw
       eb(i)  = -Wband + real(i-1,8)*de
       dosb(i)= dens_bethe(eb(i),Wband)*de
    enddo

    delta=zero
    allocate(sd(NL),sp(NL),gp(NL),gd(NL),gpd(NL))
    do i=1,NL
       iw   = xi*wm(i)
       calG0=zero
       do iorb=1,Norb
          calG0(iorb,iorb)=iw+xmu
          if(iorb==2)calG0(iorb,iorb)=calG0(iorb,iorb)-ep0
          do jorb=1,Norb
             calG0(iorb,jorb) = calG0(iorb,jorb)-delta_and(iw,bath,iorb,jorb,1)             
             Gloc(iorb,jorb)  = impGmats(iorb,jorb,1,i)
          enddo
       enddo
       call matrix_inverse(Gloc)
       Sloc = calG0 - Gloc
       Foo  = zero
       do ie=1,Nw
          Hk(1,1) = alpha*eb(ie)
          Hk(2,2) = ep0+eb(ie)
          Hk(1,2) = tpd
          Hk(2,1) = tpd
          zeta(1) = iw + xmu - Sloc(1,1)
          zeta(2) = iw + xmu - Sloc(2,2)
          Foo = Foo + inverse_gk(zeta,Hk)*dosb(ie)
       enddo
       gd(i) = Foo(1,1)
       gp(i) = Foo(2,2)
       gpd(i)= Foo(1,2)
       sd(i) = Sloc(1,1)
       sp(i) = Sloc(2,2)
       !
       delta(1,1,i) = Wband**2*gd(i)/4.d0
       delta(2,2,i) = (alpha*Wband)**2*gp(i)/4.d0
       delta(1,2,i) = alpha*Wband**2*gpd(i)/4.d0
       delta(2,1,i) = delta(1,2,i)
       !
    enddo

    !call splot("Delta_iw.ed",wm,delta)
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("Gpd_iw.ed",wm,gpd)
    call splot("Sigmapp_iw.ed",wm,sp)
    call splot("Sigmadd_iw.ed",wm,sd)

    ! zdd=1.d0+abs(sd(1))/wm(1);zdd=1.d0/zdd
    ! zpp=1.d0+abs(sp(1))/wm(1);zpp=1.d0/zpp
    ! call splot("ntot.zdd.zpp.ed",iloop,nimp(1)+nimp(2),zdd,zpp,gmu,append=.true.)
    deallocate(sd,sp,gp,gd,gpd)

    allocate(sd(Nw),sp(Nw),gp(Nw),gd(Nw),gpd(Nw))
    do i=1,Nw
       iw=cmplx(wr(i),eps)
       calG0=zero
       do iorb=1,Norb
          calG0(iorb,iorb)=wr(i)+xmu
          if(iorb==2)calG0(iorb,iorb)=calG0(iorb,iorb)-ep0
          do jorb=1,Norb
             calG0(iorb,jorb) = calG0(iorb,jorb)-delta_and(wr(i)+zero,bath,iorb,jorb,1)             
             Gloc(iorb,jorb)  = impGreal(iorb,jorb,1,i)
          enddo
       enddo
       call matrix_inverse(Gloc)
       Sloc = calG0 - Gloc
       Foo  = zero
       do ie=1,Nw
          Hk(1,1) = alpha*eb(ie)
          Hk(2,2) = ep0+eb(ie)
          Hk(1,2) = tpd
          Hk(2,1) = tpd
          zeta(1) = iw + xmu - Sloc(1,1)
          zeta(2) = iw + xmu - Sloc(2,2)
          Foo = Foo + inverse_gk(zeta,Hk)*dosb(ie)
       enddo
       gd(i) = Foo(1,1)
       gp(i) = Foo(2,2)
       gpd(i)= Foo(1,2)
       sd(i) = Sloc(1,1)
       sp(i) = Sloc(2,2)
    enddo
    ! call splot("Gdd_realw.ed",wr,gd)
    ! call splot("Gpp_realw.ed",wr,gp)
    call splot("DOSdd.ed",wr,-dimag(gd)/pi)
    call splot("DOSpp.ed",wr,-dimag(gp)/pi)
    call splot("DOSpp.ed",wr,-dimag(gpd)/pi)
    call splot("Sigmapp_realw.ed",wr,sp)
    call splot("Sigmadd_realw.ed",wr,sd)
    deallocate(sd,sp,gp,gd,gpd)

    if(ntype==1)then
       nobj=nimp(1)
    elseif(ntype==2)then
       nobj=nimp(2)
    else
       nobj=nimp(1)+nimp(2)
    endif
  end subroutine get_delta_tddpam_2imp


  function inverse_gk(zeta,hk) result(gk)
    integer                     :: i,M
    real(8),dimension(2,2)   :: hk
    complex(8),dimension(2)     :: zeta
    complex(8),dimension(2,2)   :: gk
    complex(8)                  :: delta,ppi,vmix
    gk=zero
    delta = zeta(1) - hk(1,1)
    ppi   = zeta(2) - hk(2,2)
    vmix  = -hk(1,2)
    gk(1,1) = one/(delta - abs(vmix)**2/ppi)
    gk(2,2) = one/(ppi - abs(vmix)**2/delta)
    gk(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    gk(2,1) = conjg(gk(1,2))
  end function inverse_gk

end program fullED



