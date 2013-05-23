!This program solves the DMFT equations
!using complete ED at finite (high) T 
!reading the Hamiltonian model H(k)
!from a file with the standard
!wannier90 form.
!IMPORTANT:
!The ED solver uses HFmode=.true. by default
!this means that chemical potential SHOULD NOT
!include the U/2 shift.
program fulled_lda
  USE DMFT_FULLED
  USE FFTGF
  implicit none
  integer                :: i,ik,iorb,jorb,iloop,Lk,Nb
  integer                :: Norb_d,Norb_p,Norb,Nineq
  logical                :: converged
  real(8)                :: kx,ky,foo,ntot,npimp
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  complex(8)             :: iw
  !variables for the model:
  character(len=32)      :: file
  !Bath:
  real(8),allocatable    :: Bath(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: fg0(:,:,:)


  call read_input("inputED.in")
  call parse_cmd_variable(file,"FILE",default="hkfile.in")

  !Allocate:
  allocate(wm(NL),wr(Nw),tau(0:Ltau))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)
  tau = linspace(0.d0,beta,Ltau+1)

  open(50,file=trim(file),status='old')
  read(50,*)Lk,Norb_d,Norb_p,Nineq,foo
  Norb=Norb_d+Norb_p
  allocate(Hk(Norb,Norb,Lk))
  do ik=1,Lk
     read(50,"(3(F10.7,1x))")kx,ky,foo
     do iorb=1,Norb
        read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Norb)
     enddo
  enddo

  !Check non-interacting DOS:
  ! allocate(fg0(NL,Norb,Norb))
  ! do ik=1,Lk
  !    do i=1,NL
  !       iw = cmplx(wr(i),eps,8)+xmu
  !       fg0(i,:,:)=fg0(i,:,:) - dimag(inverse_g0k(iw,Hk(:,:,ik)))/pi/real(Lk,8)
  !    enddo
  ! enddo
  ! call splot("dosU0.ed",wr,fg0(:,1,1),fg0(:,2,2))
  ! deallocate(fg0)


  !Allocate Weiss Field:
  allocate(delta(NL))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_ed_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(:),bath,ichan=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(ntot,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains

  subroutine get_delta
    integer                   :: i,j
    complex(8)                :: iw,zita(2),fg(Norb,Norb)
    complex(8),dimension(NL)  :: gp,gd,sp,sd,g0d
    complex(8),dimension(Nw)  :: gpr,gdr,spr,sdr,g0dr
    real(8),dimension(0:Ltau) :: gptau,gdtau
    delta=zero
    do i=1,NL
       iw     = xi*wm(i)
       g0d(i) = iw + xmu - delta_and(iw,bath,1)
       sd(i)  = g0d(i) - one/Giw(1,i)      
       zita(1)= iw + xmu - sd(i)
       zita(2)= iw + xmu 
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))
       enddo
       gp(i)  = fg(2,2)/dble(Lk)
       gd(i)  = fg(1,1)/dble(Lk)
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
       g0dr(i) = wr(i) + xmu - delta_and(wr(i)+zero,bath,1)
       sdr(i)  = g0dr(i) - one/Gwr(1,i)
       zita(1) = iw + xmu - sdr(i)
       zita(2) = iw + xmu 
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))
       enddo
       gpr(i)  = fg(2,2)/dble(Lk) 
       gdr(i)  = fg(1,1)/dble(Lk)
    enddo
    !Print:
    call splot("Delta_iw.ed",wm,delta(:))
    call splot("Gdd_iw.ed",wm,gd)
    call splot("Gpp_iw.ed",wm,gp)
    call splot("Sigmadd_iw.ed",wm,sd)

    call splot("G0dd_realw.ed",wr,one/g0dr)
    call splot("Gdd_realw.ed",wr,gdr)
    call splot("Gpp_realw.ed",wr,gpr)
    call splot("DOSdd.ed",wr,-dimag(gdr)/pi)
    call splot("DOSpp.ed",wr,-dimag(gpr)/pi)

    call splot("Sigmadd_realw.ed",wr,sdr)

    call splot("G_tau_ddpp.ed",tau,gdtau,gptau)
    call splot("np.ntot.ed.all",npimp,ntot,append=.true.)
    call splot("np.ntot.ed",npimp,ntot)
  end subroutine get_delta



  function inverse_g0k(iw,hk) result(g0k)
    integer                     :: i,M
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k

  function inverse_gk(zeta,hk) result(gk)
    integer                     :: i,M
    complex(8),dimension(2,2)   :: hk
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

end program fulled_lda



