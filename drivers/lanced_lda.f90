!This program solves the DMFT equations
!using complete ED at finite (high) T 
!reading the Hamiltonian model H(k)
!from a file with the standard
!wannier90 form.
!IMPORTANT: The ED solver uses HFmode=.true. 
!by default this means that chemical potential 
!SHOULD NOT include the U/2 shift.
program fulled_lda
  USE DMFT_ED
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: i,ik,iorb,jorb,iloop,Lk,Nb
  integer                :: Norb_d,Norb_p,Nopd,Nineq
  logical                :: converged
  real(8)                :: kx,ky,foo,ntotal,npimp
  real(8),allocatable    :: wm(:),wr(:),tau(:)
  complex(8)             :: iw
  !variables for the model:
  character(len=32)      :: file
  character(len=10)      :: string
  integer                :: ntype
  real(8)                :: nobj
  logical                :: rerun,bool
  !Bath:
  real(8),allocatable    :: Bath(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  complex(8),allocatable :: Dconv(:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  logical                :: fbethe
  real(8)                :: wbath

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !call parse_cmd_variable(rerun,"RERUN",default=.false.)
  !ADD HERE THE RERUN STUFF, copy from _lda1b.f90 1-band case
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  call read_input("inputED.in")
  call parse_cmd_variable(file,"FILE",default="hkfile.in")
  call parse_cmd_variable(ntype,"NTYPE",default=0)
  call parse_cmd_variable(fbethe,"FBETHE",default=.false.)
  call parse_cmd_variable(wbath,"WBATH",default=1.d0)


  !Allocate grids: 
  allocate(wm(NL),wr(Nw),tau(0:Ltau))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)
  tau = linspace(0.d0,beta,Ltau+1)

  !Read the lda-hamiltonian:
  call read_hk(trim(file))
  if(Nopd/=Norb)call warning("ACTHUNG: Nopd != Norb. I hope you know what you are doing!")


  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,NL))


  !Setup solver
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
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nobj,converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains


  subroutine read_hk(file)
    character(len=*) :: file
    integer :: ik
    real(8) :: de,e
    open(50,file=file,status='old')
    read(50,*)Lk,Norb_d,Norb_p,Nineq,foo
    Nopd=Norb_d+Norb_p
    allocate(Hk(Nopd,Nopd,Lk))
    allocate(dos_wt(Lk))
    do ik=1,Lk
       read(50,"(3(F10.7,1x))")kx,ky,foo
       do iorb=1,Nopd
          read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Nopd)
       enddo
    enddo
    dos_wt=2.d0*Wbath/dfloat(Lk)
    if(fbethe)then
       de=2.d0*wbath/dfloat(Lk)
       do ik=1,Lk
          e = -wbath + dfloat(ik-1)*de
          dos_wt(ik)=dens_bethe(e,wbath)*de
       enddo
    endif
  end subroutine read_hk



  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT


  subroutine get_delta
    integer                            :: i,j,iorb,jorb
    complex(8)                         :: iw,zita(Norb)
    complex(8),dimension(Norb,Norb)    :: fg
    complex(8),dimension(Norb,Norb,NL) :: g0d,sd,gloc,kdelta
    complex(8),dimension(Norb,Norb,Nw) :: g0dr,sdr,glocr,kdeltar



    Kdelta=0.d0;forall(iorb=1:Norb,i=1:NL)Kdelta(iorb,iorb,i)=xi*wm(i)+xmu
    do iorb=1,Norb
       do jorb=1,Norb
          do i=1,NL
             iw=xi*wm(i)
             g0d(iorb,jorb,i)= kdelta(iorb,jorb,i) - delta_and(iw,bath,iorb,jorb,1)
          enddo
       enddo
    enddo
    do i=1,NL
       fg = impGmats(:,:,1,i)
       call matrix_inverse(fg)
       sd(:,:,i) = g0d(:,:,i) - fg(:,:)
    enddo
    do i=1,NL
       iw     = xi*wm(i)
       do iorb=1,Norb
          zita(iorb)= kdelta(iorb,iorb,i) - sd(iorb,iorb,i)
       enddo
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i)  = fg(:,:)
    enddo
    !
    delta=zero
    do i=1,NL
       fg = gloc(:,:,i)
       call matrix_inverse(fg)
       delta(:,:,i) = kdelta(:,:,i) - fg(:,:) - sd(:,:,i)
    enddo
    !

    ! npimp=get_density_fromFFT(gp,beta)
    ntotal=sum(nimp)
    !write(*,*)"np  =",npimp
    write(*,*)"ntot=",ntotal


    Kdeltar=0.d0;forall(iorb=1:Norb,i=1:Nw)Kdeltar(iorb,iorb,i)=cmplx(wr(i),eps,8)+xmu
    do iorb=1,Norb
       do jorb=1,Norb
          do i=1,Nw
             iw=cmplx(wr(i),eps,8)
             g0dr(iorb,jorb,i)= kdeltar(iorb,jorb,i) - delta_and(iw,bath,iorb,jorb,1)
          enddo
       enddo
    enddo
    do i=1,Nw
       fg = impGreal(:,:,1,i)
       call matrix_inverse(fg)
       sdr(:,:,i) = g0dr(:,:,i) - fg(:,:)
    enddo
    do i=1,Nw
       do iorb=1,Norb
          zita(iorb)= kdeltar(iorb,iorb,i) - sdr(iorb,iorb,i)
       enddo
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       glocr(:,:,i)  = fg(:,:)
    enddo

    !Print:
    do iorb=1,Norb
       do jorb=1,Norb
          string=reg(txtfy(iorb))//reg(txtfy(jorb))
          call splot("Delta_"//reg(string)//"_iw.ed",wm,delta(iorb,jorb,:))
          call splot("G_"//reg(string)//"_iw.ed",wm,gloc(iorb,jorb,:))
          call splot("Sigma_"//reg(string)//"_iw.ed",wm,sd(iorb,jorb,:))
          call splot("DOS_"//reg(string)//".ed",wr,-dimag(glocr(iorb,jorb,:))/pi)
          call splot("Sigma_"//reg(string)//"_realw.ed",wr,sdr(iorb,jorb,:))
       enddo
    enddo

    ! call splot("np.ntot.ed.all",npimp,ntotal,append=.true.)
    ! call splot("np.ntot.ed",npimp,ntotal)

    if(ntype==1)then
       nobj=nimp(1)
    elseif(ntype==2)then
       nobj=nimp(2)
    else
       nobj=sum(nimp)
    endif

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



