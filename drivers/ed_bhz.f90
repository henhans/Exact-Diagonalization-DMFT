!                    MODEL Hamiltonian is:
!
! |     h^{2x2}(k)              &         hso^{2x2}(k)        |
! |      [hso^{2x2}]*(k)        &        [h^{2x2}]*(-k)       |
!
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
! hso^{2x2}(k):=
! | xi*rh*(sin(kx)-xi*sin(ky))  &         \delta              |
! |         -\delta             &             0               |
program ed_bhz
  USE DMFT_ED
  USE FFTGF
  USE TOOLS
  USE MATRIX
  implicit none
  integer                :: i,ik,iorb,jorb,ispin,iloop,Lk
  integer                :: Norb_d,Lorb
  logical                :: converged
  real(8)                :: kx,ky,foo,ntotal,npimp
  !variables for the model:
  character(len=32)      :: file
  integer                :: ntype
  real(8)                :: nobj
  logical                :: rerun,bool
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: dos_wt(:),e0(:)

  !Parse additional variables
  call parse_cmd_variable(file,"FILE",default="hkfile.in")
  call parse_cmd_variable(ntype,"NTYPE",default=0)

  !Read input file:
  call read_input("inputED.in")
  if(Nspin/=2.OR.Norb/=2)stop "wrong setup from input file: Nspin=Norb=2"

  !Read Hamiltoanian H(k)
  call read_hk(trim(file))

  !Allocate Weiss Field:
  allocate(delta(Nspin,Norb,NL))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.OR.iloop>=nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(1,1:Norb,1:NL),bath,ispin=1)
     call chi2_fitgf(delta(2,1:Norb,1:NL),bath,ispin=2)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),eps_error,nsuccess,nloop)
     call end_loop
  enddo


contains


  subroutine read_hk(file)
    character(len=*)       :: file
    integer                :: ik
    real(8)                :: de,e
    complex(8),allocatable :: fg(:,:,:)
    complex(8)             :: iw
    real(8)                :: wr(Nw)

    wr = linspace(wini,wfin,Nw)

    open(50,file=file,status='old')
    read(50,*)Lk,Norb_d,foo,foo,foo
    if(Norb_d/=2)stop "Norb_d in Hk should be 2."
    Lorb=Nspin*Norb_d
    if(Lorb/=4)stop "# of SObands should be 4"
    allocate(Hk(Lorb,Lorb,Lk),dos_wt(Lk),e0(Lorb))
    do ik=1,Lk
       read(50,"(3(F10.7,1x))")kx,ky,foo
       do iorb=1,Lorb
          read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Lorb)
       enddo
    enddo
    dos_wt=1.d0/dble(Lk)
    write(*,*)"# of k-points:",Lk
    write(*,*)"# of SO-bands :",Lorb
    write(*,*)"Evaluate non-interacting G(w):"
    allocate(fg(Nw,Lorb,Lorb))
    fg=zero
    e0=0.d0
    do ik=1,Lk
       do i=1,Nw
          iw = cmplx(wr(i),eps,8)+xmu
          fg(i,:,:)=fg(i,:,:)+inverse_g0k(iw,Hk(:,:,ik))*dos_wt(ik)
       enddo
       do i=1,Lorb
          e0(i)=e0(i)+Hk(i,i,ik)*dos_wt(ik)
       enddo
    enddo
    print*,e0
    call splot("DOS0.ed",wr,-dimag(fg(:,1,1))/pi,-dimag(fg(:,2,2))/pi,&
         -dimag(fg(:,3,3))/pi,-dimag(fg(:,4,4))/pi)
    deallocate(fg)
  end subroutine read_hk



  subroutine get_delta
    integer                                 :: i,j,iorb,ispin
    complex(8)                              :: iw,zita(Lorb),fg(Lorb,Lorb)
    complex(8),dimension(:,:,:),allocatable :: gloc
    real(8)                                 :: wm(NL),wr(Nw)

    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)

    delta=zero

    print*,"Get Gloc_iw:"
    allocate(gloc(Nspin,Norb,NL))
    do i=1,NL
       iw = xi*wm(i)
       zita(1)= iw + xmu - impSmats(1,1,1,i)
       zita(2)= iw + xmu - impSmats(1,2,2,i)
       zita(3)= iw + xmu - impSmats(2,1,1,i)
       zita(4)= iw + xmu - impSmats(2,2,2,i)

       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(1,1,i)  = fg(1,1)
       gloc(1,2,i)  = fg(2,2)
       gloc(2,1,i)  = fg(3,3)
       gloc(2,2,i)  = fg(4,4)
       call matrix_inverse(fg)
       !
       delta(1,1,i) = zita(1)-fg(1,1)
       delta(1,2,i) = zita(2)-fg(2,2)
       delta(2,1,i) = zita(3)-fg(3,3)
       delta(2,2,i) = zita(4)-fg(4,4)
       !
    enddo
    call splot("Delta_iw.ed",wm,delta(1,1,:),delta(1,2,:),delta(2,1,:),delta(2,2,:))
    call splot("Gloc_iw.ed",wm,gloc(1,1,:),gloc(1,2,:),gloc(2,1,:),gloc(2,2,:))
    deallocate(gloc)

    allocate(gloc(Nspin,Norb,Nw))
    print*,"Get Gloc_realw:"
    do i=1,Nw
       iw=cmplx(wr(i),eps)
       zita(1)= iw + xmu - impSreal(1,1,1,i)
       zita(2)= iw + xmu - impSreal(1,2,2,i)
       zita(3)= iw + xmu - impSreal(2,1,1,i)
       zita(4)= iw + xmu - impSreal(2,2,2,i)
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(1,1,i)  = fg(1,1)
       gloc(1,2,i)  = fg(2,2)
       gloc(2,1,i)  = fg(3,3)
       gloc(2,2,i)  = fg(4,4)
    enddo
    call splot("Gloc_realw.ed",wr,gloc(1,1,:),gloc(1,2,:),gloc(2,1,:),gloc(2,2,:))
    call splot("DOS.ed",wr,-dimag(gloc(1,1,:))/pi,-dimag(gloc(1,2,:))/pi,&
         -dimag(gloc(2,1,:))/pi,-dimag(gloc(2,2,:))/pi)
    deallocate(gloc)


    ntotal=sum(nimp)
    write(*,*)"ntot=",ntotal
    if(ntype==1)then
       nobj=nimp(1)
    else
       nobj=ntotal
    endif

  end subroutine get_delta



  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k=zero
    g0k(1:2,1:2) = inverse_g0k2x2(iw,hk(1:2,1:2))
    g0k(3:4,3:4) = inverse_g0k2x2(iw,hk(3:4,3:4))
  end function inverse_g0k
  !
  function inverse_g0k2x2(iw,hk) result(g0k)
    integer                     :: i
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
  end function inverse_g0k2x2


  function inverse_gk(zeta,hk) result(gk)
    complex(8)                  :: zeta(4)
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: gk
    gk=zero
    gk(1:2,1:2) = inverse_gk2x2(zeta(1:2),hk(1:2,1:2))
    gk(3:4,3:4) = inverse_gk2x2(zeta(3:4),hk(3:4,3:4))
  end function inverse_gk
  !
  function inverse_gk2x2(zeta,hk) result(gk)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8),dimension(2)     :: zeta
    complex(8),dimension(2,2)   :: gk
    Complex(8)                  :: delta,ppi,vmix
    gk=zero
    delta = zeta(1) - hk(1,1)
    ppi   = zeta(2) - hk(2,2)
    vmix  = -hk(1,2)
    gk(1,1) = one/(delta - abs(vmix)**2/ppi)
    gk(2,2) = one/(ppi - abs(vmix)**2/delta)
    gk(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    gk(2,1) = conjg(gk(1,2))
  end function inverse_gk2x2


  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT

end program ed_bhz



