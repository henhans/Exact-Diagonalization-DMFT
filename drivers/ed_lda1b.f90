!This program solves the DMFT equations
!using complete ED at finite (high) T 
!reading the Hamiltonian model H(k)
!from a file with the standard
!wannier90 form.
!IMPORTANT: The ED solver uses HFmode=.true. 
!by default this means that chemical potential 
!SHOULD NOT include the U/2 shift.
program ed_lda1b
  USE DMFT_ED
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: i,ik,iorb,jorb,iloop,Lk
  integer                :: Norb_d,Norb_p,Nopd,Nineq
  logical                :: converged
  real(8)                :: kx,ky,foo,ntotal,npimp
  !Bath:
  integer,dimension(2)   :: Nb
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  character(len=32)      :: file
  integer                :: ntype
  real(8)                :: nobj
  logical                :: rerun,bool
  logical                :: fbethe
  real(8)                :: whband


  !parse additional variables
  call parse_cmd_variable(rerun,"RERUN",default=.false.)
  call parse_cmd_variable(file,"FILE",default="hkfile.in")
  call parse_cmd_variable(ntype,"NTYPE",default=0)
  call parse_cmd_variable(fbethe,"FBETHE",default=.false.)
  call parse_cmd_variable(whband,"WHBAND",default=1.d0)

  !+++++++++++++++++RERUN MODE+++++++++++++++++++++++++++++++++++++++
  if(rerun)then
     write(*,*)"+RERUN-mode: solve one and exit"
     call read_input("used.inputED.in")
     print*,"mu=",xmu
     inquire(file=trim(Hfile),exist=BOOL)
     if(.not.BOOL)then
        print*,"can't find ",trim(Hfile),". Exiting.."
        stop
     endif
     call read_hk('hkfile.in')
     !setup solver
     Nb=get_bath_size()
     allocate(bath(Nb(1),Nb(2)))
     call init_ed_solver(bath)
     call ed_solver(bath) 
     allocate(delta(Norb,NL))
     call get_delta
     stop
  endif
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  call read_input("inputED.in")

  !Read Hamiltoanian H(k)
  call read_hk(trim(file))

  !Allocate Weiss Field:
  allocate(delta(Norb,NL))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.OR.iloop>nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta,bath,ispin=1,iverbose=.true.)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nobj,converged)
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
    dos_wt=1.d0/dble(Lk)
    if(fbethe)then
       de=2.d0*WHband/dble(Lk)
       do ik=1,Lk
          e = -WHband + dfloat(ik-1)*de
          dos_wt(ik)=dens_bethe(e,WHband)*de
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
    integer                   :: i,j
    complex(8)                :: iw,zita(2),fg(Nopd,Nopd)
    complex(8),dimension(Nopd,NL)  :: gloc
    complex(8),dimension(Nopd,Nw)  :: grloc
    real(8)                   :: wm(NL),wr(Nw)

    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)

    delta=zero
    do i=1,NL
       iw     = xi*wm(i)
       zita(1)= iw + xmu - impSmats(1,1,1,i)
       zita(2)= iw + xmu 
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(1,i)  = fg(1,1)
       gloc(2,i)  = fg(2,2)
       !
       delta(1,i) = zita(1) - one/gloc(1,i)
       !
    enddo
    !Print:
    call splot("Delta_iw.ed",wm,delta(1,:))
    call splot("Gloc_iw.ed",wm,gloc(1,:),gloc(2,:))

    npimp=get_density_fromFFT(gloc(2,:),beta)
    ntotal=nimp(1)+npimp
    write(*,"(A,F25.18)")"np  =",npimp
    write(*,"(A,F25.18)")"ntot=",ntotal
    !Print:
    call splot("np.ntot_all.ed",npimp,ntotal,append=.true.)
    call splot("np.ntot.ed",npimp,ntotal)


    do i=1,Nw
       iw=dcmplx(wr(i),eps)
       zita(1) = iw + xmu - impSreal(1,1,1,i)
       zita(2) = iw + xmu 
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       grloc(1,i)  = fg(1,1)
       grloc(2,i)  = fg(2,2)
    enddo
    !Print:
    call splot("Gloc_realw.ed",wr,grloc(1,:),grloc(2,:))
    call splot("DOS.ed",wr,-dimag(grloc(1,:))/pi,-dimag(grloc(2,:))/pi)

    if(ntype==1)then
       nobj=nimp(1)
    elseif(ntype==2)then
       nobj=npimp
    else
       nobj=ntotal
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

  include 'search_mu.f90'
end program ed_lda1b



