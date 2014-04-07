program ed_tddpam_lattice
  USE DMFT_ED
  USE COMMON_VARS
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  USE ERROR
  USE ARRAYS
  USE MATRIX
  USE IOTOOLS
  USE PARSE_INPUT
  implicit none
  integer                :: iloop
  logical                :: converged
  integer                :: Npd
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  complex(8),allocatable :: Delta(:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  character(len=32)      :: hkfile,finput
  integer                :: Nx,Lk,ntype
  real(8)                :: nobj
  real(8)                :: alpha,tpp,ep0,tpd,v0,gzero,gzerop,gzerom,gmu
  logical                :: bool

  !parse additional variables
  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(ntype,"NTYPE",trim(finput),default=0)
  call parse_input_variable(Nx,"NX",trim(finput),default=100)
  call parse_input_variable(tpp,"TPP",trim(finput),default=0.25d0)
  call parse_input_variable(alpha,"ALPHA",trim(finput),default=0.d0)
  call parse_input_variable(tpd,"TPD",trim(finput),default=0.d0)
  call parse_input_variable(v0,"V0",trim(finput),default=0.d0)
  call parse_input_variable(ep0,"EP0",trim(finput),default=0.d0)
  call ed_read_input(trim(finput))


  inquire(file="last_mu.restart",exist=bool)
  if(bool.AND.nread/=0.d0)then
     open(100,file="last_mu.restart")
     read(100,*)xmu
     close(100)
     print*,"UPDATE XMU:",xmu
  endif

  !this shift contain |ep0-ed0|
  gmu=xmu
  gzerop=0.5d0*(ep0 + sqrt(ep0**2 + 4.d0*tpd**2))
  gzerom=0.5d0*(ep0 - sqrt(ep0**2 + 4.d0*tpd**2))
  gzero=0.d0
  if(ep0 < 0.d0)gzero=gzerop
  if(ep0 > 0.d0)gzero=gzerom
  if(ep0/= 0.d0)xmu=gmu+gzero
  write(*,*)'shift mu to (from) = ',xmu,'(',gmu,')'
  write(*,*)'shift is           = ',gzero

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)

  !Read/Build the H(k)
  call build_hk()

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop)
     if(nread/=0.d0)call search_chemical_potential(nobj,converged)
     call end_loop
  enddo


  if(nread/=0.d0)then
     open(100,file='last_mu.restart')
     write(100,*)xmu-gzero
     close(100)
  endif


contains


  subroutine build_hk()
    integer                :: ix,iy,i,j,ik,iorb,jorb,unit
    real(8)                :: kx,ky
    real(8),allocatable    :: peloc(:)
    complex(8),allocatable :: pHloc(:,:)
    Npd=2
    Lk=Nx*Nx
    allocate(Hk(Npd,Npd,Lk),pHloc(Npd,Npd),peloc(Npd))
    allocate(dos_wt(Lk))
    unit=free_unit()
    open(unit,file="Eigenbands.ed")
    ik=0
    do ix=1,Nx
       kx=-pi + 2.d0*pi*dble(ix-1)/dble(nx)
       do iy=1,Nx
          ky=-pi + 2.d0*pi*dble(iy-1)/dble(Nx)
          ik=ik+1
          Hk(:,:,ik)=Hk_model(kx,ky)
       enddo
    enddo
    dos_wt = 1.d0/dble(Lk)
    ik=0
    do ix=1,100
       ik=ik+1
       kx = 0.d0 + pi*real(ix-1,8)/100.d0
       ky = 0.d0
       pHloc=Hk_model(kx,ky) 
       call matrix_diagonalize(pHloc,peloc)
       write(unit,*)ik,peloc(1),peloc(2)
    enddo
    !From X=(pi,0) to M=(pi,pi): 100 steps
    do iy=1,100
       ik=ik+1
       kx = pi
       ky = 0.d0 + pi*real(iy-1,8)/100.d0
       pHloc=Hk_model(kx,ky) 
       call matrix_diagonalize(pHloc,peloc)
       write(unit,*)ik,peloc(1),peloc(2)
    enddo
    !From M=(pi,pi) to \Gamma=(0,0): 100 steps
    do ix=1,100
       ik=ik+1
       iy=ix
       kx = pi - pi*real(ix-1,8)/100.d0
       ky = pi - pi*real(iy-1,8)/100.d0
       pHloc=Hk_model(kx,ky) 
       call matrix_diagonalize(pHloc,peloc)
       write(unit,*)ik,peloc(1),peloc(2)
    enddo
  end subroutine build_hk
  function Hk_model(kx,ky) result(Hk)
    real(8)                   :: kx,ky,epsik,vpsik
    complex(8),dimension(2,2) :: Hk
    epsik = cos(kx)+cos(ky)
    vpsik = sin(kx)*sin(ky)
    Hk(1,1) = Hloc(1,1,1,1) - 2.d0*alpha*tpp*epsik
    Hk(2,2) = ep0           - 2.d0*tpp*epsik
    Hk(1,2) = tpd - 4.d0*v0*vpsik
    Hk(2,1) = tpd - 4.d0*v0*vpsik
  end function Hk_model


  !+----------------------------------------+



  subroutine get_delta
    integer                                 :: i,j,ik,iorb,jorb
    complex(8)                              :: iw,zita(2),fg(Npd,Npd)
    complex(8),dimension(:,:,:),allocatable :: gloc
    real(8)                                 :: wm(Lmats),wr(Lreal),npimp,ntotal
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    delta=zero
    allocate(gloc(Npd,Npd,Lmats))
    do i=1,Lmats
       zita(1)= xi*wm(i)+xmu  - impSmats(1,1,1,1,i)
       zita(2)= xi*wm(i)+xmu
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i)  = fg
       !
       delta(1,1,i) = zita(1) - Hloc(1,1,1,1) - one/fg(1,1)
       !
    enddo
    !Print:
    call splot("Delta_iw.ed",wm,delta(1,1,:))
    do iorb=1,Npd
       do jorb=1,Npd
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
       enddo
    enddo
    npimp=get_density_fromFFT(gloc(2,2,:),beta)
    ntotal=ed_dens(1)+npimp
    write(*,"(A,F25.18)")"np  =",npimp
    write(*,"(A,F25.18)")"ntot=",ntotal
    call splot("np.ntot_all.ed",npimp,ntotal,append=.true.)
    call splot("np.ntot.ed",npimp,ntotal)
    deallocate(gloc)


    allocate(gloc(Npd,Npd,Lreal))
    do i=1,Lreal
       zita(1) = dcmplx(wr(i),eps)+xmu - impSreal(1,1,1,1,i)
       zita(2) = dcmplx(wr(i),eps)+xmu
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i)  = fg
    enddo
    !Print:
    do iorb=1,Npd
       do jorb=1,Npd
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,gloc(iorb,jorb,:))
          call splot("DOS_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,-dimag(gloc(iorb,jorb,:))/pi)
       enddo
    enddo
    deallocate(gloc)

    if(ntype==1)then
       nobj=ed_dens(1)
    elseif(ntype==2)then
       nobj=npimp
    else
       nobj=ntotal
    endif

  end subroutine get_delta


  !+----------------------------------------+


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

  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT

  function eplus(hk)
    complex(8),dimension(2,2) :: hk
    real(8)                   :: eplus
    eplus = hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eplus = eplus/2.d0
  end function eplus

  function eminus(hk)
    complex(8),dimension(2,2) :: hk
    real(8)                   :: eminus
    eminus = hk(1,1)+hk(2,2) -sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eminus = eminus/2.d0
  end function eminus

end program ed_tddpam_lattice



