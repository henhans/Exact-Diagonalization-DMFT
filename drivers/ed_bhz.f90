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
  USE TOOLS
  USE ARRAYS
  USE ERROR
  USE MATRIX
  implicit none
  integer                :: iloop,Lk,Nsporb
  logical                :: converged
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:),Fold(:,:,:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  character(len=16)      :: finput,fhloc
  character(len=32)      :: hkfile
  integer                :: ntype
  real(8)                :: nobj,wmixing
  logical                :: fnomag

#ifdef _MPI
  call ed_init_mpi()
#endif

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(hkfile,"HKFILE",default="hkfile.in")
  call parse_cmd_variable(ntype,"NTYPE",default=0)
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_cmd_variable(fhloc,"FHLOC",default='inputHLOC.in')
  call parse_cmd_variable(fnomag,"FNOMAG",default=.true.)
  call parse_cmd_variable(wmixing,"WMIXING",default=0.5d0)
  !
  call ed_read_input(trim(finput),trim(fhloc))


  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,NL))
  allocate(Fold(Nspin,Nspin,Norb,Norb,NL))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)

  !
  call read_hk(trim(hkfile))


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     if(iloop>1)delta=wmixing*delta + (1.d0-wmixing)*Fold
     Fold=Delta

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(1,1,:,:,:),bath,ispin=1,iverbose=.true.)
     if(fnomag)then
        call spin_symmetrize_bath(bath)
     else
        call chi2_fitgf(delta(2,2,:,:,:),bath,ispin=2)
     endif
     converged = check_convergence(delta(1,1,1,1,:)+delta(1,1,2,2,:),dmft_error,nsuccess,nloop)
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif

     call end_loop
  enddo

#ifdef _MPI
  call ed_finalize_mpi()
#endif


contains


  subroutine read_hk(file)
    character(len=*)       :: file
    integer                :: i,j,ik,iorb,jorb,Norb_d,isporb,jsporb,ispin,jspin
    real(8)                :: de,e,kx,ky,kz,foo
    real(8)                :: ntotal,npimp
    complex(8),allocatable :: fg(:,:,:)
    complex(8)             :: iw
    real(8)                :: wr(Nw)
    !
    open(50,file=file,status='old')
    read(50,*)Lk,Norb_d,foo,foo,foo
    if(Nspin/=2.OR.Norb/=2)stop "wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
    if(Norb_d/=Norb)stop "Norb_d in Hk should be 2."
    Nsporb=Nspin*Norb
    if(Nsporb/=4)stop "# of SObands should be 4"
    allocate(Hk(Nsporb,Nsporb,Lk),dos_wt(Lk),bhzHloc(Nsporb,Nsporb))
    do ik=1,Lk
       read(50,"(3(F10.7,1x))")kx,ky,foo
       do isporb=1,Nsporb
          read(50,"(10(2F10.7,1x))")(Hk(isporb,jsporb,ik),jsporb=1,Nsporb)
       enddo
    enddo
    dos_wt=1.d0/dble(Lk)
    bhzHloc = sum(Hk(:,:,:),dim=3)/dble(Lk)
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    write(*,*)"# of k-points :",Lk
    write(*,*)"# of SO-bands :",Nsporb
    Hloc = j2ls_func(bhzHloc)
    write(*,*)"Updated Hloc:"
    call print_Hloc(Hloc)
  end subroutine read_hk







  subroutine get_delta
    integer                                     :: i,j,ik,iorb,jorb,ispin,jspin
    complex(8),dimension(Nsporb,Nsporb)         :: zeta,fg,gdelta,self,gkinv
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc
    complex(8)                                  :: iw
    real(8)                                     :: wm(NL),wr(Nw)
    character(len=20)                           :: suffix
    !
    write(*,*)"WARNING: this version of the code is spin-symmetric"
    write(*,*)"WARNING: do NOT use it with Rashba coupling yet"
    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)
    delta=zero

    print*,"Get Gloc_iw:"
    allocate(gloc(Nspin,Nspin,Norb,Norb,NL))
    do i=1,NL
       iw = xi*wm(i)
       zeta = -ls2j_func(impSmats(:,:,:,:,i))
       forall(iorb=1:Nsporb)zeta(iorb,iorb)=zeta(iorb,iorb)+(iw+xmu)
       fg=zero
       do ik=1,Lk         
          gkinv = inverse_gk(zeta,Hk(:,:,ik))
          fg = fg + gkinv*dos_wt(ik)
       enddo
       gloc(:,:,:,:,i) = j2ls_func(fg)
       !
       !Get Delta=\Delta or G_0
       call matrix_inverse(fg)
       if(cg_scheme=='weiss')then
          self   = ls2j_func(impSmats(:,:,:,:,i))
          gdelta = fg + self
          call matrix_inverse(gdelta)
       else
          gdelta = zeta - bhzHloc - fg
       endif
       !
       delta(:,:,:,:,i) = j2ls_func(gdelta)
       !
    enddo
    ispin=1
    do iorb=1,Norb
       do jorb=iorb,Norb
          suffix="l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw.ed"
          call splot("Delta_"//reg(suffix),wm,delta(ispin,ispin,iorb,jorb,:))
          call splot("Gloc_"//reg(suffix),wm,gloc(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    deallocate(gloc)



    !REAL AXIS
    allocate(gloc(Nspin,Nspin,Norb,Norb,Nw))
    print*,"Get Gloc_realw:"
    do i=1,Nw
       iw=dcmplx(wr(i),eps)
       zeta = -ls2j_func(impSreal(:,:,:,:,i))
       forall(iorb=1:Nsporb)zeta(iorb,iorb)=zeta(iorb,iorb)+(iw+xmu)
       fg=zero
       do ik=1,Lk         
          gkinv = inverse_gk(zeta,Hk(:,:,ik))
          fg=fg+gkinv*dos_wt(ik)
       enddo
       !
       gloc(:,:,:,:,i) = j2ls_func(fg)
    enddo
    ispin=1
    do iorb=1,Norb
       do jorb=iorb,Norb
          suffix="l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_realw.ed"
          call splot("Gloc_"//reg(suffix),wr,gloc(ispin,ispin,iorb,jorb,:))
          call splot("DOS_"//reg(suffix),wr,-dimag(gloc(ispin,ispin,iorb,jorb,:))/pi)
       enddo
    enddo
    deallocate(gloc)


    write(*,*)"ntot=",sum(nimp)
    if(ntype==1)then
       nobj=nimp(1)
    else
       nobj=sum(nimp)
    endif

  end subroutine get_delta



  function inverse_gk(zeta,hk) result(gk)
    complex(8)                  :: zita(2)
    complex(8),dimension(4,4)   :: zeta,hk
    complex(8),dimension(4,4)   :: gk
    gk=zero
    zita(1)=zeta(1,1);zita(2)=zeta(2,2)
    gk(1:2,1:2) = inverse_gk2x2(zita,hk(1:2,1:2))
    zita(1)=zeta(3,3);zita(2)=zeta(4,4)
    gk(3:4,3:4) = inverse_gk2x2(zita,hk(3:4,3:4))
    ! gk=zeta-Hk                  !
    ! call matrix_inverse(gk)
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



  function ls2j(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop"error ls2j: iorb>Norb"
    if(ispin>Nspin)stop"error ls2j: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function ls2j


  function ls2j_func(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nsporb,Nsporb)         :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=ls2j(ispin,iorb)
                j=ls2j(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function ls2j_func

  function j2ls_func(fg) result(g)
    complex(8),dimension(Nsporb,Nsporb)         :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=ls2j(ispin,iorb)
                j=ls2j(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2ls_func




end program



