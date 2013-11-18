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
  integer                :: iloop,Lk,Lorb
  logical                :: converged
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:),Fold(:,:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  real(8),allocatable    :: e0(:)
  !variables for the model:
  character(len=16)      :: finput
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
  call parse_cmd_variable(fnomag,"FNOMAG",default=.true.)
  call parse_cmd_variable(wmixing,"WMIXING",default=0.5d0)
  !
  call ed_read_input(trim(finput))
  !
  call read_hk(trim(hkfile))


  !Allocate Weiss Field:
  allocate(delta(Nspin,Norb,Norb,NL))
  allocate(Fold(Nspin,Norb,Norb,NL))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)


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
     call chi2_fitgf(delta(1,:,:,:),bath,ispin=1)
     if(fnomag)then
        call spin_symmetrize_bath(bath)
     else
        call chi2_fitgf(delta(2,:,:,:),bath,ispin=2)
     endif

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,1,:)+delta(1,2,2,:)+&
          delta(1,3,3,:)+delta(1,4,4,:),dmft_error,nsuccess,nloop)
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
    integer                :: i,j,ik,iorb,jorb,Norb_d
    real(8)                :: de,e,kx,ky,kz,foo
    real(8)                :: ntotal,npimp
    complex(8),allocatable :: fg(:,:,:)
    complex(8)             :: iw
    real(8)                :: wr(Nw)
    !
    open(50,file=file,status='old')
    read(50,*)Lk,Norb_d,foo,foo,foo
    if(Nspin/=2.OR.Norb/=2)stop "wrong setup from input file: Nspin=Norb=2"
    if(Norb_d/=Norb)stop "Norb_d in Hk should be 2."
    Lorb=Nspin*Norb
    if(Lorb/=4)stop "# of SObands should be 4"
    allocate(Hk(Lorb,Lorb,Lk),dos_wt(Lk),e0(Lorb),bhzHloc(Lorb,Lorb))
    do ik=1,Lk
       read(50,"(3(F10.7,1x))")kx,ky,foo
       do iorb=1,Lorb
          read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Lorb)
       enddo
    enddo
    dos_wt=1.d0/dble(Lk)
    bhzHloc(:,:) = sum(Hk(:,:,:),dim=3)/dble(Lk)
    where(dreal(bhzHloc) < 1.d-15)bhzHloc=0.d0
    write(*,*)"# of k-points:",Lk
    write(*,*)"# of SO-bands :",Lorb
    write(*,*)"  Hloc:",bhzHloc
    print*,""
    ! write(*,*)"Evaluate non-interacting G(w):"
    ! wr = linspace(wini,wfin,Nw)
    ! allocate(fg(Nw,Lorb,Lorb))
    ! fg=zero
    ! e0=0.d0
    ! do ik=1,Lk
    !    do i=1,Nw
    !       iw = cmplx(wr(i),eps,8)+xmu
    !       fg(i,:,:)=fg(i,:,:)+inverse_g0k(iw,Hk(:,:,ik))*dos_wt(ik)
    !    enddo
    !    do i=1,Lorb
    !       e0(i)=e0(i)+Hk(i,i,ik)*dos_wt(ik)
    !    enddo
    ! enddo
    ! write(*,*)" Centers of mass:",e0
    ! call splot("DOS0.ed",wr,-dimag(fg(:,1,1))/pi,-dimag(fg(:,2,2))/pi,&
    !      -dimag(fg(:,3,3))/pi,-dimag(fg(:,4,4))/pi)
    ! deallocate(fg)
  end subroutine read_hk



  subroutine get_delta
    integer                                   :: i,j,ik,iorb,jorb,ispin
    integer                                   :: so2j(Nspin,Norb)
    complex(8)                                :: iw
    complex(8),dimension(Lorb,Lorb)           :: zeta,fg,gdelta
    complex(8),dimension(:,:,:,:),allocatable :: gloc
    real(8)                                   :: wm(NL),wr(Nw)
    !
    do ispin=1,Nspin
       do iorb=1,Norb
          so2j(ispin,iorb)=(ispin-1)*Nspin + iorb
       enddo
    enddo
    !
    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)
    !
    delta=zero
    !
    print*,"Get Gloc_iw:"
    allocate(gloc(Nspin,Norb,Norb,NL));gloc=zero
    do i=1,NL
       zeta=zero
       iw = xi*wm(i)
       do ispin=1,Nspin
          do iorb=1,Norb
             zeta(so2j(ispin,iorb),so2j(ispin,iorb))= (iw + xmu)
             do jorb=1,Norb
                zeta(so2j(ispin,iorb),so2j(ispin,jorb)) = zeta(so2j(ispin,iorb),so2j(ispin,jorb)) -&
                     impSmats(ispin,iorb,jorb,i)
             enddo
          enddo
       enddo
       !
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo

       do ispin=1,Nspin
          do iorb=1,Norb
             gloc(ispin,iorb,iorb,i)  = fg(so2j(ispin,iorb),so2j(ispin,iorb))
          enddo
       enddo
       !
       call matrix_inverse(fg)
       !
       gdelta = zeta - bhzHloc - fg
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                delta(ispin,iorb,jorb,i)  = gdelta(so2j(ispin,iorb),so2j(ispin,iorb))
             enddo
          enddo
       enddo
       !
    enddo
    do ispin=1,Norb
       do iorb=1,Norb
          call splot("Delta_so"//reg(txtfy(so2j(ispin,iorb)))//"_iw.ed",wm,delta(ispin,iorb,iorb,:))
          call splot("Gloc_so"//reg(txtfy(so2j(ispin,iorb)))//"_iw.ed",wm,gloc(ispin,iorb,iorb,:))
       enddo
    enddo
    deallocate(gloc)



    !REAL AXIS
    allocate(gloc(Nspin,Norb,Norb,Nw))
    print*,"Get Gloc_realw:"
    do i=1,Nw
       iw=cmplx(wr(i),eps)
       do ispin=1,Nspin
          do iorb=1,Norb
             zeta(so2j(ispin,iorb),so2j(ispin,iorb))= (iw + xmu)
             do jorb=1,Norb
                zeta(so2j(ispin,iorb),so2j(ispin,jorb)) = zeta(so2j(ispin,iorb),so2j(ispin,jorb)) -&
                     impSreal(ispin,iorb,jorb,i)
             enddo
          enddo
       enddo
       !
       fg=zero
       do ik=1,Lk         
          fg=fg+inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             gloc(ispin,iorb,iorb,i)  = fg(so2j(ispin,iorb),so2j(ispin,iorb))
          enddo
       enddo
    enddo
    do ispin=1,Norb
       do iorb=1,Norb
          call splot("Gloc_so"//reg(txtfy(so2j(ispin,iorb)))//"_realw.ed",wr,gloc(ispin,iorb,iorb,:))
          call splot("DOS_so"//reg(txtfy(so2j(ispin,iorb)))//".ed",wr,-dimag(gloc(ispin,iorb,iorb,:))/pi)
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



  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k=zero
    g0k(1:2,1:2) = inverse_g0k2x2(iw,hk(1:2,1:2))
    g0k(3:4,3:4) = inverse_g0k2x2(iw,hk(3:4,3:4))
    !g0k=-Hk
    !forall(i=1:4)g0k(i,i)=iw-Hk
    !call matrix_inverse(g0k)
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
    complex(8)                  :: zita(2)
    complex(8),dimension(4,4)   :: zeta,hk
    complex(8),dimension(4,4)   :: gk
    gk=zero
    zita(1)=zeta(1,1);zita(2)=zeta(2,2)
    gk(1:2,1:2) = inverse_gk2x2(zita,hk(1:2,1:2))
    zita(1)=zeta(3,3);zita(2)=zeta(4,4)
    gk(3:4,3:4) = inverse_gk2x2(zita,hk(3:4,3:4))
    !
    ! gk=zeta-Hk                  !
    ! call matrix_inverse(gk)
    !
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

end program ed_bhz



