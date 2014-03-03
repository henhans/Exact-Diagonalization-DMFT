program ed_bhz_afm
  USE DMFT_ED
  USE PARSE_INPUT
  USE COMMON_VARS
  USE TOOLS
  USE ARRAYS
  USE ERROR
  USE MATRIX
  USE IOTOOLS
  implicit none
  integer                :: ip,iloop,Lk,Nsporb,afmNsporb,unit
  logical                :: converged
  integer                :: Nindep
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:,:)
  complex(8),allocatable :: Sreal(:,:,:,:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  character(len=16)      :: finput
  character(len=32)      :: hkfile

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")

  !
  call ed_read_input(trim(finput))
  !
  Nindep=4                      !number of independent sites, 4 for AFM ordering
  if(Nindep/=4.OR.Nspin/=2.OR.Norb/=2)stop "wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nsporb=Nspin*Norb
  afmNsporb=Nindep*Nsporb
  if(Nsporb/=4)stop "# of SObands should be 4"
  if(afmNsporb/=4)stop "# of AFM-SObands should be 16"

  !Allocate Weiss Field:
  allocate(delta(Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nindep,Nspin,Nspin,Norb,Norb,Lreal))

  !
  call read_hk(trim(hkfile))
  !

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nindep,Nb(1),Nb(2)))
  do ip=1,Nindep
     Hloc = bhzHloc_site(ip)
     write(*,*)"Updated Hloc:"
     call print_Hloc(Hloc)
     call init_ed_solver(bath(ip,:,:))
  enddo

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     do ip=1,Nindep
        call ed_solver(bath(ip,:,:))
        Smats(ip,:,:,:,:,:) = impSmats
        Sreal(ip,:,:,:,:,:) = impSreal
     enddo

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     do ip=1,Nindep
        call chi2_fitgf(delta(ip,1,1,:,:,:),bath(ip,:,:),ispin=1,iverbose=.true.)
        call chi2_fitgf(delta(ip,2,2,:,:,:),bath(ip,:,:),ispin=2,iverbose=.true.)
     enddo

     converged = check_convergence(sum(delta(:,1,1,1,1,:),1),dmft_error,nsuccess,nloop)

     call end_loop
  enddo



contains


  subroutine read_hk(file)
    character(len=*)       :: file
    integer                :: i,j,ik,iorb,jorb,isporb,jsporb,ispin,jspin
    real(8)                :: kx,ky,foo
    complex(8),allocatable :: fg(:,:,:)
    !
    open(50,file=file,status='old')
    read(50,*)Lk,foo,foo,foo,foo
    allocate(Hk(afmNsporb,afmNsporb,Lk),dos_wt(Lk),bhzHloc(afmNsporb,afmNsporb))
    do ik=1,Lk
       read(50,"(3(F10.7,1x))")kx,ky,foo
       do isporb=1,afmNsporb
          read(50,"(100(2F10.7,1x))")(Hk(isporb,jsporb,ik),jsporb=1,afmNsporb)
       enddo
    enddo
    dos_wt=1.d0/dble(Lk)
    bhzHloc = sum(Hk(:,:,:),dim=3)/dble(Lk)
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nsporb
    write(*,*)"# of AFM-SO-bands :",afmNsporb
  end subroutine read_hk




  subroutine get_delta
    integer                                       :: i,j,ik,iorb,jorb,ispin,jspin
    complex(8),dimension(Nindep,Nsporb,Nsporb)    :: zeta4x4,fg4x4
    complex(8),dimension(afmNsporb,afmNsporb)     :: zeta,gkinv,fg
    complex(8),dimension(Nsporb,Nsporb)           :: gdelta,self
    complex(8),dimension(:,:,:,:,:,:),allocatable :: gloc
    complex(8)                                    :: iw
    real(8)                                       :: wm(Lmats),wr(Lreal)
    character(len=20)                             :: suffix
    !
    write(*,*)"WARNING: do NOT use it with Rashba coupling yet"
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    delta=zero

    print*,"Get Gloc_iw:"
    allocate(gloc(Nindep,Nspin,Nspin,Norb,Norb,Lmats))
    do i=1,Lmats
       iw = xi*wm(i)
       do ip=1,Nindep
          forall(iorb=1:Nsporb)zeta4x4(ip,iorb,iorb)=iw+xmu
          zeta4x4(ip,:,:) = zeta4x4(ip,:,:)-ls2j(Smats(ip,:,:,:,:,i))
       enddo
       zeta = blocks4x4_to_matrix(zeta4x4)
       fg   = zero
       do ik=1,Lk         
          gkinv = inverse_gk(zeta,Hk(:,:,ik))
          fg = fg + gkinv*dos_wt(ik)
       enddo
       fg4x4 = matrix_to_blocks4x4(fg)
       do ip=1,Nindep
          gloc(ip,:,:,:,:,i) = j2ls(fg4x4(ip,:,:))
          !
          !Get Delta=\Delta or G_0
          call matrix_inverse(fg4x4(ip,:,:))
          if(cg_scheme=='weiss')then
             gdelta = fg4x4(ip,:,:) + ls2j(Smats(ip,:,:,:,:,i))
             call matrix_inverse(gdelta)
          else
             gdelta = zeta4x4(ip,:,:) - ls2j(bhzHloc_site(ip)) - fg4x4(ip,:,:)
          endif
          !
          delta(ip,:,:,:,:,i) = j2ls(gdelta)
          !
       enddo
    enddo
    do ispin=1,Nspin
       do iorb=1,Norb
          do jorb=iorb,Norb
             do ip=1,Nindep
                suffix="_site"//reg(txtfy(ip))//&
                     "_l"//reg(txtfy(iorb))//&
                     "_m"//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//"_iw.ed"
                call splot("Delta"//reg(suffix),wm,delta(ip,ispin,ispin,iorb,jorb,:))
                call splot("Gloc"//reg(suffix),wm,gloc(ip,ispin,ispin,iorb,jorb,:))
             enddo
          enddo
       enddo
    enddo
    deallocate(gloc)



    !REAL AXIS
    allocate(gloc(Nindep,Nspin,Nspin,Norb,Norb,Lreal))
    print*,"Get Gloc_realw:"
    do i=1,Lreal
       iw=dcmplx(wr(i),eps)
       do ip=1,Nindep
          forall(iorb=1:Nsporb)zeta4x4(ip,iorb,iorb)=iw+xmu
          zeta4x4(ip,:,:) = zeta4x4(ip,:,:)-ls2j(Sreal(ip,:,:,:,:,i))
       enddo

       zeta = blocks4x4_to_matrix(zeta4x4)
       fg   = zero
       do ik=1,Lk         
          gkinv = inverse_gk(zeta,Hk(:,:,ik))
          fg = fg + gkinv*dos_wt(ik)
       enddo
       fg4x4 = matrix_to_blocks4x4(fg)
       do ip=1,Nindep
          gloc(ip,:,:,:,:,i) = j2ls(fg4x4(ip,:,:))
       enddo
    enddo

    do ispin=1,Nspin
       do iorb=1,Norb
          do jorb=iorb,Norb
             do ip=1,Nindep
                suffix="_site"//reg(txtfy(ip))//&
                     "_l"//reg(txtfy(iorb))//&
                     "_m"//reg(txtfy(jorb))//&
                     "_s"//reg(txtfy(ispin))//"_realw.ed"
                call splot("Gloc_"//reg(suffix),wr,gloc(ip,ispin,ispin,iorb,jorb,:))
                call splot("DOS_"//reg(suffix),wr,-dimag(gloc(ip,ispin,ispin,iorb,jorb,:))/pi)
             enddo
          enddo
       enddo
    enddo
    deallocate(gloc)

  end subroutine get_delta




  function inverse_gk(zeta,hk) result(gk)
    complex(8),dimension(afmNsporb,afmNsporb)   :: zeta,hk
    complex(8),dimension(afmNsporb,afmNsporb)   :: gk
    gk=zeta-Hk
    call matrix_inverse(gk)
  end function inverse_gk





  !TRANSFORMATION BETWEEN DIFFERENT BASIS:
  function bhzHloc_site(isite) result(H)
    integer    :: isite,i,j
    complex(8) :: H(Nspin,Nspin,Norb,Norb)
    i=1+(isite-1)*Nsporb
    j=isite*Nsporb
    H = j2ls(bhzHloc(i:j,i:j))
  end function bhzHloc_site



  function blocks4x4_to_matrix(Vblocks) result(Matrix)
    complex(8),dimension(Nindep,Nsporb,Nsporb) :: Vblocks
    complex(8),dimension(afmNsporb,afmNsporb)  :: Matrix
    integer                                    :: i,j,ip
    Matrix=zero
    do ip=1,Nindep
       i = (ip-1)*Nsporb + 1
       j = ip*Nsporb
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks4x4_to_matrix


  function matrix_to_blocks4x4(Matrix) result(Vblocks)
    complex(8),dimension(Nindep,Nsporb,Nsporb) :: Vblocks
    complex(8),dimension(afmNsporb,afmNsporb)  :: Matrix
    integer                                    :: i,j,ip
    Vblocks=zero
    do ip=1,Nindep
       i = (ip-1)*Nsporb + 1
       j = ip*Nsporb
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks4x4


  function ls2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop"error ls2j_index: iorb>Norb"
    if(ispin>Nspin)stop"error ls2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function ls2j_index


  function ls2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nsporb,Nsporb)         :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=ls2j_index(ispin,iorb)
                j=ls2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function ls2j

  function j2ls(fg) result(g)
    complex(8),dimension(Nsporb,Nsporb)         :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=ls2j_index(ispin,iorb)
                j=ls2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2ls


end program ed_bhz_afm



