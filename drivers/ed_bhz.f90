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
  USE SCIFOR
  implicit none
  integer                :: iloop,Lk,Nso
  logical                :: converged
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:),Bath_(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  integer                :: Nk
  real(8)                :: mh,lambda,wmixing,akrange
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  logical                :: spinsym,getak

  call MPI_INIT(ED_MPI_ERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ED_MPI_ID,ED_MPI_ERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ED_MPI_SIZE,ED_MPI_ERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',ED_MPI_ID,' of ',ED_MPI_SIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,ED_MPI_ERR)

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(getak,"GETAK",finput,default=.false.)
  call parse_input_variable(akrange,"AKRANGE",finput,default=3.d0)
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"spinsym",finput,default=.true.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  !
  call ed_read_input(trim(finput))

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(delta(Nspin,Nspin,Norb,Norb,Lmats))


  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))

  !Setup solver
  Nb=get_bath_size()
  allocate(Bath(Nb(1),Nb(2)))
  allocate(Bath_(Nb(1),Nb(2)))
  call init_ed_solver(bath)
  Hloc = j2so(bhzHloc)
  write(LOGfile,*)"Updated Hloc:"
  call print_Hloc(Hloc)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(ED_MPI_ID==0)call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta
     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(1,1,:,:,:),bath,ispin=1)
     if(.not.spinsym)then
        call chi2_fitgf(delta(2,2,:,:,:),bath(:,:),ispin=2)
     else
        call spin_symmetrize_bath(bath(:,:))
     endif

     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     Bath_=Bath
     if(ED_MPI_ID==0)converged = check_convergence(delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ED_MPI_ERR)
     if(ED_MPI_ID==0)call end_loop
  enddo


  call MPI_FINALIZE(ED_MPI_ERR)

contains




  !---------------------------------------------------------------------
  !PURPOSE: GET DELTA FUNCTION
  !---------------------------------------------------------------------
  subroutine get_delta
    integer                                     :: i,j,ik,iorb,jorb,ispin,jspin,iso,unit
    complex(8),dimension(Nso,Nso)               :: zeta,fg,gdelta,fgk
    complex(8),dimension(:,:,:,:,:),allocatable :: gloc,Sreal
    complex(8),dimension(:,:,:,:),allocatable :: gk
    complex(8),dimension(:,:,:,:),allocatable   :: gfoo
    complex(8)                                  :: iw
    real(8)                                     :: wm(Lmats),wr(Lreal),reS(Nspin),imS(Nspin),ww
    character(len=20)                           :: suffix
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    delta=zero

    if(getak)then
       print*,"Get A(k,w):"
       allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
       do iorb=1,Norb
          unit=free_unit()
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_realw.ed"
          open(unit,file="impSigma"//reg(suffix),status='old')
          do i=1,Lreal
             read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
             Sreal(ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
          enddo
          close(unit)
       enddo

       allocate(gk(Lk,Nspin,Norb,Lreal))
       allocate(gfoo(Nspin,Nspin,Norb,Norb))
       call start_progress(LOGfile)
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          zeta=zero
          forall(iso=1:Nso)zeta(iso,iso)=iw+xmu
          zeta(:,:) = zeta(:,:)-so2j(Sreal(:,:,:,:,i))
          do ik=1,Lk
             fgk = inverse_gk(zeta,Hk(:,:,ik))
             gfoo(:,:,:,:) = j2so(fgk(:,:))
             do ispin=1,Nspin
                do iorb=1,Norb
                   gk(ik,ispin,iorb,i) = gfoo(ispin,ispin,iorb,iorb)
                enddo
             enddo
          enddo
          call progress(i,Lreal)
       enddo
       call stop_progress()

       !PRINT
       do ispin=1,Nspin
          do iorb=1,Norb
             unit=free_unit()
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call splot3d("Ak"//reg(suffix),(/(dble(ik),ik=1,Lk)/),wr,-dimag(gk(:,ispin,iorb,:))/pi,ymin=-akrange,ymax=akrange,nosurface=.true.)
             print*,"printing: Ak"//reg(suffix)
          enddo
       enddo
       deallocate(gk,gfoo)
       return
    endif





    !MATSUBARA AXIS
    print*,"Get Gloc_iw:"
    allocate(gloc(Nspin,Nspin,Norb,Norb,Lmats))
    do i=1,Lmats
       iw = xi*wm(i)
       forall(iorb=1:Nso)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - so2j(impSmats(:,:,:,:,i))
       fg=zero
       do ik=1,Lk
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,:,:,i) = j2so(fg)
       !
       !Get Delta=\Delta or G_0
       call matrix_inverse(fg)
       if(cg_scheme=='weiss')then
          gdelta = fg + so2j(impSmats(:,:,:,:,i))
          call matrix_inverse(gdelta)
       else
          gdelta = zeta(:,:) - bhzHloc - fg(:,:)
       endif
       !
       delta(:,:,:,:,i) = j2so(gdelta(:,:))
       !
    enddo
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
          call splot("Delta"//reg(suffix),wm,delta(ispin,ispin,iorb,iorb,:))
          call splot("Gloc"//reg(suffix),wm,gloc(ispin,ispin,iorb,iorb,:))
       enddo
    enddo
    deallocate(gloc)






    !REAL AXIS
    allocate(gloc(Nspin,Nspin,Norb,Norb,Lreal))
    print*,"Get Gloc_realw:"
    do i=1,Lreal
       iw=dcmplx(wr(i),eps)
       forall(iorb=1:Nso)zeta(iorb,iorb)=iw+xmu
       zeta(:,:) = zeta(:,:) - so2j(impSreal(:,:,:,:,i))
       fg=zero
       do ik=1,Lk         
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,:,:,i) = j2so(fg)
    enddo
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call splot("Gloc"//reg(suffix),wr,-dimag(gloc(ispin,ispin,iorb,iorb,:))/pi,dreal(gloc(ispin,ispin,iorb,iorb,:)))
       enddo
    enddo
    deallocate(gloc)

  end subroutine get_delta















  subroutine build_hk(file)
    character(len=*)                    :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Lmats,Nso,Nso) :: fg
    complex(8),dimension(Lreal,Nso,Nso) :: fgr
    real(8)                             :: wm(Lmats),wr(Lreal),dw,n0(Nso),eig(Nso)
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    call start_progress(LOGfile)
    if(getak)then
       write(LOGfile,*)"Build H(k) BHZ on the path GXMG:"
       unit=free_unit() 
       open(unit,file="Eigenbands.dat")
       Lk=3*Nk
       ik = 0
       allocate(Hk(Nso,Nso,Lk))
       !From \Gamma=(0,0) to X=(pi,0): Nk steps
       do ix=1,Nk
          ik=ik+1
          kx = 0.d0 + pi*real(ix-1,8)/dble(Nk)
          ky = 0.d0
          Hk(:,:,ik) = hk_bhz(kx,ky)
          eig        = Eigk(Hk(:,:,ik))
          call progress(ik,Lk)
          write(unit,"(I,16F25.12)")ik,(eig(i),i=1,Nso)
       enddo
       !From X=(pi,0) to M=(pi,pi): Nk steps
       do iy=1,Nk
          ik=ik+1
          kx = pi
          ky = 0.d0 + pi*real(iy-1,8)/dble(Nk)
          Hk(:,:,ik) = hk_bhz(kx,ky)
          eig        = Eigk(Hk(:,:,ik))
          call progress(ik,Lk)
          write(unit,"(I,16F25.12)")ik,(eig(i),i=1,Nso)
       enddo
       !From M=(pi,pi) to \Gamma=(0,0): Nk steps
       do ix=1,Nk
          ik=ik+1
          iy=ix
          kx = pi - pi*real(ix-1,8)/dble(Nk)
          ky = pi - pi*real(iy-1,8)/dble(Nk)
          Hk(:,:,ik) = hk_bhz(kx,ky)
          eig        = Eigk(Hk(:,:,ik))
          call progress(ik,Lk)
          write(unit,"(I,16F25.12)")ik,(eig(i),i=1,Nso)
       enddo
       close(unit)
       !
    else
       !
       write(LOGfile,*)"Build H(k) for BHZ:"
       Lk=Nk**2
       allocate(Hk(Nso,Nso,Lk))
       unit=free_unit()
       open(unit,file=file)
       fg=zero
       fgr=zero
       do ix=1,Nk
          kx = -pi + 2.d0*pi*dble(ix-1)/dble(Nk)
          do iy=1,Nk
             ky = -pi + 2.d0*pi*dble(iy-1)/dble(Nk)
             ik=ik+1
             Hk(:,:,ik) = hk_bhz(kx,ky)
             write(unit,"(3(F10.7,1x))")kx,ky,pi
             do i=1,Nso
                write(unit,"(100(2F10.7,1x))")(Hk(i,j,ik),j=1,Nso)
             enddo
             do i=1,Lreal
                fgr(i,:,:)=fgr(i,:,:) + inverse_g0k(dcmplx(wr(i),eps)+xmu,Hk(:,:,ik))
             enddo
             do i=1,Lmats
                fg(i,:,:) =fg(i,:,:)  + inverse_g0k(xi*wm(i)+xmu,Hk(:,:,ik))
             enddo
             call progress(ik,Lk)
          enddo
       enddo
       call stop_progress()
       write(unit,*)""
       allocate(dos_wt(Lk))
       dos_wt=1.d0/dble(Lk)
       fgr= fgr/dble(Lk)
       fg = fg/dble(Lk)
       do i=1,Nso
          n0(i) = -2.d0*sum(dimag(fgr(:,i,i))*fermi(wr(:),beta))*dw/pi
       enddo
       write(unit,"(24F20.12)")mh,lambda,xmu,(n0(i),i=1,Nso),sum(n0)
       write(LOGfile,"(24F20.12)")mh,lambda,xmu,(n0(i),i=1,Nso),sum(n0)
       open(10,file="U0_DOS.ed")
       do i=1,Lreal
          write(10,"(100(F25.12))") wr(i),(-dimag(fgr(i,iorb,iorb))/pi,iorb=1,Nso)
       enddo
       close(10)
       open(11,file="U0_Gloc_iw.ed")
       do i=1,Lmats
          write(11,"(20(2F20.12))") wm(i),(fg(i,iorb,iorb),iorb=1,Nso)
       enddo
       close(11)
       allocate(bhzHloc(Nso,Nso))
       bhzHloc = sum(Hk(:,:,:),dim=3)/dble(Lk)
       where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    endif
    !
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nso
  end subroutine build_hk







  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kx,ky) result(hk)
    real(8)                   :: kx,ky
    complex(8),dimension(4,4) :: hk
    Hk          = zero
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
    ! Hk(1,4) = -delta ; Hk(4,1)=-delta
    ! Hk(2,3) =  delta ; Hk(3,2)= delta
    ! Hk(1,3) = xi*rh*(sin(kx)-xi*sin(ky))
    ! Hk(3,1) =-xi*rh*(sin(kx)+xi*sin(ky))
  end function hk_bhz

  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky,epsik
    complex(8),dimension(2,2) :: hk
    epsik   = cos(kx)+cos(ky)
    hk(1,1) = mh - epsik
    hk(2,2) =-mh + epsik
    hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
    hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  end function hk_bhz2x2

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


  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k=zero
    g0k(1:2,1:2) = inverse_g0k2x2(iw,hk(1:2,1:2))
    g0k(3:4,3:4) = inverse_g0k2x2(iw,hk(3:4,3:4))
    ! else
    !    g0k = -hk
    !    forall(i=1:4)g0k(i,i) = iw + xmu + g0k(i,i)
    !    call matrix_inverse(g0k)
    ! endif
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


  function Eigk(hk) result(eig)
    complex(8),dimension(4,4) :: hk
    real(8),dimension(4)      :: eig
    call matrix_diagonalize(hk,eig)
  end function Eigk

  function Eigk2x2(hk) result(eig)
    complex(8),dimension(2,2) :: hk
    real(8),dimension(2)      :: eig
    call matrix_diagonalize(hk,eig)
    eig(1)=hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig(2)=hk(1,1)+hk(2,2) - sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig = eig/2.d0
  end function Eigk2x2





  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS:
  !--------------------------------------------------------------------!
  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop"error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop"error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nso,Nso)         :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nso,Nso)         :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so




end program ed_bhz



