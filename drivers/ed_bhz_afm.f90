program ed_bhz_afm
  USE DMFT_ED
  USE SCIFOR
  implicit none
  integer                :: ip,iloop,Lk,Nso,afmNso
  logical                :: converged
  integer                :: Nindep
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:,:),Bath_(:,:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:,:,:)
  complex(8),allocatable :: Smats(:,:,:,:,:,:)
  complex(8),allocatable :: Sreal(:,:,:,:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:),bhzHloc(:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  integer                :: Nk
  real(8)                :: mh,lambda,wmixing
  character(len=16)      :: finput
  character(len=32)      :: hkfile
  character(len=32)      :: hamfile
  logical :: waverage
  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(waverage,"WAVERAGE",finput,default=.false.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  !
  call ed_read_input(trim(finput))
  !
  Nindep=4                      !number of independent sites, 4 for AFM ordering
  if(Nspin/=2.OR.Norb/=2)stop "wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb
  afmNso=Nindep*Nso!=16


  !Allocate Weiss Field:
  allocate(delta(Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nindep,Nspin,Nspin,Norb,Norb,Lreal))

  !
  call build_hk(trim(hkfile))
  hamfile=trim(Hfile)
  !

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nindep,Nb(1),Nb(2)))
  allocate(Bath_(Nindep,Nb(1),Nb(2)))
  do ip=1,Nindep
     ed_file_suffix="_is"//reg(txtfy(ip))
     call init_ed_solver(bath(ip,:,:))
     call break_symmetry_bath(bath(ip,:,:),sb_field,(-1.d0)**(ip+1))
     Hloc = bhzHloc_site(ip)
     write(LOGfile,*)"Updated Hloc, site:",ip
     call print_Hloc(Hloc)
  enddo

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     do ip=1,Nindep
        write(LOGfile,*)"Solving site:",ip
        ed_file_suffix="_is"//reg(txtfy(ip))
        Hloc = bhzHloc_site(ip)
        call ed_solver(bath(ip,:,:))
        Smats(ip,:,:,:,:,:) = impSmats
        Sreal(ip,:,:,:,:,:) = impSreal
     enddo

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Fit the new bath, starting from the old bath + the supplied delta
     do ip=1,Nindep
        ed_file_suffix="_is"//reg(txtfy(ip))
        Hloc = bhzHloc_site(ip)
        call chi2_fitgf(delta(ip,1,1,:,:,:),bath(ip,:,:),ispin=1)
        call chi2_fitgf(delta(ip,2,2,:,:,:),bath(ip,:,:),ispin=2)
     enddo


     !AVERAGE OUT 1-->4, 2-->3
     if(waverage)then
        bath(1,:,:)=(bath(1,:,:)+bath(4,:,:))/2.d0
        bath(2,:,:)=(bath(2,:,:)+bath(3,:,:))/2.d0
        bath(3,:,:)=bath(2,:,:)
        bath(4,:,:)=bath(1,:,:)
     endif
     !MIXING:
     do ip=1,Nindep
        if(iloop>1)bath(ip,:,:) = wmixing*bath(ip,:,:) + (1.d0-wmixing)*Bath_(ip,:,:)
        Bath_(ip,:,:)=bath(ip,:,:)
     enddo
     converged = check_convergence(delta(1,1,1,1,1,:)+delta(4,1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo



contains



  !---------------------------------------------------------------------
  !PURPOSE: GET DELTA FUNCTION
  !---------------------------------------------------------------------
  subroutine get_delta
    integer                                       :: i,j,ik,iorb,jorb,ispin,jspin,iso,jso
    complex(8),dimension(Nindep,Nso,Nso)          :: zeta_site,fg_site
    complex(8),dimension(afmNso,afmNso)           :: zeta,fg
    complex(8),dimension(Nso,Nso)                 :: gdelta,self
    complex(8),dimension(:,:,:,:,:,:),allocatable :: gloc
    complex(8)                                    :: iw
    complex(8),dimension(Lmats,afmNso,afmNso) :: gmats
    complex(8),dimension(Lreal,afmNso,afmNso) :: greal
    real(8)                                       :: wm(Lmats),wr(Lreal)
    character(len=32)                             :: suffix
    !
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    delta=zero

    print*,"Get Gloc_iw:"
    allocate(gloc(Nindep,Nspin,Nspin,Norb,Norb,Lmats))
    do i=1,Lmats
       iw = xi*wm(i)
       zeta_site=zero
       do ip=1,Nindep
          forall(iso=1:Nso)zeta_site(ip,iso,iso)=iw+xmu
          zeta_site(ip,:,:) = zeta_site(ip,:,:)-so2j(Smats(ip,:,:,:,:,i))
       enddo
       zeta = blocks_to_matrix(zeta_site)
       fg   = zero
       do ik=1,Lk
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       fg_site = matrix_to_blocks(fg)
       do ip=1,Nindep
          gloc(ip,:,:,:,:,i) = j2so(fg_site(ip,:,:))
          call matrix_inverse(fg_site(ip,:,:)) !Get G_loc^{-1} used later
          !
          if(cg_scheme=='weiss')then
             gdelta = fg_site(ip,:,:) + so2j(Smats(ip,:,:,:,:,i))
             call matrix_inverse(gdelta)
          else
             gdelta = zeta_site(ip,:,:) - so2j(bhzHloc_site(ip)) - fg_site(ip,:,:)
          endif
          !
          delta(ip,:,:,:,:,i) = j2so(gdelta)
          !
       enddo
    enddo
    !PRINT
    do ip=1,Nindep
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_is"//reg(txtfy(ip))//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call splot("Delta"//reg(suffix),wm,delta(ip,ispin,ispin,iorb,iorb,:))
             call splot("Gloc"//reg(suffix),wm,gloc(ip,ispin,ispin,iorb,iorb,:))
          enddo
       enddo
    enddo
    deallocate(gloc)



    !REAL AXIS
    allocate(gloc(Nindep,Nspin,Nspin,Norb,Norb,Lreal))
    print*,"Get Gloc_realw:"
    do i=1,Lreal
       iw=dcmplx(wr(i),eps)
       zeta_site=zero
       do ip=1,Nindep
          forall(iso=1:Nso)zeta_site(ip,iso,iso)=iw+xmu
          zeta_site(ip,:,:) = zeta_site(ip,:,:)-so2j(Sreal(ip,:,:,:,:,i))
       enddo
       zeta = blocks_to_matrix(zeta_site)
       fg   = zero
       do ik=1,Lk
          fg = fg + inverse_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       fg_site = matrix_to_blocks(fg)
       do ip=1,Nindep
          gloc(ip,:,:,:,:,i) = j2so(fg_site(ip,:,:))
       enddo
    enddo
    !PRINT
    do ip=1,Nindep
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_is"//reg(txtfy(ip))//"_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call splot("Gloc"//reg(suffix),wr,gloc(ip,ispin,ispin,iorb,iorb,:))
             call splot("DOS"//reg(suffix),wr,-dimag(gloc(ip,ispin,ispin,iorb,iorb,:))/pi)
          enddo
       enddo
    enddo
    deallocate(gloc)

  end subroutine get_delta








  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
  !--------------------------------------------------------------------!
  subroutine build_hk(file)
    character(len=*)                          :: file
    integer                                   :: i,j,ik=0
    integer                                   :: ix,iy
    real(8)                                   :: kx,ky    
    integer                                   :: iorb,jorb
    integer                                   :: isporb,jsporb
    integer                                   :: ispin,jspin
    real(8)                                   :: foo
    integer                                   :: unit
    complex(8),dimension(afmNso,afmNso)       :: Hkmix
    complex(8),dimension(Lmats,afmNso,afmNso) :: fg
    complex(8),dimension(Lreal,afmNso,afmNso) :: fgr
    real(8)                                   :: wm(Lmats),wr(Lreal),dw,n0(afmNso)
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal,mesh=dw)
    Lk=Nk**2
    write(LOGfile,*)"Build H(k) AFM-BHZ:"
    allocate(Hk(afmNso,afmNso,Lk))
    unit=free_unit()
    open(unit,file=file)
    call start_progress(LOGfile)
    fg=zero
    fgr=zero
    do ix=1,Nk
       kx = -pi + 2.d0*pi*dble(ix-1)/dble(Nk)
       do iy=1,Nk
          ky = -pi + 2.d0*pi*dble(iy-1)/dble(Nk)
          ik=ik+1
          Hkmix         = hk_bhz_afm(kx,ky)
          Hk(:,:,ik)    = iso2site(Hkmix)
          write(unit,"(3(F10.7,1x))")kx,ky,pi
          do i=1,afmNso
             write(unit,"(100(2F10.7,1x))")(Hk(i,j,ik),j=1,Norb)
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
    do i=1,afmNso
       n0(i) = -2.d0*sum(dimag(fgr(:,i,i))*fermi(wr(:),beta))*dw/pi
    enddo
    write(unit,"(24F20.12)")mh,lambda,xmu,(n0(i),i=1,afmNso),sum(n0)
    write(LOGfile,"(24F20.12)")mh,lambda,xmu,(n0(i),i=1,afmNso),sum(n0)
    open(10,file="U0_DOS.ed")
    do i=1,Lreal
       write(10,"(100(F25.12))") wr(i),(-dimag(fgr(i,iorb,iorb))/pi,iorb=1,afmNso)
    enddo
    close(10)
    open(11,file="U0_Gloc_iw.ed")
    do i=1,Lmats
       write(11,"(20(2F20.12))") wm(i),(fg(i,iorb,iorb),iorb=1,afmNso)
    enddo
    close(11)
    allocate(bhzHloc(afmNso,afmNso))
    bhzHloc = sum(Hk(:,:,:),dim=3)/dble(Lk)
    where(abs(dreal(bhzHloc))<1.d-9)bhzHloc=0.d0
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nso
    write(*,*)"# of AFM-SO-bands :",afmNso
  end subroutine build_hk



  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function inverse_gk(zeta,hk) result(gk)
    complex(8),dimension(afmNso,afmNso)   :: zeta,hk
    complex(8),dimension(afmNso,afmNso)   :: gk
    integer :: i
    ! gk=zeta-Hk
    gk = -Hk
    forall(i=1:afmNso)gk(i,i)=zeta(i,i) - hk(i,i)
    call matrix_inverse(gk)
  end function inverse_gk

  function inverse_g0k(zeta,hk,type) result(g0k)
    integer :: i
    complex(8)                  :: zeta
    complex(8),dimension(16,16) :: hk
    complex(8),dimension(16,16) :: g0k,g0ktmp
    integer,optional            :: type
    integer                     :: type_
    type_=0 ; if(present(type))type_=type
    if(type_==0)then
       g0k = -hk
       forall(i=1:16)g0k(i,i)=zeta - hk(i,i)
       call matrix_inverse(g0k)
    else
       g0k=zero
       g0k(1:8,1:8)   = inverse_g0k_8x8(zeta,hk(1:8,1:8))
       g0k(9:16,9:16) = inverse_g0k_8x8(zeta,hk(9:16,9:16))
    endif
  end function inverse_g0k

  function inverse_g0k_8x8(zeta,hk) result(g0k)
    complex(8)                :: zeta
    complex(8),dimension(8,8) :: hk
    complex(8),dimension(8,8) :: g0k
    integer :: i
    g0k = -hk
    forall(i=1:8)g0k(i,i)=zeta - hk(i,i)
    call matrix_inverse(g0k)
  end function inverse_g0k_8x8



  !--------------------------------------------------------------------!
  !BHZ - AFM HAMILTONIAN:
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function hk_bhz_afm(kx,ky) result(hk)
    real(8)                     :: kx,ky
    complex(8),dimension(16,16) :: hk
    hk            = zero
    hk(1:8,1:8)   = hk_bhz_afm8x8(kx,ky)
    hk(9:16,9:16) = conjg(hk_bhz_afm8x8(-kx,-ky))
  end function hk_bhz_afm

  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function hk_bhz_afm8x8(kx,ky) result(hk)
    real(8)                     :: kx,ky
    complex(8),dimension(8,8)   :: hk
    hk(1,1:8)=[one*mh,zero,-one/2.d0,xi*lambda/2.d0,-one/2.d0,-one*lambda/2.d0,zero,zero]
    hk(2,1:8)=[zero,-one*mh,xi*lambda/2.d0,one/2.d0,one*lambda/2.d0,one/2.d0,zero,zero]
    hk(3,1:8)=[-one/2.d0,-xi*lambda/2.d0,one*mh,zero,zero,zero,-one/2.d0,-one*lambda/2.d0]
    hk(4,1:8)=[-xi*lambda/2.d0,one/2.d0,zero,-one*mh,zero,zero,one*lambda/2.d0,one/2.d0]
    hk(5,1:8)=[-one/2.d0,one*lambda/2.d0,zero,zero,one*mh,zero,-one/2.d0,xi*lambda/2.d0]
    hk(6,1:8)=[-one*lambda/2.d0,one/2.d0,zero,zero,zero,-one*mh,xi*lambda/2.d0,one/2.d0]
    hk(7,1:8)=[zero,zero,-one/2.d0,one*lambda/2.d0,-one/2.d0,-xi*lambda/2.d0,one*mh,zero]
    hk(8,1:8)=[zero,zero,-one*lambda/2.d0,one/2.d0,-xi*lambda/2.d0,one/2.d0,zero,-one*mh]
    !
    hk(1,1:8)=hk(1,1:8)+[zero,zero,exp(-xi*2.d0*kx)*(-0.5d0),exp(-xi*2.d0*kx)*(-xi*lambda/2.d0),exp(xi*2.d0*ky)*(-0.5d0),exp(xi*2.d0*ky)*(lambda/2.d0),zero,zero]
    hk(2,1:8)=hk(2,1:8)+[zero,zero,exp(-xi*2.d0*kx)*(-xi*lambda/2.d0),exp(-xi*2.d0*kx)*(0.5d0),exp(xi*2.d0*ky)*(-lambda/2.d0),exp(xi*2.d0*ky)*(0.5d0),zero,zero]
    hk(3,1:8)=hk(3,1:8)+[exp(xi*2.d0*kx)*(-0.5d0),exp(xi*2.d0*kx)*(xi*lambda/2.d0),zero,zero,zero,zero,exp(xi*2.d0*ky)*(-0.5d0),exp(xi*2.d0*ky)*(lambda/2.d0)]
    hk(4,1:8)=hk(4,1:8)+[exp(xi*2.d0*kx)*(xi*lambda/2.d0),exp(xi*2.d0*kx)*(0.5d0),zero,zero,zero,zero,exp(xi*2.d0*ky)*(-lambda/2.d0),exp(xi*2.d0*ky)*(0.5d0)]
    hk(5,1:8)=hk(5,1:8)+[exp(-xi*2.d0*ky)*(-0.5d0),exp(-xi*2.d0*ky)*(-lambda/2.d0),zero,zero,zero,zero,exp(-xi*2.d0*kx)*(-0.5d0),exp(-xi*2.d0*kx)*(-xi*lambda/2.d0)]
    hk(6,1:8)=hk(6,1:8)+[exp(-xi*2.d0*ky)*(lambda/2.d0),exp(-xi*2.d0*ky)*(0.5d0),zero,zero,zero,zero,exp(-xi*2.d0*kx)*(-xi*lambda/2.d0),exp(-xi*2.d0*kx)*(0.5d0)]
    hk(7,1:8)=hk(7,1:8)+[zero,zero,exp(-xi*2.d0*ky)*(-0.5d0),exp(-xi*2.d0*ky)*(-lambda/2.d0),exp(xi*2.d0*kx)*(-0.5d0),exp(xi*2.d0*kx)*(xi*lambda/2.d0),zero,zero]
    hk(8,1:8)=hk(8,1:8)+[zero,zero,exp(-xi*2.d0*ky)*(lambda/2.d0),exp(-xi*2.d0*ky)*(0.5d0),exp(xi*2.d0*kx)*(xi*lambda/2.d0),exp(xi*2.d0*kx)*(0.5d0),zero,zero]
  end function hk_bhz_afm8x8









  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS:
  !--------------------------------------------------------------------!
  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function bhzHloc_site(isite) result(H)
    integer                              :: isite
    complex(8),dimension(Nindep,Nso,Nso) :: Vblocks
    complex(8)                           :: H(Nspin,Nspin,Norb,Norb)
    Vblocks = matrix_to_blocks(bhzHloc)
    H = j2so(Vblocks(isite,:,:))
  end function bhzHloc_site


  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function blocks_to_matrix(Vblocks) result(Matrix)
    complex(8),dimension(Nindep,Nso,Nso) :: Vblocks
    complex(8),dimension(afmNso,afmNso)  :: Matrix
    integer                                    :: i,j,ip
    Matrix=zero
    do ip=1,Nindep
       i = (ip-1)*Nso + 1
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix


  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function matrix_to_blocks(Matrix) result(Vblocks)
    complex(8),dimension(Nindep,Nso,Nso) :: Vblocks
    complex(8),dimension(afmNso,afmNso)  :: Matrix
    integer                                    :: i,j,ip
    Vblocks=zero
    do ip=1,Nindep
       i = (ip-1)*Nso + 1
       j = ip*Nso
       Vblocks(ip,:,:) = Matrix(i:j,i:j)
    enddo
  end function matrix_to_blocks



  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function matrix_to_block(ip,Matrix) result(Vblock)
    complex(8),dimension(Nso,Nso) :: Vblock
    complex(8),dimension(afmNso,afmNso)  :: Matrix
    integer                              :: i,j,ip
    Vblock=zero
    i = (ip-1)*Nso + 1
    j = ip*Nso
    Vblock(:,:) = Matrix(i:j,i:j)
  end function matrix_to_block




  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop"error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop"error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index

  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
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


  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
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


  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function iso2site(hktmp) result(hk)
    complex(8),dimension(16,16) :: hktmp,hk
    integer                     :: icount,isite,ispin,iband,iorb
    integer                     :: jcount,jsite,jspin,jband,jorb
    integer                     :: row,col
    icount=0
    do isite=1,4
       do ispin=1,2
          do iorb=1,2
             icount=icount+1
             row = (ispin-1)*2*4 + (isite-1)*2 + iorb
             jcount=0
             do jsite=1,4
                do jspin=1,2
                   do jorb=1,2
                      jcount=jcount+1
                      col = (jspin-1)*2*4 + (jsite-1)*2 + jorb
                      hk(icount,jcount) = hktmp(row,col)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function iso2site


  !--------------------------------------------------------------------!
  !PURPOSE: 
  !--------------------------------------------------------------------!
  function site2iso(hk) result(hktmp)
    complex(8),dimension(16,16) :: hktmp,hk
    integer                     :: icount,isite,ispin,iband,iorb
    integer                     :: jcount,jsite,jspin,jband,jorb
    integer                     :: row,col
    icount=0
    do isite=1,4
       do ispin=1,2
          do iorb=1,2
             icount=icount+1
             row = (ispin-1)*2*4 + (isite-1)*2 + iorb
             jcount=0
             do jsite=1,4
                do jspin=1,2
                   do jorb=1,2
                      jcount=jcount+1
                      col = (jspin-1)*2*4 + (jsite-1)*2 + jorb
                      hktmp(row,col) = hk(icount,jcount)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function site2iso

end program ed_bhz_afm



