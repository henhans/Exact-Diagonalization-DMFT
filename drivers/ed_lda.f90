!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE COMMON_VARS
  USE IOTOOLS
  USE FUNCTIONS
  USE TOOLS
  USE MATRIX
  USE ERROR
  USE ARRAYS
  USE PARSE_INPUT
  implicit none
  integer                :: iloop,Lk
  logical                :: converged
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: fg0(:,:,:)
  real(8),allocatable    :: dos_wt(:)
  !variables for the model:
  character(len=32)      :: hkfile
  integer                :: ntype
  real(8)                :: nobj
  logical                :: fbethe
  real(8)                :: alpha(2),ts,thop(3)
  character(len=16)      :: finput,fhloc

  !parse additional variables && read input file && read Hk
  call parse_cmd_variable(finput,'FINPUT',default='inputED.in')
  !call parse_input_variable(hkfile,"HKFILE",finput,default='hkfile.in')
  call parse_cmd_variable(ntype,"NTYPE",default=0)
  call parse_cmd_variable(fbethe,"FBETHE",default=.true.)
  call parse_cmd_variable(Lk,"Lk",default=1000)
  call parse_cmd_variable(ts,"TS",default=0.5d0)
  call parse_cmd_variable(alpha,"alpha",default=[1.d0,1.d0])
  !
  call ed_read_input(trim(finput))
  !



  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,Lmats))


  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)

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

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(nobj,converged)
     call end_loop
  enddo


contains


  subroutine read_hk(file)
    character(len=*) :: file
    integer          :: i,j,ik,iorb,jorb,Norb_d,Norb_p
    real(8)          :: de,e,kx,ky,kz,foo

    if(fbethe)then
       !
       allocate(Hk(Norb,Norb,Lk))
       allocate(dos_wt(Lk))
       thop(1)=ts
       do iorb=2,Norb
          thop(iorb)=alpha(iorb-1)*ts
       enddo
       de=2.d0/dble(Lk)
       do ik=1,Lk
          e = -1.d0 + dble(ik-1)*de
          Hk(:,:,ik) = Hloc(1,1,:,:) !set local part
          do iorb=1,Norb
             Hk(iorb,iorb,ik) = Hk(iorb,iorb,ik) - 2.d0*thop(iorb)*e
          enddo
          dos_wt(ik)=dens_bethe(e,2.d0*ts)*de
       enddo
       !
    else
       !
       open(50,file=file,status='old')
       read(50,*)Lk,Norb_d,Norb_p,foo,foo
       if((Norb_d+Norb_p)/=Norb)stop "Can not read Hk.file: check Norb_d,Norb_p or FBETHE flag"
       allocate(Hk(Norb,Norb,Lk))
       allocate(dos_wt(Lk))
       do ik=1,Lk
          read(50,"(3(F10.7,1x))")kx,ky,kz
          do iorb=1,Norb
             read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Norb)
          enddo
       enddo
       close(50)
       dos_wt=1.d0/dble(Lk)
       Hloc(1,1,:Norb,:Norb) = sum(Hk(:,:,:),dim=3)/dble(Lk)
       !
    endif
    print*,"Hloc:"
    call print_Hloc(Hloc)
  end subroutine read_hk


  !+----------------------------------------+


  subroutine get_delta
    integer                                 :: i,j,ik,iorb,jorb
    complex(8)                              :: iw
    complex(8),dimension(Norb,Norb)         :: zeta,fg,self,gdelta
    complex(8),allocatable,dimension(:,:,:) :: gloc
    real(8)                                 :: wm(Lmats),wr(Lreal)

    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    !
    delta=zero
    !
    allocate(gloc(Norb,Norb,Lmats))
    do i=1,Lmats
       zeta=zero
       do iorb=1,Norb
          zeta(iorb,iorb) = (xi*wm(i) + xmu)
          do jorb=1,Norb
             zeta(iorb,jorb) = zeta(iorb,jorb)  - impSmats(1,1,iorb,jorb,i)
          enddo
       enddo
       !
       fg=zero
       do ik=1,Lk
          fg = fg + invert_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
       !
       if(Norb>1)then
          call matrix_inverse(fg)
       else
          fg(1,1)=one/fg(1,1)
       endif

       if(cg_scheme=='weiss')then
          gdelta = fg + impSmats(1,1,:,:,i)
          call matrix_inverse(gdelta)
       else
          gdelta = zeta - hloc(1,1,:,:) - fg
       endif
       delta(:,:,i)=gdelta

    enddo
    !print
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
          call splot("Delta_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,delta(iorb,jorb,:))
       enddo
    enddo
    deallocate(gloc)


    allocate(gloc(Norb,Norb,Lreal))
    do i=1,Lreal
       zeta=zero
       do iorb=1,Norb
          zeta(iorb,iorb) = dcmplx(wr(i),eps) + xmu
          do jorb=1,Norb
             zeta(iorb,jorb) = zeta(iorb,jorb)  - impSreal(1,1,iorb,jorb,i)
          enddo
       enddo
       !
       fg=zero
       do ik=1,Lk
          fg = fg + invert_gk(zeta,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gloc(:,:,i) = fg
    enddo
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,gloc(iorb,jorb,:))
          call splot("DOS_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,-dimag(gloc(iorb,jorb,:))/pi)
       enddo
    enddo
    deallocate(gloc)

    nobj=ed_dens(1)
    if(Norb>1)then
       if(ntype==1)then
          nobj=ed_dens(1)
       elseif(ntype==2)then
          nobj=ed_dens(2)
       else
          nobj=sum(ed_dens)
       endif
    endif
  end subroutine get_delta



  !+----------------------------------------+

  ! function inverse_g0k(iw,hk) result(g0k)
  !   integer                     :: i,M
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8)                  :: iw
  !   complex(8),dimension(2,2)   :: g0k
  !   complex(8)                  :: delta,ppi,vmix
  !   g0k=zero
  !   delta = iw - hk(1,1)
  !   ppi   = iw - hk(2,2)
  !   vmix  = -hk(1,2)
  !   g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
  !   g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
  !   g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
  !   g0k(2,1) = conjg(g0k(1,2))
  ! end function inverse_g0k

  ! function inverse_gk(zeta,hk) result(gk)
  !   integer                     :: i,M
  !   complex(8),dimension(2,2)   :: hk
  !   complex(8),dimension(2)     :: zeta
  !   complex(8),dimension(2,2)   :: gk
  !   complex(8)                  :: delta,ppi,vmix
  !   gk=zero
  !   delta = zeta(1) - hk(1,1)
  !   ppi   = zeta(2) - hk(2,2)
  !   vmix  = -hk(1,2)
  !   gk(1,1) = one/(delta - abs(vmix)**2/ppi)
  !   gk(2,2) = one/(ppi - abs(vmix)**2/delta)
  !   gk(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
  !   gk(2,1) = conjg(gk(1,2))
  ! end function inverse_gk

  function invert_gk(zeta,Hk) result(invg)
    complex(8),dimension(Norb,Norb) :: zeta,Hk
    complex(8),dimension(Norb,Norb) :: invg
    invg=zeta-Hk
    if(Norb>1)then
       call matrix_inverse(invg)
    else
       invg(1,1)=one/invg(1,1)
    endif
  end function invert_gk

  !include 'search_mu.f90'

end program lancED



