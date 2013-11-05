!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE TIMER
  USE IOTOOLS, only:free_unit,reg
  USE TOOLS, only: arange,linspace
  USE MATRIX, only: matrix_inverse
  USE PLAIN_LANCZOS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  !
  implicit none
  private 

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer :: state_vec
  real(8)                      :: state_e

  !Impurity GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:) :: impGmats
  complex(8),allocatable,dimension(:,:,:,:) :: impGreal


  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)    :: Chitau
  real(8),allocatable,dimension(:)      :: Chitautot
  complex(8),allocatable,dimension(:,:) :: Chiw,Chiiw
  complex(8),allocatable,dimension(:)   :: Chiwtot,Chiiwtot

  public :: lanc_ed_getgf
  public :: full_ed_getgf
  public :: lanc_ed_getchi
  public :: full_ed_getchi

contains



  !####################################################################
  !                    LANCZOS DIAGONALIZATION 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions using Lanczos algorithm
  !+------------------------------------------------------------------+
  subroutine lanc_ed_getgf()
    integer :: izero,iorb,jorb,ispin,i
    integer :: isect0,numstates
    real(8) :: norm0
    call allocate_grids

    if(.not.allocated(impGmats))allocate(impGmats(Nspin,Norb,Norb,NL))
    if(.not.allocated(impGreal))allocate(impGreal(Nspin,Norb,Norb,Nw))
    impGmats=zero
    impGreal=zero

    do ispin=1,Nspin
       do iorb=1,Norb
          if(mpiID==0)write(LOGfile,"(A)")" Get G_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))
          call lanc_ed_buildgf(iorb,ispin,.false.)
       enddo
    enddo

    if(bath_type=='hybrid')then
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                if(mpiID==0)write(LOGfile,"(A)")" Get G_l"//&
                     reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))
                call lanc_ed_buildgf_mix(iorb,jorb,ispin,.false.)
             enddo
          enddo
       enddo

       !Put here off-diagonal manipulation by symmetry:
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             impGmats(:,iorb,jorb,:) = 0.5d0*(impGmats(:,iorb,jorb,:) + &
                  (xi-1.d0)*(impGmats(:,iorb,iorb,:)+impGmats(:,jorb,jorb,:)))
             !
             impGreal(:,iorb,jorb,:) = 0.5d0*(impGreal(:,iorb,jorb,:) + &
                  (xi-1.d0)*(impGreal(:,iorb,iorb,:)+impGreal(:,jorb,jorb,:)))
          enddo
       enddo

       forall(iorb=1:Norb,jorb=1:Norb,jorb<iorb)
          impGmats(:,iorb,jorb,:)=impGmats(:,jorb,iorb,:)
          impGreal(:,iorb,jorb,:)=impGreal(:,jorb,iorb,:)
       end forall
    endif

    ! impGmats=impGmats/zeta_function
    ! impGreal=impGreal/zeta_function

    !Print impurity functions:
    if(mpiID==0)call print_imp_gf
    deallocate(wm,wr,tau,vm)
    deallocate(impGmats,impGreal)
  end subroutine lanc_ed_getgf







  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_buildgf(iorb,ispin,iverbose)
    integer             :: iorb,ispin,isite,isect0,izero
    integer             :: idim0,jsect0
    integer             :: jup0,jdw0,jdim0
    integer             :: ib(Ntot)
    integer             :: m,i,j,r,numstates
    real(8)             :: sgn,norm2,norm0
    complex(8)          :: cnorm2
    real(8),allocatable :: vvinit(:),alfa_(:),beta_(:)
    integer             :: Nitermax,Nlanc
    logical,optional    :: iverbose
    logical             :: iverbose_
    integer,allocatable,dimension(:)     :: HImap,HJmap    !map of the Sector S to Hilbert space H
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    !
    Nitermax=lanc_nGFiter
    allocate(alfa_(Nitermax),beta_(Nitermax))
    !
    isite=impIndex(iorb,ispin)
    !
    numstates=numgs
    if(finiteT)numstates=state_list%size
    !   
    if(mpiID==0)call start_progress
    do izero=1,numstates
       if(mpiID==0)call progress(izero,numstates)
       isect0     =  es_return_sector(state_list,izero)
       state_e    =  es_return_energy(state_list,izero)
       state_vec  => es_return_vector(state_list,izero)
       norm0=sqrt(dot_product(state_vec,state_vec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"

       idim0  = getdim(isect0)
       !allocate map from isect0 to HS
       allocate(HImap(idim0))
       call build_sector(isect0,HImap)
       !
       !ADD ONE PARTICLE:
       jsect0 = getCDGsector(ispin,isect0)
       if(jsect0/=0)then 
          jdim0  = getdim(jsect0)
          jup0   = getnup(jsect0)
          jdw0   = getndw(jsect0)
          if(mpiID==0.AND.iverbose_)write(LOGfile,"(A,2I3,I15)")'add particle:',jup0,jdw0,jdim0
          allocate(HJmap(jdim0),vvinit(jdim0))
          call build_sector(jsect0,HJmap) !note that here you are doing twice the map building...
          vvinit=0.d0
          do m=1,idim0                     !loop over |gs> components m
             i=HImap(m)                    !map m to Hilbert space state i
             call bdecomp(i,ib)            !i into binary representation
             if(ib(isite)==0)then          !if impurity is empty: proceed
                call cdg(isite,i,r,sgn)
                j=binary_search(HJmap,r)      !map r back to  jsect0
                vvinit(j) = sgn*state_vec(m)  !build the cdg_up|gs> state
             endif
          enddo
          deallocate(HJmap)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)

          call lanczos_plain_set_htimesv_d(lanc_spHtimesV_d)
          call setup_Hv_sector(jsect0)
          call ed_geth(jsect0)
          alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
          call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
          call delete_Hv_sector()
          call lanczos_plain_delete_htimesv
          cnorm2=cmplx(norm2,0.d0,8)
          call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,iorb,ispin)
          deallocate(vvinit)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       !
       !REMOVE ONE PARTICLE:
       jsect0 = getCsector(ispin,isect0)
       if(jsect0/=0)then
          jdim0  = getdim(jsect0)
          jup0    = getnup(jsect0)
          jdw0    = getndw(jsect0)
          if(mpiID==0.AND.iverbose_)write(LOGfile,"(A,2I3,I15)")'del particle:',&
               jup0,jdw0,jdim0
          allocate(HJmap(jdim0),vvinit(jdim0))
          call build_sector(jsect0,HJmap)
          vvinit=0.d0
          do m=1,idim0
             i=HImap(m)
             call bdecomp(i,ib)
             if(ib(isite)==1)then
                call c(isite,i,r,sgn)
                j=binary_search(HJmap,r)
                vvinit(j) = sgn*state_vec(m)
             endif
          enddo
          deallocate(HJmap)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          call lanczos_plain_set_htimesv_d(lanc_spHtimesV_d)
          call setup_Hv_sector(jsect0)
          call ed_geth(jsect0)
          alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
          call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
          call delete_Hv_sector()
          call lanczos_plain_delete_htimesv
          cnorm2=cmplx(norm2,0.d0,8)
          call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,iorb,ispin)
          deallocate(vvinit)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       nullify(state_vec)
       deallocate(HImap)
       !
    enddo
    if(mpiID==0)call stop_progress
    deallocate(alfa_,beta_)
  end subroutine lanc_ed_buildgf




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_buildgf_mix(iorb,jorb,ispin,iverbose)
    integer                          :: iorb,jorb,ispin,isite,jsite,isect0,izero
    integer                          :: idim0,jsect0
    integer                          :: jup0,jdw0,jdim0
    integer                          :: ib(Ntot)
    integer                          :: m,i,j,r,numstates
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    real(8),allocatable              :: vvinit(:)
    complex(8),allocatable           :: cvinit(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: Nitermax,Nlanc
    logical,optional                 :: iverbose
    logical                          :: iverbose_
    integer,allocatable,dimension(:) :: HImap,HJmap    !map of the Sector S to Hilbert space H
    !
    iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
    !
    Nitermax=lanc_nGFiter
    allocate(alfa_(Nitermax),beta_(Nitermax))
    isite=impIndex(iorb,ispin)  !orbital 1
    jsite=impIndex(jorb,ispin)  !orbital 2
    !
    numstates=numgs
    if(finiteT)numstates=state_list%size
    !   
    if(mpiID==0)call start_progress
    do izero=1,numstates
       if(mpiID==0)call progress(izero,numstates)
       isect0     =  es_return_sector(state_list,izero)
       state_e    =  es_return_energy(state_list,izero)
       state_vec  => es_return_vector(state_list,izero)
       norm0=sqrt(dot_product(state_vec,state_vec))
       if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
       !
       idim0  = getdim(isect0)
       allocate(HImap(idim0))
       call build_sector(isect0,HImap)
       !
       !EVALUATE (c^+_iorb + c^+_jorb)|gs>
       jsect0 = getCDGsector(ispin,isect0)
       if(jsect0/=0)then 
          jdim0  = getdim(jsect0)
          if(mpiID==0.AND.iverbose_)write(*,"(A,2I3,I15)")'add particle:',getnup(jsect0),getndw(jsect0),jdim0
          allocate(HJmap(jdim0),vvinit(jdim0))
          call build_sector(jsect0,HJmap)
          vvinit=0.d0
          do m=1,idim0
             i=HImap(m)
             call bdecomp(i,ib)
             if(ib(isite)==0)then
                call cdg(isite,i,r,sgn)
                j=binary_search(HJmap,r)
                vvinit(j) = sgn*state_vec(m)
             endif
             if(ib(jsite)==0)then
                call cdg(jsite,i,r,sgn)
                j=binary_search(HJmap,r)
                vvinit(j) = vvinit(j) + sgn*state_vec(m)
             endif
          enddo
          deallocate(HJmap)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          call lanczos_plain_set_htimesv_d(lanc_spHtimesV_d)
          call setup_Hv_sector(jsect0)
          call ed_geth(jsect0)
          alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
          call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
          call delete_Hv_sector()
          call lanczos_plain_delete_htimesv
          cnorm2=dcmplx(norm2,0.d0)
          call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
          deallocate(vvinit)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif

       !EVALUATE (c_iorb + c_jorb)|gs>
       jsect0 = getCsector(ispin,isect0)
       if(jsect0/=0)then
          jdim0   = getdim(jsect0)
          if(mpiID==0.AND.iverbose_)write(*,"(A,2I3,I15)")'del particle:',getnup(jsect0),getndw(jsect0),jdim0
          allocate(HJmap(jdim0),vvinit(jdim0))
          call build_sector(jsect0,HJmap)
          vvinit=0.d0
          do m=1,idim0
             i=HImap(m)
             call bdecomp(i,ib)
             if(ib(isite)==1)then
                call c(isite,i,r,sgn)
                j=binary_search(HJmap,r)
                vvinit(j) = sgn*state_vec(m)
             endif
             if(ib(jsite)==1)then
                call c(jsite,i,r,sgn)
                j=binary_search(HJmap,r)
                vvinit(j) = vvinit(j) + sgn*state_vec(m)
             endif
          enddo
          deallocate(HJmap)
          norm2=dot_product(vvinit,vvinit)
          vvinit=vvinit/sqrt(norm2)
          call lanczos_plain_set_htimesv_d(lanc_spHtimesV_d)
          call setup_Hv_sector(jsect0)
          call ed_geth(jsect0)
          alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
          call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nlanc)
          call delete_Hv_sector()
          call lanczos_plain_delete_htimesv
          cnorm2=dcmplx(norm2,0.d0)
          call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
          deallocate(vvinit)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif



       !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
       jsect0 = getCDGsector(ispin,isect0)
       if(jsect0/=0)then 
          jdim0  = getdim(jsect0)
          if(mpiID==0.AND.iverbose_)write(*,"(A,2I3,I15)")'add particle:',getnup(jsect0),getndw(jsect0),jdim0
          allocate(HJmap(jdim0),cvinit(jdim0))
          call build_sector(jsect0,HJmap)
          cvinit=zero
          do m=1,idim0
             i=HImap(m)
             call bdecomp(i,ib)
             if(ib(isite)==0)then
                call cdg(isite,i,r,sgn)
                j=binary_search(HJmap,r)
                cvinit(j) = sgn*state_vec(m)
             endif
             if(ib(jsite)==0)then
                call cdg(jsite,i,r,sgn)
                j=binary_search(HJmap,r)
                cvinit(j) = cvinit(j) + xi*sgn*state_vec(m)
             endif
          enddo
          deallocate(HJmap)
          norm2=dot_product(cvinit,cvinit)
          cvinit=cvinit/sqrt(norm2)
          call lanczos_plain_set_htimesv_c(lanc_spHtimesV_c)
          call setup_Hv_sector(jsect0)
          call ed_geth(jsect0)
          alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
          call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc)
          call delete_Hv_sector()
          call lanczos_plain_delete_htimesv
          cnorm2=dcmplx(0.d0,-norm2)
          call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,1,iorb,jorb,ispin)
          deallocate(cvinit)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !

       !EVALUATE (c_iorb - xi*c_jorb)|gs>
       jsect0 = getCsector(ispin,isect0)
       if(jsect0/=0)then
          jdim0   = getdim(jsect0)
          if(mpiID==0.AND.iverbose_)write(*,"(A,2I3,I15)")'del particle:',getnup(jsect0),getndw(jsect0),jdim0
          allocate(HJmap(jdim0),cvinit(jdim0))
          call build_sector(jsect0,HJmap)
          cvinit=zero
          do m=1,idim0
             i=HImap(m)
             call bdecomp(i,ib)
             if(ib(isite)==1)then
                call c(isite,i,r,sgn)
                j=binary_search(HJmap,r)
                cvinit(j) = sgn*state_vec(m)
             endif
             if(ib(jsite)==1)then
                call c(jsite,i,r,sgn)
                j=binary_search(HJmap,r)
                cvinit(j) = cvinit(j) - xi*sgn*state_vec(m)
             endif
          enddo
          deallocate(HJmap)
          norm2=dot_product(cvinit,cvinit)
          cvinit=cvinit/sqrt(norm2)
          call lanczos_plain_set_htimesv_c(lanc_spHtimesV_c)
          call setup_Hv_sector(jsect0)
          call ed_geth(jsect0)
          alfa_=0.d0 ; beta_=0.d0 ; nlanc=nitermax
          call lanczos_plain_tridiag_c(cvinit,alfa_,beta_,nlanc)
          call delete_Hv_sector()
          call lanczos_plain_delete_htimesv
          cnorm2=dcmplx(0.d0,-norm2)
          call add_to_lanczos_gf(cnorm2,state_e,nlanc,alfa_,beta_,-1,iorb,jorb,ispin)
          deallocate(cvinit)
          if(spH0%status)call sp_delete_matrix(spH0)
       endif
       !
       nullify(state_vec)
       deallocate(HImap)
       !
    enddo
    if(mpiID==0)call stop_progress
    deallocate(alfa_,beta_)
  end subroutine lanc_ed_buildgf_mix





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine add_to_lanczos_gf(vnorm2,Ei,nlanc,alanc,blanc,isign,iorb,jorb,ispin)
    complex(8)                                 :: vnorm2
    real(8)                                    :: Ei,Egs,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(nlanc)                   :: alanc,blanc 
    integer                                    :: isign,iorb,jorb,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    pesoBZ = vnorm2/zeta_function
    if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,NL
          iw=xi*wm(i)
          impGmats(ispin,iorb,jorb,i)=impGmats(ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Nw
          iw=cmplx(wr(i),eps,8)
          impGreal(ispin,iorb,jorb,i)=impGreal(ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf





  !####################################################################
  !                    LANC GET SUSCPTIBILITY
  !####################################################################
  include 'ed_lanc_getchi.f90'






  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  include 'ed_full_getgf_chi.f90'






  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    allocate(wm(NL))
    wm     = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(vm(0:NL))          !bosonic frequencies
    do i=0,NL
       vm(i) = pi/beta*2.d0*real(i,8)
    enddo
    allocate(wr(Nw))
    wr     = linspace(wini,wfin,Nw)
    allocate(tau(0:Ltau))
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_gf
    integer                                  :: i,j,ispin,unit(6),iorb,jorb
    complex(8)                               :: iw
    complex(8),dimension(Nspin,Norb,Norb,NL) :: impG0iw
    complex(8),dimension(Nspin,Norb,Norb,Nw) :: impG0wr
    complex(8),dimension(Norb,Norb)          :: invGimp,impG0
    character(len=20)                        :: suffix
    write(LOGfile,"(A)")"Printing the impurity GF"
    impSmats = zero
    impSreal = zero
    select case(bath_type)
    case default
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,NL
                iw=xi*wm(i)
                impG0iw(ispin,iorb,iorb,i) = iw+xmu-eloc(iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
                impSmats(ispin,iorb,iorb,i)= impG0iw(ispin,iorb,iorb,i) - one/impGmats(ispin,iorb,iorb,i)
                impG0iw(ispin,iorb,iorb,i) = delta_bath(ispin,iorb,iw,dmft_bath)
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                impG0wr(ispin,iorb,iorb,i) = wr(i)+xmu-eloc(iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
                impSreal(ispin,iorb,iorb,i)= impG0wr(ispin,iorb,iorb,i) - one/impGreal(ispin,iorb,iorb,i)
                impG0wr(ispin,iorb,iorb,i) = delta_bath(ispin,iorb,iw,dmft_bath)
             enddo
          enddo
       enddo
       !Print the impurity functions:
       do iorb=1,Norb
          suffix="_orb"//reg(txtfy(iorb))
          call open_units(reg(suffix))
          do i=1,NL
             write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impGmats(ispin,iorb,iorb,i)),dreal(impGmats(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(3),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impSmats(ispin,iorb,iorb,i)),dreal(impSmats(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impG0iw(ispin,iorb,iorb,i)),dreal(impG0iw(ispin,iorb,iorb,i)),ispin=1,Nspin)
          enddo
          do i=1,Nw
             write(unit(2),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impGreal(ispin,iorb,iorb,i)),dreal(impGreal(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impSreal(ispin,iorb,iorb,i)),dreal(impSreal(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(6),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impG0wr(ispin,iorb,iorb,i)),dreal(impG0wr(ispin,iorb,iorb,i)),ispin=1,Nspin)
          enddo
          call close_units
       enddo

       !=====================================================================================
    case ('hybrid')
       do ispin=1,Nspin
          do iorb=1,Norb
             !Get WF diagonals:
             do i=1,NL
                iw=xi*wm(i)
                impG0iw(ispin,iorb,iorb,i)= iw + xmu-eloc(iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                impG0wr(ispin,iorb,iorb,i)= wr(i)+xmu-eloc(iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
             enddo
             !Get WF off-diagonals
             do jorb=iorb+1,Norb
                do i=1,NL
                   iw=xi*wm(i)
                   impG0iw(ispin,iorb,jorb,i)= -delta_bath(ispin,iorb,jorb,iw,dmft_bath)
                enddo
                do i=1,Nw
                   iw=cmplx(wr(i),eps)
                   impG0wr(ispin,iorb,jorb,i)= -delta_bath(ispin,iorb,jorb,iw,dmft_bath)
                enddo
             enddo
          enddo
       enddo
       !
       forall(iorb=1:Norb,jorb=1:Norb,jorb<iorb)
          impG0iw(:,iorb,jorb,:)=impG0iw(:,jorb,iorb,:)
          impG0wr(:,iorb,jorb,:)=impG0wr(:,jorb,iorb,:)
       end forall
       !
       !get Sigma and WF by matrix inversions:
       do ispin=1,Nspin
          do i=1,NL
             invGimp = impGmats(ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSmats(ispin,:,:,i) = impG0iw(ispin,:,:,i) - invGimp
          enddo
          do i=1,Nw
             invGimp = impGreal(ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSreal(ispin,:,:,i) = impG0wr(ispin,:,:,i) - invGimp
          enddo
          do i=1,NL
             impG0=impG0iw(ispin,:,:,i)
             call matrix_inverse(impG0)
             impG0iw(ispin,:,:,i)=impG0
          enddo
          do i=1,Nw
             impG0=impG0wr(ispin,:,:,i)
             call matrix_inverse(impG0)
             impG0wr(ispin,:,:,i)=impG0
          enddo
       enddo
       !
       !Print the impurity functions:
       do iorb=1,Norb
          do jorb=iorb,Norb
             suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
             call open_units(reg(suffix))
             do i=1,NL
                write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impGmats(ispin,iorb,jorb,i)),dreal(impGmats(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(3),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSmats(ispin,iorb,jorb,i)),dreal(impSmats(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impG0iw(ispin,iorb,jorb,i)),dreal(impG0iw(ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Nw
                write(unit(2),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impGreal(ispin,iorb,jorb,i)),dreal(impGreal(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSreal(ispin,iorb,jorb,i)),dreal(impSreal(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(6),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impG0wr(ispin,iorb,jorb,i)),dreal(impG0wr(ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             call close_units()
          enddo
       enddo
    end select

  contains

    subroutine open_units(string)
      character(len=*) :: string
      unit(1)=free_unit()
      open(unit(1),file="impG"//string//"_iw.ed")
      unit(2)=free_unit()
      open(unit(2),file="impG"//string//"_realw.ed")
      unit(3)=free_unit()
      open(unit(3),file="impSigma"//string//"_iw.ed")
      unit(4)=free_unit()
      open(unit(4),file="impSigma"//string//"_realw.ed")
      unit(5)=free_unit()
      open(unit(5),file="impDelta"//string//"_iw.ed")
      unit(6)=free_unit()
      open(unit(6),file="impDelta"//string//"_realw.ed")
    end subroutine open_units

    subroutine close_units()
      do i=1,6
         close(unit(i))
      enddo
    end subroutine close_units

  end subroutine print_imp_gf




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_chi
    integer                               :: i,j,iorb
    integer                               :: unit(3)
    do iorb=1,Norb
       unit(1)=free_unit()
       open(unit(1),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_tau.ed")
       unit(2)=free_unit()
       open(unit(2),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_realw.ed")
       unit(3)=free_unit()
       open(unit(3),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_iw.ed")
       do i=0,Ltau/2
          write(unit(1),*)tau(i),chitau(iorb,i)
       enddo
       do i=1,Nw
          if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(iorb,i)),dreal(chiw(iorb,i))
       enddo
       do i=0,NL
          write(unit(3),*)vm(i),dimag(chiiw(iorb,i)),dreal(chiiw(iorb,i))
       enddo
       close(unit(1))
       close(unit(2))
       close(unit(3))
    enddo
    ! unit(1)=free_unit()
    ! open(unit(1),file=trim(CHIfile)//"_tot_tau.ed")
    ! unit(2)=free_unit()
    ! open(unit(2),file=trim(CHIfile)//"_tot_realw.ed")
    ! unit(3)=free_unit()
    ! open(unit(3),file=trim(CHIfile)//"_tot_iw.ed")
    ! do i=0,Ltau
    !    write(unit(1),*)tau(i),chitautot(i)
    ! enddo
    ! do i=1,Nw
    !    if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiwtot(i)),dreal(chiwtot(i))
    ! enddo
    ! do i=0,NL
    !    write(unit(3),*)vm(i),dimag(chiiwtot(i)),dreal(chiiwtot(i))
    ! enddo
    ! close(unit(1))
    ! close(unit(2))
    ! close(unit(3))
  end subroutine print_imp_chi

end MODULE ED_GREENS_FUNCTIONS
