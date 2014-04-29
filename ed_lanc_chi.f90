!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Spin Susceptibility using Lanczos algorithm
!+------------------------------------------------------------------+
subroutine lanc_ed_getchi()
  integer :: iorb,jorb,ispin
  call allocate_grids
  allocate(Chitau(Norb,0:Ltau),Chiw(Norb,Lreal),Chiiw(Norb,0:Lmats))
  Chitau=0.d0
  Chiw=zero
  Chiiw=zero
  do iorb=1,Norb
     if(ED_MPI_ID==0)write(LOGfile,"(A)")"Evaluating Chi_Orb"//reg(txtfy(iorb))
     select case(ed_type)
     case default
        call lanc_ed_buildchi_d(iorb)
     case ('c')
        call lanc_ed_buildchi_c(iorb)
     end select
  enddo
  Chitau = Chitau/zeta_function
  Chiw   = Chiw/zeta_function
  Chiiw  = Chiiw/zeta_function
  if(ED_MPI_ID==0)call print_imp_chi()
  deallocate(Chitau,Chiw,Chiiw)
  deallocate(wm,wr,tau,vm)
end subroutine lanc_ed_getchi


!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine lanc_ed_buildchi_d(iorb,iverbose)
  integer                          :: iorb,isite,isect0,izero
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  real(8),allocatable              :: vvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !
  if(ED_MPI_ID==0)call start_progress
  do izero=1,numstates
     call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_vec  => es_return_vector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim0  = getdim(isect0)
     allocate(HImap(idim0),vvinit(idim0))
     call build_sector(isect0,HImap)
     vvinit=0.d0
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = real(ib(iorb),8)-real(ib(iorb+Ns),8)
        vvinit(m) = sgn*state_vec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HImap)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call ed_buildH_d(isect0)
     call lanczos_plain_set_htimesv_d(lanc_spHtimesV_dd)
     call lanczos_plain_tridiag_d(vvinit,alfa_,beta_,nitermax)
     call lanczos_plain_delete_htimesv
     call add_to_lanczos_chi(norm0,state_e,nitermax,alfa_,beta_,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_vec)
  enddo
  if(ED_MPI_ID==0)call stop_progress
  deallocate(alfa_,beta_)
end subroutine lanc_ed_buildchi_d


subroutine lanc_ed_buildchi_c(iorb,iverbose)
  integer                          :: iorb,isite,isect0,izero
  integer                          :: numstates
  integer                          :: nlanc,idim0
  integer                          :: iup0,idw0
  integer                          :: ib(Ntot)
  integer                          :: m,i,j,r
  real(8)                          :: norm0,sgn
  real(8),allocatable              :: alfa_(:),beta_(:)
  complex(8),allocatable           :: vvinit(:)
  integer                          :: Nitermax
  logical,optional                 :: iverbose
  logical                          :: iverbose_
  integer,allocatable,dimension(:) :: HImap    !map of the Sector S to Hilbert space H
  !
  iverbose_=.false.;if(present(iverbose))iverbose_=iverbose
  !
  Nitermax=lanc_nGFiter
  allocate(alfa_(Nitermax),beta_(Nitermax))
  idim0  = getdim(isect0)
  !
  numstates=numgs
  if(finiteT)numstates=state_list%size
  !
  if(ED_MPI_ID==0)call start_progress
  do izero=1,numstates
     call progress(izero,numstates)
     isect0     =  es_return_sector(state_list,izero)
     state_e    =  es_return_energy(state_list,izero)
     state_cvec => es_return_cvector(state_list,izero)
     norm0=sqrt(dot_product(state_vec,state_vec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     allocate(HImap(idim0),vvinit(idim0))
     call build_sector(isect0,HImap)
     vvinit=0.d0
     do m=1,idim0                     !loop over |gs> components m
        i=HImap(m)
        call bdecomp(i,ib)
        sgn = dble(ib(iorb))-dble(ib(iorb+Ns))
        vvinit(m) = sgn*state_cvec(m)   !build the cdg_up|gs> state
     enddo
     deallocate(HImap)
     norm0=sqrt(dot_product(vvinit,vvinit))
     vvinit=vvinit/norm0
     alfa_=0.d0 ; beta_=0.d0 ; nlanc=0
     call lanczos_plain_set_htimesv_c(lanc_spHtimesV_cc)
     call ed_buildH_c(isect0)
     call lanczos_plain_tridiag_c(vvinit,alfa_,beta_,nitermax)
     call lanczos_plain_delete_htimesv
     call add_to_lanczos_chi(norm0,state_e,nitermax,alfa_,beta_,iorb)
     deallocate(vvinit)
     if(spH0%status)call sp_delete_matrix(spH0)
     nullify(state_cvec)
  enddo
  if(ED_MPI_ID==0)call stop_progress
  deallocate(alfa_,beta_)
end subroutine lanc_ed_buildchi_c





!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_chi(vnorm,Ei,nlanc,alanc,blanc,iorb)
  real(8)                                    :: vnorm,Ei,Egs,pesoBZ,de,peso
  integer                                    :: nlanc
  real(8),dimension(nlanc)                   :: alanc,blanc 
  integer                                    :: isign,iorb
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw,chisp
  !
  Egs = state_list%emin
  pesoBZ = vnorm**2/zeta_function 
  if(finiteT)pesoBZ = pesoBZ*exp(-beta*(Ei-Egs))
  diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
  forall(i=1:Nlanc)Z(i,i)=1.d0
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
  do j=1,nlanc
     de = diag(j)-Ei
     peso = pesoBZ*Z(1,j)*Z(1,j)
     if(de>cutoff)chiiw(iorb,0)=chiiw(iorb,0) - peso*(exp(-beta*de)-1.d0)/de
     do i=1,Lmats
        iw=xi*vm(i)
        chiiw(iorb,i)=chiiw(iorb,i) + peso*(exp(-beta*de)-1.d0)/(iw-de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        chiw(iorb,i)=chiw(iorb,i) + peso*(exp(-beta*de)-1.d0)/(iw-de)
     enddo
     !
     do i=0,Ltau 
        chitau(iorb,i)=chitau(iorb,i) + peso*exp(-tau(i)*de)
     enddo
  enddo
end subroutine add_to_lanczos_chi
