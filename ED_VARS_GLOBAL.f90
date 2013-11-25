MODULE ED_BATH_TYPE
  implicit none
  type effective_bath
     real(8),dimension(:,:,:),allocatable :: e
     real(8),dimension(:,:,:),allocatable :: v
     logical                              :: status=.false.
  end type effective_bath
END MODULE ED_BATH_TYPE
!
MODULE ED_VARS_GLOBAL
  USE SCIFOR_VERSION
  USE COMMON_VARS
  USE PARSE_CMD
  USE IOTOOLS, only:free_unit,reg,splot,file_length
  !LOCAL
  USE ED_BATH_TYPE
  USE MATRIX_SPARSE
  USE EIGEN_SPACE
#ifdef _MPI
  USE MPI
#endif
  implicit none

  !GIT VERSION
  !this file is generated at compilation time in the Makefile
  include "revision.inc"

  !SIZE OF THE PROBLEM
  !Ns   = # of levels per spin
  !Ntot = 2*Ns = total #  of levels
  !NN   = 2**Ntot = 2**(2Ns) max size of the Hilbert space
  !Nbath=# of bath sites (per orbital or not depending on bath_type)
  !Norb =# of impurity orbitals
  !Nbo  =# number of bath sites (all sites - impurity sites)
  !Nsect=# of sectors
  !Nspin=# spin degeneracy (max 2)
  !=========================================================
  integer                                     :: Ns,Ntot,NN
  integer                                     :: Nbath,Nsect
  integer                                     :: Norb,Nspin,Nbo


  !Global variables
  !=========================================================
  integer                                     :: nloop          !max dmft loop variables
  real(8)                                     :: Ust,Jh         !intra-orbitals interactions
  real(8),dimension(3)                        :: Uloc           !local interactions
  real(8),dimension(:,:,:,:),allocatable      :: reHloc         !local hamiltonian, real part 
  real(8),dimension(:,:,:,:),allocatable      :: imHloc         !local hamiltonian, imag part
  complex(8),dimension(:,:,:,:),allocatable   :: Hloc           !local hamiltonian
  real(8)                                     :: xmu            !chemical potential
  real(8)                                     :: beta           !inverse temperature
  real(8)                                     :: eps            !broadening
  real(8)                                     :: wini,wfin      !
  integer                                     :: Nsuccess       !
  logical                                     :: Jhflag         !spin-exchange and pair-hopping flag.
  logical                                     :: chiflag        !
  logical                                     :: HFmode         !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                                     :: cutoff         !cutoff for spectral summation
  real(8)                                     :: dmft_error     !dmft convergence threshold
  integer                                     :: lanc_niter     !Max number of Lanczos iterations
  integer                                     :: lanc_neigen    !Max number of required eigenvalues per sector
  integer                                     :: lanc_ngfiter   !Max number of iteration in resolvant tri-diagonalization
  integer                                     :: lanc_nstates   !Max number of states hold in the finite T calculation
  integer                                     :: cg_Niter       !Max number of iteration in the fit
  real(8)                                     :: cg_Ftol        !Tolerance in the cg fit
  integer                                     :: cg_Weight      !CGfit mode 0=normal,1=1/n weight, 2=1/w weight
  character(len=5)                            :: cg_Scheme      !fit scheme: delta (default), weiss for G0^
  logical                                     :: finiteT        !flag for finite temperature calculation
  character(len=4)                            :: ed_method      !flag to set ed method solution: lanc=lanczos method, full=full diagonalization
  character(len=1)                            :: ed_type        !flag to set real or complex Ham: d=symmetric H (real), c=hermitian H (cmplx)
  character(len=7)                            :: bath_type      !flag to set bath type: irreducible (1bath/imp), reducible(1bath)
  real(8)                                     :: nread          !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)                                     :: nerr           !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)                                     :: ndelta         !initial chemical potential step
  integer                                     :: niter


  !Dimension of the functions:
  !=========================================================
  integer                                     :: NL,Ltau,Nw,Nfit


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)          :: getsector
  integer,allocatable,dimension(:,:)          :: getCsector
  integer,allocatable,dimension(:,:)          :: getCDGsector
  integer,allocatable,dimension(:,:)          :: getBathStride
  integer,allocatable,dimension(:,:)          :: impIndex
  integer,allocatable,dimension(:)            :: getdim,getnup,getndw


  !Effective Bath used in the ED code (this is opaque to user)
  !=========================================================
  type(effective_bath)                        :: dmft_bath


  !Eigenvalues,Eigenvectors FULL DIAGONALIZATION
  !=========================================================
  type(full_espace),dimension(:),allocatable  :: espace


  !Variables for DIAGONALIZATION
  !=========================================================
  integer                                     :: numgs
  type(sparse_espace)                         :: state_list
  type(sparse_matrix)                         :: spH0
  integer,allocatable,dimension(:)            :: neigen_sector

  !Partition function
  !=========================================================
  real(8)                                     :: zeta_function

  !Local Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSreal


  !Global observables
  !=========================================================
  real(8),dimension(:),allocatable            ::  nimp,dimp


  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=32)                           :: Hfile
  integer                                     :: LOGfile


  namelist/EDvars/Norb,Nbath,Nspin,&       
       beta,xmu,nloop,Uloc,Ust,Jh,& 
       eps,wini,wfin,      &
       NL,Nw,Nfit,         &
       nread,nerr,ndelta,       &
       chiflag,Jhflag,cutoff,HFmode,   &
       dmft_error,Nsuccess,      &
       ed_method,ed_type,bath_type,&
       lanc_neigen,lanc_niter,lanc_ngfiter,lanc_nstates,&
       cg_niter,cg_ftol,cg_weight,cg_scheme,  &
       Hfile,LOGfile

END MODULE ED_VARS_GLOBAL
