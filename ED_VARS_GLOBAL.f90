MODULE ED_VARS_GLOBAL  
  USE ED_BATH_TYPE
  USE MATRIX_SPARSE
  USE EIGEN_SPACE
#ifdef _MPI
  USE MPI
#endif
  implicit none

  !SIZE OF THE PROBLEM
  !Ns   = # of levels per spin
  !Ntot = 2*Ns = total #  of levels
  !NN   = 2**Ntot = 2**(2Ns) max size of the Hilbert space
  !Nbo  =# number of bath sites (all sites - impurity sites)
  !Nsect=# of sectors
  !=========================================================
  integer                                     :: Ns,Ntot,NN
  integer                                     :: Nsect
  integer                                     :: Nbo


  !local part of the Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable   :: Hloc           !local hamiltonian


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)          :: getsector
  integer,allocatable,dimension(:,:)          :: getCsector
  integer,allocatable,dimension(:,:)          :: getCDGsector
  integer,allocatable,dimension(:,:)          :: getBathStride
  integer,allocatable,dimension(:,:)          :: impIndex
  integer,allocatable,dimension(:)            :: getdim,getnup,getndw,getsz


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
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSmats,impSAmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impSreal,impSAreal


  !Density and double occupancy
  !=========================================================
  real(8),dimension(:),allocatable            ::  ed_dens,ed_docc,ed_phisc


#ifdef _MPI
  !MPI Parallel environment variables
  !=========================================================
  integer                                     :: ED_MPI_ID=0
  integer                                     :: ED_MPI_SIZE=1
  integer                                     :: ED_MPI_ERR
#endif  

END MODULE ED_VARS_GLOBAL
