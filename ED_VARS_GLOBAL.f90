! MODULE ED_BATH_TYPE
!   implicit none
!   type effective_bath
!      real(8),dimension(:,:,:),allocatable :: e
!      real(8),dimension(:,:,:),allocatable :: v
!      logical                              :: status=.false.
!   end type effective_bath
! END MODULE ED_BATH_TYPE
!
MODULE ED_VARS_GLOBAL
  USE SCIFOR_VERSION
  USE COMMON_VARS
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
  integer                                     :: Nsect
  integer                                     :: Nbo


  !local part of the Hamiltonian
  !=========================================================
  real(8),dimension(:,:,:,:),allocatable      :: reHloc         !local hamiltonian, real part 
  real(8),dimension(:,:,:,:),allocatable      :: imHloc         !local hamiltonian, imag part
  complex(8),dimension(:,:,:,:),allocatable   :: Hloc           !local hamiltonian


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


  !Density and double occupancy
  !=========================================================
  real(8),dimension(:),allocatable            ::  nimp,dimp


END MODULE ED_VARS_GLOBAL
