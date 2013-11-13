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
  USE IOTOOLS, only:free_unit,reg,splot
  !LOCAL
  USE ED_BATH_TYPE
  USE MATRIX_SPARSE
  USE EIGEN_SPACE
#ifdef _MPI
  USE MPI
#endif
  implicit none

  !GIT VERSION
  include "revision.inc"

  !SIZE OF THE PROBLEM (as PARAMETERS to allocate in the stack)
  !VARIABLES and DEFINITIONS in common to every subroutine
  !Ns=numero di siti del sistema (ogni sito puo' contenere 2e^-)
  !Ntot=2*Ns=numero totale di particelle
  !Lo spazio di Hilbert e' dato dal prodotto tensore 2^Ns * 2^Ns=2^N =NN (LARGE)
  !Nbath=Ns-1=numero di siti del bagno (sistema - impurezza)
  !NP=dimensione del sottospazio piu' grande (settore).
  !=========================================================
  integer :: Ns,Nbo,Norb,Nspin,Nbath,Ntot,NN,Nsect


  !Global variables
  !=========================================================
  integer                   :: nloop          !max dmft loop variables
  real(8),dimension(3)      :: Uloc           !local interactions
  real(8)                   :: Ust,Jh         !intra-orbitals interactions
  real(8),dimension(3)      :: eloc           !local energies
  complex(8),dimension(3,3) :: hybrd          !hybridizations
  real(8),dimension(3,3)    :: rehybrd        !hybridizations
  real(8),dimension(3,3)    :: imhybrd        !hybridizations
  real(8)                   :: xmu            !chemical potential
  real(8)                   :: beta           !inverse temperature
  real(8)                   :: eps            !broadening
  real(8)                   :: wini,wfin      !
  integer                   :: Nsuccess       !
  logical                   :: Jhflag         !spin-exchange and pair-hopping flag.
  logical                   :: chiflag        !
  logical                   :: HFmode         !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                   :: cutoff         !cutoff for spectral summation
  real(8)                   :: eps_error      !
  integer                   :: lanc_niter     !Max number of Lanczos iterations
  integer                   :: lanc_neigen    !Max number of required eigenvalues per sector
  integer                   :: lanc_ngfiter   !Max number of iteration in resolvant tri-diagonalization
  integer                   :: lanc_nstates   !Max number of states hold in the finite T calculation
  integer                   :: cg_Niter       !Max number of iteration in the fit
  real(8)                   :: cg_Ftol        !Tolerance in the cg fit
  integer                   :: cg_Weight      !CGfit mode 0=normal,1=1/n weight, 2=1/w weight
  character(len=5)          :: cg_Scheme      !fit scheme: delta (default), weiss for G0^
  logical                   :: finiteT        !flag for finite temperature calculation
  character(len=4)          :: ed_method      !flag to set ed method solution: lanc=lanczos method, full=full diagonalization
  character(len=1)          :: ed_type        !flag to set real or complex Ham: d=symmetric H (real), c=hermitian H (cmplx)
  character(len=7)          :: bath_type      !flag to set bath type: irreducible (1bath/imp), reducible(1bath)

  !Dimension of the functions:
  !=========================================================
  integer :: NL,Ltau,Nw,Nfit


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)   :: getsector
  integer,allocatable,dimension(:,:)   :: getCsector
  integer,allocatable,dimension(:,:)   :: getCDGsector
  integer,allocatable,dimension(:,:)   :: getBathStride
  integer,allocatable,dimension(:,:)   :: impIndex
  integer,allocatable,dimension(:)     :: getdim,getnup,getndw


  !Effective Bath used in the ED code (this is opaque to user)
  !=========================================================
  type(effective_bath) :: dmft_bath


  !Eigenvalues,Eigenvectors FULL DIAGONALIZATION
  !=========================================================
  type(full_espace),dimension(:),allocatable :: espace


  !Variables for DIAGONALIZATION
  !=========================================================
  integer                                :: numgs
  type(sparse_espace)                    :: state_list
  type(sparse_matrix)                    :: spH0
  integer,allocatable,dimension(:)       :: neigen_sector

  !Partition function
  !=========================================================
  real(8) :: zeta_function


  !Functions for GETGFUNX
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:) :: impSmats
  complex(8),allocatable,dimension(:,:,:,:) :: impSreal


  !Variables for fixed density mu-loop 
  !=========================================================
  real(8) :: nread,nerr,ndelta

  !Qties needed to get energy
  !=========================================================

  real(8),dimension(:),allocatable   ::  nimp,dimp,nupimp,ndwimp,magimp
  real(8),dimension(:,:),allocatable ::  m2imp


  !NML READ/WRITE UNITS
  character(len=32) :: Hfile,Ofile,GFfile,CHIfile
  integer           :: LOGfile


  !Writing/Reading Units:
  character(len=32) :: USEDinput


  namelist/EDvars/Norb,Nbath,Nspin,&       
       beta,xmu,nloop,uloc,Ust,Jh,Eloc,rehybrd,imhybrd,   &
       eps,wini,wfin,      &
       NL,Nw,Ltau,Nfit,         &
       nread,nerr,ndelta,       &
       chiflag,Jhflag,cutoff,HFmode,   &
       eps_error,Nsuccess,      &
       ed_method,ed_type,bath_type,&
       lanc_neigen,lanc_niter,lanc_ngfiter,lanc_nstates,&
       cg_niter,cg_ftol,cg_weight,cg_scheme,  &
       Hfile,Ofile,GFfile,CHIfile,LOGfile

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine read_input(INPUTunit)
    character(len=*) :: INPUTunit
    logical          :: control

    if(mpiID==0)call version(revision)

    !DEFAULT VALUES OF THE PARAMETERS:
    !ModelConf
    Norb       = 1
    Nbath      = 4
    Nspin      = 1
    Uloc       = [ 2.d0, 0.d0, 0.d0 ]
    Ust        = 0.d0
    Jh         = 0.d0
    reHybrd    = 0.d0
    imHybrd    = 0.d0
    Eloc       = 0.d0
    xmu        = 0.d0
    beta       = 500.d0
    !Loops
    nloop      = 100
    chiflag    =.true.
    Jhflag     =.false.
    hfmode     =.true.
    !parameters
    NL         = 2000
    Nw         = 2000
    Ltau       = 1000
    Nfit       = 1000
    eps        = 0.01d0
    nread      = 0.d0
    nerr       = 1.d-4
    ndelta     = 0.1d0
    wini       =-4.d0
    wfin       = 4.d0
    cutoff     = 1.d-9
    eps_error  = 1.d-5
    nsuccess   = 2
    lanc_niter = 512
    lanc_neigen = 1
    lanc_ngfiter = 100
    lanc_nstates = 1            !set to T=0 calculation
    cg_niter   = 200
    cg_Ftol     = 1.d-9
    cg_weight     = 0
    cg_scheme     = 'delta'
    ed_method    = 'lanc'
    ed_type = 'd'
    bath_type='normal' !hybrid,superc

    !ReadUnits
    Hfile  ="hamiltonian.restart"
    GFfile ="impG"
    CHIfile ="Chi"
    Ofile  ="observables"
    LOGfile=6

    inquire(file=INPUTunit,exist=control)    
    if(control)then
       open(50,file=INPUTunit,status='old')
       read(50,nml=EDvars)
       close(50)
    else
#ifdef _MPI
       if(mpiID==0)then
#endif
          print*,"Can not find INPUT file"
          print*,"Printing a default version in default."//INPUTunit
          open(50,file="default."//INPUTunit)
          write(50,nml=EDvars)
          write(50,*)""
#ifdef _MPI
       endif
#endif
       stop
    endif

    call parse_cmd_variable(Norb,"NORB")    
    call parse_cmd_variable(Nbath,"NBATH")
    call parse_cmd_variable(Nspin,"NSPIN")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(uloc,"ULOC")
    call parse_cmd_variable(ust,"UST")
    call parse_cmd_variable(Jh,"JH")
    call parse_cmd_variable(eloc,"ELOC")
    call parse_cmd_variable(reHybrd,"reHybrd")
    call parse_cmd_variable(imHybrd,"imHybrd")
    call parse_cmd_variable(nloop,"NLOOP")
    call parse_cmd_variable(eps_error,"EPS_ERROR")
    call parse_cmd_variable(nsuccess,"NSUCCESS")
    call parse_cmd_variable(NL,"NL")
    call parse_cmd_variable(Nw,"NW")
    call parse_cmd_variable(Ltau,"LTAU")
    call parse_cmd_variable(Nfit,"NFIT")
    call parse_cmd_variable(nread,"NREAD")
    call parse_cmd_variable(nerr,"NERR")
    call parse_cmd_variable(ndelta,"NDELTA")
    call parse_cmd_variable(wini,"WINI")
    call parse_cmd_variable(wfin,"WFIN")
    call parse_cmd_variable(chiflag,"CHIFLAG")
    call parse_cmd_variable(hfmode,"HFMODE")
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(cutoff,"CUTOFF")
    call parse_cmd_variable(lanc_neigen,"LANC_NEIGEN")
    call parse_cmd_variable(lanc_niter,"LANC_NITER")
    call parse_cmd_variable(lanc_nstates,"LANC_NSTATES")
    call parse_cmd_variable(lanc_ngfiter,"LANC_NGFITER")
    call parse_cmd_variable(cg_niter,"CG_NITER")
    call parse_cmd_variable(cg_scheme,"CG_SCHEME")
    call parse_cmd_variable(cg_ftol,"CG_FTOL")
    call parse_cmd_variable(cg_weight,"CG_WEIGHT")
    call parse_cmd_variable(ed_Type,"ED_TYPE")
    call parse_cmd_variable(ed_Method,"ED_METHOD")
    call parse_cmd_variable(bath_type,"BATH_TYPE")
    call parse_cmd_variable(Hfile,"HFILE")
    call parse_cmd_variable(Ofile,"OFILE")
    call parse_cmd_variable(GFfile,"GFFILE")
    call parse_cmd_variable(CHIfile,"CHIFILE")
    call parse_cmd_variable(LOGfile,"LOGFILE")

#ifdef _MPI
    if(mpiID==0)then
#endif
       open(50,file="used."//INPUTunit)
       write(50,nml=EDvars)
       close(50)
       write(*,*)"CONTROL PARAMETERS"
       write(*,nml=EDvars)
#ifdef _MPI
    endif
#endif

  end subroutine read_input

END MODULE ED_VARS_GLOBAL
