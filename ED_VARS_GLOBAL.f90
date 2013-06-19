MODULE ED_VARS_GLOBAL
  USE SCIFOR_VERSION
  USE COMMON_VARS
  USE TIMER
  USE PARSE_CMD
  USE IOTOOLS
  USE MATRIX, only: matrix_diagonalize,matrix_inverse
  USE OPTIMIZE
  USE TOOLS, only: arange,linspace
  !LOCAL
  USE MATRIX_SPARSE
  USE EIGEN_SPACE
  USE PLAIN_LANCZOS
  USE ARPACK_LANCZOS
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
  integer :: Ns,Norb,Nspin,Nbath,Ntot,NN,Nsect


  !Global variables
  !=========================================================
  integer :: nloop          !max dmft loop variables
  real(8) :: u              !local,non-local interaction
  real(8) :: xmu            !chemical potential
  real(8) :: beta           !inverse temperature
  real(8) :: eps            !broadening
  real(8) :: tpd,ep0        !Nimp=2 variables
  real(8) :: wini,wfin      !
  integer :: Nsuccess       !
  real(8) :: weight         !
  real(8) :: heff           !
  logical :: chiflag        !
  logical :: HFmode         !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8) :: cutoff         !cutoff for spectral summation
  real(8) :: eps_error      !
  integer :: nLancitermax   !Max number of Lanczos iterations
  integer :: nGFitermax     !Max number of iteration in resolvant tri-diagonalization
  integer :: cgNitmax       !Max number of iteration in the fit
  real(8) :: cgFtol         !Tolerance in the cg fit
  integer :: cgType         !CGfit mode 0=normal,1=1/n weight, 2=1/w weight

  !Dimension of the functions:
  !=========================================================
  integer :: NL,Ltau,Nw,Nfit



  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  type HSmap
     integer,dimension(:),pointer :: map
  end type HSmap
  type(HSmap),dimension(:),allocatable :: Hmap
  integer,allocatable,dimension(:,:) :: invHmap
  integer,allocatable,dimension(:,:) :: getsector
  integer,allocatable,dimension(:,:) :: getCsector
  integer,allocatable,dimension(:,:) :: getCDGsector
  integer,allocatable,dimension(:,:) :: impIndex
  integer,allocatable,dimension(:)   :: getdim,getnup,getndw
  integer                            :: startloop,lastloop



  !Eigenvalues,Eigenvectors FULL DIAGONALIZATION
  !=========================================================
  type(eigenspace),dimension(:),allocatable :: espace


  !Ground state variables LANCZOS DIAGONALIZATINO
  !=========================================================
  integer                                 :: numzero
  type(eig_space)                         :: groundstate
  type(sparse_matrix)                     :: spH0


  !Partition function
  !=========================================================
  real(8) :: zeta_function


  !Functions for GETGFUNX (names must be changed)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:) :: Giw,Siw
  complex(8),allocatable,dimension(:,:,:,:) :: Gwr,Swr
  real(8),allocatable,dimension(:,:,:)      :: Chitau
  real(8),allocatable,dimension(:)          :: Chitautot
  complex(8),allocatable,dimension(:,:,:)   :: Chiw
  complex(8),allocatable,dimension(:)       :: Chiwtot


  !Variables for fixed density mu-loop 
  !=========================================================
  real(8) :: nread,nerr,ndelta

  !Qties needed to get energy
  !=========================================================
  real(8),dimension(:),allocatable   ::  nimp,dimp,nupimp,ndwimp,magimp
  real(8),dimension(:,:),allocatable ::  m2imp


  !NML READ/WRITE UNITS
  character(len=32) :: Hfile,Ofile,GMfile,GRfile,CTfile,CWfile
  integer           :: LOGfile


  !Writing/Reading Units:
  character(len=32) :: USEDinput


  namelist/EDvars/Ns,Norb,Nspin,&       
       beta,xmu,nloop,u,tpd,ep0,&
       eps,wini,wfin,heff,      &
       NL,Nw,Ltau,Nfit,         &
       nread,nerr,ndelta,       &
       chiflag,cutoff,HFmode,   &
       eps_error,Nsuccess,      &
       nLancitermax,nGFitermax,&
       cgNitmax,cgFtol,cgType,   &
       Hfile,Ofile,GMfile,GRfile,CTfile,CWfile,LOGfile

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine read_input(INPUTunit)
    character(len=*) :: INPUTunit
    integer          :: i,NP,nup,ndw
    logical          :: control

    !ModelConf
    Ns         = 5
    Norb       = 1
    Nspin      = 1
    u          = 2.d0
    tpd        = 0.d0
    ep0        = 0.d0
    xmu        = 0.d0
    beta       = 50.d0
    !Loops
    nloop      = 100
    chiflag    =.true.
    hfmode     =.true.
    !parameters
    NL         = 1024
    Nw         = 1024
    Ltau       = 500
    Nfit       = 1024
    eps        = 0.005d0
    nread      = 0.d0
    nerr       = 1.d-4
    ndelta     = 0.1d0
    wini       =-4.d0
    wfin       = 4.d0
    heff       = 0.d0
    cutoff     = 1.d-9
    eps_error  = 1.d-5
    nsuccess   = 2
    nLancitermax = 512
    nGFitermax = 100
    cgNitmax   = 1000
    cgFtol     = 1.d-9
    cgType     = 0
    !ReadUnits
    Hfile  ="hamiltonian.restart"
    GMfile ="impG_iw"
    GRfile ="impG_realw"
    CTfile ="Chi_tau"
    CWfile ="Chi_realw"    
    Ofile  ="observables"
    LOGfile=6

    inquire(file=INPUTunit,exist=control)    
    if(control)then
       open(50,file=INPUTunit,status='old')
       read(50,nml=EDvars)
       close(50)
    else
       print*,"Can not find INPUT file"
       print*,"Printing a default version in default."//INPUTunit
       open(50,file="default."//INPUTunit)
       write(50,nml=EDvars)
       write(50,*)""
       stop
    endif

    call parse_cmd_variable(Ns,"NS")
    call parse_cmd_variable(Norb,"NORB")
    call parse_cmd_variable(Nspin,"NSPIN")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(u,"U")
    call parse_cmd_variable(tpd,"TPD")
    call parse_cmd_variable(ep0,"EP0")
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
    call parse_cmd_variable(wini,"WINI","WMIN")
    call parse_cmd_variable(wfin,"WFIN","WMAX")
    call parse_cmd_variable(chiflag,"CHIFLAG")
    call parse_cmd_variable(hfmode,"HFMODE")
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(cutoff,"CUTOFF")
    call parse_cmd_variable(heff,"HEFF")
    call parse_cmd_variable(nGFitermax,"NGFITERMAX")
    call parse_cmd_variable(nLancitermax,"NLANCITERMAX")
    call parse_cmd_variable(cgNitmax,"CGNITMAX")
    call parse_cmd_variable(cgFtol,"CGFTOL")
    call parse_cmd_variable(cgType,"CGTYPE")
    call parse_cmd_variable(Hfile,"HFILE")
    call parse_cmd_variable(Ofile,"OFILE")
    call parse_cmd_variable(GMfile,"GMFILE")
    call parse_cmd_variable(GRfile,"GRFILE")
    call parse_cmd_variable(CTfile,"CTFILE")
    call parse_cmd_variable(CWfile,"CWFILE")
    call parse_cmd_variable(LOGfile,"LOGFILE")

    call version(revision)

    call allocate_system_structure()

    nup=Ns/2
    ndw=Ns-nup
    NP=(factorial(Ns)/factorial(nup)/factorial(Ns-nup))
    NP=NP*(factorial(Ns)/factorial(ndw)/factorial(Ns-ndw))
    write(*,*)"CONTROL PARAMETERS"
    write(*,nml=EDvars)
    write(*,*)"--------------------------------------------"
    write(*,*)'| Total number of sites/spin   = ',Ns
    write(*,*)'| Maximum dimension            = ',NP
    write(*,*)'| Number of impurities         = ',Norb
    write(*,*)'| Bath`s number of sites/spin  = ',Nbath
    write(*,*)'| Total size, Hilber space dim.= ',Ntot,NN
    write(*,*)'| Number of sectors            = ',Nsect
    write(*,*)"--------------------------------------------"
    print*,''
    USEDinput="used."//INPUTunit
    open(50,file=trim(adjustl(trim(USEDinput))))
    write(50,nml=EDvars)
    close(50)

    !Some check:
    if(Nfit>NL)Nfit=NL
    if(Norb>2)call abort("Norb > 2 is not yet supported!")
    if(nerr > eps_error) nerr=eps_error    

    !allocate functions
    allocate(Giw(Norb,Norb,Nspin,NL),Siw(Norb,Norb,Nspin,NL))
    allocate(Gwr(Norb,Norb,Nspin,Nw),Swr(Norb,Norb,Nspin,Nw))
    if(chiflag)allocate(Chitau(Norb,Norb,0:Ltau),Chiw(Norb,Norb,Nw))

    allocate(nimp(Norb),dimp(Norb),nupimp(Norb),ndwimp(Norb),magimp(Norb))
    allocate(m2imp(Norb,Norb))
  end subroutine read_input



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine allocate_system_structure()
    Nbath = Ns-Norb
    Ntot  = 2*Ns
    NN    = 2**Ntot
    Nsect = (Ns+1)*(Ns+1)
    allocate(impIndex(Norb,2))
    allocate(Hmap(Nsect),invHmap(Nsect,NN))
    allocate(getdim(Nsect),getnup(Nsect),getndw(Nsect))
    allocate(getsector(0:Ns,0:Ns))
    allocate(getCsector(2,Nsect))
    allocate(getCDGsector(2,Nsect))
  end subroutine allocate_system_structure


  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the Heaviside  function
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial

END MODULE ED_VARS_GLOBAL
