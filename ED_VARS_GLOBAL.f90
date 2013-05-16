MODULE ED_VARS_GLOBAL
  USE COMMON_VARS
  USE TIMER
  USE PARSE_CMD
  USE SLPLOT
  USE FUNCTIONS
  USE MATRIX
  USE TOOLS
  implicit none

  !GIT VERSION
  include "revision.inc"

  !SIZE OF THE PROBLEM (as PARAMETERS to allocate in the stack)
  !VARIABLES and DEFINITIONS in common to every subroutine
  !Ns=numero di siti del sistema (ogni sito puo' contenere 2e^-)
  !N=2*Ns=numero totale di particelle
  !Lo spazio di Hilbert e' dato dal prodotto tensore 2^Ns * 2^Ns=2^N =NN (LARGE)
  !Nbath=Ns-1=numero di siti del bagno (sistema - impurezza)
  !NP=dimensione del sottospazio piu' grande (settore).
  !=========================================================
  integer :: Ns,Nimp,Nspin,Nbath,NP,N,NN,Nsect


  !Global variables
  !=========================================================
  integer :: nloop    !dmft loop variables
  real(8) :: d              !bandwidth
  real(8) :: ts,tsp,tpp,tdd !n.n./n.n.n. hopping amplitude
  real(8) :: u,v            !local,non-local interaction
  real(8) :: tpd,vpd        !hybridization,band-band coupling
  real(8) :: ed0,ep0        !orbital energies
  real(8) :: xmu            !chemical potential
  real(8) :: beta      !inverse temperature
  real(8) :: dtau      
  real(8) :: fmesh
  real(8) :: eps            !broadening
  logical :: printf
  integer :: Nx,Ny
  real(8) :: wini,wfin
  integer :: Nsuccess
  real(8) :: weight
  real(8) :: heff
  logical :: chiflag
  real(8) :: cutoff !cutoff for the energy of the states for spectral summation
  real(8) :: eps_error

  !Dimension of the functions:
  !=========================================================
  integer :: NL,Ltau,Nw,Nfit

  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:) :: nmap,invnmap
  integer,allocatable,dimension(:,:) :: getloop
  integer,allocatable,dimension(:)   :: getCUPloop,getCDWloop
  integer,allocatable,dimension(:)   :: deg,getin,getis
  integer                            :: startloop,lastloop


  !Eigenvalues,Eigenvectors 
  !=========================================================
  type eigenspace
     real(8),dimension(:),pointer   :: e
     real(8),dimension(:,:),pointer :: M
  end type eigenspace
  type(eigenspace),dimension(:),allocatable :: espace


  !Partition function
  !=========================================================
  real(8) :: zeta_function


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm
  real(8),dimension(:),allocatable :: tau
  real(8),dimension(:),allocatable :: wr


  !Bath parameters (to be used in H)
  !=========================================================
  real(8),allocatable,dimension(:,:) :: ebath!epsiup,epsidw
  real(8),allocatable,dimension(:,:) :: vbath!vup,vdw


  !Functions for GETGFUNX 
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: Delta
  complex(8),allocatable,dimension(:,:) :: Giw,Siw
  complex(8),allocatable,dimension(:,:) :: Gwr,Swr
  real(8),allocatable,dimension(:)      :: Chitau

  !Variables for fixed density mu-loop 
  !=========================================================
  real(8) :: nread,nerr,ndelta,ndelta1,nindex,nindex1

  !Qties needed to get energy
  !=========================================================
  real(8) ::  nsimp,dimp,nupimp,ndwimp,magimp,m2imp

  !NML READ/WRITE UNITS
  character(len=32) :: Hfile,Ofile,GMimp1file,GRimp1file,CTimp1file

  !Writing/Reading Units:
  character(len=32) :: USEDinput

  namelist/EDvars/Ns,Nimp,Nspin,d,beta,xmu,u,v,tpd,tpp,ep0,ed0,ts,tsp,&
       nloop,chiflag,NL,Nw,Ltau,eps,weight,nread,nerr,ndelta,&
       wini,wfin,heff,Nx,Ny,Nfit,cutoff,eps_error,Nsuccess,&
       Hfile,Ofile,GMimp1file,GRimp1file,CTimp1file

contains

  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine read_input(INPUTunit)
    character(len=*) :: INPUTunit
    integer          :: i
    logical          :: control

    !ModelConf
    Ns         = 5
    Nimp       = 1
    Nspin      = 1
    !Hvars    
    d          = 1.d0
    ts         = 0.5d0
    tsp        = 0.d0
    beta       = 50.d0
    xmu        = 0.d0
    u          = 2.d0
    tpd        = 0.d0
    tpp        = 0.d0
    ed0        = 0.d0
    ep0        = 0.d0
    v          = 0.d0
    !Loops
    nloop      = 100
    ! iantif   =0
    chiflag    =.false.
    !parameters
    NL         = 1024
    Nw         = 1024
    Ltau       = 256
    Nfit       = 1024
    eps        = 0.01d0
    weight     = 0.8d0
    nread      = 0.d0
    nerr       = 1.d-4
    ndelta     = 0.1d0
    wini       =-4.d0
    wfin       = 4.d0
    heff       = 0.d0
    Nx         = 100
    Ny         = 100
    cutoff     =1.d-9
    eps_error  =1.d-5
    nsuccess   =2
    !ReadUnits
    Hfile      ="hamiltonian.restart"
    GMimp1file ="Gimp1_iw.ed"
    GRimp1file ="Gimp1_realw.ed"
    CTimp1file ="Chi_imp1_tau.ed"
    Ofile      ="observables.ed"

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
       write(50,*)" Ns   [5]    -- Number of bath sites per spin."
       write(50,*)" Norb [1]    -- Number of local orbitals"
       write(50,*)" Nspin [1]   -- Number of spin channels (max=2)"
       write(50,*)" nloop [50]  -- Max. number of iterations."
       write(50,*)" d [1.0]     -- Half-bandwidth."
       write(50,*)" beta [50.0] -- Inverse temperature."
       write(50,*)" xmu [0.0]   -- Chemical potential."
       write(50,*)" u [2.0]     -- Local interaction (Hubbard term)."
       write(50,*)" v [0.0]     -- N.N. interaction."
       write(50,*)" ts [0.5]    -- N.N. Hopping parameter."
       write(50,*)" tsp [0.0]   -- N.N.N. Hopping parameter."
       write(50,*)" tpd [0.0]   -- Local hybridization."
       write(50,*)" tpp [0.0]   -- N.N. diagonal hopping."
       write(50,*)" ep0 [0.0]   -- Impurity electronic level 1."
       write(50,*)" ed0 [0.0]   -- Impurity electronic level 2."
       write(50,*)" NL [2048]   -- Number of Matsubara frequencies."
       write(50,*)" Nw [1024]   -- Number of real frequencies."
       write(50,*)" Ltau [512]  -- Number of imaginary time points."
       write(50,*)" Nfit [2048] -- Number of fitted frequncies."
       write(50,*)" eps_error=[1.d-4] -- Convergence tolerance"
       write(50,*)" Nsuccess =[2]     -- Number of successive convergence treshold"
       write(50,*)" chiflag [.false.] -- Evaluation flag for spin susceptibility."
       write(50,*)" eps [0.035]       -- Broadening constant."
       write(50,*)" weigth [0.8]      -- Mixing weight (mispelled)."
       write(50,*)" nread [0.0]       -- Target density for chemical potential search."
       write(50,*)" nerror [1.d-4]    -- Tolerance in chemical potential search."
       write(50,*)" ndelta [0.10]     -- Starting delta for chemical potential search."
       write(50,*)" wmin [-4.0]       -- Lower bound frequency interval."
       write(50,*)" wmax [4.0]        -- Upper bound frequency interval."
       write(50,*)" heff [0.0]        -- Symmetry Breaking field."
       write(50,*)" Nx [100]          -- Number of points on a lattice mesh (x-axis)."
       write(50,*)" Ny [100]          -- Number of points on a lattice mesh (y-axis)."
       write(50,*)" cutoff [1.e-9]    -- Cutoff parameter for the spectrum contributing to GF calculation."
       write(50,*)" Hfile [Hamiltonian.restart] -- Store bath hamiltonian ."
       write(50,*)" Ofile [Observables.data]    -- Store observables."
       write(50,*)" Dfile [Delta.restart]       -- Store delta function."
       write(50,*)" GMfile [gfMimp.data]   -- Store GF Matsubara."
       write(50,*)" GRfile [gfRimp.data]   -- Store GF Real-axis."
       write(50,*)" CTfile [chiTimp.data]  -- Store Chi Im. time."
       stop
    endif

    call parse_cmd_variable(Ns,"NS")
    call parse_cmd_variable(Nimp,"NIMP")
    call parse_cmd_variable(Nspin,"NSPIN")
    call parse_cmd_variable(d,"D")
    call parse_cmd_variable(beta,"BETA")
    call parse_cmd_variable(xmu,"XMU")
    call parse_cmd_variable(u,"U")
    call parse_cmd_variable(v,"V")
    call parse_cmd_variable(ts,"TS")
    call parse_cmd_variable(tsp,"TSP")
    call parse_cmd_variable(tpp,"TPP")
    call parse_cmd_variable(tpd,"TPD")
    call parse_cmd_variable(ep0,"EP0")
    call parse_cmd_variable(ed0,"ED0")
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
    call parse_cmd_variable(eps,"EPS")
    call parse_cmd_variable(weight,"WEIGTH","WEIGHT")
    call parse_cmd_variable(cutoff,"CUTOFF")
    call parse_cmd_variable(heff,"HEFF")
    call parse_cmd_variable(Nx,"NX")
    call parse_cmd_variable(Ny,"NY")
    call parse_cmd_variable(Hfile,"HFILE")
    call parse_cmd_variable(Ofile,"OFILE")

    call version(revision)

    call allocate_system_structure()

    write(*,*)"CONTROL PARAMETERS"
    write(*,*)"--------------------------------------------"
    write(*,*)'| Total number of sites/spin   = ',Ns
    write(*,*)'| Maximum degeneracy           = ',NP
    write(*,*)'| Number of impurities         = ',Nimp
    write(*,*)'| Bath`s number of sites/spin  = ',Nbath
    write(*,*)'| Total size, Hilber space dim.= ',N,NN
    write(*,*)'| Number of sectors            = ',Nsect
    write(*,*)"--------------------------------------------"
    write(*,*)'| mu       = ',xmu
    write(*,*)'| T        = ',1.d0/beta
    write(*,*)"--------------------------------------------"
    write(*,nml=EDvars)
    write(*,*)"--------------------------------------------"
    print*,''
    USEDinput="used."//INPUTunit
    open(50,file=trim(adjustl(trim(USEDinput))))
    write(50,nml=EDvars)
    close(50)

    call global_check
    call global_allocation
  end subroutine read_input



  !+-------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine allocate_system_structure()
    select case(Ns)
    case (3)
       NP=9
    case (4)
       NP=36
    case (5)
       NP=100
    case (6)
       NP=400
    case (7)
       NP=1225
    case (8)
       NP=4900
    end select
    Nbath = Ns-Nimp
    N     = 2*Ns
    NN    = 2**N
    Nsect = ((Ns+1)*(Ns+1) - 1)
    allocate(nmap(Nsect,NP),invnmap(Nsect,NN))
    allocate(deg(Nsect),getin(Nsect),getis(Nsect))
    allocate(getloop(N,-N:N),getCUPloop(Nsect),getCDWloop(Nsect))

  end subroutine allocate_system_structure



  !+-------------------------------------------------------------------+
  !PURPOSE  : perform a sequence of check before calculation start
  ! if something wrong either EXIT with some error msg or fix on the fly 
  !+-------------------------------------------------------------------+
  subroutine global_check
    if(Nfit>NL)Nfit=NL
    if(Nimp>1)call abort("Norb > 1 is not yet supported!")
    if(Ns>8)call abort("Ns > 8 is too big!")
    if(nerr > eps_error) nerr=eps_error
  end subroutine global_check




  !+-------------------------------------------------------------------+
  !PURPOSE  : setup common arrays
  !+-------------------------------------------------------------------+
  subroutine global_allocation
    !Freq. arrays
    allocate(wm(NL),tau(0:Ltau))
    wm    = pi/beta*real(2*arange(1,NL)-1,8)
    tau   = linspace(0.d0,beta,Ltau+1,mesh=dtau)
    allocate(wr(Nw))
    wr    = linspace(wini,wfin,Nw,mesh=fmesh)

    !Bath:
    allocate(ebath(Nspin,Nbath),vbath(Nspin,Nbath))
    ! allocate(epsiup(Nbath))
    ! allocate(epsidw(Nbath))
    ! allocate(Vup(Nbath))
    ! allocate(Vdw(Nbath))

    !Functions
    allocate(Delta(Nspin,NL))
    allocate(Giw(Nspin,NL),Siw(Nspin,NL))
    allocate(Gwr(Nspin,Nw),Swr(Nspin,Nw))

    if(chiflag)allocate(Chitau(0:Ltau))
  end subroutine global_allocation


END MODULE ED_VARS_GLOBAL
