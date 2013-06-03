MODULE ED_VARS_GLOBAL
  USE COMMON_VARS
  USE SCIFOR_VERSION
  USE TIMER
  USE PARSE_CMD
  USE IOTOOLS
  USE FUNCTIONS
  USE MATRIX
  USE OPTIMIZE
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
  integer :: nGFitermax     !Max number of iteration in resolvant tri-diagonalization
  integer :: cgNitmax       !Max number of iteration in the fit
  real(8) :: cgFtol         !Tolerance in the cg fit
  integer :: cgType         !CGfit mode 0=normal,1=1/n weight, 2=1/w weight

  !Dimension of the functions:
  !=========================================================
  integer :: NL,Ltau,Nw,Nfit


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:) :: nmap,invnmap
  integer,allocatable,dimension(:,:) :: getloop
  integer,allocatable,dimension(:)   :: getCUPloop,getCDWloop
  integer,allocatable,dimension(:)   :: getCDGUPloop,getCDGDWloop
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


  !Functions for GETGFUNX (names must be changed)
  !=========================================================
  complex(8),allocatable,dimension(:,:) :: Giw,Siw
  complex(8),allocatable,dimension(:,:) :: Gwr,Swr
  real(8),allocatable,dimension(:)      :: Chitau
  complex(8),allocatable,dimension(:)   :: Chiw
  !2imp
  complex(8),allocatable,dimension(:,:) :: G2iw
  complex(8),allocatable,dimension(:,:) :: G2wr
  real(8),allocatable,dimension(:)      :: Chi2tau,Chi12tau
  complex(8),allocatable,dimension(:)   :: Chi2w,Chi12w


  !Variables for fixed density mu-loop 
  !=========================================================
  real(8) :: nread,nerr,ndelta

  !Qties needed to get energy
  !=========================================================
  real(8) ::  nsimp,dimp,nupimp,ndwimp,magimp,m2imp
  real(8) ::  nimp2,dimp2,nupimp2,ndwimp2,magimp2,m2imp2
  real(8) ::  m2imp12


  !NML READ/WRITE UNITS
  character(len=32) :: Hfile,Ofile,&
       GMfile,GRfile,  &
       GM2file,GR2file,&
       CTfile,CT2file,CT12file, &
       CWfile,CW2file,CW12file
  integer           :: LOGfile


  !Writing/Reading Units:
  character(len=32) :: USEDinput


  namelist/EDvars/Ns,Nimp,Nspin,&       
       beta,xmu,nloop,u,tpd,ep0,&
       eps,wini,wfin,heff,      &
       NL,Nw,Ltau,Nfit,         &
       nread,nerr,ndelta,       &
       chiflag,cutoff,HFmode,   &
       eps_error,Nsuccess,      &
       cgNitmax,cgFtol,cgType,  &
       Hfile,Ofile,LOGfile,     &
       GMfile,GRfile,           &  
       GM2file,GR2file,         &
       CTfile,CT2file,CT12file,  &
       CWfile,CW2file,CW12file

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
    nGFitermax = 100
    cgNitmax   = 1000
    cgFtol     = 1.d-9
    cgType     = 0
    !ReadUnits
    Hfile  ="hamiltonian.restart"
    GMfile="Gimp1_iw.ed"
    GRfile="Gimp1_realw.ed"
    GM2file="Gimp2_iw.ed"
    GR2file="Gimp2_realw.ed"
    CTfile="Chi1_tau.ed"
    CWfile="Chi1_realw.ed"
    CT2file="Chi2_tau.ed"
    CW2file="Chi2_realw.ed"
    CT12file="Chi12_tau.ed"  
    CW12file="Chi12_realw.ed"  
    Ofile  ="observables.ed"
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
    call parse_cmd_variable(Nimp,"NIMP")
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
    call parse_cmd_variable(cgNitmax,"CGNITMAX")
    call parse_cmd_variable(cgFtol,"CGFTOL")
    call parse_cmd_variable(cgType,"CGTYPE")
    call parse_cmd_variable(Hfile,"HFILE")
    call parse_cmd_variable(Ofile,"OFILE")
    call parse_cmd_variable(GMfile,"GMFILE")
    call parse_cmd_variable(GRfile,"GRFILE")
    call parse_cmd_variable(CTfile,"CTFILE")
    call parse_cmd_variable(CWfile,"CWFILE")
    call parse_cmd_variable(GM2file,"GM2FILE")
    call parse_cmd_variable(GR2file,"GR2FILE")
    call parse_cmd_variable(CT2file,"CT2FILE")
    call parse_cmd_variable(CW2file,"CW2FILE")
    call parse_cmd_variable(CT12file,"CT12FILE")
    call parse_cmd_variable(CW12file,"CW12FILE")
    call parse_cmd_variable(LOGfile,"LOGFILE")

    call version(revision)

    call allocate_system_structure()

    write(*,*)"CONTROL PARAMETERS"
    write(*,nml=EDvars)
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
    print*,''
    USEDinput="used."//INPUTunit
    open(50,file=trim(adjustl(trim(USEDinput))))
    write(50,nml=EDvars)
    close(50)

    !Some check:
    if(Nfit>NL)Nfit=NL
    if(Nimp>2)call abort("Norb > 2 is not well supported yet!")
    if(Ns>8)call abort("Ns > 8 is too big!")
    if(nerr > eps_error) nerr=eps_error    

    !allocate functions
    allocate(Giw(Nspin,NL),Siw(Nspin,NL))
    allocate(Gwr(Nspin,Nw),Swr(Nspin,Nw))
    if(chiflag)allocate(Chitau(0:Ltau),Chiw(Nw))
    if(Nimp==2)then
       allocate(G2iw(Nspin,NL))
       allocate(G2wr(Nspin,Nw))
       if(chiflag)allocate(Chi2tau(0:Ltau),Chi12tau(0:Ltau),Chi2w(Nw),Chi12w(Nw))
    endif
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
    allocate(getloop(N,-N:N))
    allocate(getCUPloop(Nsect),getCDWloop(Nsect))
    allocate(getCDGUPloop(Nsect),getCDGDWloop(Nsect))
  end subroutine allocate_system_structure



END MODULE ED_VARS_GLOBAL
