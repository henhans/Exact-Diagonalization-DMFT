Module ED_VARS_GLOBAL
  !########################################################################
  !PROGRAM  : VARS_GLOBAL
  !TYPE     : Module
  !PURPOSE  : Declare the variables/arrays/matrices in common to DIAG_MODULE
  !AUTHORS  : Adriano Amaricci
  !########################################################################
  USE COMMON_VARS
  USE SLPLOT
  USE TOOLS
  USE CHRONOBAR
  implicit none

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


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:) :: nmap,invnmap
  integer,allocatable,dimension(:,:) :: getloop
  integer,allocatable,dimension(:)   :: getCUPloop,getCDWloop
  integer,allocatable,dimension(:)   :: deg,getin,getis
  integer                            :: startloop,lastloop
  integer                            :: startloop_mod,lastloop_mod


  !Dimension of the functions:
  !=========================================================
  integer :: NL,Ltau,Nw,Nfit
  integer :: Nx,Ny


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
  real(8),dimension(:),allocatable :: wm,tau
  real(8),dimension(:),allocatable :: wr


  !Bath parameters (to be used in H)
  !=========================================================
  !Spin order UP/DW
  real(8),allocatable,dimension(:) :: epsiup,epsidw
  real(8),allocatable,dimension(:) :: vup,vdw


  !Functions for GETGFUNX 
  !=========================================================
  !1imp
  complex(8),allocatable,dimension(:,:) :: Delta
  complex(8),allocatable,dimension(:,:) :: Giw,Siw
  complex(8),allocatable,dimension(:,:) :: Gwr,Swr
  real(8),allocatable,dimension(:)      :: Chitau
  !2imp
  complex(8),allocatable,dimension(:,:) :: G2iw
  complex(8),allocatable,dimension(:,:) :: G2wr
  real(8),allocatable,dimension(:)      :: Chi2tau,Chi12tau


  !Variables for fixed density mu-loop 
  !=========================================================
  real(8) :: nread,nerr,ndelta,ndelta1,nindex,nindex1

  !Global variables
  !=========================================================
  real(8) :: wini,wfin
  real(8) :: heff
  real(8) :: weigth
  real(8) :: pt
  logical :: chiflag
  real(8) :: cutoff !cutoff for the energy of the states for spectral summation
  integer :: Nsuccess
  real(8) :: eps_error

  !Qties needed to get energy
  !=========================================================
  real(8) ::  nimp1,dimp1,nupimp1,ndwimp1,magimp1,m2imp1
  real(8) ::  nimp2,dimp2,nupimp2,ndwimp2,magimp2,m2imp2
  real(8) ::  m2imp12


  !NML READ/WRITE UNITS
  character(len=32) :: Hfile,Ofile,&
       GMimp1file,GRimp1file,GMimp2file,GRimp2file,&
       CTimp1file,CTimp2file,CTimpAfile

  !Writing/Reading Units:
  character(len=32),parameter :: USEDinput

  namelist/EDvars/Ns,Nimp,Nspin,d,beta,xmu,u,v,tpd,tpp,ep0,ed0,ts,tsp,&
       nloop,chiflag,NL,Nw,Ltau,eps,weigth,nread,nerr,ndelta,&
       wini,wfin,heff,Nx,Ny,Nfit,cutoff,eps_error,Nsuccess,&
       Hfile,Ofile,GMimp1file,GRimp1file,GMimp2file,GRimp2file,&
       CTimp1file,CTimp2file,CTimpAfile



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

    include "help_buffer.f90"

    include "nml_default_values.f90"
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
       close(50)
       stop
    endif
    include "nml_read_cml.f90"

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

    temp=1.d0/beta  ; pt = pi/beta
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
    if(Nimp>2)call abort("Norb > 2 is not well supported yet!")
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
    allocate(epsiup(Nbath))
    allocate(epsidw(Nbath))
    allocate(Vup(Nbath))
    allocate(Vdw(Nbath))

    !Functions
    allocate(Delta(Nspin,NL))
    allocate(Giw(Nspin,NL),Siw(Nspin,NL))
    allocate(Gwr(Nspin,Nw),Swr(Nspin,Nw))

    if(chiflag)allocate(Chitau(0:Ltau))

    if(Nimp==2)then
       allocate(G2iw(Nspin,NL))
       allocate(G2wr(Nspin,Nw))
       if(chiflag)allocate(Chi2tau(0:Ltau),Chi12tau(0:Ltau))
    endif
  end subroutine global_allocation


END MODULE ED_VARS_GLOBAL
