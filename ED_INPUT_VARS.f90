MODULE ED_INPUT_VARS
  USE SCIFOR_VERSION
  USE PARSE_INPUT
  implicit none

  !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile


  !input variables
  !=========================================================
  integer                                     :: Nbath               !Nbath=# of bath sites (per orbital or not depending on bath_type)
  integer                                     :: Norb                !Norb =# of impurity orbitals
  integer                                     :: Nspin               !Nspin=# spin degeneracy (max 2)
  integer                                     :: nloop               !max dmft loop variables
  real(8)                                     :: Ust,Jh              !intra-orbitals interactions
  real(8),dimension(3)                        :: Uloc                !local interactions
  real(8)                                     :: xmu                 !chemical potential
  real(8)                                     :: deltasc             !breaking symmetry field
  real(8)                                     :: beta                !inverse temperature
  real(8)                                     :: eps                 !broadening
  real(8)                                     :: wini,wfin           !
  integer                                     :: Nsuccess            !
  logical                                     :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                     :: chiflag             !
  logical                                     :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                                     :: cutoff              !cutoff for spectral summation
  real(8)                                     :: dmft_error          !dmft convergence threshold
  real(8)                                     :: sb_field            !symmetry breaking field
  integer                                     :: lanc_niter          !Max number of Lanczos iterations
  integer                                     :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer                                     :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer                                     :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer                                     :: cg_Niter            !Max number of iteration in the fit
  real(8)                                     :: cg_Ftol             !Tolerance in the cg fit
  integer                                     :: cg_Weight           !CGfit mode 0=normal,1=1/n weight, 2=1/w weight
  character(len=5)                            :: cg_Scheme           !fit scheme: delta (default), weiss for G0^
  logical                                     :: finiteT             !flag for finite temperature calculation
  character(len=4)                            :: ed_method           !flag to set ed method solution: lanc=lanczos method, full=full diagonalization
  character(len=1)                            :: ed_type             !flag to set real or complex Ham: d=symmetric H (real), c=hermitian H (cmplx)
  logical                                     :: ed_supercond        !flag to set ed symmetry type: F=normal (default), T=superc=superconductive
  character(len=7)                            :: bath_type           !flag to set bath type: irreducible (1bath/imp), reducible(1bath)
  character(len=100)                          :: ed_file_suffix      !suffix string attached to the output files.
  real(8)                                     :: nread               !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)                                     :: nerr                !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)                                     :: ndelta              !initial chemical potential step
  integer                                     :: niter
  !logical                                     :: ed_verbose !verbose flag:
  integer                                     :: ed_verbose
  !Some parameters for function dimension:
  !=========================================================
  integer                                     :: Lmats
  integer                                     :: Lreal
  integer                                     :: Lfit
  integer                                     :: Ltau

  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)                          :: Hfile
  integer                                     :: LOGfile


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit)
    character(len=*) :: INPUTunit
    logical          :: control
    integer          :: iorb,jorb,ispin,jspin
    call version(revision)
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1)
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=4)
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1)
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2.d0,0.d0,0.d0])
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0)
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0)
    call parse_input_variable(beta,"BETA",INPUTunit,default=500.d0)
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0)
    call parse_input_variable(deltasc,"DELTASC",INPUTunit,default=0.02d0)
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100)
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0)
    call parse_input_variable(sb_field,"SB_FIELD",INPUTunit,default=0.1d0)
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1)
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=2000)
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=2000)
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000)
    call parse_input_variable(Lfit,"LFIT",INPUTunit,default=1000)
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0)
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4)
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0)
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0)
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0)
    call parse_input_variable(chiflag,"CHIFLAG",INPUTunit,default=.false.)
    call parse_input_variable(jhflag,"JHFLAG",INPUTunit,default=.false.)
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.)
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0)
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9)
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=1)
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1)
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512)
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200)
    call parse_input_variable(cg_niter,"CG_NITER",INPUTunit,default=200)
    call parse_input_variable(cg_scheme,"CG_SCHEME",INPUTunit,default='delta')
    call parse_input_variable(cg_ftol,"CG_FTOL",INPUTunit,default=0.00001d0)
    call parse_input_variable(cg_weight,"CG_WEIGHT",INPUTunit,default=0)
    call parse_input_variable(ed_Type,"ED_TYPE",INPUTunit,default='d')
    call parse_input_variable(ed_Supercond,"ED_SUPERCOND",INPUTunit,default=.false.)
    call parse_input_variable(ed_Method,"ED_METHOD",INPUTunit,default='lanc')
    call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal')
    call parse_input_variable(Hfile,"HFILE",INPUTunit,default="hamiltonian")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6)
    call parse_input_variable(ed_verbose,"ED_VERBOSE",INPUTunit,default=0)
    call parse_input_variable(ed_file_suffix,"ED_FILE_SUFFIX",INPUTunit,default=".ed")
    call substring_delete(ed_file_suffix,".ed")
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".ed")
    call check_input_file(INPUTunit)
    Ltau=max(int(beta),Ltau)
  end subroutine ed_read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE ED_INPUT_VARS
