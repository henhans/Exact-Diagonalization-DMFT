!########################################################################
!PROGRAM  : ED_AUX_FUNX
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_AUX_FUNX
  USE TIMER
  USE ED_VARS_GLOBAL
  implicit none
  private

  public :: ed_read_input
  public :: print_Hloc
#ifdef _MPI
  public :: ed_init_mpi
  public :: ed_finalize_mpi
#endif
  !
  public :: init_ed_structure
  public :: search_chemical_potential
  !
  public :: setup_pointers
  public :: build_sector
  public :: bdecomp
  public :: c,cdg
  public :: binary_search


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit,Hunit)
    character(len=*) :: INPUTunit
    character(len=*) :: Hunit
    logical          :: control
    integer          :: iorb,jorb,ispin,jspin
    if(mpiID==0)call version(revision)
    !DEFAULT VALUES OF THE PARAMETERS:
    !ModelConf
    Norb       = 1
    Nbath      = 4
    Nspin      = 1
    Uloc       = 0.d0;Uloc(1)=2.d0
    Ust        = 0.d0
    Jh         = 0.d0
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
    dmft_error  = 1.d-5
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
    LOGfile=6

    inquire(file=INPUTunit,exist=control)    
    if(control)then
       open(50,file=INPUTunit,status='old')
       read(50,nml=EDvars)
       close(50)
    else
       if(mpiID==0)then
          print*,"Can not find INPUT file"
          print*,"Printing a default version in default."//INPUTunit
          open(50,file="default."//INPUTunit)
          write(50,nml=EDvars)
          write(50,*)""
       endif
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
    call parse_cmd_variable(nloop,"NLOOP")
    call parse_cmd_variable(dmft_error,"DMFT_ERROR")
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
    call parse_cmd_variable(LOGfile,"LOGFILE")
    Ltau=max(int(beta),100)
    !
    if(mpiID==0)then
       open(50,file="used."//INPUTunit)
       write(50,nml=EDvars)
       close(50)
       write(*,*)"CONTROL PARAMETERS"
       write(*,nml=EDvars)
    endif


    write(LOGfile,"(A)")"U_local:"
    write(LOGfile,"(90F12.6,1x)")(Uloc(iorb),iorb=1,Norb)


    allocate(reHloc(Nspin,Nspin,Norb,Norb))
    allocate(imHloc(Nspin,Nspin,Norb,Norb))
    allocate(Hloc(Nspin,Nspin,Norb,Norb))
    reHloc = 0.d0
    imHloc = 0.d0

    inquire(file=Hunit,exist=control)    
    if(control)then
       open(50,file=Hunit,status='old')
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       close(50)
    else
       if(mpiID==0)then
          print*,"Can not find Uloc/Hloc file"
          print*,"Printing a default version in default."//Hunit
          open(50,file="default."//Hunit)
          do ispin=1,Nspin
             do iorb=1,Norb
                write(50,"(90F12.6)")((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
             enddo
          enddo
          write(50,*)""
          do ispin=1,Nspin
             do iorb=1,Norb
                write(50,"(90F12.6)")((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
             enddo
          enddo
          write(50,*)""
          close(50)
       endif
       stop
    endif

    hloc = dcmplx(reHloc,imHloc)

    write(LOGfile,"(A)")"H_local:"
    call print_Hloc(Hloc)
  end subroutine ed_read_input


  subroutine print_Hloc(hloc)
    integer                                     :: iorb,jorb,ispin,jspin
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(20(A1,F7.3,A1,F7.3,A1,2x))")&
               (&
               (&
               '(',dreal(Hloc(ispin,jspin,iorb,jorb)),',',dimag(Hloc(ispin,jspin,iorb,jorb)),')',&
               jorb =1,Norb),&
               jspin=1,Nspin)
       enddo
    enddo
  end subroutine print_Hloc


#ifdef _MPI
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_init_mpi
    call MPI_INIT(mpiERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
    write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
    call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  end subroutine ed_init_mpi

  subroutine ed_finalize_mpi
    call MPI_FINALIZE(mpiERR)
  end subroutine ed_finalize_mpi
#endif


  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure
    integer          :: i,NP,nup,ndw,iorb,jorb,ispin,jspin
    !Norb=# of impurity orbitals
    !Nbath=# of bath sites (per orbital or not depending on bath_type)
    !Ns=total number of sites
    !Nbo=total number of bath sites (all sites - impurity sites)
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
    end select
    Nbo   = Ns-Norb
    Ntot  = 2*Ns
    NN    = 2**Ntot
    Nsect = (Ns+1)*(Ns+1)
    !
    nup=Ns/2
    ndw=Ns-nup
    NP=(factorial(Ns)/factorial(nup)/factorial(Ns-nup))
    NP=NP*(factorial(Ns)/factorial(ndw)/factorial(Ns-ndw))
    if(mpiID==0)then
       write(*,*)"Summary:"
       write(*,*)"--------------------------------------------"
       write(*,*)'Number of impurities         = ',Norb
       write(*,*)'Number of bath/impurity      = ',Nbath
       write(*,*)'Total # of Bath sites/spin   = ',Nbo
       write(*,*)'Total # of sites/spin        = ',Ns
       write(*,*)'Maximum dimension            = ',NP
       write(*,*)'Total size, Hilber space dim.= ',Ntot,NN
       write(*,*)'Number of sectors            = ',Nsect
       write(*,*)"--------------------------------------------"
       print*,''
    endif

    allocate(impIndex(Norb,2))
    allocate(getdim(Nsect),getnup(Nsect),getndw(Nsect))
    allocate(getsector(0:Ns,0:Ns))
    allocate(getCsector(2,Nsect))
    allocate(getCDGsector(2,Nsect))
    allocate(getBathStride(Norb,Nbath))
    allocate(neigen_sector(Nsect))


    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates==1)then     !is you only want to keep 1 state
       lanc_neigen=1            !set the required eigen per sector to 1 see later for neigen_sector
       finiteT=.false.          !set to do zero temperature calculations
       if(mpiID==0)then
          write(LOGfile,"(A)")"Required Lanc_Nstates=1 => set T=0 calculation"
       endif
    endif


    !check whether lanc_neigen and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_neigen,2)/=0)then
          lanc_neigen=lanc_neigen+1
          if(mpiID==0)then
             write(LOGfile,"(A,I10)")"Increased Lanc_Neigen:",lanc_neigen
          endif
       endif
       if(mod(lanc_nstates,2)/=0)then
          lanc_nstates=lanc_nstates+1
          if(mpiID==0)then
             write(LOGfile,"(A,I10)")"Increased Lanc_Nstates:",lanc_nstates
          endif
       endif

    endif

    if(finiteT)then
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
    endif

    !Some check:
    if(Nfit>NL)Nfit=NL
    if(Nspin>2)stop "Nspin > 2 ERROR. ask developer or develop your own on separate branch"
    if(Norb>3)stop "Norb > 3 ERROR. ask developer or develop your own on separate branch" 
    if(nerr < dmft_error) nerr=dmft_error
    if(ed_method=='full'.AND.bath_type=='hybrid')stop "FULL ED & HYBRID not implemented yet:ask developer..."
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop
       nloop=(i-1)*niter                !increase the max number of dmft loop allowed so to do threshold loop
       write(LOGfile,"(A,I10)")"Increased Nloop to:",nloop
    endif

    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,NL))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Nw))

    !allocate observables
    allocate(nimp(Norb),dimp(Norb))
  end subroutine init_ed_structure




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: i,in,dim,isector,jsector,dimup,dimdw
    integer                          :: nup,ndw,jup,jdw,iorb
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    if(mpiID==0)write(LOGfile,"(A)")"Setting up pointers:"
    call start_timer
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getsector(nup,ndw)=isector
          getnup(isector)=nup
          getndw(isector)=ndw
          dimup=(factorial(Ns)/factorial(nup)/factorial(Ns-nup))
          dimdw=(factorial(Ns)/factorial(ndw)/factorial(Ns-ndw))
          dim=dimup*dimdw
          getdim(isector)=dim
          neigen_sector(isector) = min(dim,lanc_neigen)   !init every sector to required eigenstates
       enddo
    enddo
    call stop_timer

    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo

    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)      = Norb + i
       enddo
    end select

    getCsector=0
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo

    getCDGsector=0
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsect
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers



  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map from the 
  !states i\in Hilbert_space to the states count in H_sector.
  !+------------------------------------------------------------------+
  !|ImpUP,BathUP>|ImpDW,BathDW >
  subroutine build_sector(isector,map)
    integer              :: i,j,isector,iup,idw,count
    integer              :: nup,ndw
    integer              :: ivec(Ntot)
    integer,dimension(:) :: map
    nup = getnup(isector)
    ndw = getndw(isector)
    count=0
    do i=1,NN
       call bdecomp(i,ivec)
       iup = sum(ivec(1:Ns))
       idw = sum(ivec(Ns+1:2*Ns))
       if(iup==nup.AND.idw==ndw)then
          count             = count+1 !count the states in the sector (n_up,n_dw)
          map(count)        = i       !build the map to full space states
       endif
    enddo
  end subroutine build_sector


  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Ntot)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  subroutine bdecomp(i,ivec)
    integer :: ivec(Ntot)         
    integer :: l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end subroutine bdecomp




  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(m,i,j,sgn)
    integer :: ib(Ntot)
    integer :: i,j,m,km,k
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==0)then
       j=0
    else
       if(m==1)then
          j=i-1
       else
          km=0
          do k=1,m-1
             km=km+ib(k)
          enddo
          !km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i-2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine c



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j,sgn)
    integer :: ib(Ntot)
    integer :: i,j,m,km,k
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==1)then
       j=0
    else
       if(m==1)then
          j=i+1
       else
          km=0
          do k=1,m-1
             km=km+ib(k)
          enddo
          !km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i+2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine cdg




  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
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




  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine search_chemical_potential(ntmp,niter,converged)
    real(8),intent(in)    :: ntmp
    integer,intent(in)    :: niter
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0
    integer,save          :: nindex=0
    integer               :: nindex1
    real(8)               :: ndelta1,nratio
    integer,save          :: nth_magnitude=-1,nth_magnitude_old=-1
    real(8),save          :: nth=1.d-1
    logical,save          :: ireduce=.true.
    integer :: unit
    !
    ndiff=ntmp-nread
    !nratio = 0.5d0
    nratio = 1.d0/(6.d0/11.d0*pi)
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    nindex1=nindex
    ndelta1=ndelta
    if(ndiff >= nth)then      !if((ntmp >= nread+nth))then
       nindex=-1
    elseif(ndiff <= -nth)then !elseif(ntmp <= nread-nth)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       ndelta=ndelta1*nratio !decreasing the step
    else
       ndelta=ndelta1
    endif
    !
    !update chemical potential
    xmu=xmu+dble(nindex)*ndelta
    !
    !Print information
    write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
    if(nindex>0)then
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," ==>"
    elseif(nindex<0)then
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," <=="
    else
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," == "
    endif
    write(LOGfile,"(A,f15.9)")"xmu  = ",xmu
    write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
    unit=free_unit()
    open(unit,file="search_mu_iteration.ed",position="append")
    write(unit,*)xmu,ntmp,ndiff
    close(unit)
    !
    !check convergence within actual threshold
    !if reduce is activetd
    !if density is in the actual threshold
    !if DMFT is converged
    !if threshold is larger than nerror (i.e. this is not last loop)
    bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
    if(bool)then
       nth_magnitude_old=nth_magnitude        !save old threshold magnitude
       nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
       nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
       count=0                                !reset the counter
       converged=.false.                      !reset convergence
       !experimental
       ndelta=ndelta1*nratio
       !
    endif
    !
    !if density is not converged set convergence to .false.
    if(abs(ntmp-nread)>nth)converged=.false.
    !
    !check convergence for the smallest threshold
    !if smallest threshold
    !if reduce is active
    !if # iterations > max number
    !if not yet converged
    !set threshold back to the previous larger one.
    bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
    if(bool)then
       ireduce=.false.
       nth=10.d0**(nth_magnitude_old)
    endif
    !
    write(LOGfile,"(A,I5)")"count= ",count
    write(LOGfile,"(A,L2)"),"Converged=",converged
    print*,""
    !
  end subroutine search_chemical_potential



  ! subroutine search_mu(ntmp,convergence)
  !   logical,intent(inout) :: convergence
  !   real(8)               :: ntmp
  !   logical               :: check
  !   integer,save          :: count=0
  !   integer,save          :: nindex=0
  !   real(8)               :: ndelta1,nindex1
  !   if(count==0)then
  !      inquire(file="searchmu_file.restart",exist=check)
  !      if(check)then
  !         open(10,file="searchmu_file.restart")
  !         read(10,*)ndelta,nindex
  !         close(10)
  !      endif
  !   endif
  !   count=count+1
  !   nindex1=nindex
  !   ndelta1=ndelta
  !   if((ntmp >= nread+nerr))then
  !      nindex=-1
  !   elseif(ntmp <= nread-nerr)then
  !      nindex=1
  !   else
  !      nindex=0
  !   endif
  !   if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
  !      ndelta=ndelta1/2.d0 !decreasing the step       
  !   else
  !      ndelta=ndelta1
  !   endif
  !   xmu=xmu+real(nindex,8)*ndelta
  !   if(abs(ntmp-nread)>nerr)convergence=.false.
  !   write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",ntmp," /",nread,&
  !        "| shift=",nindex*ndelta,"| xmu=",xmu
  !   write(*,"(A,f15.12)")"dn=",abs(ntmp-nread)
  !   print*,""
  !   print*,"Convergence:",convergence
  !   print*,""
  !   open(10,file="searchmu_file.restart.new")
  !   write(10,*)ndelta,nindex,xmu
  !   close(10)
  ! end subroutine search_mu



END MODULE ED_AUX_FUNX
