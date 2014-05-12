!########################################################################
!PROGRAM  : ED_AUX_FUNX
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_AUX_FUNX
  USE MPI_VARS
  USE TIMER
  USE IOTOOLS, only:free_unit,reg
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  implicit none
  private

  public :: print_Hloc
  !
  public :: init_ed_structure
  public :: search_chemical_potential
  !
  public :: setup_pointers
  public :: setup_pointers_sc
  public :: build_sector
  public :: bdecomp
  public :: c,cdg
  public :: binary_search

contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure(Hunit)
    character(len=64)                        :: Hunit
    logical                                  :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: imHloc         !local hamiltonian, imag part
    integer                                  :: i,NP,nup,ndw,iorb,jorb,ispin,jspin
    !
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
    !
    nup=Ns/2
    ndw=Ns-nup
    if(.not.ed_supercond)then
       Nsect = (Ns+1)*(Ns+1)
       NP=get_sector_dimension(nup,ndw)
    else
       Nsect = Ntot+1
       NP=get_sc_sector_dimension(0)
    endif
    !
    if(ED_MPI_ID==0)then
       write(LOGfile,*)"Summary:"
       write(LOGfile,*)"--------------------------------------------"
       write(LOGfile,*)'Number of impurities         = ',Norb
       write(LOGfile,*)'Number of bath/impurity      = ',Nbath
       write(LOGfile,*)'Total # of Bath sites/spin   = ',Nbo
       write(LOGfile,*)'Total # of sites/spin        = ',Ns
       write(LOGfile,*)'Maximum dimension            = ',NP
       write(LOGfile,*)'Total size, Hilber space dim.= ',Ntot,NN
       write(LOGfile,*)'Number of sectors            = ',Nsect
       write(LOGfile,*)"--------------------------------------------"
    endif

    allocate(Hloc(Nspin,Nspin,Norb,Norb))
    reHloc = 0.d0
    imHloc = 0.d0

    inquire(file=Hunit,exist=control)
    if(control)then
       if(ED_MPI_ID==0)write(LOGfile,*)"Reading Hloc from file: "//Hunit
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
       if(ED_MPI_ID==0)then
          write(LOGfile,*)"Hloc file not found."
          write(LOGfile,*)"Hloc should be defined elsewhere..."
       endif
    endif
    Hloc = dcmplx(reHloc,imHloc)
    if(ED_MPI_ID==0)then
       write(LOGfile,"(A)")"H_local:"
       call print_Hloc(Hloc)
    endif



    allocate(impIndex(Norb,2))
    allocate(getdim(Nsect),getnup(Nsect),getndw(Nsect),getsz(Nsect))
    if(.not.ed_supercond)then
       allocate(getsector(0:Ns,0:Ns))
    else
       allocate(getsector(-Ns:Ns,1))
    endif
    allocate(getCsector(2,Nsect))
    allocate(getCDGsector(2,Nsect))
    allocate(getBathStride(Norb,Nbath))
    allocate(neigen_sector(Nsect))


    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       lanc_nstates_sector=1            !set the required eigen per sector to 1 see later for neigen_sector
       finiteT=.false.          !set to do zero temperature calculations
       if(ED_MPI_ID==0)then
          write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
       endif
    endif


    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          if(ED_MPI_ID==0)write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          if(ED_MPI_ID==0)write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif

    endif

    if(finiteT)then
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
    else
       if(ED_MPI_ID==0)write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
    endif

    !Some check:
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "Nspin > 2 ERROR. ask developer or develop your own on separate branch"
    if(Norb>3)stop "Norb > 3 ERROR. ask developer or develop your own on separate branch" 
    if(nerr < dmft_error) nerr=dmft_error
    if(ed_method=='full'.AND.bath_type=='hybrid')stop "FULL ED & HYBRID not implemented yet:ask developer..."
    if(ed_supercond)then
       if(Nspin>1)stop "SC+AFM ERROR. ask developer or develop your own on separate branch" 
       if(ed_method=='full')stop "FULL ED & SUPERC is not implemented yet:ask developer..."
       if(Norb>1)stop "SC Multi-Band not yet implemented. Wait for the developer to understand what to do..."
       if(ed_type=='c')stop "SC with Hermitian H not yet implemented. Wait for the developer to code it..."
    endif
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
       !nloop=(i-1)*niter                !increase the max number of dmft loop allowed so to do threshold loop
       !write(LOGfile,"(A,I10)")"Increased Nloop to:",nloop
    endif

    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    if(ed_supercond)then
       allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats))
       allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal))
    endif

    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb))
    if(ed_supercond)allocate(ed_phisc(Norb))
  end subroutine init_ed_structure





  subroutine print_Hloc(hloc,unit)
    integer,optional                            :: unit
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
    if(present(unit))then
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dreal(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,"(90F12.6)")((dimag(Hloc(ispin,jspin,iorb,jorb)),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       write(unit,*)""
    endif
  end subroutine print_Hloc






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: i,in,dim,isector,jsector,dimup,dimdw
    integer                          :: nup,ndw,jup,jdw,iorb
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getsector(nup,ndw)=isector
          getnup(isector)=nup
          getndw(isector)=ndw
          dim = get_sector_dimension(nup,ndw)
          getdim(isector)=dim
          neigen_sector(isector) = min(dim,lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    enddo
    if(ED_MPI_ID==0)call stop_timer

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


  subroutine setup_pointers_sc
    integer                          :: i,isz,in,dim,isector,jsector
    integer                          :: sz,iorb,dim2,jsz
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    if(ED_MPI_ID==0)write(LOGfile,"(A)")"Setting up pointers:"
    if(ED_MPI_ID==0)call start_timer
    isector=0
    do isz=-Ns,Ns
       sz=abs(isz)
       isector=isector+1
       getsector(isz,1)=isector
       getsz(isector)=isz
       dim = get_sc_sector_dimension(isz)
       getdim(isector)=dim
       neigen_sector(isector) = min(dim,lanc_nstates_sector)   !init every sector to required eigenstates
       ! !<DEBUG
       ! allocate(imap(dim))
       ! call build_sector(isector,imap,dim2)
       ! print*,isz,dim,dim2
       ! deallocate(imap)
       ! !>DEBUG
    enddo
    if(ED_MPI_ID==0)call stop_timer

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
    !c_up
    do isector=1,Nsect
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCsector(1,isector)=jsector
    enddo
    !c_dw
    do isector=1,Nsect
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCsector(2,isector)=jsector
    enddo

    getCDGsector=0
    !cdg_up
    do isector=1,Nsect
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCDGsector(1,isector)=jsector
    enddo
    !cdg_dw
    do isector=1,Nsect
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_sc

















  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !+------------------------------------------------------------------+
  !|ImpUP,BathUP>|ImpDW,BathDW >
  subroutine build_sector(isector,map,dim2)
    integer              :: i,j,isector,iup,idw,mz,dim
    integer,optional     :: dim2
    integer              :: nup,ndw,sz
    integer              :: ivec(Ntot)
    integer,dimension(:) :: map
    !if(size(map)/=getdim(isector)stop "error in build_sector: wrong dimension of map"
    dim=0
    if(.not.ed_supercond)then
       nup = getnup(isector)
       ndw = getndw(isector)
       do i=1,NN
          call bdecomp(i,ivec)
          iup = sum(ivec(1:Ns))
          idw = sum(ivec(Ns+1:2*Ns))
          if(iup==nup.AND.idw==ndw)then
             dim           = dim+1 !count the states in the sector (n_up,n_dw)
             map(dim)      = i       !build the map to full space states
          endif
       enddo
    else
       sz = getsz(isector)
       do i=1,NN
          call bdecomp(i,ivec)
          mz = sum(ivec(1:Ns)) - sum(ivec(Ns+1:2*Ns))
          if(mz==sz)then
             dim             = dim+1 !count the states in the sector (n_up,n_dw)
             map(dim)        = i       !build the map to full space states
          endif
       enddo
    endif
    if(present(dim2))dim2=dim
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
  !PURPOSE  : calculate the factorial
  !+------------------------------------------------------------------+
  function get_sector_dimension(nup,ndw) result(dim)
    integer :: nup,ndw,dim,dimup,dimdw
    dimup=(factorial(Ns)/factorial(nup)/factorial(Ns-nup))
    dimdw=(factorial(Ns)/factorial(ndw)/factorial(Ns-ndw))
    dim=dimup*dimdw
  end function get_sector_dimension




  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial
  !+------------------------------------------------------------------+
  function get_sc_sector_dimension(mz) result(dim)
    integer :: mz
    integer :: i,dim,Nb
    dim=0
    Nb=Ns-mz
    do i=0,Nb/2 
       dim=dim + 2**(Nb-2*i)*nchoos(ns,Nb-2*i)*nchoos(ns-Nb+2*i,i)
    enddo
  end function get_sc_sector_dimension


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
  !PURPOSE  : calculate the binomial factor
  !+------------------------------------------------------------------+
  function nchoos(n1,n2)
    real(8) :: xh
    integer :: n1,n2,i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*real(n1+1-i,8)/real(i,8)
    enddo
    nchoos = int(xh + 0.5d0)
  end function nchoos




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
  subroutine search_chemical_potential(ntmp,converged)
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2
    logical,save          :: ireduce=.true.
    integer               :: unit
    !
    if(ED_MPI_ID==0)then
       ndiff=ntmp-nread
       nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
       !
       !check actual value of the density *ntmp* with respect to goal value *nread*
       count=count+1
       totcount=totcount+1
       if(count>2)then
          do i=1,2
             nindex_old(i+1)=nindex_old(i)
          enddo
       endif
       nindex_old(1)=nindex
       !
       if(ndiff >= nth)then
          nindex=-1
       elseif(ndiff <= -nth)then
          nindex=1
       else
          nindex=0
       endif
       !
       ndelta_old=ndelta
       bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
       !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       if(bool)then
          ndelta=ndelta_old*nratio !decreasing the step
       else
          ndelta=ndelta_old
       endif
       !
       if(ndelta_old<1.d-9)then
          ndelta_old=0.d0
          nindex=0
       endif
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
       open(unit,file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
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
          ndelta=ndelta_old*nratio                  !reduce the delta step
          !
       endif
       !
       !if density is not converged set convergence to .false.
       if(abs(ntmp-nread)>nth)converged=.false.
       !
       !check convergence for this threshold
       !!---if smallest threshold-- NO MORE
       !if reduce is active (you reduced the treshold at least once)
       !if # iterations > max number
       !if not yet converged
       !set threshold back to the previous larger one.
       !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
       bool=ireduce.AND.(count>niter).AND.(.not.converged)
       if(bool)then
          ireduce=.false.
          nth=10.d0**(nth_magnitude_old)
       endif
       !
       write(LOGfile,"(A,I5)")"count= ",count
       write(LOGfile,"(A,L2)"),"Converged=",converged
       print*,""
       !
    endif
#ifdef _MPI
    call MPI_BCAST(xmu,1,MPI_Double_Precision,0,MPI_COMM_WORLD,mpiERR)
#endif
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
