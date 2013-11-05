!########################################################################
!PROGRAM  : ED_AUX_FUNX
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_AUX_FUNX
  USE TIMER
  USE ED_VARS_GLOBAL
  implicit none
  private

  public :: init_ed_structure
  public :: setup_pointers
  public :: build_sector
  public :: bdecomp
  public :: c,cdg
  public :: binary_search

contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure
    integer          :: i,NP,nup,ndw
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
          write(*,*)"--------------------------------------------"
          write(LOGfile,"(A)")"Required Lanc_Nstates=1 => set T=0 calculation"
          write(*,*)"--------------------------------------------"
          write(*,*)""
       endif
    endif


    !check whether lanc_neigen and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_neigen,2)/=0)then
          lanc_neigen=lanc_neigen+1
          if(mpiID==0)then
             write(*,*)"--------------------------------------------"
             write(LOGfile,*)"Increased Lanc_Neigen:",lanc_neigen
             write(*,*)"--------------------------------------------"
             write(*,*)""
          endif
       endif
       if(mod(lanc_nstates,2)/=0)then
          lanc_nstates=lanc_nstates+1
          if(mpiID==0)then
             write(*,*)"--------------------------------------------"
             write(LOGfile,*)"Increased Lanc_Nstates:",lanc_nstates
             write(*,*)"--------------------------------------------"
             write(*,*)""
          endif
       endif
    endif

    !Some check:
    if(Nfit>NL)Nfit=NL
    if(Nspin>2)stop "Nspin > 2 ERROR. I guess you need to open the code at this point..."
    if(Norb>3)stop "Norb > 3 ERROR. I guess more than 3 bands is not feasing at this stage... " 
    if(nerr > eps_error) nerr=eps_error    

    !allocate functions
    allocate(impSmats(Nspin,Norb,Norb,NL))
    allocate(impSreal(Nspin,Norb,Norb,Nw))

    !allocate observables
    allocate(nimp(Norb),dimp(Norb),nupimp(Norb),ndwimp(Norb),magimp(Norb))
    allocate(m2imp(Norb,Norb))
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


END MODULE ED_AUX_FUNX
