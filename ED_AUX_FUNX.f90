!########################################################################
!PROGRAM  : ED_AUX_FUNX
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_AUX_FUNX
  USE ED_VARS_GLOBAL
  implicit none
  private

  public :: init_ed_structure
  public :: setup_pointers
  public :: build_sector
  public :: bdecomp
  public :: c,cdg
  public :: search_mu
  public :: binary_search

contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Init calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure
    integer          :: i,NP,nup,ndw
    !!<MPI
    ! if(mpiID==0)then
    !!>MPI
    !Norb=# of impurity orbitals
    !Nbath=# of bath sites per orbital
    !Ns=total number of sites
    !Nbo=total number of bath sites (all sites - impurity sites)
    Ns=(Nbath+1)*Norb
    Nbo   = Ns-Norb
    Ntot  = 2*Ns
    NN    = 2**Ntot
    Nsect = (Ns+1)*(Ns+1)
    !
    nup=Ns/2
    ndw=Ns-nup
    NP=(factorial(Ns)/factorial(nup)/factorial(Ns-nup))
    NP=NP*(factorial(Ns)/factorial(ndw)/factorial(Ns-ndw))
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
    !!<MPI
    ! endif
    !!>MPI

    allocate(impIndex(Norb,2))
    !allocate(Hmap(Nsect),invHmap(Nsect,NN))
    allocate(getdim(Nsect),getnup(Nsect),getndw(Nsect))
    allocate(getsector(0:Ns,0:Ns))
    allocate(getCsector(2,Nsect))
    allocate(getCDGsector(2,Nsect))
    allocate(neigen_sector(Nsect))

    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates==1)then     !is you only want to keep 1 state
       lanc_neigen=1            !set the required eigen per sector to 1 see later for neigen_sector
       finiteT=.false.          !set to do zero temperature calculations
       write(*,*)"--------------------------------------------"
       write(LOGfile,"(A)")"Required Lanc_Nstates=1 => set T=0 calculation"
       write(*,*)"--------------------------------------------"
       write(*,*)""
    endif

    !check whether lanc_neigen and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_neigen,2)/=0)then
          lanc_neigen=lanc_neigen+1
          write(*,*)"--------------------------------------------"
          write(LOGfile,*)"Increased Lanc_Neigen:",lanc_neigen
          write(*,*)"--------------------------------------------"
          write(*,*)""
       endif
       if(mod(lanc_nstates,2)/=0)then
          lanc_nstates=lanc_nstates+1
          write(*,*)"--------------------------------------------"
          write(LOGfile,*)"Increased Lanc_Nstates:",lanc_nstates
          write(*,*)"--------------------------------------------"
          write(*,*)""
       endif
    endif

    !Some check:
    if(Nfit>NL)Nfit=NL
    if(Norb>3)stop "Norb > 3 ERROR. I guess you need to open the code at this point..." 
    if(nerr > eps_error) nerr=eps_error    

    !allocate functions
    allocate(impGmats(Nspin,Norb,NL),impSmats(Nspin,Norb,NL))
    allocate(impGreal(Nspin,Norb,Nw),impSreal(Nspin,Norb,Nw))

    !allocate observables
    allocate(nimp(Norb),dimp(Norb),nupimp(Norb),ndwimp(Norb),magimp(Norb))
    allocate(m2imp(Norb,Norb))
  end subroutine init_ed_structure


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: i,in,dim,isector,jsector,dimup,dimdw
    integer                          :: nup,ndw,jup,jdw
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    ! !<DEBUG
    integer,dimension(:),allocatable :: Ninv,Nstock
    ! integer :: ivec(Ntot),kcp,is
    ! !>DEBUG

    write(LOGfile,"(A)")"Setting up pointers:"
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



    !!>DEBUG
    ! !check that invHmap is always filled with zeros, what a weist!!
    ! do isector=1,Nsect
    !    print*,isector
    !    do i=1,NN
    !       print*,invHmap(isector,i)
    !    enddo
    !    print*,""
    ! end do
    ! nup=Ns/2;ndw=nup
    ! isector=getsector(nup,ndw) !half-filling
    ! dim    =getdim(isector)     !Nstates(nup,ndw)
    ! allocate(Nstock(dim))
    ! allocate(Ninv(NN))
    ! call build_sector(isector,Nstock,Ninv)
    ! ! Ninv=0
    ! ! do i=1,dim
    ! !    Ninv(Hmap(isector)%map(i))=i
    ! ! enddo
    ! print*,size(Ninv)!,size(invHmap(isector,:))
    ! do i=1,NN
    !    print*,i,Ninv(i)!,invHmap(isector,i)
    ! enddo
    ! ! for Nbath=3, state 90 in HS is related to state 10 in Sector(2,2)
    ! print*,binary_search(Nstock,90)
    ! !THIS WAS TO TEST THE F77 build_basis code.
    ! allocate(Nstock(dim))
    ! kcp=0
    ! do i=1,NN
    !    call bdecomp(i,ivec) !|0.1.1.1.0.0>
    !    jup = sum(ivec(1:Ns))    !conta quanti up in questa configurazione
    !    jdw = sum(ivec(Ns+1:2*Ns)) !conta quanti dw " " " 
    !    if(jup==nup.AND.jdw==ndw)then
    !       kcp=kcp+1           !count the states in the sector (n_up,n_dw)//kcp
    !       do is=1,Ntot
    !          if(ivec(is)==1)Nstock(kcp)=ibset(Nstock(kcp),is-1)
    !       enddo
    !    endif
    ! enddo
    ! print*,""
    ! do i=1,dim
    !    print*,Nstock(i),Hmap(isector)%map(i)
    ! enddo
    ! stop
    !call build_sector(isector)
    !stop 
    ! !!<DEBUG



    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo

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
  !PURPOSE  : constructs the pointers for the different sectors and
  !the vectors isrt and jsrt with the corresponding
  !ordering definition of the sub-basis within each sector
  !+------------------------------------------------------------------+
  !|ImpUP,BathUP;,ImpDW,BathDW >
  subroutine build_sector(isector,map)!,invmap)
    integer :: i,j,isector,iup,idw,count,dim
    integer :: nup,ndw
    integer :: ivec(Ntot)
    integer,dimension(:)  :: map
    !integer,dimension(NN),optional :: invmap
    nup = getnup(isector)
    ndw = getndw(isector)
    dim = getdim(isector)
    ! if(allocated(Hmap))deallocate(Hmap)
    ! if(allocated(invHmap))deallocate(invHmap)
    ! allocate(Hmap(dim),invHmap(NN))
    !if(present(invmap))invmap=0
    count=0
    do i=1,NN
       call bdecomp(i,ivec)
       iup = sum(ivec(1:Ns))
       idw = sum(ivec(Ns+1:2*Ns))
       if(iup==nup.AND.idw==ndw)then
          count             = count+1 !count the states in the sector (n_up,n_dw)
          map(count)        = i       !build the map to full space states
          !if(present(invmap))invmap(i)         = count
       endif
    enddo
  end subroutine build_sector
  !


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
    do l=0,Ntot-1                  !loop sul numero di "siti"
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
    return 
  end subroutine bdecomp




  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(m,i,j)!,sgn)
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
          ! sgn=(-1.d0)**km
          ! j=(i-2**(m-1))!*isg
       endif
    endif
    return
  end subroutine c



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j)!,sgn)
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
          ! sgn=(-1.d0)**km
          ! j=(i+2**(m-1))!*isg
       endif
    endif
    return
  end subroutine cdg





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine search_mu(ntmp,convergence)
    logical,intent(inout) :: convergence
    real(8)               :: ntmp
    logical               :: check
    integer,save          :: count=0
    integer,save          :: nindex=0
    real(8)               :: ndelta1,nindex1
    if(count==0)then
       inquire(file="searchmu_file.restart",exist=check)
       if(check)then
          open(10,file="searchmu_file.restart")
          read(10,*)ndelta,nindex          
          close(10)
       endif
    endif
    count=count+1
    nindex1=nindex
    ndelta1=ndelta
    if((ntmp >= nread+nerr))then
       nindex=-1
    elseif(ntmp <= nread-nerr)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
       ndelta=ndelta1/2.d0 !decreasing the step       
    else
       ndelta=ndelta1
    endif
    xmu=xmu+real(nindex,8)*ndelta
    write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",ntmp," /",nread,&
         "| shift=",nindex*ndelta,"| xmu=",xmu
    write(*,"(A,f15.12)")"dn=",abs(ntmp-nread)
    print*,""
    if(abs(ntmp-nread)>nerr)then
       convergence=.false.
       ! else
       !    convergence=.true.
    endif
    print*,""
    print*,"Convergence:",convergence
    print*,""
    open(10,file="searchmu_file.restart.new")
    write(10,*)ndelta,nindex
    close(10)
  end subroutine search_mu




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
