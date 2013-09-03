!########################################################################
!PROGRAM  : ED_AUX_FUNX
!AUTHORS  : Adriano Amaricci
!########################################################################
MODULE ED_AUX_FUNX
  USE ED_VARS_GLOBAL
  implicit none
  private
  public :: build_sector
  public :: bdecomp
  public :: c,cdg
  public :: setup_pointers
  public :: search_mu

contains




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_pointers
    integer                          :: in,dim,isector,jsector,dimup,dimdw
    integer                          :: nup,ndw,jup,jdw
    integer,dimension(:),allocatable :: imap
    integer,dimension(:),allocatable :: invmap
    call msg("Setting up pointers:")
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
          allocate(Hmap(isector)%map(dim))
          call build_sector(nup,ndw,dim,&
               Hmap(isector)%map(:),invHmap(isector,:))
       enddo
    enddo
    call stop_timer

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
  subroutine build_sector(iup,idw,dim,imap,invmap)
    integer :: i,j,iup,idw,dim
    integer :: nup,ndw
    integer :: imap(:),invmap(:)
    integer :: ivec(Ntot)
    dim=0
    imap=0
    invmap=0
    if(size(invmap)/=NN)stop "ERROR1 in imp_sectorns" 
    do i=1,NN
       call bdecomp(i,ivec)
       nup = sum(ivec(1:Ns))
       ndw = sum(ivec(Ns+1:2*Ns))
       if(nup==iup.AND.ndw==idw)then
          dim=dim+1           !count the states in the sector (n_up,n_dw)
          imap(dim)=i         !build the map to full space states
          invmap(i)=dim       !and the inverse map
       endif
    enddo
    if(size(imap)/=dim)stop "ERROR2 in imp_sectorns" 
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
  subroutine c(m,i,j)
    integer :: ib(Ntot)
    integer :: i,j,m,km,k
    integer :: isg
    call bdecomp(i,ib)
    if (ib(m).eq.0)then
       j=0
    else
       if(m.eq.1)then
          j=i-1
       else
          km=0
          do k=1,m-1
             km=km+ib(k)
          enddo
          isg=(-1)**km
          j=(i-2**(m-1))*isg
       endif
    endif
    return
  end subroutine c



  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j)
    integer :: ib(Ntot)
    integer :: i,j,m,km,k
    integer :: isg
    call bdecomp(i,ib)
    if (ib(m)==1)then
       j=0
    else
       if(m.eq.1)then
          j=i+1
       else
          km=0
          do k=1,m-1
             km=km+ib(k)
          enddo
          isg=(-1)**km
          j=(i+2**(m-1))*isg
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






END MODULE ED_AUX_FUNX
