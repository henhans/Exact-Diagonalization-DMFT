MODULE ED_AUX_FUNX
  !########################################################################
  !PROGRAM  : FUNCS_GLOBAL
  !TYPE     : Module
  !PURPOSE  : Constructs some functions used in other places. 
  !AUTHORS  : Adriano Amaricci
  !LAST UPDATE: 10/2009
  !########################################################################
  USE ED_VARS_GLOBAL
  implicit none
  private
  public :: imp_sectorns,bdecomp,cdg,c,delta_and
  public :: init_bath_ed,search_mu

contains

  subroutine search_mu(ntmp,convergence)
    logical,intent(inout) :: convergence
    real(8)               :: ntmp
    nindex1=nindex
    ndelta1=ndelta
    if((ntmp >= nread+nerr))then
       nindex=-1
    elseif(ntmp <= nread-nerr)then
       nindex=1
    else
       nindex=0
    endif
    if(nindex1+nindex==0)then !avoid loop forth and back
       ndelta=ndelta1/2.d0 !decreasing the step
       xmu=xmu+real(nindex,8)*ndelta
    else
       ndelta=ndelta1
       xmu=xmu+real(nindex,8)*ndelta
    endif
    write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",ntmp," /",nread,&
         "| shift=",nindex*ndelta,"| xmu=",xmu
    write(*,"(A,f15.12)")"dn=",abs(ntmp-nread)
    print*,""
    if(abs(ntmp-nread)>nerr)convergence=.false.
  end subroutine search_mu

  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_bath_ed
    integer :: i
    logical :: IOfile
    !Initialize the parameter for every mu-loop                  
    inquire(file=trim(Hfile),exist=IOfile)
    if(.NOT.IOfile)then
       call guess_bath_params
    else
       write(*,*)'Reading the seed from file'       
       open(51,file=trim(Hfile))
       read(51,*)xmu
       do i=1,Nbath
          read(51,*)epsiup(i),epsidw(i),vup(i),vdw(i)
       enddo
       close(51)
    endif
  end subroutine init_bath_ed



  !+------------------------------------------------------------------+
  !PURPOSE  : Build the parameters for the Hamiltonian
  !+------------------------------------------------------------------+
  subroutine guess_bath_params
    integer :: i,n2
    write(*,*)'Generating the seed'
    n2=Nbath/2;if(n2==0)n2=1
    do i=0,Nbath-1            !(1->NC=Nbath/2=#sites in each Bath)
       epsiup(i+1)=2.d0*dfloat(i-1-n2)/dfloat(n2) !d/2.d0+heff
       vup(i+1)=dsqrt(1.d0/dfloat(Nbath)) !d**2/4.
       epsidw(i+1)=epsiup(i+1)
       vdw(i+1)=vup(i+1)
    enddo
  end subroutine guess_bath_params



  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the pointers for the different sectors and
  !the vectors isrt and jsrt with the corresponding
  !ordering definition of the sub-basis within each sector
  !+------------------------------------------------------------------+
  !|ImpUP,BathUP;,ImpDW,BathDW >
  subroutine imp_sectorns(in,is,idg,imap,invmap)
    integer :: i,j,in,is,idg,NR
    integer :: ibn,ibs
    integer :: imap(:),invmap(:)
    integer :: ib2(N)
    idg=0
    imap=0
    invmap=0
    if(size(invmap)/=NN)stop "ERROR1 in imp_sectorns"
    if(size(imap)/=NP)stop "ERROR2 in imp_sectorns"
    NR=Nimp+1
    do i=1,NN
       call bdecomp(i,ib2)
       ibs=0
       ibn=0
       ibs=ib2(1) - ib2(1+Ns) 
       ibn=ib2(1) + ib2(1+Ns) 
       if(Nimp==2)then
          ibs=ibs + ib2(2) - ib2(2+Ns)
          ibn=ibn + ib2(2) + ib2(2+Ns)
       endif
       !add more if Nimp>2
       do j=NR,Ns
          ibs=ibs+ib2(j)-ib2(j+Ns)
          ibn=ibn+ib2(j)+ib2(j+Ns)
       enddo
       if(ibn==in.AND.ibs==is)then
          idg=idg+1           !count the states in the sector (n,s)
          imap(idg)=i         !build the map to full space states
          invmap(i)=idg       !and the inverse map
       endif
    enddo
    return
  end subroutine imp_sectorns
  !==================================================================
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************





  !+------------------------------------------------------------------+
  !PROGRAM  : bdecomp
  !TYPE     : subroutine
  !PURPOSE  : input a state |i> and output a vector ivec(N)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  subroutine bdecomp(i,ivec)
    integer :: ivec(N)         
    integer :: l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,N>
    !obtained from binary decomposition of the state/number i\in 2^N
    do l=0,N-1                  !loop sul numero di "siti"
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
    return 
  end subroutine bdecomp
  !=======================================================================
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : C
  !TYPE     : subroutine
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(m,i,j)
    integer :: ib(N)
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
  !=======================================================================
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************




  !+-------------------------------------------------------------------+
  !PROGRAM  : CDG
  !TYPE     : subroutine
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j)
    integer :: ib(N)
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
  !=======================================================================
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************




  pure function delta_and(x,epsi,vi) result(fg)
    complex(8),intent(in)            :: x
    real(8),dimension(:),intent(in)  :: epsi,vi
    complex(8)                       :: fg
    integer                          :: i
    fg=zero
    do i=1,size(epsi)
       fg=fg+vi(i)**2/(x-epsi(i))
    enddo
  end function delta_and



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************



  ! !+------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  : 
  ! !+------------------------------------------------------------------+
  ! subroutine getloop_range(initloop)
  !   integer                        :: initloop
  !   integer                        :: sectchunk,sum,id,isloop     
  !   integer,dimension(0:mpiSIZE-1) :: lloop
  !   startloop=initloop
  !   sectchunk=NN/mpiSIZE;sum=0;id=0
  !   do isloop=startloop,Nsect
  !      sum=sum+deg(isloop)
  !      if(sum > sectchunk)then
  !         lloop(id)=isloop
  !         sum=0;id=id+1
  !      endif
  !   enddo
  !   lloop(mpiSIZE-1)=Nsect
  !   if(mpiID==0)lastloop=lloop(0)
  !   do id=1,mpiSIZE-1
  !      if(mpiID==id)then
  !         startloop=lloop(id-1)+1
  !         lastloop=lloop(id)
  !      endif
  !   enddo
  ! end subroutine getloop_range
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************

END MODULE ED_AUX_FUNX
