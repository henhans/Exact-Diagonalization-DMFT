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
  if(abs(ntmp-nread)>nerr)convergence=.false.
  write(*,"(A,f15.12,A,f15.12,A,f15.12,A,f15.12)")" n=",ntmp," /",nread,&
       "| shift=",nindex*ndelta,"| xmu=",xmu
  write(*,"(A,f15.12)")"dn=",abs(ntmp-nread)
  print*,""
  print*,"Convergence:",convergence
  print*,""
  open(10,file="searchmu_file.restart.new")
  write(10,*)ndelta,nindex
  close(10)
end subroutine search_mu
