!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
!subroutine search_chemical_potential(ntmp,converged)
subroutine search_mu(ntmp,niter,converged)
  real(8),intent(in)    :: ntmp
  integer,intent(in)    :: niter
  logical,intent(inout) :: converged
  logical               :: check
  integer,save          :: count=0
  integer,save          :: nindex=0
  integer               :: nindex1
  real(8)               :: ndelta1
  integer,save          :: nth_magnitude=-1,nth_magnitude_old=-1
  real(8),save          :: nth=1.d-1
  logical,save          :: ireduce=.true.
  !
  !check actual value of the density *ntmp* with respect to goal value *nread*
  count=count+1
  nindex1=nindex
  ndelta1=ndelta
  if((ntmp >= nread+nth))then
     nindex=-1
  elseif(ntmp <= nread-nth)then
     nindex=1
  else
     nindex=0
  endif
  if(nindex1+nindex==0.AND.nindex/=0)then !avoid loop forth and back
     ndelta=ndelta1/2.d0 !decreasing the step       
  else
     ndelta=ndelta1
  endif
  !
  !update chemical potential
  xmu=xmu+dble(nindex)*ndelta
  !
  !Print information
  write(*,"(A,f16.9,A,f15.9)")" n    =",ntmp," /",nread
  if(nindex>0)then
     write(*,"(A,es16.9,A)")" shift=",nindex*ndelta," --->"
  elseif(nindex<0)then
     write(*,"(A,es16.9,A)")" shift=",nindex*ndelta," <---"
  else
     write(*,"(A,es16.9,A)")" shift=",nindex*ndelta," -- "
  endif
  write(*,"(A,f15.9)")" xmu  =",xmu
  write(*,"(A,ES16.9,A,ES16.9)")" dn   =",abs(ntmp-nread),"/",nth
  !
  !check convergence within actual threshold
  !if reduce is activetd
  !if density is in the actual threshold
  !if DMFT is converged
  if(ireduce.AND.abs(ntmp-nread)<nth.AND.converged)then
     nth_magnitude_old=nth_magnitude        !save old threshold magnitude
     nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude
     nth=max(nerr,10.d0**(nth_magnitude)) !set the new threshold 
     count=0                                !reset the counter
  endif
  !
  !if density is not converged set convergence to .false.
  if(abs(ntmp-nread)>nth)converged=.false.
  !
  !if you can not converge dmft for the smallest threshold
  !set it back to the previous larger one.
  if(count>niter.AND.nth==nerr.AND.ireduce.AND..not.converged)then
     ireduce=.false.
     nth=10.d0**(nth_magnitude_old) !previous converging threshold
  endif
  print*,""
  print*,"Converged:",converged,count
  print*,""
  !
end subroutine search_mu



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

