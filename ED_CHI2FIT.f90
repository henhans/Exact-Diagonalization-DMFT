!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE CGFIT
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, ONLY:delta_and
  implicit none
  private

  public :: chi2_fitgf

contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf(fdelta,epsi,vi)
    complex(8),dimension(:) :: fdelta
    real(8),dimension(:)    :: epsi,vi
    call fitgreen(wm(1:Nfit),fdelta(1:Nfit),epsi,vi)
    call dump_fit_result(fdelta(1:Nfit),epsi,vi)
  end subroutine chi2_fitgf


  !********************************************************************
  !********************************************************************
  !********************************************************************


  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine dump_fit_result(fg,epsi,vi)
    complex(8),dimension(:)        :: fg
    complex(8),dimension(size(fg)) :: fgand
    real(8),dimension(:)    :: epsi,vi
    integer                 :: i,j
    fgand=zero
    do i=1,Nfit
       fgand(i) = delta_and(xi*wm(i),epsi,vi)
    enddo
    call splot("fit_delta.ed",wm(1:Nfit),aimag(fg(1:Nfit)),aimag(fgand(1:Nfit)),real(fg(1:Nfit)),real(fgand(1:Nfit)))
  end subroutine dump_fit_result



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


  ! !+-------------------------------------------------------------+
  ! !PURPOSE  : 
  ! !+-------------------------------------------------------------+
  ! function scc_fit(loop)
  !   integer                   :: loop,scc_fit
  !   call Mconvergence(loop) !this writes on success
  !   call fit_gf(loop)
  !   scc_fit=success
  !   return
  ! end function scc_fit



  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine Mconvergence(loop)
  !   real(8) :: mismatch
  !   integer :: loop
  !   real(8) :: sum
  !   integer :: i
  !   if(loop==nloop)return
  !   mismatch=0.d0
  !   if(loop==1)then
  !      mismatch=1.d0
  !   else
  !      sum=0.d0
  !      select case(order)
  !      case default
  !         do i=1,NL
  !            mismatch=mismatch+abs(delta(i)-gold(i))
  !            sum=sum+abs(delta(i))
  !         end do
  !      case("co")
  !         if(lat_label == "A")then
  !            do i=1,NL
  !               mismatch=mismatch+abs(delta(i)-gaold(i))
  !               sum=sum+abs(delta(i))
  !            end do
  !         elseif(lat_label == "B")then
  !            do i=1,NL
  !               mismatch=mismatch+abs(delta(i)-gbold(i))
  !               sum=sum+abs(delta(i))
  !            end do
  !         endif
  !      case("lro")
  !         do i=1,NL
  !            mismatch=mismatch+abs(deltaup(i)-goldu(i))+abs(deltadw(i)-goldd(i))
  !            sum=sum+abs(deltaup(i)+deltadw(i))
  !         enddo
  !      end select
  !      mismatch=mismatch/sum
  !      write(*,"(A,f15.12)") "DMFT-Error=",mismatch
  !      if(mismatch < ConvErr)success=success+1
  !   end if
  ! end subroutine  Mconvergence
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************







  ! !+-------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  : 
  ! !+-------------------------------------------------------------+
  ! subroutine fit_gf(loop)
  !   integer :: loop
  !   call mix_gf_std(loop)
  !   if(trim(order)=="lro")then
  !      select case(lat_label)
  !      case default
  !         call fit_gf_lro()
  !      case("A")
  !         call fit_gf_co_lro()
  !      case("B")
  !         call fit_gf_co_lro()
  !      end select
  !      return
  !   else
  !      call fit_gf_std()
  !      return
  !   endif
  ! end subroutine fit_gf
  ! include "fit_gf.h"
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************




  ! !+-------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  :On INPUT Delta = Impurity Green function 
  ! !+-------------------------------------------------------------+
  ! subroutine mix_gf_std(loop)
  !   integer :: loop
  !   if(trim(order)=="lro")then    !LRO ordering:
  !      select case(lat_label)
  !      case default
  !         if(.not.allocated(goldu))allocate(goldu(NL))
  !         if(.not.allocated(goldd))allocate(goldd(NL))
  !         if(loop > 1)then
  !            deltaup=weigth*deltaup+(1.d0-weigth)*goldu(1:size(deltaup))
  !            deltadw=weigth*deltadw+(1.d0-weigth)*goldd(1:size(deltadw))
  !         endif
  !         goldu(1:size(deltaup))=deltaup
  !         goldd(1:size(deltadw))=deltadw
  !      case("A")
  !         if(.not.allocated(gaoldu))allocate(gaoldu(NL))
  !         if(.not.allocated(gaoldd))allocate(gaoldd(NL))
  !         if(iloop > 1)then
  !            deltaup=weigth*deltaup+(1.d0-weigth)*gaoldu(1:size(deltaup))
  !            deltadw=weigth*deltadw+(1.d0-weigth)*gaoldd(1:size(deltadw))
  !         endif
  !         gaoldu(1:size(deltaup))=deltaup
  !         gaoldd(1:size(deltadw))=deltadw
  !      case("B")
  !         if(.not.allocated(gboldu))allocate(gboldu(NL))
  !         if(.not.allocated(gboldd))allocate(gboldd(NL))
  !         if(iloop > 1)then
  !            deltaup=weigth*deltaup+(1.d0-weigth)*gboldu(1:size(deltaup))
  !            deltadw=weigth*deltadw+(1.d0-weigth)*gboldd(1:size(deltadw))
  !         endif
  !         gboldu(1:size(deltaup))=deltaup
  !         gboldd(1:size(deltadw))=deltadw
  !      end select
  !   else                          !NO LRO ordering
  !      select case(lat_label)
  !      case default
  !         if(.not.allocated(gold))allocate(gold(NL))
  !         if(loop > 1)delta=weigth*delta+(1.d0-weigth)*gold(1:size(delta))
  !         gold(1:size(delta))=delta
  !      case("A")
  !         if(.not.allocated(gaold))allocate(gaold(NL))     
  !         if(loop > 1)delta=weigth*delta+(1.d0-weigth)*gaold(1:size(delta)) !mix delta functions
  !         gaold(1:size(delta))=delta
  !      case("B")
  !         if(.not.allocated(gbold))allocate(gbold(NL))
  !         if(loop > 1)delta=weigth*delta+(1.d0-weigth)*gbold(1:size(delta)) !mix 
  !         gbold(1:size(delta))=delta
  !      end select
  !   endif
  ! end subroutine mix_gf_std
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************



  ! !+-------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  : 
  ! !+-------------------------------------------------------------+
  ! subroutine dump_delta(x,fg,filename)
  !   complex(8),dimension(:) :: x,fg
  !   character(len = *) :: filename
  !   integer :: i
  !   open(70,file=trim(filename))
  !   do i=1,NL
  !      write(70,*)aimag(x(i)) , aimag(fg(i)), real(fg(i))
  !   enddo
  !   close(70)
  !   return
  ! end subroutine dump_delta
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************





  ! !+-------------------------------------------------------------------+
  ! !PROGRAM  : 
  ! !TYPE     : subroutine
  ! !PURPOSE  : 
  ! !+-------------------------------------------------------------------+
  ! subroutine get_g0and(eps,vps)
  !   complex(8),dimension(NL) :: gg,ome
  !   real(8),dimension(:) :: eps,vps
  !   integer :: i,j
  !   complex(8) :: iw
  !   gg=zero
  !   do i=1,NL
  !      iw=xi*pt*dble(2*i-1);ome(i)=iw
  !      do j=1,Nbath
  !         gg(i)=gg(i)+vps(j)**2/(iw-eps(j))
  !      enddo
  !   enddo
  !   call dump_delta(ome,gg,trim(Dimpfile))
  ! end subroutine get_g0and
  ! !*********************************************************************
  ! !*********************************************************************
  ! !*********************************************************************



end MODULE ED_CHI2FIT
