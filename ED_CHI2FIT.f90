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





end MODULE ED_CHI2FIT
