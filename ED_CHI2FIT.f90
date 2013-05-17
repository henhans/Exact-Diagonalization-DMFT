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
    complex(8),dimension(:)         :: fdelta
    real(8),dimension(:)            :: epsi,vi
    real(8),dimension(size(fdelta)) :: wm
    integer  :: Lw
    Lw=Nfit
    if(Nfit>size(fdelta))then
       Lw=size(fdelta)
       call msg("Fitting with "//txtfy(Lw))
    endif
    wm = pi/beta*dble(2*arange(1,Lw)-1)
    call fitgreen(wm(1:Lw),fdelta(1:Lw),epsi,vi)
    call dump_fit_result(wm(1:Lw),fdelta(1:Lw),epsi,vi)
  end subroutine chi2_fitgf


  !********************************************************************
  !********************************************************************
  !********************************************************************


  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine dump_fit_result(wm,fg,epsi,vi)
    complex(8),dimension(:)        :: fg
    real(8),dimension(size(fg))    :: wm
    complex(8),dimension(size(fg)) :: fgand
    real(8),dimension(:)           :: epsi,vi
    integer                        :: i,j,Lw
    Lw=size(fg)
    fgand=zero
    do i=1,Lw
       fgand(i) = delta_and(xi*wm(i),epsi,vi)
    enddo
    call splot("fit_delta.ed",wm,dimag(fg),dimag(fgand),dreal(fg),dreal(fgand))
  end subroutine dump_fit_result



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_CHI2FIT
