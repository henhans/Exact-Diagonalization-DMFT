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


end MODULE ED_CHI2FIT
