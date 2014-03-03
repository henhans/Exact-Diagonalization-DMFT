MODULE ED_BATH_TYPE
  implicit none
  type effective_bath
     real(8),dimension(:,:,:),allocatable :: e
     real(8),dimension(:,:,:),allocatable :: v
     real(8),dimension(:,:,:),allocatable :: d
     logical                              :: status=.false.
  end type effective_bath
END MODULE ED_BATH_TYPE
