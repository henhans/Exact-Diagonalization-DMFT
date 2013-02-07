program CuO2fulled
  !########################################################################
  !PROGRAM  : HMED
  !TYPE     : main code
  !PURPOSE  : Complete ED with STAR GEOMETRY for CuO2 model (Emery-OKA)
  !AUTHORS  : Adriano Amaricci, L. de'Medici,G.Sordi, M.Rozenberg
  !########################################################################
  !LOCAL:
  USE VARS_GLOBAL
  USE DIAG
  USE TOFITGF
  USE MPI
  implicit none

  !LOCAL VARIABLES
  integer :: success

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)

  !Fix the size of the problem && read params:
  !----------------------------------------------------------------  
  call readval(INunit,.false.)
  call size_alloc()   

  !Store initial parameters:
  xmu0=xmu+ed0;xmu=xmu0

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  print*,'Processor ',mpiID,' of ',mpiSIZE,' is alive'

  success=0
  !Starts DMFT iteration
  do iloop=1,nloop
     !Solve the EFFECTIVE IMPURITY PROBLEM
     !################################################################
     call ed_solver(model,order,searchmode,lat_label)

     !Perform SELF-CONSISTENCY and Check convergency
     !################################################################
     if(mpiID == 0)then
        call get_delta(trim(model))  !Obtain the \Delta_loc
        call Mconvergence(success)
        call fit_gf()
     endif
  enddo

  !Save the BATH and finalize the calculation
  !################################################################
  call finalize(success)
  CALL MPI_FINALIZE(mpiERR)
end program CuO2fulled



