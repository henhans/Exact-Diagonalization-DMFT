program fullED2DEG
  !########################################################################
  !PROGRAM  : HMED
  !TYPE     : main code
  !PURPOSE  : Complete ED with STAR GEOMETRY for the Hubbard model
  !AUTHORS  : Adriano Amaricci, G.Sordi, M.Rozenberg
  !COMMENTS : F95
  !########################################################################
  USE VARS_GLOBAL
  USE DIAG
  USE TOFITGF
  USE TOOLS
  USE MPI

  implicit none

  !LOCAL VARIABLES
  integer :: success

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)

  !Fix the size of the problem && read params:
  !----------------------------------------------------------------  
  call readval(INunit,.FALSE.)
  call size_alloc()   

  !Store initial parameters:
  nread=1.d0+doping2deg
  extension="v.ed"
  xmu0=xmu
  xmuA=xmu0
  xmuB=xmu0

  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  print*,'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  success=0
  !Starts DMFT iteration
  do iloop=1,nloop

     !SOLVE LATTICE A
     !################################################################
     call dump("Solve lattice A")
     xmu=xmuA
     lat_label="A"
     muflag=0
     call ed_solver(model,order,searchmode,lat_label,TT)
     if(mpiID == 0)then
        call get_delta(deltaMODE)        
        call Mconvergence(success,lat_label)
        call fit_gf(order,lat_label)
     endif
     zna=npimp
     xmuB=xmu0-v*zna

     !SOLVE LATTICE B
     !################################################################
     call dump("",3)
     call dump("Solve lattice B")
     xmu=xmuB
     lat_label="B"
     muflag=1
     call ed_solver(model,order,searchmode,lat_label,FF)
     if(mpiID == 0)then
        call get_delta(deltaMODE)        
        call Mconvergence(success,lat_label)
        call fit_gf(order,lat_label)
     endif
     znb=npimp
     xmuA=xmu0-v*znb
  enddo

  call finalize(success)
  CALL MPI_FINALIZE(mpiERR)
end program fullED2DEG



