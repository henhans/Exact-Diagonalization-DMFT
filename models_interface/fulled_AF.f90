program magHMED
  !########################################################################
  !PROGRAM  : HMED
  !TYPE     : main code
  !PURPOSE  : Complete ED with STAR GEOMETRY for the Hubbard model
  !AUTHORS  : Adriano Amaricci, G.Sordi, M.Rozenberg
  !COMMENTS : F95
  !########################################################################
  USE VARS_GLOBAL
  USE DIAG
  USE SEARCHMU
  USE GETFUNX
  USE GETOUT
  USE VEOGF
  USE TOFITGF
  USE MPI
  implicit none

  !LOCAL VARIABLES
  real(8) :: DMFTErr
  integer :: success

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)

  !Fix the size of the problem && read params:
  !----------------------------------------------------------------  
  call readval(INunit)
  if(mpiID==0)call dumpinival(.false.)
  call size_alloc()   
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  print*,'Processor ',mpiID,' of ',mpiSIZE,' is alive'


  !Store initial parameters:
  if(heff < 0.d0)heff=-heff !to always havea positive field
  heff=0.5d0*heff           !spin multiplication

  DMFTerr=1.d4 
  success=0
  !Starts DMFT iteration
  do iloop=1,nloop
     !fix some parameters - INIT the bath
     !################################################################
     include "init_dmft.f90"	

     !Solve the EFFECTIVE IMPURITY PROBLEM
     !################################################################
     select case (muflag)
     case default
        call imp_diag(muflag)
     case (1)
        call imp_searchmu
     end select
     call imp_getfunx_LRO(Nimp,0,"lro")
     if(mpiID==0)call imp_getout
     deallocate(storvec)
     if(last .AND. mpiID==0 .AND. veoflag)call hm_veogf_LRO

     !Perform SELF-CONSISTENCY and Check convergence
     !################################################################
     if(mpiID == 0)then
        call get_delta_BLHM_LRO()
        if(iloop < nloop)then 
           call Mconvergence(DMFTerr,"LRO");print*, "DMFT-Error=",DMFTerr
           if(DMFTerr < ConvErr)success=success+1
           call mix_gf_LRO()
           call fit_LRO() 
        endif
     endif
  enddo

  !Save the BATH and finalize the calculation
  !################################################################
  if(mpiID==0)then
     Print*,"Bravo hai totalizzato",success," successi"
     cmd = 'mv '//trim(Hfile)//'.new '//trim(Hfile)
     call system(cmd)
     call writeval(USEDinput)
  endif

  CALL MPI_FINALIZE(mpiERR)
end program magHMED



