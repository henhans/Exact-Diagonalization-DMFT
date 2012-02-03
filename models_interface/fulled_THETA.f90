program THETAed
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
  real(8) :: beta0,u0,v0,nread0
  real(8) :: DMFTErr, ConvErr
  integer :: success
  character(len=128) :: cmd
  character(len=1)   :: namedir

  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  !Fix the BATH size
  Nbath=Ns-1

  !Fix DMFT convergence error
  ConvErr = 1.d-6

  !read SOME init  parameters
  call readval(INunit)
  if(mpiID==0)call dumpinival(Nbath,.false.)
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
  print*,'Processor ',mpiID,' of ',mpiSIZE,' is alive'


  !Store initial parameters:
  beta0=beta;temp=1.d0/beta
  u0 =u
  v0 =v
  xmu0 =xmu
  nread=1.d0
  nread0=nread
  extension="v.ed"
  muflag=.FALSE.;if(nread /= 0.d0)muflag=.TRUE.
  lm=max(lm,Nfreq);lm=int(lm/(pi2*temp)) !Get maxNf at lowest T


  gmu0=xmu0
  xmu=xmu0
  xmuA=xmu0
  xmuB=xmu0
  DMFTerr=1.d4 
  success=0
  !Starts DMFT iteration
  do iloop=1,nloop
     include "init_dmft.f90"	
     if(iloop==1)then
        allocate(epsiAup(Nbath),epsiBup(Nbath),&
             vAup(Nbath),vBup(Nbath))

        allocate(epsiAdw(Nbath),epsiBdw(Nbath),&
             vAdw(Nbath),vBdw(Nbath))

        epsiAup=epsiup-heff;epsiBup=epsiup-heff
        vAup=vup;vBup=vup

        epsiAdw=epsiup+heff;epsiBdw=epsiup+heff
        vAdw=vdw;vBdw=vdw
     endif

     !SOLVE LATTICE A
     !-------------------------------------------
     epsiup=epsiAup;epsidw=epsiAdw
     vup=vAup;vdw=vAdw
     xmu=xmuA
     call imp1_diag(0)
     call imp1_getfunx_LRO(0)
     !Calcola le osservabili:
     if(mpiID==0)call imp1_getout_latA
     deallocate(storvec)
     if(mpiID == 0)then
        if(last)call hm_veogf_LRO
        call mix_gf_latA_LRO()
        !Get the Delta_loc 
        call get_delta_BLHM_LRO()
        !Perform the chi^2 fit to update the bath
        call fit_latA_LRO() 
     endif
     zna=npimp
     xmuB=xmu0-v*zna


     !SOLVE LATTICE B
     !-------------------------------------------
     epsiup=epsiBup;epsidw=epsiBdw
     vup=vBup;vdw=vBdw
     xmu=xmuB
     call imp1_searchmu_CO
     call imp1_getfunx_LRO(0)
     !Calcola le osservabili:
     if(mpiID==0)call imp1_getout_latB
     deallocate(storvec)
     if(last .AND. mpiID==0)call hm_veogf_LRO
     if(iloop < nloop)then 
        if(mpiID == 0)then
           call mix_gf_latB_LRO()
           !Get the Delta_loc 
           call get_delta_BLHM_LRO()
           !Perform the chi^2 fit to update the bath
           call fit_latB_LRO() 
        endif
     endif
     znb=npimp
     xmuA=xmu0-v*znb
  enddo


  if(mpiID==0)call writeval(USEDinput)
  CALL MPI_FINALIZE(mpiERR)
end program THETAed



