  include 'eispack.f90'
  module HLZPACK_LANCZOS
    USE ED_VARS_GLOBAL
    USE ED_AUX_FUNX, only: HtimesV
    implicit none
    private
    public :: lanczos_hlzpack
    public :: lanczos_gs
    public :: lanczos_simple
    public :: lanczos_one_step
  contains


    subroutine lanczos_simple(vin,vout,alfa,beta,Niter)
      real(8),dimension(:),intent(inout) :: vin,vout
      real(8),dimension(:),intent(inout) :: alfa,beta
      integer,intent(inout)              :: Niter
      integer                            :: iter,nlanc
      real(8)                            :: dummy,a_,b_
      a_=zero ; b_=zero ; nlanc=0
      do iter=1,size(alfa)
         nlanc=nlanc+1
         call lanczos_one_step(vin,vout,a_,b_,iter)
         alfa(iter)=a_
         if(iter/=size(beta))beta(iter+1)=b_
         if(abs(b_)<1.d-10)exit
      enddo
      Niter=Nlanc
    end subroutine lanczos_simple


    subroutine lanczos_one_step(vin,vout,a,b,iter)
      real(8),dimension(:),intent(inout) :: vin,vout
      real(8),intent(inout)              :: a,b
      integer                            :: iter,i
      real(8)                            :: dummy
      if(iter>1)then
         do i=1,size(vin)
            dummy=vin(i)
            vin(i)=vout(i)/b
            vout(i)=-b*dummy
         end do
      endif
      call HtimesV(vin,vout)
      a  = dot_product(vout,vin)
      vout  = vout - a*vin
      b = sqrt(dot_product(vout,vout))
    end subroutine lanczos_one_step


    function lanczos_hlzpack(N,NREIG,NSTEPS,LPRNT,LN,EVAL,EVEC,NECONV,NLANC,V,BFLAG) result(status)
      ! C* from file hlzdrd.f in hlzpack/src/
      ! C*  - ARGUMENTS :                                                      *
      ! C*    LFLAG : flag for the Hermitian-Lanczos algorithm                 *
      ! C*            < 0, non standard exit (error condition)                 *
      ! C*            = 0, must be used on the first call to HLZDRD            *
      ! C*            = 1, the user must compute R=H*Q and call HLZDRD again   *
      ! C*            = 2, NREIG, NSTEPS or an invariant subspace reached      *
      ! C*            = 3, compute the eigenvectors                            *
      ! C*    LPSET : holds basic data for the Hermitian-Lanczos algorithm     *
      ! C*            LPSET(1)=N     , .gt.0, dimension of the matrix          *
      ! C*            LPSET(2)=NREIG , .gt.0, number of eigenpairs required    *
      ! C*            LPSET(3)=NSTEPS, .gt.0, maximum number of steps          *
      ! C*            LPSET(4)=LPRNT , .ge.0, level of printing                *
      ! C*            LPSET(5)=LN    , .ge.N, first dimension of X and BASIS   *
      ! C*    INFO  : holds information about the execution                    *
      ! C*            INFO(1)=J    , dimension of the basis                    *
      ! C*            INFO(2)=NCEIG, number of eigenpairs converged with J     *
      ! C*            INFO(3)=LWARN, code for warning messages                 *
      ! C*            INFO(4)=LRERR, code for error messages                   *
      ! C*    HNRM  : norm of H                                                *
      ! C*    EVAL   : eigenvalue approximations                               *
      ! C*    RNRM  : residual norms of the eigenpair approximations           *
      ! C*    W     : work array                                               *
      ! C*    Q     : Lanczos vector                                           *
      ! C*    R     : H*Q                                                      *
      ! C*    BASIS : stores the Lanczos vectors                               *
      ! C*    EVEC  : eigenvector approximations                               *
      !Control parameter:
      INTEGER                 :: LFLAG
      INTEGER,DIMENSION(5)    :: LPSET
      INTEGER,INTENT(IN)      :: N      !LPSET(1), Dimension of the problem (Hamiltonian and vector)
      INTEGER,INTENT(IN)      :: NREIG  !LPSET(2), Number of required eigenpairs
      INTEGER,INTENT(IN)      :: NSTEPS !LPSET(3), Max number of Lanczos steps
      INTEGER,INTENT(IN)      :: LPRNT  !LPSET(4), Printing level
      INTEGER,INTENT(IN)      :: LN     !LPSET(5), first dimension of X and Basis (usually = N)
      !Eigenvalues, Eigenvectors:
      REAL(8),INTENT(INOUT),DIMENSION(NSTEPS)      :: EVAL
      COMPLEX(8),INTENT(INOUT),DIMENSION(LN,NSTEPS) :: EVEC
      !Info array
      INTEGER,DIMENSION(5)    :: INFO(4)
      !Norm of H, residual norm:
      REAL(8)                 :: HNRM
      REAL(8),ALLOCATABLE     :: RNRM(:)
      !Work array:
      INTEGER                 :: MAXW,MAXW_MAX
      REAL(8),ALLOCATABLE     :: W(:)
      !Lanczos vector:
      COMPLEX(8),ALLOCATABLE  :: Q(:)
      COMPLEX(8),ALLOCATABLE  :: R(:) !H*Q
      !Eigenpairs:
      COMPLEX(8),ALLOCATABLE  :: BASIS(:,:)
      !INFO:
      INTEGER,OPTIONAL,INTENT(OUT) :: NECONV,NLANC
      !CONSTANTS
      INTEGER                 :: I,J
      REAL(8)                 :: IPART,RPART,Qnorm
      INTEGER                 :: WARN,IBIT,ICONV
      LOGICAL                 :: BOOL
      CHARACTER(LEN=64)       :: WMSG(0:5),EMSG(0:7)
      REAL(8),DIMENSION(:),OPTIONAL   :: V
      LOGICAL,OPTIONAL                :: BFLAG
      LOGICAL                         :: BFLAG_
      !
      INTEGER     :: STATUS
      !
      BFLAG_=.false.;if(present(BFLAG))BFLAG_=BFLAG
      LFLAG = 0      
      HNRM  = 0.0D0
      !NSTEPS=min(NSTEPS,N)
      !if(NSTEPS>N)call error("Error in HLZPACK_LANCZOS/lanczos_hlzpack: wrong number of Nsteps > N.")
      LPSET(1) = N
      LPSET(2) = NREIG
      LPSET(3) = NSTEPS
      LPSET(4) = LPRNT
      LPSET(5) = LN
      MAXW=(NSTEPS+10)*NSTEPS

      WMSG(0)="normal execution"
      WMSG(1)="nothing to be done: N=1"
      WMSG(2)="invariant subspace found at the j-th steps"
      WMSG(3)="no solution converged with NSTEPS"
      WMSG(4)="eigenvectors can not be computed"
      WMSG(5)="not enough solutions converged"

      EMSG(0)="normal execution"
      EMSG(1)="invalid specification of LFLAG"
      EMSG(2)="invalid specification of LPSET(1):N"
      EMSG(3)="invalid specification of LPSET(2):NREIG"
      EMSG(4)="invalid specification of LPSET(3):NSTEPS"
      EMSG(5)="invalid specification of LPSET(4):LPRNT"
      EMSG(5)="invalid specification of LPSET(5):LN"
      EMSG(5)="reduced eigenproblem could not be solved"


      allocate(RNRM(N),W(MAXW),Q(N),R(N),BASIS(LN,NSTEPS))

      if(LPRNT/=0)write(*,"(A28)")"LANCZOS                        "
      if(present(V))then
         if(size(V)/=N)call error("Error HLZPACK_LANCZOS/lanczos_hlzpack: optional vectore V has wrong dimension")         
         Q = V
      else         
         Q = cmplx(0.d0,0.d0)
      endif

      if(LPRNT/=0)then
         write(*,"(A28,I6)")"Required states:            ",NREIG
         write(*,"(A28,I6)")"Required Lanczos iterations:",NSTEPS
         write(*,"(A28,I6)")"Dimension of the problem:   ",N
      end if
      ICONV = 0
      HZLCYCLE: DO
         CALL HLZDRD (LFLAG,LPSET,INFO,HNRM,EVAL,RNRM,W,Q,R,BASIS,EVEC)
         IF (LFLAG < 0) THEN
            WRITE(*,"(A)") "INFO(1), INFO(2), INFO(3), INFO(4)"
            WRITE(*,"(4I8)") INFO(1),INFO(2),INFO(3),INFO(4)
            WRITE(*,"(A,I2,A)")"Warning! LANCZOS/lanczos_hzlpack_sparse: INFO(4)=",INFO(4),":"         
            do ibit=0,7
               if(btest(INFO(4),ibit))write(*,"(A,A1,1x)",advance='no')trim(adjustl(trim(EMSG(ibit)))),", "
            enddo
            write(*,*)
            IF ( INFO(4) .NE. 0 ) call error("Error in LANCZOS/lanczos_hzlpack: info(4)/=0")
         ELSE IF (LFLAG==1) THEN
            R=cmplx(0.d0,0.d0)
            call HtimesV(Q,R)
         ELSE
            exit HZLCYCLE
         END IF
      ENDDO HZLCYCLE
      CALL HLZDRD (3,LPSET,INFO,HNRM,EVAL,RNRM,W,Q,R,BASIS,EVEC)
      if(BFLAG_)EVEC=BASIS

      !OUTPUT:
      ICONV= INFO(2)
      WARN = INFO(3)      
      if(present(NLANC))NLANC=INFO(1)
      if(present(NECONV))NECONV=INFO(2)
      if(WARN/=0 .AND. LPRNT/=0)then
         write(*,"(A,I2,A)",advance='no')"Warning! LANCZOS/lanczos_hzlpack_sparse: INFO(3)=",WARN,":"
         do ibit=0,5
            if(btest(WARN,ibit))write(*,"(A,A1,1x)",advance='no')trim(adjustl(trim(WMSG(ibit)))),", "
         enddo
         write(*,*)
      endif
      if(ICONV<NREIG)then
         if(LPRNT/=0)write(*,"(A,I3,A,I3,A)")"Warning! converged EIGENSTATES:",ICONV,' over',NREIG,' required'
         STATUS=1
      else
         if(LPRNT/=0)write(*,"(A,I3,A,I3,A)")"Lanczos converged ",ICONV," eigen-states, with ",NLANC," iterations"
         STATUS=0
      endif
      if(LPRNT/=0)write(*,*)""

      deallocate(RNRM,W,Q,R,BASIS)

      return
    end function lanczos_hlzpack


    subroutine lanczos_gs(egs,vect,itermax)
      integer                            :: i,j,ierr,iter
      real(8),dimension(:)               :: Vect
      real(8)                            :: Egs
      real(8),dimension(size(Vect))      :: Vin,Vout
      integer                            :: Itermax,Nlanc,Ndim
      real(8),dimension(Itermax)         :: alfa,beta
      real(8),dimension(:,:),allocatable :: Z
      real(8),dimension(:),allocatable   :: diag,subdiag
      real(8)                            :: a_,b_,xnorm
      Ndim = size(vect)     
      Nlanc= min(Ndim,Itermax)
      !Copy arrays:
      vin = vect ; vout= zero
      !Perform lanczos, tri-diagonalize the Hamiltonian on a sufficiently large basis
      !Nlanc=Nitermax on input, # of actual iterations on output
      alfa=zero ; beta=zero
      call lanczos_simple(vin,vout,alfa,beta,Nlanc)

      !Build Z, tridiagonal matrix. Diagonalize and get ground-state energy.
      allocate(Z(Nlanc,Nlanc))
      allocate(diag(Nlanc),subdiag(Nlanc))
      Z=zero ; forall(i=1:Nlanc)Z(i,i)=one
      diag    = alfa(1:Nlanc)
      subdiag = beta(2:Nlanc)
      call tql2(Nlanc,diag,subdiag,Z,ierr)
      egs=diag(1)

      !Here one should put a check to improve convergence of the ground state, see Massimo.
      vin =vect
      vout=zero
      vect=zero
      do iter=1,nlanc
         call lanczos_one_step(vin,vout,alfa(iter),beta(iter),iter)
         vect = vect + vin*Z(iter,1)
      end do
      xnorm=sqrt(dot_product(vect,vect))
      vect=vect/xnorm
    end subroutine lanczos_gs

  END MODULE HLZPACK_LANCZOS
