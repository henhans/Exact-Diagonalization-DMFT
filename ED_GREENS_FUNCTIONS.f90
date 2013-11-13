!###################################################################
!PURPOSE  : Build the impurity Green's function using spectral sum 
!NOTE: in the MPI implementation we may require all the nodes to 
!evaluate the GF, this is safer, simpler and works for both Lanc &
!Ed. For Lanc we can indeed assign the contribution from each state 
!to different node and accumulate the result at the end.
!AUTHORS  : Adriano Amaricci
!###################################################################
MODULE ED_GREENS_FUNCTIONS
  USE TIMER
  USE IOTOOLS, only: free_unit,reg
  USE TOOLS,   only: arange,linspace
  USE MATRIX,  only: matrix_inverse
  USE PLAIN_LANCZOS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  USE ED_HAMILTONIAN
  !
  implicit none
  private 

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable          :: wm,tau,wr,vm

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer              :: state_vec
  complex(8),dimension(:),pointer           :: state_cvec
  real(8)                                   :: state_e

  !Impurity GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:) :: impGmats
  complex(8),allocatable,dimension(:,:,:,:) :: impGreal


  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)        :: Chitau
  complex(8),allocatable,dimension(:,:)     :: Chiw,Chiiw

  public                                    :: lanc_ed_getgf
  public                                    :: full_ed_getgf
  public                                    :: lanc_ed_getchi
  public                                    :: full_ed_getchi


contains


  !                    LANCZOS DIAGONALIZATION 
  !+------------------------------------------------------------------+
  include 'ed_lanc_gf.f90'


  !                    LANC SUSCPTIBILITY
  !+------------------------------------------------------------------+
  include 'ed_lanc_chi.f90'

  !                    FULL DIAGONALIZATION
  !+------------------------------------------------------------------+
  include 'ed_full_gf.f90'


  !                    FULL SUSCEPTIBILITY
  !+------------------------------------------------------------------+
  include 'ed_full_chi.f90'


  !                    COMPUTATIONAL ROUTINES
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    allocate(wm(NL))
    wm     = pi/beta*real(2*arange(1,NL)-1,8)
    allocate(vm(0:NL))
    do i=0,NL
       vm(i) = pi/beta*2.d0*real(i,8)
    enddo
    allocate(wr(Nw))
    wr     = linspace(wini,wfin,Nw)
    allocate(tau(0:Ltau))
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_gf
    integer                                  :: i,j,ispin,unit(6),iorb,jorb
    complex(8)                               :: iw
    complex(8),dimension(Nspin,Norb,Norb,NL) :: impG0iw
    complex(8),dimension(Nspin,Norb,Norb,Nw) :: impG0wr
    complex(8),dimension(Norb,Norb)          :: invGimp,impG0
    character(len=20)                        :: suffix
    !
    write(LOGfile,"(A)")"Printing the impurity GF"
    impSmats = zero
    impSreal = zero
    !
    select case(bath_type)
    case default
       do ispin=1,Nspin
          do iorb=1,Norb
             do i=1,NL
                iw=xi*wm(i)
                impG0iw(ispin,iorb,iorb,i) = iw+xmu-eloc(iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
                impSmats(ispin,iorb,iorb,i)= impG0iw(ispin,iorb,iorb,i) - one/impGmats(ispin,iorb,iorb,i)
                impG0iw(ispin,iorb,iorb,i) = delta_bath(ispin,iorb,iw,dmft_bath)
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                impG0wr(ispin,iorb,iorb,i) = wr(i)+xmu-eloc(iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
                impSreal(ispin,iorb,iorb,i)= impG0wr(ispin,iorb,iorb,i) - one/impGreal(ispin,iorb,iorb,i)
                impG0wr(ispin,iorb,iorb,i) = delta_bath(ispin,iorb,iw,dmft_bath)
             enddo
          enddo
       enddo
       !Print the impurity functions:
       do iorb=1,Norb
          suffix="_orb"//reg(txtfy(iorb))
          call open_units(reg(suffix))
          do i=1,NL
             write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impGmats(ispin,iorb,iorb,i)),dreal(impGmats(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(3),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impSmats(ispin,iorb,iorb,i)),dreal(impSmats(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                  (dimag(impG0iw(ispin,iorb,iorb,i)),dreal(impG0iw(ispin,iorb,iorb,i)),ispin=1,Nspin)
          enddo
          do i=1,Nw
             write(unit(2),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impGreal(ispin,iorb,iorb,i)),dreal(impGreal(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impSreal(ispin,iorb,iorb,i)),dreal(impSreal(ispin,iorb,iorb,i)),ispin=1,Nspin)
             write(unit(6),"(F26.15,6(F26.15))")wr(i),&
                  (dimag(impG0wr(ispin,iorb,iorb,i)),dreal(impG0wr(ispin,iorb,iorb,i)),ispin=1,Nspin)
          enddo
          call close_units
       enddo



    case ('hybrid')
       do ispin=1,Nspin
          do iorb=1,Norb
             !Get WF diagonals:
             do i=1,NL
                iw=xi*wm(i)
                impG0iw(ispin,iorb,iorb,i)= iw + xmu-eloc(iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
             enddo
             do i=1,Nw
                iw=cmplx(wr(i),eps)
                !impG0wr(ispin,iorb,iorb,i)= wr(i)+xmu-eloc(iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
                impG0wr(ispin,iorb,iorb,i)= iw+xmu-eloc(iorb)-delta_bath(ispin,iorb,iorb,iw,dmft_bath)
             enddo
          enddo
          !Get WF off-diagonals
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                do i=1,NL
                   iw=xi*wm(i)
                   impG0iw(ispin,iorb,jorb,i)= -delta_bath(ispin,iorb,jorb,iw,dmft_bath)
                   impG0iw(ispin,jorb,iorb,i)= -delta_bath(ispin,jorb,iorb,iw,dmft_bath)
                enddo
                do i=1,Nw
                   iw=cmplx(wr(i),eps)
                   impG0wr(ispin,iorb,jorb,i)= -delta_bath(ispin,iorb,jorb,iw,dmft_bath)
                   impG0wr(ispin,jorb,iorb,i)= -delta_bath(ispin,jorb,iorb,iw,dmft_bath)
                enddo
             enddo
          enddo
       enddo

       !get Sigma and WF by matrix inversions:
       do ispin=1,Nspin
          do i=1,NL
             invGimp = impGmats(ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSmats(ispin,:,:,i) = impG0iw(ispin,:,:,i) - invGimp
          enddo
          do i=1,Nw
             invGimp = impGreal(ispin,:,:,i)
             call matrix_inverse(invGimp)
             impSreal(ispin,:,:,i) = impG0wr(ispin,:,:,i) - invGimp
          enddo
          do i=1,NL
             impG0=impG0iw(ispin,:,:,i)
             call matrix_inverse(impG0)
             impG0iw(ispin,:,:,i)=impG0
          enddo
          do i=1,Nw
             impG0=impG0wr(ispin,:,:,i)
             call matrix_inverse(impG0)
             impG0wr(ispin,:,:,i)=impG0
          enddo
       enddo
       !
       !Print the impurity functions:
       do iorb=1,Norb
          do jorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))
             call open_units(reg(suffix))
             do i=1,NL
                write(unit(1),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impGmats(ispin,iorb,jorb,i)),dreal(impGmats(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(3),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impSmats(ispin,iorb,jorb,i)),dreal(impSmats(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(5),"(F26.15,6(F26.15))")wm(i),&
                     (dimag(impG0iw(ispin,iorb,jorb,i)),dreal(impG0iw(ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             do i=1,Nw
                write(unit(2),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impGreal(ispin,iorb,jorb,i)),dreal(impGreal(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(4),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impSreal(ispin,iorb,jorb,i)),dreal(impSreal(ispin,iorb,jorb,i)),ispin=1,Nspin)
                write(unit(6),"(F26.15,6(F26.15))")wr(i),&
                     (dimag(impG0wr(ispin,iorb,jorb,i)),dreal(impG0wr(ispin,iorb,jorb,i)),ispin=1,Nspin)
             enddo
             call close_units()
          enddo
       enddo
    end select

  contains

    subroutine open_units(string)
      character(len=*) :: string
      unit(1)=free_unit()
      open(unit(1),file="impG"//string//"_iw.ed")
      unit(2)=free_unit()
      open(unit(2),file="impG"//string//"_realw.ed")
      unit(3)=free_unit()
      open(unit(3),file="impSigma"//string//"_iw.ed")
      unit(4)=free_unit()
      open(unit(4),file="impSigma"//string//"_realw.ed")
      unit(5)=free_unit()
      open(unit(5),file="impDelta"//string//"_iw.ed")
      unit(6)=free_unit()
      open(unit(6),file="impDelta"//string//"_realw.ed")
    end subroutine open_units

    subroutine close_units()
      do i=1,6
         close(unit(i))
      enddo
    end subroutine close_units

  end subroutine print_imp_gf




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_chi
    integer                               :: i,j,iorb
    integer                               :: unit(3)
    do iorb=1,Norb
       unit(1)=free_unit()
       open(unit(1),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_tau.ed")
       unit(2)=free_unit()
       open(unit(2),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_realw.ed")
       unit(3)=free_unit()
       open(unit(3),file=reg(CHIfile)//"_orb"//reg(txtfy(iorb))//"_iw.ed")
       do i=0,Ltau/2
          write(unit(1),*)tau(i),chitau(iorb,i)
       enddo
       do i=1,Nw
          if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(iorb,i)),dreal(chiw(iorb,i))
       enddo
       do i=0,NL
          write(unit(3),*)vm(i),dimag(chiiw(iorb,i)),dreal(chiiw(iorb,i))
       enddo
       close(unit(1))
       close(unit(2))
       close(unit(3))
    enddo
  end subroutine print_imp_chi


end MODULE ED_GREENS_FUNCTIONS
