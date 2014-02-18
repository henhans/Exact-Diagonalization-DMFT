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
  USE IOTOOLS, only: free_unit,reg,free_units
  USE ARRAYS,   only: arange,linspace
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
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm

  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer             :: state_cvec
  real(8)                                     :: state_e

  !Impurity GF
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:) :: impGmats,impFmats
  complex(8),allocatable,dimension(:,:,:,:,:) :: impGreal,impFreal
  complex(8),allocatable,dimension(:,:) :: Gaux_mats,Gaux_real

  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)          :: Chitau
  complex(8),allocatable,dimension(:,:)       :: Chiw,Chiiw

  public                                      :: lanc_ed_getgf
  public                                      :: full_ed_getgf
  public                                      :: lanc_ed_getchi
  public                                      :: full_ed_getchi

contains


  subroutine lanc_ed_getgf()
    if(.not.ed_supercond)then
       call lanc_ed_getgf_normal()
    else
       call lanc_ed_getgf_superc()
    end if
  end subroutine lanc_ed_getgf



  !                    LANCZOS DIAGONALIZATION 
  !+------------------------------------------------------------------+
  include 'ed_lanc_gf_normal.f90'
  include 'ed_lanc_gf_superc.f90'
  

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
    allocate(vm(0:NL))          !bosonic frequencies
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
    integer                                        :: i,j,ispin,unit(6),iorb,jorb
    complex(8)                                     :: iw,fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,NL) :: impG0iw
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nw) :: impG0wr
    complex(8),dimension(Norb,Norb)                :: invGimp,impG0
    character(len=20)                              :: suffix
    !
    write(LOGfile,"(A)")"Printing the impurity GF"
    impSmats = zero
    impSreal = zero
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       !                        !this is ensured by the special *per impurity" bath structure
       !                        !no intra-orbital hoopings
       include 'ed_print_gf_normal.f90'

    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       !                        !intra-orbital hopping allow for mixed _ab GF
       include 'ed_print_gf_hybrid.f90'

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
      open(unit(5),file="impG0"//string//"_iw.ed")
      unit(6)=free_unit()
      open(unit(6),file="impG0"//string//"_realw.ed")
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
  subroutine print_imp_gf_sc
    integer                                        :: i,j,ispin,unit(12),iorb,jorb
    complex(8)                                     :: iw
    complex(8),allocatable,dimension(:)            :: det
    complex(8),allocatable,dimension(:,:)          :: fg0,fg,sigma
    complex(8),dimension(Nspin,Nspin,Norb,Norb,NL) :: impG0mats,impF0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nw) :: impG0real,impF0real
    character(len=20)                              :: suffix
    !
    write(LOGfile,"(A)")"Printing the impurity GF"
    impSmats = zero
    impSreal = zero
    impSAmats = zero
    impSAreal = zero
    !
    !Diagonal in both spin and orbital
    !this is ensured by the special *per impurity" bath structure
    !no intra-orbital hoopings
    !THIS IS SUPERCONDUCTING CASE

    allocate(fg0(2,NL),fg(2,NL),sigma(2,NL),det(NL))
    do ispin=1,Nspin
       do iorb=1,Norb
          det     =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
          fg(1,:) =  conjg(impGmats(ispin,ispin,iorb,iorb,:))/det
          fg(2,:) =  impFmats(ispin,ispin,iorb,iorb,:)/det
          do i=1,NL
             iw = xi*wm(i)
             fg0(1,i) = iw+xmu-hloc(ispin,ispin,iorb,iorb)-delta_bath(ispin,iorb,iw,dmft_bath)
             fg0(2,i) = -fdelta_bath(ispin,iorb,iw,dmft_bath)
          enddo
          impSmats(ispin,ispin,iorb,iorb,:)= fg0(1,:) - fg(1,:)
          impSAmats(ispin,ispin,iorb,iorb,:)= fg0(2,:) - fg(2,:)
          det     =  abs(fg0(1,:))**2 + (fg0(2,:))**2
          impG0mats(ispin,ispin,iorb,iorb,:) = conjg(fg0(1,:))/det
          impF0mats(ispin,ispin,iorb,iorb,:) = fg0(2,:)/det
       enddo
    enddo
    deallocate(fg0,fg,sigma,det)

    
    allocate(fg0(2,Nw),fg(2,Nw),sigma(2,Nw),det(Nw))
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Nw
             iw=cmplx(wr(i),eps)
             det(i)  = impGreal(ispin,ispin,iorb,iorb,i)*conjg(impGreal(ispin,ispin,iorb,iorb,Nw+1-i)) + &
                  impFreal(ispin,ispin,iorb,iorb,i)*conjg(impFreal(ispin,ispin,iorb,iorb,Nw+1-i))
             fg(1,i) =  conjg(impGreal(ispin,ispin,iorb,iorb,Nw+1-i))/det(i)
             fg(2,i) =  conjg(impFreal(ispin,ispin,iorb,iorb,Nw+1-i))/det(i)
             fg0(1,i) = iw+xmu-hloc(ispin,ispin,iorb,iorb)-delta_bath(ispin,iorb,wr(i),eps,dmft_bath)
             fg0(2,i) = -fdelta_bath(ispin,iorb,wr(i),eps,dmft_bath)
          enddo
          impSreal(ispin,ispin,iorb,iorb,:)= fg0(1,:) - fg(1,:)
          impSAreal(ispin,ispin,iorb,iorb,:)= fg0(2,:) - fg(2,:)
          ! surely wrong...
          do i=1,Nw          
             det(i)     =  fg0(1,i)*conjg(fg0(2,Nw+1-i)) + fg0(2,i)*conjg(fg0(2,Nw+1-i))
             impG0real(ispin,ispin,iorb,iorb,i) = conjg(fg0(1,Nw+1-i))/det(i)
             impF0real(ispin,ispin,iorb,iorb,i) = conjg(fg0(2,Nw+1-i))/det(i)
          enddo
       enddo
    enddo
    deallocate(fg0,fg,sigma,det)

    !
    do iorb=1,Norb
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))
       call open_units(reg(suffix))
       do i=1,NL
          write(unit(1),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impGmats(ispin,ispin,iorb,iorb,i)),dreal(impGmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(2),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impFmats(ispin,ispin,iorb,iorb,i)),dreal(impFmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(3),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impSmats(ispin,ispin,iorb,iorb,i)),dreal(impSmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(4),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impSAmats(ispin,ispin,iorb,iorb,i)),dreal(impSAmats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(5),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impG0mats(ispin,ispin,iorb,iorb,i)),dreal(impG0mats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(6),"(F26.15,6(F26.15))")wm(i),&
               (dimag(impF0mats(ispin,ispin,iorb,iorb,i)),dreal(impF0mats(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
       enddo
       do i=1,Nw
          write(unit(7),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impGreal(ispin,ispin,iorb,iorb,i)),dreal(impGreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(8),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impFreal(ispin,ispin,iorb,iorb,i)),dreal(impFreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(9),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impSreal(ispin,ispin,iorb,iorb,i)),dreal(impSreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(10),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impSAreal(ispin,ispin,iorb,iorb,i)),dreal(impSAreal(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(11),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impG0real(ispin,ispin,iorb,iorb,i)),dreal(impG0real(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
          write(unit(12),"(F26.15,6(F26.15))")wr(i),&
               (dimag(impF0real(ispin,ispin,iorb,iorb,i)),dreal(impF0real(ispin,ispin,iorb,iorb,i)),ispin=1,Nspin)
       enddo
       call close_units
    enddo



  contains

    subroutine open_units(string)
      character(len=*) :: string
      unit=free_units(12)
      open(unit(1),file="impG"//string//"_iw.ed")
      open(unit(2),file="impF"//string//"_iw.ed")
      open(unit(3),file="impSigma"//string//"_iw.ed")
      open(unit(4),file="impSelf"//string//"_iw.ed")
      open(unit(5),file="impG0"//string//"_iw.ed")
      open(unit(6),file="impF0"//string//"_iw.ed")
      !
      open(unit(7),file="impG"//string//"_realw.ed")
      open(unit(8),file="impF"//string//"_realw.ed")
      open(unit(9),file="impSigma"//string//"_realw.ed")
      open(unit(10),file="impSelf"//string//"_realw.ed")
      open(unit(11),file="impG0"//string//"_realw.ed")
      open(unit(12),file="impF0"//string//"_realw.ed")
    end subroutine open_units

    subroutine close_units()
      do i=1,size(unit)
         close(unit(i))
      enddo
    end subroutine close_units

  end subroutine print_imp_gf_sc










  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_imp_chi
    integer                               :: i,j,iorb
    integer                               :: unit(3)
    write(LOGfile,"(A)")"Printing the spin Chi:"
    do iorb=1,Norb
       unit(1)=free_unit()
       open(unit(1),file="Chi_orb"//reg(txtfy(iorb))//"_tau.ed")
       unit(2)=free_unit()
       open(unit(2),file="Chi_orb"//reg(txtfy(iorb))//"_realw.ed")
       unit(3)=free_unit()
       open(unit(3),file="Chi_orb"//reg(txtfy(iorb))//"_iw.ed")
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
    print*,""
  end subroutine print_imp_chi


end MODULE ED_GREENS_FUNCTIONS
