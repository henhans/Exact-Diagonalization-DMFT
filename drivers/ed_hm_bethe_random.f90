!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE COMMON_VARS
  USE FUNCTIONS
  USE IOTOOLS
  USE TOOLS
  USE MATRIX
  USE ERROR
  USE ARRAYS
  USE FFTGF
  USE PARSE_INPUT
  implicit none
  integer                :: iloop,Nb(2),is
  logical                :: converged
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:,:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:),Gmats(:,:),Smats(:,:),Greal(:,:),Sreal(:,:)
  character(len=16)      :: finput,fhloc
  integer                :: nsample,idum
  character(len=5)       :: tmp_suffix
  real(8),allocatable    :: eloc(:)

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')

  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(nsample,"nsample",finput,default=20)
  call ed_read_input(trim(finput))
  Hloc=zero

  !Allocate Weiss Field:
  allocate(delta(Nsample,Norb,Norb,Lmats))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(nsample,Nb(1),Nb(2)))
  allocate(eloc(nsample))
  idum=123456
  call random_seed(idum)
  call random_number(eloc)
  eloc = (2.d0*eloc-1.d0)*wband/2.d0
  do is=1,nsample
     write(tmp_suffix,'(I4.4)') is
     ed_file_suffix="_site"//trim(tmp_suffix)
     call init_ed_solver(bath(is,:,:),bcenter=eloc(is))
  end do
  allocate(Smats(nsample,Lmats),Sreal(Nsample,Lreal))


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     do is=1,nsample
        !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
        write(tmp_suffix,'(I4.4)') is
        ed_file_suffix="_site"//trim(tmp_suffix)
        call ed_solver(bath(is,:,:)) 
        Smats(is,:) = impSmats(1,1,1,1,:)
        Sreal(is,:) = impSreal(1,1,1,1,:)
     enddo

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe

     !Perform the SELF-CONSISTENCY by fitting the new bath
     do is=1,nsample
        write(tmp_suffix,'(I4.4)') is
        ed_file_suffix="_site"//trim(tmp_suffix)
        call chi2_fitgf(delta(is,:,:,:),bath(is,:,:),ispin=1)
     enddo

     !Check convergence (if required change chemical potential)
     converged = check_convergence(sum(delta(:,1,1,:),1)/dble(nsample),dmft_error,nsuccess,nloop,reset=.false.)
     call end_loop
  enddo


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                     :: i,j,iorb
    complex(8)                  :: iw,zita,g0loc
    complex(8),dimension(nsample,Lmats) :: gloc
    complex(8),dimension(nsample,Lreal) :: grloc
    real(8)                     :: wm(Lmats),wr(Lreal)
    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    iorb=1
    do is=1,nsample
       do i=1,Lmats
          iw = xi*wm(i)
          zita    = iw + xmu - Smats(is,i)
          gloc(is,i) = gfbethe(wm(i),zita,Wband)
          if(cg_scheme=='weiss')then
             delta(is,1,1,i)= one/(one/gloc(is,i) + Smats(is,i))
          else
             delta(is,1,1,i)= iw + xmu - Smats(is,i) - one/gloc(is,i)
          endif
       enddo

       do i=1,Lreal
          iw=cmplx(wr(i),eps)
          zita     = iw + xmu - Sreal(is,i)
          grloc(is,i) = gfbether(wr(i),zita,Wband)
       enddo
    enddo
    call splot("LGloc_iw.ed",wm,gloc(:,:))
    call splot("aGloc_iw.ed",wm,sum(gloc(:,:),1)/nsample)
    call splot("LGloc_realw.ed",wr,-dimag(grloc(:,:))/pi,dreal(grloc(:,:)))
    call splot("aGloc_realw.ed",wr,-dimag(sum(grloc(:,:),1))/pi/nsample)
    call splot("Delta_iw.ed",wm(1:Lmats),Delta(1:Nsample,1,1,1:Lmats))


  end subroutine get_delta_bethe
  !+----------------------------------------+










end program lancED



