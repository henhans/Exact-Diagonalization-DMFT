!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE FUNCTIONS
  USE TOOLS
  USE INTEGRATE
  implicit none
  integer :: iloop,Nb,Ne,ie
  logical :: converged
  real(8),allocatable    :: wm(:),wr(:)
  real(8)                :: wband,ts,de
  !Bath:
  real(8),allocatable    :: Bath(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:)
  real(8),allocatable :: epsik(:),wt(:)

  call read_input("inputED.in")
  call parse_cmd_variable(ts,"TS",default=1.d0)
  call parse_cmd_variable(Ne,"NE",default=2000)

  allocate(wm(NL),wr(Nw))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)

  !Allocate Weiss Field:
  allocate(delta(Norb,NL))

  allocate(wt(Ne),epsik(Ne))
  wband=4.d0*ts
  epsik = linspace(-wband,wband,Ne,mesh=de)
  do ie=1,Ne
     wt(ie)=dens_2dsquare(epsik(ie),ts)
  enddo
  wt=wt/trapz(de,wt)
  call splot("DOS2d.qmc",epsik,wt)
  wt = wt*de

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_ed_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     !if(nread/=0.d0)call search_mu(nimp(1),converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains

  !+----------------------------------------+
  subroutine get_delta
    integer                   :: i,j,ie
    complex(8)                :: iw,zita,g0and,g0loc,gg
    complex(8),dimension(NL)  :: self,gloc
    complex(8),dimension(Nw)  :: selfr,grloc

    do i=1,NL
       iw = xi*wm(i)
       zita    = iw + xmu - impSmats(1,1,i)
       gloc(i) = zero
       do ie=1,Ne
          gloc(i)=gloc(i)+wt(ie)/(zita-epsik(ie))
       enddo
       delta(1,i)= iw + xmu - impSmats(1,1,i) - one/gloc(i)
    enddo

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       zita     = iw + xmu - impSreal(1,1,i)
       grloc(i) = zero
       do ie=1,Ne
          grloc(i)=grloc(i)+wt(ie)/(zita-epsik(ie))
       enddo
    enddo
    call splot("Gloc_iw.ed",wm,gloc)
    call splot("Gloc_realw.ed",wr,grloc)
    call splot("DOS.ed",wr,-dimag(grloc)/pi)
    call splot("Delta_iw.ed",wm,delta(1,:))

  end subroutine get_delta
  !+----------------------------------------+

end program lancED



