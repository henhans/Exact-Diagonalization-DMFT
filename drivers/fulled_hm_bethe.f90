!###################################################################
!PURPOSE  : Fulle ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program fullED
  USE DMFT_ED
  USE FUNCTIONS
  USE TOOLS
  implicit none
  integer :: iloop,Nb
  logical :: converged
  real(8),allocatable    :: wm(:),wr(:)
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:)
  !The local hybridization function:
<<<<<<< HEAD
  complex(8),allocatable :: Delta(:)
=======
  complex(8),allocatable :: Delta(:,:,:)
>>>>>>> devel_multi-orbital_nomix


  call read_input("inputED.in")
  call parse_cmd_variable(wband,"wband","D",default=1.d0)


  allocate(wm(NL),wr(Nw))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)

  !Allocate Weiss Field:
<<<<<<< HEAD
  allocate(delta(NL))
=======
  allocate(delta(Norb,Norb,NL))
>>>>>>> devel_multi-orbital_nomix

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_full_ed_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call full_ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe

     !Perform the SELF-CONSISTENCY by fitting the new bath
<<<<<<< HEAD
     call chi2_fitgf(delta(:),bath,ichan=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nimp,converged)
=======
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nimp(1),converged)
>>>>>>> devel_multi-orbital_nomix
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                   :: i,j
    complex(8)                :: iw,zita,g0and,g0loc
    complex(8),dimension(NL)  :: self,gloc
    complex(8),dimension(Nw)  :: selfr,grloc

    do i=1,NL
       iw = xi*wm(i)
<<<<<<< HEAD
       g0and   = iw + xmu - delta_and(iw,bath,1)
       self(i) = g0and - one/Giw(1,i)
       zita    = iw + xmu - self(i)
       gloc(i) = gfbethe(wm(i),zita,Wband)
       g0loc   = self(i) + one/gloc(i)
       delta(i)= iw + xmu - g0loc
=======
       g0and   = iw + xmu - delta_and(iw,bath,1,1,1)
       self(i) = g0and - one/Giw(1,1,1,i)
       zita    = iw + xmu - self(i)
       gloc(i) = gfbethe(wm(i),zita,Wband)
       g0loc   = self(i) + one/gloc(i)
       delta(1,1,i)= iw + xmu - g0loc
>>>>>>> devel_multi-orbital_nomix
    enddo

    do i=1,Nw
       iw=cmplx(wr(i),eps)
<<<<<<< HEAD
       g0and    = wr(i) + xmu - delta_and(wr(i)+zero,bath,1)
       selfr(i) = g0and - one/Gwr(1,i)    
=======
       g0and    = wr(i) + xmu - delta_and(wr(i)+zero,bath,1,1,1)
       selfr(i) = g0and - one/Gwr(1,1,1,i)    
>>>>>>> devel_multi-orbital_nomix
       zita     = iw + xmu - selfr(i)
       grloc(i) = gfbether(wr(i),zita,Wband)
    enddo
    call splot("G_iw.ed",wm,gloc)
    call splot("Sigma_iw.ed",wm,self)
    call splot("G_realw.ed",wr,grloc)
    call splot("Sigma_realw.ed",wr,selfr)
<<<<<<< HEAD
    call splot("Delta_iw.ed",wm,delta(:))
=======
    call splot("Delta_iw.ed",wm,delta(1,1,:))
>>>>>>> devel_multi-orbital_nomix
    call splot("DOS.ed",wr,-dimag(grloc)/pi)
    return    
  end subroutine get_delta_bethe
  !+----------------------------------------+

end program fullED



