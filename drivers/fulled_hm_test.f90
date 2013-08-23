!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
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
  complex(8),allocatable :: Delta(:,:)


  call read_input("inputED.in")
  call parse_cmd_variable(wband,"wband",default=1.d0)

  allocate(wm(NL),wr(Nw))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)

  !Allocate Weiss Field:
  allocate(delta(Norb,NL))

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
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nimp(1),converged)
     if(iloop>nloop)converged=.true.
     call end_loop
  enddo


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                   :: i,j,iorb
    complex(8)                :: iw,zita,g0and,g0loc
    complex(8),dimension(NL)  :: self,gloc
    complex(8),dimension(Nw)  :: selfr,grloc

    do iorb=1,Norb
       do i=1,NL
          iw = xi*wm(i)
          g0and   = iw + xmu - eloc(iorb) - delta_and(1,iorb,iw,bath)
          self(i) = g0and - one/impGmats(1,iorb,i)
          zita    = iw + xmu - self(i)
          gloc(i) = gfbethe(wm(i),zita,Wband)
          g0loc   = self(i) + one/gloc(i)
          delta(iorb,i)= iw + xmu - g0loc
       enddo

       do i=1,Nw
          iw=cmplx(wr(i),eps)
          g0and    = wr(i) + xmu - eloc(iorb) - delta_and(1,iorb,wr(i)+zero,bath)
          selfr(i) = g0and - one/impGreal(1,iorb,i)    
          zita     = iw + xmu - selfr(i)
          grloc(i) = gfbether(wr(i),zita,Wband)
       enddo
       call splot("Gloc_"//reg(txtfy(iorb))//"_iw.ed",wm,gloc)
       call splot("Sigma_"//reg(txtfy(iorb))//"_iw.ed",wm,self)
       call splot("Gloc_"//reg(txtfy(iorb))//"_realw.ed",wr,grloc)
       call splot("Sigma_"//reg(txtfy(iorb))//"_realw.ed",wr,selfr)
       call splot("DOS"//reg(txtfy(iorb))//".ed",wr,-dimag(grloc)/pi)
       call splot("Delta_"//reg(txtfy(iorb))//"_iw.ed",wm,delta(iorb,:))
    enddo
  end subroutine get_delta_bethe
  !+----------------------------------------+

end program lancED



