!###################################################################
!PURPOSE  : Complete ED solution of DMFT equations 
!AUTHORS  : A. Amaricci
!###################################################################
program FULLED_AF_HM_BETHE
  USE DMFT_FULLED
  implicit none
  logical             :: converged
  real(8),allocatable :: wt(:),epsik(:)
  integer             :: Lk

  !Read the input file
  call read_input("inputED.in")

  !Create the Bethe lattice DOS (used to build local GF)
  D=2.d0*ts                     !half-bandwidth (set ts=0.5 in inputfile)
  Lk=Nx**2                      !number of energy points in the DOS
  allocate(wt(Lk),epsik(Lk))    !wt=dos(e), epsik=e\in[-D,D]
  call bethe_lattice(wt,epsik,Lk,D)

  !Setup the ED solver on first call with statue =-2
  !Bath is generated or read, hilbers space is created w/ all pointers
  call ed_solver(status=-2)

  !Start the DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM
     call ed_solver()

     !Get the Weiss field/Delta function to be fitted (user defined)
     !This routine MUST define a suitable function call delta(Nchannel,:)
     call get_delta_bethe_af_hm

     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(1,:),epsiup,vup) !UP
     call chi2_fitgf(delta(2,:),epsidw,vdw) !DOWN

     !Check convergence of the calculation
     converged = check_convergence(delta(:,:),eps_error,nsuccess,nloop)

     !If fixed density is required change the chemical potential (a simple but working method).
     if(nread/=0.d0)call search_mu(nimp1,converged)
     if(iloop>nloop)converged=.true.

     call end_loop
  enddo

  !finalize calculation
  call ed_solver(status=-1)


contains


  !+----------------------------------------+
  subroutine get_delta_bethe_af_hm
    integer                         :: i,j,ik
    real(8)                         :: w,nup,ndw,mu(Nspin)
    complex(8)                      :: iw,zita(Nspin)
    real(8),dimension(Nspin,0:NL)   :: gtau
    complex(8),dimension(Nspin,NL)  :: g0and,self,gloc
    complex(8),dimension(Nspin,Nw)  :: rg0and,rself,rgloc
    delta=zero
    mu(1)=xmu-heff
    mu(2)=xmu+heff
    do i=1,NL
       iw=xi*wm(i)
       !Get the Anderson non-interacting GF: \cal{G}^{And}^{-1}
       g0and(1,i) = iw + mu(1) -ed0 - delta_and(iw,epsiup,vup)
       g0and(2,i) = iw + mu(2) -ed0 - delta_and(iw,epsidw,vdw)
       !Get Sigma from Dyson Equation (Giw is the output of ed_solver)
       Self(:,i)  = g0and(:,i) - one/Giw(:,i)
       !Get the local GF from energy integral w/ respect to Bethe DOS
       !Check on the RMP96 eq.97
       zita       = iw + mu(:) - Self(:,i)
       gloc(1,i)  = sum(wt(:)*zita(2)/(zita(1)*zita(2)-epsik(:)**2))
       gloc(2,i)  = sum(wt(:)*zita(1)/(zita(1)*zita(2)-epsik(:)**2))
    enddo
    !
    !Get the Delta (or Hybridization) function
    !Delta_\sigma = D^2/4.0*G_{-\sigma}
    delta(1,:) =  Gloc(2,:)*D**2/4.d0
    delta(2,:) =  Gloc(1,:)*D**2/4.d0
    !
    do i=1,Nw
       iw=cmplx(wr(i),eps)
       !Get the Anderson non-interacting GF: \cal{G}^{And}^{-1}
       rG0and(1,i) = wr(i) + mu(1) - ed0 - delta_and(wr(i)+zero,epsiup,vup)
       rG0and(2,i) = wr(i) + mu(2) - ed0 - delta_and(wr(i)+zero,epsidw,vdw)
       !Get Sigma from Dyson Equation (Giw is the output of ed_solver)
       rSelf(:,i)  = rG0and(:,i) - one/Gwr(:,i)
       !Get the local GF from energy integral w/ respect to Bethe DOS
       zita      = iw + mu(:) - rSelf(:,i)
       rgloc(1,i)  = sum(wt(:)*zita(2)/(zita(1)*zita(2)-epsik(:)**2))
       rgloc(2,i)  = sum(wt(:)*zita(1)/(zita(1)*zita(2)-epsik(:)**2))
    enddo
    call splot("Delta_iw.ed",wm,delta(1,:))
    call splot("Delta_iw.ed",wm,delta(2,:),append=.true.)
    call splot("G_iw.ed",wm,Giw(1,:))
    call splot("G_iw.ed",wm,gloc(2,:),append=.true.)
    call splot("Sigma_iw.ed",wm,self(1,:))
    call splot("Sigma_iw.ed",wm,self(2,:),append=.true.)
    ! 
    call splot("G_realw.ed",wr,rGloc(1,:))
    call splot("G_realw.ed",wr,rGloc(2,:),append=.true.)
    call splot("Sigma_realw.ed",wr,rSelf(1,:))
    call splot("Sigma_realw.ed",wr,rSelf(2,:),append=.true.)

  end subroutine get_delta_bethe_af_hm
  !+----------------------------------------+

end program FULLED_AF_HM_BETHE



