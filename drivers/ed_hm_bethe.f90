!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE FUNCTIONS
  USE TOOLS
  USE MATRIX
  USE ERROR
  USE ARRAYS
  USE FFTGF
  implicit none
  integer                :: iloop,Nb(2)
  logical                :: converged
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  character(len=16)      :: finput,fhloc

#ifdef _MPI
  call ed_init_mpi()
#endif

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_cmd_variable(fhloc,"FHLOC",default='inputHLOC.in')
  call parse_cmd_variable(wband,"wband",default=1.d0)
  !
  call ed_read_input(trim(finput),trim(fhloc))
  Hloc=zero

  !Allocate Weiss Field:
  allocate(delta(Norb,Norb,NL))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call chi2_fitgf(delta,bath,ispin=1,iverbose=.true.)

     !Check convergence (if required change chemical potential)
     if(mpiID==0)then
        converged = check_convergence(delta(1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
        if(nread/=0.d0)call search_chemical_potential(nimp(1),niter,converged)
     endif
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
     call end_loop
  enddo

#ifdef _MPI
  call ed_finalize_mpi()
#endif


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                   :: i,j,iorb
    complex(8)                :: iw,zita,g0loc
    complex(8),dimension(NL)  :: gloc,sigma
    complex(8),dimension(Nw)  :: grloc
    real(8)                   :: wm(NL),wr(Nw),tau(0:NL),C0,C1,n0
    real(8),dimension(0:NL)   :: sigt
    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)
    tau(0:) = linspace(0.d0,beta,NL+1)
    do iorb=1,Norb
       do i=1,NL
          iw = xi*wm(i)
          zita    = iw + xmu - impSmats(1,1,iorb,iorb,i)
          gloc(i) = gfbethe(wm(i),zita,Wband)
          if(cg_scheme=='weiss')then
             delta(iorb,iorb,i)= one/(one/gloc(i) + impSmats(1,1,iorb,iorb,i))
          else
             delta(iorb,iorb,i)= iw + xmu - impSmats(1,1,iorb,iorb,i) - one/gloc(i)
          endif
       enddo

       do i=1,Nw
          iw=cmplx(wr(i),eps)
          zita     = iw + xmu - impSreal(1,1,iorb,iorb,i)
          grloc(i) = gfbether(wr(i),zita,Wband)
       enddo
       call splot("Gloc_"//reg(txtfy(iorb))//"_iw.ed",wm,gloc)
       call splot("Gloc_"//reg(txtfy(iorb))//"_realw.ed",wr,grloc)
       call splot("DOS"//reg(txtfy(iorb))//".ed",wr,-dimag(grloc)/pi)
       call splot("Delta_"//reg(txtfy(iorb))//"_iw.ed",wm,delta(iorb,iorb,:))

       n0=nimp(1)/2.d0
       C0=Uloc(1)*(n0-0.5d0)
       C1=Uloc(1)**2*n0*(1.d0-n0)
       print*,n0,C0,C1
       sigma = impSmats(1,1,iorb,iorb,:)  - C1/(xi*wm) - C0
       call fftgf_iw2tau(sigma,sigt(0:),beta,notail=.true.)
       sigt=sigt-C1*0.5d0
       call splot("Sigma_"//reg(txtfy(iorb))//"_tau.ed",tau,sigt)

       call fftgf_tau2iw(sigt(0:),sigma,beta)

       call splot("Sigma_"//reg(txtfy(iorb))//"_iw.ed",wm,sigma)
    enddo



  end subroutine get_delta_bethe
  !+----------------------------------------+

end program lancED



