!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE FUNCTIONS
  USE TOOLS
  implicit none
  integer                :: iloop,Nb(2)
  logical                :: converged
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:)

#ifdef _MPI
  call ed_init_mpi()
#endif

  call read_input("inputED.in")
  call parse_cmd_variable(wband,"wband",default=1.d0)

  !Allocate Weiss Field:
  allocate(delta(Norb,NL))

  !setup solver
  Nb=get_bath_size()
  allocate(bath(Nb(1),Nb(2)))
  call init_ed_solver(bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.OR.iloop>nloop)
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
        converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
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
    complex(8),dimension(NL)  :: gloc
    complex(8),dimension(Nw)  :: grloc
    real(8)                   :: wm(NL),wr(Nw)

    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)

    do iorb=1,Norb
       do i=1,NL
          iw = xi*wm(i)
          zita    = iw + xmu - impSmats(1,iorb,iorb,i)
          gloc(i) = gfbethe(wm(i),zita,Wband)
          if(cg_scheme=='weiss')then
             delta(iorb,i)= one/(one/gloc(i) + impSmats(1,iorb,iorb,i))
          else
             delta(iorb,i)= iw + xmu - impSmats(1,iorb,iorb,i) - one/gloc(i)
          endif
       enddo

       do i=1,Nw
          iw=cmplx(wr(i),eps)
          zita     = iw + xmu - impSreal(1,iorb,iorb,i)
          grloc(i) = gfbether(wr(i),zita,Wband)
       enddo
       call splot("Gloc_"//reg(txtfy(iorb))//"_iw.ed",wm,gloc)
       call splot("Gloc_"//reg(txtfy(iorb))//"_realw.ed",wr,grloc)
       call splot("DOS"//reg(txtfy(iorb))//".ed",wr,-dimag(grloc)/pi)
       call splot("Delta_"//reg(txtfy(iorb))//"_iw.ed",wm,delta(iorb,:))
    enddo
  end subroutine get_delta_bethe
  !+----------------------------------------+

  include 'search_mu.f90'
end program lancED



