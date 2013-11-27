!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE FUNCTIONS
  USE TOOLS
  USE MATRIX
  implicit none
  integer                :: iloop,Le
  logical                :: converged
  real(8)                :: wband
  !Bath:
  integer                :: Nb(2)
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  character(len=16)      :: input_file
#ifdef _MPI
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif


  !parse additional variables && read input file
  call parse_cmd_variable(wband,"wband",default=1.d0)
  call parse_cmd_variable(Le,"LE",default=1000)
  call parse_cmd_variable(input_file,'input',default='inputED.in')
  call read_input(trim(input_file))


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
     call chi2_fitgf(delta,bath,ispin=1)

     !Check convergence (if required change chemical potential)
     if(mpiID==0)then
        converged = check_convergence(delta(1,1,:),eps_error,nsuccess,nloop)
        if(nread/=0.d0)call search_mu(nimp(1),converged)
     endif
#ifdef _MPI
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
     call end_loop
  enddo

#ifdef _MPI
  call MPI_FINALIZE(mpiERR)
#endif


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                                 :: i,j,iorb,jorb,ie
    complex(8)                              :: iw
    complex(8),dimension(Norb,Norb)         :: zeta,fg,He
    complex(8),allocatable,dimension(:,:,:) :: gloc
    real(8)                                 :: wm(NL),wr(Nw)
    real(8)                                 :: de,e(Le),dos_wt(Le)

    wm = pi/beta*real(2*arange(1,NL)-1,8)
    wr = linspace(wini,wfin,Nw)

    de=2.d0*wband/real(Le,8)
    do ie=1,Le
       e(ie) = -wband + dfloat(ie-1)*de
       dos_wt(ie)=dens_bethe(e(ie),wband)*de
    enddo

    delta=zero

    allocate(gloc(Norb,Norb,NL))
    do i=1,NL
       zeta=zero
       do iorb=1,Norb
          zeta(iorb,iorb) = (xi*wm(i) + xmu)
          do jorb=1,Norb
             zeta(iorb,jorb) = zeta(iorb,jorb) - Hloc(iorb,jorb) - impSmats(1,iorb,jorb,i)
          enddo
       enddo
       !
       fg=zero
       do ie=1,Le
          fg = fg + invert_ge(e(ie),zeta)*dos_wt(ie)
       enddo
       gloc(:,:,i) = fg
       !
       call matrix_inverse(fg)
       delta(:,:,i)=zeta - fg 
    enddo
    !print
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,gloc(iorb,jorb,:))
          call splot("Delta_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_iw.ed",wm,delta(iorb,jorb,:))
       enddo
    enddo
    deallocate(gloc)


    allocate(gloc(Norb,Norb,Nw))
    do i=1,Nw
       zeta=zero
       do iorb=1,Norb
          zeta(iorb,iorb) = (cmplx(wr(i),eps) + xmu) 
          do jorb=1,Norb
             zeta(iorb,jorb) = zeta(iorb,jorb) - Hloc(iorb,jorb) - impSreal(1,iorb,jorb,i)
          enddo
       enddo
       !
       fg=zero
       do ie=1,Le
          fg = fg + invert_ge(e(ie),zeta)*dos_wt(ie)
       enddo
       gloc(:,:,i) = fg
    enddo
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("Gloc_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,gloc(iorb,jorb,:))
          call splot("DOS_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_realw.ed",wr,-dimag(gloc(iorb,jorb,:))/pi)
       enddo
    enddo
    deallocate(gloc)
  end subroutine get_delta_bethe
  !+----------------------------------------+

  function invert_ge(ek,zeta) result(invg)
    real(8) :: ek
    complex(8),dimension(Norb,Norb) :: zeta
    complex(8),dimension(Norb,Norb) :: He,invg
    integer :: i,j
    He=zero
    He(1,1)=1.d-9*ek
    He(2,2)=ek
    invg=zeta-He
    call matrix_inverse(invg)
  end function invert_ge

  include 'search_mu.f90'

end program lancED



