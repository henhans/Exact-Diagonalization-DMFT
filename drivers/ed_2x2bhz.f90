!This program solves the DMFT equations
!using complete ED at finite (high) T 
!reading the Hamiltonian model H(k)
!from a file with the standard
!wannier90 form.
!IMPORTANT: The ED solver uses HFmode=.true. 
!by default this means that chemical potential 
!SHOULD NOT include the U/2 shift.
program ed_lda1b
  USE DMFT_ED
  USE FFTGF
  USE TOOLS
  USE FUNCTIONS
  implicit none
  integer                :: i,ik,iorb,jorb,ispin,iloop,Lk,Nb
  integer                :: Norb_d,Lorb
  logical                :: converged
  real(8)                :: kx,ky,foo,ntotal,npimp
  real(8),allocatable    :: wm(:),wr(:)
  complex(8)             :: iw
  !variables for the model:
  character(len=32)      :: file
  integer                :: ntype
  real(8)                :: nobj
  logical                :: rerun,bool
  !Bath:
  real(8),allocatable    :: Bath(:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:)
  !Hamiltonian input:
  complex(8),allocatable :: Hk(:,:,:)
  real(8),allocatable    :: dos_wt(:),e0(:)
  logical                :: fbethe
  real(8)                :: wbath

  call read_input("inputED.in")
  !parse additional variables
  call parse_cmd_variable(file,"FILE",default="hkfile.in")
  call parse_cmd_variable(ntype,"NTYPE",default=0)
  call parse_cmd_variable(fbethe,"FBETHE",default=.false.)
  call parse_cmd_variable(wbath,"WBATH",default=1.d0)


  !Allocate:
  allocate(wm(NL),wr(Nw))
  wm = pi/beta*real(2*arange(1,NL)-1,8)
  wr = linspace(wini,wfin,Nw)

  !Read Hamiltoanian H(k)
  call read_hk(trim(file))

  !Allocate Weiss Field:
  allocate(delta(Nspin,Norb,NL))

  !Setup solver
  Nb=get_bath_size()
  allocate(bath(Nb))
  call init_ed_solver(bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver(bath) 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta


     !Fit the new bath, starting from the old bath + the supplied delta
     call chi2_fitgf(delta(1,1:Norb,1:NL),bath,ispin=1)

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,1,:),eps_error,nsuccess,nloop)
     !if(nread/=0.d0)call search_mu(nobj,converged)
     if(iloop>=nloop)converged=.true.
     call end_loop
  enddo


contains



  subroutine read_hk(file)
    character(len=*)       :: file
    integer                :: ik
    real(8)                :: de,e
    complex(8),allocatable :: fg(:,:,:)
    complex(8) :: w
    open(50,file=file,status='old')
    read(50,*)Lk,foo,Norb_d,foo,foo
    Lorb=Norb_d
    if(Lorb/=2)stop "Norb_d!=2"
    allocate(Hk(Lorb,Lorb,Lk))
    allocate(dos_wt(Lk))
    allocate(e0(Lorb))
    do ik=1,Lk
       read(50,"(3(F10.7,1x))")kx,ky,foo
       do iorb=1,Lorb
          read(50,"(10(2F10.7,1x))")(Hk(iorb,jorb,ik),jorb=1,Lorb)
       enddo
    enddo
    write(*,*)"# of k-points:",Lk
    write(*,*)"# of d-bands :",Norb_d
    dos_wt=1.d0/dfloat(Lk)
    write(*,*)"Evaluate non-interacting G(w):"
    allocate(fg(Nw,Lorb,Lorb))
    fg=zero
    e0=0.d0
    do ik=1,Lk
       do i=1,Nw
          w = cmplx(wr(i),eps,8)+xmu
          fg(i,:,:)=fg(i,:,:)+inverse_g0k(w,Hk(:,:,ik))*dos_wt(ik)
       enddo
       do i=1,Lorb
          e0(i)=e0(i)+Hk(i,i,ik)*dos_wt(ik)
       enddo
    enddo
    print*,e0
    call splot("DOS0.ed",wr,-dimag(fg(:,1,1))/pi,-dimag(fg(:,2,2))/pi)
    deallocate(fg)
  end subroutine read_hk



  function inverse_g0k(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    Complex(8)                  :: delta,ppi,vmix12,vmix21
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix12= -hk(1,2)
    vmix21= -hk(2,1)
    g0k(1,1) = one/(delta - vmix21*vmix12/ppi)
    g0k(2,2) = one/(ppi - vmix21*vmix12/delta)
    g0k(1,2) =-vmix21/(ppi*delta - vmix21*vmix12)
    g0k(2,1) =-vmix12/(ppi*delta - vmix21*vmix12)
  end function inverse_g0k

  function inverse_gk(zeta,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8),dimension(2)     :: zeta
    complex(8),dimension(2,2)   :: g0k
    Complex(8)                  :: delta,ppi,vmix12,vmix21
    g0k=zero
    delta = zeta(1) - hk(1,1)
    ppi   = zeta(2) - hk(2,2)
    vmix12= -hk(1,2)
    vmix21= -hk(2,1)
    g0k(1,1) = one/(delta - vmix21*vmix12/ppi)
    g0k(2,2) = one/(ppi - vmix21*vmix12/delta)
    g0k(1,2) =-vmix21/(ppi*delta - vmix21*vmix12)
    g0k(2,1) =-vmix12/(ppi*delta - vmix21*vmix12)
  end function inverse_gk




  subroutine get_delta
    integer                                 :: i,j,iorb,ispin
    complex(8)                              :: iw,zita(Lorb),fg(Lorb,Lorb)
    complex(8),dimension(:,:,:),allocatable :: gd,sd,g0d
    delta=zero

    ispin=1
    allocate(gd(1,Norb,NL),sd(1,Norb,NL),g0d(1,Norb,NL))
    print*,"Get Gloc_iw:"
    open(100,file="Delta_iw.ed")
    open(101,file="G_iw.ed")
    open(102,file="Sigma_iw.ed")
    do i=1,NL
       iw = xi*wm(i)
       do iorb=1,Norb
          g0d(ispin,iorb,i) = iw + xmu - delta_and(ispin,iorb,iw,bath)
          sd(ispin,iorb,i)  = g0d(ispin,iorb,i) - one/impGmats(ispin,iorb,i)
          zita(iorb)= iw + xmu - sd(ispin,iorb,i)
       enddo

       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gd(1,1,i)  = fg(1,1)
       gd(1,2,i)  = fg(2,2)
       !
       delta(1,1,i) = iw+xmu-one/gd(1,1,i)-sd(1,1,i)
       delta(1,2,i) = iw+xmu-one/gd(1,2,i)-sd(1,2,i)
       !
       write(100,"(5F25.12)")wm(i),dimag(delta(1,1,i)),dimag(delta(1,2,i)),dreal(delta(1,1,i)),dreal(delta(1,2,i))
       write(101,"(5F25.12)")wm(i),dimag(gd(1,1,i)),dimag(gd(1,2,i)),dreal(gd(1,1,i)),dreal(gd(1,2,i))
       write(102,"(5F25.12)")wm(i),dimag(sd(1,1,i)),dimag(sd(1,2,i)),dreal(sd(1,1,i)),dreal(sd(1,2,i))
    enddo
    deallocate(gd,sd,g0d)
    close(100);close(101);close(102)


    allocate(gd(1,Norb,Nw),sd(1,Norb,Nw),g0d(1,Norb,Nw))
    print*,"Get Gloc_realw:"
    open(100,file="DOS.ed")
    open(101,file="G_realw.ed")
    open(102,file="Sigma_realw.ed")
    do i=1,Nw
       iw=cmplx(wr(i),eps)
       do iorb=1,Norb
          g0d(ispin,iorb,i) = iw + xmu - delta_and(ispin,iorb,iw,bath)
          sd(ispin,iorb,i)  = g0d(ispin,iorb,i) - one/impGreal(ispin,iorb,i)
          zita(iorb)= iw + xmu - sd(iorb,iorb,i)
       enddo
       !
       fg=zero
       do ik=1,Lk
          fg=fg+inverse_gk(zita,Hk(:,:,ik))*dos_wt(ik)
       enddo
       gd(1,1,i)  = fg(1,1)
       gd(1,2,i)  = fg(2,2)

       write(100,"(3F25.12)")wr(i),-dimag(gd(1,1,i))/pi,-dimag(gd(1,2,i))/pi
       write(101,"(5F25.12)")wr(i),dimag(gd(1,1,i)),dimag(gd(1,2,i)),dreal(gd(1,1,i)),dreal(gd(1,2,i))
       write(102,"(5F25.12)")wr(i),dimag(sd(1,1,i)),dimag(sd(1,2,i)),dreal(sd(1,1,i)),dreal(sd(1,2,i))
    enddo
    deallocate(gd,sd,g0d)
    close(100);close(101);close(102)

    ntotal=sum(nimp)
    write(*,*)"ntot=",ntotal
  end subroutine get_delta




  function get_density_fromFFT(giw,beta) result(n)
    complex(8),dimension(:) :: giw
    real(8)                 :: gtau(0:size(giw))
    real(8)                 :: beta,n
    call fftgf_iw2tau(giw,gtau,beta)
    n = -2.d0*gtau(size(giw))
  end function get_density_fromFFT




end program ed_lda1b



