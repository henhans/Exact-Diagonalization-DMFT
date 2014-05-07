!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE CONSTANTS
  USE IOTOOLS
  USE FUNCTIONS
  USE TOOLS
  USE MATRIX
  USE ERROR
  USE ARRAYS
  USE FFTGF
  implicit none
  integer                :: iloop,Nb(2),Lk
  logical                :: converged
  real(8)                :: wband,ts
  !Bath:
  real(8),allocatable    :: Bath(:,:)
  !The local hybridization function:
  complex(8),allocatable :: Delta(:,:,:,:)
  character(len=16)      :: finput,fhloc

  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')
  call parse_input_variable(wband,"wband",finput,default=1.d0,comment="Bethe Lattice bandwidth")
  call parse_input_variable(Lk,"Lk",finput,default=500,comment="Number of energy levels for Bethe DOS integration")
  !
  call ed_read_input(trim(finput))

  !Allocate Weiss Field:
  allocate(delta(2,Norb,Norb,Lmats))

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
     converged = check_convergence(delta(1,1,1,:)+delta(2,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
     if(nread/=0.d0)call search_chemical_potential(ed_dens(1),converged)
     call end_loop
  enddo

  call get_sc_internal_energy(Lmats)

contains

  !+----------------------------------------+
  subroutine get_delta_bethe
    integer                    :: i,j,iorb,ik
    complex(8)                 :: iw,zita,g0loc,cdet,zita1,zita2
    complex(8),dimension(Lmats)   :: zeta
    complex(8),dimension(2,Lmats) :: gloc,calG
    complex(8),dimension(2,Lreal) :: grloc
    real(8)                    :: wm(Lmats),wr(Lreal),tau(0:Lmats)
    real(8),dimension(Lk)       :: epsik,wt

    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr = linspace(wini,wfin,Lreal)
    call bethe_lattice(wt,epsik,Lk,1.d0)

    delta=zero
    do i=1,Lmats
       iw = xi*wm(i)
       zita    = iw + xmu - impSmats(1,1,1,1,i)
       gloc(:,i)=zero
       do ik=1,Lk
          cdet = abs(zita-epsik(ik))**2 + impSAmats(1,1,1,1,i)**2
          gloc(1,i)=gloc(1,i) + wt(ik)*(conjg(zita)-epsik(ik))/cdet
          gloc(2,i)=gloc(2,i) - wt(ik)*impSAmats(1,1,1,1,i)/cdet
       enddo
       if(cg_scheme=='weiss')then
          !Get G0^{-1} matrix components:
          cdet      =  abs(gloc(1,i))**2 + (gloc(2,i))**2
          calG(1,i) =  conjg(gloc(1,i))/cdet + impSmats(1,1,1,1,i)
          calG(2,i) =  gloc(2,i)/cdet        + impSAmats(1,1,1,1,i) 
          !Get Weiss field G0 components:
          cdet            =  abs(calG(1,i))**2 + (calG(2,i))**2
          delta(1,1,1,i)  =  conjg(calG(1,i))/cdet
          delta(2,1,1,i)  =  calG(2,i)/cdet
          write(200,*)wm(i),dimag(delta(1,1,1,i)),dreal(delta(2,1,1,i))
       else
          cdet            = abs(gloc(1,i))**2 + (gloc(2,i))**2
          delta(1,1,1,i)  = iw + xmu - impSmats(1,1,1,1,i)  - conjg(gloc(1,i))/cdet 
          delta(2,1,1,i)  =          - impSAmats(1,1,1,1,i) - gloc(2,i)/cdet 
       endif
    enddo
    !
    zeta(:) = cmplx(wr(:),eps,8) + xmu - impSreal(1,1,1,1,:)
    do i=1,Lreal
       zita1 = zeta(i)
       zita2 = conjg(zeta(Lreal+1-i))
       grloc(:,i) = zero
       do ik=1,Lk
          cdet = (zita1-epsik(ik))*(zita2-epsik(ik)) + impSAreal(1,1,1,1,i)*impSAreal(1,1,1,1,i)
          grloc(1,i) = grloc(1,i) + wt(ik)*(zita2-epsik(ik))/cdet
          grloc(2,i) = grloc(2,i) + wt(ik)*impSAreal(1,1,1,1,i)/cdet
       enddo
    enddo

    call splot("Gloc_iw.ed",wm,gloc(1,:))
    call splot("Floc_iw.ed",wm,gloc(2,:))
    call splot("Gloc_realw.ed",wr,grloc(1,:))
    call splot("Floc_realw.ed",wr,grloc(2,:))
    call splot("DOS.ed",wr,-dimag(grloc(1,:))/pi)
    call splot("Delta_iw.ed",wm,delta(1,1,1,:),delta(2,1,1,:))


  end subroutine get_delta_bethe
  !+----------------------------------------+




  subroutine get_sc_internal_energy(L)
    integer                       :: L
    real(8)                       :: wm(Lmats)
    complex(8)                    :: fg(2,L),sigma(2,L)
    real(8)                       :: matssum,fmatssum,checkP,checkdens,vertex,Dssum
    complex(8)                    :: iw,gkw,fkw,g0kw,f0kw
    real(8)                       :: Epot,Etot,Eint,kin,kinsim,Ds,docc
    real(8)                       :: Sigma_infty,S_infty,det,det_infty,csi,Ei,thermal_factor
    real(8)                       :: free(Lk),Ffree(Lk),n_k(Lk),n,delta,u,ts
    integer                       :: i,j,iorb,ik
    complex(8)                    :: zita,g0loc,cdet,zita1,zita2
    complex(8),dimension(Lmats)   :: zeta
    real(8),dimension(Lk)         :: epsik,wt

    wm = pi/beta*real(2*arange(1,Lmats)-1,8)
    call bethe_lattice(wt,epsik,Lk,1.d0)
    ts=0.5d0
    sigma(1,:)=impSmats(1,1,1,1,:)
    sigma(2,:)=impSAmats(1,1,1,1,:)!-conjg(impSAmats(1,1,1,1,:))
    fg=zero
    do i=1,L
       iw   = xi*wm(i)
       zita = iw + xmu - sigma(1,i)
       do ik=1,Lk
          cdet = abs(zita-epsik(ik))**2 + sigma(2,i)**2
          fg(1,i)=fg(1,i) + wt(ik)*(conjg(zita)-epsik(ik))/cdet
          fg(2,i)=fg(2,i) - wt(ik)*sigma(2,i)/cdet
       enddo
    enddo
    !fg(2,:)=-conjg(fg(2,:))
    u    = uloc(1)
    n    = ed_dens(1)/2.d0
    delta= ed_phisc(1)*u

    !Get asymptotic self-energies
    Sigma_infty =   dreal(sigma(1,L))
    S_infty     =   dreal(sigma(2,L))

    checkP=0.d0 ; checkdens=0.d0 ;          ! test variables

    kin=0.d0                      ! kinetic energy (generic)
    Ds=0.d0                       ! superfluid stiffness (Bethe)
    do ik=1,Lk
       csi            = epsik(ik)-(xmu-Sigma_infty)
       Ei             = dsqrt(csi**2 + S_infty**2)
       thermal_factor = dtanh(0.5d0*beta*Ei)
       free(ik)        = 0.5d0*(1.d0 - csi/Ei)*thermal_factor
       Ffree(ik)       =-(0.5d0*S_infty)/Ei*thermal_factor
       fmatssum= 0.d0
       matssum = 0.d0
       Dssum   = 0.d0
       vertex=(4.d0*ts**2-epsik(ik)**2)/3.d0
       do i=1,L
          iw       = xi*wm(i)
          det      = abs(iw+xmu-epsik(ik)-sigma(1,i))**2 + dreal(sigma(2,i))**2
          det_infty= wm(i)**2 + (epsik(ik)-(xmu-Sigma_infty))**2 + S_infty**2
          gkw = (-iw+xmu - epsik(ik) - conjg(sigma(1,i)) )/det
          fkw = -sigma(2,i)/det
          g0kw= (-iw - (epsik(ik)-(xmu-Sigma_infty)))/det_infty
          f0kw=-S_infty/det_infty
          matssum =  matssum +  dreal(gkw)-dreal(g0kw)
          fmatssum= fmatssum +  dreal(fkw)-dreal(f0kw)
          Dssum   = Dssum    +  fkw*fkw
       enddo
       n_k(ik)   = 4.d0/beta*matssum + 2.d0*free(ik)
       checkP    = checkP    - wt(ik)*(2.d0/Beta*fmatssum+Ffree(ik))
       checkdens = checkdens + wt(ik)*n_k(ik)
       kin    = kin    + wt(ik)*n_k(ik)*epsik(ik)
       Ds=Ds + 8.d0/beta* wt(ik)*vertex*Dssum
    enddo

    kinsim=0.d0
    kinsim = sum(fg(1,:)*fg(1,:)+conjg(fg(1,:)*fg(1,:))-2.d0*fg(2,:)*fg(2,:))*2.d0*ts**2/beta

    Epot=zero
    Epot = sum(fg(1,:)*sigma(1,:) + fg(2,:)*sigma(2,:))/beta*2.d0

    docc = 0.5d0*n**2
    if(u > 0.01d0)docc=-Epot/u + n - 0.25d0

    Eint=kin+Epot

    Ds=zero
    Ds = sum(fg(2,:)*fg(2,:))/beta*2.d0

    write(*,*)"Asymptotic Self-Energies",Sigma_infty, S_infty
    write(*,*)"n,delta",n,delta
    write(*,*)"Dn% ,Ddelta%",(n-0.5d0*checkdens)/n,(delta + u*checkP)/delta ! u is positive
    write(*,*)'========================================='
    write(*,*)"Kinetic energy",kin
    write(*,*)'========================================='
    write(*,*)"double occupancy   =",docc
    write(*,*)'========================================='
    write(*,*) 'Kinetic Energy TEST (simple formula)'
    write(*,*) '###ACTHUNG: FOR BETHE ONLY####',kinsim
    write(*,*) 'Dkin%',(kin-kinsim)/kin
    write(*,*)'========================================='
    write(*,*) 'Superfluid stiffness',Ds
    write(*,*) 'Potential Energy U(n_up-1/2)(n_do-1/2)',Epot
    write(*,*) 'Internal Energy',Eint
    write(*,*)'========================================='
    call splot("nk_distribution.ipt",epsik,n_k/2.d0,free)
    open(100,file="columns.ipt")
    write(100,"(11A21)")"1vbias","2u","3beta","4n","5kin","6docc","7Ds","8Epot","9Eint"
    close(100)
    open(200,file="thermodynamics.ipt")
    write(200,"(11F21.12)")0.d0,u,beta,n,kinsim,docc,Ds,Epot,Eint
    close(200)
    return 
  end subroutine get_sc_internal_energy




end program lancED



