!###################################################################
!PURPOSE  : LANCZOS ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program lancED
  USE DMFT_ED
  USE COMMON_VARS
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
  call parse_input_variable(wband,"wband",finput,default=1.d0)
  call parse_input_variable(Lk,"Lk",finput,default=500)
  !
  call ed_read_input(trim(finput))!,trim(fhloc))

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
     if(mpiID==0)then
        converged = check_convergence(delta(1,1,1,:)+delta(2,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
        if(nread/=0.d0)call search_chemical_potential(ed_dens(1),converged)
     endif
     call end_loop
  enddo

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

end program lancED



