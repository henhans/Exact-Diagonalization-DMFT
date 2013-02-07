!###################################################################
!PURPOSE  : Fulle ED solution of DMFT problem for Hubbard model.
!AUTHORS  : A. Amaricci
!###################################################################
program fullED
  USE DMFT_FULLED
  implicit none
  logical :: converged

  call read_input("inputED.in")

  !Setup solver:
  call ed_solver(status=-2)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solver() 

     !Get the Weiss field/Delta function to be fitted (user defined)
     call get_delta_bethe

     !Perform the SELF-CONSISTENCY by fitting the new bath
     call chi2_fitgf(delta(1,:),epsiup,vup)
     epsidw=epsiup;vdw=vup

     !Check convergence (if required change chemical potential)
     converged = check_convergence(delta(1,:),eps_error,nsuccess,nloop)
     if(nread/=0.d0)call search_mu(nimp1,converged)
     call end_loop
  enddo

  !Finalize calculation
  call ed_solver(status=-1)


contains


  !+----------------------------------------+
  subroutine get_delta_bethe
    integer    :: i,j
    real(8)    :: gtau(0:Ltau)
    complex(8) :: iw,deltaAnd,g0and,zetan,g0loc,deltaLoc
    complex(8) :: gloc(NL),grloc(Nw),self(NL),selfr(Nw)

    gloc=zero;grloc=zero
    do i=1,NL
       iw=xi*wm(i)
       g0and= iw + xmu -ed0 -delta_and(iw,epsiup,vup)
       self(i)  = g0and - one/Giw(1,i)
       zetan = iw + xmu - ed0 - self(i)
       gloc(i)=gfbethe(wm(i),zetan,D)
       g0loc=self(i) + one/gloc(i)
       delta(1,i)= iw+xmu-g0loc
    enddo

    do i=1,Nw
       iw=cmplx(wr(i),eps)
       g0and = wr(i) + xmu - ed0 - delta_and(wr(i)+zero,epsiup,vup)
       selfr(i) = g0and - one/Gwr(1,i)    
       zetan=iw + xmu - ed0 - selfr(i)
       grloc(i)=gfbether(wr(i),zetan,D)
    enddo
    call splot("G_iw.ed",wm,gloc)
    call splot("Sigma_iw.ed",wm,self)
    call splot("G_realw.ed",wr,grloc)
    call splot("Sigma_realw.ed",wr,selfr)
    call splot("Delta_iw.ed",wm,delta(1,:))
    return    
  end subroutine get_delta_bethe
  !+----------------------------------------+

end program fullED



