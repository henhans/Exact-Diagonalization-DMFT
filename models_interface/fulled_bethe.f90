!###################################################################
!PROGRAM  : FULLED
!TYPE     : main code
!PURPOSE  : Complete ED solution of DMFT equations (STD interface HM/PAM)
!AUTHORS  : A. Amaricci
!COMMENTS : USE THIS INTERFACE TO BUILD YOUR OWN CODE 
!###################################################################
program fullED
  USE VARS_GLOBAL
  USE DIAG
  USE TOFITGF
  USE TOOLS
  implicit none
  integer :: success

  call initialize("inputED.in")

  store_treshold=100
  do iloop=1,nloop
     call ed_solver(iloop,.true.) !Solve the EFFECTIVE IMPURITY PROBLEM
     call get_delta_bethe_integral()
     success = scc_fit(iloop)     !Perform SELF-CONSISTENCY and Check convergency
  enddo
  call finalize(success)

contains

  !+----------------------------------------+
  subroutine get_delta_bethe_integral
    integer :: i,j
    real(8) :: w
    complex(8) :: iw,deltaAnd,g0and,self,zetan,g0loc,deltaLoc
    complex(8) :: gloc(NL),grloc(-Nw:Nw)
    if(allocated(delta))deallocate(delta)
    allocate(delta(NL))

    gloc=zero;grloc=zero
    do i=1,NL
       w=wm(i);iw=xi*w
       deltaAnd=sum(vup(1:Nbath)**2/(iw-epsiup(1:Nbath)))
       g0and = iw + xmu - ed0 -deltaAnd
       self  = g0and - one/Giw(i)
       zetan = iw + xmu - ed0 - self
       gloc(i)=gfbethe(w,zetan,D)
       g0loc=self + one/gloc(i)
       delta(i)= iw+xmu-g0loc
    enddo

    do i=-Nw,Nw
       w=wr(i);iw=cmplx(w,eps)
       deltaAnd=sum(vup(1:Nbath)**2/(iw-epsiup(1:Nbath)))
       g0and = iw + xmu - ed0 -deltaAnd
       self = g0and - one/Gwr(i)       
       zetan=iw + xmu - ed0 - self
       grloc(i)=gfbethe(w,zetan,D)
    enddo
    call splot("locG_iw.data",wm(1:NL0),gloc(1:NL0))
    call splot("locG_realw.data",wr,grloc)
    return    
  end subroutine get_delta_bethe_integral
  !+----------------------------------------+

end program fullED



