!+------------------------------------------------------+
!PROGRAM  : 
!TYPE     : subroutine
!PURPOSE  :On OUTPUT \Delta = g
!STRATEGY : Build the imp Sigma function on the fly
!           Use Sigma to get Glocal=sum_k 1/iw+mu-Sigma
!           Use Gloc and Sigma to get G0loc
!           Use G0loc to get Delta_loc
!+------------------------------------------------------+
subroutine getdelta_bethe_integral
  integer :: i,j
  real(8) :: w,en0
  complex(8) :: iw,deltaAnd,g0and,self,zetan,g0loc,deltaLoc
  complex(8) :: gloc(NL),grloc(-Nw:Nw)
  if(allocated(delta))deallocate(delta)
  allocate(delta(NL))
  en0=ed0
  gloc=zero;delta=zero;grloc=zero
  do i=1,NL
     w=wm(i);iw=xi*w
     deltaAnd=zero
     do j=1,Nbath
        deltaAnd=deltaAnd+vup(j)**2/(iw-epsiup(j))
     enddo
     g0and = iw + xmu - en0 -deltaAnd
     self  = g0and - one/Giw(i)
     zetan = iw + xmu - en0 - self
     do j=1,Lk
        gloc(i)=gloc(i)+wt(j)/(zetan-epsik(j))
     enddo
     g0loc=self + one/gloc(i)
     delta(i)= iw+xmu-g0loc
  enddo
  do i=-Nw,Nw
     w=wr(i);iw=cmplx(w,eps)
     deltaAnd=zero
     do j=1,Nbath
        deltaAnd=deltaAnd+vup(j)**2/(iw-epsiup(j))
     enddo
     g0and = iw + xmu - en0 -deltaAnd
     self = g0and - one/Gwr(i)       
     zetan=iw + xmu - en0 - self
     do j=1,Lk
        grloc(i)=grloc(i)+wt(j)/(zetan-epsik(j))
     enddo
  enddo
  call splot("locG_iw.data",wm(1:NL0),gloc(1:NL0))
  call splot("locG_realw.data",wr,grloc)
  return    
end subroutine getdelta_bethe_integral
!*********************************************************************
!*********************************************************************
!*********************************************************************
