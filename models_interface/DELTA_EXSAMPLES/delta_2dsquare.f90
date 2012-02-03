!+-------------------------------------------------------------------+
!PROGRAM  : 
!TYPE     : Subroutine
!PURPOSE  : 
!COMMENT  : 
!+-------------------------------------------------------------------+
subroutine getdelta_2dsquare()
  integer          :: i,j,ik
  real(8)          :: w,deltaw
  complex(8)       :: iw,zeta,deltaAnd,g0and,self,gloc,g0loc,deltaLoc
  complex(8),dimension(NL)   :: Glociw,Siw
  real(8),dimension(0:NL)    :: dummyG
  real(8),dimension(0:Ltau)    :: Gtau
  complex(8),dimension(-Nw:Nw) :: Glocwr,Swr
  if(allocated(delta))deallocate(delta)
  allocate(delta(NL))
  do i=1,NL
     w=wm(i); iw=xi*w
     deltaAnd=zero
     do j=1,Nbath
        deltaAnd=deltaAnd+vup(j)**2/(iw-epsiup(j))
     enddo
     g0and = iw + xmu - ed0 -deltaAnd     
     self = g0and - one/Giw(i)       
     zeta=iw+xmu-self
     gloc=zero
     do ik=1,Lk
        gloc=gloc + wt(ik)/(zeta - epsik(ik))
     enddo
     g0loc=self + one/gloc
     deltaLoc= iw+xmu-g0loc
     delta(i)=deltaLoc

     Glociw(i)=gloc
     Siw(i)=self
  enddo
  call cfft_iw2tau(Glociw,dummyG,beta)
  call extract(dummyG,Gtau)

  do i=-Nw,Nw
     w=wr(i);iw=cmplx(w,eps)
     deltaAnd=zero
     do j=1,Nbath
        deltaAnd=deltaAnd+vup(j)**2/(iw-epsiup(j))
     enddo
     g0and= iw + xmu - ed0 -deltaAnd     
     self = g0and - one/Gwr(i)       
     zeta = iw + xmu - self
     gloc = zero
     do ik=1,Lk
        gloc=gloc + wt(ik)/(zeta - epsik(ik))
     enddo
     Glocwr(i)=gloc
     Swr(i)=self
  enddo

  call splot("locG_iw.data",wm(1:NL),Glociw(1:NL))
  call splot("locS_iw.data",wm(1:NL),Siw(1:NL))

  call splot("locG_realw.data",wr,Glocwr)
  call splot("locS_realw.data",wr,Swr)

  call splot("locG_tau.data",tau(0:Ltau),Gtau(0:Ltau))

  return
end subroutine getdelta_2dsquare
!*******************************************************************
!*******************************************************************
!*******************************************************************
