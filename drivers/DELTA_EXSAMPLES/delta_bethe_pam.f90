!+------------------------------------------------------+
!PROGRAM  : 
!TYPE     : subroutine
!PURPOSE  :On OUTPUT \Delta = g
!+------------------------------------------------------+
subroutine getdelta_bethe_pam
  integer    :: i,j
  real(8)    :: w
  complex(8) :: iw,deltaAnd,g0d,sd,sp,zita,gp
  complex(8) :: gploc(NL),gdloc(NL)
  complex(8) :: gprloc(-Nw:Nw),gdrloc(-Nw:Nw)

  if(allocated(delta))deallocate(delta)
  allocate(delta(NL))

  do i=1,NL
     w=wm(i);iw=xi*w
     deltaAnd=zero
     do j=1,Nbath
        deltaAnd=deltaAnd+vup(j)**2/(iw-epsiup(j))
     enddo
     g0d = iw + xmu - ed0 - deltaAnd
     sd  = g0d - one/Giw(i)
     sp  = tpd**2/(iw + xmu - ed0 - sd)
     zita= iw + xmu -ep0 - sp
     gp=zero
     do j=1,Lk
        gp=gp+wt(j)/(zita-epsik(j))
     enddo
     gploc(i)=gp
     ! gploc(i)  = gfbethe(w,zita,d)
     gdloc(i)  = one/(iw+xmu-ed0-sd) + tpd**2/(iw+xmu-ed0-sd)**2*gploc(i)
     delta(i)=tpd**2/(iw+xmu-d**2*gploc(i)/4.d0)
  enddo

  do i=-Nw,Nw
     w=wr(i);iw=cmplx(w,eps)
     deltaAnd=zero
     do j=1,Nbath
        deltaAnd=deltaAnd+vup(j)**2/(iw-epsiup(j))
     enddo
     g0d = iw + xmu - ed0 - deltaAnd     
     sd  = g0d - one/Gwr(i)
     sp  = tpd**2/(iw + xmu - ed0 - sd)
     zita= iw + xmu -ep0 - sp
     gp=zero
     do j=1,Lk
        gp=gp+wt(j)/(zita-epsik(j))
     enddo
     gprloc(i)=gp
     ! gprloc(i)  = gfbether(w,zita,D)
     gdrloc(i)  = one/(iw+xmu-ed0-sd) + tpd**2/(iw+xmu-ed0-sd)**2*gprloc(i)
  enddo

  call splot("locGdd_iw.data",wm(1:NL0),gdloc(1:NL0))
  call splot("locGpp_iw.data",wm(1:NL0),gploc(1:NL0))
  call splot("locGdd_realw.data",wr,gdrloc)
  call splot("locGpp_realw.data",wr,gprloc)

  return    
end subroutine getdelta_bethe_pam
!*********************************************************************
!*********************************************************************
!*********************************************************************
