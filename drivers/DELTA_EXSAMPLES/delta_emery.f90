!+------------------------------------------------------+
!PROGRAM  : 
!TYPE     : subroutine
!PURPOSE  :On OUTPUT \Delta = g
!+------------------------------------------------------+
subroutine getdelta_EMERY()
  implicit none
  integer :: i,j,ix,iy,ik
  real(8) :: w,sx,sy,ndd,npx,npy,xrr,normx,normy,delt
  real(8) :: xk,yk,num1,num2,tp,deltaw,dosx,dosy
  real(8),dimension(0:NL) :: Gddt,Gpxt,Gpyt
  complex(8) :: iw,ggp,tpdx,tpdy,frac1,frac2,g0
  complex(8) :: detG,alpha,gamma
  complex(8),dimension(-Nw:Nw) :: Gdd_r,Sigma_r,Gpx_r,Gpy_r
  complex(8),dimension(NL) :: Gdd,Sigma,Gpx,Gpy

  if(allocated(delta))deallocate(delta)
  allocate(delta(NL))

  delt=0.01

  Gdd=zero
  Gpx=zero
  Gpy=zero
  do i=1,NL
     w=wm(i);iw=xi*w
     g0=zero
     do j=1,Nbath
        g0=g0+vup(j)**2/(iw-epsiup(j))
     enddo
     g0 = iw + xmu - g0  - ed0
     sigma(i) = g0 - one/Giw(i)       
  enddo

  Gdd_r=zero
  Gpx_r=zero
  Gpy_r=zero
  do i=1,Nw
     w=wr(i)
     g0=zero
     do j=1,Nbath
        g0=g0+vup(j)**2/(w-xi*eps-epsiup(j))
     enddo
     g0 = w + xmu  - g0 - ed0
     sigma_r(i) = g0 - one/Gwr(i)
  enddo

  !get Glocal
  do ik=1,Lk
     ix=ik2ix(ik)
     iy=ik2iy(ik)
     xk=kgrid(ix,iy)%x
     yk=kgrid(ix,iy)%y

     sx=sin(xk/2) ; sy=sin(yk/2)
     tpdx=2*xi*tpd*sx
     tpdy=2*xi*tpd*sy
     tp = 4*tpp*sx*sy 
     num1=(tpdx**2 + tpdy**2)
     num2=2.d0*tp*tpdx*tpdy

     do i=1,NL
        w=wm(i)
        iw=xi*w + xmu
        alpha=iw-sigma(i) -ed0
        gamma=iw-ep0
        frac1=num1/(gamma - tp**2/gamma)
        frac2=num2/(gamma**2 -tp**2)
        Gdd(i)=Gdd(i)+wt(ik)/(alpha+frac1+frac2)
        detG=num2+alpha*(gamma**2-tp**2)+gamma*(tpdx**2+tpdy**2)
        Gpx(i)=Gpx(i)+wt(ik)*(tpdy**2+alpha*gamma)/detG
        Gpy(i)=Gpy(i)+wt(ik)*(tpdx**2+alpha*gamma)/detG   
     enddo
     do i=-Nw,Nw
        w=wr(i)
        iw=cmplx(w+xmu,eps)
        alpha=iw-sigma_r(i)-ed0
        gamma=iw-ep0
        frac1=num1/(gamma - tp**2/gamma)
        frac2=num2/(gamma**2 -tp**2)
        Gdd_r(i)=Gdd_r(i)+wt(ik)/(alpha+frac1+frac2)
        detG=num2+alpha*(gamma**2-tp**2)+gamma*(tpdx**2+tpdy**2)
        Gpx_r(i)=Gpx_r(i)+wt(ik)*(tpdy**2+alpha*gamma)/detG
        Gpy_r(i)=Gpy_r(i)+wt(ik)*(tpdx**2+alpha*gamma)/detG
     enddo
  enddo

  open(90,file="Gddloc_iw.data")
  open(91,file="Gpxloc_iw.data")
  open(92,file="Gpyloc_iw.data")
  open(93,file="avGpploc_iw.data")
  open(61,file='Deltaloc_iw.data')
  do i=1,NL
     w=pt*dble(2*i-1)
     write(90,*)w,aimag(Gdd(i)),real(Gdd(i))
     write(91,*)w,aimag(Gpx(i)),real(Gpx(i))
     write(92,*)w,aimag(Gpy(i)),real(Gpy(i))
     ggp=(Gpx(i) + Gpy(i))/2.d0
     write(93,*)w,aimag(ggp),real(ggp)
     !Get Delta
     g0 = sigma(i) + one/Gdd(i)
     delta(i)=xi*w + xmu - g0 - ed0
     write(61,*)w,aimag(delta(i)),real(delta(i))
  enddo
  close(90);close(91);close(92);close(93);close(61)

  open(90,file="Gddloc_realw.data")
  open(91,file="Gpxloc_realw.data")
  open(92,file="Gpyloc_realw.data")
  open(93,file="avGpploc_realw.data")
  do i=-Nw,Nw
     w=wr(i)
     write(90,*)w,-aimag(Gdd_r(i))/pi
     write(91,*)w,-aimag(Gpx_r(i))/pi
     write(92,*)w,-aimag(Gpy_r(i))/pi
     write(93,*)w,-(aimag(Gpx_r(i))+aimag(Gpy_r(i)))/2.d0/pi
  enddo
  close(90);close(91);close(92);close(93)


  open(20,file="nddVSiloop.data",access="append")
  open(21,file="npxVSiloop.data",access="append")
  open(22,file="npyVSiloop.data",access="append")
  open(23,file="ntotVSiloop.data",access="append")
  npx=0.d0
  npy=0.d0
  normx=0.d0
  normy=0.d0
  do i=-Nw,Nw-1
     w=wr(i)
     dosx=-dimag(Gpx_r(i))/pi
     dosy=-dimag(Gpy_r(i))/pi
     npx=npx+deltaw*fermi(w,beta)*dosx
     npy=npy+deltaw*fermi(w,beta)*dosy
     normx=normx+deltaw*dosx
     normy=normy+deltaw*dosy
  enddo
  npx=2.d0*npx/normx
  npy=2.d0*npy/normy

  call cfft_iw2tau(Gdd,Gddt,beta)
  ndd=2.d0*(1.d0-Gddt(0))

  write(20,*)iloop,ndd
  write(21,*)iloop,npx
  write(22,*)iloop,npy
  write(23,*)iloop,npimp+npx+npy
  close(20);close(21);close(22);close(23)


  if(last)then
     call cfft_iw2tau(Gpx,Gpxt,beta)
     call cfft_iw2tau(Gpy,Gpyt,beta)
     open(90,file="Gddloc_tau.data")
     open(91,file="Gpxloc_tau.data")
     open(92,file="Gpyloc_tau.data")
     do i=0,NL
        write(90,*)dble(i)*beta/dble(NL),Gddt(i)
        write(91,*)dble(i)*beta/dble(NL),Gpxt(i)
        write(92,*)dble(i)*beta/dble(NL),Gpyt(i)
     enddo
     close(90);close(91);close(92)

     if(itu.eq.'m')xrr=xmu
     if(itu.eq.'t')xrr=temp
     if(itu.eq.'u')xrr=u
     open(20,file='nddVS'//extension,access='append')
     open(21,file='npxVS'//extension,access='append')
     open(22,file='npyVS'//extension,access='append')
     open(23,file='ntotVS'//extension,access='append')
     write(20,*)xrr,ndd
     write(21,*)xrr,npx,2.d0*(1.d0-Gpxt(0))
     write(22,*)xrr,npy,2.d0*(1.d0-Gpyt(0))
     write(23,*)xrr,npimp+npx+npy,&
          2.d0*(1.d0-Gpxt(0))+2.d0*(1.d0-Gpyt(0))+npimp
     close(20);close(21);close(22);close(23)
  endif
end subroutine getdelta_EMERY
!*********************************************************************
!*********************************************************************
!*********************************************************************
