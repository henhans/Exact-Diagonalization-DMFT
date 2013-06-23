!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_geth(isector,h)
  real(8)                  :: h(:,:)
  integer                  :: ib(Ntot)
  integer                  :: dim
  integer                  :: i,j,k,r,m,ms
  integer                  :: kp,isector
  real(8),dimension(Nbath) :: eup,edw,vup,vdw
  real(8)                  :: ndup,nddw,npup,npdw,nup,ndw,sg1,sg2,tef

  dim=getdim(isector)
  if(size(h,1)/=dim)call error("FULL_ED_GETH: wrong dimension 1 of H")
  if(size(h,2)/=dim)call error("FULL_ED_GETH: wrong dimension 2 of H")

  h=0.d0

  eup=ebath(1,:)
  vup=vbath(1,:)
  edw=eup
  vdw=vup
  if(Nspin==2)then
     edw=ebath(2,:)
     vdw=vbath(2,:)
  endif

  do i=1,dim
     m=Hmap(isector)%map(i)!Hmap(isector,i)
     call bdecomp(m,ib)

     nup=real(ib(1),8)
     ndw=real(ib(1+Ns),8)

     !Diagonal part
     select case(hfmode)
     case(.true.)
        h(i,i)= -xmu*(nup+ndw) + U*(nup-0.5d0)*(ndw-0.5d0) + heff*(nup-ndw)
     case (.false.)
        h(i,i)= -(xmu+U/2d0)*(nup+ndw) + U*nup*ndw + heff*(nup-ndw)
     end select

     !energy of the bath=\sum_{n=1,Ntot}\e_l n_l
     do kp=2,Ns
        h(i,i)=h(i,i)+eup(kp-1)*real(ib(kp),8)
        h(i,i)=h(i,i)+edw(kp-1)*real(ib(kp+Ns),8)
     enddo

     !non-diagonal part
     do ms=2,Ns
        if(ib(1)==1 .AND. ib(ms)==0)then
           call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
           call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
           j=invHmap(isector,r)
           tef=vup(ms-1)
           h(i,j)=tef*sg1*sg2
           h(j,i)=h(i,j)!Hermitian conjugate 
        endif
        if(ib(1+Ns)==1 .AND. ib(ms+Ns)==0)then
           call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
           call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
           j=invHmap(isector,r)
           tef=vdw(ms-1)   
           h(i,j)=tef*sg1*sg2
           h(j,i)=h(i,j)
        endif
     enddo
  enddo
  return
end subroutine full_ed_geth
