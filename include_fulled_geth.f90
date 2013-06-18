!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_geth(isector,h)
  real(8)                       :: h(:,:)
  integer                       :: ib(Ntot)
  integer                       :: dim
  integer                       :: i,j,k,r,m,ms,iorb,ispin,iup,idw
  integer                       :: kp,isector
  real(8),dimension(Nbath)      :: eup,edw
  real(8),dimension(Norb,Nbath) :: vup,vdw
  real(8)                       :: nup(Norb),ndw(Norb),sg1,sg2,tef,htmp
  dim=getdim(isector)
  if(size(h,1)/=dim)call error("FULL_ED_GETH: wrong dimension 1 of H")
  if(size(h,2)/=dim)call error("FULL_ED_GETH: wrong dimension 2 of H")
  h=0.d0
  eup=ebath(1,:)   ; edw=ebath(Nspin,:)
  vup=vbath(:,1,:) ; vdw=vbath(:,Nspin,:)
  do i=1,dim
     m=Hmap(isector)%map(i)
     call bdecomp(m,ib)
     do iorb=1,Norb
        nup(iorb)=real(ib(iorb),8)
        ndw(iorb)=real(ib(iorb+Ns),8)
     enddo
     !LOCAL HAMILTONIAN PART:
     !diagonal terms:
     htmp = 0.d0
     htmp = -xmu*(sum(nup)+sum(ndw)) + heff*(sum(nup)-sum(ndw))
     if(Norb > 1)htmp =  htmp + ep0*(nup(2)+ndw(2)) !we should find a better way here:
     select case(hfmode)
     case(.true.)
        htmp =htmp + u*(nup(1)-0.5d0)*(ndw(1)-0.5d0)
     case (.false.)
        htmp =htmp + u*nup(1)*ndw(1)
     end select
     !Hbath:
     do kp=Norb+1,Ns
        htmp =htmp + eup(kp-Norb)*real(ib(kp),8)
        htmp =htmp + edw(kp-Norb)*real(ib(kp+Ns),8)
     enddo
     h(i,i)=htmp
     !
     !Hybridization terms
     if(tpd/=0.d0.AND.Norb>1)then
        if(ib(1) == 1 .AND. ib(2) == 0)then
           call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
           call cdg(2,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
           j=invHmap(isector,r)
           tef=tpd
           h(i,j)=tef*sg1*sg2
           h(j,i)=h(i,j)
        endif
        if(ib(1+Ns) == 1 .AND. ib(2+Ns) == 0)then
           call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
           call cdg(2+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
           j=invHmap(isector,r)
           tef=tpd
           h(i,j)=tef*sg1*sg2
           h(j,i)=h(i,j)
        endif
     endif
     !
     !NON-LOCAL PART:
     do iorb=1,Norb
        do ms=Norb+1,Ns
           !UP
           if(ib(iorb) == 1 .AND. ib(ms) == 0)then
              call c(iorb,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
              call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
              j=invHmap(isector,r)
              tef=vup(iorb,ms-Norb)
              h(i,j)=tef*sg1*sg2
              h(j,i)=h(i,j)
           endif
           !DW
           if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
              call c(iorb+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
              call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
              j=invHmap(isector,r)
              tef=vdw(iorb,ms-Norb)
              h(i,j)=tef*sg1*sg2
              h(j,i)=h(i,j)
           endif
        enddo
     enddo
  enddo
  return
end subroutine full_ed_geth

