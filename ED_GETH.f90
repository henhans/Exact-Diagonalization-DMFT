!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_GETH
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, only: cdg,c,bdecomp
  implicit none
  private
  public :: imp_geth
  save

contains

  !+------------------------------------------------------------------+
  !PROGRAM  : 
  !TYPE     : subroutine
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_geth(isloop)
    !LOCAL VARIABLES
    !Configuration vector
    integer :: ib(N)
    integer :: idg
    integer :: i,j,k,r,m,ms
    integer :: kp,isloop
    integer :: NR
    real(8) :: ndup,nddw,npup,npdw,nup,ndw,sg1,sg2,tef
    real(8),pointer :: h(:,:)

    idg=deg(isloop)
    allocate(h(idg,idg))
    h => espace(isloop)%M
    h=0.d0

    NR=Nimp + 1
    do i=1,idg
       m=nmap(isloop,i)
       call bdecomp(m,ib)

       if(Nimp==1)then
          nup=real(ib(1),8)
          ndw=real(ib(1+Ns),8)
          !Diagonal part
          h(i,i)=(-xmu+ed0)*(nup+ndw) + u*(nup-0.5d0)*(ndw-0.5d0) + heff*(nup-ndw)
          !energy of the bath=\sum_{n=1,N}\e_l n_l
          do kp=2,Ns
             h(i,i)=h(i,i)+epsiup(kp-1)*real(ib(kp),8)
             h(i,i)=h(i,i)+epsidw(kp-1)*real(ib(kp+Ns),8)
          enddo
          !NON-Diagonal part
          do ms=2,Ns
             if(ib(1) == 1 .AND. ib(ms) == 0)then
                call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
                j=invnmap(isloop,r)
                tef=vup(ms-1)
                h(i,j)=tef*sg1*sg2
                h(j,i)=h(i,j)!Hermitian conjugate 
             endif
             if(ib(1+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
                j=invnmap(isloop,r)
                tef=vdw(ms-1)   
                h(i,j)=tef*sg1*sg2
                h(j,i)=h(i,j)
             endif
          enddo

       elseif(Nimp==2)then
          ndup=dble(ib(1))
          nddw=dble(ib(1+Ns))
          npup=dble(ib(2))
          npdw=dble(ib(2+Ns))
          !Diagonal part
          h(i,i)=(-xmu+ed0)*(ndup+nddw)    &
               +u*(ndup-0.5d0)*(nddw-0.5d0)&
               +heff*(ndup-nddw)           &  
               +(-xmu)*(npup+npdw)         &
               + heff*(npup-npdw)
          !energy of the bath=\sum_{n=1,N}\e_l n_l
          do kp=3,Ns
             h(i,i)=h(i,i)+epsiup(kp-2)*(dfloat(ib(kp)))
             h(i,i)=h(i,i)+epsidw(kp-2)*(dfloat(ib(kp+Ns)))
          enddo
          !NON-Diagonal part       
          !UP-SPIN PART
          do ms=3,Ns
             if(ib(2) == 1 .AND. ib(ms) == 0)then
                call c(2,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
                j=invnmap(isloop,r)
                tef=vup(ms-2)
                h(i,j)=tef*sg1*sg2
                h(j,i)=h(i,j)!Hermitian conjugate 
             endif
          enddo
          !DW-SPIN PART
          do ms=3,Ns
             if(ib(2+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(Nimp+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
                j=invnmap(isloop,r)
                tef=vdw(ms-2)   
                h(i,j)=tef*sg1*sg2
                h(j,i)=h(i,j)
             endif
          enddo
          !Hybridization part
          if(ib(1) == 1 .AND. ib(2) == 0)then
             call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
             call cdg(2,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=invnmap(isloop,r)
             tef=tpd
             h(i,j)=tef*sg1*sg2
             h(j,i)=h(i,j)
          endif
          if(ib(1+Ns) == 1 .AND. ib(2+Ns) == 0)then
             call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
             call cdg(2+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=invnmap(isloop,r)
             tef=tpd
             h(i,j)=tef*sg1*sg2
             h(j,i)=h(i,j)
          endif
       endif
    enddo
    return
  end subroutine imp_geth
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_GETH
