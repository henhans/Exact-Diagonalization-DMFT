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
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine imp_geth(isloop)
    integer                  :: ib(N)
    integer                  :: idg
    integer                  :: i,j,k,r,m,ms
    integer                  :: kp,isloop
    real(8),dimension(Nbath) :: eup,edw,vup,vdw
    real(8)                  :: ndup,nddw,npup,npdw,nup,ndw,sg1,sg2,tef
    real(8),pointer          :: h(:,:)

    idg=deg(isloop)
    allocate(h(idg,idg))
    h => espace(isloop)%M
    h=0.d0

    eup=ebath(1,:)
    vup=vbath(1,:)
    edw=eup
    vdw=vup
    if(Nspin==2)then
       edw=ebath(2,:)
       vdw=vbath(2,:)
    endif

    do i=1,idg
       m=nmap(isloop,i)
       call bdecomp(m,ib)

       nup=real(ib(1),8)
       ndw=real(ib(1+Ns),8)

       !Diagonal part
       h(i,i)=(-xmu+ed0)*(nup+ndw) + u*(nup-0.5d0)*(ndw-0.5d0) + heff*(nup-ndw)

       !energy of the bath=\sum_{n=1,N}\e_l n_l
       do kp=2,Ns
          h(i,i)=h(i,i)+eup(kp-1)*real(ib(kp),8)
          h(i,i)=h(i,i)+edw(kp-1)*real(ib(kp+Ns),8)
       enddo

       !non-diagonal part
       do ms=2,Ns
          if(ib(1)==1 .AND. ib(ms)==0)then
             call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
             call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=invnmap(isloop,r)
             tef=vup(ms-1)
             h(i,j)=tef*sg1*sg2
             h(j,i)=h(i,j)!Hermitian conjugate 
          endif
          if(ib(1+Ns)==1 .AND. ib(ms+Ns)==0)then
             call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
             call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
             j=invnmap(isloop,r)
             tef=vdw(ms-1)   
             h(i,j)=tef*sg1*sg2
             h(j,i)=h(i,j)
          endif
       enddo
    enddo
    return
  end subroutine imp_geth
  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_GETH
