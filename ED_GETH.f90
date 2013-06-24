!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_GETH
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none
  private
  public :: full_ed_geth
  public :: lanc_ed_geth
  public :: spHtimesV_d,spHtimesV_c
  public :: set_Hsector

  integer :: Hsector

contains

  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  include 'include_fulled_geth.f90'



  !####################################################################
  !                    LANCZOS DIAGONALIZATION (T=0, GS only)
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine  set_Hsector(isector)
    integer :: isector
    Hsector=isector
  end subroutine set_Hsector


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine spHtimesV_d(N,v,Hv)
    integer              :: N
    real(8),dimension(N) :: v
    real(8),dimension(N) :: Hv
    Hv=zero
    call sp_matrix_vector_product(N,spH0,v,Hv)
  end subroutine SpHtimesV_d

  subroutine spHtimesV_c(N,v,Hv)
    integer              :: N
    complex(8),dimension(N) :: v
    complex(8),dimension(N) :: Hv
    Hv=zero
    call sp_matrix_vector_product(N,spH0,v,Hv)
  end subroutine SpHtimesV_c



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine lanc_ed_geth(isector)
    integer                       :: isector
    integer                       :: ib(Ntot)
    integer                       :: dim,iup,idw
    integer                       :: i,j,k,r,m,ms,iorb,ispin
    integer                       :: kp
    real(8),dimension(Nbath)      :: eup,edw
    real(8),dimension(Norb,Nbath) :: vup,vdw
    real(8)                       :: nup(Norb),ndw(Norb),sg1,sg2,tef,htmp

    dim=getdim(isector)
    if(.not.spH0%status)&
         call error("LANC_ED_GETH: spH0 not initialized at sector:"//txtfy(isector))

    eup=ebath(1,:)   ; edw=ebath(Nspin,:)
    vup=vbath(:,1,:) ; vdw=vbath(:,Nspin,:)

    do i=1,dim
       m=Hmap(isector)%map(i)!Hmap(isector,i)
       call bdecomp(m,ib)
       do iorb=1,Norb
          nup(iorb)=real(ib(iorb),8)
          ndw(iorb)=real(ib(iorb+Ns),8)
       enddo
       !LOCAL HAMILTONIAN PART:
       !diagonal terms:
       htmp=0.d0
       htmp = -xmu*(sum(nup)+sum(ndw)) + heff*(sum(nup)-sum(ndw))
       if(Norb > 1)htmp = htmp + ep0*(nup(2)+ndw(2)) !find a better way to include Crystal field
       select case(hfmode)
       case(.true.)
          htmp =htmp + u*(nup(1)-0.5d0)*(ndw(1)-0.5d0)
       case (.false.)
          htmp =htmp + u*nup(1)*ndw(1)
       end select
       !+energy of the bath=\sum_{n=1,N}\e_l n_l
       do kp=Norb+1,Ns
          htmp =htmp + eup(kp-Norb)*real(ib(kp),8) + edw(kp-Norb)*real(ib(kp+Ns),8)
       enddo
       call sp_insert_element(spH0,htmp,i,i)
       !
       !Hybridization terms
       if(tpd/=0.d0.AND.Norb>1)then
          if(ib(1) == 1 .AND. ib(2) == 0)then
             call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
             call cdg(2,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=invHmap(isector,r)
             tef=tpd
             htmp=tef*sg1*sg2
             call sp_insert_element(spH0,htmp,i,j)
             call sp_insert_element(spH0,htmp,j,i)
          endif
          if(ib(1+Ns) == 1 .AND. ib(2+Ns) == 0)then
             call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
             call cdg(2+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
             j=invHmap(isector,r)
             tef=tpd
             htmp=tef*sg1*sg2
             call sp_insert_element(spH0,htmp,i,j)
             call sp_insert_element(spH0,htmp,j,i)
          endif
       endif
       !
       !NON-LOCAL PART
       do iorb=1,Norb
          do ms=Norb+1,Ns
             !UP
             if(ib(iorb) == 1 .AND. ib(ms) == 0)then
                call c(iorb,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
                j=invHmap(isector,r)
                tef=vup(iorb,ms-Norb)
                htmp = tef*sg1*sg2
                call sp_insert_element(spH0,htmp,i,j)
                call sp_insert_element(spH0,htmp,j,i)
             endif

             !DW
             if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(iorb+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
                j=invHmap(isector,r)
                tef=vdw(iorb,ms-Norb)
                htmp=tef*sg1*sg2
                call sp_insert_element(spH0,htmp,i,j)
                call sp_insert_element(spH0,htmp,j,i)
             endif
          enddo
       enddo


    enddo
    return
  end subroutine lanc_ed_geth


  !*********************************************************************
  !*********************************************************************
  !*********************************************************************



  ! subroutine HtimesV(Nv,v,Hv)
  !   integer                  :: Nv
  !   real(8),dimension(Nv)    :: v
  !   real(8),dimension(Nv)    :: Hv
  !   integer                  :: isector
  !   integer                  :: ib(Ntot)
  !   integer                  :: dim
  !   integer                  :: i,j,k,r,m,ms,ispin
  !   integer                  :: kp
  !   integer                  :: iimp,ibath
  !   real(8),dimension(Nbath) :: eup,edw,vup,vdw
  !   real(8)                  :: ndup,nddw,npup,npdw,nup,ndw,sg1,sg2,tef,htmp

  !   isector=Hsector
  !   dim=getdim(isector)
  !   eup=ebath(1,:)
  !   vup=vbath(1,:)
  !   edw=eup
  !   vdw=vup
  !   if(Nspin==2)then
  !      edw=ebath(2,:)
  !      vdw=vbath(2,:)
  !   endif

  !   if(Nv/=dim)call error("HtimesV error in dimensions")
  !   Hv=0.d0

  !   do i=1,dim
  !      m=Hmap(isector)%map(i)
  !      call bdecomp(m,ib)
  !      nup=real(ib(1),8)
  !      ndw=real(ib(1+Ns),8)

  !      !Diagonal part
  !      !local part of the impurity Hamiltonian: (-mu+\e0)*n + U*(n_up-0.5)*(n_dw-0.5) + heff*mag
  !      !+ energy of the bath=\sum_{n=1,N}\e_l n_l
  !      htmp=0.d0
  !      select case(hfmode)
  !      case(.true.)
  !         htmp = -xmu*(nup+ndw) + U*(nup-0.5d0)*(ndw-0.5d0) + heff*(nup-ndw)
  !      case (.false.)
  !         htmp = -(xmu+U/2d0)*(nup+ndw) + U*nup*ndw + heff*(nup-ndw)
  !      end select
  !      !energy of the bath=\sum_{n=1,N}\e_l n_l
  !      do kp=2,Ns
  !         htmp=htmp + eup(kp-1)*real(ib(kp),8)
  !         htmp=htmp + edw(kp-1)*real(ib(kp+Ns),8)
  !      enddo

  !      Hv(i) = Hv(i) + htmp*v(i)

  !      !Non-Diagonal part
  !      do ms=2,Ns
  !         if(ib(1)==1 .AND. ib(ms)==0)then
  !            call c(1,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
  !            call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
  !            j=invHmap(isector,r)
  !            tef=vup(ms-1)
  !            htmp = tef*sg1*sg2
  !            Hv(i) = Hv(i) + htmp*v(j)
  !            Hv(j) = Hv(j) + htmp*v(i)
  !         endif
  !         if(ib(1+Ns)==1 .AND. ib(ms+Ns)==0)then
  !            call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
  !            call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
  !            j=invHmap(isector,r)
  !            tef=vdw(ms-1)   
  !            htmp=tef*sg1*sg2
  !            Hv(i) = Hv(i) + htmp*v(j)
  !            Hv(j) = Hv(j) + htmp*v(i)
  !         endif
  !      enddo
  !   enddo
  !   return
  ! end subroutine HtimesV


end MODULE ED_GETH
