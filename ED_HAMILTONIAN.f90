!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none
  private
  public :: ed_geth
  public :: spHtimesV_d,spHtimesV_c
  !!<MPI
  !public :: spHtimesV_mpi 
  !!>MPI
  public :: set_Hsector

  integer :: Hsector

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_geth(isector,h)
    real(8),optional,dimension(:,:)  :: h
    integer                          :: isector
    integer,dimension(Ntot)          :: ib
    integer                          :: dim,iup,idw
    integer                          :: i,j,k,r,m,ms,iorb,jorb,ispin
    integer                          :: kp,k1,k2,k3,k4
    real(8)                          :: sg1,sg2,sg3,sg4
    real(8)                          :: tef,htmp
    real(8),dimension(Norb,Nbath)    :: eup,edw
    real(8),dimension(Norb,Nbath)    :: vup,vdw
    real(8),dimension(Norb)          :: nup,ndw
    logical                          :: Jcondition,flanc
    integer,allocatable,dimension(:) :: Hmap    !map of the Sector S to Hilbert space H


    dim=getdim(isector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)

    flanc=.true. ; if(present(h))flanc=.false.

    if(flanc)then
       if(spH0%status)then
          print*,"ED_GETH: spH0 was already initialized in sector:"//txtfy(isector)
          call sp_delete_matrix(spH0) 
       endif
       call sp_init_matrix(spH0,dim)
    else
       if(size(h,1)/=dim)stop "ED_GETH: wrong dimension 1 of H"
       if(size(h,2)/=dim)stop "ED_GETH: wrong dimension 2 of H"
       h=0.d0
    endif

    eup=ebath(1,:,:)   ; edw=ebath(Nspin,:,:)
    vup=vbath(1,:,:)   ; vdw=vbath(Nspin,:,:)

    do i=1,dim
       m=Hmap(i)
       call bdecomp(m,ib)

       htmp=0.d0

       do iorb=1,Norb
          nup(iorb)=real(ib(iorb),8)
          ndw(iorb)=real(ib(iorb+Ns),8)
       enddo

       !LOCAL HAMILTONIAN PART:
       htmp = -xmu*(sum(nup)+sum(ndw))  + dot_product(eloc,nup+ndw) !+ heff*(sum(nup)-sum(ndw))

       !Density-density interaction: same orbital, opposite spins
       htmp = htmp + dot_product(uloc,nup*ndw)!=\sum=i U_i*(n_u*n_d)_i
       if(hfmode)htmp=htmp - 0.5d0*dot_product(uloc,nup+ndw) + 0.25d0*sum(uloc)

       if(Norb>1)then
          !density-density interaction: different orbitals, opposite spins
          do iorb=1,Norb         ! n_up_i*n_dn_j
             do jorb=iorb+1,Norb ! n_up_j*n_dn_i
                htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
             enddo
          enddo

          !density-density interaction: different orbitals, parallel spins
          !Jhund effect: U``=U`-J smallest of the interactions
          do iorb=1,Norb         ! n_up_i*n_up_j
             do jorb=iorb+1,Norb ! n_dn_i*n_dn_j
                htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
             enddo
          enddo
       endif

       !Hbath: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
       do iorb=1,Norb
          do kp=1,Nbath
             ms=Norb+(iorb-1)*Nbath + kp
             htmp =htmp + eup(iorb,kp)*real(ib(ms),8) + edw(iorb,kp)*real(ib(ms+Ns),8)
          enddo
       enddo

       if(flanc)then
          call sp_insert_element(spH0,htmp,i,i)
       else
          h(i,i)=h(i,i)+htmp
       endif


       if(Norb>1.AND.Jhflag)then
          !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
          !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
          !it shoud rather be (not ordered product):
          !S-E: J c^+_iorb_up c_iorb_dw   c^+_jorb_dw    c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb}  c_{iorb+Ns} c^+_{jorb+Ns}  c_{jorb}
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ib(jorb)==1).AND.&
                     (ib(iorb+Ns)==1).AND.&
                     (ib(jorb+Ns)==0).AND.&
                     (ib(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1)      ;sg1=dfloat(k1)/dfloat(abs(k1));k1=abs(k1)
                   call c(iorb+Ns,k1,k2)  ;sg2=dfloat(k2)/dfloat(abs(k2));k2=abs(k2)
                   call cdg(jorb+Ns,k2,k3);sg3=dfloat(k3)/dfloat(abs(k3));k3=abs(k3)
                   call cdg(iorb,k3,k4)   ;sg4=dfloat(k4)/dfloat(abs(k4));k4=abs(k4)
                   ! call c(jorb,m,k1)      ;sg1=dfloat(k1)/dfloat(abs(k1));k1=abs(k1)
                   ! call cdg(jorb+Ns,k1,k2);sg2=dfloat(k2)/dfloat(abs(k2));k2=abs(k2)
                   ! call c(iorb+Ns,k2,k3)  ;sg3=dfloat(k3)/dfloat(abs(k3));k3=abs(k3)
                   ! call cdg(iorb,k3,k4)   ;sg4=dfloat(k4)/dfloat(abs(k4));k4=abs(k4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   if(flanc)then
                      call sp_insert_element(spH0,htmp,i,j)
                      call sp_insert_element(spH0,htmp,j,i)
                   else
                      h(i,j)=h(i,j)+htmp
                      h(j,i)=h(i,j)
                   endif
                endif
             enddo
          enddo
          !
          !PAIR-HOPPING (P-H) TERMS
          !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
          !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
          do iorb=1,Norb
             do jorb=1,Norb
                Jcondition=(&
                     (iorb/=jorb).AND.&
                     (ib(jorb)==1).AND.&
                     (ib(jorb+Ns)==1).AND.&
                     (ib(iorb+Ns)==0).AND.&
                     (ib(iorb)==0))
                if(Jcondition)then
                   call c(jorb,m,k1)      ;sg1=dfloat(k1)/dfloat(abs(k1));k1=abs(k1)
                   call c(jorb+Ns,k1,k2)  ;sg2=dfloat(k2)/dfloat(abs(k2));k2=abs(k2)
                   call cdg(iorb+Ns,k2,k3);sg3=dfloat(k3)/dfloat(abs(k3));k3=abs(k3)
                   call cdg(iorb,k3,k4)   ;sg4=dfloat(k4)/dfloat(abs(k4));k4=abs(k4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   if(flanc)then
                      call sp_insert_element(spH0,htmp,i,j)
                      call sp_insert_element(spH0,htmp,j,i)
                   else
                      h(i,j)=h(i,j)+htmp
                      h(j,i)=h(i,j)
                   endif
                endif
             enddo
          enddo
       endif


       !NON-LOCAL PART
       do iorb=1,Norb
          do kp=1,Nbath!Norb+1,Ns
             ms=Norb+(iorb-1)*Nbath + kp
             !UP
             if(ib(iorb) == 1 .AND. ib(ms) == 0)then
                call c(iorb,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
                j=binary_search(Hmap,r)
                tef=vup(iorb,kp)
                htmp = tef*sg1*sg2
                if(flanc)then
                   call sp_insert_element(spH0,htmp,i,j)
                   call sp_insert_element(spH0,htmp,j,i)
                else
                   h(i,j)=h(i,j)+htmp
                   h(j,i)=h(i,j)
                endif
             endif
             !DW
             if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(iorb+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
                call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
                j=binary_search(Hmap,r)
                tef=vdw(iorb,kp)
                htmp=tef*sg1*sg2
                if(flanc)then
                   call sp_insert_element(spH0,htmp,i,j)
                   call sp_insert_element(spH0,htmp,j,i)
                else
                   h(i,j)=h(i,j)+htmp
                   h(j,i)=h(i,j)
                endif
             endif
          enddo
       enddo
    enddo
    deallocate(Hmap)
  end subroutine ed_geth









  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine  set_Hsector(isector)
    integer :: isector
    Hsector=isector
  end subroutine set_Hsector


  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! Lanczos algorithm using serial double real, complex, and MPI (commented)
  !+------------------------------------------------------------------+
  subroutine spHtimesV_d(N,v,Hv)
    integer              :: N
    real(8),dimension(N) :: v
    real(8),dimension(N) :: Hv
    Hv=zero
    call sp_matrix_vector_product(N,spH0,v,Hv)
  end subroutine SpHtimesV_d
  !---------------------------------!
  subroutine spHtimesV_c(N,v,Hv)
    integer              :: N
    complex(8),dimension(N) :: v
    complex(8),dimension(N) :: Hv
    Hv=zero
    call sp_matrix_vector_product(N,spH0,v,Hv)
  end subroutine SpHtimesV_c

  !!<MPI
  ! subroutine spHtimesV_mpi(Q,R,Nloc,N,v,Hv)
  !   integer                 :: Q,R,Nloc,N
  !   real(8),dimension(Nloc) :: v,Hv
  !   real(8),dimension(N)    :: vin,vout
  !   integer                 :: i,j
  !   vout=0.d0
  !   do i=mpiID*Q+1,(mpiID+1)*Q+R
  !      vout(i)=v(i-mpiID*Q)
  !   enddo
  !   call MPI_AllReduce(vout,vin,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,mpiErr)
  !   vout=0.d0
  !   call sp_matrix_vector_product_mpi(Q,R,N,spH0,vin,vout)
  !   Hv=0.d0
  !   do i=mpiID*Q+1,(mpiID+1)*Q+R
  !      Hv(i-mpiID*Q)=vout(i)
  !   enddo
  ! end subroutine SpHtimesV_mpi
  !!>MPI


  ! !+------------------------------------------------------------------+
  ! !PURPOSE  : Direct Matrix-vector multiplication H*v used in 
  ! !Lanczos algorithm. this DOES NOT store the H-matrix (slower but 
  ! !more memory efficient)
  ! !+------------------------------------------------------------------+
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
  !      m=Hmap(i)
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
  !            j=invHmap(r)
  !            tef=vup(ms-1)
  !            htmp = tef*sg1*sg2
  !            Hv(i) = Hv(i) + htmp*v(j)
  !            Hv(j) = Hv(j) + htmp*v(i)
  !         endif
  !         if(ib(1+Ns)==1 .AND. ib(ms+Ns)==0)then
  !            call c(1+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
  !            call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
  !            j=invHmap(r)
  !            tef=vdw(ms-1)   
  !            htmp=tef*sg1*sg2
  !            Hv(i) = Hv(i) + htmp*v(j)
  !            Hv(j) = Hv(j) + htmp*v(i)
  !         endif
  !      enddo
  !   enddo
  !   return
  ! end subroutine HtimesV



  ! !####################################################################
  ! !                    FULL DIAGONALIZATION
  ! !####################################################################
  ! !+------------------------------------------------------------------+
  ! !PURPOSE  : 
  ! !+------------------------------------------------------------------+
  ! subroutine full_ed_geth(isector,h)
  !   real(8)                       :: h(:,:)
  !   integer                       :: ib(Ntot)
  !   integer                       :: dim
  !   integer                       :: i,j,k,r,m,ms,iorb,jorb,ispin,iup,idw
  !   integer                       :: kp,isector,k1,k2,k3,k4
  !   real(8)                       :: sg1,sg2,sg3,sg4
  !   real(8),dimension(Norb,Nbath) :: eup,edw
  !   real(8),dimension(Norb,Nbath) :: vup,vdw
  !   real(8)                       :: nup(Norb),ndw(Norb),No(Norb),Nno(Norb),tef,htmp
  !   logical                       :: Jcondition
  !   dim=getdim(isector)
  !   if(size(h,1)/=dim)call error("FULL_ED_GETH: wrong dimension 1 of H")
  !   if(size(h,2)/=dim)call error("FULL_ED_GETH: wrong dimension 2 of H")
  !   h=0.d0
  !   eup=ebath(1,:,:)   ; edw=ebath(Nspin,:,:)
  !   vup=vbath(1,:,:)   ; vdw=vbath(Nspin,:,:)
  !   do i=1,dim
  !      m=Hmap(i)
  !      call bdecomp(m,ib)
  !      htmp=0.d0
  !      do iorb=1,Norb
  !         nup(iorb)=real(ib(iorb),8)
  !         ndw(iorb)=real(ib(iorb+Ns),8)
  !      enddo
  !      !LOCAL HAMILTONIAN PART:
  !      htmp = -xmu*(sum(nup)+sum(ndw))  + dot_product(eloc,nup+ndw) !+ heff*(sum(nup)-sum(ndw))
  !      !+
  !      !Density-density interaction: same orbital, opposite spins
  !      htmp = htmp + dot_product(uloc,nup*ndw)!=\sum=i U_i*(n_u*n_d)_i
  !      if(hfmode)htmp=htmp - 0.5d0*dot_product(uloc,nup+ndw) + 0.25d0*sum(uloc)
  !      !+
  !      if(Norb>1)then
  !         !density-density interaction: different orbitals, opposite spins
  !         do iorb=1,Norb         ! n_up_i*n_dn_j
  !            do jorb=iorb+1,Norb ! n_up_j*n_dn_i
  !               htmp = htmp + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
  !            enddo
  !         enddo
  !         !density-density interaction: different orbitals, parallel spins
  !         !Jhund effect: U``=U`-J smallest of the interactions
  !         do iorb=1,Norb         ! n_up_i*n_up_j
  !            do jorb=iorb+1,Norb ! n_dn_i*n_dn_j
  !               htmp = htmp + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
  !            enddo
  !         enddo
  !      endif
  !      !+
  !      !Hbath: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
  !      do iorb=1,Norb
  !         do kp=1,Nbath
  !            ms=Norb+(iorb-1)*Nbath + kp
  !            htmp =htmp + eup(iorb,kp)*real(ib(ms),8) + edw(iorb,kp)*real(ib(ms+Ns),8)
  !         enddo
  !      enddo
  !      !=
  !      h(i,i)=htmp
  !      !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
  !      if(Jhflag.AND.Norb>1)then
  !         !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
  !         !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
  !         !it shoud rather be (not ordered product):
  !         !S-E: J c^+_iorb_up c_iorb_dw   c^+_jorb_dw    c_jorb_up  (i.ne.j) 
  !         !S-E: J c^+_{iorb}  c_{iorb+Ns} c^+_{jorb+Ns}  c_{jorb}
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               Jcondition=(&
  !                    (iorb/=jorb).AND.&
  !                    (ib(jorb)==1).AND.&
  !                    (ib(iorb+Ns)==1).AND.&
  !                    (ib(jorb+Ns)==0).AND.&
  !                    (ib(iorb)==0))
  !               if(Jcondition)then
  !                  call c(jorb,m,k1)      ;sg1=dfloat(k1)/dfloat(abs(k1));k1=abs(k1)
  !                  call c(iorb+Ns,k1,k2)  ;sg2=dfloat(k2)/dfloat(abs(k2));k2=abs(k2)
  !                  call cdg(jorb+Ns,k2,k3);sg3=dfloat(k3)/dfloat(abs(k3));k3=abs(k3)
  !                  call cdg(iorb,k3,k4)   ;sg4=dfloat(k4)/dfloat(abs(k4));k4=abs(k4)
  !                  ! call c(jorb,m,k1)      ;sg1=dfloat(k1)/dfloat(abs(k1));k1=abs(k1)
  !                  ! call cdg(jorb+Ns,k1,k2);sg2=dfloat(k2)/dfloat(abs(k2));k2=abs(k2)
  !                  ! call c(iorb+Ns,k2,k3)  ;sg3=dfloat(k3)/dfloat(abs(k3));k3=abs(k3)
  !                  ! call cdg(iorb,k3,k4)   ;sg4=dfloat(k4)/dfloat(abs(k4));k4=abs(k4)
  !                  j=invHmap(k4)
  !                  htmp = Jh*sg1*sg2*sg3*sg4
  !                  call sp_insert_element(spH0,htmp,i,j)
  !                  call sp_insert_element(spH0,htmp,j,i)
  !               endif
  !            enddo
  !         enddo
  !         !
  !         !PAIR-HOPPING (P-H) TERMS
  !         !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
  !         !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
  !         do iorb=1,Norb
  !            do jorb=1,Norb
  !               Jcondition=(&
  !                    (iorb/=jorb).AND.&
  !                    (ib(jorb)==1).AND.&
  !                    (ib(jorb+Ns)==1).AND.&
  !                    (ib(iorb+Ns)==0).AND.&
  !                    (ib(iorb)==0))
  !               if(Jcondition)then
  !                  call c(jorb,m,k1)      ;sg1=dfloat(k1)/dfloat(abs(k1));k1=abs(k1)
  !                  call c(jorb+Ns,k1,k2)  ;sg2=dfloat(k2)/dfloat(abs(k2));k2=abs(k2)
  !                  call cdg(iorb+Ns,k2,k3);sg3=dfloat(k3)/dfloat(abs(k3));k3=abs(k3)
  !                  call cdg(iorb,k3,k4)   ;sg4=dfloat(k4)/dfloat(abs(k4));k4=abs(k4)
  !                  j=invHmap(k4)
  !                  htmp = Jh*sg1*sg2*sg3*sg4
  !                  call sp_insert_element(spH0,htmp,i,j)
  !                  call sp_insert_element(spH0,htmp,j,i)
  !               endif
  !            enddo
  !         enddo
  !      endif
  !      !NON-LOCAL PART:
  !      do iorb=1,Norb
  !         do kp=1,Nbath!Norb+1,Ns
  !            ms=Norb+(iorb-1)*Nbath + kp
  !            !UP
  !            if(ib(iorb) == 1 .AND. ib(ms) == 0)then
  !               call c(iorb,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
  !               call cdg(ms,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)
  !               j=invHmap(r)
  !               tef=vup(iorb,kp)!ms-Norb)
  !               h(i,j)=tef*sg1*sg2
  !               h(j,i)=h(i,j)
  !            endif
  !            !DW
  !            if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
  !               call c(iorb+Ns,m,k);sg1=dfloat(k)/dfloat(abs(k));k=abs(k)
  !               call cdg(ms+Ns,k,r);sg2=dfloat(r)/dfloat(abs(r));r=abs(r)           
  !               j=invHmap(r)
  !               tef=vdw(iorb,kp)!ms-Norb)
  !               h(i,j)=tef*sg1*sg2
  !               h(j,i)=h(i,j)
  !            endif
  !         enddo
  !      enddo
  !   enddo
  !   return
  ! end subroutine full_ed_geth

end MODULE ED_HAMILTONIAN
