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

  !Get sparse sector Hamiltonian
  public :: ed_geth

  !Sparse Matrix-vector product using stored sparse matrix 
  public :: spHtimesV_d,spHtimesV_c
#ifdef _MPI
  public :: spHtimesV_mpi 
#endif

  !Direct Matrix-vector product (no allocation of H)
  public :: setup_Hv_sector
  public :: delete_Hv_sector
  public :: HtimesV
#ifdef _MPI
  public :: HtimesV_mpi
#endif


  integer                      :: Hsector
  integer,dimension(:),pointer :: Hmap    !map of the Sector S to Hilbert space H

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_geth(isector,h)
    real(8),optional,dimension(:,:)  :: h
    integer                          :: isector
    integer,dimension(Ntot)          :: ib
    integer                          :: dim,iup,idw
    integer                          :: i,j,m,ms,iorb,jorb,ispin
    integer                          :: kp,k1,k2,k3,k4
    real(8)                          :: sg1,sg2,sg3,sg4
    real(8)                          :: tef,htmp
    real(8),dimension(Norb,Nbath)    :: eup,edw
    real(8),dimension(Norb,Nbath)    :: vup,vdw
    real(8),dimension(Norb)          :: nup,ndw
    logical                          :: Jcondition,flanc

    dim=getdim(isector)
    call setup_Hv_sector(isector)
    flanc=.true. ; if(present(h))flanc=.false.

    if(flanc)then
       if(spH0%status)then
#ifdef _MPI
          if(mpiID==0)then
#endif
             print*,"ED_GETH: spH0 was already initialized in sector:"//txtfy(isector)
#ifdef _MPI
          endif
#endif
          call sp_delete_matrix(spH0) 
       endif
       call sp_init_matrix(spH0,dim)
    else
       if(size(h,1)/=dim)stop "ED_GETH: wrong dimension 1 of H"
       if(size(h,2)/=dim)stop "ED_GETH: wrong dimension 2 of H"
       h=0.d0
    endif
    !
    eup=ebath(1,:,:)   ; edw=ebath(Nspin,:,:)
    vup=vbath(1,:,:)   ; vdw=vbath(Nspin,:,:)
    !
    do i=1,dim
       m=Hmap(i)
       call bdecomp(m,ib)
       htmp=0.d0
       do iorb=1,Norb
          nup(iorb)=real(ib(iorb),8)
          ndw(iorb)=real(ib(iorb+Ns),8)
       enddo
       !
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
       !
       !Hbath: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
       do iorb=1,Norb
          do kp=1,Nbath
             ms=Norb+(iorb-1)*Nbath + kp
             htmp =htmp + eup(iorb,kp)*real(ib(ms),8) + edw(iorb,kp)*real(ib(ms+Ns),8)
          enddo
       enddo
       !
       !
       if(flanc)then
          call sp_insert_element(spH0,htmp,i,i)
       else
          h(i,i)=h(i,i)+htmp
       endif
       !
       !
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
                   call c(jorb,m,k1,sg1)
                   call c(iorb+Ns,k1,k2,sg2)
                   call cdg(jorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   ! call c(jorb,m,k1,sg1)
                   ! call cdg(jorb+Ns,k1,k2,sg2)
                   ! call c(iorb+Ns,k2,k3,sg3)
                   ! call cdg(iorb,k3,k4,sg4)
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
                   call c(jorb,m,k1,sg1)
                   call c(jorb+Ns,k1,k2,sg2)
                   call cdg(iorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
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
       !
       !
       !NON-LOCAL PART
       do iorb=1,Norb
          do kp=1,Nbath!Norb+1,Ns
             ms=Norb+(iorb-1)*Nbath + kp
             !UP
             if(ib(iorb) == 1 .AND. ib(ms) == 0)then
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j=binary_search(Hmap,k2)
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
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
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
    call delete_Hv_sector()
  end subroutine ed_geth



  !####################################################################
  !               RELATED COMPUTATIONAL ROUTINES
  !####################################################################
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

#ifdef _MPI
  subroutine spHtimesV_mpi(Q,R,Nloc,N,v,Hv)
    integer                 :: Q,R,Nloc,N
    real(8),dimension(Nloc) :: v,Hv
    real(8),dimension(N)    :: vin,vout
    integer                 :: i,j
    vout=0.d0
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       vout(i)=v(i-mpiID*Q)
    enddo
    call MPI_AllReduce(vout,vin,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,mpiErr)
    vout=0.d0
    call sp_matrix_vector_product_mpi(Q,R,N,spH0,vin,vout)
    Hv=0.d0
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       Hv(i-mpiID*Q)=vout(i)
    enddo
  end subroutine SpHtimesV_mpi
#endif








  !+------------------------------------------------------------------+
  !PURPOSE  : Direct Matrix-vector multiplication H*v used in 
  !Lanczos algorithm. this DOES NOT store the H-matrix (slower but 
  !more memory efficient)
  !+------------------------------------------------------------------+
  subroutine HtimesV(Nv,v,Hv)
    integer                       :: Nv
    real(8),dimension(Nv)         :: v
    real(8),dimension(Nv)         :: Hv
    integer                       :: isector
    integer,dimension(Ntot)       :: ib
    integer                       :: dim,iup,idw
    integer                       :: i,j,m,ms,iorb,jorb,ispin
    integer                       :: kp,k1,k2,k3,k4
    real(8)                       :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb,Nbath) :: eup,edw
    real(8),dimension(Norb,Nbath) :: vup,vdw
    real(8),dimension(Norb)       :: nup,ndw
    real(8)                       :: tef,htmp
    logical                       :: Jcondition,flanc
    isector=Hsector
    dim=getdim(isector)
    if(.not.associated(Hmap).AND.size(Hmap)/=dim)stop "HtimesV: wrong allocation of Hmap"
    !
    eup=ebath(1,:,:)   ; edw=ebath(Nspin,:,:)
    vup=vbath(1,:,:)   ; vdw=vbath(Nspin,:,:)
    !
    if(Nv/=dim)stop "HtimesV error in dimensions"
    Hv=0.d0
    do i=1,dim
       m=Hmap(i)
       call bdecomp(m,ib)
       do iorb=1,Norb
          nup(iorb)=real(ib(iorb),8)
          ndw(iorb)=real(ib(iorb+Ns),8)
       enddo
       !
       !Diagonal part
       !local part of the impurity Hamiltonian: (-mu+\e0)*n + U*(n_up-0.5)*(n_dw-0.5) + heff*mag
       !+ energy of the bath=\sum_{n=1,N}\e_l n_l
       htmp=0.d0
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
       !
       Hv(i) = Hv(i) + htmp*v(i)
       !
       !
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
                   call c(jorb,m,k1,sg1)
                   call c(iorb+Ns,k1,k2,sg2)
                   call cdg(jorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   ! call c(jorb,m,k1,sg1)
                   ! call cdg(jorb+Ns,k1,k2,sg2)
                   ! call c(iorb+Ns,k2,k3,sg3)
                   ! call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   !
                   Hv(i) = Hv(i) + htmp*v(j)
                   Hv(j) = Hv(j) + htmp*v(i)
                   !
                endif
             enddo
          enddo
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
                   call c(jorb,m,k1,sg1)
                   call c(jorb+Ns,k1,k2,sg2)
                   call cdg(iorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   !
                   Hv(i) = Hv(i) + htmp*v(j)
                   Hv(j) = Hv(j) + htmp*v(i)
                   !
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
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                tef=vup(iorb,kp)
                htmp = tef*sg1*sg2
                !
                Hv(i) = Hv(i) + htmp*v(j)
                Hv(j) = Hv(j) + htmp*v(i)
                !
             endif
             !DW
             if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                tef=vdw(iorb,kp)
                htmp=tef*sg1*sg2
                !
                Hv(i) = Hv(i) + htmp*v(j)
                Hv(j) = Hv(j) + htmp*v(i)
                !
             endif
          enddo
       enddo
    enddo
  end subroutine HtimesV

#ifdef _MPI
  subroutine HtimesV_mpi(Nchunk,NRest,Nloc,Nv,v,Hv)
    integer                       :: Nchunk,Nrest
    integer                       :: Nloc !the small dimension (LDV in parpack)
    integer                       :: Nv   !the large dimension (Ns in parpack)
    real(8),dimension(Nloc)       :: v    !this required by parpack and is small
    real(8),dimension(Nloc)       :: Hv   !this required by parpack and is small
    real(8),dimension(Nv)         :: vin       !this used here and is large
    real(8),dimension(Nv)         :: vtmp !
    integer                       :: isector
    integer,dimension(Ntot)       :: ib
    integer                       :: dim,iup,idw
    integer                       :: i,j,k,m,ms,iorb,jorb,ispin
    integer                       :: kp,k1,k2,k3,k4
    real(8)                       :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb,Nbath) :: eup,edw
    real(8),dimension(Norb,Nbath) :: vup,vdw
    real(8),dimension(Norb)       :: nup,ndw
    real(8)                       :: tef,htmp
    logical                       :: Jcondition,flanc
    isector=Hsector
    dim=getdim(isector)
    if(.not.associated(Hmap).AND.size(Hmap)/=dim)stop "HtimesV: wrong allocation of Hmap"
    !
    eup=ebath(1,:,:)   ; edw=ebath(Nspin,:,:)
    vup=vbath(1,:,:)   ; vdw=vbath(Nspin,:,:)
    !
    if(Nv/=dim)stop "HtimesV error in dimensions"

    !each processor dump its small piece of vector v(i) into the large vtmp vector
    !that is reduced so to have it shared among nodes
    vtmp=0.d0
    do i=mpiID*Nchunk+1,(mpiID+1)*Nchunk+Nrest
       vtmp(i)=v(i-mpiID*Nchunk)
    enddo
    call MPI_AllReduce(vtmp,vin,Nv,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,mpiErr)

    !each node perform a part of the matrix vector product and store it in vtmp
    vtmp=0.d0
    do i=mpiID*Nchunk+1,(mpiID+1)*Nchunk+Nrest
       m=Hmap(i)
       call bdecomp(m,ib)
       do iorb=1,Norb
          nup(iorb)=real(ib(iorb),8)
          ndw(iorb)=real(ib(iorb+Ns),8)
       enddo
       !
       !Diagonal part
       !local part of the impurity Hamiltonian: (-mu+\e0)*n + U*(n_up-0.5)*(n_dw-0.5) + heff*mag
       !+ energy of the bath=\sum_{n=1,N}\e_l n_l
       htmp=0.d0
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
       !
       vtmp(i) = vtmp(i) + htmp*vin(i)
       !
       !
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
                   call c(jorb,m,k1,sg1)
                   call c(iorb+Ns,k1,k2,sg2)
                   call cdg(jorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   ! call c(jorb,m,k1,sg1)
                   ! call cdg(jorb+Ns,k1,k2,sg2)
                   ! call c(iorb+Ns,k2,k3,sg3)
                   ! call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   !
                   vtmp(i) = vtmp(i) + htmp*vin(j)
                   vtmp(j) = vtmp(j) + htmp*vin(i)
                   !
                endif
             enddo
          enddo
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
                   call c(jorb,m,k1,sg1)
                   call c(jorb+Ns,k1,k2,sg2)
                   call cdg(iorb+Ns,k2,k3,sg3)
                   call cdg(iorb,k3,k4,sg4)
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   !
                   vtmp(i) = vtmp(i) + htmp*vin(j)
                   vtmp(j) = vtmp(j) + htmp*vin(i)
                   !
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
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                tef=vup(iorb,kp)
                htmp = tef*sg1*sg2
                !
                vtmp(i) = vtmp(i) + htmp*vin(j)
                vtmp(j) = vtmp(j) + htmp*vin(i)
                !
             endif
             !DW
             if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                tef=vdw(iorb,kp)
                htmp=tef*sg1*sg2
                !
                vtmp(i) = vtmp(i) + htmp*vin(j)
                vtmp(j) = vtmp(j) + htmp*vin(i)
                !
             endif
          enddo
       enddo
    enddo

    ! !the tmp array vtmp is now reduced to all nodes
    ! vout=0.d0
    ! call MPI_ALLREDUCE(vtmp,vout,Nloc,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpiERR)

    !and each piece is dumped back to the each copy of the small vector on each node
    Hv=0.d0
    do i=mpiID*Nchunk+1,(mpiID+1)*Nchunk+Nrest
       Hv(i-mpiID*Nchunk)=vtmp(i)
    enddo

    !note I am not sure that the last ALLREDUCE is strictly necessary, because
    !each node should have its piece of vector that must be copied to Hv...
  end subroutine HtimesV_mpi
#endif



  !####################################################################
  !               RELATED COMPUTATIONAL ROUTINES
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE : 
  !+------------------------------------------------------------------+
  subroutine  setup_Hv_sector(isector)
    integer                              :: isector
    integer                              :: dim
    Hsector=isector
    dim = getdim(Hsector)
    allocate(Hmap(dim))
    call build_sector(isector,Hmap)
  end subroutine setup_Hv_sector

  subroutine delete_Hv_sector()
    deallocate(Hmap)
  end subroutine delete_Hv_sector



end MODULE ED_HAMILTONIAN
