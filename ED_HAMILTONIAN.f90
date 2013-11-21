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
  public :: ed_buildH_d,ed_buildH_c

  !Sparse Matrix-vector product using stored sparse matrix 
  public :: spHtimesV_dd,spHtimesV_dc,spHtimesV_cc
  public :: lanc_spHtimesV_dd,lanc_spHtimesV_dc,lanc_spHtimesV_cc

  !Direct Matrix-vector product (no allocation of H)
  public :: setup_Hv_sector
  public :: delete_Hv_sector
  !   public :: HtimesV
  ! #ifdef _MPI
  !   public :: HtimesV_mpi
  ! #endif

  integer                      :: Hsector
  integer,dimension(:),pointer :: Hmap    !map of the Sector S to Hilbert space H

contains



  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE PRECISION
  !+------------------------------------------------------------------+
  subroutine ed_buildH_d(isector,h)
    real(8),optional,dimension(:,:)  :: h
    integer                          :: isector
    integer,dimension(Ntot)          :: ib
    integer                          :: mpiQ,mpiR                
    integer                          :: dim,iup,idw
    integer                          :: i,j,m,ms,iorb,jorb,ispin
    integer                          :: kp,k1,k2,k3,k4
    real(8)                          :: sg1,sg2,sg3,sg4
    real(8)                          :: htmp
    real(8),dimension(Norb)          :: nup,ndw,eloc
    logical                          :: Jcondition,flanc
    integer                          :: first_state,last_state

    !
    dim=getdim(isector)
    flanc=.true. ; if(present(h))flanc=.false.
    !
    first_state= 1
    last_state = dim
    !
    if(flanc)then
       if(spH0%status)call sp_delete_matrix(spH0) 
#ifdef _MPI
       mpiQ = dim/mpiSIZE
       mpiR = 0
       if(mpiID==(mpiSIZE-1))mpiR=mod(dim,mpiSIZE)
       call sp_init_matrix(spH0,mpiQ+mpiR)
       first_state= mpiID*mpiQ+1
       last_state = (mpiID+1)*mpiQ+mpiR
#else
       call sp_init_matrix(spH0,dim)
#endif
    else
       if(size(h,1)/=dim)stop "ED_GETH: wrong dimension 1 of H"
       if(size(h,2)/=dim)stop "ED_GETH: wrong dimension 2 of H"
       h=0.d0
    endif
    !
    !dump Hloc
    do iorb=1,Norb
       eloc(iorb)=dreal(Hloc(iorb,iorb))
    enddo
    !
    do i=first_state,last_state
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
       do iorb=1,size(dmft_bath%e,2)
          do kp=1,Nbath
             ! ms=Norb+(iorb-1)*Nbath + kp
             ! if(bath_type=='hybrid')ms=Norb+kp
             ms=getBathStride(iorb,kp)
             htmp =htmp + dmft_bath%e(1,iorb,kp)*real(ib(ms),8) + &
                  dmft_bath%e(Nspin,iorb,kp)*real(ib(ms+Ns),8)
          enddo
       enddo
       !
       !
       if(flanc)then
#ifdef _MPI
          call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,i)
#else
          call sp_insert_element(spH0,htmp,i,i)
#endif
       else
          h(i,i)=h(i,i)+htmp
       endif
       !
       !
       if(Norb>1.AND.Jhflag)then
          !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
          !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
          !it shoud rather be (not ordered product) but this changes sign in the code,
          !so it is more natura to NORMAL order the products of operators:
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
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   if(flanc)then
#ifdef _MPI
                      call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                      call sp_insert_element(spH0,htmp,i,j)
#endif
                   else
                      h(i,j)=h(i,j)+htmp
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
#ifdef _MPI
                      call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                      call sp_insert_element(spH0,htmp,i,j)
#endif
                   else
                      h(i,j)=h(i,j)+htmp
                   endif
                endif
             enddo
          enddo
       endif
       !

       !LOCAL HYBRIDIZATION
       do iorb=1,Norb
          do jorb=1,Norb
             !SPIN UP
             if((ib(iorb)==0).AND.(ib(jorb)==1))then
                call c(jorb,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(iorb,jorb)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !SPIN DW
             if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
                call c(jorb+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(iorb,jorb)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
          enddo
       enddo

       !IMP-BATH HYBRIDIZATION
       do iorb=1,Norb
          do kp=1,Nbath
             ms=getBathStride(iorb,kp)
             if(ib(iorb) == 1 .AND. ib(ms) == 0)then
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j = binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !
             if(ib(iorb) == 0 .AND. ib(ms) == 1)then
                call c(ms,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !
             if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !
             if(ib(iorb+Ns) == 0 .AND. ib(ms+Ns) == 1)then
                call c(ms+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
          enddo
       enddo
    enddo
  end subroutine ed_buildH_d








  !+------------------------------------------------------------------+
  !PURPOSE  : Build Hamiltonian sparse matrix DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine ed_buildH_c(isector,h)
    complex(8),optional,dimension(:,:) :: h
    integer                            :: isector
    integer,dimension(Ntot)            :: ib
    integer                            :: mpiQ,mpiR                
    integer                            :: dim,iup,idw
    integer                            :: i,j,m,ms,iorb,jorb,ispin
    integer                            :: kp,k1,k2,k3,k4
    real(8)                            :: sg1,sg2,sg3,sg4
    complex(8)                         :: htmp
    real(8),dimension(Norb)            :: nup,ndw,eloc
    logical                            :: Jcondition,flanc
    integer                            :: first_state,last_state
    !
    dim=getdim(isector)
    flanc=.true. ; if(present(h))flanc=.false.
    !
    first_state= 1
    last_state = dim
    !
    if(flanc)then
       if(spH0%status)call sp_delete_matrix(spH0) 
#ifdef _MPI
       mpiQ = dim/mpiSIZE
       mpiR = 0
       if(mpiID==(mpiSIZE-1))mpiR=mod(dim,mpiSIZE)
       call sp_init_matrix(spH0,mpiQ+mpiR)
       first_state= mpiID*mpiQ+1
       last_state = (mpiID+1)*mpiQ+mpiR
#else
       call sp_init_matrix(spH0,dim)
#endif
    else
       if(size(h,1)/=dim)stop "ED_GETH: wrong dimension 1 of H"
       if(size(h,2)/=dim)stop "ED_GETH: wrong dimension 2 of H"
       h=0.d0
    endif
    !
    !dump Hloc
    do iorb=1,Norb
       eloc(iorb)=dreal(Hloc(iorb,iorb))
    enddo
    !
    do i=first_state,last_state
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
       do iorb=1,size(dmft_bath%e,2)
          do kp=1,Nbath
             ! ms=Norb+(iorb-1)*Nbath + kp
             ! if(bath_type=='hybrid')ms=Norb+kp
             ms=getBathStride(iorb,kp)
             htmp =htmp + dmft_bath%e(1,iorb,kp)*real(ib(ms),8) + &
                  dmft_bath%e(Nspin,iorb,kp)*real(ib(ms+Ns),8)
          enddo
       enddo
       !
       !
       if(flanc)then
#ifdef _MPI
          call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,i)
#else
          call sp_insert_element(spH0,htmp,i,i)
#endif
       else
          h(i,i)=h(i,i)+htmp
       endif
       !
       !
       if(Norb>1.AND.Jhflag)then
          !SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
          !S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
          !S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
          !it shoud rather be (not ordered product) but this changes sign in the code,
          !so it is more natura to NORMAL order the products of operators:
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
                   j=binary_search(Hmap,k4)
                   htmp = Jh*sg1*sg2*sg3*sg4
                   if(flanc)then
#ifdef _MPI
                      call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                      call sp_insert_element(spH0,htmp,i,j)
#endif
                   else
                      h(i,j)=h(i,j)+htmp
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
#ifdef _MPI
                      call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                      call sp_insert_element(spH0,htmp,i,j)
#endif
                   else
                      h(i,j)=h(i,j)+htmp
                   endif
                endif
             enddo
          enddo
       endif
       !

       !LOCAL HYBRIDIZATION
       do iorb=1,Norb
          do jorb=1,Norb
             !SPIN UP
             if((ib(iorb)==0).AND.(ib(jorb)==1))then
                call c(jorb,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(iorb,jorb)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !SPIN DW
             if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
                call c(jorb+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = Hloc(iorb,jorb)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
          enddo
       enddo

       !IMP-BATH HYBRIDIZATION
       do iorb=1,Norb
          do kp=1,Nbath
             ms=getBathStride(iorb,kp)
             if(ib(iorb) == 1 .AND. ib(ms) == 0)then
                call c(iorb,m,k1,sg1)
                call cdg(ms,k1,k2,sg2)
                j = binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !
             if(ib(iorb) == 0 .AND. ib(ms) == 1)then
                call c(ms,m,k1,sg1)
                call cdg(iorb,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !
             if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
                call c(iorb+Ns,m,k1,sg1)
                call cdg(ms+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
             !
             if(ib(iorb+Ns) == 0 .AND. ib(ms+Ns) == 1)then
                call c(ms+Ns,m,k1,sg1)
                call cdg(iorb+Ns,k1,k2,sg2)
                j=binary_search(Hmap,k2)
                htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2
                if(flanc)then
#ifdef _MPI
                   call sp_insert_element(spH0,htmp,i-mpiID*mpiQ,j)
#else
                   call sp_insert_element(spH0,htmp,i,j)
#endif
                else
                   h(i,j)=h(i,j)+htmp
                endif
             endif
          enddo
       enddo
    enddo
  end subroutine ed_buildH_c









  !               RELATED COMPUTATIONAL ROUTINES
  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! Lanczos algorithm using serial double real, complex, and MPI
  !+------------------------------------------------------------------+
  subroutine spHtimesV_dd(N,Nloc,v,Hv)
    integer                    :: N,Nloc
    real(8),dimension(Nloc)    :: v
    real(8),dimension(Nloc)    :: Hv
    integer                    :: Q,R
    real(8),dimension(N)       :: vin,vtmp
    integer                    :: i
#ifdef _MPI
    Q = N/mpiSIZE ; R = 0
    if(mpiID==(mpiSIZE-1))R=mod(N,mpiSIZE)
    vtmp=0.d0
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       vtmp(i)=v(i-mpiID*Q)
    enddo
    call MPI_AllReduce(vtmp,vin,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,mpiErr)
    call sp_matrix_vector_product_mpi_dd(spH0,Q,R,N,vin,Nloc,Hv)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_dd

  subroutine spHtimesV_dc(N,Nloc,v,Hv)
    integer                    :: N,Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    integer                    :: Q,R
    complex(8),dimension(N)    :: vin,vtmp
    integer                    :: i
#ifdef _MPI
    Q = N/mpiSIZE ; R = 0
    if(mpiID==(mpiSIZE-1))R=mod(N,mpiSIZE)
    vtmp=zero
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       vtmp(i)=v(i-mpiID*Q)
    enddo
    call MPI_AllReduce(vtmp,vin,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,mpiErr)
    call sp_matrix_vector_product_mpi_dc(spH0,Q,R,N,vin,Nloc,Hv)
#else
    Hv=zero
    call sp_matrix_vector_product_dc(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_dc

  subroutine spHtimesV_cc(N,Nloc,v,Hv)
    integer                    :: N,Nloc
    complex(8),dimension(Nloc) :: v
    complex(8),dimension(Nloc) :: Hv
    integer                    :: Q,R
    complex(8),dimension(N)    :: vin,vtmp
    integer                    :: i
#ifdef _MPI
    Q = N/mpiSIZE ; R = 0
    if(mpiID==(mpiSIZE-1))R=mod(N,mpiSIZE)
    vtmp=zero
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       vtmp(i)=v(i-mpiID*Q)
    enddo
    call MPI_AllReduce(vtmp,vin,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,mpiErr)
    call sp_matrix_vector_product_mpi_cc(spH0,Q,R,N,vin,Nloc,Hv)
#else
    Hv=zero
    call sp_matrix_vector_product_cc(spH0,Nloc,v,Hv)
#endif
  end subroutine spHtimesV_cc



  !+------------------------------------------------------------------+
  !PURPOSE  : Perform the matrix-vector product H*v used in the
  ! Plain Lanczos algorithm for GF using MPI
  !NOTE that the integer arguments are inverted here with respect 
  ! to previous routines to respect the structure of the hamiltonian
  ! and the P_ARPACK algorithm using smaller block arrays (Nloc).
  !+------------------------------------------------------------------+
  subroutine lanc_spHtimesV_dd(Nloc,N,v,Hv)
    integer                 :: Nloc,N
    real(8),dimension(N)    :: v,Hv,Hvtmp
    real(8),dimension(Nloc) :: vout
    integer                 :: Q,R
    integer                 :: i
#ifdef _MPI
    Q = N/mpiSIZE ; R = 0
    if(mpiID==(mpiSIZE-1))R=mod(N,mpiSIZE)
    call sp_matrix_vector_product_mpi_dd(spH0,Q,R,N,v,Nloc,vout)
    Hvtmp=0.d0
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       Hvtmp(i)=vout(i-mpiID*Q)
    enddo
    Hv=0.d0
    call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Precision,MPI_Sum,MPI_Comm_World,mpiErr)
#else
    Hv=0.d0
    call sp_matrix_vector_product_dd(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_dd

  subroutine lanc_spHtimesV_dc(Nloc,N,v,Hv)
    integer                 :: Nloc,N
    complex(8),dimension(N)    :: v,Hv,Hvtmp
    complex(8),dimension(Nloc) :: vout
    integer                 :: Q,R
    integer                 :: i
#ifdef _MPI
    Q = N/mpiSIZE ; R = 0
    if(mpiID==(mpiSIZE-1))R=mod(N,mpiSIZE)
    call sp_matrix_vector_product_mpi_dc(spH0,Q,R,N,v,Nloc,vout)
    Hvtmp=0.d0
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       Hvtmp(i)=vout(i-mpiID*Q)
    enddo
    Hv=zero
    call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,mpiErr)
#else
    Hv=zero
    call sp_matrix_vector_product_dc(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_dc

  subroutine lanc_spHtimesV_cc(Nloc,N,v,Hv)
    integer                 :: Nloc,N
    complex(8),dimension(N)    :: v,Hv,Hvtmp
    complex(8),dimension(Nloc) :: vout
    integer                 :: Q,R
    integer                 :: i
#ifdef _MPI
    Q = N/mpiSIZE ; R = 0
    if(mpiID==(mpiSIZE-1))R=mod(N,mpiSIZE)
    call sp_matrix_vector_product_mpi_cc(spH0,Q,R,N,v,Nloc,vout)
    Hvtmp=0.d0
    do i=mpiID*Q+1,(mpiID+1)*Q+R
       Hvtmp(i)=vout(i-mpiID*Q)
    enddo
    Hv=zero
    call MPI_AllReduce(Hvtmp,Hv,N,MPI_Double_Complex,MPI_Sum,MPI_Comm_World,mpiErr)
#else
    Hv=zero
    call sp_matrix_vector_product_cc(spH0,N,v,Hv)
#endif
  end subroutine lanc_spHtimesV_cc


  !+------------------------------------------------------------------+
  !PURPOSE  : Direct Matrix-vector multiplication H*v used in 
  !Lanczos algorithm. this DOES NOT store the H-matrix (slower but 
  !more memory efficient)
  !+------------------------------------------------------------------+
  !include "ed_htimesv_direct.f90"
#ifdef _MPI
  !include "ed_htimesv_direct_mpi.f90"
#endif

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
