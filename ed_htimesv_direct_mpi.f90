subroutine HtimesV_mpi(Nchunk,NRest,Nloc,Nv,v,Hv)
  integer                       :: Nchunk,Nrest
  integer                       :: Nloc !the small dimension (LDV in parpack)
  integer                       :: Nv   !the large dimension (Ns in parpack)
  real(8),dimension(Nloc)       :: v    !this required by parpack and is small
  real(8),dimension(Nloc)       :: Hv   !this required by parpack and is small
  real(8),dimension(Nv)         :: vin       !this used here and is large
  real(8),dimension(Nv)         :: vtmp
  integer                       :: isector
  integer,dimension(Ntot)       :: ib
  integer                       :: dim,iup,idw
  integer                       :: i,j,m,ms,iorb,jorb,ispin
  integer                       :: kp,k1,k2,k3,k4
  real(8)                       :: sg1,sg2,sg3,sg4
  real(8),dimension(Norb)       :: nup,ndw
  real(8)                       :: tef,htmp
  logical                       :: Jcondition,flanc
  isector=Hsector
  dim=getdim(isector)
  if(.not.associated(Hmap).AND.size(Hmap)/=dim)stop "HtimesV: wrong allocation of Hmap"
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
     do iorb=1,size(dmft_bath%e,2)
        do kp=1,Nbath
           ms=Norb+(iorb-1)*Nbath + kp
           htmp =htmp + dmft_bath%e(1,iorb,kp)*real(ib(ms),8) + dmft_bath%e(Nspin,iorb,kp)*real(ib(ms+Ns),8)
           !htmp =htmp + eup(iorb,kp)*real(ib(ms),8) + edw(iorb,kp)*real(ib(ms+Ns),8)
        enddo
     enddo
     !
     vtmp(i) = vtmp(i) + htmp*vin(i)
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
                 j=binary_search(Hmap,k4)
                 htmp = Jh*sg1*sg2*sg3*sg4
                 !
                 vtmp(i) = vtmp(i) + htmp*vin(j)
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
                 !
              endif
           enddo
        enddo
     endif
     !NON-LOCAL PART
     do iorb=1,Norb
        do kp=1,Nbath
           ms=Norb+(iorb-1)*Nbath + kp
           !UP
           if(ib(iorb) == 1 .AND. ib(ms) == 0)then
              call c(iorb,m,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2 !vup(iorb,kp)
              vtmp(i) = vtmp(i) + htmp*vin(j)
           endif
           !
           if(ib(iorb) == 0 .AND. ib(ms) == 1)then
              call c(ms,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp = dmft_bath%v(1,iorb,kp)*sg1*sg2 !vup(iorb,kp)
              vtmp(i) = vtmp(i) + htmp*vin(j)
           endif
           !DW
           if(ib(iorb+Ns) == 1 .AND. ib(ms+Ns) == 0)then
              call c(iorb+Ns,m,k1,sg1)
              call cdg(ms+Ns,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2!vdw(iorb,kp)
              vtmp(i) = vtmp(i) + htmp*vin(j)
           endif
           !
           if(ib(iorb+Ns) == 0 .AND. ib(ms+Ns) == 1)then
              call c(ms+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(Hmap,k2)
              htmp=dmft_bath%v(Nspin,iorb,kp)*sg1*sg2!vdw(iorb,kp)
              vtmp(i) = vtmp(i) + htmp*vin(j)
           endif
        enddo
     enddo
  enddo

  ! !the tmp array vtmp is now reduced to all nodes
  !vout=0.d0
  !call MPI_ALLREDUCE(vtmp,vout,Nloc,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpiERR)

  !and each piece is dumped back to the each copy of the small vector on each node
  Hv=0.d0
  do i=mpiID*Nchunk+1,(mpiID+1)*Nchunk+Nrest
     Hv(i-mpiID*Nchunk)=vtmp(i)!vout(i)
  enddo

  !note I am not sure that the last ALLREDUCE is strictly necessary, because
  !each node should have its piece of vector that must be copied to Hv...
end subroutine HtimesV_mpi


