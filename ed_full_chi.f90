!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_getchi()
  real(8)                               :: cdgmat(Norb),chij(Norb),chitot,spin,spintot
  integer,dimension(N)                  :: ib(Ntot)
  integer                               :: i,j,k,r,ll,m,in,is,ispin,iorb,jorb,isector
  integer                               :: idim,ia,unit(3)
  real(8)                               :: Ei,Ej,cc,peso(Norb),pesotot
  real(8)                               :: expterm,de,w0,it
  complex(8)                            :: iw
  integer,allocatable,dimension(:)      :: HImap    !map of the Sector S to Hilbert space H
  write(LOGfile,"(A)")"Evaluating Suceptibility:"
  call allocate_grids
  allocate(Chitau(Norb,0:Ltau),Chiw(Norb,Lreal),Chiiw(Norb,0:Lmats))
  Chitau=0.d0
  Chiw=zero
  Chiiw=zero
  !Spin susceptibility \X(tau). |<i|S_z|j>|^2
  write(LOGfile,"(A)")"Evaluating Chi_Sz"
  call start_progress(LOGfile)

  do isector=1,Nsect !loop over <i| total particle number
     call progress(isector,Nsect)
     idim=getdim(isector)
     allocate(HImap(idim))
     call build_sector(isector,HImap)
     do i=1,idim 
        do j=1,idim
           chij=0.d0;chitot=0.d0
           expterm=exp(-beta*espace(isector)%e(j))
           if(expterm<cutoff)cycle
           do ll=1,idim 
              ia=HImap(ll)
              call bdecomp(ia,ib)
              spintot=0.d0
              do iorb=1,Norb
                 spin=real(ib(iorb),8)-real(ib(iorb+Ns),8) !nup - ndw
                 chij(iorb)=chij(iorb)+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
              enddo
              spintot=spintot+real(ib(iorb),8)-real(ib(iorb+Ns),8)
              chitot=chitot+espace(isector)%M(ll,i)*spintot*espace(isector)%M(ll,j)
           enddo
           Ei=espace(isector)%e(i)
           Ej=espace(isector)%e(j)
           de=Ei-Ej
           !ACTHUNG!! chij should not be squared 'cause is <j|Sz|i> and not <i|c+|j><j|c|i>
           peso=chij/zeta_function !chij**2/zeta_function
           pesotot=chitot/zeta_function !chitot**2/zeta_function
           do iorb=1,Norb
              !Matsubara (bosonic) frequency
              if(de>cutoff)chiiw(iorb,0)=chiiw(iorb,0)-peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
              do m=1,Lmats
                 iw=xi*vm(m)
                 chiiw(iorb,m)=chiiw(iorb,m)+peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
              enddo
              !Real-frequency
              do m=1,Lreal
                 w0=wr(m);iw=cmplx(w0,eps,8)
                 !Retarded = Commutator = response function
                 chiw(iorb,m)=chiw(iorb,m)+peso(iorb)*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
                 !Time-ordered
                 ! chiw(iorb,m)=chiw(iorb,m)-exp(-beta*espace(isector)%e(j))*&
                 !      (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso(iorb)
              enddo

              !Imaginary time:
              do m=0,Ltau 
                 it=tau(m)
                 chitau(iorb,m)=chitau(iorb,m) + exp(-it*espace(isector)%e(i))*&
                      exp(-(beta-it)*espace(isector)%e(j))*peso(iorb)
                 ! chitau(iorb,m)=chitau(iorb,m) + peso(iorb)*exp(-beta*Ej)*exp(-it*de)
              enddo
           enddo
           ! !Matsubara (bosonic) frequency
           ! if(de>cutoff)chiiwtot(0)=chiiwtot(0)-pesotot*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
           ! do m=1,Lmats
           !    iw=xi*vm(m)
           !    chiiwtot(m)=chiiwtot(m)+pesotot*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
           ! enddo
           ! do m=0,Ltau 
           !    it=tau(m)
           !    chitautot(m)=chitautot(m) + exp(-it*espace(isector)%e(i))*&
           !         exp(-(beta-it)*espace(isector)%e(j))*pesotot
           !    !chitautot(m)=chitautot(m) + pesotot*exp(-beta*Ej)*exp(-it*de)
           ! enddo
           ! do m=1,Lreal
           !    w0=wr(m);iw=cmplx(w0,eps,8)
           !    ! chiwtot(m)=chiwtot(m)-exp(-beta*espace(isector)%e(j))*&
           !    !      (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*pesotot
           !    chiwtot(m)=chiwtot(m)+&
           !         exp(-beta*Ej)*pesotot*(exp(-beta*de)-1.d0)/(iw-de)
           ! enddo
        enddo
     enddo
  enddo
  call stop_progress
  call print_imp_chi()
  deallocate(Chitau,Chiw,Chiiw)
  deallocate(wm,tau,wr,vm)
end subroutine full_ed_getchi
