!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_getgf()
  integer :: iorb,jorb,ispin
  !----------------------------------------------
  !<i|C^+|j>=<in,is,idim|C^+|jn,js,jdim>=C^+_{ij} |
  !----------------------------------------------
  call allocate_grids
  !Initialize some functions
  Giw   =zero
  Gwr   =zero
  call start_timer
  do ispin=1,Nspin
     do iorb=1,Norb
        call full_ed_buildgf(iorb,ispin)
     enddo
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           call full_ed_buildgf_mix(iorb,jorb,ispin)
        enddo
     enddo
  enddo
  call stop_timer
  call print_imp_gf
  deallocate(wm,tau,wr)
end subroutine full_ed_getgf



!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_buildgf(iorb,ispin)
  integer                 :: iorb,ispin,isite,jsite,nsite
  real(8),allocatable     :: cdgmat(:),cc(:)
  integer,dimension(Ntot) :: ib
  integer                 :: i,j,k,r,ll,m,in,is
  integer                 :: idim,jdim,isector,jsector,ia
  real(8)                 :: Ei,Ej,matcdg
  real(8)                 :: expterm,peso,de,w0,it,chij1
  complex(8)              :: iw
  nsite=1
  isite=impIndex(iorb,ispin)
  allocate(cdgmat(nsite),cc(nsite))
  call msg("Evaluating G_imp_Orb"//reg(txtfy(iorb))//reg(txtfy(iorb))//&
       "_Spin"//reg(txtfy(ispin)),unit=LOGfile)

  do isector=startloop,lastloop
     jsector=getCsector(1,isector);if(jsector==0)cycle
     !call eta(isector,lastloop,file="ETA_GF_Orb"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_Spin"//reg(txtfy(ispin))//".ed")
     idim=getdim(isector)     !i-th sector dimension
     jdim=getdim(jsector)     !j-th sector dimension
     do i=1,idim          !loop over the states in the i-th sect.
        do j=1,jdim       !loop over the states in the j-th sect.
           cdgmat=0.d0
           expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
           if(expterm < cutoff)cycle
           !
           do ll=1,jdim              !loop over the component of |j> (IN state!)
              m=Hmap(jsector)%map(ll)!map from IN state (j) 2 full Hilbert space
              call bdecomp(m,ib)
              if(ib(isite) == 0)then
                 call cdg(isite,m,k);cc(1)=dble(k)/dble(abs(k));k=abs(k)
                 r=invHmap(isector,k)
                 cdgmat(1)=cdgmat(1)+espace(isector)%M(r,i)*cc(1)*espace(jsector)%M(ll,j)
              endif
           enddo
           Ei=espace(isector)%e(i)
           Ej=espace(jsector)%e(j)
           de=Ej-Ei
           peso=expterm/zeta_function
           matcdg=peso*cdgmat(1)**2
           !build Matsubara GF
           do m=1,NL
              iw=xi*wm(m)
              Giw(iorb,iorb,ispin,m)=Giw(iorb,iorb,ispin,m)+matcdg/(iw+de)
           enddo
           !build Real-freq. GF
           do m=1,Nw 
              w0=wr(m);iw=cmplx(w0,eps)
              Gwr(iorb,iorb,ispin,m)=Gwr(iorb,iorb,ispin,m)+matcdg/(iw+de)
           enddo
        enddo
     enddo
  enddo
  !call stop_timer
end subroutine full_ed_buildgf



!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_buildgf_mix(iorb,jorb,ispin)
  integer                 :: iorb,jorb,ispin,isite,jsite,nsite
  real(8),allocatable     :: cdgmat(:),cc(:)
  integer,dimension(Ntot) :: ib
  integer                 :: i,j,k,r,ll,m,in,is
  integer                 :: idim,jdim,isector,jsector,ia
  real(8)                 :: Ei,Ej,matcdg
  real(8)                 :: expterm,peso,de,w0,it,chij1
  complex(8)              :: iw
  if(iorb==jorb)then
     call error("FULL_ED_BUILDGF_MIX: jorb == iorb. Here I get off-diagonal GF only")
  else
     nsite=2
     isite=impIndex(iorb,ispin)
     jsite=impIndex(jorb,ispin)
     if(isite==jsite)call error("FULL_ED_BUILDGF_MIX: isite == jsite. Here I get off-diagonal GF only")
  endif
  allocate(cdgmat(nsite),cc(nsite))
  call msg("Evaluating G_imp_Orb"//reg(txtfy(iorb))//reg(txtfy(jorb))//&
       "_Spin"//reg(txtfy(ispin)),unit=LOGfile)
  !call start_timer
  do isector=startloop,lastloop
     jsector=getCsector(1,isector);if(jsector==0)cycle
     !call eta(isector,lastloop,file="ETA_GF_Orb"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_Spin"//reg(txtfy(ispin))//".ed")
     idim=getdim(isector)     !i-th sector dimension
     jdim=getdim(jsector)     !j-th sector dimension
     do i=1,idim          !loop over the states in the i-th sect.
        do j=1,jdim       !loop over the states in the j-th sect.
           cdgmat=0.d0
           expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
           if(expterm < cutoff)cycle
           !
           do ll=1,jdim              !loop over the component of |j> (IN state!)
              m=Hmap(jsector)%map(ll)!map from IN state (j) 2 full Hilbert space
              call bdecomp(m,ib)
              if(ib(isite) == 0)then
                 call cdg(isite,m,k);cc(1)=dble(k)/dble(abs(k));k=abs(k)
                 r=invHmap(isector,k)
                 cdgmat(1)=cdgmat(1)+espace(isector)%M(r,i)*cc(1)*espace(jsector)%M(ll,j)
              endif
              if(ib(jsite) == 0)then
                 call cdg(jsite,m,k);cc(nsite)=dble(k)/dble(abs(k));k=abs(k)
                 r=invHmap(isector,k)
                 cdgmat(nsite)=cdgmat(nsite)+espace(isector)%M(r,i)*cc(nsite)*espace(jsector)%M(ll,j)
              endif
           enddo
           Ei=espace(isector)%e(i)
           Ej=espace(jsector)%e(j)
           de=Ej-Ei
           peso=expterm/zeta_function
           matcdg=peso*cdgmat(1)*cdgmat(nsite)
           !build Matsubara GF
           do m=1,NL
              iw=xi*wm(m)
              Giw(iorb,jorb,ispin,m)=Giw(iorb,jorb,ispin,m)+matcdg/(iw+de)
           enddo
           !build Real-freq. GF
           do m=1,Nw 
              w0=wr(m);iw=cmplx(w0,eps)
              Gwr(iorb,jorb,ispin,m)=Gwr(iorb,jorb,ispin,m)+matcdg/(iw+de)
           enddo
        enddo
     enddo
  enddo
  !call stop_timer
end subroutine full_ed_buildgf_mix







!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_getchi()
  real(8)                  :: cdgmat(Norb,Norb),chij(Norb,Norb),chitot,spin,spintot
  integer,dimension(N)     :: ib(Ntot)
  integer                  :: i,j,k,r,ll,m,in,is,ispin,iorb,jorb,isector
  integer                  :: idim,ia,unit(6)
  real(8)                  :: Ei,Ej,cc,peso(Norb,Norb),pesotot
  real(8)                  :: expterm,de,w0,it
  complex(8)               :: iw
  call allocate_grids

  Chitau=0.d0
  Chiw=zero

  !Spin susceptibility \X(tau). |<i|S_z|j>|^2
  call msg("Evaluating Chi_Sz",unit=LOGfile)
  call start_timer
  do isector=1,Nsect !loop over <i| total particle number
     call eta(isector,lastloop,file="ETA_chi.ed")
     idim=getdim(isector)
     do i=1,idim 
        do j=1,idim
           chij=0.d0;chitot=0.d0
           expterm=exp(-beta*espace(isector)%e(j))
           if(expterm<cutoff)cycle
           do ll=1,idim 
              ia=Hmap(isector)%map(ll)
              call bdecomp(ia,ib)
              spintot=0.d0
              do iorb=1,Norb
                 do jorb=1,Norb
                    spin=real(ib(iorb),8)-real(ib(jorb+Ns),8) !nup - ndw
                    chij(iorb,jorb)=chij(iorb,jorb)+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
                 enddo
                 spintot=spintot+real(ib(iorb),8)-real(ib(iorb+Ns),8)
                 chitot=chitot+espace(isector)%M(ll,i)*spintot*espace(isector)%M(ll,j)
              enddo
           enddo
           Ei=espace(isector)%e(i)
           Ej=espace(isector)%e(j)
           de=Ej-Ei
           peso=chij/zeta_function
           pesotot=chitot/zeta_function
           do iorb=1,Norb
              do jorb=1,Norb
                 !Imaginary-time
                 do m=0,Ltau 
                    it=tau(m)
                    chitau(iorb,jorb,m)=chitau(iorb,jorb,m) + exp(-it*espace(isector)%e(i))*&
                         exp(-(beta-it)*espace(isector)%e(j))*peso(iorb,jorb)
                 enddo
                 !Real-frequency
                 do m=1,Nw
                    w0=wr(m);iw=cmplx(w0,eps,8)
                    chiw(iorb,jorb,m)=chiw(iorb,jorb,m)-exp(-beta*espace(isector)%e(j))*&
                         (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso(iorb,jorb)
                 enddo
              enddo
           enddo
           !Imaginary-time
           do m=0,Ltau 
              it=tau(m)
              chitautot(m)=chitautot(m) + exp(-it*espace(isector)%e(i))*&
                   exp(-(beta-it)*espace(isector)%e(j))*pesotot
           enddo
           !Real-frequency
           do m=1,Nw
              w0=wr(m);iw=cmplx(w0,eps,8)
              chiwtot(m)=chiwtot(m)-exp(-beta*espace(isector)%e(j))*&
                   (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*pesotot
           enddo

        enddo
     enddo
  enddo
  call stop_timer

  !###########################PRINTING######################################
  do iorb=1,Norb
     do jorb=1,Norb
        unit(1)=free_unit()
        open(unit(1),file=trim(CTfile)//"_orb"//reg(txtfy(iorb))//reg(txtfy(jorb))//".ed")
        unit(2)=free_unit()
        open(unit(2),file=trim(CWfile)//"_orb"//reg(txtfy(iorb))//reg(txtfy(jorb))//".ed")
        do i=0,Ltau
           write(unit(1),*)tau(i),chitau(iorb,jorb,i)
        enddo
        do i=1,Nw
           if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiw(iorb,jorb,i)),dreal(chiw(iorb,jorb,i))
        enddo
        close(unit(1))
        close(unit(2))
     enddo
  enddo
  unit(1)=free_unit()
  open(unit(1),file=trim(CTfile)//"_tot.ed")
  unit(2)=free_unit()
  open(unit(2),file=trim(CWfile)//"_tot.ed")
  do i=0,Ltau
     write(unit(1),*)tau(i),chitautot(i)
  enddo
  do i=1,Nw
     if(wr(i)>=0.d0)write(unit(2),*)wr(i),dimag(chiwtot(i)),dreal(chiwtot(i))
  enddo
  close(unit(1))
  close(unit(2))
  deallocate(wm,tau,wr)
end subroutine full_ed_getchi
