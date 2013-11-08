!+------------------------------------------------------------------+
!PURPOSE  : Evalute green's function and suscpetibility
! using complete spectral decomposition 
!+------------------------------------------------------------------+
subroutine full_ed_getgf()
  integer :: iorb,jorb,ispin
  write(LOGfile,"(A)")"Evaluating Green's functions"
  call allocate_grids
  impGmats   =zero
  impGreal   =zero
  if(mpiID==0)call start_timer
  do ispin=1,Nspin
     do iorb=1,Norb
        call full_ed_buildgf(iorb,ispin)
     enddo
  enddo
  if(mpiID==0)call print_imp_gf
  if(mpiID==0)call stop_timer
  deallocate(wm,tau,wr,vm)
end subroutine full_ed_getgf



!+------------------------------------------------------------------+
!PURPOSE  :  Evalue the spectral sum for impurity Green's function. 
!<i|C^+|j>=<in,is,idim|C^+|jn,js,jdim>=C^+_{ij} |
!+------------------------------------------------------------------+
subroutine full_ed_buildgf(iorb,ispin)
  integer                 :: iorb,ispin,isite,jsite,nsite
  real(8)                 :: cdgmat,cc
  integer,dimension(Ntot) :: ib
  integer                 :: i,j,k,r,ll,m,in,is
  integer                 :: idim,jdim,isector,jsector,ia
  real(8)                 :: Ei,Ej,matcdg
  real(8)                 :: expterm,peso,de,w0,it,chij1
  complex(8)              :: iw
  integer,allocatable,dimension(:)     :: HJmap,HImap    !map of the Sector S to Hilbert space H


  nsite=1
  isite=impIndex(iorb,ispin)
  if(mpiID==0)write(LOGfile,"(A)")"Evaluating G_imp_Orb"//reg(txtfy(iorb))//"_Spin"//reg(txtfy(ispin))
  if(mpiID==0)call start_progress(LOGfile)

  do isector=1,Nsect
     jsector=getCsector(1,isector);if(jsector==0)cycle
     if(mpiID==0)call progress(isector,Nsect)
     idim=getdim(isector)     !i-th sector dimension
     jdim=getdim(jsector)     !j-th sector dimension
     allocate(HImap(idim))
     allocate(HJmap(jdim))
     call build_sector(isector,HImap)
     call build_sector(jsector,HJmap)
     do i=1,idim          !loop over the states in the i-th sect.
        do j=1,jdim       !loop over the states in the j-th sect.
           cdgmat=0.d0
           expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
           if(expterm < cutoff)cycle
           !
           do ll=1,jdim              !loop over the component of |j> (IN state!)
              !m=Hmap(jsector)%map(ll)!map from IN state (j) to full Hilbert space
              m = HJmap(ll)
              call bdecomp(m,ib)
              if(ib(isite) == 0)then
                 call cdg(isite,m,k,cc)!;cc=dble(k)/dble(abs(k));k=abs(k)
                 r = binary_search(HImap,k)
                 cdgmat=cdgmat+espace(isector)%M(r,i)*cc*espace(jsector)%M(ll,j)
              endif
           enddo
           Ei=espace(isector)%e(i)
           Ej=espace(jsector)%e(j)
           de=Ej-Ei
           peso=expterm/zeta_function
           matcdg=peso*cdgmat**2
           !build Matsubara GF
           do m=1,NL
              iw=xi*wm(m)
              impGmats(ispin,iorb,iorb,m)=impGmats(ispin,iorb,iorb,m)+matcdg/(iw+de)
           enddo
           !build Real-freq. GF
           do m=1,Nw 
              w0=wr(m);iw=cmplx(w0,eps)
              impGreal(ispin,iorb,iorb,m)=impGreal(ispin,iorb,iorb,m)+matcdg/(iw+de)
           enddo
        enddo
     enddo
     deallocate(HImap,HJmap)
  enddo
  if(mpiID==0)call stop_progress
end subroutine full_ed_buildgf


