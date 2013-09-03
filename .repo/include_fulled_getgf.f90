!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_getgf()
  real(8)                      :: cdgmat,matcdg
  integer,dimension(Ntot)      :: ib
  integer                      :: i,j,k,r,ll,m,in,is,ispin,iorb
  integer                      :: idim,jdim,isector,jsector,ia
  real(8)                      :: Ei,Ej
  real(8)                      :: cc,spin1,peso1
  real(8)                      :: expterm,peso,de,w0,it,chij1
  complex(8)                   :: iw
  !----------------------------------------------
  !<i|C^+|j>=<in,is,idim|C^+|jn,js,jdim>=C^+_{ij} |
  !----------------------------------------------
  call allocate_grids
  !Initialize some functions
  impGmats   =zero
  impGreal   =zero
  if(Nspin*Norb>1)call start_timer
  do ispin=1,Nspin
     do iorb=1,Norb
        call full_ed_buildgf(iorb,ispin)
     enddo
  enddo
  if(Nspin*Norb>1)call stop_timer
  call print_imp_gf
  deallocate(wm,tau,wr)
end subroutine full_ed_getgf




subroutine full_ed_buildgf(iorb,ispin)
  integer                 :: iorb,ispin,isite
  real(8)                 :: cdgmat,matcdg
  integer,dimension(Ntot) :: ib
  integer                 :: i,j,k,r,ll,m,in,is
  integer                 :: idim,jdim,isector,jsector,ia
  real(8)                 :: Ei,Ej
  real(8)                 :: cc,spin1,peso1
  real(8)                 :: expterm,peso,de,w0,it,chij1
  complex(8)              :: iw
  isite=impIndex(iorb,ispin)
  call msg("Evaluating G_imp_Orb"//trim(txtfy(iorb))//"_Spin"//trim(txtfy(ispin)),unit=LOGfile)
  call start_timer
  do isector=startloop,lastloop
     !if(isector < minCsector(ispin))cycle
     jsector=getCsector(1,isector);if(jsector==0)cycle
     call eta(isector,lastloop,file="ETA_GF_Orb"//trim(txtfy(iorb))//"_Spin"//trim(txtfy(ispin))//".ed")
     idim=getdim(isector)     !i-th sector dimension
     jdim=getdim(jsector)     !j-th sector dimension
     do i=1,idim          !loop over the states in the i-th sect.
        do j=1,jdim       !loop over the states in the j-th sect.
           cdgmat=0.d0
           expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
           if(expterm > cutoff)then
              do ll=1,jdim              !loop over the component of |j> (IN state!)
                 m=Hmap(jsector)%map(ll)!map from IN state (j) 2 full Hilbert space
                 call bdecomp(m,ib)
                 if(ib(isite) == 0)then
                    call cdg(isite,m,k);cc=dble(k)/dble(abs(k));k=abs(k)
                    r=invHmap(isector,k)
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
                 impGmats(ispin,m)=impGmats(ispin,m)+matcdg/(iw+de)
              enddo
              !build Real-freq. GF
              do m=1,Nw 
                 w0=wr(m);iw=cmplx(w0,eps)
                 impGreal(ispin,m)=impGreal(ispin,m)+matcdg/(iw+de)
              enddo
           endif
        enddo
     enddo
  enddo
  call stop_timer
end subroutine full_ed_buildgf



!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine full_ed_getchi()
  ! USE STATISTICS
  real(8)                  :: cdgmat,matcdg
  integer,dimension(Ntot)  :: ib,ibi,ibj
  integer                  :: i,j,k,r,ll,m,iup,idw,ispin,kk
  integer                  :: idim,jdim,isector,jsector,ia,unit(2)
  real(8)                  :: cc,spin,peso,chij,weigth
  real(8)                  :: de,w0,it,Ei,Ej,expterm,Pchi(Nsect),totPchi
  complex(8)               :: iw

  real(8),allocatable,dimension(:)      :: Chitau,vm
  complex(8),allocatable,dimension(:)   :: Chiw,Chiiw

  call allocate_grids

  allocate(vm(0:NL))
  do i=0,NL
     vm(i) = pi/beta*2.d0*real(i,8)
  enddo
  allocate(Chitau(0:Ltau),Chiw(Nw),Chiiw(0:NL))

  Chitau=0.d0
  Chiw=zero
  Chiiw=zero

  !Spin susceptibility \X(tau). |<i|S_z|j>|^2
  call msg("Evaluating Chi_Sz",unit=LOGfile)
  call start_timer
  do isector=1,Nsect !loop over <i| total particle number
     call eta(isector,lastloop,file="ETA_chi.ed")
     idim=getdim(isector)
     do i=1,idim 
        do j=1,idim
           chij=0.d0
           expterm=exp(-beta*espace(isector)%e(j))
           if(expterm<cutoff)cycle
           do ll=1,idim 
              ia=Hmap(isector)%map(ll)
              call bdecomp(ia,ib)
              spin=real(ib(1),8)-real(ib(1+Ns),8) !nup - ndw
              chij=chij+espace(isector)%M(ll,i)*spin*espace(isector)%M(ll,j)
           enddo
           Ei=espace(isector)%e(i)
           Ej=espace(isector)%e(j)
           de=Ei-Ej
           peso=chij**2/zeta_function

           !Matsubara (bosonic) frequency
           if(de>cutoff)chiiw(0)=chiiw(0)-peso*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/de
           do m=1,NL
              iw=xi*vm(m)
              chiiw(m)=chiiw(m)+peso*exp(-beta*Ej)*(exp(-beta*de)-1.d0)/(iw-de)
           enddo

           !Imaginary time:
           do m=0,Ltau
              it=tau(m)
              chitau(m)=chitau(m) + peso*exp(-beta*Ej)*exp(-it*de)
           enddo

           !Real-frequency
           do m=1,Nw
              w0=wr(m);iw=cmplx(w0,eps,8)
              !Time-ordered
              ! chiw(m)=chiw(m)-exp(-beta*espace(isector)%e(j))*&
              !      (one/(w0+xi*eps+de) + one/(w0-xi*eps-de))*peso
              !Retarded = Commutator = response function
              chiw(m)=chiw(m)+&
                   exp(-beta*Ej)*peso*(exp(-beta*de)-1.d0)/(iw-de)
           enddo
        enddo
     enddo
  enddo
  call stop_timer

  unit(1)=free_unit()
  open(unit(1),file=trim(CTfile))
  unit(2)=free_unit()
  open(unit(2),file="Chi_w.ed")!trim(CWfile))
  do i=0,Ltau/2
     write(unit(1),*)tau(i),chitau(i)
  enddo
  do i=1,Nw
     if(wr(i)>=0.d0)&
          write(unit(2),*)wr(i),dimag(chiw(i)),dreal(chiw(i))
  enddo
  close(unit(1))
  close(unit(2))

  m=free_unit()
  open(m,file="Chi_iv.ed")
  do i=0,NL
     write(m,*)vm(i),dimag(chiiw(i)),dreal(chiiw(i))
  enddo
  close(m)

  deallocate(wm,tau,wr,vm)
end subroutine full_ed_getchi



