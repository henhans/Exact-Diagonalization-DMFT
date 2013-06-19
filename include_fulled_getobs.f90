!+-------------------------------------------------------------------+
!PURPOSE  : Evaluate and print out many interesting physical qties
!+-------------------------------------------------------------------+
subroutine full_ed_getobs()
  !Configuration vector
  integer,dimension(Ntot)  :: ib
  integer                  :: i,j,ispin,ia,isector
  integer                  :: dim,iup,idw,iorb,jorb
  real(8)                  :: gs
  real(8)                  :: wm1,wm2
  real(8)                  :: Ei,nup(Norb),ndw(Norb),peso
  real(8)                  :: w
  logical                  :: last
  nimp   =0.d0
  nupimp =0.d0
  ndwimp =0.d0
  dimp   =0.d0
  magimp =0.d0
  m2imp  =0.d0
  do isector=startloop,lastloop
     dim=getdim(isector)
     do i=1,dim
        Ei=espace(isector)%e(i)
        peso=exp(-beta*Ei)/zeta_function
        if(peso < cutoff)cycle
        do j=1,dim
           ia=Hmap(isector)%map(j)
           gs=espace(isector)%M(j,i)
           call bdecomp(ia,ib)
           do iorb=1,Norb
              nup(iorb)=real(ib(iorb),8)
              ndw(iorb)=real(ib(iorb+Ns),8)
              nimp(iorb)   = nimp(iorb)    +  (nup(iorb)+ndw(iorb))*peso*gs**2
              nupimp(iorb) = nupimp(iorb)  +  (nup(iorb))*peso*gs**2
              ndwimp(iorb) = ndwimp(iorb)  +  (ndw(iorb))*peso*gs**2
              dimp(iorb)   = dimp(iorb)    +  (nup(iorb)*ndw(iorb))*peso*gs**2
              magimp(iorb) = magimp(iorb)  +  (nup(iorb)-ndw(iorb))*peso*gs**2
              do jorb=1,Norb
                 m2imp(iorb,jorb)  = m2imp(iorb,jorb)  +  gs**2*peso*(nup(iorb)-ndw(iorb))*(nup(jorb)-ndw(jorb))
              enddo
           enddo
        enddo
     enddo
  enddo

  wm1 = pi/beta ; wm2=3.d0*pi/beta
  allocate(simp(Norb,Nspin),zimp(Norb,Nspin),rimp(Norb,Nspin))
  do ispin=1,Nspin
     do iorb=1,Norb
        simp(iorb,ispin) = dimag(Siw(iorb,iorb,ispin,1)) - &
             wm1*(dimag(Siw(iorb,iorb,ispin,2))-dimag(Siw(iorb,iorb,ispin,1)))/(wm2-wm1)
        zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(Siw(iorb,iorb,ispin,1))/wm1 ))
        rimp(iorb,ispin)   = dimag(Giw(iorb,iorb,ispin,1)) - &
             wm1*(dimag(Giw(iorb,iorb,ispin,2))-dimag(Giw(iorb,iorb,ispin,1)))/(wm2-wm1)
     enddo
  enddo

  if(iolegend)call write_legend

  loop=loop+1
  call write_to_unit_column()

  call msg("Main observables:",unit=LOGfile)
  write(LOGfile,"(A,10f18.10)")"nimp=  ",(nimp(iorb),iorb=1,Norb)
  write(LOGfile,"(A,10f18.12)")"docc=  ",(dimp(iorb),iorb=1,Norb)
  write(LOGfile,"(A,10f18.12)")"mom2=  ",(m2imp(iorb,iorb),iorb=1,Norb)
  if(Nspin==2)then
     write(LOGfile,"(A,10f18.12)")"mag=   ",(magimp(iorb),iorb=1,Norb)
  endif
  write(LOGfile,*)""
  deallocate(simp,zimp,rimp)
end subroutine full_ed_getobs
