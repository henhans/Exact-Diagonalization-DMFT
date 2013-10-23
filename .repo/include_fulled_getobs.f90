!+-------------------------------------------------------------------+
!PURPOSE  : Evaluate and print out many interesting physical qties
!+-------------------------------------------------------------------+
subroutine full_ed_getobs()
  !Configuration vector
  integer,dimension(Ntot)     :: ib
  integer                  :: i,j,ia,isector,ispin
  integer                  :: dim
  real(8)                  :: gs
  real(8)                  :: wm1,wm2
  real(8)                  :: Ei,nup,ndw,peso
  real(8)                  :: w
  complex(8)               :: iw,alpha,greend0,selfd,zita
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
           ia=Hmap(isector)%map(j)!Hmap(isector,j)
           gs=espace(isector)%M(j,i)
           call bdecomp(ia,ib)
           nup=real(ib(1),8)
           ndw=real(ib(1+Ns),8)
           nimp  = nimp  +  (nup+ndw)*peso*gs**2
           nupimp = nupimp +  (nup)*peso*gs**2
           ndwimp = ndwimp +  (ndw)*peso*gs**2
           dimp   = dimp   +  (nup*ndw)*peso*gs**2
           magimp = magimp +  (nup-ndw)*peso*gs**2
           m2imp  = m2imp  +  gs**2*peso*(nup-ndw)**2
        enddo
     enddo
  enddo

  wm1 = pi/beta ; wm2=3.d0*pi/beta
  allocate(simp(Nspin),zimp(Nspin),rimp(Nspin))
  do ispin=1,Nspin
     simp(ispin) = dimag(impSmats(ispin,1)) - &
          wm1*(dimag(impSmats(ispin,2))-dimag(impSmats(ispin,1)))/(wm2-wm1)
     zimp(ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,1))/wm1 ))
     rimp(ispin)   = dimag(impGmats(ispin,1)) - wm1*(dimag(impGmats(ispin,2))-dimag(impGmats(ispin,1)))/(wm2-wm1)
  enddo

  if(iolegend)call write_legend

  loop=loop+1
  open(10,file="all_"//trim(Ofile)//".ed",access="append")
  call write_to_unit_column(10)
  close(10)
  
  open(20,file=trim(Ofile)//".ed")
  call write_to_unit_column(20)
  close(20)


  call msg("Main observables:",unit=LOGfile)
  write(LOGfile,"(A,f18.10)")"nimp=  ",nimp
  write(LOGfile,"(A,f18.12)")"docc=  ",dimp
  write(LOGfile,"(A,f18.12)")"mom2=  ",m2imp
  if(Nspin==2)then
     write(LOGfile,"(A,f18.12)")"mag=   ",magimp
  endif
  write(LOGfile,*)""
  deallocate(rimp,zimp,simp)
end subroutine full_ed_getobs
