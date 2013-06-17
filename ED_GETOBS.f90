!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_GETOBS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, only:bdecomp
  implicit none
  private

  public       :: full_ed_getobs
  public       :: lanc_ed_getobs

  logical,save        :: iolegend=.true.
  integer,save        :: loop=0
  real(8),allocatable :: zimp(:,:),simp(:,:),rimp(:,:)
  real(8)             :: freenimp

contains 

  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
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

    freenimp=0.d0

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
                iup=impIndex(iorb,1)
                idw=impIndex(iorb,2)
                nup(iorb)=real(ib(iup),8)
                ndw(iorb)=real(ib(idw),8)

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

    freenimp=-(1.d0/beta)*log(zeta_function)
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





     !####################################################################
     !                    LANCZOS DIAGONALIZATION (T=0, GS only)
     !####################################################################
     !+-------------------------------------------------------------------+
     !PURPOSE  : Evaluate and print out many interesting physical qties
     !+-------------------------------------------------------------------+
     subroutine lanc_ed_getobs()
       integer,dimension(Ntot)      :: ib
       integer                      :: i,j
       integer                      :: k,r
       integer                      :: izero,isect0,jsect0,m
       integer                      :: dim0,iup,idw,iorb,jorb,ispin
       real(8)                      :: gs
       real(8)                      :: wm1,wm2
       real(8)                      :: norm0,sgn,nup(Norb),ndw(Norb)
       real(8)                      :: factor
       real(8),dimension(:),pointer :: vec

       factor=real(numzero,8)
       nimp  = 0.d0
       nupimp = 0.d0
       ndwimp = 0.d0
       dimp   = 0.d0
       magimp = 0.d0
       m2imp  = 0.d0
       do izero=1,numzero   
          !GET THE GROUNDSTATE (make some checks)
          isect0 = es_get_sector(groundstate,izero)
          dim0   = getdim(isect0)
          vec    => es_get_vector(groundstate,izero)
          do i=1,dim0
             m=Hmap(isect0)%map(i)
             call bdecomp(m,ib)
             gs=vec(i)
             do iorb=1,Norb
                iup=impIndex(iorb,1)
                idw=impIndex(iorb,2)
                nup(iorb)=real(ib(iup),8)
                ndw(iorb)=real(ib(idw),8)

                nimp(iorb)   = nimp(iorb)    +  (nup(iorb)+ndw(iorb))*gs**2
                nupimp(iorb) = nupimp(iorb)  +  (nup(iorb))*gs**2
                ndwimp(iorb) = ndwimp(iorb)  +  (ndw(iorb))*gs**2
                dimp(iorb)   = dimp(iorb)    +  (nup(iorb)*ndw(iorb))*gs**2
                magimp(iorb) = magimp(iorb)  +  (nup(iorb)-ndw(iorb))*gs**2
                do jorb=1,Norb
                   m2imp(iorb,jorb)  = m2imp(iorb,jorb)  +  gs**2*(nup(iorb)-ndw(iorb))*(nup(jorb)-ndw(jorb))
                enddo
             enddo
          enddo
       enddo
       nimp  = nimp/factor
       nupimp = nupimp/factor
       ndwimp = ndwimp/factor
       dimp   = dimp/factor
       magimp = magimp/factor
       m2imp  = m2imp/factor

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
     end subroutine lanc_ed_getobs




     subroutine write_legend()
       integer :: unit,iorb,jorb,ispin
       unit = free_unit()
       open(unit,file="columns_occupation.ed")
       write(unit,"(A1,3A18,1A7,90A18)")"#","u","xmu","beta","loop",&
            ("nimp_"//reg(txtfy(iorb)),iorb=1,Norb),&
            ("docc_"//reg(txtfy(iorb)),iorb=1,Norb)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="columns_nspin_mag.ed")
       write(unit,"(A1,3A18,1A7,90A18)")"#","u","xmu","beta","loop",&
            ("nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
            ("ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
            ("mag_"//reg(txtfy(iorb)),iorb=1,Norb)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="columns_mom.ed")
       write(unit,"(A1,3A18,1A7,90A18)")"#","u","xmu","beta","loop",&
            (("mom2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="columns_z.ed")
       write(unit,"(A1,3A18,1A7,90A18)")"#","u","xmu","beta","loop",&
            (("z_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="columns_rho_sig.ed")
       write(unit,"(A1,3A18,1A7,90A18)")"#","u","xmu","beta","loop",&
            (("rho_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
            (("sig_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       iolegend=.false.
     end subroutine write_legend




     subroutine write_to_unit_column()
       integer :: unit,flast
       integer :: iorb,jorb,ispin
       unit = free_unit()
       open(unit,file=reg(Ofile)//"_occupation.ed.all",position='append')
       flast = free_unit()
       open(flast,file=reg(Ofile)//"_occupation.ed")
       write(unit,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            (nimp(iorb),iorb=1,Norb),&
            (dimp(iorb),iorb=1,Norb)
       write(flast,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            (nimp(iorb),iorb=1,Norb),&
            (dimp(iorb),iorb=1,Norb)
       close(unit);close(flast)

       unit = free_unit()
       open(unit,file=reg(Ofile)//"_nspin_mag.ed.all",position='append')
       flast = free_unit()
       open(flast,file=reg(Ofile)//"_nspin_mag.ed")
       write(unit,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            (nupimp(iorb),iorb=1,Norb),&
            (ndwimp(iorb),iorb=1,Norb),&
            (magimp(iorb),iorb=1,Norb)
       write(flast,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            (nupimp(iorb),iorb=1,Norb),&
            (ndwimp(iorb),iorb=1,Norb),&
            (magimp(iorb),iorb=1,Norb)
       close(unit);close(flast)
       !
       unit = free_unit()
       open(unit,file=reg(Ofile)//"_mom.ed.all",position='append')
       flast = free_unit()
       open(flast,file=reg(Ofile)//"_mom.ed")
       write(unit,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            ((m2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
       write(flast,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            ((m2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
       close(unit);close(flast)
       !
       unit = free_unit()
       open(unit,file=reg(Ofile)//"_z.ed.all",position='append')
       flast = free_unit()
       open(flast,file=reg(Ofile)//"_z.ed.all")
       write(unit,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       write(flast,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit);close(flast)
       !
       unit = free_unit()
       open(unit,file=reg(Ofile)//"_rho_sig.ed.all",position='append')
       flast = free_unit()
       open(flast,file=reg(Ofile)//"_rho_sig.ed")
       write(unit,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            ((rimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       write(flast,"(3F18.9,I7,90F18.9)")u,xmu,beta,loop,&
            ((rimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit);close(flast)

     end subroutine write_to_unit_column










   end MODULE ED_GETOBS
