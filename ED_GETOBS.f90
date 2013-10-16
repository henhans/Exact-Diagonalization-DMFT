!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_GETOBS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none
  private

  public       :: full_ed_getobs
  public       :: lanc_ed_getobs

  logical,save        :: iolegend=.true.
  integer,save        :: loop=0
  real(8),allocatable :: zimp(:,:),simp(:,:),rimp(:,:)

contains 

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
    integer                      :: dim0,iup0,idw0,iorb,jorb,ispin
    real(8)                      :: gs
    real(8)                      :: norm0,sgn,nup(Norb),ndw(Norb)
    real(8)                      :: factor
    real(8),dimension(:),pointer :: gsvec

    !!<MPI
    !if(mpiID==0)then
    !!>MPI

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
       iup0  = getnup(isect0)
       idw0  = getndw(isect0)
       gsvec  => es_get_vector(groundstate,izero)
       norm0=sqrt(dot_product(gsvec,gsvec))
       if(abs(norm0-1.d0)>1.d-9)then
          write(*,*) "GS : "//reg(txtfy(izero))//"is not normalized:"//txtfy(norm0)
          stop
       endif
       !
       do i=1,dim0
          m=Hmap(isect0)%map(i)
          call bdecomp(m,ib)
          gs=gsvec(i)
          do iorb=1,Norb
             nup(iorb)=real(ib(iorb),8)
             ndw(iorb)=real(ib(iorb+Ns),8)
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
       nullify(gsvec)
    enddo
    nimp  = nimp/factor
    nupimp = nupimp/factor
    ndwimp = ndwimp/factor
    dimp   = dimp/factor
    magimp = magimp/factor
    m2imp  = m2imp/factor

    allocate(simp(Norb,Nspin),zimp(Norb,Nspin),rimp(Norb,Nspin))
    call get_szr

    if(iolegend)call write_legend

    loop=loop+1
    call write_to_unit_column()

    call msg("Main observables:",unit=LOGfile)
    write(LOGfile,"(A,10f18.12)")"nimp=  ",(nimp(iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12)")"docc=  ",(dimp(iorb),iorb=1,Norb)
    write(LOGfile,"(A,10f18.12)")"mom2=  ",(m2imp(iorb,iorb),iorb=1,Norb)
    if(Nspin==2)then
       write(LOGfile,"(A,10f18.12)")"mag=   ",(magimp(iorb),iorb=1,Norb)
    endif
    write(LOGfile,*)""
    deallocate(simp,zimp,rimp)
    !!<MPI
    ! endif
    ! call MPI_BCAST(nimp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    !!>MPI
  end subroutine lanc_ed_getobs




  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine full_ed_getobs()
    integer,dimension(Ntot)  :: ib
    integer                  :: i,j,ispin,ia,isector
    integer                  :: dim,iup,idw,iorb,jorb
    real(8)                  :: gs
    real(8)                  :: Ei,nup(Norb),ndw(Norb),peso
    real(8)                  :: w
    logical                  :: last

    !!<MPI
    !if(mpiID==0)then
    !!>MPI

    nimp   =0.d0
    nupimp =0.d0
    ndwimp =0.d0
    dimp   =0.d0
    magimp =0.d0
    m2imp  =0.d0
    do isector=1,Nsect
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

    allocate(simp(Norb,Nspin),zimp(Norb,Nspin),rimp(Norb,Nspin))
    call get_szr
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
    !!<MPI
    ! endif
    ! call MPI_BCAST(nimp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
    !!>MPI
  end subroutine full_ed_getobs




  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: i,j,ispin,iorb,isector
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3.d0*pi/beta
    do ispin=1,Nspin
       do iorb=1,Norb
          simp(iorb,ispin) = dimag(impSmats(ispin,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,iorb,2))-dimag(impSmats(ispin,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,iorb,1))/wm1 ))
          rimp(iorb,ispin)   = dimag(impGmats(ispin,iorb,1)) - &
               wm1*(dimag(impGmats(ispin,iorb,2))-dimag(impGmats(ispin,iorb,1)))/(wm2-wm1)
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="columns_all.ed")
    write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
         ("nimp_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (("mom2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         (("z_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         (("rho_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         (("sig_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="columns_dens_docc.ed")
    write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
         ("nimp_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("docc_"//reg(txtfy(iorb)),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="columns_nspin_mag.ed")
    write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
         ("nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("mag_"//reg(txtfy(iorb)),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="columns_mom.ed")
    write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
         (("mom2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="columns_z.ed")
    write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
         (("z_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="columns_rho_sig.ed")
    write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
         (("rho_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         (("sig_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    iolegend=.false.
  end subroutine write_legend



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine write_to_unit_column()
    integer :: unit,flast
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file=reg(Ofile)//"_all.ed",position='append')
    write(unit,"(2F18.12,I7,90F18.12)")xmu,beta,loop,&
         (nimp(iorb),iorb=1,Norb),&
         (dimp(iorb),iorb=1,Norb),&
         (nupimp(iorb),iorb=1,Norb),&
         (ndwimp(iorb),iorb=1,Norb),&
         (magimp(iorb),iorb=1,Norb),&
         ((m2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((rimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
    !Last loop observables:
    flast = free_unit()
    open(flast,file=reg(Ofile)//"_dens_docc.ed")
    write(flast,"(2F18.12,I7,90F18.12)")xmu,beta,loop,&
         (nimp(iorb),iorb=1,Norb),&
         (dimp(iorb),iorb=1,Norb)
    close(flast)
    !
    flast = free_unit()
    open(flast,file=reg(Ofile)//"_nspin_mag.ed")
    write(flast,"(2F18.12,I7,90F18.12)")xmu,beta,loop,&
         (nupimp(iorb),iorb=1,Norb),&
         (ndwimp(iorb),iorb=1,Norb),&
         (magimp(iorb),iorb=1,Norb)
    close(flast)
    !
    flast = free_unit()
    open(flast,file=reg(Ofile)//"_mom.ed")
    write(flast,"(2F18.12,I7,90F18.12)")xmu,beta,loop,&
         ((m2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
    close(flast)
    !
    flast = free_unit()
    open(flast,file=reg(Ofile)//"_z.ed")
    write(flast,"(2F18.12,I7,90F18.12)")xmu,beta,loop,&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(flast)
    !
    flast = free_unit()
    open(flast,file=reg(Ofile)//"_rho_sig.ed")
    write(flast,"(2F18.12,I7,90F18.12)")xmu,beta,loop,&
         ((rimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(flast)
  end subroutine write_to_unit_column


end MODULE ED_GETOBS
