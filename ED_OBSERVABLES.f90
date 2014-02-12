!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none
  private

  public                             :: ed_getobs

  logical,save                       :: iolegend=.true.
  integer,save                       :: loop=0
  real(8),dimension(:),allocatable   :: nupimp,ndwimp,magimp
  real(8),dimension(:,:),allocatable :: sz2imp,n2imp
  real(8),dimensioN(:,:),allocatable :: zimp,simp
  real(8)                            :: s2tot
contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine ed_getobs()
    integer,dimension(Ntot)          :: ib
    integer                          :: i,j
    integer                          :: k,r
    integer                          :: ia,isector
    integer                          :: izero,isect0,jsect0,m
    integer                          :: dim,dim0,iup,idw
    integer                          :: iorb,jorb,ispin,numstates
    real(8)                          :: gs_weight
    real(8)                          :: Ei,Egs
    real(8)                          :: norm0,sgn,peso
    real(8)                          :: nup(Norb),ndw(Norb),Sz(Norb),nt(Norb)
    real(8)                          :: factor
    real(8),dimension(:),pointer     :: gsvec
    complex(8),dimension(:),pointer  :: gscvec
    integer,allocatable,dimension(:) :: Hmap
    !
    if(mpiID==0)then
       write(LOGfile,"(A)")"Evaluating Observables:"
       allocate(nupimp(Norb),ndwimp(Norb),magimp(Norb),sz2imp(Norb,Norb),n2imp(Norb,Norb))
       Egs    = state_list%emin
       nimp   = 0.d0
       nupimp = 0.d0
       ndwimp = 0.d0
       dimp   = 0.d0
       magimp = 0.d0
       sz2imp = 0.d0
       n2imp  = 0.d0
       s2tot  = 0.d0
       select case(ed_method)
       case default
          numstates=numgs
          if(finiteT)numstates=state_list%size
          do izero=1,numstates
             isect0 = es_return_sector(state_list,izero)
             Ei     = es_return_energy(state_list,izero)
             dim0   = getdim(isect0)
             if(ed_type=='d')then
                gsvec  => es_return_vector(state_list,izero)
                norm0=sqrt(dot_product(gsvec,gsvec))
             elseif(ed_type=='c')then
                gscvec  => es_return_cvector(state_list,izero)
                norm0=sqrt(dot_product(gscvec,gscvec))
             endif
             if(abs(norm0-1.d0)>1.d-9)then
                write(LOGfile,*) "GS : "//reg(txtfy(izero))//"is not normalized:"//txtfy(norm0)
                stop
             endif
             peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
             peso = peso/zeta_function
             !
             allocate(Hmap(dim0))
             call build_sector(isect0,Hmap)
             !
             do i=1,dim0
                m=Hmap(i)
                call bdecomp(m,ib)
                if(ed_type=='d')then
                   gs_weight=gsvec(i)**2
                elseif(ed_type=='c')then
                   gs_weight=abs(gscvec(i))**2
                endif
                !Get operators:
                do iorb=1,Norb
                   nup(iorb)= dble(ib(iorb))
                   ndw(iorb)= dble(ib(iorb+Ns))
                   sz(iorb) =(nup(iorb)-ndw(iorb))/2.d0
                   nt(iorb) = nup(iorb)+ndw(iorb)
                enddo
                !Evaluate averages of observables:
                do iorb=1,Norb
                   nimp(iorb)   = nimp(iorb)    +  nt(iorb)*gs_weight*peso
                   nupimp(iorb) = nupimp(iorb)  +  nup(iorb)*gs_weight*peso
                   ndwimp(iorb) = ndwimp(iorb)  +  ndw(iorb)*gs_weight*peso
                   dimp(iorb)   = dimp(iorb)    +  nup(iorb)*ndw(iorb)*gs_weight*peso
                   magimp(iorb) = magimp(iorb)  +  (nup(iorb)-ndw(iorb))*gs_weight*peso
                   sz2imp(iorb,iorb) = sz2imp(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight*peso
                   n2imp(iorb,iorb)  = n2imp(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight*peso
                   do jorb=iorb+1,Norb
                      sz2imp(iorb,jorb) = sz2imp(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight*peso
                      sz2imp(jorb,iorb) = sz2imp(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight*peso
                      n2imp(iorb,jorb)  = n2imp(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight*peso
                      n2imp(jorb,iorb)  = n2imp(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight*peso
                   enddo
                enddo
                s2tot = s2tot  + (sum(sz))**2*gs_weight*peso
             enddo
             if(associated(gsvec))nullify(gsvec)
             if(associated(gscvec))nullify(gscvec)
             deallocate(Hmap)
          enddo
       case ('full')
          do isector=1,Nsect
             dim=getdim(isector)
             allocate(Hmap(dim))
             call build_sector(isector,Hmap)
             do i=1,dim
                Ei=espace(isector)%e(i)
                peso=exp(-beta*Ei)/zeta_function
                if(peso < cutoff)cycle
                do j=1,dim
                   ia=Hmap(j)
                   gs_weight=espace(isector)%M(j,i)**2
                   call bdecomp(ia,ib)
                   !Get operators:
                   do iorb=1,Norb
                      nup(iorb)= dble(ib(iorb))
                      ndw(iorb)= dble(ib(iorb+Ns))
                      sz(iorb) =(nup(iorb)-ndw(iorb))/2.d0
                      nt(iorb) = nup(iorb)+ndw(iorb)
                   enddo
                   !Evaluate averages of observables:
                   do iorb=1,Norb
                      nimp(iorb)   = nimp(iorb)    +  nt(iorb)*gs_weight*peso
                      nupimp(iorb) = nupimp(iorb)  +  nup(iorb)*gs_weight*peso
                      ndwimp(iorb) = ndwimp(iorb)  +  ndw(iorb)*gs_weight*peso
                      dimp(iorb)   = dimp(iorb)    +  nup(iorb)*ndw(iorb)*gs_weight*peso
                      magimp(iorb) = magimp(iorb)  +  (nup(iorb)-ndw(iorb))*gs_weight*peso
                      sz2imp(iorb,iorb) = sz2imp(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight*peso
                      n2imp(iorb,iorb)  = n2imp(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight*peso
                      do jorb=iorb+1,Norb
                         sz2imp(iorb,jorb) = sz2imp(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight*peso
                         sz2imp(jorb,iorb) = sz2imp(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight*peso
                         n2imp(iorb,jorb)  = n2imp(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight*peso
                         n2imp(jorb,iorb)  = n2imp(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight*peso
                      enddo
                   enddo
                   s2tot = s2tot  + (sum(sz))**2*gs_weight*peso
                enddo
             enddo
             deallocate(Hmap)
          enddo

       end select

       allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
       call get_szr
       if(iolegend)call write_legend
       loop=loop+1
       call write_to_unit_column()
       write(LOGfile,"(A,10f18.12,f18.12)")"nimp=  ",(nimp(iorb),iorb=1,Norb),sum(nimp)
       write(LOGfile,"(A,10f18.12)")"docc=  ",(dimp(iorb),iorb=1,Norb)
       write(LOGfile,"(A,20f18.12)")"sz2 =  ",((sz2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb)
       if(Nspin==2)then
          write(LOGfile,"(A,10f18.12)")"mag=   ",(magimp(iorb),iorb=1,Norb)
       endif
       write(LOGfile,*)""
       deallocate(nupimp,ndwimp,magimp,sz2imp,n2imp)
       deallocate(simp,zimp)
    endif

#ifdef _MPI
    call MPI_BCAST(nimp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpiERR)
#endif

  end subroutine ed_getobs



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
          simp(iorb,ispin) = dimag(impSmats(ispin,ispin,iorb,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,ispin,iorb,iorb,2))-dimag(impSmats(ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,ispin,iorb,iorb,1))/wm1 ))
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="columns_info.ed")
    write(unit,"(A1,1A7,90(A10,5X))")"#","loop","xmu",&
         (reg(txtfy(2+iorb))//"nimp_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+2+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+2+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+2+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+2+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+3))//"s2",&
         ((reg(txtfy(5*Norb+3+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy(7*Norb+3+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy(9*Norb+3+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy(10*Norb+3+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend



  !+-------------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------------+
  subroutine write_to_unit_column()
    integer :: unit,flast
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all.ed",position='append')
    write(unit,"(I7,90F15.9)")loop,xmu,&
         (nimp(iorb),iorb=1,Norb),&
         (dimp(iorb),iorb=1,Norb),&
         (nupimp(iorb),iorb=1,Norb),&
         (ndwimp(iorb),iorb=1,Norb),&
         (magimp(iorb),iorb=1,Norb),&
         s2tot,&
         ((sz2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         

    unit = free_unit()
    open(unit,file="observables_last.ed")
    write(unit,"(I7,90F15.9)")loop,xmu,&
         (nimp(iorb),iorb=1,Norb),&
         (dimp(iorb),iorb=1,Norb),&
         (nupimp(iorb),iorb=1,Norb),&
         (ndwimp(iorb),iorb=1,Norb),&
         (magimp(iorb),iorb=1,Norb),&
         s2tot,&
         ((sz2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((n2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
  end subroutine write_to_unit_column


end MODULE ED_OBSERVABLES
