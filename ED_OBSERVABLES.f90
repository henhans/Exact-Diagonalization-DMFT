!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none
  private

  public              :: ed_getobs

  logical,save        :: iolegend=.true.
  integer,save        :: loop=0
  real(8),allocatable :: zimp(:,:),simp(:,:),rimp(:,:)

contains 

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine ed_getobs()
    integer,dimension(Ntot)      :: ib
    integer                      :: i,j
    integer                      :: k,r
    integer                      :: ia,isector
    integer                      :: izero,isect0,jsect0,m
    integer                      :: dim,dim0,iup,idw
    integer                      :: iorb,jorb,ispin,numstates
    real(8)                      :: gs_weight
    real(8)                      :: Ei,Egs,norm0,sgn,nup(Norb),ndw(Norb),peso
    real(8)                      :: factor
    real(8),dimension(:),pointer :: gsvec
    complex(8),dimension(:),pointer :: gscvec
    integer,allocatable,dimension(:)     :: Hmap
    !
#ifdef _MPI
    if(mpiID==0)then
#endif
       write(LOGfile,"(A)")"Evaluating Observables:"
       Egs = state_list%emin       !es_return_energy(state_list,1)
       nimp  = 0.d0
       nupimp = 0.d0
       ndwimp = 0.d0
       dimp   = 0.d0
       magimp = 0.d0
       m2imp  = 0.d0
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
             !
             peso   = 1.d0
             if(finiteT)peso=exp(-beta*(Ei-Egs))
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
                do iorb=1,Norb
                   nup(iorb)=real(ib(iorb),8)
                   ndw(iorb)=real(ib(iorb+Ns),8)
                   nimp(iorb)   = nimp(iorb)    +  (nup(iorb)+ndw(iorb))*gs_weight*peso
                   nupimp(iorb) = nupimp(iorb)  +  (nup(iorb))*gs_weight*peso
                   ndwimp(iorb) = ndwimp(iorb)  +  (ndw(iorb))*gs_weight*peso
                   dimp(iorb)   = dimp(iorb)    +  (nup(iorb)*ndw(iorb))*gs_weight*peso
                   magimp(iorb) = magimp(iorb)  +  (nup(iorb)-ndw(iorb))*gs_weight*peso
                   do jorb=1,Norb
                      m2imp(iorb,jorb)  = m2imp(iorb,jorb)  +  gs_weight*(nup(iorb)-ndw(iorb))*(nup(jorb)-ndw(jorb))*peso
                   enddo
                enddo
             enddo
             nullify(gsvec)
             deallocate(Hmap)
          enddo
          nimp  = nimp/zeta_function
          nupimp = nupimp/zeta_function
          ndwimp = ndwimp/zeta_function
          dimp   = dimp/zeta_function
          magimp = magimp/zeta_function
          m2imp  = m2imp/zeta_function
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
                   do iorb=1,Norb
                      nup(iorb)=real(ib(iorb),8)
                      ndw(iorb)=real(ib(iorb+Ns),8)
                      nimp(iorb)   = nimp(iorb)    +  (nup(iorb)+ndw(iorb))*peso*gs_weight
                      nupimp(iorb) = nupimp(iorb)  +  (nup(iorb))*peso*gs_weight
                      ndwimp(iorb) = ndwimp(iorb)  +  (ndw(iorb))*peso*gs_weight
                      dimp(iorb)   = dimp(iorb)    +  (nup(iorb)*ndw(iorb))*peso*gs_weight
                      magimp(iorb) = magimp(iorb)  +  (nup(iorb)-ndw(iorb))*peso*gs_weight
                      do jorb=1,Norb
                         m2imp(iorb,jorb)  = m2imp(iorb,jorb)  +  gs_weight*peso*(nup(iorb)-ndw(iorb))*(nup(jorb)-ndw(jorb))
                      enddo
                   enddo
                enddo
             enddo
             deallocate(Hmap)
          enddo
       end select
       allocate(simp(Norb,Nspin),zimp(Norb,Nspin))!,rimp(Norb,Nspin))
       call get_szr
       !
       if(iolegend)call write_legend
       !
       loop=loop+1
       call write_to_unit_column()
       write(LOGfile,"(A,10f18.12)")"nimp=  ",(nimp(iorb),iorb=1,Norb)
       write(LOGfile,"(A,10f18.12)")"docc=  ",(dimp(iorb),iorb=1,Norb)
       write(LOGfile,"(A,10f18.12)")"mom2=  ",(m2imp(iorb,iorb),iorb=1,Norb)
       if(Nspin==2)then
          write(LOGfile,"(A,10f18.12)")"mag=   ",(magimp(iorb),iorb=1,Norb)
       endif
       write(LOGfile,*)""
       deallocate(simp,zimp)!,rimp)

#ifdef _MPI
    endif
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
          simp(iorb,ispin) = dimag(impSmats(ispin,iorb,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,iorb,iorb,2))-dimag(impSmats(ispin,iorb,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,iorb,iorb,1))/wm1 ))
          ! rimp(iorb,ispin)   = dimag(impGmats(ispin,iorb,iorb,1)) - &
          !      wm1*(dimag(impGmats(ispin,iorb,iorb,2))-dimag(impGmats(ispin,iorb,iorb,1)))/(wm2-wm1)
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
    write(unit,"(A1,1A7,90A18)")"#","loop",&
         ("nimp_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         ("mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (("mom2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         (("z_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
                                ! (("rho_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         (("sig_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    ! unit = free_unit()
    ! open(unit,file="columns_dens.ed")
    ! write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
    !      ("nimp_"//reg(txtfy(iorb)),iorb=1,Norb),&
    !      ("nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
    !      ("ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
    !      ("mag_"//reg(txtfy(iorb)),iorb=1,Norb)
    ! close(unit)
    ! !
    ! unit = free_unit()
    ! open(unit,file="columns_docc_mom_z.ed")
    ! write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
    !      ("docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
    !      (("mom2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
    !      (("z_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    ! close(unit)
    ! !
    ! unit = free_unit()
    ! open(unit,file="columns_rho_sig.ed")
    ! write(unit,"(A1,2A18,1A7,90A18)")"#","xmu","beta","loop",&
    !      (("rho_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
    !      (("sig_"//reg(txtfy(iorb))//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    ! close(unit)
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
    write(unit,"(I7,90F18.12)")loop,&
         (nimp(iorb),iorb=1,Norb),&
         (dimp(iorb),iorb=1,Norb),&
         (nupimp(iorb),iorb=1,Norb),&
         (ndwimp(iorb),iorb=1,Norb),&
         (magimp(iorb),iorb=1,Norb),&
         ((m2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
                                ! ((rimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         

    unit = free_unit()
    open(unit,file=reg(Ofile)//"_last.ed")
    write(unit,"(I7,90F18.12)")loop,&
         (nimp(iorb),iorb=1,Norb),&
         (dimp(iorb),iorb=1,Norb),&
         (nupimp(iorb),iorb=1,Norb),&
         (ndwimp(iorb),iorb=1,Norb),&
         (magimp(iorb),iorb=1,Norb),&
         ((m2imp(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
         ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
                                ! ((rimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
         ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    close(unit)         
  end subroutine write_to_unit_column


end MODULE ED_OBSERVABLES
