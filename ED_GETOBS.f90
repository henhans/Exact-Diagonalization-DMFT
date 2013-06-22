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

  logical,save :: iolegend=.true.
  integer,save :: loop=0
  real(8)      :: zupimp,zdwimp
  real(8)      :: supimp,sdwimp
  real(8)      :: rupimp,rdwimp
  real(8)      :: freenimp

contains 

  !####################################################################
  !                    FULL DIAGONALIZATION
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine full_ed_getobs()
    !Configuration vector
    integer,dimension(Ntot)     :: ib
    integer                  :: i,j,ia,isector
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

    freenimp=0.d0

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

    freenimp=-(1.d0/beta)*log(zeta_function)
    wm1 = pi/beta ; wm2=3.d0*pi/beta
    supimp   = dimag(Siw(1,1)) - wm1*(dimag(Siw(1,2))-dimag(Siw(1,1)))/(wm2-wm1)
    zupimp   = 1.d0/( 1.d0 + abs( dimag(Siw(1,1))/wm1 ))
    rupimp   = dimag(Giw(1,1)) - wm1*(dimag(Giw(1,2))-dimag(Giw(1,1)))/(wm2-wm1)
    if(Nspin==2)then
       sdwimp   = dimag(Siw(2,1)) - wm1*(dimag(Siw(2,2))-dimag(Siw(2,1)))/(wm2-wm1)
       zdwimp   = 1.d0/( 1.d0 + abs( dimag(Siw(2,1))/wm1 ))
       rdwimp   = dimag(Giw(2,1)) - wm1*(dimag(Giw(2,2))-dimag(Giw(2,1)))/(wm2-wm1)
    endif

    if(iolegend)call write_legend

    loop=loop+1
    open(10,file=trim(Ofile)//".all",access="append")
    call write_to_unit_column(10)
    close(10)

    open(20,file=trim(Ofile))
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
    integer                      :: dim0,ispin
    real(8)                      :: gs
    real(8)                      :: wm1,wm2
    real(8)                      :: norm0,sgn,nup,ndw
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
          m=Hmap(isect0)%map(i)!Hmap(isect0,i)
          call bdecomp(m,ib)
          nup=real(ib(1),8)
          ndw=real(ib(1+Ns),8)
          gs=vec(i)
          nimp  = nimp  +  (nup+ndw)*gs**2
          nupimp = nupimp +  (nup)*gs**2
          ndwimp = ndwimp +  (ndw)*gs**2
          dimp   = dimp   +  (nup*ndw)*gs**2
          magimp = magimp +  (nup-ndw)*gs**2
          m2imp  = m2imp  +  gs**2*(nup-ndw)**2
       enddo
    enddo
    nimp  = nimp/factor
    nupimp = nupimp/factor
    ndwimp = ndwimp/factor
    dimp   = dimp/factor
    magimp = magimp/factor
    m2imp  = m2imp/factor

    wm1 = pi/beta ; wm2=3.d0*pi/beta
    supimp   = dimag(Siw(1,1)) - wm1*(dimag(Siw(1,2))-dimag(Siw(1,1)))/(wm2-wm1)
    zupimp   = 1.d0/( 1.d0 + abs( dimag(Siw(1,1))/wm1 ))
    rupimp   = dimag(Giw(1,1)) - wm1*(dimag(Giw(1,2))-dimag(Giw(1,1)))/(wm2-wm1)

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

  end subroutine lanc_ed_getobs







  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  subroutine write_legend()
    if(Nspin==1)then            
       open(unit=50,file="columns_info.ed")
       write(50,"(A1,3A18,1A7,8A18)")"#","1u","2xmu","3beta","4loop","5nimp","6docc","7mag","8mom2","9z","10sig","11rho","12freeE"
       close(50)
    else
       open(unit=50,file="columns_info.ed")
       write(50,"(A1,3A18,1A7,13A18)")"#","1u","2xmu","3beta","4loop","5nimp","6n_up","7n_dw","8docc","9mag","10mom2","11z_up",&
            "12_dw","13sig_up","14sig_dw","15rho_up","16rho_dw","17freeE"
       close(50)
    endif
    iolegend=.false.
  end subroutine write_legend

  subroutine write_to_unit_column(unit)
    integer :: unit
    if(Nspin==1)then
       write(unit,"(3f18.12,I7,8f18.12)")u,xmu,beta,loop,nimp,dimp,magimp,m2imp,zupimp,supimp,rupimp,freenimp
    else
       write(unit,"(3f18.12,I7,14f18.12)")u,xmu,beta,loop,nimp,nupimp,ndwimp,dimp,magimp,m2imp,zupimp,zdwimp,supimp,sdwimp,&
            rupimp,rdwimp,freenimp
    endif
  end subroutine write_to_unit_column

  subroutine write_to_unit_list(unit)
    integer :: unit
    write(unit,"(A,f18.12)")"nimp=  ",nimp
    write(unit,"(A,f18.12)")"docc=  ",dimp
    write(unit,"(A,f18.12)")"mom2=  ",m2imp
    if(Nspin==1)then
       write(unit,"(A,f18.12)")"z=     ",zupimp
       write(unit,"(A,f18.12)")"sig=   ",supimp
       write(unit,"(A,f18.12)")"rho=   ",rupimp
    else
       write(unit,"(A,f18.12)")"n_up=  ",nupimp
       write(unit,"(A,f18.12)")"n_dw=  ",ndwimp
       write(unit,"(A,f18.12)")"mag=   ",magimp
       write(unit,"(A,f18.12)")"z_up=  ",zupimp
       write(unit,"(A,f18.12)")"z_dw=  ",zdwimp
       write(unit,"(A,f18.12)")"sig_up=",supimp
       write(unit,"(A,f18.12)")"sig_dw=",sdwimp
       write(unit,"(A,f18.12)")"rho_up=",rupimp
       write(unit,"(A,f18.12)")"rho_dw=",rdwimp
    endif
    write(unit,"(A,f18.12)")"freeE= ",freenimp
    write(unit,*)""
  end subroutine write_to_unit_list




end MODULE ED_GETOBS
