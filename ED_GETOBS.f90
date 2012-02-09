!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_GETOBS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX, only:bdecomp
  implicit none
  private

  public :: imp_getobs

contains 



  !+-------------------------------------------------------------------+
  !PROGRAM  : GETOUTNS
  !TYPE     : subroutine
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine imp_getobs
    !Configuration vector
    integer,dimension(N) :: ib
    integer :: i,j,ia,isloop
    integer :: idg
    real(8) :: gs
    real(8) :: zupimp,zdwimp
    real(8) :: supimp,sdwimp
    real(8) :: rupimp1,rdwimp1,rupimp2,rdwimp2
    real(8) :: Ei,nup1,ndw1,nup2,ndw2,peso
    real(8) :: Etot,freenimp,w
    complex(8) :: iw,alpha,greend0,selfd,zita
    complex(8),dimension(NL) :: dummy
    real(8),dimension(0:NL) :: dummyt

    nimp1   =0.d0
    nupimp1 =0.d0
    ndwimp1 =0.d0
    dimp1   =0.d0
    magimp1 =0.d0
    m2imp1  =0.d0


    nimp2   =0.d0
    nupimp2 =0.d0
    ndwimp2 =0.d0
    dimp2   =0.d0
    magimp2 =0.d0
    m2imp2  =0.d0

    m2imp12 =0.d0

    freenimp=0.d0
    Etot=zero

    do isloop=startloop,lastloop
       idg=deg(isloop)
       do i=1,idg
          Ei=espace(isloop)%e(i)
          peso=exp(-beta*Ei)/zeta_function
          if(peso < 1.d-12)cycle
          do j=1,idg
             ia=nmap(isloop,j)
             gs=espace(isloop)%M(j,i)
             call bdecomp(ia,ib)

             Etot    =Etot    +  peso*Ei

             nup1=dble(ib(1))
             ndw1=dble(ib(1+Ns))

             nimp1   = nimp1   +  (nup1+ndw1)*peso*gs**2
             nupimp1 = nupimp1 +  (nup1)*peso*gs**2
             ndwimp1 = ndwimp1 +  (ndw1)*peso*gs**2
             dimp1   = dimp1   +  (nup1*ndw1)*peso*gs**2
             magimp1 = magimp1 +  (nup1-ndw1)*peso*gs**2
             m2imp1  = m2imp1  +  gs**2*peso*(nup1-ndw1)**2

             if(Nimp==2)then
                nup2=dble(ib(2))
                ndw2=dble(ib(2+Ns))
                nimp2   = nimp2   +  (nup2+ndw2)*peso*gs**2
                nupimp2 = nupimp2 +  (nup2)*peso*gs**2
                ndwimp2 = ndwimp2 +  (ndw2)*peso*gs**2
                dimp2   = dimp2   +  (nup2*ndw2)*peso*gs**2
                magimp2 = magimp2 +  (nup2-ndw2)*peso*gs**2
                m2imp2  = m2imp2  +  gs**2*peso*(nup2-ndw2)**2

                m2imp12 = m2imp12  + gs**2*peso*(nup2-ndw2)*(nup1-ndw1)
             endif
          enddo
       enddo
    enddo

    freenimp=-(1.d0/beta)*log(zeta_function)

    

    open(10,file=trim(Ofile),access="append")
    call write_to_unit_column(10)
    close(10)


  contains


    subroutine write_to_unit_column(unit)
      integer :: unit
      logical,save :: head=.true.
      if(head)write(unit,*)
      select case(Nimp)
      case default
         !
         supimp   = dimag(Siw(1,1)) - wm(1)*(dimag(Siw(1,2))-dimag(Siw(1,1)))/(wm(2)-wm(1))
         sdwimp   = dimag(Siw(2,1)) - wm(1)*(dimag(Siw(2,2))-dimag(Siw(2,1)))/(wm(2)-wm(1))
         !
         zupimp   = 1.d0/( 1.d0 + abs( dimag(Siw(1,1))/wm(1) ))
         zdwimp   = 1.d0/( 1.d0 + abs( dimag(Siw(2,1))/wm(1) ))
         !
         rupimp1   = dimag(Giw(1,1)) - wm(1)*(dimag(Giw(1,2))-dimag(Giw(1,1)))/(wm(2)-wm(1))
         rdwimp1   = dimag(Giw(2,1)) - wm(1)*(dimag(Giw(2,2))-dimag(Giw(2,1)))/(wm(2)-wm(1))
         !
         if(Nspin==1)then
            if(head)write(unit,"(11A18)")"1nimp","2docc","3mag","4mom2","5z","6sig","7rho","8freeE","9u","10xmu","11beta"
            write(unit,"(13f18.12)")nimp1,dimp1,magimp1,m2imp1,zupimp,supimp,rupimp1,freenimp,u,xmu,beta
         else
            if(head)write(unit,"(16A18)")"1nimp","2n_up","3n_dw","4docc","5mag","6mom2","7z_up","8z_dw","9sig_up","10sig_dw","11rho_up","12rho_dw","13freeE","14u","15xmu","16beta"
            write(unit,"(18f18.12)")nimp1,nupimp1,ndwimp1,dimp1,magimp1,m2imp1,zupimp,zdwimp,supimp,sdwimp,rupimp1,rdwimp1,freenimp,u,xmu,beta
         endif

      case (2)
         rupimp1   = dimag(Giw(1,1)) - wm(1)*(dimag(Giw(1,2))-dimag(Giw(1,1)))/(wm(2)-wm(1))
         rdwimp1   = dimag(Giw(2,1)) - wm(1)*(dimag(Giw(2,2))-dimag(Giw(2,1)))/(wm(2)-wm(1))
         rupimp2   = dimag(G2iw(1,1)) - wm(1)*(dimag(G2iw(1,2))-dimag(G2iw(1,1)))/(wm(2)-wm(1))
         rdwimp2   = dimag(G2iw(2,1)) - wm(1)*(dimag(G2iw(2,2))-dimag(G2iw(2,1)))/(wm(2)-wm(1))
         if(Nspin==1)then
            if(head)write(unit,"(17A18)")"1nimp_1","2nimp_2","3docc_1","4docc_2","5mag_1","6mag_2","7mom2_1","8mom2_2","9mom12","10rho_1","11rho_2","12freeE","13u","14xmu","15beta"
            write(unit,"(18f18.12)")nimp1,nimp2,dimp1,dimp2,magimp1,magimp2,m2imp1,m2imp2,m2imp12,rupimp1,rupimp2,freenimp,u,xmu,beta
         else
            if(head)write(unit,"(26A18)")"1nimp_1","2n_1up","3n_1dw","4nimp_2","5n_2up","6n_2dw","7docc_1","8docc_2","9mag_1","10mag_2","11mom2_1","12mom2_2","13mom2_12","14rho_1up","15rho_1dw","16rho_2up","17rho_2dw","19freeE","20u","21xmu","22beta"
            write(unit,"(26f18.12)")nimp1,nupimp1,ndwimp1,nimp2,nupimp2,ndwimp2,dimp1,dimp2,magimp1,magimp2,m2imp1,m2imp2,m2imp12,rupimp1,rdwimp1,rupimp2,rdwimp2,freenimp,u,xmu,beta
         endif
      end select
      head=.false.
    end subroutine write_to_unit_column



    subroutine write_to_unit_list(unit)
      integer :: unit
      select case(Nimp)
      case default
         rupimp1   = dimag(Giw(1,1)) - wm(1)*(dimag(Giw(1,2))-dimag(Giw(1,1)))/(wm(2)-wm(1))
         rdwimp1   = dimag(Giw(2,1)) - wm(1)*(dimag(Giw(2,2))-dimag(Giw(2,1)))/(wm(2)-wm(1))
         write(unit,"(A,f18.12)")"nimp=  ",nimp1
         write(unit,"(A,f18.12)")"n_up=  ",nupimp1
         write(unit,"(A,f18.12)")"n_dw=  ",ndwimp1
         write(unit,"(A,f18.12)")"docc=  ",dimp1
         write(unit,"(A,f18.12)")"mag=   ",magimp1
         write(unit,"(A,f18.12)")"mom2=  ",m2imp1
         if(Nspin==1)then
            write(unit,"(A,f18.12)")"z=     ",zupimp
            write(unit,"(A,f18.12)")"sig=   ",supimp
            write(unit,"(A,f18.12)")"rho=   ",rupimp1
         else
            write(unit,"(A,f18.12)")"z_up=  ",zupimp
            write(unit,"(A,f18.12)")"z_dw=  ",zdwimp
            write(unit,"(A,f18.12)")"sig_up=",supimp
            write(unit,"(A,f18.12)")"sig_dw=",sdwimp
            write(unit,"(A,f18.12)")"rho_up=",rupimp1
            write(unit,"(A,f18.12)")"rho_dw=",rdwimp1
         endif
         write(unit,"(A,f18.12)")"freeE= ",freenimp
         write(unit,*)""
      case (2)
         rupimp1   = dimag(Giw(1,1)) - wm(1)*(dimag(Giw(1,2))-dimag(Giw(1,1)))/(wm(2)-wm(1))
         rdwimp1   = dimag(Giw(2,1)) - wm(1)*(dimag(Giw(2,2))-dimag(Giw(2,1)))/(wm(2)-wm(1))
         rupimp2   = dimag(G2iw(1,1)) - wm(1)*(dimag(G2iw(1,2))-dimag(G2iw(1,1)))/(wm(2)-wm(1))
         rdwimp2   = dimag(G2iw(2,1)) - wm(1)*(dimag(G2iw(2,2))-dimag(G2iw(2,1)))/(wm(2)-wm(1))
         write(unit,"(A,2f18.12)")"nimp=  ",nimp1,nimp2
         write(unit,"(A,2f18.12)")"n_up=  ",nupimp1,nupimp2
         write(unit,"(A,2f18.12)")"n_dw=  ",ndwimp1,ndwimp2
         write(unit,"(A,2f18.12)")"docc=  ",dimp1,dimp2
         write(unit,"(A,2f18.12)")"mag=   ",magimp1,magimp2
         write(unit,"(A,2f18.12)")"mom2=  ",m2imp1,m2imp2
         if(Nspin==1)then
            write(unit,"(A,f18.12)")"z=     ",zupimp
            write(unit,"(A,f18.12)")"sig=   ",supimp
            write(unit,"(A,2f18.12)")"rho=   ",rupimp1,rupimp2
         else
            write(unit,"(A,f18.12)")"z_up=  ",zupimp
            write(unit,"(A,f18.12)")"z_dw=  ",zdwimp
            write(unit,"(A,f18.12)")"sig_up=",supimp
            write(unit,"(A,f18.12)")"sig_dw=",sdwimp
            write(unit,"(A,2f18.12)")"rho_up=",rupimp1,rupimp2
            write(unit,"(A,2f18.12)")"rho_dw=",rdwimp1,rdwimp2
         endif
         write(unit,"(A,f18.12)")"mom12=  ",m2imp12
         write(unit,"(A,f18.12)")"freeEn= ",freenimp
         write(unit,*)""
      end select
    end subroutine write_to_unit_list


  end subroutine imp_getobs




  !*********************************************************************
  !*********************************************************************
  !*********************************************************************







end MODULE ED_GETOBS
