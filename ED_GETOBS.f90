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
          peso=exp(-beta*Ei)/zeta
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

    freenimp=-(1.d0/beta)*log(zeta)



    open(10,file=trim(Ofile),access="append")
    select case(Nimp)
    case default
       write(10,"(A,f12.9)")"nimp = ",nimp1
       write(10,"(A,f12.9)")"n_up = ",nupimp1
       write(10,"(A,f12.9)")"n_dw = ",ndwimp1
       write(10,"(A,f12.9)")"docc = ",dimp1
       write(10,"(A,f12.9)")"mag  =",magimp1
       write(10,"(A,f12.9)")"mom2 =",m2imp1
       write(10,"(A,f12.9)")"freeE=",freenimp
       write(10,*)""

    case (2)
       write(10,"(A,2f12.9)")"nimp  = ",nimp1,nimp2
       write(10,"(A,2f12.9)")"n_up  = ",nupimp1,nupimp2
       write(10,"(A,2f12.9)")"n_dw  = ",ndwimp1,ndwimp2
       write(10,"(A,2f12.9)")"docc  = ",dimp1,dimp2
       write(10,"(A,2f12.9)")"mag   =",magimp1,magimp2
       write(10,"(A,2f12.9)")"mom2  =",m2imp1,m2imp2
       write(10,"(A,f12.9)")"mom12= ",m2imp12
       write(10,"(A,f12.9)")"freeEn=",freenimp
       write(10,*)""
    end select
    close(10)

  end subroutine imp_getobs


  !*********************************************************************
  !*********************************************************************
  !*********************************************************************







end MODULE ED_GETOBS
