!########################################################################
!PURPOSE  : Perform the \Chi^2 fit procedure on the Delta function
!########################################################################
MODULE ED_CHI2FIT
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none
  private

  public :: chi2_fitgf

  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb
  integer,dimension(:,:),allocatable    :: IndxOrb
  integer,dimension(:),allocatable      :: getIorb,getJorb
  real(8),dimension(:,:),allocatable    :: DeltaOrb
contains

  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf(fg,bath,ispin)
    complex(8),dimension(:,:,:)        :: fg
    real(8),dimension(:),intent(inout) :: bath
    integer                            :: ispin
    real(8),dimension((Norb+1)*Nbath)  :: a
    integer                            :: iter,stride,ifirst,ilast,i,j,corb
    real(8)                            :: chi
    !
    call msg("FIT Delta function:",unit=LOGfile)
    if(size(fg,1)/=Norb)call error("CHI2_FITGF: wrong dimension 1 in Delta_input")
    if(size(fg,2)/=Norb)call error("CHI2_FITGF: wrong dimension 2 in Delta_input")
    call check_bath_dimension(bath)
    allocate(IndxOrb(Norb,Norb),DeltaOrb(Norb,Norb),&
         getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
    corb=0
    DeltaOrb=0.d0
    do i=1,Norb
       DeltaOrb(i,i)=1.d0
       do j=i,Norb
          corb=corb+1
          IndxOrb(i,j)=corb
          getIorb(corb)=i
          getJorb(corb)=j
       enddo
    enddo
    totNorb=corb
    !
    Ldelta = size(fg,3)    
    allocate(Fdelta(Norb*(Norb+1)/2,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    do i=1,totNorb
       Fdelta(i,:) = fg(getIorb(i),getJorb(i),:)
    enddo
    Xdelta = pi/beta*real(2*arange(1,Ldelta)-1,8)
    select case(CGtype)
    case(0)
       Wdelta=(/(1.d0,i=1,Ldelta)/)
    case(1)
       Wdelta=(/(real(i,8),i=1,Ldelta)/)
    case(2)
       Wdelta=Xdelta
    end select
    !
    stride=(ispin-1)*(Norb+1)*Nbath
    ifirst=stride+1 ; ilast =stride+(Norb+1)*Nbath
    a(:) = bath(ifirst:ilast)
    call fmin_cg(a,chi2,dchi2,iter,chi,itmax=cgNitmax,ftol=cgFtol)
    bath(ifirst:ilast) = a(:)
    !
    call dump_fit_result(bath,ispin)
    write(*,"(A,ES18.9,A,I5)") 'chi^2|iter = ',chi," | ",iter
    print*," "
    deallocate(Fdelta,Xdelta,Wdelta)
    deallocate(IndxOrb,DeltaOrb,getIorb,getJorb)
  end subroutine chi2_fitgf


  !********************************************************************
  !********************************************************************
  !********************************************************************



  !FUNCTION CHI^2 to be minimized:
  function chi2(a)
    real(8),dimension(:)  ::  a
    real(8)               ::  chi2
    real(8)               ::  chig(Norb*(Norb+1)/2)
    integer               ::  i,corb
    complex(8)            ::  g0(Norb*(Norb+1)/2,Ldelta)
    chig=0.d0
    do corb=1,totNorb
       do i=1,Ldelta   !Number of freq. in common to the module
          g0(corb,i)   = gand(xdelta(i),getIorb(corb),getJorb(corb),a)
       enddo
       chig(corb) = sum(abs(Fdelta(corb,:)-g0(corb,:))**2/Wdelta(:))
    enddo
    chi2=sum(chig)
  end function chi2
  !GRADIENT OF CHI^2 used in the CG minimization:
  function dchi2(a)
    real(8),dimension(:)                                 :: a
    real(8),dimension(size(a))                           :: dchi2
    real(8),dimension(Norb*(Norb+1)/2,size(a))           :: df
    integer                                              :: i,j,corb
    complex(8),dimension(Norb*(Norb+1)/2,Ldelta)         :: g0
    complex(8),dimension(Norb*(Norb+1)/2,Ldelta,size(a)) :: dg0
    df  =0.d0
    do corb=1,totNorb
       do i=1,Ldelta
          g0(corb,i)    = gand(xdelta(i),getIorb(corb),getJorb(corb),a)
          dg0(corb,i,:) = grad_gand(xdelta(i),getIorb(corb),getJorb(corb),a)
       enddo
       do j=1,size(a)
          df(corb,j)=&
               sum(dreal(Fdelta(corb,:)-g0(corb,:))*dreal(dg0(corb,:,j))/Wdelta(:) ) + &
               sum(dimag(Fdelta(corb,:)-g0(corb,:))*dimag(dg0(corb,:,j))/Wdelta(:) )
       enddo
    enddo
    dchi2 = -2.d0*sum(df,dim=1)
  end function dchi2


  function gand(w,orb1,orb2,a) result(gg)
    real(8)                       :: w
    real(8),dimension(:)          :: a
    real(8),dimension(Nbath)      :: epsk
    real(8),dimension(Norb,Nbath) :: vpsk
    complex(8)                    :: gg
    integer                       :: i,orb1,orb2,j
    epsk=a(1:Nbath)
    do j=1,Norb
       vpsk(j,:)=a((j*Nbath+1):((j+1)*Nbath))
    enddo
    gg=zero
    do i=1,Nbath
       gg=gg + vpsk(orb1,i)*vpsk(orb2,i)/(xi*w-epsk(i))
    enddo
  end function gand



  function grad_gand(w,orb1,orb2,a) result(dgz)
    real(8)                       :: w
    real(8),dimension(:)          :: a
    real(8),dimension(Nbath)      :: epsk
    real(8),dimension(Norb,Nbath) :: vpsk
    complex(8),dimension(size(a)) :: dgz
    integer                       :: i,orb1,orb2,j
    dgz =zero
    epsk=a(1:Nbath)
    do j=1,Norb
       vpsk(j,:)=a((j*Nbath+1):((j+1)*Nbath))
    enddo
    do i=1,Nbath
       dgz(i)    = Vpsk(orb1,i)*Vpsk(orb2,i)/(xi*w-epsk(i))**2
       do j=1,Norb
          dgz(i+j*Nbath) = &
               (DeltaOrb(j,orb1)*Vpsk(orb2,i) + DeltaOrb(j,orb2)*Vpsk(orb1,i))&
               /(xi*w-epsk(i))
       enddo
    enddo
  end function grad_gand






  !+-------------------------------------------------------------+
  !PURPOSE  : 
  !+-------------------------------------------------------------+
  subroutine dump_fit_result(bath,ispin)
    real(8),dimension(:)         :: bath
    complex(8),dimension(Ldelta) :: fgand
    integer                      :: i,j,k,iorb,jorb,ispin
    real(8)                      :: w
    character(len=20)            :: suffix
    do k=1,totNorb
       iorb=getIorb(k)
       jorb=getJorb(k)
       suffix="_orb"//reg(txtfy(iorb))//reg(txtfy(jorb))//".ed"
       fgand=zero
       do i=1,Ldelta
          w=Xdelta(i)
          fgand(i) = delta_and(xi*w,bath,iorb,jorb,ispin)
       enddo
       call splot("fit_delta"//reg(suffix),&
            Xdelta,dimag(Fdelta(k,:)),dimag(fgand),dreal(Fdelta(k,:)),dreal(fgand))
    enddo
  end subroutine dump_fit_result



  !*********************************************************************
  !*********************************************************************
  !*********************************************************************


end MODULE ED_CHI2FIT
