!Solve a self-consistent problem using different 
!iteration mixing schemes.
MODULE MIXING
  USE CONSTANTS
  USE MATRIX, only: matrix_inverse
  implicit none
  private

  interface broyden_mix
     module procedure d_broyden_mix,c_broyden_mix
  end interface broyden_mix

  interface broyden_mix_
     module procedure d_broyden_mix_,c_broyden_mix_
  end interface broyden_mix_

  public :: broyden_mix
  public :: broyden_mix_

contains

  subroutine d_broyden_mix(vout,vin,alpha,M,iter,w0)
    real(8)                                 :: vout(:)
    real(8)                                 :: vin(size(vout))
    real(8),intent(in)                      :: alpha
    integer,intent(in)                      :: M
    integer                                 :: N
    real(8),optional                        :: w0
    real(8),save                            :: w0_
    integer                                 :: iter
    integer                                 :: iter_used,ipos,inext
    real(8),allocatable,dimension(:,:),save :: Df,Dv,Beta
    real(8),allocatable,dimension(:),save   :: Curv,Cm
    real(8)                                 :: norm,gamma,curvature
    integer                                 :: i,j
    N = size(vout)
    !vout=v[n+1]
    !vin=v[n]
    !Get F[n+1]=v[n+1]-v[n]
    vout = vout-vin
    !linear mixing if M=0 or first iteration
    if(M==0)then
       vout=vin+alpha*vout
       return
    endif
    !Get broyden mixing pointer
    iter_used = min(iter-1,M)
    ipos      = iter-1-((iter-2)/M)*M
    inext     = iter-((iter-1)/M)*M
    !allocate aux arrays and define the 
    !DeltaF^(n) = F^(n+1)-F^(n)/|F^(n+1)-F^(n)| 
    !DeltaV^(n) = V^(n+1)-V^(n)/|F^(n+1)-F^(n)| 
    if(iter==1)then
       vout=vin+alpha*vout
       w0_=0.01d0;if(present(w0))w0_=w0
       allocate(Df(M,N),Dv(M,N),Beta(M,M),Cm(M),Curv(N))
       Df=zero;Dv=zero;Beta=zero;Cm=zero;Curv=zero
       return
    else
       Df(ipos,:) = Vout - Df(ipos,:)
       Dv(ipos,:) = Vin  - Dv(ipos,:)
       norm = 1.d0/sqrt(dot_product(Df(ipos,:),Df(ipos,:)))
       Df(ipos,:)=Df(ipos,:)/norm
       Dv(ipos,:)=Dv(ipos,:)/norm
    endif
    !Build 
    !we are assuming w(i)=w(j)=1
    !beta  = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
    !cm(i) = Df(i)*f
    beta=zero
    do i=1,iter_used
       do j=1,iter_used
          beta(i,j) = dot_product(Df(i,:),Df(j,:))
       enddo
       beta(i,i) = beta(i,i)+w0_*w0_
       cm(i) = dot_product(Df(i,:),vout)
    enddo
    call matrix_inverse(beta(:iter_used,:iter_used))
    !v^{i+1} = v^i + alpha*f - sum_i sum_j cm(j)*b(j,i)*w(j)*u(i)*w(i)
    !u(i) = (alpha*DeltaF(i) + DeltaV(i))
    curv = alpha*vout !alpha*f
    do i=1,iter_used
       gamma=dot_product(beta(:,i),Cm)
       curv(:) = curv(:) - gamma*(Dv(i,:)+alpha*Df(i,:))
    enddo
    !Store the next entries for DeltaF and DeltaV
    Df(inext,:)=vout
    Dv(inext,:)=vin
    curvature=dot_product(vout,curv)
    if(curvature>-1.d0)then
       vout = vin + curv
    else
       vout = vin + alpha*0.5d0*vout
    endif
  end subroutine d_broyden_mix


  subroutine c_broyden_mix(vout,vin,alpha,M,iter,w0)
    complex(8)                                 :: vout(:)
    complex(8)                                 :: vin(size(vout))
    real(8),intent(in)                         :: alpha
    integer,intent(in)                         :: M
    integer                                    :: N
    real(8),optional                           :: w0
    real(8),save                               :: w0_
    integer                                    :: iter
    integer                                    :: iter_used,ipos,inext
    complex(8),allocatable,dimension(:,:),save :: Df,Dv,Beta
    complex(8),allocatable,dimension(:),save   :: Curv,Cm
    real(8)                                    :: norm,gamma,curvature
    integer                                    :: i,j
    N=size(vout)
    !vout=v[n+1]
    !vin=v[n]
    !Get F[n+1]=v[n+1]-v[n]
    vout = vout-vin
    !linear mixing if M=0 or first iteration
    if(M==0)then
       vout=vin+alpha*vout
       return
    endif
    !Get broyden mixing pointer
    iter_used = min(iter-1,M)
    ipos      = iter-1-((iter-2)/M)*M
    inext     = iter-((iter-1)/M)*M
    !allocate aux arrays and define the 
    !DeltaF^(n) = F^(n+1)-F^(n)/|F^(n+1)-F^(n)| 
    !DeltaV^(n) = V^(n+1)-V^(n)/|F^(n+1)-F^(n)| 
    if(iter==1)then
       vout=vin+alpha*vout
       w0_=0.01d0;if(present(w0))w0_=w0
       allocate(Df(M,N),Dv(M,N),Beta(M,M),Cm(M),Curv(N))
       Df=zero;Dv=zero;Beta=zero;Cm=zero;Curv=zero
       return
    else
       Df(ipos,:) = Vout - Df(ipos,:)
       Dv(ipos,:) = Vin  - Dv(ipos,:)
       norm = 1.d0/sqrt(dot_product(Df(ipos,:),Df(ipos,:)))
       Df(ipos,:)=Df(ipos,:)/norm
       Dv(ipos,:)=Dv(ipos,:)/norm
    endif
    !Build 
    !we are assuming w(i)=w(j)=1
    !beta  = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
    !cm(i) = Df(i)*f
    beta=zero
    do i=1,iter_used
       do j=1,iter_used
          beta(i,j) = dot_product(Df(i,:),Df(j,:))
       enddo
       beta(i,i) = beta(i,i)+w0_*w0_
       cm(i) = dot_product(Df(i,:),vout)
    enddo
    call matrix_inverse(beta(:iter_used,:iter_used))
    !v^{i+1} = v^i + alpha*f - sum_i sum_j cm(j)*b(j,i)*w(j)*u(i)*w(i)
    !u(i) = (alpha*DeltaF(i) + DeltaV(i))
    curv = alpha*vout !alpha*f
    do i=1,iter_used
       gamma= dot_product(beta(:,i),Cm)
       curv = curv - gamma*(Dv(i,:)+alpha*Df(i,:))
    enddo
    !Store the next entries for DeltaF and DeltaV
    Df(inext,:)=vout
    Dv(inext,:)=vin
    curvature=dot_product(vout,curv)
    if(curvature>-1.d0)then
       vout = vin + curv
    else
       vout = vin + alpha*0.5d0*vout
    endif
  end subroutine c_broyden_mix





  !*********************************************************************
  !*********************************************************************
  !*********************************************************************





  subroutine d_broyden_mix_(N,v2,v1,alpha,M,iter)
    integer                  :: N
    real(8)                  :: v1(N)
    real(8)                  :: v2(N)
    real(8)                  :: alpha
    integer                  :: M
    integer                  :: iter
    !allocatable array to be saved here
    real(8),allocatable,save :: U(:,:)
    real(8),allocatable,save :: Vt(:,:)
    real(8),allocatable,save :: F(:)
    real(8),allocatable,save :: DF(:)
    real(8),allocatable,save :: Vold(:)
    integer,save             :: last_iter
    real(8)                  :: b(M,M)
    real(8)                  :: cm(M)
    !other variables
    real(8)                  :: w0
    real(8)                  :: gmi
    real(8)                  :: wtmp,dfnorm,fnorm,rms
    integer                  :: i,j,k,iter_broyden
    !the weights w(k=1,M) are set to be 1 so all the terms
    !w(i)*w(j) disappear from equations.
    if(M==0)then
       v2 = v1 + alpha*(v2-v1)
       return
    endif
    if(iter==1)then
       ! first iteration:allocate saved arrays 
       allocate(U(N,M),Vt(N,M),F(N),DF(N),Vold(N))
       ! preform linear mixing
       last_iter = 0
       f    = v2 - v1
       vold = v1
       v2   = v1 + alpha*f
    else
       !iter > 1: this is where the non-linear mixing is done
       ! last_iter = last_iter+1      ! update pointers
       ! if( iter > M) then     ! set current lenght of broyden cycle
       !    iter_broyden = M
       ! else
       !    iter_broyden = last_iter         !lastm1
       ! endif
       iter_broyden=min(iter-1,M)
       !
       w0=0.01d0               ! set weighting factor for the zeroth iteration
       !
       !df   := f[i] - f[i-1]
       !f[i] := vector(2)[i] - vector(1)[i]
       Df = v2 - v1 - f
       f  = v2 - v1
       ! vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
       !      vold := vector(1) 
       ! vector(1) := df/|df|
       dfnorm=sqrt(dot_product(df,df))
       v2   = alpha*df/dfnorm + (v1-vold)/dfnorm
       vold = v1
       v1   = df/dfnorm
       !store vector(1) and vector(2) in the stacks u and vt restpectively
       if( iter-1 < M ) then
          U(:,iter-1)  = v2
          Vt(:,iter-1) = v1
       else
          do j = 1,M - 1
             U(:,j)  = U(:,j+1)
             Vt(:,j) = Vt(:,j+1)
          enddo
          U(:,M)  = v2
          Vt(:,M) = v1
       endif
       !calculate coefficient matrices, beta(i,j) and sum cm(i) for corrections:
       !beta(i,j) = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
       !a(i,j) = v(i)*v(j)
       !cm(i)  = sum_{l=1,N} v(l)*f(l)
       do i=1,iter_broyden
          do j = i+1,iter_broyden
             b(i,j) = dot_product(vt(:,j),vt(:,i))
             b(j,i) = b(i,j)
          enddo
          b(i,i) = dot_product(vt(:,i),vt(:,i))+w0**2
          cm(i)  = dot_product(vt(:,i),f)
       enddo
       call matrix_inverse(b(1:iter_broyden,1:iter_broyden))
       !mix vectors: v2 = vold + alpha*f - sum_i sum_j cm(j)*beta(j,i)*w(j)*u(i)*w(i)
       v2 = vold + alpha*f
       do i=1,iter_broyden
          gmi = dot_product(cm,b(:,i))
          v2 = v2 - gmi*u(:,i)
       enddo
    endif
  end subroutine d_broyden_mix_


  subroutine c_broyden_mix_(N,vnew,vprv,alpha,M,iter)
    integer                     :: N
    complex(8)                  :: vnew(N)
    complex(8)                  :: vprv(N)
    real(8)                     :: alpha
    integer                     :: M
    integer                     :: iter
    !allocatable array to be saved here
    complex(8),allocatable,save :: F(:)
    complex(8),allocatable,save :: DF(:,:),DV(:,:)
    complex(8),allocatable,save :: Vold(:)
    complex(8)                  :: b(M,M)
    complex(8)                  :: cm(M)
    complex(8)                  :: df_(N),dv_(N)
    !other variables
    real(8)                     :: w0
    complex(8)                  :: gmi
    real(8)                     :: nrm
    integer                     :: i,j,k,iter_broyden
    !the weights w(k=1,M) are set to be 1.
    ! set weighting factor for the zeroth iteration
    w0=0.01d0               
    if(M==0)then
       vnew = vprv + alpha*(vnew-vprv)
       return
    endif
    if(iter==1)then
       ! first iteration:allocate saved arrays and perform linear mixing
       allocate(F(N),DF(N,M),DV(N,M),Vold(N))
       F    = vnew - vprv
       vold = vprv
       vnew = vprv + alpha*F
    else
       !iter > 1: this is where the non-linear mixing is done
       iter_broyden = min(iter-1,M)
       !Df[i]:= F[i] - F[i-1]/|f[i]-f[i-1]|
       !Dv[i]:= V[i] - V[i-1]/|f[i]-f[i-1]|
       !F[i] := V[i] - V[i-1]       
       df_  = vnew - vprv - f
       dv_  = vprv - vold
       nrm = sqrt(dot_product(df_,df_))
       df_  = df_/nrm
       dv_  = dv_/nrm
       f    = vnew - vprv
       vold = vprv
       !store DeltaF(i), DeltaV(i)
       if(iter_broyden < M) then
          Dv(:,iter_broyden) = dv_
          Df(:,iter_broyden) = df_
       else
          do j=1,M-1
             Dv(:,j) = Dv(:,j+1)
             Df(:,j) = Df(:,j+1)
          enddo
          Dv(:,M) = dv_
          Df(:,M) = df_
       endif
       !
       !calculate coefficient matrices, beta(i,j) and cm(i) for corrections:
       !b(i,j) = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*(v(i)*v(j))]^-1
       !cm(i)  = v(i)*f
       forall(i=1:iter_broyden,j=1:iter_broyden)b(i,j) = dot_product(Df(:,i),Df(:,j))
       forall(i=1:iter_broyden)
          b(i,i) = b(i,i)+w0**2
          cm(i)  = dot_product(Df(:,i),f)
       end forall
       call matrix_inverse(b(1:iter_broyden,1:iter_broyden))
       !mix vectors: vnew = vold + alpha*F - sum_i sum_j cm(j)*beta(j,i)*w(j)*u(i)*w(i)
       vnew = vold + alpha*f
       do i=1,iter_broyden
          gmi = dot_product(cm,b(:,i))
          vnew = vnew - gmi*(alpha*Df(:,i) + Dv(:,i))
       enddo
    endif
  end subroutine c_broyden_mix_





END MODULE MIXING
