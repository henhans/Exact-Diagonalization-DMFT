MODULE LANCZOS_HTIMEV
  implicit none
  abstract interface 
     subroutine lanc_htimesv(n,vin,vout)
       integer :: n
       real(8) :: vin(n)
       real(8) :: vout(n)
     end subroutine lanc_htimesv
  end interface
END MODULE LANCZOS_HTIMEV



MODULE LANCZOS_SIMPLE
  USE LANCZOS_HTIMEV
  implicit none
  private
  procedure(lanc_htimesv),pointer :: p_hprod
  logical :: verb=.false.
  real(8) :: threshold_=1.d-12
  real(8) :: ncheck_=10
  public :: plain_lanczos_set_htimesv
  public :: plain_lanczos_get_groundstate
  !public :: plain_lanczos_tridiag
  public :: plain_lanczos_iteration
  public :: tql2

contains

  subroutine plain_lanczos_set_htimesv(Hprod)
    procedure(lanc_htimesv)  :: Hprod
    if(associated(p_hprod))nullify(p_hprod)
    p_hprod=>Hprod
  end subroutine plain_lanczos_set_htimesv



  subroutine plain_lanczos_get_groundstate(ns,nitermax,egs,vect,Nlanc,iverbose,threshold,ncheck)
    integer                              :: ns,nitermax
    real(8),dimension(nitermax)          :: egs
    real(8),dimension(ns)                :: vect
    real(8),dimension(ns)                :: vin,vout
    integer                              :: iter,nlanc
    real(8),dimension(nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag,esave
    real(8)                              :: a_,b_,norm,diff
    integer                              :: i,j,ierr
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    if(.not.associated(p_hprod))then
       print*,"PLAIN_LANCZOS: p_hprod is not set. call plain_lanczos_set_htimesv"
       stop
    endif
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       call random_number(vect)
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb)write(*,*)"PLAIN_LANCZOS: random initial vector generated:"
    endif
    !
    !============= LANCZOS ITERATION =====================
    vin = vect
    vout= 0.d0
    alanc=0.d0
    blanc=0.d0
    nlanc=0
    lanc_loop: do iter=1,Nitermax
       if(verb)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       call plain_lanczos_iteration(iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       if(verb)print*,"alanc,blanc=",a_,b_
       nlanc=nlanc+1
       alanc(iter)  =a_
       blanc(iter+1)=b_   
       diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
       forall(i=1:Nlanc)Z(i,i)=1.d0
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call tql2(Nlanc,diag,subdiag,Z,ierr)
       if(verb)then
          print *,'---> lowest eigenvalue  <---'
          write(*,*)"E_lowest    = ",diag(1)
          open(10,file="lanc_eigenvals.dat")
          do i=1,Nlanc
             write(10,*)i,diag(i)
          enddo
          close(10)
       endif
       if(nlanc >= ncheck_)then
          esave(nlanc-(Ncheck_-1))=diag(1)
          if(nlanc >= (Ncheck_+1))then
             diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
             if(verb)write(*,*)'test deltaE = ',diff
             if(abs(diff).le.threshold_)exit lanc_loop
          endif
       endif
    enddo lanc_loop
    if(.not.verb)write(*,*)'test deltaE lanczos = ',diff
    if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
    !
    !============== END LANCZOS LOOP ======================
    !
    diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
    forall(i=1:Nlanc)Z(i,i)=1.d0
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    !Get the Eigenvalues:
    egs(1:Nlanc)=diag(1:Nlanc)
    !
    !Get the Eigenvector:
    vin =vect
    vout=0.d0
    vect=0.d0
    do iter=1,nlanc
       call plain_lanczos_iteration(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine plain_lanczos_get_groundstate



  subroutine plain_lanczos_tridiag(vin,alfa,beta,Nitermax,Nlanc,iverbose,threshold)
    real(8),dimension(:),intent(inout)          :: vin
    real(8),dimension(size(vin))                :: vout
    integer                                     :: i,nitermax,ierr
    real(8),dimension(nitermax+1),intent(inout) :: alfa,beta 
    integer                                     :: iter,nlanc
    real(8)                                     :: a_,b_
    real(8),optional                            :: threshold
    logical,optional                            :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(.not.associated(p_hprod))then
       print*,"PLAIN_LANCZOS: p_hprod is not set. call plain_lanczos_set_htimesv"
       stop
    endif
    a_=0.d0
    b_=0.d0
    nlanc=0
    lanc_loop: do iter=1,nitermax
       call plain_lanczos_iteration(iter,vin,vout,a_,b_)
       if(verb)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       if(abs(b_)<threshold_)exit lanc_loop
       nlanc=nlanc+1
       if(verb)print*,"alfa,beta=",a_,b_
       alfa(iter)=a_
       beta(iter+1)=b_
    enddo lanc_loop
  end subroutine plain_lanczos_tridiag


  ! subroutine plain_lanczos_step(vin,vout,alfa,beta,Nitermax,Nlanc,iverbose,threshold,ncheck)
  !   real(8),dimension(:),intent(inout)        :: vin,vout
  !   integer                                   :: i,nitermax,ierr
  !   real(8),dimension(nitermax),intent(inout) :: alfa,beta 
  !   real(8),optional                          :: threshold
  !   integer,optional                          :: ncheck
  !   logical,optional                          :: iverbose
  !   integer                                   :: iter,nlanc
  !   real(8)                                   :: a_,b_,diff
  !   real(8),dimension(Nitermax,Nitermax)      :: Z
  !   real(8),dimension(Nitermax)               :: diag,subdiag,esave
  !   !
  !   if(present(iverbose))verb=iverbose
  !   if(present(threshold))threshold_=threshold
  !   if(present(ncheck))ncheck_=ncheck
  !   !
  !   a_=0.d0
  !   b_=0.d0
  !   nlanc=0
  !   lanc_loop: do iter=1,nitermax
  !      call plain_lanczos_iteration(iter,vin,vout,a_,b_)
  !      if(verb)then
  !         print*,""
  !         write(*,*)"Lanczos iteration:",iter
  !      endif
  !      if(abs(b_)<threshold_)exit lanc_loop
  !      nlanc=nlanc+1
  !      if(verb)print*,"alfa=",a_
  !      if(verb)print*,"beta=",b_
  !      alfa(iter)=a_
  !      beta(iter+1)=b_   
  !      diag=0.d0 ; subdiag=0.d0 ; Z=0.d0
  !      forall(i=1:Nlanc)Z(i,i)=1.d0
  !      diag(1:Nlanc)    = alfa(1:Nlanc)
  !      subdiag(2:Nlanc) = beta(2:Nlanc)
  !      call tql2(Nlanc,diag,subdiag,Z,ierr)
  !      if(verb)then
  !         print *,'---> lowest eigenvalue  <---'
  !         write(*,*)"E_lowest    = ",diag(1)
  !         open(10,file="lanc_eigenvals.dat")
  !         do i=1,Nlanc
  !            write(10,*)i,diag(i)
  !         enddo
  !         close(10)
  !      endif
  !      if(nlanc >= ncheck_)then
  !         esave(nlanc-(Ncheck_-1))=diag(1)
  !         if(nlanc >= (Ncheck_+1))then
  !            diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
  !            if(verb)write(*,*)'test deltaE = ',diff
  !            if(abs(diff).le.threshold_)exit lanc_loop
  !         endif
  !      endif
  !   enddo lanc_loop
  !   if(.not.verb)write(*,*)'test deltaE lanczos = ',diff
  !   if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
  ! end subroutine plain_lanczos_step


  subroutine plain_lanczos_iteration(iter,vin,vout,a,b)
    real(8),dimension(:),intent(inout)         :: vin
    real(8),dimension(size(vin)),intent(inout) :: vout
    real(8),dimension(size(vin))               :: dummy,tmp
    real(8),intent(inout)              :: a,b
    integer                            :: i,iter,ns_
    real(8)                            :: norm
    ns_=size(vin)
    if(iter==1)then
       norm=sqrt(dot_product(vin,vin))
       vin=vin/norm
       b=0.d0
    end if
    call p_Hprod(ns_,vin,tmp)
    tmp=tmp-b*vout
    a = dot_product(vin,tmp)
    dummy=tmp-a*vin
    b = sqrt(dot_product(dummy,dummy))
    vout = vin
    vin = dummy/b
  end subroutine plain_lanczos_iteration


  subroutine tql2 ( n, d, e, z, ierr )
    !! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
    !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
    !    tridiagonal matrix by the QL method.  The eigenvectors of a full
    !    symmetric matrix can also be found if TRED2 has been used to reduce this
    !    full matrix to tridiagonal form.
    !  Parameters:
    !    Input, integer ( kind = 4 ) N, the order of the matrix.
    !
    !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
    !    the matrix.  On output, the eigenvalues in ascending order.  If an error
    !    exit is made, the eigenvalues are correct but unordered for indices
    !    1,2,...,IERR-1.
    !
    !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
    !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
    !    On output, E has been destroyed.
    !
    !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
    !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
    !    the tridiagonal matrix are desired, Z must contain the identity matrix.
    !    On output, Z contains the orthonormal eigenvectors of the symmetric
    !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
    !    the eigenvectors associated with the stored eigenvalues.
    !
    !    Output, integer ( kind = 4 ) IERR, error flag.
    !    0, normal return,
    !    J, if the J-th eigenvalue has not been determined after
    !    30 iterations.
    !
    implicit none
    integer ( kind = 4 ) n
    real ( kind = 8 ) c
    real ( kind = 8 ) c2
    real ( kind = 8 ) c3
    real ( kind = 8 ) d(n)
    real ( kind = 8 ) dl1
    real ( kind = 8 ) e(n)
    real ( kind = 8 ) el1
    real ( kind = 8 ) f
    real ( kind = 8 ) g
    real ( kind = 8 ) h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) l1
    integer ( kind = 4 ) l2
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real ( kind = 8 ) p
    real ( kind = 8 ) r
    real ( kind = 8 ) s
    real ( kind = 8 ) s2
    real ( kind = 8 ) tst1
    real ( kind = 8 ) tst2
    real ( kind = 8 ) z(n,n)
    ierr = 0
    if ( n == 1 ) then
       return
    end if
    do i = 2, n
       e(i-1) = e(i)
    end do
    f = 0.0D+00
    tst1 = 0.0D+00
    e(n) = 0.0D+00
    do l = 1, n
       j = 0
       h = abs ( d(l) ) + abs ( e(l) )
       tst1 = max ( tst1, h )
       !
       !  Look for a small sub-diagonal element.
       !
       do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
             exit
          end if
       end do
       if ( m == l ) then
          go to 220
       end if
130    continue
       if ( 30 <= j ) then
          ierr = l
          return
       end if
       j = j + 1
       !
       !  Form shift.
       !
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
       r = pythag ( p, 1.0D+00 )
       d(l) = e(l) / ( p + sign ( r, p ) )
       d(l1) = e(l) * ( p + sign ( r, p ) )
       dl1 = d(l1)
       h = g - d(l)
       d(l2:n) = d(l2:n) - h
       f = f + h
       !
       !  QL transformation.
       !
       p = d(m)
       c = 1.0D+00
       c2 = c
       el1 = e(l1)
       s = 0.0D+00
       mml = m - l
       do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
          !
          !  Form vector.
          !
          do k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
          end do
       end do
       p = - s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + abs ( e(l) )
       if ( tst2 > tst1 ) then
          go to 130
       end if
220    continue
       d(l) = d(l) + f
    end do
    !
    !  Order eigenvalues and eigenvectors.
    !
    do ii = 2, n
       i = ii - 1
       k = i
       p = d(i)

       do j = ii, n
          if ( d(j) < p ) then
             k = j
             p = d(j)
          end if
       end do
       if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             call r8_swap ( z(j,i), z(j,k) )
          end do
       end if
    end do
    return
  end subroutine tql2

  function pythag ( a, b )
    !! PYTHAG computes SQRT ( A * A + B * B ) carefully.
    !    The formula
    !      PYTHAG = sqrt ( A * A + B * B )
    !    is reasonably accurate, but can fail if, for example, A**2 is larger
    !    than the machine overflow.  The formula can lose most of its accuracy
    !    if the sum of the squares is very large or very small.
    !  Parameters:
    !    Input, real ( kind = 8 ) A, B, the two legs of a right triangle.
    !    Output, real ( kind = 8 ) PYTHAG, the length of the hypotenuse.
    implicit none
    real ( kind = 8 ) a
    real ( kind = 8 ) b
    real ( kind = 8 ) p
    real ( kind = 8 ) pythag
    real ( kind = 8 ) r
    real ( kind = 8 ) s
    real ( kind = 8 ) t
    real ( kind = 8 ) u
    p = max ( abs ( a ), abs ( b ) )
    if ( p /= 0.0D+00 ) then
       r = ( min ( abs ( a ), abs ( b ) ) / p )**2
       do
          t = 4.0D+00 + r
          if ( t == 4.0D+00 ) then
             exit
          end if
          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u )**2 * r
       end do
    end if
    pythag = p
    return
  end function pythag

  subroutine r8_swap ( x, y )
    !! R8_SWAP swaps two R8's.
    !  Parameters:
    !    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
    !    Y have been interchanged.
    implicit none
    real ( kind = 8 ) x
    real ( kind = 8 ) y
    real ( kind = 8 ) z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap

END MODULE LANCZOS_SIMPLE
