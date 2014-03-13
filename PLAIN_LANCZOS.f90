MODULE PLAIN_LANCZOS_HTIMESV_INTERFACE
  implicit none
  abstract interface 
     subroutine lanc_htimesv_d(nloc,n,vin,vout)
       integer :: n,nloc
       real(8) :: vin(n)
       real(8) :: vout(n)
     end subroutine lanc_htimesv_d
     !
     subroutine lanc_htimesv_c(nloc,n,vin,vout)
       integer    :: n,nloc
       complex(8) :: vin(n)
       complex(8) :: vout(n)
     end subroutine lanc_htimesv_c
  end interface
END MODULE PLAIN_LANCZOS_HTIMESV_INTERFACE


MODULE PLAIN_LANCZOS
  USE COMMON_VARS
  USE IOFILE, only: reg
  USE PLAIN_LANCZOS_HTIMESV_INTERFACE
  USE ED_VARS_GLOBAL
  USE ED_INPUT_VARS, only:ed_file_suffix
  implicit none
  private
  procedure(lanc_htimesv_d),pointer     :: dp_hprod
  procedure(lanc_htimesv_c),pointer     :: cp_hprod
  logical                               :: verb=.false.
  real(8)                               :: threshold_=1.d-12
  real(8)                               :: ncheck_=10

  public :: lanczos_plain_set_htimesv_d
  public :: lanczos_plain_set_htimesv_c
  public :: lanczos_plain_delete_htimesv
  public :: lanczos_plain_iteration_d
  public :: lanczos_plain_iteration_c
  public :: lanczos_plain_tridiag_d
  public :: lanczos_plain_tridiag_c
  public :: lanczos_plain_d
  public :: lanczos_plain_c
  !
  public :: tql2

contains


  !---------------------------------------------------------------------
  !Purpose: set the H*v function operator (use abstract interface)
  !---------------------------------------------------------------------
  subroutine lanczos_plain_set_htimesv_d(Hprod)
    procedure(lanc_htimesv_d)  :: Hprod
    if(associated(dp_hprod))nullify(dp_hprod)
    dp_hprod=>Hprod
  end subroutine lanczos_plain_set_htimesv_d
  !
  subroutine lanczos_plain_set_htimesv_c(Hprod)
    procedure(lanc_htimesv_c)  :: Hprod
    if(associated(cp_hprod))nullify(cp_hprod)
    cp_hprod=>Hprod
  end subroutine lanczos_plain_set_htimesv_c


  !---------------------------------------------------------------------
  !Purpose: delete the H*v function operator (use abstract interface)
  !---------------------------------------------------------------------
  subroutine lanczos_plain_delete_htimesv()
    if(associated(dp_hprod))nullify(dp_hprod)
    if(associated(cp_hprod))nullify(cp_hprod)
  end subroutine lanczos_plain_delete_htimesv



  !---------------------------------------------------------------------
  !Purpose: plain homebrew lanczos iteration (no orthogonalization)
  !note: the a,b variables are real, even in the complex matrix case
  !to understand why check out the Gollub-Van Loan textbook.
  !a it is easy: hermiticity->diag\in\RRR
  !b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
  !the absolute value. A lemma shows that the phase can be chosen 
  !identically zero
  !---------------------------------------------------------------------
  subroutine lanczos_plain_iteration_d(iter,vin,vout,a,b)
    real(8),dimension(:),intent(inout)         :: vin
    real(8),dimension(size(vin)),intent(inout) :: vout
    real(8),dimension(size(vin))               :: dummy,tmp
    real(8),intent(inout)                      :: a,b
    integer                                    :: i,iter,ns_,nloc
    real(8)                                    :: norm
    integer                                    :: mpiQ,mpiR
    ns_=size(vin)
    if(iter==1)then
       norm=sqrt(dot_product(vin,vin))
       if(norm==0.d0)stop "lanczos_plain_iteration: norm =0!!"
       vin=vin/norm
       b=0.d0
    end if
    nloc=1
#ifdef _MPI
    mpiQ = ns_/mpiSIZE
    mpiR = 0
    if(mpiID==(mpiSIZE-1))mpiR=mod(ns_,mpiSIZE)
    nloc=mpiQ+mpiR
#endif
    call dp_Hprod(nloc,ns_,vin,tmp)
    tmp=tmp-b*vout
    a = dot_product(vin,tmp)
    dummy=tmp-a*vin
    b = sqrt(dot_product(dummy,dummy))
    vout = vin
    vin = dummy/b
  end subroutine lanczos_plain_iteration_d
  !
  subroutine lanczos_plain_iteration_c(iter,vin,vout,a,b)
    complex(8),dimension(:),intent(inout)         :: vin
    complex(8),dimension(size(vin)),intent(inout) :: vout
    complex(8),dimension(size(vin))               :: dummy,tmp
    real(8),intent(inout)                         :: a,b
    integer                                       :: i,iter,ns_,nloc
    real(8)                                       :: norm
    integer                                       :: mpiQ,mpiR
    ns_=size(vin)
    if(iter==1)then
       norm=sqrt(dot_product(vin,vin))
       if(norm==0.d0)stop "lanczos_plain_iteration: norm =0!!"
       vin=vin/norm
       b=0.d0
    end if
    nloc = 1 
#ifdef _MPI
    mpiQ = ns_/mpiSIZE
    mpiR = 0
    if(mpiID == mpiSIZE-1)mpiR=mod(Ns_,mpiSIZE)
    nloc=mpiQ+mpiR
#endif
    call cp_Hprod(nloc,ns_,vin,tmp)
    tmp=tmp-b*vout
    a = dot_product(vin,tmp)
    dummy=tmp-a*vin
    b = sqrt(dot_product(dummy,dummy))
    vout = vin
    vin = dummy/b
  end subroutine lanczos_plain_iteration_c




  !---------------------------------------------------------------------
  !Purpose: use simple Lanczos to tri-diagonalize a matrix H (defined 
  ! in the H*v function).
  !---------------------------------------------------------------------
  subroutine lanczos_plain_tridiag_d(vin,alanc,blanc,nitermax,iverbose,threshold)
    real(8),dimension(:),intent(inout)        :: vin
    real(8),dimension(size(vin))              :: vout
    real(8),dimension(nitermax),intent(inout) :: alanc
    real(8),dimension(nitermax),intent(inout) :: blanc
    integer                                   :: i,nitermax,ierr
    integer                                   :: iter,nlanc
    real(8)                                   :: a_,b_,diff
    real(8),optional                          :: threshold
    logical,optional                          :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    a_=0.d0
    b_=0.d0
    vout=0.d0
    do iter=1,nitermax
       call lanczos_plain_iteration_d(iter,vin,vout,a_,b_)
       if(verb)print*,iter,a_,b_
       if(abs(b_)<threshold_)exit
       alanc(iter)=a_
       if(iter<nitermax)blanc(iter+1)=b_
    enddo
  end subroutine lanczos_plain_tridiag_d
  !
  subroutine lanczos_plain_tridiag_c(vin,alanc,blanc,nitermax,iverbose,threshold)
    complex(8),dimension(:),intent(inout)     :: vin
    complex(8),dimension(size(vin))           :: vout
    real(8),dimension(nitermax),intent(inout) :: alanc
    real(8),dimension(nitermax),intent(inout) :: blanc
    integer                                   :: i,nitermax,ierr
    integer                                   :: iter,nlanc
    real(8)                                   :: a_,b_,diff
    real(8),optional                          :: threshold
    logical,optional                          :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    a_=0.d0
    b_=0.d0
    vout=zero
    do iter=1,nitermax
       call lanczos_plain_iteration_c(iter,vin,vout,a_,b_)
       if(verb)print*,iter,a_,b_
       if(abs(b_)<threshold_)exit
       alanc(iter)=a_
       if(iter<nitermax)blanc(iter+1)=b_
    enddo
    if(iter==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
  end subroutine lanczos_plain_tridiag_c




  !---------------------------------------------------------------------
  !Purpose: use plain lanczos to get the groundstate energy
  !---------------------------------------------------------------------
  subroutine lanczos_plain_d(ns,nitermax,egs,vect,Nlanc,iverbose,threshold,ncheck)
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
    if(.not.associated(dp_hprod))then
       print*,"LANCZOS_PLAIN: dp_hprod is not set. call lanczos_plain_set_htimesv"
       stop
    endif
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       call random_number(vect)
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb)write(*,*)"LANCZOS_PLAIN: random initial vector generated:"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
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
       call lanczos_plain_iteration_d(iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       if(verb)print*,"alanc,blanc=",a_,b_
       nlanc=nlanc+1
       alanc(iter) = a_ ; blanc(iter+1) = b_
       diag = 0.d0 ; subdiag = 0.d0 ; Z = 0.d0
       forall(i=1:Nlanc)Z(i,i)=1.d0
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call tql2(Nlanc,diag,subdiag,Z,ierr)
       if(verb)then
          print *,'---> lowest eigenvalue  <---'
          write(*,*)"E_lowest    = ",diag(1)
          open(10,file="lanc_eigenvals"//reg(ed_file_suffix)//".ed")
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
       call lanczos_plain_iteration_d(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine lanczos_plain_d

  subroutine lanczos_plain_c(ns,nitermax,egs,vect,Nlanc,iverbose,threshold,ncheck)
    integer                              :: ns,nitermax
    real(8),dimension(nitermax)          :: egs
    complex(8),dimension(ns)             :: vect
    complex(8),dimension(ns)             :: vin,vout
    integer                              :: iter,nlanc
    real(8),dimension(nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag,esave
    real(8)                              :: a_,b_,norm,diff,ran(2)
    integer                              :: i,j,ierr
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    if(.not.associated(cp_hprod))then
       print*,"LANCZOS_PLAIN: cp_hprod is not set. call lanczos_plain_set_htimesv"
       stop
    endif
    norm=dot_product(vect,vect)
    if(norm==0.d0)then
       do i=1,ns
          call random_number(ran)
          vect(i)=cmplx(ran(1),ran(2),8)
       enddo
       vect=vect/sqrt(dot_product(vect,vect))
       if(verb)write(*,*)"LANCZOS_PLAIN: random initial vector generated:"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
    vin = vect
    vout= zero
    alanc=0.d0
    blanc=0.d0
    nlanc=0
    lanc_loop: do iter=1,Nitermax
       if(verb)then
          print*,""
          write(*,*)"Lanczos iteration:",iter
       endif
       call lanczos_plain_iteration_c(iter,vin,vout,a_,b_)
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
          open(10,file="lanc_eigenvals"//reg(ed_file_suffix)//".ed")
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
       call lanczos_plain_iteration_c(iter,vin,vout,alanc(iter),blanc(iter))
       vect = vect + vin*Z(iter,1)
    end do
    norm=sqrt(dot_product(vect,vect))
    vect=vect/norm
  end subroutine lanczos_plain_c







  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---------------------------------------------------------------------
  ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
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
  !---------------------------------------------------------------------
  subroutine tql2 ( n, d, e, z, ierr )
    integer :: n
    real(8) :: c
    real(8) :: c2
    real(8) :: c3
    real(8) :: d(n)
    real(8) :: dl1
    real(8) :: e(n)
    real(8) :: el1
    real(8) :: f
    real(8) :: g
    real(8) :: h
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
    real(8) :: p
    real(8) :: r
    real(8) :: s
    real(8) :: s2
    real(8) :: tst1
    real(8) :: tst2
    real(8) :: z(n,n)
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


  !---------------------------------------------------------------------
  ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
  !    The formula
  !    PYTHAG = sqrt ( A * A + B * B )
  !    is reasonably accurate, but can fail if, for example, A**2 is larger
  !    than the machine overflow.  The formula can lose most of its accuracy
  !    if the sum of the squares is very large or very small.
  !  Parameters:
  !    Input, real(8) :: A, B, the two legs of a right triangle.
  !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
  !---------------------------------------------------------------------
  function pythag ( a, b )
    implicit none
    real(8) :: a
    real(8) :: b
    real(8) :: p
    real(8) :: pythag
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: u
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

  !---------------------------------------------------------------------
  ! PURPOSE: swaps two R8's.
  !  Parameters:
  !    Input/output, real(8) :: X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !---------------------------------------------------------------------
  subroutine r8_swap ( x, y )
    real(8) :: x
    real(8) :: y
    real(8) :: z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap

END MODULE PLAIN_LANCZOS
