!!Some notes here:
!!this moule provides an interface to sparse matrix in ll format
!!the routines perform some generic action on the sparse_matrix object
!!but we still miss some elementary action such as +update_value, +check_value_exist, +delete_value
!!that I am not gonna use in ED code (for which this module is developed).
MODULE MATRIX_SPARSE
  USE COMMON_VARS
  !!<MPI
  !USE MPI
  !!>MPI
  implicit none
  private


  type sparse_element
     private
     real(8)                               :: val  !value of the entry: double precision
     complex(8)                            :: cval !value of the entry: double complex
     integer                               :: col  !col connected to this compress value
     type(sparse_element),pointer          :: next !link to next entry in the row
  end type sparse_element

  type sparse_row
     private
     integer                               :: size    !size of the list
     type(sparse_element),pointer          :: root    !head/root of the list\== list itself
  end type sparse_row


  type sparse_matrix
     integer                               :: size
     logical                               :: status=.false.
     type(sparse_row),dimension(:),pointer :: row
  end type sparse_matrix

  !Here we interface some routines that behave differently for real and compelx numbers
  interface sp_insert_element
     module procedure sp_insert_element_d,sp_insert_element_c
  end interface sp_insert_element

  interface sp_load_matrix
     module procedure sp_load_matrix_d,sp_load_matrix_c
  end interface sp_load_matrix

  interface sp_dump_matrix
     module procedure sp_dump_matrix_d,sp_dump_matrix_c
  end interface sp_dump_matrix

  interface sp_matrix_vector_product
     module procedure sp_matrix_vector_product_dd, sp_matrix_vector_product_dc
  end interface sp_matrix_vector_product

  !!<MPI
  ! interface sp_matrix_vector_product_mpi
  !    module procedure sp_matrix_vector_product_d_mpi, sp_matrix_vector_product_c_mpi
  ! end interface sp_matrix_vector_product_mpi
  !!>MPI
  
  public :: sparse_matrix
  !
  public :: sp_init_matrix      !checked
  public :: sp_delete_matrix    !checked
  !
  public :: sp_insert_element   !checked
  public :: sp_get_element_d    !checked
  public :: sp_get_element_c    !checked
  public :: sp_delete_element   !checked
  public :: sp_inquire_element  !checked
  !
  public :: sp_load_matrix      !checked
  public :: sp_dump_matrix      !checked
  public :: sp_print_matrix     !checked
  !
  public :: sp_matrix_vector_product !checked
  public :: sp_matrix_vector_product_cc !checked

  !!<MPI
  ! public :: sp_matrix_vector_product_mpi
  !!>MPI

contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix(sparse,N)
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i,N
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix: alreay allocate can not init"
    sparse%size=N
    sparse%status=.true.
    allocate(sparse%row(N))
    do i=1,N
       allocate(sparse%row(i)%root)
       sparse%row(i)%root%next => null()
       sparse%row(i)%size=0
    end do
  end subroutine sp_init_matrix




  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_d(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    real(8),intent(in)                :: value
    integer,intent(in)                :: i,j
    call insert_element_in_row_d(sparse%row(i),value,j)
  end subroutine sp_insert_element_d
  subroutine sp_insert_element_c(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),intent(in)             :: value
    integer,intent(in)                :: i,j
    call insert_element_in_row_c(sparse%row(i),value,j)
  end subroutine sp_insert_element_c


  !+------------------------------------------------------------------+
  !PURPOSE: insert an element in a given row (private) 
  !+------------------------------------------------------------------+
  subroutine insert_element_in_row_d(row,value,column)
    type(sparse_row),intent(inout)    :: row
    real(8) ,intent(in)               :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: p,c
    integer                           :: i
    p => row%root
    c => p%next
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    p%next%val = value
    p%next%col = column
    row%size   = row%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine insert_element_in_row_d

  subroutine insert_element_in_row_c(row,value,column)
    type(sparse_row),intent(inout)    :: row
    complex(8) ,intent(in)            :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: p,c
    integer                           :: i
    p => row%root
    c => p%next
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    allocate(p%next)                !Create a new element in the list
    p%next%cval= value
    p%next%col = column
    row%size   = row%size+1
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine insert_element_in_row_c






  !+------------------------------------------------------------------+
  !PURPOSE: get an element from position (i,j) of the sparse matrix
  !+------------------------------------------------------------------+
  function sp_get_element_d(sparse,i,j) result(value)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    real(8)                           :: value
    call get_element_from_row_d(sparse%row(i),value,j)
  end function sp_get_element_d

  function sp_get_element_c(sparse,i,j) result(value)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    complex(8)                        :: value
    call get_element_from_row_c(sparse%row(i),value,j)
  end function sp_get_element_c


  !+------------------------------------------------------------------+
  !PURPOSE: get an element from a given row of the matrix (private)
  !+------------------------------------------------------------------+
  subroutine get_element_from_row_d(row,value,column)
    type(sparse_row),intent(inout)    :: row
    real(8)                           :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    value=0.d0
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    value = c%val
  end subroutine get_element_from_row_d

  subroutine get_element_from_row_c(row,value,column)
    type(sparse_row),intent(inout)    :: row
    complex(8)                        :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    value=cmplx(0.d0,0.d0,8)
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    !
    value = c%cval
  end subroutine get_element_from_row_c




  !+------------------------------------------------------------------+
  !PURPOSE: check if a given element exists
  !+------------------------------------------------------------------+
  function sp_inquire_element(sparse,i,j) result(exist)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    logical                           :: exist
    exist = inquire_element_from_row(sparse%row(i),j)
  end function sp_inquire_element

  !+------------------------------------------------------------------+
  !PURPOSE: check if an element in a given row of the matrix exist (private)
  !+------------------------------------------------------------------+
  function inquire_element_from_row(row,column) result(exist)
    type(sparse_row),intent(inout)    :: row
    logical                           :: exist
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    exist=.false.
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    exist=.true.
  end function inquire_element_from_row




  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix(sparse)    
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i,Ndim
    do i=1,sparse%size
       call delete_row(sparse%row(i))
       deallocate(sparse%row(i)%root)
    end do
    deallocate(sparse%row)
    sparse%size=0
    sparse%status=.false.
  end subroutine sp_delete_matrix



  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire row from the sparse matrix (private)
  !+------------------------------------------------------------------+
  subroutine delete_row(row)
    type(sparse_row),intent(inout) :: row
    type(sparse_element),pointer   :: p,c
    do
       p => row%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       deallocate(c)
    end do
  end subroutine delete_row







  !+------------------------------------------------------------------+
  !PURPOSE: delete a single element at (i,j) from the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_element(matrix,i,j)
    type(sparse_matrix),intent(inout) :: matrix
    integer,intent(in)                :: i,j
    logical :: delete
    delete = delete_element_from_row(matrix%row(i),col=j)
    if(.not.delete)write(*,"(A,I3,I3)")"sp_delete_element: can not delete element in",i,j
  end subroutine sp_delete_element


  !+------------------------------------------------------------------+
  !PURPOSE: delete a given element from a row of the sparse matrix (private)
  !+------------------------------------------------------------------+
  function delete_element_from_row(row,n,col) result(found)
    type(sparse_row),intent(inout)    :: row
    integer,optional                  :: n
    integer,optional                  :: col
    integer                           :: i,pos
    type(sparse_element),pointer      :: p,c
    logical                           :: found
    pos= row%size ; if(present(n))pos=n
    p => row%root
    c => p%next
    found = .false.
    if(present(col))then
       do 
          if(found .OR. .not.associated(c))return
          if(col == c%col)then
             found=.true.
             exit
          else
             p => c
             c => c%next
          endif
       end do
       if(found)then
          p%next => c%next !reallocate skipping the deleted link
          deallocate(c)           !free link
          row%size=row%size-1
       endif
    else
       do i=1,pos 
          if(.not.associated(c))return !empty list
          p => c
          c => c%next
       end do
       found=.true.
       p%next => c%next !reallocate skipping the deleted link
       deallocate(c)           !free link
       row%size=row%size-1
    endif
  end function delete_element_from_row









  !+------------------------------------------------------------------+
  !PURPOSE: load a regular matrix (2dim array) into a sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_d(matrix,sparse)
    real(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_matrix),intent(inout)  :: sparse    
    integer                            :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%size /= Ndim1)stop"Warning SPARSE/load_matrix: dimensions error"
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call sp_insert_element_d(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_d

  subroutine sp_load_matrix_c(matrix,sparse)
    complex(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_matrix),intent(inout)     :: sparse    
    integer                               :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%size /= Ndim1)stop "Warning SPARSE/load_matrix: dimensions error"
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=cmplx(0.d0,0.d0,8))call sp_insert_element_c(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_c





  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_d(sparse,matrix)
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(:,:),intent(inout)  :: matrix
    type(sparse_element),pointer          :: c
    integer                               :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%size /= Ndim1)stop "Warning SPARSE/load_matrix: dimensions error"
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do 
          if(.not.associated(c))exit
          matrix(i,c%col) = c%val
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_d

  subroutine sp_dump_matrix_c(sparse,matrix)
    type(sparse_matrix),intent(in)          :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element),pointer            :: c
    integer                                 :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%size /= Ndim1)stop "Warning SPARSE/load_matrix: dimensions error"
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do 
          if(.not.associated(c))exit
          matrix(i,c%col) = c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_c





  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+
  subroutine sp_print_matrix(sparse,unit,fmt,type,full)
    type(sparse_matrix)            :: sparse
    integer,optional               :: unit
    integer                        :: i,j,unit_,Ns
    character(len=*),optional      :: fmt
    character(len=64)              :: fmt_
    character(len=1),optional      :: type
    character(len=1)               :: type_
    logical,optional               :: full
    logical                        :: full_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F8.3';if(present(fmt))fmt_=fmt
    type_='d';if(present(type))type_=type
    full_=.false.;if(present(full))full_=full
    Ns=sparse%size
    select case(type_)
    case('d')
       if(full_)then
          write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
          do i=1,Ns
             write(*,"(100"//trim(fmt_)//",1X)")(sp_get_element_d(sparse,i,j),j=1,Ns)
          enddo
       else
          write(*,*)"Print sparse matrix (compact mode) ->",unit_
          do i=1,Ns
             call print_row_d(sparse%row(i),unit_,fmt_)
          enddo
       endif
    case('c')
       if(full_)then
          write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
          do i=1,Ns
             write(*,"(100("//trim(fmt_)//",A1,"//trim(fmt_)//",2X))")(&
                  real(sp_get_element_c(sparse,i,j)),",",imag(sp_get_element_c(sparse,i,j)),j=1,Ns)
          enddo
       else
          write(*,*)"Print sparse matrix (compact mode) ->",unit_
          do i=1,Ns
             call print_row_c(sparse%row(i),unit_,fmt_)
          enddo
       endif
    end select
    write(unit_,*)
  end subroutine sp_print_matrix

  !+------------------------------------------------------------------+
  !PURPOSE: print an entire row of the sparse matrix (private)
  !+------------------------------------------------------------------+
  subroutine print_row_d(row,unit,fmt)
    type(sparse_row),intent(in)  :: row
    type(sparse_element),pointer :: c
    integer                      :: count=0,i
    integer                      :: unit
    character(len=*)             :: fmt
    c => row%root%next   !assume is associated,ie list exists
    do
       if(.not.associated(c))exit
       count=count+1       
       write(unit,"("//trim(fmt)//",A1,I3,1X)",advance='no')c%val,',',c%col
       c => c%next  !traverse list
    end do
    write(unit,*)
  end subroutine print_row_d

  subroutine print_row_c(row,unit,fmt)
    type(sparse_row),intent(in)   :: row
    type(sparse_element),pointer  :: c
    integer                       :: count=0
    integer,optional :: unit
    integer          :: unit_
    character(len=*),optional :: fmt
    character(len=64)         :: fmt_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F15.9';if(present(fmt))fmt_=fmt
    c => row%root%next   !assume is associated,ie list exists
    do
       if(.not.associated(c))exit
       count=count+1
       write(unit_,"(2"//trim(fmt_)//",A1,I3,3X)",advance='no')c%cval,',',c%col
       c => c%next  !traverse list
    end do
    write(unit_,*)
  end subroutine print_row_c






  !+------------------------------------------------------------------+
  !PURPOSE: given a vector vin, perform the matrix-vector multiplication
  ! H_sparse * vin and put the result in vout.
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_dd(Ndim,sparse,vin,vout)
    integer                               :: Ndim
    type(sparse_matrix),intent(in)        :: sparse
    real(8),dimension(Ndim),intent(in)    :: vin
    real(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer          :: c
    integer                               :: i
    vout=0.d0
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_dd
  !+------------------------------------------------------------------+
  subroutine sp_matrix_vector_product_dc(Ndim,sparse,vin,vout)
    integer                                  :: Ndim
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    ! if(.not.sparse%status)then
    !    print*,"Error SPARSE/matrix_vector_product: sparse matrix not allocated."
    !    stop
    ! endif
    ! if(Ndim/=sparse%size)then
    !    print*,"Error SPARSE/matrix_vector_product: wrong dimensions matrix vector."
    !    stop
    ! endif
    vout=cmplx(0.d0,0.d0,8)
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do  
          if(.not.associated(c))exit matmul
          ! if(c%col > ndim)then
          !    print*,"Error SPARSE_MATRIX/sp_matrix_vector_product: column index > size(vector)!"
          !    stop
          ! endif
          vout(i) = vout(i) + c%val*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_dc

  subroutine sp_matrix_vector_product_cc(Ndim,sparse,vin,vout)
    integer                                  :: Ndim
    type(sparse_matrix),intent(in)           :: sparse
    complex(8),dimension(Ndim),intent(in)    :: vin
    complex(8),dimension(Ndim),intent(inout) :: vout
    type(sparse_element),pointer             :: c
    integer                                  :: i
    vout=cmplx(0.d0,0.d0,8)
    do i=1,Ndim
       c => sparse%row(i)%root%next       
       matmul: do  
          if(.not.associated(c))exit matmul
          vout(i) = vout(i) + c%cval*vin(c%col)
          c => c%next  !traverse list
       end do matmul
    end do
  end subroutine sp_matrix_vector_product_cc


  !!<MPI
  ! subroutine sp_matrix_vector_product_d_mpi(Q,R,Ndim,sparse,vin,vout)
  !   integer                               :: Ndim
  !   type(sparse_matrix),intent(in)        :: sparse
  !   real(8),dimension(Ndim),intent(in)    :: vin
  !   real(8),dimension(Ndim),intent(inout) :: vout
  !   real(8),dimension(Ndim)               :: vtmp
  !   type(sparse_element),pointer          :: c
  !   integer                               :: i
  !   integer                               :: R,Q,Nini,Nfin
  !   vtmp=0.d0
  !   vout=0.d0
  !   ! Nini=mpiID*R + 1            !1,R+1,2R+1,...
  !   ! Nfin=(mpiID+1)*R            !R,2R,3R
  !   ! if(Q/=0 .AND. mpiID==mpiSIZE-1)Nfin=Nfin+Q
  !   do i=mpiID*Q+1,(mpiID+1)*Q+R!Nini,Nfin
  !      c => sparse%row(i)%root%next       
  !      matmul: do
  !         if(.not.associated(c))exit matmul
  !         vtmp(i) = vtmp(i) + c%val*vin(c%col)
  !         c => c%next
  !      end do matmul
  !   end do
  !   call MPI_ALLREDUCE(vtmp,vout,Ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpiERR)
  ! end subroutine sp_matrix_vector_product_d_mpi
  ! !+------------------------------------------------------------------+
  ! subroutine sp_matrix_vector_product_c_mpi(Q,R,Ndim,sparse,vin,vout)
  !   integer                                  :: Ndim
  !   type(sparse_matrix),intent(in)           :: sparse
  !   complex(8),dimension(Ndim),intent(in)    :: vin
  !   complex(8),dimension(Ndim),intent(inout) :: vout
  !   complex(8),dimension(Ndim)               :: vtmp
  !   type(sparse_element),pointer             :: c
  !   integer                                  :: i
  !   integer                                  :: mpiID,mpiSIZE,mpiERR
  !   integer                                  :: R,Q,Nini,Nfin
  !   call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  !   call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  !   call MPI_BCAST(vin,Ndim,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,mpiERR)
  !   vout=cmplx(0.d0,0.d0,8)
  !   vtmp=cmplx(0.d0,0.d0,8)
  !   ! R=Ndim/mpiSIZE ; Q=mod(Ndim,mpiSIZE)
  !   ! Nini=mpiID*R + 1            !1,R+1,2R+1,...
  !   ! Nfin=(mpiID+1)*R            !R,2R,3R
  !   ! if(Q/=0 .AND. mpiID==mpiSIZE-1)Nfin=Nfin+Q
  !   do i=mpiID*Q+1,(mpiID+1)*Q+R!Nini,Nfin
  !      c => sparse%row(i)%root%next       
  !      matmul: do
  !         if(.not.associated(c))exit matmul
  !         vtmp(i) = vtmp(i) + c%val*vin(c%col)
  !         c => c%next
  !      end do matmul
  !   end do
  !   call MPI_ALLREDUCE(vtmp,vout,Ndim,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mpiERR)
  ! end subroutine sp_matrix_vector_product_c_mpi
  !!>MPI

end module MATRIX_SPARSE
