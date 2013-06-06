module EIGEN_SPACE
  implicit none
  private

  type,public :: eigenspace
     real(8),dimension(:),pointer   :: e
     real(8),dimension(:,:),pointer :: M
  end type eigenspace

  type,public :: eig_object
     integer                      :: sector
     real(8)                      :: e
     real(8),dimension(:),pointer :: vec
  end type eig_object

  type,public :: eig_state
     private
     type(eig_object)         :: obj  !object value of the node (content of the box)
     type(eig_state),pointer :: next !link to prev box (chain)
  end type eig_state

  type,public :: eig_space     
     integer                 :: size
     logical                 :: status=.false.
     type(eig_state),pointer :: root !head/root of the list\== list itself
  end type eig_space

  interface es_insert_state
     module procedure insert_element_obj,insert_element_components
  end interface es_insert_state

  interface es_remove_state
     module procedure remove_element_obj,remove_element_components
  end interface es_remove_state

  interface operator(<)
     module procedure less_than_obj
  end interface operator(<)

  interface operator(<=)
     module procedure less_or_equal_than_obj
  end interface operator(<=)

  interface operator(>)
     module procedure greater_than_obj
  end interface operator(>)

  interface operator(>=)
     module procedure greater_or_equal_than_obj
  end interface operator(>=)

  interface operator(==)
     module procedure equal_to_obj
  end interface operator(==)


  public :: es_init_espace
  public :: es_destroy_espace
  public :: es_free_espace
  public :: es_print_espace
  !
  public :: es_insert_state
  public :: es_remove_state
  public :: es_get_size_state
  !
  public :: es_get_component
  public :: es_get_sector
  public :: es_get_energy
  public :: es_get_vector
  !
  public :: operator(<)
  public :: operator(<=)
  public :: operator(>)
  public :: operator(>=)
  public :: operator(==)

contains        !some routine to perform simple operation on the lists


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_init_espace() result(space)
    type(eig_space) :: space
    allocate(space%root)
    space%status=.true.
    space%root%next => null()
    space%size=0
  end function es_init_espace


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine es_destroy_espace(space)
    type(eig_space),intent(inout) :: space
    type(eig_state),pointer       :: p,c
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       deallocate(c)
    end do
    deallocate(space%root)
  end subroutine es_destroy_espace



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine es_free_espace(space)
    type(eig_space),intent(inout) :: space
    type(eig_state),pointer       :: p,c
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       deallocate(c)
    end do
    space%size=0
  end subroutine es_free_espace






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine insert_element_obj(space,obj)
    type(eig_space),intent(inout) :: space
    type(eig_object),intent(in)   :: obj
    type(eig_state),pointer       :: p,c
    p => space%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list of beginning of the list
       if(obj <= c%obj) exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    p%next%obj = obj
    space%size = space%size+1
    !
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine insert_element_obj





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine insert_element_components(space,e,vec,sector)
    type(eig_space),intent(inout)  :: space
    real(8),intent(in)             :: e
    real(8),dimension(:),intent(in):: vec
    integer,intent(in)             :: sector
    type(eig_object)               :: obj
    type(eig_state),pointer        :: p,c
    p => space%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list of beginning of the list
       if(e <= c%obj%e) exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    obj%e=e
    allocate(obj%vec(size(vec)))
    obj%vec=vec
    obj%sector=sector
    p%next%obj = obj
    space%size = space%size+1
    !
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine insert_element_components





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine remove_element_obj(space,obj,found)
    type(eig_space),intent(inout)        :: space
    type(eig_object),optional,intent(in) :: obj
    !
    integer                           :: i
    type(eig_state),pointer           :: p,c
    logical                           :: found
    p => space%root
    c => p%next
    found = .false.
    do 
       if(found .OR. .not.associated(c))return
       if(obj == c%obj)then
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
       space%size=space%size-1
    endif
  end subroutine remove_element_obj






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine remove_element_components(space,n,e,found)
    type(eig_space),intent(inout)        :: space
    integer,optional,intent(in)          :: n
    real(8),optional,intent(in)          :: e
    integer                           :: i,pos
    type(eig_state),pointer           :: p,c
    logical                           :: found
    pos= space%size ; if(present(n))pos=n
    p => space%root
    c => p%next
    found = .false.
    if(present(e))then
       do 
          if(found .OR. .not.associated(c))return
          if(e == c%obj%e)then
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
          space%size=space%size-1
       endif
    else
       do i=1,pos 
          if(.not.associated(c))return
          p => c
          c => c%next
       end do
       p%next => c%next !reallocate skipping the deleted link
       deallocate(c)           !free link
       space%size=space%size-1
    endif
  end subroutine remove_element_components





  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine es_print_espace(space)
    type(eig_space),intent(in) :: space
    type(eig_state),pointer      :: c
    integer                        :: counter
    c => space%root%next   !assume is associated,ie list exists
    counter = 0
    if(space%size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          write(*,"(A,I5,A)",advance='no')"element: ",counter," |"
          call print_obj(c%obj)
          c => c%next  !traverse list
       end do
    else
       write(*,*)"Empty space"
       return
    endif
  end subroutine es_print_espace






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_get_size_state(space,n) result(nsize)
    type(eig_space),intent(in) :: space
    integer,intent(in)         :: n
    type(eig_state),pointer    :: p,c
    integer                    :: i,nsize
    c => space%root%next
    nsize=0
    do i=1,n
       if(.not.associated(c))return
       nsize=size(c%obj%vec)
       c => c%next
    end do
  end function es_get_size_state



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_get_component(space,n) result(obj)
    type(eig_space),intent(in)   :: space
    integer,intent(in)           :: n
    type(eig_object)            :: obj
    type(eig_state),pointer      :: c
    integer                      :: i
    c => space%root%next   !assume is associated,ie list exists
    do i=1,n
       if(.not.associated(c))exit
       obj = c%obj
       c => c%next  !traverse list
    end do
  end function es_get_component



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_get_sector(space,n) result(sector)
    type(eig_space),intent(in)   :: space
    integer,intent(in)           :: n
    integer                      :: sector
    type(eig_state),pointer      :: c
    integer                      :: i
    c => space%root%next   !assume is associated,ie list exists
    do i=1,n
       if(.not.associated(c))exit
       sector = c%obj%sector
       c => c%next  !traverse list
    end do
  end function es_get_sector




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_get_energy(space,n) result(egs)
    type(eig_space),intent(in)   :: space
    integer,intent(in)           :: n
    real(8)                      :: egs
    type(eig_state),pointer      :: c
    integer                      :: i
    c => space%root%next   !assume is associated,ie list exists
    do i=1,n
       if(.not.associated(c))exit
       egs = c%obj%e
       c => c%next  !traverse list
    end do
  end function es_get_energy



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_get_vector(space,n) result(vector)
    type(eig_space),intent(in)   :: space
    integer,intent(in)           :: n
    real(8),dimension(:),pointer :: vector
    type(eig_state),pointer      :: c
    integer                      :: i
    c => space%root   !assume is associated,ie list exists
    do i=1,n
       c => c%next  !traverse list
       if(.not.associated(c))exit
    end do
    vector => c%obj%vec
  end function es_get_vector




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function less_than_obj(O1,O2) result(boolean)
    type(eig_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%e < O2%e
  end function less_than_obj
  !
  function less_or_equal_than_obj(O1,O2) result(boolean)
    type(eig_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%e <= O2%e
  end function less_or_equal_than_obj
  !
  function greater_than_obj(O1,O2) result(boolean)
    type(eig_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%e > O2%e
  end function greater_than_obj
  !
  function greater_or_equal_than_obj(O1,O2) result(boolean)
    type(eig_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%e >= O2%e
  end function greater_or_equal_than_obj
  !
  function equal_to_obj(O1,O2) result(boolean)
    type(eig_object),intent(in) :: O1,O2
    logical                      :: boolean
    boolean = O1%e == O2%e
  end function equal_to_obj
  !
  subroutine print_obj(obj)
    type(eig_object) :: obj
    integer          :: i
    write(*,"(A10,f18.9)")"Energy: ",obj%e
    do i=1,size(obj%vec)
       write(*,"(A25,f18.9)")"Vec:",obj%vec(i)
    enddo
  end subroutine print_obj


end module EIGEN_SPACE
