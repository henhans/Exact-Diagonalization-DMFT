module EIGEN_SPACE
  implicit none
  private

  type full_espace
     real(8),dimension(:),pointer   :: e
     real(8),dimension(:,:),pointer :: M
  end type full_espace

  type sparse_estate
     integer                         :: sector
     real(8)                         :: e
     real(8),dimension(:),pointer    :: vec
     complex(8),dimension(:),pointer :: cvec
     logical                         :: isreal=.true.
     type(sparse_estate),pointer     :: next !link to next box (chain)
  end type sparse_estate

  type sparse_espace
     integer                     :: size
     logical                     :: status=.false.
     type(sparse_estate),pointer :: root !head/root of the list\== list itself
  end type sparse_espace


  interface es_insert_state
     module procedure es_insert_state_d,es_insert_state_c
  end interface es_insert_state

  public :: full_espace
  public :: sparse_estate
  public :: sparse_espace
  !
  public :: es_init_espace      !init the espace                 !checked
  public :: es_delete_espace    !del the espace                  !checked
  public :: es_free_espace      !free the espace                 !checked
  public :: es_print_espace     !print the espace                !checked
  !
  public :: es_insert_state     !insert a state                  !checked
  public :: es_delete_state     !remove a state                  !checked
  !
  public :: es_pop_size         !get the size of a state         !checked
  public :: es_pop_type         !get the type of a state         !checked
  public :: es_pop_sector       !get the sector of a state       !checked
  public :: es_pop_energy       !get the energy of a state       !checked
  public :: es_pop_vector       !get the vector of a state       !checked
  public :: es_pop_cvector      !get the vector of a state       !checked
  !

contains        !some routine to perform simple operation on the lists


  !+------------------------------------------------------------------+
  !PURPOSE  : initialize the list of states
  !+------------------------------------------------------------------+
  function es_init_espace() result(space)
    type(sparse_espace) :: space
    allocate(space%root)
    space%status=.true.
    space%root%next => null()
    space%size=0
  end function es_init_espace



  !+------------------------------------------------------------------+
  !PURPOSE  : destroy the list of states
  !+------------------------------------------------------------------+
  subroutine es_delete_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       deallocate(c)
    end do
    deallocate(space%root)
  end subroutine es_delete_espace




  !+------------------------------------------------------------------+
  !PURPOSE  : empty the list of states
  !+------------------------------------------------------------------+
  subroutine es_free_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
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
  !PURPOSE  : pretty print the list of states
  !+------------------------------------------------------------------+
  subroutine es_print_espace(space,unit)
    type(sparse_espace),intent(in) :: space
    type(sparse_estate),pointer    :: c
    integer                        :: counter,i
    integer,optional               :: unit
    integer                        :: unit_
    unit_=6;if(present(unit))unit_=unit
    write(*,*)"Print sparse espace ->",unit_
    c => space%root%next   !assume is associated,ie list exists
    counter = 0
    if(space%size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          write(unit_,"(A10,I5)")"Element : ",counter
          write(unit_,"(A10,I5)")"Sector  : ",c%sector
          if(c%isreal)then
             write(unit_,"(A10,I5)")"Size    : ",size(c%vec)
             write(unit_,"(A10,f18.9)")"Energy  : ",c%e
             write(unit_,"(A10)")"Vec     : "
             write(unit_,*)c%vec
          else
             write(unit_,"(A10,I5)")"size    : ",size(c%vec)
             write(unit_,"(A10,f18.9)")"Energy: ",c%e
             write(unit_,"(A10)")"Vec     : "
             write(unit_,*)c%cvec
          endif
          c => c%next  !traverse list
          write(unit_,*)""
       end do
    else
       write(unit_,*)"Empty space"
       return
    endif
  end subroutine es_print_espace



  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_insert_state_d(space,e,vec,sector)
    type(sparse_espace),intent(inout) :: space
    real(8),intent(in)                :: e
    real(8),dimension(:),intent(in)   :: vec
    integer,intent(in)                :: sector
    type(sparse_estate),pointer       :: p,c
    p => space%root
    c => p%next
    do                            !traverse the list until obj < value (ordered list)
       if(.not.associated(c))exit !empty list of beginning of the list
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    p%next%e = e
    allocate(p%next%vec(size(vec)))
    p%next%vec = vec
    p%next%isreal=.true.
    p%next%sector=sector
    space%size = space%size+1
    !
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine es_insert_state_d

  subroutine es_insert_state_c(space,e,vec,sector)
    type(sparse_espace),intent(inout)  :: space
    real(8),intent(in)                 :: e
    complex(8),dimension(:),intent(in) :: vec
    integer,intent(in)                 :: sector
    type(sparse_estate),pointer        :: p,c
    p => space%root
    c => p%next
    do                            !traverse the list until e < value (ordered list)
       if(.not.associated(c))exit
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    p%next%e = e
    allocate(p%next%cvec(size(vec)))
    p%next%cvec = vec
    p%next%isreal=.false.
    p%next%sector=sector
    space%size = space%size+1
    !
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       p%next%next  => null()
    else
       p%next%next  => c      !the %next of the new node come to current
    end if
  end subroutine es_insert_state_c





  !+------------------------------------------------------------------+
  !PURPOSE  : remove last element from the list, if +n is given remove 
  ! the n-th element, if +e is given remove the state with state%e=e
  !+------------------------------------------------------------------+
  subroutine es_delete_state(space,n)
    type(sparse_espace),intent(inout) :: space
    integer,optional,intent(in)       :: n
    integer                           :: i,pos
    type(sparse_estate),pointer       :: p,c
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)stop "es_delete_state: pos > espace.size"
    c => space%root
    do i=1,pos 
       p => c
       c => c%next
       if(.not.associated(c))return !empty or end of the list
    end do
    p%next => c%next !reallocate skipping the deleted link
    if(c%isreal)then
       deallocate(c%vec)
    else
       deallocate(c%cvec)
    endif
    deallocate(c)           !free link
    space%size=space%size-1
  end subroutine es_delete_state




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_pop_size(space,n) result(nsize)
    type(sparse_espace),intent(in) :: space
    integer,intent(in)             :: n
    type(sparse_estate),pointer    :: p,c
    integer                        :: i,nsize
    if(.not.space%status) stop "es_pop_size: espace not allocated"
    if(n>space%size)      stop "es_pop_size: n > espace.size"
    c => space%root
    do i=1,n
       c => c%next
       if(.not.associated(c))exit
    end do
    if(c%isreal)then
       nsize=size(c%vec)
    else
       nsize=size(c%cvec)
    endif
  end function es_pop_size



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_pop_sector(space,n) result(sector)
    type(sparse_espace),intent(in) :: space
    integer,intent(in)             :: n
    integer                        :: sector
    type(sparse_estate),pointer    :: c
    integer                        :: i
    if(.not.space%status) stop "es_pop_sector: espace not allocated"
    if(n>space%size)      stop "es_pop_sector: n > espace.size"
    c => space%root
    do i=1,n
       c => c%next
       if(.not.associated(c))exit
    end do
    sector = c%sector
  end function es_pop_sector



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_pop_type(space,n) result(ctype)
    type(sparse_espace),intent(in) :: space
    integer,intent(in)             :: n
    logical                        :: ctype
    type(sparse_estate),pointer    :: c
    integer                        :: i
    if(.not.space%status) stop "es_pop_type: espace not allocated"
    if(n>space%size)      stop "es_pop_type: n > espace.size"
    c => space%root
    do i=1,n
       c => c%next
       if(.not.associated(c))exit
    end do
    ctype = c%isreal
  end function es_pop_type



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_pop_energy(space,n) result(egs)
    type(sparse_espace),intent(in) :: space
    integer,intent(in)             :: n
    real(8)                        :: egs
    type(sparse_estate),pointer    :: c
    integer                        :: i
    c => space%root             !assume is associated,ie list exists
    egs=0.d0
    do i=1,n
       c => c%next  !traverse list
       if(.not.associated(c))exit
    end do
    egs = c%e
  end function es_pop_energy



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_pop_vector(space,n) result(vector)
    type(sparse_espace),intent(in) :: space
    integer,intent(in)             :: n
    real(8),dimension(:),pointer   :: vector
    type(sparse_estate),pointer    :: c
    integer                        :: i
    c => space%root   !assume is associated,ie list exists
    do i=1,n
       c => c%next  !traverse list
       if(.not.associated(c))exit
    end do
    if(.not.c%isreal)stop "pop vector isreal=F: can not associate to complex vector"
    vector => c%vec
  end function es_pop_vector

  function es_pop_cvector(space,n) result(vector)
    type(sparse_espace),intent(in)  :: space
    integer,intent(in)              :: n
    complex(8),dimension(:),pointer :: vector
    type(sparse_estate),pointer     :: c
    integer                         :: i
    c => space%root   !assume is associated,ie list exists
    do i=1,n
       c => c%next  !traverse list
       if(.not.associated(c))exit
    end do
    if(c%isreal)stop "pop cvector isreal=T: can not associate to real vector"
    vector => c%cvec
  end function es_pop_cvector





end module EIGEN_SPACE
