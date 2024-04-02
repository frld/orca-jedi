module pm_linked_list_mod
! Contains class for linked list where the links contain key-value pairs. The value
! is of abstract type (can be anything).Key is integer returned when link is added.
use pm_link_mod
implicit none
private
integer :: global_count=0 ! Global count to make sure key for any list within executable
                          ! is unique (unless you use insert function)
public :: pm_linked_list
type pm_linked_list
  private
  integer                :: count  = 0  ! Incremented when link is added to list
  integer                :: nelem  = 0  ! Current number of links in list
  type(pm_link), pointer :: head   => NULL()
  logical                :: ldebug=.false.
  contains
  procedure :: add
  procedure :: get
  procedure :: remove
  procedure :: replace
  procedure :: finalize
  procedure :: insert
  procedure :: iterate_keys
  procedure :: set_debug
  procedure :: print
end type pm_linked_list
contains
!=========================================================================
subroutine set_debug(self,ldon)
! Return pointer to object in list given key
class(pm_linked_list), intent(inout) :: self
logical,               intent(in)    :: ldon

self%ldebug = ldon
end subroutine set_debug
!=========================================================================
subroutine add(self,key,object)
! Add link to linked list
class(pm_linked_list), intent(inout) :: self
integer              , intent(out)   :: key
class(*) ::  object

type(pm_link),pointer :: current

if(global_count == huge(global_count)) global_count = 0 ! this should hopefully not happen
global_count = global_count+1
write(6,*) 'DJL pm_linked_list global_count ',global_count
self%count = global_count
self%nelem = self%nelem+1
key = self%count ! Provides unique key as always incremented

allocate(current)
call current%make_link(key,object,self%head)
self%head => current
if(self%ldebug) call self%print('add',key)
end subroutine add
!=========================================================================
subroutine get(self,key,obj_ptr)
! Return pointer to object in list given key
class(pm_linked_list), intent(inout) :: self
integer              , intent(in)    :: key
class (*),pointer                    :: obj_ptr
type(pm_link), pointer :: current

current => self%head
obj_ptr => NULL()

 !sweep the linked list to find matching key
do while(associated(current))
  if (key == current%get_key()) then
    obj_ptr => current%point_object()
    exit
  endif
  current=>current%get_next()
enddo

if (.not.associated(obj_ptr)) then
  write(0,*) 'pm_linked_list:get key not found DJL1 ',key,self%nelem
  call abor1_ftn('pm_linked_list:get key not found DJL2 ')
endif
end subroutine get
!=========================================================================
subroutine print(self,cdstring,key)
! Return pointer to object in list given key
class(pm_linked_list), intent(inout) :: self
character(len=*) ,     intent(in)    :: cdstring
integer          ,     intent(in)    :: key
type(pm_link), pointer :: current
integer :: ii

current => self%head

 !sweep the linked list to find matching key
if(associated(current)) then
  write(20,*)'pm_linked_list:print ',cdstring,' key=',key,' HEAD-key=',current%get_key(),self%nelem
else
  write(20,*)'pm_linked_list:print ',cdstring, 'EMPTY',' key=',key,self%nelem
endif
ii = 0
do while(associated(current))
  ii = ii+1
  write(20,*) 'pm_linked_list:print ',ii,current%get_key()
  current=>current%get_next()
enddo

end subroutine print
!=========================================================================
subroutine remove(self,key,ldealloc)
! Remove link with matching key from list
class(pm_linked_list), intent(inout) :: self
integer          , intent(in) :: key
logical,optional , intent(in) :: ldealloc
type(pm_link), pointer :: elem
type(pm_link), pointer :: prev
logical :: lldealloc

lldealloc = .true.
if(present(ldealloc)) lldealloc = ldealloc
elem => self%head
prev => null()
!sweep the linked list to find matching key,
do while(associated(elem))
  if (key == elem%get_key()) exit
  prev => elem
  elem => elem%get_next()
enddo

if (associated(elem)) then
  call elem%delete_link(prev,self%head,lldealloc)
  if(lldealloc) deallocate(elem)
  self%nelem = self%nelem-1
endif

if(self%ldebug) call self%print('rem',key)


end subroutine remove
!=========================================================================
subroutine replace(self,key,obj_ptr)
! Replace object in list matching key
class(pm_linked_list), intent(inout) :: self
integer, intent(in) :: key
class(*), pointer, intent(in) :: obj_ptr

type(pm_link), pointer :: elem

elem => self%head
do while(associated(elem))
  if (key == elem%get_key()) exit
  elem => elem%get_next()
enddo

if (associated(elem)) then
  call elem%replace_link(obj_ptr)
else
  call abor1_ftn('pm_linked_list:replace key not found')
endif

end subroutine replace
!=========================================================================
subroutine finalize(self)
! Empty list
class(pm_linked_list), intent(inout) :: self
type(pm_link), pointer :: current


current => self%head
do while(associated(self%head))
  current => self%head
  self%head => self%head%get_next()
  call current%delete_link()
  deallocate(current)
enddo
self%nelem = 0
self%count = 0

end subroutine finalize
!=========================================================================
subroutine insert(self,key,object)
! Add link with given key to linked list
class(pm_linked_list), intent(inout) :: self
integer              , intent(in)   :: key
class(*), pointer, intent(in) ::  object

type(pm_link),pointer :: current

current => self%head
 !sweep the linked list to find matching key
do while(associated(current))
  if (key == current%get_key()) then
    call abor1_ftn('pm_linked_list:insert - key already in use')
  endif
  current=>current%get_next()
enddo
current => null()
self%count = self%count+1
self%nelem = self%nelem+1

allocate(current)
call current%insert_link(key,object,self%head)
self%head => current
end subroutine insert
!=========================================================================
logical function  iterate_keys(self,key)
! Iterate over keys in linked list
class(pm_linked_list), intent(inout) :: self
integer              , intent(inout) :: key
type(pm_link), pointer :: current

current => self%head
if(key /= 0) then
 !sweep the linked list to find matching key
  do while(associated(current))
    if (key == current%get_key())then
      current=>current%get_next()
      exit
    endif
    current=>current%get_next()
  enddo
endif

if(associated(current)) then
  key = current%get_key()
  if(key == 0) then
    iterate_keys = .false.
  else
    iterate_keys = .true.
  endif
else
  iterate_keys = .false.
  key = 0
endif

end function iterate_keys
!=========================================================================

end module pm_linked_list_mod
