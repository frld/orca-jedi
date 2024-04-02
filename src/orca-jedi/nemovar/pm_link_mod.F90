module pm_link_mod
! Contains pm_link class and methods - link in linked list of key-value pairs
! where values are of type unlimited polymorphic
! Link object used in pm_linked_list_mod.F90
implicit none

private
public :: pm_link
type pm_link
  private
  integer :: key=0
  class(*), pointer :: object => null()
  type(pm_link),pointer :: next => null()
  contains
  procedure :: make_link
  procedure :: insert_link
  procedure :: get_key
  procedure :: get_next
  procedure :: point_object
  procedure :: delete_link
  procedure :: replace_link
end type pm_link
contains
!=========================================================================
subroutine make_link(self,key,object,next)
! Create a link type
class(pm_link),intent(inout) :: self
integer,intent(in) :: key
class(*),intent(in) :: object
type(pm_link),pointer,intent(in) :: next

self%key = key
self%next => next
allocate(self%object,source=object)

end subroutine make_link
!=========================================================================
subroutine insert_link(self,key,object,next)
! Create a link type
class(pm_link),intent(inout) :: self
integer,intent(in) :: key
class(*),pointer,intent(in) :: object
type(pm_link),pointer,intent(in) :: next

self%key = key
self%next => next
self%object => object

end subroutine insert_link
!=========================================================================
function get_key(self)
! Return Key of link
class(pm_link),intent(in) :: self
integer :: get_key
get_key = self%key
end function get_key
!=========================================================================
function get_next(self)
! Return self%next
class(pm_link),intent(in) :: self
type(pm_link),pointer :: get_next
get_next => self%next
end function get_next
!=========================================================================
function point_object(self)
! Return pointer to object in link
class(pm_link),intent(in) :: self
class(*),pointer :: point_object
point_object => self%object
end function point_object
!=========================================================================
subroutine delete_link(self,prev,head,ldealloc)
! Remove it self from linked list
class(pm_link),intent(inout) :: self
type(pm_link),pointer,optional,intent(inout) :: prev
type(pm_link),pointer,optional,intent(inout) :: head
logical,optional , intent(in) :: ldealloc
logical :: lldealloc

if(present(ldealloc)) lldealloc = ldealloc

!reconnect the list
if(present(prev)) then
  if (associated(prev)) then
    prev%next => self%next
  else
    head => self%next
  endif
endif

!remove the node and set key to 0
self%key = 0
self%next => null()
if(lldealloc) then
  if(associated(self%object)) deallocate(self%object)
endif

end subroutine delete_link
!=========================================================================
subroutine replace_link(self,object)
! Replace object in link
class(pm_link),intent(inout) :: self
class(*), pointer, intent(in) :: object

!YT cannot deallocate here, should return pointer and decide higher up
!yt if(associated(self%object)) deallocate(self%object)
self%object => object

end subroutine replace_link
!=========================================================================

end module pm_link_mod
