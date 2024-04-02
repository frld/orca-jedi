! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module fieldtlad_mod
   use iso_c_binding
   use pm_linked_list_mod

   use nv_init_wrapper

   use geom_mod
   use vars_mod

   use grid_def
   use nv_oce_def
   use nv_par_kind
   use par_def
   use variables

   implicit none

   private

   public &
      & get_fieldtlad, &          !: Access function to objects stored in fieldtlad_list
      & c_nv_fieldtlad_create, &
      & c_nv_fieldtlad_delete

   integer, parameter :: max_string = 200

   ! ------------------------------------------------------------------------------

   !> Global registry

   type(pm_linked_list) :: fieldtlad_list

contains
   ! ------------------------------------------------------------------------------

   subroutine c_nv_fieldtlad_create( &
      & c_key_self, &
      & c_key_geom, &
      & c_key_vars ) bind(c,name='nv_fieldtlad_create_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_fieldtlad_create  ***
      !!
      !! ** Purpose : construct field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in)    :: c_key_geom
      integer(c_int), intent(in)    :: c_key_vars

      type(foce_typ),       pointer :: self
      type(foce_typ)                :: template
      type(fgrid_typ),      pointer :: geom
      type(fvariables_typ), pointer :: vars

      ! Allocate field object

      call fieldtlad_list%add(c_key_self, template)
      self => get_fieldtlad(c_key_self)

      ! Retrieve objects

      geom => get_geom(c_key_geom)
      vars => get_vars(c_key_vars)

      ! Allocate a control vector

      call self%alloc( &
         & piom_ctl,   &
         & geom%fmpp,  &
         & vars )

      ! Zero out

      call self%zero( &
         & piom_ctl )

   end subroutine c_nv_fieldtlad_create

   ! ------------------------------------------------------------------------------

   subroutine c_nv_fieldtlad_delete(c_key_self) bind(c,name='nv_fieldtlad_delete_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_fieldtlad_delete  ***
      !!
      !! ** Purpose : delete a field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self

      type(foce_typ), pointer  :: self

      ! Retreive field object

      self => get_fieldtlad(c_key_self)

      ! Deallocate a control vector

      call self%dealloc( &
         & piom_ctl )

      ! Remove from list

      call fieldtlad_list%remove(c_key_self)

   end subroutine c_nv_fieldtlad_delete

   ! ------------------------------------------------------------------------------

   subroutine c_nv_fieldtlad_zero( &
      & c_key_self ) bind(c,name='nv_fieldtlad_zero_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_fieldtlad_zero  ***
      !!
      !! ** Purpose : set to zero a field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self

      type(foce_typ), pointer :: self

      ! Retrieve a field object

      self => get_fieldtlad(c_key_self)

      ! Set to zero

      call self%zero( &
         & piom_ctl )

   end subroutine c_nv_fieldtlad_zero

   ! ------------------------------------------------------------------------------

   function get_fieldtlad(c_self)
   integer(kind=c_int), intent(in) :: c_self
   type(foce_typ), pointer :: get_fieldtlad
   class(*), pointer :: this=>null()

   call fieldtlad_list%get(c_self, this)

   select type(self => this)
   type is (foce_typ)
      get_fieldtlad=>self
   class default
      call abor1_ftn('get_fieldtlad : unexpected type ')
   end select

   end function get_fieldtlad

   ! ------------------------------------------------------------------------------

end module fieldtlad_mod
