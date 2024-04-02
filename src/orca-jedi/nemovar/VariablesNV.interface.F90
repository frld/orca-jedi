! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module vars_mod

   use iso_c_binding
   use fckit_configuration_module, only: fckit_configuration
   use pm_linked_list_mod

   use nv_init_wrapper
   use variables

   implicit none

   private

   public &
      & get_vars   !: !: Access function to objects stored in vars_list

   ! ------------------------------------------------------------------------------

   !> Global registry

   type(pm_linked_list) :: vars_list

contains
   ! ------------------------------------------------------------------------------
   subroutine c_nv_variables_create( &
      & c_key_self, &
      & c_conf) bind(c,name='nv_variables_create_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_variables_create  ***
      !!
      !! ** Purpose : create variables object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      type(c_ptr),    intent(in)    :: c_conf

      type(fvariables_typ), pointer :: self
      type(fvariables_typ)          :: template
      integer :: ivar

      ! Local variables
      type(fckit_configuration) :: f_conf

      ! Interface
      f_conf = fckit_configuration(c_conf)

      ! Allocate variables object

      call vars_list%add(c_key_self, template)
      self => get_vars(c_key_self)

      ! Construct variables object

      !ivar = config_get_int(c_conf, "mode")
      call f_conf%get_or_die("mode",ivar)
      !!ivar = 0           ! DJL   (mode = 0 bg, mode=2 cov/cost fun?)
      !!ivar = 2
      write(6,*) 'c_nv_variables_create ivar / mode = ',ivar
      if (ivar < 0 .or. ivar > 4) call abor1_ftn("nv_variables_create: variables invalid")

      call self%set( piom_ctl, ivar)

   end subroutine c_nv_variables_create

   ! ------------------------------------------------------------------------------

   subroutine c_nv_variables_setup( &
      & c_key_self, &
      & c_var) bind(c,name='nv_variables_setup_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_variables_setup  ***
      !!
      !! ** Purpose : setup variables object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in)    :: c_var

      type(fvariables_typ), pointer :: self
      type(fvariables_typ)          :: template
      integer :: ivar

      ivar = c_var

      ! Allocate variables object

      call vars_list%add(c_key_self, template)
      self => get_vars(c_key_self)

      ! Construct variables object

      call self%set( piom_ctl, ivar)

   end subroutine c_nv_variables_setup

   ! ------------------------------------------------------------------------------

   subroutine c_nv_variables_clone( &
      & c_key_self, &
      & c_key_other ) bind(c,name='nv_variables_clone_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_variables_clone  ***
      !!
      !! ** Purpose : clone variables object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in)    :: c_key_other

      type(fvariables_typ), pointer :: self
      type(fvariables_typ)          :: template
      type(fvariables_typ), pointer :: other

      ! Allocate variables object

      call vars_list%add(c_key_self, template)
      self => get_vars(c_key_self)

      ! Retrieve variables object

      other => get_vars(c_key_other)

      ! Clone

      self = other

   end subroutine c_nv_variables_clone

   ! ------------------------------------------------------------------------------

   subroutine c_nv_variables_delete( &
      & c_key_self ) bind(c,name='nv_variables_delete_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_variables_delete  ***
      !!
      !! ** Purpose : delete variables object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self

      ! Remove variables object from the list

      call vars_list%remove(c_key_self)

   end subroutine c_nv_variables_delete

   ! ------------------------------------------------------------------------------

   function get_vars(c_self)
   integer(kind=c_int), intent(in) :: c_self
   type(fvariables_typ), pointer :: get_vars
   class(*), pointer :: this=>null()

   call vars_list%get(c_self, this)

   select type(self => this)
   type is (fvariables_typ)
      get_vars=>self
   class default
      call abor1_ftn('get_vars : unexpected type ')
   end select

   end function get_vars

   ! ------------------------------------------------------------------------------

end module vars_mod
