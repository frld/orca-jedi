! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module field_mod

   use iso_c_binding
   use fckit_configuration_module, only: fckit_configuration
   use pm_linked_list_mod

   use nv_init_wrapper

   use fieldtlad_mod
   use geom_mod
!   use goms_mod
!   use locs_mod
!   use obsspace_mod
   use vars_mod

   use bal_opt
   use control_vectors
!   use obs_inter_opt
   use gen_opt
   use inc_def
   use locations
   use modelAtLocations
   use nv_asminc_tam
   use grid_def
   use nv_lbclnk
   use nv_lbclnk_tam
   use nv_oce_def
   use nv_par_kind
!   use obsstore
   use par_def
   use variables
   use variables_base
   
   implicit none

   private

   public &
      & get_field       !: Access function to objects stored in field_list

   ! ------------------------------------------------------------------------------

   integer, parameter :: max_string = 200

   !> Global registry

   type(pm_linked_list) :: field_list

contains
   ! ------------------------------------------------------------------------------
   subroutine c_nv_field_create( &
      & c_key_self, &
      & c_key_geom, &
      & c_key_vars ) bind(c,name='nv_field_create_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_create  ***
      !!
      !! ** Purpose : construct field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in)    :: c_key_geom
      integer(c_int), intent(in)    :: c_key_vars

      type(ctlvec),         pointer :: self
      type(ctlvec)                  :: template
      type(fgrid_typ),      pointer :: geom
      type(fvariables_typ), pointer :: vars

      ! Allocate field object

      call field_list%add(c_key_self, template)
      self => get_field(c_key_self)

      ! Retrieve objects

      geom => get_geom(c_key_geom)
      vars => get_vars(c_key_vars)

      ! Allocate a control vector

      write(6,*) 'DJL c_nv_field_create: alloc_ctlvec start' 

      call alloc_ctlvec( &
         & piom_ctl,  &
         & geom%fmpp, &
         & self,      &
         & vars,      &
         & 'dx' )        !: this needs to be made more flexible

      write(6,*) 'DJL c_nv_field_create: alloc_ctlvec done' 


   end subroutine c_nv_field_create

   ! ------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------
   subroutine c_nv_field_set( &
      & c_key_self, &
      & len, &
      & array ) bind(c,name='nv_field_set_f90')
      !!
      !!      DJL
      !!
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in) :: len
!!      real(c_double), intent(in) :: array(:)
      real(c_double), intent(in) :: array(len)

      type(ctlvec),         pointer :: self
          
      ! Retrieve a field object

      write(6,*) 'retrieve a field object c_key_self ',c_key_self
      call flush(6)

      self => get_field(c_key_self)

      write(6,*) 'set field to array'
      call flush(6)
      write(6,*) 'size(self%pdata(:) ', size(self%pdata(:))
      call flush(6)
      write(6,*) 'size(array(:) ', size(array(:))
      call flush(6)

      write(6,*) 'array(:) max min', maxval(array(:)), minval(array(:))
      call flush(6)

!!      self%pdata(:) = array(:)
      self%pdata(1:len) = array(1:len)
      write(6,*) 'set field to array complete'
      call flush(6)

   end subroutine c_nv_field_set

   ! ------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------
   subroutine c_nv_field_get( &
      & c_key_self, &
      & len, &
      & array ) bind(c,name='nv_field_get_f90')
      !!
      !!      DJL
      !!
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(inout) :: len
!!      real(c_double), intent(out) :: array(:)
      real(c_double), intent(out) :: array(len)

      type(ctlvec),         pointer :: self
            
      ! Retrieve a field object

      write(6,*) 'retrieve a field object c_key_self ',c_key_self
      call flush(6)

      self => get_field(c_key_self)

      write(6,*) 'get array from field'
      call flush(6)

      write(6,*) 'size(self%pdata(:) ', size(self%pdata(:))
      call flush(6)
      write(6,*) 'size(array(:) ', size(array(:))
      call flush(6)

      write(6,*) 'self%pdata(1:len) max min ',maxval(self%pdata(1:len)), minval(self%pdata(1:len))
      call flush(6)
      write(6,*) 'array(1:len) b4 max min ',maxval(array(1:len)), minval(array(1:len))
      call flush(6)

!!      array(:) = self%pdata(:)
      array(1:len) = self%pdata(1:len)
      write(6,*) 'set array from field complete'
      call flush(6)

   end subroutine c_nv_field_get
       

   subroutine c_nv_field_delete(c_key_self) bind(c,name='nv_field_delete_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_delete  ***
      !!
      !! ** Purpose : delete a field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self

      type(ctlvec), pointer  :: self

      ! Retreive field object

      self => get_field(c_key_self)

      ! Deallocate a control vector

      call dealloc_ctlvec( &
         & piom_ctl, &
         & self )

      ! Remove from list

      call field_list%remove(c_key_self)

   end subroutine c_nv_field_delete

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_zero( &
      & c_key_self ) bind(c,name='nv_field_zero_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_zero  ***
      !!
      !! ** Purpose : set to zero a field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self

      type(ctlvec), pointer :: self

      ! Retrieve a field object

      write(6,*) 'retrieve a field object c_key_self ',c_key_self
      call flush(6)

      self => get_field(c_key_self)

      ! Set to zero

      write(6,*) 'set field to zero'
      call flush(6)

      self%pdata(:) = 0.0_wp
      write(6,*) 'set field to zero complete'
      call flush(6)

   end subroutine c_nv_field_zero

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_copy( &
      & c_key_self, &
      & c_key_rhs ) bind(c,name='nv_field_copy_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_copy  ***
      !!
      !! ** Purpose : copy a field object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_rhs

      type(ctlvec), pointer  :: self
      type(ctlvec), pointer  :: rhs

      ! Retrieve field objects

      self => get_field(c_key_self)
      rhs => get_field(c_key_rhs)

      ! Copy

      self = rhs

   end subroutine c_nv_field_copy

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_self_add( &
      & c_key_self, &
      & c_key_other ) bind(c,name='nv_field_self_add_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_add  ***
      !!
      !! ** Purpose : add field objects
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in)    :: c_key_other

      type(ctlvec), pointer  :: self
      type(ctlvec), pointer  :: other

      ! Retrieve field objects

      self => get_field(c_key_self)
      other => get_field(c_key_other)

      ! Add

      self%pdata(:) = self%pdata(:) + other%pdata(:)   !: implement a procedure in ctlvec

   end subroutine c_nv_field_self_add

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_self_sub( &
      & c_key_self, &
      & c_key_other ) bind(c,name='nv_field_self_sub_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_sub  ***
      !!
      !! ** Purpose : subtract field objects
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_other

      type(ctlvec), pointer  :: self
      type(ctlvec), pointer  :: other

      ! Retrieve field objects

      self => get_field(c_key_self)
      other => get_field(c_key_other)

      ! Subtract

      self%pdata(:) = self%pdata(:) - other%pdata(:)   !: implement a procedure in ctlvec

   end subroutine c_nv_field_self_sub

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_self_mul( &
      & c_key_self, &
      & c_zz ) bind(c,name='nv_field_self_mul_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_mul  ***
      !!
      !! ** Purpose : multiply a field objects by a scalar
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      real(c_double), intent(in) :: c_zz

      type(ctlvec), pointer  :: self
      real(kind=wp) :: zz

      ! Retrieve a field object

      self => get_field(c_key_self)

      !
      zz = c_zz

      ! Multiply

      self%pdata(:) = zz * self%pdata(:)         !: implement a procedure in ctlvec

   end subroutine c_nv_field_self_mul

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_axpy( &
      & c_key_self,  &
      & c_key_other, &
      & c_zz ) bind(c,name='nv_field_axpy_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_axpy  ***
      !!
      !! ** Purpose : axpy
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_other
      real(c_double), intent(in) :: c_zz

      type(ctlvec), pointer :: self
      type(ctlvec), pointer :: other
      real(kind=wp) :: zz

      ! Retrieve objects

      self => get_field(c_key_self)
      other => get_field(c_key_other)
      zz =  c_zz

      ! AXPY

      self%pdata(:) = self%pdata(:) + zz * other%pdata(:)  !: use blas axpy

   end subroutine c_nv_field_axpy

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_dot_prod( &
      & c_key_self,  &
      & c_key_other, &
      & c_key_geom, &
      & c_zz ) bind(c,name='nv_field_dot_prod_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_dot_prod  ***
      !!
      !! ** Purpose : dot product
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in)  :: c_key_self
      integer(c_int), intent(in)  :: c_key_other
      integer(c_int), intent(in)  :: c_key_geom
      real(c_double), intent(out) :: c_zz

      type(ctlvec), pointer :: self
      type(ctlvec), pointer :: other
      type(fgrid_typ), pointer :: geom
      real(kind=wp) :: zz

      ! Retrieve objects

      self => get_field(c_key_self)
      other => get_field(c_key_other)

      geom => get_geom(c_key_geom)

      ! Calculate dot product

      zz = dot_product( geom%fmpp, self, other )

      c_zz = zz

   end subroutine c_nv_field_dot_prod

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_norm( &
      & c_key_self, &
      & c_key_geom, &
      & c_zz ) bind(c,name='nv_field_norm_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_norm  ***
      !!
      !! ** Purpose : calculate norm
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in)  :: c_key_self
      integer(c_int), intent(in)  :: c_key_geom
      real(c_double), intent(out) :: c_zz

      type(ctlvec), pointer :: self
      type(fgrid_typ), pointer :: geom
      real(kind=wp) :: zz

      ! Retrieve a field object

      self => get_field(c_key_self)

      geom => get_geom(c_key_geom)

      ! Calculate norm

      zz = norm( geom%fmpp, self )

      c_zz = zz

   end subroutine c_nv_field_norm

   ! ------------------------------------------------------------------------------

!   subroutine c_nv_field_interp_tl( &
!      & c_key_self_tl, &
!      & c_key_geom, &
!      & c_key_locs, &
!      & c_key_goms, &
!      & c_iter ) bind(c,name='nv_field_interp_tl_f90')
!      !!----------------------------------------------------------------------
!      !!                    ***  ROUTINE c_nv_interp_tl  ***
!      !!
!      !! ** Purpose : interpolation operator
!      !!
!      !!----------------------------------------------------------------------
!      implicit none
!
!      integer(c_int), intent(in)    :: c_key_self_tl
!      integer(c_int), intent(in)    :: c_key_geom
!      integer(c_int), intent(in)    :: c_key_locs
!      integer(c_int), intent(in)    :: c_key_goms
!      integer(c_int), intent(in)    :: c_iter
!
!      type(foce_typ),              pointer :: self_tl
!      type(fgrid_typ),             pointer :: geom
!      type(flocations_typ),        pointer :: locs
!      type(fmodelAtLocations_typ), pointer :: goms
!
!      integer :: j_stp, j_ext
!      integer :: j_iter
!
!      ! Retrieve objects
!
!      self_tl => get_fieldtlad(c_key_self_tl)
!      geom => get_geom(c_key_geom)
!      locs => get_locs(c_key_locs)
!      goms => get_goms(c_key_goms)
!
!      j_iter = c_iter
!
!      do j_stp = j_iter, j_iter
!
!          if( j_stp == 0 ) then
!
!              call obs_count_tan_init( &
!                 & pobshand, &
!                 & locs%locs ) ! TEMPORARY
!
!          endif
!
!          call obs_inter_tan( &
!             & locs,     &
!             & goms,     &
!             & pvar_ctl, &
!             & self_tl, &
!             & pobshand, &
!             & piom_ctl, &
!             & geom,     &
!             & j_stp )
!
!          if( j_stp == piom_ctl%nn_itend ) THEN
!
!              call goms%getFromODB( &
!                 & piom_ctl,  &
!                 & locs%locs, &
!                 & pobshand,  &
!                 & 'tan' )  ! TEMPORARY
!
!          endif
!
!      end do
!
!   end subroutine c_nv_field_interp_tl

   ! ------------------------------------------------------------------------------

!   subroutine c_nv_field_interp_ad( &
!      & c_key_self_ad, &
!      & c_key_geom,    &
!      & c_key_locs,    &
!      & c_key_goms,    &
!      & c_iter ) bind(c,name='nv_field_interp_ad_f90')
!      !!----------------------------------------------------------------------
!      !!                    ***  ROUTINE c_nv_interp_ad  ***
!      !!
!      !! ** Purpose : adjoint of the interpolation operator
!      !!
!      !!----------------------------------------------------------------------
!      implicit none
!
!      integer(c_int), intent(in)    :: c_key_self_ad
!      integer(c_int), intent(in)    :: c_key_geom
!      integer(c_int), intent(in)    :: c_key_locs
!      integer(c_int), intent(in)    :: c_key_goms
!      integer(c_int), intent(in)    :: c_iter
!
!      type(foce_typ),              pointer :: self_ad
!      type(fgrid_typ),             pointer :: geom
!      type(flocations_typ),        pointer :: locs
!      type(fmodelAtLocations_typ), pointer :: goms
!
!      integer :: j_stp, j_ext
!      integer :: j_iter
!
!      ! Retrieve objects
!
!      self_ad => get_fieldtlad(c_key_self_ad)
!      geom => get_geom(c_key_geom)
!      locs => get_locs(c_key_locs)
!      goms => get_goms(c_key_goms)
!
!      j_iter = c_iter
!
!      ! Allocate working types
!
!      do j_stp = j_iter, j_iter
!
!          if( j_stp == piom_ctl%nn_itend ) THEN
!
!              call obs_count_adj_init( &
!                 & pobshand, &
!                 & locs%locs ) ! TEMPORARY
!
!          endif
!
!          call obs_inter_adj( &
!             & locs,     &
!             & goms,     &
!             & pvar_ctl, &
!             & self_ad,  &
!             & pobshand, &
!             & piom_ctl, &
!             & geom,     &
!             & j_stp )
!
!      end do
!
!   end subroutine c_nv_field_interp_ad

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_totlfields( &
      & c_key_self,    &
      & c_key_self_tl, &
      & c_key_geom,    &
      & c_key_vars,    &
      & c_skip_bal ) bind(c,name='nv_field_totlfields_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_totlfields  ***
      !!
      !! ** Purpose : transform from incr fields to tl fields
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(inout) :: c_key_self_tl
      integer(c_int), intent(in)    :: c_key_geom
      integer(c_int), intent(in)    :: c_key_vars
      integer(c_int), intent(in)    :: c_skip_bal

      type(ctlvec),         pointer :: self
      type(foce_typ),       pointer :: self_tl
      type(fgrid_typ),      pointer :: geom
      type(fvariables_typ), pointer :: vars
      type(fvar_typ),       pointer :: z_var => NULL()
      logical                       :: ll_skip_bal

      type(finc_typ) :: z_inc_tl
      integer :: j_ext
      integer :: i_var
      integer :: j_var
      REAL(KIND=wp) :: z_sgn


      ! Retrieve objects

      self => get_field(c_key_self)
      geom => get_geom(c_key_geom)
      call c_nv_fieldtlad_create( &
         & c_key_self_tl, &
         & c_key_geom,    &
         & c_key_vars )
      self_tl => get_fieldtlad(c_key_self_tl)

      if (c_skip_bal == 1) then
        ll_skip_bal = .true.
      else
        ll_skip_bal = .false.
      endif

      ! Reset the ocean tangent variables

      call z_inc_tl%reset( &
         &  piom_ctl,   &
         &  geom%fmpp, &
         &  self%vars )

      if (.NOT. ll_skip_bal) then

         call z_inc_tl%getFromCtlvec( &
            & self,     &
            & piom_ctl, &
            & geom%fmpp )

         ! Apply the balance operator

         call bal_lbc_tan_opt( &
            & piom_ctl,          &
            & self%vars,         &
            & pbal_ctl,          &
            & pbal,              &
            & geom,              &
            & z_inc_tl%vars3d(vars%nid_t)%var(:,:,:), &
            & z_inc_tl%vars3d(vars%nid_s)%var(:,:,:), &
            & z_inc_tl%vars3d(vars%nid_u)%var(:,:,:), &
            & z_inc_tl%vars3d(vars%nid_v)%var(:,:,:), &
            & z_inc_tl%vars2d(vars%nid_ssh)%var(:,:)  )

         ! Apply the lateral boundary conditions, if not applied during the balance

         DO j_var = 1, vars%n2d

            i_var = vars%getVar2d(j_var)

            IF ( vars%nid_ssh == i_var ) CYCLE  ! lbc_lnk already applied

            IF ( vars%ln_in_ctlvec_2d(i_var) ) THEN

               z_var => get_var_2d( i_var )

               CALL nv_lbc_lnk( &
                  &  piom_ctl, &
                  &  geom%fmpp, &
                  &  z_inc_tl%vars2d(i_var)%var(:,:), &
                  &  z_var%cgrid, &
                  &  1.0_wp )

            END IF

         END DO

         DO j_var = 1, vars%n3d

            i_var = vars%getVar3d(j_var)

            IF ( vars%nid_t == i_var ) CYCLE  ! lbc_lnk already applied
            IF ( vars%nid_s == i_var ) CYCLE
            IF ( vars%nid_u == i_var ) CYCLE
            IF ( vars%nid_v == i_var ) CYCLE

            IF ( vars%ln_in_ctlvec_3d(i_var) ) THEN

               z_sgn = 1.0_wp

               IF (      TRIM( vars%c3d_fldname(j_var) ) == 'u' &
                  & .OR. TRIM( vars%c3d_fldname(j_var) ) == 'v' ) z_sgn = -1.0_wp

               z_var => get_var_3d( i_var )

               CALL nv_lbc_lnk( &
                  &  piom_ctl, &
                  &  geom%fmpp, &
                  &  z_inc_tl%vars3d(i_var)%var(:,:,:), &
                  &  z_var%cgrid, &
                  &  z_sgn )

            END IF

         END DO

         if ( piom_ctl%lwp ) THEN
            write(piom_ctl%numout,*)
            write(piom_ctl%numout,*) ' Max/min of the initial increment'
            write(piom_ctl%numout,*) ' -------------------------------- '
         endif

         call gen_diag_tan( &
            &  z_inc_tl, &
            &  piom_ctl, &
            &  vars, &
            &  geom )

      endif

      if ( pasm_ctl%ln_asmdin .OR. ( pvar_ctl%nvarex == jp_inc_3dvar ) ) THEN

         ! Direct initialization

         call nv_asm_inc_tan( &
            &  pasm_ctl, &
            &  z_inc_tl, &
            &  self_tl,  &
            &  piom_ctl, &
            &  vars, &
            &  pbal_ctl, &
            &  piom_ctl%nn_it000 - 1 )

      endif

   end subroutine c_nv_field_totlfields

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_toincrfields( &
      & c_key_self,    &
      & c_key_self_ad, &
      & c_key_geom,    &
      & c_key_vars ) bind(c,name='nv_field_toincrfields_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_incfields  ***
      !!
      !! ** Purpose : transform from tl fields to increment fields
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(inout) :: c_key_self_ad
      integer(c_int), intent(in)    :: c_key_geom
      integer(c_int), intent(in)    :: c_key_vars

      type(ctlvec),         pointer :: self
      type(foce_typ),       pointer :: self_ad
      type(fgrid_typ),      pointer :: geom
      type(fvariables_typ), pointer :: vars
      type(fvar_typ),       pointer :: z_var => NULL()

      type(finc_typ) :: z_inc_ad
      integer :: j_ext
      integer :: i_var
      integer :: j_var
      REAL(KIND=wp) :: z_sgn
      
      ! Retrieve objects

      self => get_field(c_key_self)
      self_ad => get_fieldtlad(c_key_self_ad)
      geom => get_geom(c_key_geom)

      ! Reset the adjoint variables

      call z_inc_ad%reset( &
         &  piom_ctl, &
         &  geom%fmpp, &
         &  self%vars )

     if ( pasm_ctl%ln_asmdin .OR. ( pvar_ctl%nvarex == jp_inc_3dvar ) ) THEN

         ! Direct initialization

         call nv_asm_inc_adj( &
            &  pasm_ctl, &
            &  z_inc_ad, &
            &  self_ad,  &
            &  piom_ctl, &
            &  vars, &
            &  pbal_ctl, &
            &  piom_ctl%nn_it000 - 1 )

      endif

      ! Apply the lateral boundary conditions, if they will not be applied during the balance

      DO j_var = 1, vars%n2d

         i_var = vars%getVar2d(j_var)

         IF ( vars%nid_ssh == i_var ) CYCLE  ! lbc_lnk will be applied

         IF ( vars%ln_in_ctlvec_2d(i_var) ) THEN

            z_var => get_var_2d( i_var )

            CALL nv_lbc_lnk( &
               &  piom_ctl, &
               &  geom%fmpp, &
               &  z_inc_ad%vars2d(i_var)%var(:,:), &
               &  z_var%cgrid, &
               &  1.0_wp )

         END IF

      END DO

      DO j_var = 1, vars%n3d

         i_var = vars%getVar3d(j_var)

         IF ( vars%nid_t == i_var ) CYCLE  ! lbc_lnk will be applied
         IF ( vars%nid_s == i_var ) CYCLE
         IF ( vars%nid_u == i_var ) CYCLE
         IF ( vars%nid_v == i_var ) CYCLE

         IF ( vars%ln_in_ctlvec_3d(i_var) ) THEN

            z_sgn = 1.0_wp

            IF (      TRIM( vars%c3d_fldname(j_var) ) == 'u' &
               & .OR. TRIM( vars%c3d_fldname(j_var) ) == 'v' ) z_sgn = -1.0_wp

            z_var => get_var_3d( i_var )

            CALL nv_lbc_lnk( &
               &  piom_ctl, &
               &  geom%fmpp, &
               &  z_inc_ad%vars3d(i_var)%var(:,:,:), &
               &  z_var%cgrid, &
               &  z_sgn )

         END IF

      END DO

      ! Apply the balance operator

      call bal_lbc_tan_opt( &
         & piom_ctl,          &
         & self%vars,         &
         & pbal_ctl,          &
         & pbal,              &
         & geom,              &
         & z_inc_ad%vars3d(vars%nid_t)%var(:,:,:), &
         & z_inc_ad%vars3d(vars%nid_s)%var(:,:,:), &
         & z_inc_ad%vars3d(vars%nid_u)%var(:,:,:), &
         & z_inc_ad%vars3d(vars%nid_v)%var(:,:,:), &
         & z_inc_ad%vars2d(vars%nid_ssh)%var(:,:)  )

      call z_inc_ad%setCtlvec( &
         & self,     &
         & piom_ctl, &
         & geom%fmpp )

      call c_nv_fieldtlad_delete( &
         & c_key_self_ad )

   end subroutine c_nv_field_toincrfields

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_add_incr( &
      & c_key_self, &
      & c_key_other ) bind(c,name='nv_field_add_incr_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_add_incr  ***
      !!
      !! ** Purpose : add increment to the state vector
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(in)    :: c_key_other

      type(ctlvec), pointer :: self
      type(ctlvec), pointer :: other

      ! Retrieve objects

      self => get_field(c_key_self)
      other => get_field(c_key_other)

      ! Add increment

!!      self%pdata(:) = self%pdata(:) + other%pdata(:)

   end subroutine c_nv_field_add_incr

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_diff_incr( &
      & c_key_self, &
      & c_key_x1,   &
      & c_key_x2 ) bind(c,name='nv_field_diff_incr_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_diff_incr  ***
      !!
      !! ** Purpose : subtract increments
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(in)    :: c_key_x1
      integer(c_int), intent(in)    :: c_key_x2

      type(ctlvec), pointer :: self
      type(ctlvec), pointer :: x1
      type(ctlvec), pointer :: x2

      ! Retrieve objects

      self => get_field(c_key_self)
      x1 => get_field(c_key_x1)
      x2 => get_field(c_key_x2)

      !self%pdata(:) = x1%pdata(:) - x2%pdata(:)

   end subroutine c_nv_field_diff_incr

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_read_file( &
      & c_key_self, &
      & c_key_geom, &
      & c_conf ) bind(c,name='nv_field_read_file_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_read_file  ***
      !!
      !! ** Purpose : read from file
      !!
      !!----------------------------------------------------------------------
      use string_f_c_mod

      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_geom
      type(c_ptr), intent(in) :: c_conf

      character(len=max_string) :: file_name
      character(len=:),allocatable :: str

      type(ctlvec),    pointer :: self
      type(fgrid_typ), pointer :: geom

      ! Local variables
      type(fckit_configuration) :: f_conf

      ! Interface
      f_conf = fckit_configuration(c_conf)

      ! Retrieve objects

      !file_name = config_get_string(c_conf, len(file_name), "filename")
      call f_conf%get_or_die("filename",str)
      file_name = str
      self => get_field(c_key_self)
      geom => get_geom(c_key_geom)

      ! Read

      call self%read( &
         & piom_ctl,  &
         & piom_file, &
         & geom%fctl, &
         & geom%fmpp, &
         & trim(file_name) )

   end subroutine c_nv_field_read_file

   ! ------------------------------------------------------------------------------

   subroutine c_nv_field_write_file( &
      & c_key_self, &
      & c_key_geom, &
      & c_conf ) bind(c,name='nv_field_write_file_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_field_write_file  ***
      !!
      !! ** Purpose : write to a file
      !!
      !!----------------------------------------------------------------------
      use string_f_c_mod

      implicit none

      integer(c_int), intent(in)    :: c_key_self
      integer(c_int), intent(inout) :: c_key_geom
      type(c_ptr)   , intent(in)    :: c_conf

      character(len=max_string) :: file_name
      character(len=:),allocatable :: str

      type(ctlvec),    pointer :: self
      type(fgrid_typ), pointer :: geom

      ! Local variables
      type(fckit_configuration) :: f_conf

      ! Interface
      f_conf = fckit_configuration(c_conf)

      ! Retrieve objects

      !file_name = config_get_string(c_conf, len(file_name), "name")
      call f_conf%get_or_die("name",str)
      file_name = str
      self => get_field(c_key_self)
      geom => get_geom(c_key_geom)

      ! write to file

      call self%write( &
         & piom_ctl,  &
         & piom_file, &
         & geom%fctl, &
         & geom%fmpp, &
         & geom%f2dh, &
         & geom%f1dz, &
         & trim(file_name) )

   end subroutine c_nv_field_write_file

   ! ------------------------------------------------------------------------------

   function get_field(c_self)
   integer(kind=c_int), intent(in) :: c_self
   type(ctlvec), pointer :: get_field
   class(*), pointer :: this=>null()

   call field_list%get(c_self, this)

   select type(self => this)
   type is (ctlvec)
      get_field=>self
   class default
      call abor1_ftn('get_field : unexpected type ')
   end select

   end function get_field

   ! ------------------------------------------------------------------------------

end module field_mod
