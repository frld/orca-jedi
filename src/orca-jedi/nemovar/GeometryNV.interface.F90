! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module geom_mod

   use iso_c_binding
   use fckit_configuration_module, only: fckit_configuration
   use pm_linked_list_mod

   use grid_mppini      , only : grid_mpp_init
   use nemovar_common
   use nv_asminc_tam
   use grid_def
   use metrics
   use nv_init_wrapper
   use nv_in_out_manager, only : fiom_ctl_typ, nv_ctl_warn
   use nv_lib_mpp       , only : lbc_map_setup
   use nv_mpp_map       , only : nv_mppmap_init
   use nv_timing        , only : nv_timing_init

   use fckit_log_module

   implicit none

   private

   public &
      & get_geom    !: Access function to objects stored in geom_list

   ! ------------------------------------------------------------------------------

   !> global registry

   type(pm_linked_list) :: geom_list

contains
   ! ------------------------------------------------------------------------------
   subroutine c_nv_geo_setup( &
      & c_key_self, &
      & c_conf ) bind(c,name='nv_geo_setup_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_geo_setup  ***
      !!
      !! ** Purpose : construct geometry object
      !!
      !!----------------------------------------------------------------------
      use nv_lib_mpp, only : lbc_map_setup
      use nv_domcfg, only : nv_dom_cfg
      use nv_domain, only : nv_dom_init
      use coast, only : coast_init, coast_bkg

      implicit none

      integer(c_int), intent(inout) :: c_key_self
      type(c_ptr)   , intent(in)    :: c_conf

      type(fgrid_typ), pointer :: self
      type(fgrid_typ)          :: template

      ! Local variables
      !type(fckit_configuration) :: f_conf
      !
      ! Interface
      !f_conf = fckit_configuration(c_conf)

      ! Allocate geometry object

      write(6,*) 'nv_geo_setup start'
      call flush(6)

      call geom_list%add(c_key_self, template)
      self => get_geom(c_key_self)

      write(6,*) 'nv_geo_setup 1'
      call flush(6)
      
      ! mc> temporary hack
      pgrid => self
      ! mc< temporary hack

      write(6,*) 'nv_geo_setup 2'
      call flush(6)

      call grid_read_ctl_nml( &
         & piom_ctl, &
         & self )

      write(6,*) 'nv_geo_setup 3'
      call flush(6)
      
      call grid_mpp_nml_init( &
         & piom_ctl, &
         & self )

      write(6,*) 'nv_geo_setup 4'
      call flush(6)

      if( piom_ctl%nn_timing == 1 ) then

         call nv_timing_init( &
            & piom_ctl, &
            & self%fmpp )

      endif

      write(6,*) 'nv_geo_setup 5'
      call flush(6)

      call alloc_glob_geom( &
         & self%fmpp, &
         & self%fmsk, &
         & self%fctl, &
         & piom_ctl )

      write(6,*) 'nv_geo_setup 6'
      call flush(6)
      
      call grid_mpp_init( &
         & piom_ctl,  &
         & piom_file, &
         & self )

      write(6,*) 'nv_geo_setup 7'
      call flush(6)
      
      call lbc_map_setup( &
         & piom_ctl, &
         & self%fmpp )

!      call nv_asm_inc_tam_init( &
!         & pasm_ctl, &
!         & piom_ctl, &
!         & self )

      write(6,*) 'nv_geo_setup 8'
      call flush(6)
      
      if (.not. piom_ctl%lwp) then
         call log%reset()
      endif

      write(6,*) 'nv_geo_setup 9'
      call flush(6)
      
      call geom_alloc( &
         & self,    &
         & piom_ctl )

      write(6,*) 'nv_geo_setup 10'
      call flush(6)
      
      call nv_dom_cfg( &
         & piom_ctl, &
         & self )

      write(6,*) 'nv_geo_setup 11'
      call flush(6)
      
      call nv_dom_init( &
         & piom_ctl,  &
         & piom_file, &
         & self )

      write(6,*) 'nv_geo_setup 11a'
      call flush(6)

      CALL msk3d_ipt_init( &
         &  self%fmetrics, &
         &  piom_ctl, &
         &  self )

      write(6,*) 'nv_geo_setup 12'
      call flush(6)

      call coast_init( &
         & piom_ctl,  &
         & piom_file, &
         & self,      &
         & .false. )

      write(6,*) 'nv_geo_setup 13'
      call flush(6)
      
      call coast_bkg( &
         & piom_ctl,  &
         & piom_file, &
         & self,      &
         & .true. )

      write(6,*) 'nv_geo_setup 14'
      call flush(6)
      
      call nv_mppmap_init( &
         & piom_ctl, &
         & self%fctl, &
         & self%fmpp )

      write(6,*) 'nv_geo_setup end'
      call flush(6)

   end subroutine c_nv_geo_setup

   ! ------------------------------------------------------------------------------

   subroutine c_nv_geo_delete( &
      & c_key_self ) bind(c,name='nv_geo_delete_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_geo_delete  ***
      !!
      !! ** Purpose : delete a geometry object
      !!
      !!----------------------------------------------------------------------
      use coast, only : coast_dealloc, coast_bkg_dealloc

      implicit none

      integer(c_int), intent(inout) :: c_key_self

      type(fgrid_typ), pointer :: self

      ! Retrieve geometry object

      write(6,*) 'c_nv_geo_delete get geometry object'

      self => get_geom(c_key_self)

      ! Deallocate

      write(6,*) 'c_nv_geo_delete dealloc coast'

      call coast_dealloc( &
         & self%fcoast, &
         & piom_ctl )

      write(6,*) 'c_nv_geo_delete dealloc bkg'

      call coast_bkg_dealloc( &
         & self%fcoast, &
         & piom_ctl )

      write(6,*) 'c_nv_geo_delete dealloc geom'

      call geom_dealloc( &
         & self, &
         & piom_ctl )

      write(6,*) 'c_nv_geo_delete dealloc glob_geom'

      call dealloc_glob_geom( &
         & self%fmpp,     &
         & self%fmsk,     &
         & self%fmetrics, &
         & piom_ctl )

      ! Remove from list

      write(6,*) 'c_nv_geo_delete remove from list'

      call geom_list%remove( c_key_self )

      write(6,*) 'c_nv_geo_delete complete'

   end subroutine c_nv_geo_delete

   ! ------------------------------------------------------------------------------

   subroutine geom_alloc( &
      & self, &
      & piom_ctl )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE geom_alloc  ***
      !!
      !! ** Purpose : allocate geometry object
      !!
      !!----------------------------------------------------------------------
      type(fgrid_typ) , intent(inout) :: self
      type(fiom_ctl_typ), intent(inout) :: piom_ctl
      !!----------------------------------------------------------------------
      integer :: jpi, jpj, jpk
      !!
      if( self%ll_alloc .eqv. .true. ) then
        call nv_ctl_warn( piom_ctl, 'geom_alloc: fields already allocated')
        return
      endif
      !
      if( self%f2dh%ll_alloc .eqv. .true. ) then
        call nv_ctl_warn( piom_ctl, 'geom_alloc: f2dh: fields already allocated')
        return
      endif
      !
      if( self%f1dz%ll_alloc .eqv. .true. ) then
        call nv_ctl_warn( piom_ctl, 'geom_alloc: f1dz: fields already allocated')
        return
      endif
      !
      if( self%f1dt%ll_alloc .eqv. .true. ) then
        call nv_ctl_warn( piom_ctl, 'geom_alloc: f1dt: fields already allocated')
        return
      endif
      !

      jpi = self%fmpp%jpi
      jpj = self%fmpp%jpj
      jpk = self%fmpp%jpk
      !
      call self%f2dh%alloc( piom_ctl, jpi, jpj, jpk )
      call self%f1dz%alloc( piom_ctl, jpi, jpj, jpk )

      write(6,*) 'DJL geom_alloc f1dz ', jpi, jpj, jpk

      call self%f1dt%alloc( piom_ctl,           jpk )
      !
      self%ll_alloc = .true.
      write(6,*) 'DJL geom_alloc complete ll_alloc ',self%ll_alloc
      !
   end subroutine geom_alloc

   ! ------------------------------------------------------------------------------

   subroutine geom_dealloc( &
      & self, &
      & piom_ctl )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE geom_alloc  ***
      !!
      !! ** Purpose : allocate geometry object
      !!
      !!----------------------------------------------------------------------
      type(fgrid_typ) , intent(inout) :: self
      type(fiom_ctl_typ), intent(inout) :: piom_ctl
      !!----------------------------------------------------------------------
      !!
      if( self%ll_alloc .eqv. .false. ) then
        call nv_ctl_warn( piom_ctl, 'geom_dealloc: fields already deallocated')
        return
      endif
      !
      if( self%f2dh%ll_alloc .eqv. .false. ) then
        call nv_ctl_warn( piom_ctl, 'geom_dealloc: f2dh: fields already deallocated')
        return
      endif
      !
      if( self%f1dz%ll_alloc .eqv. .false. ) then
        call nv_ctl_warn( piom_ctl, 'geom_dealloc: f1dz: fields already deallocated')
        return
      endif
      !
      if( self%f1dt%ll_alloc .eqv. .false. ) then
        call nv_ctl_warn( piom_ctl, 'geom_dealloc: f1dt: fields already deallocated')
        return
      endif
      !
      call self%f2dh%dealloc( piom_ctl )
      call self%f1dz%dealloc( piom_ctl )
      call self%f1dt%dealloc( piom_ctl )
      !
      self%ll_alloc = .false.
      !
   end subroutine geom_dealloc

   ! ------------------------------------------------------------------------------

    subroutine alloc_glob_geom( &
      & pmpp, &
      & pmsk, &
      & pctl, &
      & piom_ctl )
      !!----------------------------------------------------------------------
      !!                   ***  subroutine alloc_glob_geom ***
      !!----------------------------------------------------------------------

      type(fmpp_typ), intent(inout)     :: pmpp
      type(fmsk_typ), intent(inout)     :: pmsk
      type(fctl_typ), intent(inout)     :: pctl
      type(fiom_ctl_typ), intent(inout) :: piom_ctl
      !!----------------------------------------------------------------------
      integer :: jpi, jpj, jpk, jpnij, jpidta, jpjdta, jpiglo, jpjglo
      !!
      if( pmpp%ll_alloc .eqv. .true. ) then
         call nv_ctl_warn( piom_ctl, 'alloc_glob_geom: fmpp: fields already allocated')
         return
      endif
      !
      if( pmsk%ll_alloc .eqv. .true. ) then
         call nv_ctl_warn( piom_ctl, 'alloc_glob_geom: fmsk: fields already allocated')
         return
      endif

      jpi = pmpp%jpi
      jpj = pmpp%jpj
      jpk = pmpp%jpk
      jpnij = pmpp%jpnij
      jpidta = pctl%jpidta
      jpjdta = pctl%jpjdta
      jpiglo = pctl%jpiglo
      jpjglo = pctl%jpjglo
      !

      write(6,*) 'DJL alloc pmmp jpidta, jpjdta, jpiglo, jpjglo ',jpidta, jpjdta, jpiglo, jpjglo
      call pmpp%alloc( piom_ctl, jpidta, jpjdta, jpiglo, jpjglo )
      write(6,*) 'DJL alloc pmsk jpi, jpj, jpk, jpiglo ',jpi, jpj, jpk, jpiglo
      call pmsk%alloc( piom_ctl, jpi, jpj, jpk, jpiglo )

      write(6,*) 'DJL alloc_glob_geom complete'

      !
    end subroutine alloc_glob_geom

   ! ------------------------------------------------------------------------------

    subroutine dealloc_glob_geom( &
      & pmpp, &
      & pmsk, &
      & pmetrics, &
      & piom_ctl )
      !!----------------------------------------------------------------------
      !!                   ***  subroutine alloc_glob_geom ***
      !!----------------------------------------------------------------------

      type(fmpp_typ), intent(inout)     :: pmpp
      type(fmsk_typ), intent(inout)     :: pmsk
      type(fmetrics_typ), intent(inout) :: pmetrics
      type(fiom_ctl_typ), intent(inout) :: piom_ctl

      !!
      if( pmpp%ll_alloc .eqv. .false. ) then
         call nv_ctl_warn( piom_ctl, 'dealloc_glob_geom: fmpp: fields already deallocated')
         return
      endif
      !
      if( pmsk%ll_alloc .eqv. .false. ) then
         call nv_ctl_warn( piom_ctl, 'dealloc_glob_geom: fmsk: fields already deallocated')
         return
      endif
      !
      call pmpp%dealloc( piom_ctl )
      call pmsk%dealloc( piom_ctl )
      call pmetrics%dealloc( piom_ctl )
      !
    end subroutine dealloc_glob_geom

   ! ------------------------------------------------------------------------------


   function get_geom(c_self)
   integer(kind=c_int), intent(in) :: c_self
   type(fgrid_typ), pointer :: get_geom
   class(*), pointer :: this=>null()

   call geom_list%get(c_self, this)

   select type(self => this)
   type is (fgrid_typ)
      get_geom=>self
   class default
      call abor1_ftn('get_geom : unexpected type ')
   end select

   end function get_geom

   ! ------------------------------------------------------------------------------

end module geom_mod

