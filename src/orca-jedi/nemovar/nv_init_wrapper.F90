! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module nv_init_wrapper
!!==========================================================================
!!                       ***  module nv_init_wrapper ***
!!
!!==========================================================================
!!
!! Interface to call NV_INIT from C++
!!
!! This is only temporary to get going.
!! ************************************
!! **** NV_INIT should NOT EXIST *****
!! ************************************
!!
!! Hold instances of the following data types:
!!  * fwesp_typ
!!  * fvar_ctl_typ
!!  * fdfi_ctl_typ
!!  * fdfi_obs_typ
!!  * fiom_ctl_typ
!!  * fiom_file_dsc_typ
!!
!!==========================================================================
    !! * modules used

    use nv_asm_ctl_def   , only : fasm_ctl_typ
    use bal_ctl_def      , only : fbal_ctl_typ
    use bal_typ          , only : fbal_typ
    use bge_ctl_def      , only : fbge_ctl_typ
    use bkg_def          , only : fbkg_typ
    use dfi_ctl          , only : fdfi_ctl_typ
    use dfi_def          , only : fdfi_obs_typ
    use grid_def      , only : fgrid_typ, fmpp_typ, fclndr_typ, fctl_typ, fmsk_typ, fmetrics_typ
    use nv_in_out_manager, only : fiom_ctl_typ, nv_ctl_warn
    use nv_iom_def       , only : fiom_file_dsc_typ, nv_jpmax_files
    use nv_paresp        , only : fwesp_typ
!    use obsstore         , only : fobsstore_typ
    use obs_hand_def      , only : fobshand_typ
    use obe_sdv_ctl_def  , only : fobe_sdv_ctl_typ
    use var_ctl          , only : fvar_ctl_typ
    use variables_base

    use fckit_log_module

    implicit none

    !! * Global instances

    type(fasm_ctl_typ      )                                      :: pasm_ctl
    type(fbal_typ          )                            , pointer :: pbal
    type(fbal_ctl_typ      )                            , pointer :: pbal_ctl
    type(fgrid_typ         )                            , pointer :: pgrid
    type(fwesp_typ         )                                      :: pwesp
    type(fvar_ctl_typ      )                                      :: pvar_ctl
    type(fiom_ctl_typ      )                                      :: piom_ctl
    type(fiom_file_dsc_typ ), DIMENSION(nv_jpmax_files)           :: piom_file

    type(fobshand_typ      )                                      :: pobshand
    type(fobe_sdv_ctl_typ  )                                      :: pobe_sdv_ctl
    type(fbge_ctl_typ      )                                      :: pbge_ctl
    type(fbkg_typ          )                                      :: pbkg

    !! * Accessibility

    private

    public &
       & pasm_ctl,     &
       & pbal,         &
       & pbal_ctl,     &
       & pbge_ctl,     &
       & pbkg,         &
       & piom_ctl,     &
       & piom_file,    &
       & pgrid,        &
       & pobshand,     &
       & pobe_sdv_ctl, &
       & pwesp,        &
       & pvar_ctl

contains
! --------------------------------------------------------------------------

    subroutine nv_init_c( &
       & c_conf) bind(c,name='nv_init_f90')
       !!----------------------------------------------------------------------
       !!                    ***  ROUTINE nv_init_c  ***
       !!
       !! ** Purpose : initialize oops-nemovar
       !!
       !!----------------------------------------------------------------------
       use iso_c_binding
       use fckit_configuration_module

       implicit none

       type(c_ptr), intent(in) :: c_conf
       character(len=:),allocatable :: expv
       type(fckit_configuration) :: yconfig

       write(6,*) 'nv_init_c 1'
       call flush(6)
       yconfig = fckit_configuration(c_conf)

       write(6,*) 'nv_init_c 2'
       call flush(6)

       expv = "dan1";
       !if (.not. yconfig%get("expver",expv)) call exit

       write(6,*) 'nv_init_c 3'
       call flush(6)

       call nv_init(expv)

       return

    end subroutine nv_init_c

! --------------------------------------------------------------------------

    subroutine nv_finalize_c() bind(c,name='nv_finalize_f90')
       !!----------------------------------------------------------------------
       !!                    ***  ROUTINE nv_init_c  ***
       !!
       !! ** Purpose : finalize oops-nemovar
       !!
       !!----------------------------------------------------------------------
       use iso_c_binding

       implicit none

       call nv_finalize()
       call mpi_finalize()

       return

    end subroutine nv_finalize_c

! --------------------------------------------------------------------------

    subroutine nv_init( &
       & expv )
       !!-------------------------------------------------------------------
       !!                    ***  ROUTINE NV_INIT  ***
       !!
       !! ** Purpose : Initialize NEMOVAR
       !!
       !! ** Method  :
       !!
       !! ** Action  :
       !!
       !!-------------------------------------------------------------------
       use grid_mppini      , only : grid_mpp_init
       use nemovar_common
       use nv_in_out_manager, only : nv_iom_nml
       use nv_iom_def       , only : fiom_file_dsc_typ, nv_jpmax_files
       use nv_lib_mpp       , only : lbc_map_setup
       use nv_phycst        , only : nv_phy_cst
       use nv_timing        , only : nv_timing_init

       character(len=4), intent(in) :: expv

       !! * Local declarations

       ! * Initialize mpi

       write(6,*) 'nv_init 1'
       call flush(6)

       call mpp_init( &
          & piom_ctl )

       ! * Initialize logging/namelis

       write(6,*) 'nv_init 2'
       call flush(6)

       call nemovar_log_init( &
          & piom_ctl )

       ! * Print NEMOVAR banner

       write(6,*) 'nv_init 3'
       call flush(6)

       call nemovar_banner( &
          & piom_ctl )

       ! * Initialize physical constants

       write(6,*) 'nv_init 4'
       call flush(6)

       call nv_phy_cst( &
          & piom_ctl )

       ! * Initialize IO namelist

       write(6,*) 'nv_init 5'
       call flush(6)

       call nv_iom_nml( &
          & piom_ctl )

      ! * Initialise variables

       write(6,*) 'nv_init 5a'
       call flush(6)

       CALL gvars_metadata%init_metadata( &
          &  piom_ctl )

       ! * Initialize pointers to variables

       write(6,*) 'nv_init 6'
       call flush(6)

       call gvars_metadata%set_metadata( &
          & piom_ctl )

       write(6,*) 'nv_init 7'
       call flush(6)

       call pvar_ctl%init( &
          & piom_ctl )

       ! * Hacks

       write(6,*) 'nv_init 8'
       call flush(6)

       if (.not. piom_ctl%lwp) then
          call log%reset()
       endif

       return

    end subroutine nv_init

    subroutine nv_finalize()
       !!-------------------------------------------------------------------
       !!                    ***  ROUTINE NV_FINALIZE  ***
       !!
       !! ** Purpose : Finalize NEMOVAR
       !!
       !! ** Method  :
       !!
       !! ** Action  :
       !!
       !!-------------------------------------------------------------------

    end subroutine nv_finalize

! --------------------------------------------------------------------------


end module

