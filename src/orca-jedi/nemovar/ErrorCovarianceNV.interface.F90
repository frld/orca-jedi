! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module bmat_mod

   use iso_c_binding
   use fckit_configuration_module, only: fckit_configuration
   use pm_linked_list_mod

   use field_mod
   use geom_mod
   use nv_init_wrapper

   use control_vectors
   use bal
   use bge
   use bge_cov
   use bge_ran
   use cov_def  ! DJL
   use cor_def  ! DJL
   use eof_def  ! DJL
   use dif_def  ! DJL
   use sol_def  ! DJL
   use jb_def
   use grid_def
   use par_def

   implicit none

   private

   public &
      & get_bmat    !: Access function to objects stored in bmat_list

   public dan_typ

   ! ------------------------------------------------------------------------------

   !> Global registry

   type(pm_linked_list) :: bmat_list

   TYPE dan_fdif_3d

      INTEGER :: &
         & npi,             &               !: Number of 3D grid points
         & npj,             &               !:
         & npk,             &               !:
         & ncov     = 0,    &               !: Total number of covariance models
         & nloc     = 0,    &               !: Total number of localization scales
                                            !:
         & ndiftyp,         &               !: Type of diffusion model
                                            !:   jp_dif_exp_2x1D_exp_1D = 2 x 1D explicit x 1D explicit
                                            !:   jp_dif_exp_2D_exp_1D   = 2D explicit     x 1D explicit
                                            !:   jp_dif_exp_2x1D_imp_1D = 2 x 1D explicit x 1D implicit
                                            !:   jp_dif_exp_2D_imp_1D   = 2D explicit     x 1D implicit
                                            !:   jp_dif_imp_2D_imp_1D   = 2D implicit     x 1D implicit
                                            !:   jp_dif_exp_3D          = 3D explicit
                                            !:   jp_dif_imp_3D          = 3D implicit
                                            !:   jp_dif_exp_2x1D_noz    = 2 x 1D explicit
                                            !:   jp_dif_exp_2D_noz      = 2D explicit
                                            !:   jp_dif_imp_2D_noz      = 2D implicit
                                            !:   jp_dif_noh_exp_1D      =                   1D explicit
                                            !:   jp_dif_noh_imp_1D      =                   1D implicit
                                            !:   jp_dif_exp_2x1D_pro_1D = 2 x 1D explicit x 1D projection
                                            !:   jp_dif_exp_2D_pro_1D   = 2D explicit     x 1D projection
                                            !:   jp_dif_imp_2D_pro_1D   = 2D implicit     x 1D projection
                                            !:   jp_dif_imp_1D_imp_2D   = 1D implicit     x 2D implicit
                                            !:
         & ndiften,         &               !: Type of diffusion tensor
                                            !:   jp_dif_ten_dia       = Diagonal diffusion tensor
                                            !:   jp_dif_ten_nondia_2D = Non-diagonal 2D (horizontal) diffusion tensor
                                            !:
         & ndifbcs,         &               !: Type of boundary conditions
                                            !:   jp_dif_bcs_Nh_Nz   = Neumann   for both horizontal and vertical diffusion
                                            !:   jp_dif_bcs_Dh_Dz   = Dirichlet for both horizontal and vertical diffusion
                                            !:   jp_dif_bcs_NDh_NDz = Average of Neumann and Dirichlet for both horizontal
                                            !:                        and vertical diffusion
                                            !:   jp_dif_bcs_Nh_NDz  = Neumann   for horizontal diffusion and
                                            !:                        average of Neumann and Dirichlet for vertical diffusion
                                            !:   jp_dif_bcs_Dh_NDz  = Dirichlet for horizontal diffusion and
                                            !:                        average of Neumann and Dirichlet for vertical diffusion
                                            !:
         & ndif,            &               !: Number of diffusion operations (depends on values of ndiftyp and ndifbcs)
         & ndifit_h,        &               !: Number of h-direction diffusion iterations
         & ndifit_z,        &               !: Number of z-direction diffusion iterations
         & nintstp,         &               !: Interleaving steps for 2x1D x 1D explicit diffusion: ( ndifit_h / 2  ) / nintstp must be an integer
         & nparstp,         &               !: Parallel steps for 2D and 3D implicit diffusion
         & nseqstp,         &               !: Sequential step indicator for 2D and 3D implicit diffusion
                                            !:
         & nten_h,          &               !: Type of horizontal diffusion tensor
                                            !:   jp_dif_ten_hor_cst     = Constant
                                            !:   jp_dif_ten_hor_grd     = Factor times the local grid elements
                                            !:   jp_dif_ten_hor_equ     = Equatorial parameterization
                                            !:   jp_dif_ten_hor_ros     = Rossby radius
                                            !:   jp_dif_ten_hor_rea     = Read tensor elements from file
                                            !:   jp_dif_ten_hor_lsc_rea = Read length scales from file
                                            !:
         & nten_z,          &               !: Type of vertical diffusion tensor for temperature
                                            !:   jp_dif_ten_ver_grd     = Factor times the local grid elements
                                            !:   jp_dif_ten_ver_mxl     = Mixed layer parameterization
                                            !:   jp_dif_ten_ver_ens_rea = Read from file (ens)
                                            !:   jp_dif_ten_ver_rea     = Read from file
                                            !:
         & nmxlten,         &               !: Type of mixed layer depth for vertical length scale parameterization
                                            !:   jp_dif_ten_ver_mxl_cst = Fixed depth (p_depfix_mxld)
                                            !:   jp_dif_ten_ver_mxl_bkg = Background mixed layer depth
                                            !:
         & nmxlup,          &               !: Type trigering either the writing or reading of the fixed mixed layer
                                            !: depth lookup table with normalizationfactors
                                            !:   jp_mxlup_dis = disabled = 0
                                            !:   jp_mxlup_bld = writing
                                            !:   jp_mxlup_rea = reading only normalization factors
                                            !:
         & mxl_fixlev,      &               !: Fixed mixed layer depth level to build the lookup table
                                            !: Only relevant in fgrid%fctl%ln_zps and zco vertical coordinates.
                                            !:
         & nrosten,         &               !: Type of Rossby radius for the horizontal length scale parameterization
                                            !:   jp_dif_ten_hor_ros_cst = Fixed phase speed
                                            !:   jp_dif_ten_hor_ros_bkg = Background dependent phase speed
                                            !:   jp_dif_ten_hor_ros_rea = Read in from file
                                            !:
         & nensten_h,       &               !: Type of horizontal diffusion tensor with ensembles
                                            !:   jp_dif_ten_hor_ens_dal     = Use directly the elements of the Daley tensor
                                            !:   jp_dif_ten_hor_ens_lps     = Diagonal: square of the length of the principal axes of the Daley tensor
                                            !:   jp_dif_ten_hor_ens_lps_avg = Diagonal: square of the geometrical average of the length
                                            !:                                of the principal axes of the Daley tensor
                                            !:   jp_dif_ten_hor_ens_dal_sca = Diagonal: Daley tensor scaled to preserve area of the local ellipse
                                            !:
         & nbctot_h,        &               !: Total number of BCs for horizontal diffusion
                                            !:   1 =    if Neumann   BC only        (ndifbcs = jp_dif_bcs_Nh_Nz, jp_dif_bcs_Nh_NDz)
                                            !:       or if Dirichlet BC only        (ndifbcs = jp_dif_bcs_Dh_Dz, jp_dif_bcs_Dh_NDz)
                                            !:   2 =    if Neumman and Dirichlet BC (ndifbcs = jp_dif_bcs_NDh_NDz)
                                            !:
         & nbctot_z,        &               !: Total number of BCs for vertical diffusion
                                            !:   1 =    if Neumann   BC only        (ndifbcs = jp_dif_bcs_Nh_Nz)
                                            !:       or if Dirichlet BC only        (ndifbcs = jp_dif_bcs_Dh_Dz)
                                            !:   2 =    if Neumman and Dirichlet BC (ndifbcs = jp_dif_bcs_NDh_NDz, jp_dif_bcs_Nh_NDz, jp_dif_bcs_Dh_NDz)
                                            !:
         & nintsize3d,      &               !: Size of the interior domain
         & nranvec,         &               !: Random vector number (used with ndiftyp = jp_dif_imp_1D_imp_2D)
         & n_iom_unit,      &               !: IOM unit number (used with ndiftyp = jp_dif_imp_1D_imp_2D)
         & nitsha,          &               !: Number of Shapiro iterations on the fine grid in the transfer operator
         & nsamstp,         &               !: Sampling step on the fine grid used for testing purposes only
         & nten_h_clim_hyb, &               !: Type of horizontal diffusion tensor combined with climatological tensor
                                            !:   jp_dif_ten_hor_cst     = Constant
                                            !:   jp_dif_ten_hor_grd     = Factor times the local grid elements
                                            !:   jp_dif_ten_hor_equ     = Equatorial parameterization
                                            !:   jp_dif_ten_hor_ros     = Rossby radius
         & nten_z_clim_hyb, &               !: Type of vertical diffusion tensor combined with climatological tensor
                                            !:   jp_dif_ten_ver_grd = Factor times the local grid elements
                                            !:   jp_dif_ten_ver_mxl = Mixed layer parameterization
                                            !:
         & nlstyp                           !: Type of Matern/Gaussian function length-scale for diffusion in R^d
                                            !:   jp_dif_len_sca_dal = From Weaver and Mirouze (2013) derived from Daley (1991) general relation
                                            !:      D =           L sqrt( 2M - d - 2 )  with implicit scheme (Matern)
                                            !:      D =           L sqrt( 2M         )  with explicit scheme (Gaussian)
                                            !:   jp_dif_len_sca_ras = From Rasumssen and Williams (2006)
                                            !:      D =           L sqrt( 2M - d     )  with implicit scheme (Matern)
                                            !:      D =           L sqrt( 2M         )  with explicit scheme (Gaussian)
                                            !:   jp_dif_len_sca_ste = From Stein (1999)
                                            !:      D = sqrt( 2 ) L sqrt( 2M - d     )  with implicit scheme (Matern)
                                            !:      D = sqrt( 2 ) L sqrt( 2M         )  with explicit scheme (Gaussian)
                                            !:   jp_dif_len_sca_lin = From Lindgren, Rue and Lindstrom (2011)
                                            !:      D =       2   L sqrt( 2M - d     )  with implicit scheme (Matern)
                                            !:      D =       2   L sqrt( 2M         )  with explicit scheme (Gaussian)

      ! Coarse grid

      INTEGER :: &
         & ngrid = 1                        !: Number of grids for diffusion filtering of ensemble-estimated parameters

      INTEGER, DIMENSION(:), ALLOCATABLE :: &
         & nsolit_k                         !: Total number of solver iterations

      INTEGER, DIMENSION(jp_nimpmax) :: &
         & npartot                          !: Number of time-parallel diffusion steps for 2D or 3D implicit diffusion

      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & mxlup_mskt                       !: Lookup table mixed layer level mask

      CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: &
         & cbccomb                          !: Boundary condition combinations:
                                            !:  'NN', 'ND', 'DN' or 'DD'                                 if 2D x 1D
                                            !:  'NNN', 'NND', 'NDN', 'NDD', 'DNN', 'DND', 'DDN' or 'DDD' if 3 x 1D

      CHARACTER(LEN=jp_fld_name_len) :: &
         & cfld                             !: Field type identifier

      CHARACTER(LEN=1) :: &
         & cgrid                            !: Grid-point indicator: 'T', 'U', 'V', 'F'

      CHARACTER(LEN=jp_var_comments_max_len) :: &
         & cfldname                         !: Field name

      CHARACTER(LEN=14) :: &
         & cfile_prefix                     !: Prefix of filename for diffusion tensor I/O

      CHARACTER(LEN=11) :: &
         & cdiftyp                          !: Diffusion data structure type

      REAL(KIND=wp) :: &
         & sqtm1ndif,           &           !: Inverse of square root of number of diffusion applications
         & rm1ndif,             &           !: Inverse of number of diffusion applications
         & sgn_bc,              &           !: Sign change at lateral boundary condition: +1 or -1
         & sgn_lap_hz,          &           !: Sign of the 3D Laplacian operator
                                            !:  = +1 for explicit scheme
                                            !:  = -1 for implicit scheme
         & sgn_lap_h,           &           !: Sign of the 2D Laplacian operator
                                            !:  = +1 for explicit scheme
                                            !:  = -1 for implicit scheme
         & sgn_lap_z,           &           !: Sign of the 1D Laplacian operator
                                            !:  = +1 for explicit scheme
                                            !:  = -1 for implicit scheme
                                            !:
         & DtoK_hz,             &           !: Iteration number dependent scale factor needed to convert
                                            !: the scale tensor D into a diffusion tensor K
                                            !:  = 1 /   2 M_hz       for 3D explicit diffusion (horizontal-vertical)
                                            !:  = 1 / ( 2 M_hz - 5 ) for 3D implicit diffusion (horizontal-vertical)
         & KtoD_hz,             &           !: Inverse of DtoK_hz
                                            !:
         & DtoK_h,              &           !: Iteration number dependent scale factor needed to convert
                                            !: the scale tensor D into a diffusion tensor K
                                            !:  = 1 /   2 M_h       for 2D explicit diffusion (horizontal)
                                            !:  = 1 / ( 2 M_h - 3 ) for 1D implicit diffusion (horizontal)
                                            !:  = 1 / ( 2 M_h - 4 ) for 2D implicit diffusion (horizontal)
         & KtoD_h,              &           !: Inverse of DtoK_h
                                            !:
         & DtoK_z,              &           !: Iteration number dependent scale factor needed to convert
                                            !: the scale tensor D into a diffusion tensor K
                                            !:  = 1 /   2 M_z       for 1D explicit diffusion (vertical)
                                            !:  = 1 / ( 2 M_z - 3 ) for 1D implicit diffusion (vertical)
         & KtoD_z,              &           !: Inverse of DtoK_z
                                            !:
         & scagrd_z,            &           !: Factor times local vertical grid element
         & scagrd_z_mld,        &           !: Factor times local MLD
         & scamax_h,            &           !: Maximum allowed horizontal length scale (in degrees)
         & scamax_z,            &           !: Maximum allowed vertical length scale (in metres)
         & scamin_z,            &           !: Minimum allowed vertical length scale (in metres)
         & scaten_h,            &           !: Scale factor for horizontal diffusion tensor with localization
                                            !: and normalization factor filtering
         & scaten_z,            &           !: Scale factor for vertical diffusion tensor with localization
                                            !: and normalization factor filtering
         & rcoast,              &           !: Distance to the coastline from which horizontal tensor
                                            !: elements are reduced (in km)
         & depbnd_mxl,          &           !: Maximum mixed layer depth allowed for the vertical length scale bounds check
         & depfix_mxl,          &           !: Fixed mixed layer depth
         & depmax_mxl,          &           !: Maximum mixed layer depth allowed for the vertical length scale parameterization
         & wgtten_h,            &           !: Weighting parameter for horizontal diffusion tensor
         & wgtten_z,            &           !: Weighting parameter for vertical diffusion tensor
         & fraten_h,            &           !: Fraction of maximum of K11 (K22) for defining the amplitude
                                            !: of the random perturbations to the horizontal tensor elements
         & rcoast_msk,          &           !: Distance from coastline (in km) from which to reduce mask in optimal filter
         & grdmin_fac,          &           !: Factor for multiplying minimum grid space in coarse-grid check
         & rten_hyb_cutlat,     &           !: Cutoff latitude beyond which the parameterized tensor is linearly ramped to the
                                            !: full climatological tensor
         & rten_hyb_cutlat_bnd, &           !: Linear weighting band
         & rten_hyb_cutlat_wgt              !: Weight for the climatological tensor between the equator and the cutoff latitude

      REAL(KIND=wp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
         & difchokk                         !: Cholesky factor elements of the 1D vertical implicit diffusion matrix
                                            !:  difchokk(:,:,:,1,nbc) = Inverse of the diagonal elements
                                            !:  difchokk(:,:,:,2,nbc) = Off-diagonal elements

      REAL(KIND=wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & cho1, &                          !: Constant arrays (k,i,j) involving the Cholesky factors which are used
         & cho2, &                          !: in the forward elimination / back substitution 1D vertical diffusion solver
         & cho3                             !:

      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & K11,            &                !: 11-element of the diffusion tensor at (i+1/2,j    ,k    ) point if     diagonal tensor
                                            !: 11-element of the diffusion tensor at (i+1/2,j+1/2,k    ) point if non-diagonal tensor
         & K12,            &                !: 12-element of the diffusion tensor at (i    ,j+1/2,k    ) point if     diagonal tensor
                                            !: 12-element of the diffusion tensor at (i+1/2,j+1/2,k    ) point if non-diagonal tensor
         & K22,            &                !: 22-element of the diffusion tensor at (i    ,j+1/2,k    ) point if     diagonal tensor
                                            !: 22-element of the diffusion tensor at (i+1/2,j+1/2,k    ) point if non-diagonal tensor
         & K33,            &                !: 33-element of the diffusion tensor at (i    ,j    ,k+1/2) point
         & K11_fix,        &                !: K11 associated with the correlation matrix used to perturb the diffusion coefficients
                                            !: and standard deviations
         & K22_fix,        &                !: K22 associated with the correlation matrix used to perturb the diffusion coefficients
                                            !: and standard deviations
         & K33_min,        &                !: Minimum allowed K33 element of the diffusion tensor
         & detsqt_Kh_col,  &                !: Square root of the determinant of the horizontal diffusion tensor at the colocated field point
         & KP1,            &                !: Parallel      diffusion coefficient of the rotated non-diagonal tensor at (i+1/2,j+1/2,k) point
         & KP2,            &                !: Perpendicular diffusion coefficient of the rotated non-diagonal tensor at (i+1/2,j+1/2,k) point
         & L11,            &                !: 11-element of the localization tensor at (i,j,k) point
         & L22,            &                !: 22-element of the localization tensor at (i,j,k) point
         & L33,            &                !: 33-element of the localization tensor at (i,j,k) point
         & cofmatii,       &                !: Constant diffusion matrix coefficient for the ii-derivative at (i+1/2,j    ,k    ) point
         & cofmatjj,       &                !: Constant diffusion matrix coefficient for the jj-derivative at (i    ,j+1/2,k    ) point
         & cofmatkk,       &                !: Constant diffusion matrix coefficient for the kk-derivative at (i    ,j+1/2,k+1/2) point
         & cofmatii_nd,    &                !: Constant diffusion matrix coefficient for the ii-derivative at (i+1/2,j+1/2,k    ) point with non-diagonal tensor
         & cofmatij_nd,    &                !: Constant diffusion matrix coefficient for the ij-derivative at (i+1/2,j+1/2,k    ) point with non-diagonal tensor
                                            !:        = diffusion matrix coefficient for the ji-derivative
         & cofmatjj_nd,    &                !: Constant diffusion matrix coefficient for the jj-derivative at (i+1/2,j+1/2,k    ) point with non-diagonal tensor
                                            !: Constant symmetric diffusion matrix coefficients for a 5- or 7-point stencil
         & cofsym_ijk,     &                !:  - Coefficient at (i  ,j  ,k  ) point
         & cofsym_imjk,    &                !:  - Coefficient at (i-1,j  ,k  ) point
         & cofsym_ipjk,    &                !:  - Coefficient at (i+1,j  ,k  ) point
         & cofsym_ijmk,    &                !:  - Coefficient at (i  ,j-1,k  ) point
         & cofsym_ijpk,    &                !:  - Coefficient at (i  ,j+1,k  ) point
         & cofsym_ijkm,    &                !:  - Coefficient at (i  ,j  ,k-1) point used with full 3D implicit diffusion
         & cofsym_ijkp,    &                !:  - Coefficient at (i  ,j  ,k+1) point used with full 3D implicit diffusion
         & cofsym_imjmk,   &                !:  - Coefficient at (i-1,j-1,k  ) point used with a non-diagonal horizontal tensor
         & cofsym_imjpk,   &                !:  - Coefficient at (i-1,j+1,k  ) point used with a non-diagonal horizontal tensor
         & cofsym_ipjmk,   &                !:  - Coefficient at (i+1,j-1,k  ) point used with a non-diagonal horizontal tensor
         & cofsym_ipjpk                     !:  - Coefficient at (i+1,j+1,k  ) point used with a non-diagonal horizontal tensor

      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & K11_min,         &               !: Minimum allowed K11 element of the diffusion tensor
         & K11_col_K22_min, &               !: Minimum allowed K11 element of the diffusion tensor, colocated at the K22 point
         & K22_min,         &               !: Minimum allowed K22 element of the diffusion tensor
         & K22_col_K11_min                  !: Minimum allowed K11 element of the diffusion tensor, colocated at the K11 point

      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & wrkupd1                          !: Work array used with the symmetric diffusion matrix for 2Dx1D implicit diffusion

      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & wrkupd2                          !: Work array used with the symmetric diffusion matrix for 3D implicit diffusion

      REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: &
         & soltol_k                         !: Final relative residual norm

      REAL(KIND=wp), DIMENSION(jp_jpk) :: &
         & scalen_h, &                      !: Constant horizontal length scale (in degrees) in each model level
         & scalen_z, &                      !: Constant vertical   length scale (in metres)  in each model level
         & scagrd_h                         !: Factor times local horizontal grid element in each model level

      REAL(KIND=wp), DIMENSION(:,:,:), POINTER :: &
         & difmski     => NULL(), &         !: Masking coefficient multiplying the i-derivative at (i+1/2,j    ,k    ) point
         & difmskj     => NULL(), &         !: Masking coefficient multiplying the j-derivative at (i    ,j+1/2,k    ) point
         & difmskk     => NULL(), &         !: Masking coefficient multiplying the k-derivative at (i    ,j    ,k+1/2) point
         & difmskint   => NULL(), &         !: Interior point mask at (i,j,k) point
         & grd_vp1     => NULL(), &         !: Volume elements at (i,j,k) point
         & grd_sqt_vp1 => NULL(), &         !: Square root of volume elements at (i,j,k) point
         & grd_vm1     => NULL(), &         !: Inverse of volume elements at (i,j,k) point
         & grd_sqt_vm1 => NULL(), &         !: Square root of inverse of volume elements at (i,j,k) point
         & grd_dp1     => NULL()            !: Depth elements at (i,j,k) point

      REAL(KIND=wp), DIMENSION(:,:), POINTER :: &
         & grd_sqt_ap1 => NULL(), &         !: Square root of area elements at (i,j) point
         & grd_dp1_tot => NULL(), &         !: Total depth at interior points at (i,j) point
         & grd_dm1_tot => NULL()            !: Inverse of the total depth at (i,j) point

      REAL(KIND=wp), DIMENSION(:), POINTER :: &
         & grdmin_ap1  => NULL(), &         !: Minimum geometric mean of the horizontal scale factors at interior points
         & grdmax_ap1  => NULL(), &         !: Maximum geometric mean of the horizontal scale factors at interior points
         & grd_vm1_tot => NULL(), &         !: Inverse of the total volume of interior points
         & difmsksum   => NULL()            !: Sum of interior point masks in each level

      LOGICAL :: &
         & ltrfgrd,                       & !: Switch for transferring to a coarse grid in the diffusion operator
         & lsamgrd,                       & !: Switch for subsampling the fine grid for testing the transfer operator
         & lcoast_bnd,                    & !: Namelist switch for bounding the maximum horizontal diffusion tensor elements
                                            !: by the distance to the nearest coastline
         & lmxl_bnd,                      & !: Namelist switch for bounding the maximum vertical diffusion tensor elements
                                            !: in the mixed layer by the maximum allowed mixed layer depth (depmax_mxl)
         & ldep_bnd,                      & !: Namelist switch for bounding the maximum vertical diffusion tensor elements
                                            !: by the local ocean depth
         & lcoast,                        & !: Namelist switch for reducing horizontal diffusion tensor elements near coastlines
         & lranten,                       & !: Namelist switch for randomly perturbing the horizontal diffusion tensor elements
         & lens_ten_h,                    & !: Switch for ensemble-estimated background error horizontal tensor
         & lens_ten_z,                    & !: Switch for ensemble-estimated background error vertical   tensor
         & ldif_ten_rea        = .FALSE., & !: Switch indicating if           diffusion tensor has been read from correlation restart file
         & ldif_ten_per_rea    = .FALSE., & !: Switch indicating if perturbed diffusion tensor has been read from correlation restart file
         & ldif_h_ran_wri      = .FALSE., & !: Switch for writing to   file randomization vectors of the horizontal correlation model
         & ldif_h_ran_rea      = .FALSE., & !: Switch for reading from file randomization vectors of the horizontal correlation model
         & ldif_noh            = .FALSE., & !: Switch for skipping the horizontal diffusion in the randomization computation
         & lnorran_wr4         = .FALSE., & !: Switch for writing randomization vectors to file in single precision
         & lnor_nozxnoh        = .FALSE., & !: Switch for normalization for 3D fields approximated as the product of the normalization factors
                                            !: for no horizontal diffusion with the normalization factors for no vertical diffusion
         & ldif_ten_h_clim_hyb = .FALSE., & !: Switch for combining climatological horizontal diffusion tensor with the parameterized one
         & ldif_ten_z_clim_hyb = .FALSE., & !: Switch for combining climatological verticl    diffusion tensor with the parameterized one
         & lcofmat             = .FALSE., & !: Diffusion matrix coefficient initialization indicator
         & lcofsymmat          = .FALSE.    !: Symmetric diffusion matrix coefficient initialization indicator

      TYPE( fsol_3d ) :: &
         & fsol                             !: 2D elliptic solver type

      TYPE( fgrid_typ ), POINTER :: & 
         & fgrid => NULL()                  !: Grid definitions

   CONTAINS
!
      PROCEDURE :: assign_fdif_3d
      GENERIC   :: assignment(=) => assign_fdif_3d
!      FINAL     :: final_fdif_3d

   END TYPE dan_fdif_3d

   TYPE dan_typ                         ! DJL
      TYPE(dan_fdif_3d), DIMENSION(jp_var_max) :: &
         & array                !:  model arrays
!      TYPE(fdif_3d), DIMENSION(jp_var_max) :: &
!         & array                !:  model arrays
   END TYPE dan_typ

contains
   SUBROUTINE assign_fdif_3d( &
      & self, &
      & other )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE assign_fdif_3d  ***
      !!
      !! ** Purpose : Assign fdif_3d    DJL debugging!!
      !!
      !!----------------------------------------------------------------------
      !! * Arguments

      CLASS(dan_fdif_3d), INTENT(INOUT) :: &
         & self

      TYPE(dan_fdif_3d), INTENT(IN   ) :: &
         & other

      !! *

      INTEGER :: &
         & irank, &
         & ierror

      write(6,*) 'assign_fdif_3d DJL ',self%npi, self%npj, self%npk, self%ncov, self%nloc, self%ndiftyp, self%ndiften

      !!* Non-pointer members

      self%npi                 = other%npi
      self%npj                 = other%npj
      self%npk                 = other%npk
      self%ncov                = other%ncov
      self%nloc                = other%nloc
      self%ndiftyp             = other%ndiftyp
      self%ndiften             = other%ndiften
      self%ndifbcs             = other%ndifbcs
      self%ndif                = other%ndif
      self%ndifit_h            = other%ndifit_h
      self%ndifit_z            = other%ndifit_z
      self%nintstp             = other%nintstp
      self%nparstp             = other%nparstp
      self%nseqstp             = other%nseqstp
      self%nten_h              = other%nten_h
      self%nten_z              = other%nten_z
      self%nmxlten             = other%nmxlten
      self%nmxlup              = other%nmxlup
      self%mxl_fixlev          = other%mxl_fixlev
      self%nrosten             = other%nrosten
      self%nensten_h           = other%nensten_h
      self%nbctot_h            = other%nbctot_h
      self%nbctot_z            = other%nbctot_z
      self%nintsize3d          = other%nintsize3d
      self%nranvec             = other%nranvec
      self%n_iom_unit          = other%n_iom_unit
      self%nitsha              = other%nitsha
      self%nsamstp             = other%nsamstp
      self%nten_h_clim_hyb     = other%nten_h_clim_hyb
      self%nten_z_clim_hyb     = other%nten_z_clim_hyb
      self%nlstyp              = other%nlstyp
      self%npartot             = other%npartot

      self%cfld                = other%cfld
      self%cgrid               = other%cgrid
      self%cfldname            = other%cfldname
      self%cfile_prefix        = other%cfile_prefix
      self%cdiftyp             = other%cdiftyp

      self%sqtm1ndif           = other%sqtm1ndif
      self%rm1ndif             = other%rm1ndif
      self%sgn_bc              = other%sgn_bc
      self%sgn_lap_hz          = other%sgn_lap_hz
      self%sgn_lap_h           = other%sgn_lap_h
      self%sgn_lap_z           = other%sgn_lap_z
      self%DtoK_hz             = other%DtoK_hz
      self%KtoD_hz             = other%KtoD_hz
      self%DtoK_h              = other%DtoK_h
      self%KtoD_h              = other%KtoD_h
      self%DtoK_z              = other%DtoK_z
      self%KtoD_z              = other%KtoD_z
      self%scagrd_z            = other%scagrd_z
      self%scagrd_z_mld        = other%scagrd_z_mld
      self%scamax_h            = other%scamax_h
      self%scamax_z            = other%scamax_z
      self%scamin_z            = other%scamin_z
      self%scaten_h            = other%scaten_h
      self%scaten_z            = other%scaten_z
      self%rcoast              = other%rcoast
      self%depbnd_mxl          = other%depbnd_mxl
      self%depfix_mxl          = other%depfix_mxl
      self%depmax_mxl          = other%depmax_mxl
      self%wgtten_h            = other%wgtten_h
      self%wgtten_z            = other%wgtten_z
      self%fraten_h            = other%fraten_h
      self%rcoast_msk          = other%rcoast_msk
      self%grdmin_fac          = other%grdmin_fac
      self%rten_hyb_cutlat     = other%rten_hyb_cutlat
      self%rten_hyb_cutlat_bnd = other%rten_hyb_cutlat_bnd
      self%rten_hyb_cutlat_wgt = other%rten_hyb_cutlat_wgt

      self%scalen_h            = other%scalen_h
      self%scalen_z            = other%scalen_z
      self%scagrd_h            = other%scagrd_h

      self%ltrfgrd             = other%ltrfgrd
      self%lsamgrd             = other%lsamgrd
      self%lcoast_bnd          = other%lcoast_bnd
      self%lmxl_bnd            = other%lmxl_bnd
      self%ldep_bnd            = other%ldep_bnd
      self%lcoast              = other%lcoast
      self%lranten             = other%lranten
      self%lens_ten_h          = other%lens_ten_h
      self%lens_ten_z          = other%lens_ten_z
      self%ldif_ten_rea        = other%ldif_ten_rea
      self%ldif_ten_per_rea    = other%ldif_ten_per_rea
      self%ldif_h_ran_wri      = other%ldif_h_ran_wri
      self%ldif_h_ran_rea      = other%ldif_h_ran_rea
      self%ldif_noh            = other%ldif_noh
      self%lnorran_wr4         = other%lnorran_wr4
      self%lnor_nozxnoh        = other%lnor_nozxnoh
      self%ldif_ten_h_clim_hyb = other%ldif_ten_h_clim_hyb
      self%ldif_ten_z_clim_hyb = other%ldif_ten_z_clim_hyb
      self%lcofmat             = other%lcofmat
      self%lcofsymmat          = other%lcofsymmat

      self%fsol                = other%fsol

      !!* Allocatable arrays

      IF ( ALLOCATED( other%cbccomb         ) ) self%cbccomb         = other%cbccomb
      IF ( ALLOCATED( other%soltol_k        ) ) self%soltol_k        = other%soltol_k
      IF ( ALLOCATED( other%difchokk        ) ) self%difchokk        = other%difchokk
      IF ( ALLOCATED( other%cho1            ) ) self%cho1            = other%cho1
      IF ( ALLOCATED( other%cho2            ) ) self%cho2            = other%cho2
      IF ( ALLOCATED( other%cho3            ) ) self%cho3            = other%cho3
      IF ( ALLOCATED( other%K11             ) ) self%K11             = other%K11
      IF ( ALLOCATED( other%K12             ) ) self%K12             = other%K12
      IF ( ALLOCATED( other%K22             ) ) self%K22             = other%K22
      IF ( ALLOCATED( other%K33             ) ) self%K33             = other%K33
      IF ( ALLOCATED( other%K11_fix         ) ) self%K11_fix         = other%K11_fix
      IF ( ALLOCATED( other%K22_fix         ) ) self%K22_fix         = other%K22_fix
      IF ( ALLOCATED( other%K11_min         ) ) self%K11_min         = other%K11_min
      IF ( ALLOCATED( other%K22_min         ) ) self%K22_min         = other%K22_min
      IF ( ALLOCATED( other%K33_min         ) ) self%K33_min         = other%K33_min
      IF ( ALLOCATED( other%K11_col_K22_min ) ) self%K11_col_K22_min = other%K11_col_K22_min
      IF ( ALLOCATED( other%K22_col_K11_min ) ) self%K22_col_K11_min = other%K22_col_K11_min
      IF ( ALLOCATED( other%detsqt_Kh_col   ) ) self%detsqt_Kh_col   = other%detsqt_Kh_col
      IF ( ALLOCATED( other%KP1             ) ) self%KP1             = other%KP1
      IF ( ALLOCATED( other%KP2             ) ) self%KP2             = other%KP2
      IF ( ALLOCATED( other%L11             ) ) self%L11             = other%L11
      IF ( ALLOCATED( other%L22             ) ) self%L22             = other%L22
      IF ( ALLOCATED( other%L33             ) ) self%L33             = other%L33
      IF ( ALLOCATED( other%wrkupd1         ) ) self%wrkupd1         = other%wrkupd1
      IF ( ALLOCATED( other%wrkupd2         ) ) self%wrkupd2         = other%wrkupd2
      IF ( ALLOCATED( other%cofmatii        ) ) self%cofmatii        = other%cofmatii
      IF ( ALLOCATED( other%cofmatjj        ) ) self%cofmatjj        = other%cofmatjj
      IF ( ALLOCATED( other%cofmatkk        ) ) self%cofmatkk        = other%cofmatkk
      IF ( ALLOCATED( other%cofmatii_nd     ) ) self%cofmatii_nd     = other%cofmatii_nd
      IF ( ALLOCATED( other%cofmatij_nd     ) ) self%cofmatij_nd     = other%cofmatij_nd
      IF ( ALLOCATED( other%cofmatjj_nd     ) ) self%cofmatjj_nd     = other%cofmatjj_nd
      IF ( ALLOCATED( other%cofsym_ijk      ) ) self%cofsym_ijk      = other%cofsym_ijk
      IF ( ALLOCATED( other%cofsym_imjk     ) ) self%cofsym_imjk     = other%cofsym_imjk
      IF ( ALLOCATED( other%cofsym_ipjk     ) ) self%cofsym_ipjk     = other%cofsym_ipjk
      IF ( ALLOCATED( other%cofsym_ijmk     ) ) self%cofsym_ijmk     = other%cofsym_ijmk
      IF ( ALLOCATED( other%cofsym_ijpk     ) ) self%cofsym_ijpk     = other%cofsym_ijpk
      IF ( ALLOCATED( other%cofsym_ijkm     ) ) self%cofsym_ijkm     = other%cofsym_ijkm
      IF ( ALLOCATED( other%cofsym_ijkp     ) ) self%cofsym_ijkp     = other%cofsym_ijkp
      IF ( ALLOCATED( other%cofsym_imjmk    ) ) self%cofsym_imjmk    = other%cofsym_imjmk
      IF ( ALLOCATED( other%cofsym_imjpk    ) ) self%cofsym_imjpk    = other%cofsym_imjpk
      IF ( ALLOCATED( other%cofsym_ipjmk    ) ) self%cofsym_ipjmk    = other%cofsym_ipjmk
      IF ( ALLOCATED( other%cofsym_ipjpk    ) ) self%cofsym_ipjpk    = other%cofsym_ipjpk
      IF ( ALLOCATED( other%nsolit_k        ) ) self%nsolit_k        = other%nsolit_k
      IF ( ALLOCATED( other%mxlup_mskt      ) ) self%mxlup_mskt      = other%mxlup_mskt

      !!* Pointer members

      self%difmski     => other%difmski
      self%difmskj     => other%difmskj
      self%difmskk     => other%difmskk
      self%difmskint   => other%difmskint
      self%grd_vp1     => other%grd_vp1
      self%grd_sqt_vp1 => other%grd_sqt_vp1
      self%grd_vm1     => other%grd_vm1
      self%grd_sqt_vm1 => other%grd_sqt_vm1
      self%grd_dp1     => other%grd_dp1
      self%grd_sqt_ap1 => other%grd_sqt_ap1
      self%grd_dp1_tot => other%grd_dp1_tot
      self%grd_dm1_tot => other%grd_dm1_tot
      self%grd_vm1_tot => other%grd_vm1_tot
      self%difmsksum   => other%difmsksum
      self%grdmin_ap1  => other%grdmin_ap1
      self%grdmax_ap1  => other%grdmax_ap1
      self%fgrid       => other%fgrid

      !! *

   END SUBROUTINE assign_fdif_3d


   SUBROUTINE final_fdif_3d( &
      & self )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE final_fdif_3d  ***
      !!
      !! ** Purpose : Finalize fdif_3d
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      TYPE(dan_fdif_3d), INTENT(INOUT) :: &
         & self

      !! *

!      self%difmski     => NULL()
!      self%difmskj     => NULL()
!      self%difmskk     => NULL()
!      self%difmskint   => NULL()
!      self%grd_vp1     => NULL()
!      self%grd_sqt_vp1 => NULL()
!      self%grd_vm1     => NULL()
!      self%grd_sqt_vm1 => NULL()
!      self%grd_dp1     => NULL()
!      self%grd_sqt_ap1 => NULL()
!      self%grd_dp1_tot => NULL()
!      self%grd_dm1_tot => NULL()
!      self%grd_vm1_tot => NULL()
!      self%difmsksum   => NULL()
!      self%grdmin_ap1  => NULL()
!      self%grdmax_ap1  => NULL()
!      self%fgrid       => NULL()

      !! *

   END SUBROUTINE final_fdif_3d

   ! ------------------------------------------------------------------------------
   subroutine c_nv_error_cov_3d_setup_jb( &
      & c_key_self ) bind(c,name='nv_error_cov_3d_setup_jb_f90')                !, &
!      & c_conf ) 
      
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_error_cov_3d_setup_jb  ***
      !!
      !! ** Purpose : construct background error covariance object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
!      type(c_ptr),    intent(in)    :: c_conf

      type(fjb_typ),   pointer :: self
      type(fjb_typ)            :: template

!      type(dan_typ),   pointer :: self             ! DJL test
!      type(dan_typ)            :: template        ! DJL test

      ! Local variables
      !type(fckit_configuration) :: f_conf

      ! Interface
      !f_conf = fckit_configuration(c_conf)

      write(6,*) 'started c_nv_error_cov_3d_setup_jb'

      write(6,*) 'DJL jp_var_max ',jp_var_max

      write(6,*) 'c_key_self initial: ',c_key_self

!!      write(6,*) 'template%array%npi npj npk ',template%array(1)%npi, template%array(1)%npj, template%array(1)%npk

      call flush(6)

      ! Allocate background error covariance object

      call bmat_list%add(c_key_self, template)        ! DJL causing/revealing a Segmentation fault
! DJL fbge_typ looks like the problem
! DJL cov_typ -> cor_typ -> dif_typ problem
      self => get_bmat(c_key_self)

      write(6,*) 'c_nv_error_cov_3d_setup 1'
      write(6,*) 'c_key_self get_bmat: ',c_key_self
      call flush(6)

! DJL added >>> should it go elsewhere?
! (pjb -> self)

      CALL self%fbge_vars%set( &
         & piom_ctl, &
         & jp_var_bge_nam, &
         & ld_setnamelist = .TRUE. )

      ! Construct background error covariance
      ! currently in c_nv_errorcov_3d_linearize
      ! there should be an fbge_vars set in self/x_bkg
      write(6,*) 'DJL self%fbge_vars%n2d ',self%fbge_vars%n2d,' self%fbge_vars%n3d ',self%fbge_vars%n3d
      write(6,*) 'DJL self%fbkg_vars%n2d ',self%fbkg_vars%n2d,' self%fbkg_vars%n3d ',self%fbkg_vars%n3d
      
! <<<
      call flush(piom_ctl%numout)

      write(6,*) 'finished c_nv_error_cov_3d_setup'
      call flush(6)

   end subroutine c_nv_error_cov_3d_setup_jb

   ! ------------------------------------------------------------------------------
   

   ! ------------------------------------------------------------------------------
   subroutine c_nv_error_cov_3d_setup( &
      & c_key_self,  &
      & c_key_geom,  &
      & c_key_x_bkg ) bind(c,name='nv_error_cov_3d_setup_f90')                !, &
!      & c_conf ) 
      
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_error_cov_3d_setup  ***
      !!
      !! ** Purpose : construct background error covariance object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in)    :: c_key_geom
      integer(c_int), intent(in)    :: c_key_x_bkg
!      type(c_ptr),    intent(in)    :: c_conf

      type(fjb_typ),   pointer :: self
      type(fjb_typ)            :: template

!      type(dan_typ),   pointer :: self             ! DJL test
!      type(dan_typ)            :: template        ! DJL test

      type(fgrid_typ), pointer :: geom      ! DJL
      type(ctlvec),    pointer :: x_bkg

      ! Local variables
      !type(fckit_configuration) :: f_conf

      ! Interface
      !f_conf = fckit_configuration(c_conf)

      write(6,*) 'started c_nv_error_cov_3d_setup'

      write(6,*) 'DJL jp_var_max ',jp_var_max

      write(6,*) 'c_key_geom ',c_key_geom
      write(6,*) 'c_key_x_bkg ',c_key_x_bkg
      write(6,*) 'c_key_self initial: ',c_key_self

!!      write(6,*) 'template%array%npi npj npk ',template%array(1)%npi, template%array(1)%npj, template%array(1)%npk

      call flush(6)

      ! Allocate background error covariance object

      call bmat_list%add(c_key_self, template)        ! DJL causing/revealing a Segmentation fault
! DJL fbge_typ looks like the problem
! DJL cov_typ -> cor_typ -> dif_typ problem
      self => get_bmat(c_key_self)

      write(6,*) 'c_nv_error_cov_3d_setup 1'
      write(6,*) 'c_key_self get_bmat: ',c_key_self
      call flush(6)

      ! Retrieve objects

      x_bkg => get_field(c_key_x_bkg)
      geom => get_geom(c_key_geom)

      write(6,*) 'c_nv_error_cov_3d_setup 2'

      ! Allocate fjb_typ pointers

      !allocate( self%x_bkg       )

      ! Set background pointer to the routine argument

      self%x_bkg => x_bkg
      pbal_ctl => self%fbal_ctl
      pbal => self%fbal

! DJL added >>> should it go elsewhere?
! (pjb -> self)
!
      CALL self%fbge_vars%set( &
         & piom_ctl, &
         & jp_var_bge_nam, &
         & ld_setnamelist = .TRUE. )

      ! Construct background error covariance
      ! currently in c_nv_errorcov_3d_linearize
      ! there should be an fbge_vars set in self/x_bkg
!      write(6,*) 'DJL self%fbge_vars%n2d ',self%fbge_vars%n2d,' self%fbge_vars%n3d ',self%fbge_vars%n3d
!      write(6,*) 'DJL self%fbkg_vars%n2d ',self%fbkg_vars%n2d,' self%fbkg_vars%n3d ',self%fbkg_vars%n3d
!      
      call bge_init( &
         & pvar_ctl,  &
         & pwesp,     &
         & piom_ctl,  &
         & piom_file, &
         & self,      &
         & geom,      &
         & .TRUE.,    &
         & .FALSE. )                                  ! Don't read Bkg
!! <<<
!      call flush(piom_ctl%numout)

      write(6,*) 'finished c_nv_error_cov_3d_setup'
      call flush(6)

   end subroutine c_nv_error_cov_3d_setup

   ! ------------------------------------------------------------------------------

   subroutine c_nv_error_cov_3d_delete( &
      & c_key_self, &
      & c_key_geom ) bind(c,name='nv_error_cov_3d_delete_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_error_cov_3d_delete  ***
      !!
      !! ** Purpose : delete background error covariance object
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(inout) :: c_key_self
      integer(c_int), intent(in) :: c_key_geom

      type(fjb_typ), pointer :: self
      type(fgrid_typ), pointer :: geom

      ! Retrieve objects

      self => get_bmat(c_key_self)
      geom => get_geom(c_key_geom)
      
      ! Deallocate background-error covariance data type

      call bge_final( &
         & piom_ctl,      &
         & self%fbge_vars, &
         & self%fbge,     &
         & self%fens )

      call bal_final( &
         & piom_ctl,      &
         & geom, &
         & self%fbal_ctl, &
         & self%fbal )

      self%x_bkg => NULL()

      ! Deallocate fjb_typ pointers

      !deallocate( self%x_bkg       )

      call bmat_list%remove(c_key_self)

   end subroutine c_nv_error_cov_3d_delete

   ! ------------------------------------------------------------------------------

   subroutine c_nv_error_cov_3d_linearize( &
      & c_key_self, &
      & c_key_geom, &
      & c_key_x_fg ) bind(c,name='nv_error_cov_3d_linearize_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_error_cov_3d_linearize  ***
      !!
      !! ** Purpose : linearize background error covariance
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_geom
      integer(c_int), intent(in) :: c_key_x_fg

      type(fjb_typ),   pointer :: self
      type(fgrid_typ), pointer :: geom
      type(ctlvec),    pointer :: x_fg

      ! Retrieve objects

      self => get_bmat(c_key_self)
      geom => get_geom(c_key_geom)
      x_fg => get_field(c_key_x_fg)

      ! Construct background error covariance

      self%x_bkg => x_fg

      call bge_init( &
         & pvar_ctl,  &
         & pwesp,     &
         & piom_ctl,  &
         & piom_file, &
         & self,      &
         & geom )

   end subroutine c_nv_error_cov_3d_linearize

   ! ------------------------------------------------------------------------------

   subroutine c_nv_error_cov_3d_mult( &
      & c_key_self,  &
      & c_key_geom, &
      & c_key_dx_in, &
      & c_key_dx_out ) bind(c,name='nv_error_cov_3d_mult_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_error_cov_3d_mult  ***
      !!
      !! ** Purpose : multiply by background error covariance
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_geom
      integer(c_int), intent(in) :: c_key_dx_in
      integer(c_int), intent(in) :: c_key_dx_out

      type(fjb_typ), pointer :: self
      type(fgrid_typ), pointer :: geom
      type(ctlvec), pointer :: dx_in
      type(ctlvec), pointer :: dx_out

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 1'

      self => get_bmat(c_key_self)

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 2'

      geom => get_geom(c_key_geom)
      
      write(6,*) 'DJL c_nv_error_Cov_3d_mult 3'

      dx_in => get_field(c_key_dx_in)

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 4'

      dx_out => get_field(c_key_dx_out)

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 5'

      ! Call B

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 6'

      dx_out =  dx_in

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 7'

      call bge_cov_opt( &
         & piom_ctl,       &
         & self%fbge_vars, &
         & self%fbge_ctl,  &
         & self%fbge,      &
         & self%fens_ctl,  &
         & self%fens,      &
         & geom,           &
         & dx_out )

      write(6,*) 'DJL c_nv_error_Cov_3d_mult 8'

   end subroutine c_nv_error_cov_3d_mult

   ! ------------------------------------------------------------------------------

   subroutine c_nv_error_cov_3d_randomize( &
      & c_key_self, &
      & c_key_geom, &
      & c_key_dx ) bind(c,name='nv_error_cov_3d_randomize_f90')
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE c_nv_error_cov_3d_randomize  ***
      !!
      !! ** Purpose : randomize with B^{1/2}
      !!
      !!----------------------------------------------------------------------
      implicit none

      integer(c_int), intent(in) :: c_key_self
      integer(c_int), intent(in) :: c_key_geom
      integer(c_int), intent(in) :: c_key_dx

      type(fjb_typ), pointer :: self
      type(fgrid_typ), pointer :: geom
      type(ctlvec), pointer :: dx

      type(ctlvec) :: z_v ! Control vector: v

      real(kind=wp) :: z_norm

      ! Retrieve objects

      self => get_bmat(c_key_self)
      geom => get_geom(c_key_geom)
      dx => get_field(c_key_dx)

      ! Allocate a control vector

      write(6,*) 'DJL randomize: allocate control vector v'

      call alloc_ctlvec( &
         & piom_ctl,       &
         & geom%fmpp,      &
         & z_v,            &
         & self%fbge_vars, &
         & 'v' )

      ! Generate uncorrelated Gaussian random noise: v

      call bge_covsqt_ran( &
         & piom_ctl,       &
         & self%fbge_vars, &
         & self%fbge_ctl,  &
         & self%fbge,      &
         & self%fens_ctl,  &
         & self%fens,      &
         & geom%fmpp,      &
         & geom%fmsk,      &
         & 0.0_wp,         &
         & 1.0_wp,         &
         & z_v )

      z_norm = norm( geom%fmpp, z_v)

      if( piom_ctl%lwp ) write(piom_ctl%numout,*)
      if( piom_ctl%lwp ) write(piom_ctl%numout,*) 'Norm of v: ', z_norm
      if( piom_ctl%lwp ) write(piom_ctl%numout,*)

      FLUSH(piom_ctl%numout)

      ! Compute dx_l = B^1/2 v

      call bge_covsqt_opt( &
         & piom_ctl,       &
         & self%fbge_vars, &
         & self%fbge_ctl,  &
         & self%fbge,      &
         & self%fens_ctl,  &
         & self%fens,      &
         & geom,           &
         & dx,             &
         & z_v )

      z_norm = norm( geom%fmpp, dx)

      if( piom_ctl%lwp ) write(piom_ctl%numout,*)
      if( piom_ctl%lwp ) write(piom_ctl%numout,*) 'Norm of B^{1/2}.v: ', z_norm
      if( piom_ctl%lwp ) write(piom_ctl%numout,*)

      FLUSH(piom_ctl%numout)

      ! Deallocate control vector

      call dealloc_ctlvec( &
         & piom_ctl, &
         & z_v )

   end subroutine c_nv_error_cov_3d_randomize

   ! ------------------------------------------------------------------------------

   function get_bmat(c_self)
   integer(kind=c_int), intent(in) :: c_self
   type(fjb_typ), pointer :: get_bmat
   class(*), pointer :: this=>null()

   call bmat_list%get(c_self, this)

   select type(self => this)
   type is (fjb_typ)
      get_bmat=>self
   class default
      call abor1_ftn('get_bmat : unexpected type ')
   end select

   end function get_bmat

   ! ------------------------------------------------------------------------------

end module bmat_mod
