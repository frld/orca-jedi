#ifndef NV_FORTRANNV_H
#define NV_FORTRANNV_H

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace nv {

// Geometry key type
typedef int F90geom;
// Config key type
typedef int F90conf;
// Variables key type
typedef int F90vars;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// B cov key type
typedef int F90bmat;
// Fields key type
typedef int F90flds;
// Fields TL/AD key type
typedef int F90fldstl;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
//typedef int F90bmat;
// Observation vector key type
typedef int F90ovec;
// Allobs key type
typedef int F90obop;
// Observation data base type
typedef int F90odb;
// R matrix key type
typedef int F90rmat;
// Obs bias correction key type
typedef int F90obias;
// Obs bias increment key type
typedef int F90obiasinc;
// Mod clt vec key type
typedef int F90modctlvec;
// Ens clt vec key type
typedef int F90ensctlvec;

/// Interface to Fortran NEMO/NEMOVAR model
/*!
  * The core of the NEMO/NEMOVAR model is coded in Fortran.
  * Here we define the interfaces to the Fortran code.
  */

extern "C" {
// ---------------------------------------------------------------------------
//  Geometry
// ---------------------------------------------------------------------------
void nv_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
void nv_geo_delete_f90(F90geom &);
// ---------------------------------------------------------------------------
//  Model
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//  Fields
// ---------------------------------------------------------------------------
void nv_field_create_f90(F90flds &, const F90geom &, const F90vars &);
void nv_field_delete_f90(F90flds &);
void nv_field_set_f90(const F90flds &, int64_t &, const double[]);
void nv_field_get_f90(const F90flds &, int64_t &, double[]);
void nv_field_zero_f90(const F90flds &);
void nv_field_copy_f90(const F90flds &, const F90flds &);
void nv_field_self_add_f90(const F90flds &, const F90flds &);
void nv_field_self_sub_f90(const F90flds &, const F90flds &);
void nv_field_self_mul_f90(const F90flds &, const double &);
void nv_field_axpy_f90(const F90flds &, const F90flds &, const double &);
void nv_field_dot_prod_f90(const F90flds &, const F90flds &, const F90geom &, double &);
void nv_field_norm_f90(const F90flds &, const F90geom &, double &);
//void nv_field_interp_tl_f90(const F90fldstl &, const F90geom &,
//                            const F90locs &, const F90goms &,
//                            const int &);
//void nv_field_interp_ad_f90(const F90fldstl &, const F90geom &,
//                            const F90locs &, const F90goms &,
//                            const int &);
void nv_field_totlfields_f90(const F90flds &, F90fldstl &,
                             const F90geom &, const F90vars &,
                             const int &);
void nv_field_toincrfields_f90(const F90flds &, F90fldstl &,
                               const F90geom &, const F90vars &);
void nv_field_add_incr_f90(const F90flds &, const F90flds &);
void nv_field_diff_incr_f90(const F90flds &, const F90flds &,
                            const F90flds &);
void nv_field_read_file_f90(const F90flds &, const F90geom &,
                            const eckit::Configuration * const *);
void nv_field_write_file_f90(const F90flds &, const F90geom &,
                             const eckit::Configuration * const *);

// ---------------------------------------------------------------------------
//  Fields TL/AD
// ---------------------------------------------------------------------------
void nv_fieldtlad_create_f90(F90fldstl &, const F90geom &, const F90vars &);
void nv_fieldtlad_delete_f90(F90fldstl &);
void nv_fieldtlad_zero_f90(const F90fldstl &);

// ---------------------------------------------------------------------------
//  Background error
// ---------------------------------------------------------------------------
void nv_error_cov_3d_setup_f90(F90bmat &, const F90geom &, const F90flds &);  // ,
//                               const eckit::Configuration * const *);
void nv_error_cov_3d_setup_jb_f90(F90bmat &);  // ,
//                               const eckit::Configuration * const *);
void nv_error_cov_3d_delete_f90(F90bmat &, const F90geom &);
void nv_error_cov_3d_linearize_f90(const F90bmat &, const F90geom &,
                                   const F90flds &);
void nv_error_cov_3d_mult_f90(const F90bmat &, const F90geom &, const F90flds &,
                              const F90flds &);
void nv_error_cov_3d_randomize_f90(const F90bmat &, const F90geom &, const F90flds &);

// -----------------------------------------------------------------------------
//  Localization matrix
// -----------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//  Variables
// ---------------------------------------------------------------------------
void nv_variables_create_f90(F90vars &, const eckit::Configuration * const *);
void nv_variables_setup_f90(F90vars &, const int &);
void nv_variables_clone_f90(F90vars &, const F90vars &);
void nv_variables_delete_f90(F90vars &);
// ---------------------------------------------------------------------------
//  Locations
// ---------------------------------------------------------------------------
void nv_obs_locations_f90(F90locs &, const F90odb &, const int &, const int &);
void nv_obs_locations_delete_f90(F90locs &);

// ---------------------------------------------------------------------------
//  Local Values (GOM)
// ---------------------------------------------------------------------------
void nv_gom_create_f90(F90goms &, const F90odb &,
                       const F90vars &, const F90geom &);
void nv_gom_delete_f90(F90goms &);
void nv_gom_zero_f90(F90goms &);

// ---------------------------------------------------------------------------
//  Observation Handler
// ---------------------------------------------------------------------------
  void nv_odb_setup_f90(F90odb &, const F90geom &);
  void nv_odb_delete_f90(F90odb &);
  void nv_print_jo_table_f90(const F90odb &, const F90geom &,
                             const F90ovec &, const F90ovec &);

// ---------------------------------------------------------------------------
//  Observations Operators
// ---------------------------------------------------------------------------
  void nv_allobs_equiv_f90(const F90odb &, const F90goms &, F90ovec &);

// ---------------------------------------------------------------------------
//  Observations Operators TL/AD
// ---------------------------------------------------------------------------
  void nv_allobs_equiv_tl_f90(const F90odb &, const F90goms &, F90ovec &);
  void nv_allobs_equiv_ad_f90(const F90odb &, F90goms &, const F90ovec &);

// ---------------------------------------------------------------------------
//  Observation Vectors
// ---------------------------------------------------------------------------
  void nv_obsvec_setup_f90(F90ovec &, const F90odb &);
  void nv_obsvec_clone_f90(F90ovec &, const F90ovec &);
  void nv_obsvec_delete_f90(F90ovec &);
  void nv_obsvec_assign_f90(const F90ovec &, const F90ovec &);
  void nv_obsvec_mul_scal_f90(const F90ovec &, const double &);
  void nv_obsvec_add_f90(const F90ovec &, const F90ovec &);
  void nv_obsvec_sub_f90(const F90ovec &, const F90ovec &);
  void nv_obsvec_mul_f90(const F90ovec &, const F90ovec &);
  void nv_obsvec_div_f90(const F90ovec &, const F90ovec &);
  void nv_obsvec_zero_f90(const F90ovec &);
  void nv_obsvec_axpy_f90(const F90ovec &, const F90ovec &, const double &);
  void nv_obsvec_invert_f90(const F90ovec &);
  void nv_obsvec_random_f90(const F90ovec &);
  void nv_obsvec_dotprod_f90(const F90ovec &, const F90ovec &, double &);
  void nv_obsvec_nobs_f90(const F90ovec &, int &);
  void nv_obsvec_minmaxavg_f90(const F90ovec &, double &, double &, double &);
  void nv_obsvec_read_f90(const F90odb &, const F90ovec &,
                          const char *, const int &);
  void nv_obsvec_write_f90(const F90odb &, const F90ovec &,
                           const char *, const int &);
// ---------------------------------------------------------------------------
//  Observation Errors
// ---------------------------------------------------------------------------
  void nv_r_setup_f90(F90rmat &, const F90odb &, const F90geom &,
                      const eckit::Configuration * const *);
  void nv_r_delete_f90(F90rmat &);
  void nv_r_mult_f90(const F90rmat &, const F90ovec &, const F90ovec &);
  void nv_r_imult_f90(const F90rmat &, const F90ovec &, const F90ovec &);

// ---------------------------------------------------------------------------
//  System and other libraries
// ---------------------------------------------------------------------------
  void nv_init_f90(const eckit::Configuration * const *);
  void nv_finalize_f90();

}

}  // namespace nv

#endif // NV_FORTRANNV_H
