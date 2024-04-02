/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "ErrorCovarianceNV.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "nemovar/FieldsNV.h"
#include "nemovar/FortranNV.h"
#include "nemovar/GeometryNV.h"
#include "nemovar/IncrementNV.h"
//include "nemovar/IncrModCtlVecNV.h"
#include "nemovar/StateNV.h"
//include "nemovar/VariablesNV.h"

// -----------------------------------------------------------------------------
namespace nv {
// -----------------------------------------------------------------------------

ErrorCovarianceNV::ErrorCovarianceNV(const GeometryNV & resol, const VariablesNV & cvars,
                                         const eckit::Configuration & conf, const StateNV & bg)
//    : geom_(), ctrlvars_(cvars), time_(util::DateTime(conf.getString("date"))), keyBmat_(0)
    : geom_(), ctrlvars_(cvars), time_(util::DateTime("2021-06-30T12:00:00Z")), keyBmat_(0)                 // DJL hardwire date
{
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 1a " << std::endl;
    const eckit::Configuration * configc = &conf;
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 1a2 setting up geometry needed for other methods" << std::endl;
    geom_.reset(new GeometryNV(resol));
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 1b " << std::endl;
    nv_error_cov_3d_setup_f90(keyBmat_, resol.toFortran(), bg.toFortran());    // , &configc);   // DJL
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 1c " << std::endl;
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV keyBmat_ " << keyBmat_ << std::endl;
//    ASSERT(keyBmat_ != 0);  ! DJL
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV created" << std::endl;
}

// -----

ErrorCovarianceNV::ErrorCovarianceNV(const GeometryNV & resol, const VariablesNV & cvars,
                                         const eckit::Configuration & conf)
//    : geom_(), ctrlvars_(cvars), time_(util::DateTime(conf.getString("date"))), keyBmat_(0)
    : geom_(), ctrlvars_(cvars), time_(util::DateTime("2021-06-30T12:00:00Z")), keyBmat_(0)                // DJL hardwire date
{
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 2a " << std::endl;
    StateNV bg(resol, cvars, conf);   // create empty state
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 2a2 setting up geometry needed for other methods" << std::endl;
    geom_.reset(new GeometryNV(resol));
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 2b " << std::endl;
    const eckit::Configuration * configc = &conf;
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 2c " << std::endl;
    nv_error_cov_3d_setup_f90(keyBmat_, resol.toFortran(), bg.toFortran()); // , &configc);
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 2d " << std::endl;
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV keyBmat_ " << keyBmat_ << std::endl;
    ASSERT(keyBmat_ != 0);
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV created" << std::endl;
}

// -----

ErrorCovarianceNV::ErrorCovarianceNV(const VariablesNV & cvars,
                                         const eckit::Configuration & conf)
    : ctrlvars_(cvars), keyBmat_(0)
{
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 3a " << std::endl;
    nv_error_cov_3d_setup_jb_f90(keyBmat_);
    oops::Log::trace() << "ErrorCovarianceNV::ErrorCovarianceNV 3 end " << std::endl;
}


// -----------------------------------------------------------------------------

ErrorCovarianceNV::~ErrorCovarianceNV() {
//    nv_error_cov_3d_delete_f90(keyBmat_, geom_->toFortran());   // DJL
    oops::Log::trace() << "ErrorCovarianceNV::~ErrorCovarianceNV deleted" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovarianceNV::linearize(const StateNV & fg, const GeometryNV & resol) {
    geom_.reset(new GeometryNV(resol));
    nv_error_cov_3d_linearize_f90(keyBmat_, geom_->toFortran(), fg.toFortran());
    oops::Log::trace() << "ErrorCovarianceNV::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovarianceNV::multiply(const IncrementNV & dxin, IncrementNV & dxout) const {
    oops::Log::trace() << "ErrorCovarianceNV::multiply start" << std::endl;
    
    oops::Log::trace() << "ErrorCovarianceNV::multiply 1" << std::endl;
    oops::Log::trace() << "ErrorCovarianceNV::multiply keyBmat_ " << keyBmat_ << std::endl;
    oops::Log::trace() << "ErrorCovarianceNV::multiply geom_->geomnv_int " << geom_->geomnv_int << std::endl;
//    oops::Log::trace() << "ErrorCovarianceNV::multiply keys " << geom_->toFortran() << std::endl; 
    nv_error_cov_3d_mult_f90(keyBmat_, geom_->toFortran(), dxin.toFortran(), dxout.toFortran());   // DJL
    
    oops::Log::trace() << "ErrorCovarianceNV::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovarianceNV::inverseMultiply(const IncrementNV & dxin, IncrementNV & dxout) const {
    oops::Log::trace() << "ErrorCovarianceNV::inverseMultiply not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

//void ErrorCovarianceNV::multiplySqrt(const IncrModCtlVecNV & dv,
//                                     IncrementNV & dx) const {
//    oops::Log::trace() << "ErrorCovarianceNV multiplySqrt not implemented" << std::endl;
//}

// -----------------------------------------------------------------------------

//void ErrorCovarianceNV::multiplySqrtTrans(const IncrementNV & dx,
//                                          IncrModCtlVecNV & dv) const {
//    oops::Log::trace() << "ErrorCovarianceNV multiplySqrtTrans not implemented" << std::endl;
//}


void ErrorCovarianceNV::randomize(IncrementNV & dx) const {
    nv_error_cov_3d_randomize_f90(keyBmat_, geom_->toFortran(), dx.toFortran());
    oops::Log::trace() << "ErrorCovarianceNV::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------
void ErrorCovarianceNV::print(std::ostream & os) const {
    os << "ErrorCovarianceNV::print not implemented yet." << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace nv
