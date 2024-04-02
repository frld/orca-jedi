/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "FieldsNV.h"

#include <iostream>
#include <map>
#include <math.h>
#include <string>

#include "eckit/config/Configuration.h"

//include "nemovar/GomsNV.h"
//include "nemovar/LocationsNV.h"

// -----------------------------------------------------------------------------
namespace nv {
// -----------------------------------------------------------------------------

FieldsNV::FieldsNV(const GeometryNV & geom, const VariablesNV & vars, 
                   const util::DateTime & vt)
    : tstep_(0), reftime_(vt), time_(vt), geom_(new GeometryNV(geom)), 
      vars_(new VariablesNV(vars)), keyFlds_(0), keyFldstl_(0)
{
    oops::Log::trace() << "FieldsNV::FieldsNV begin" << std::endl;    
    oops::Log::trace() << "FieldsNV::FieldsNV geom geomnv_int " << (*geom_).geomnv_int << std::endl;
    oops::Log::trace() << "FieldsNV::FieldsNV vars variablesnv_int " << (*vars_).variablesnv_int << std::endl;

    ASSERT(geom_->toFortran() != 0);
    ASSERT(vars_->toFortran() != 0);
    nv_field_create_f90(keyFlds_, geom_->toFortran(), vars_->toFortran());
    ASSERT(keyFlds_ != 0);
    oops::Log::trace() << "FieldsNV::FieldsNV done" << std::endl;
}

// -----------------------------------------------------------------------------

FieldsNV::FieldsNV(const FieldsNV & other, const bool copyFields)
    : tstep_(other.tstep_), reftime_(other.reftime_), time_(other.time_), 
      geom_(other.geom_), vars_(other.vars_), keyFlds_(0), keyFldstl_(0)
{
    ASSERT(geom_->toFortran() != 0);
    ASSERT(vars_->toFortran() != 0);
    nv_field_create_f90(keyFlds_, geom_->toFortran(), vars_->toFortran());
    ASSERT(keyFlds_ != 0);
    if (copyFields)  {
        ASSERT(other.keyFlds_ != 0);
        nv_field_copy_f90(keyFlds_, other.toFortran());
    }
    oops::Log::trace() << "FieldsNV::FieldsNV copy ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

FieldsNV::~FieldsNV() {
    nv_field_delete_f90(keyFlds_);
    oops::Log::trace() << "FieldsNV::~FieldsNV done" << std::endl;
}

// -----------------------------------------------------------------------------

FieldsNV & FieldsNV::operator=(const FieldsNV & rhs) {
    ASSERT(rhs.keyFlds_ != 0);
    nv_field_copy_f90(keyFlds_, rhs.toFortran());
    ASSERT(keyFlds_ != 0);
    oops::Log::trace() << "FieldsNV::operator= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

FieldsNV & FieldsNV::operator+=(const FieldsNV & rhs) {
    nv_field_self_add_f90(keyFlds_, rhs.toFortran());
    oops::Log::trace() << "FieldsNV::operator+= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

FieldsNV & FieldsNV::operator-=(const FieldsNV & rhs) {
    nv_field_self_sub_f90(keyFlds_, rhs.toFortran());
    oops::Log::trace() << "FieldsNV::operator+= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

FieldsNV & FieldsNV::operator*=(const double & zz) {
    nv_field_self_mul_f90(keyFlds_, zz);
    oops::Log::trace() << "FieldsNV::operator*= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

void FieldsNV::zero() {
    oops::Log::trace() << "FieldsNV::zero start" << std::endl;
    nv_field_zero_f90(keyFlds_);
    oops::Log::trace() << "FieldsNV::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::zeroTlFields() {
    ASSERT(keyFldstl_ != 0);
    nv_fieldtlad_zero_f90(keyFldstl_);
    oops::Log::trace() << "FieldsNV::zeroTlFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::zero(const util::DateTime & time) {
    oops::Log::trace() << "FieldsNV::zero time start" << std::endl;
    time_ = time;
    oops::Log::trace() << "FieldsNV::zero time 1" << std::endl;
    nv_field_zero_f90(keyFlds_);
    
// c++ set field to 1?
    // F90flds
    // pdata    
    
    oops::Log::trace() << "FieldsNV::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::axpy(const double & zz, const FieldsNV & rhs) {
    nv_field_axpy_f90(keyFlds_, rhs.toFortran(), zz);
    oops::Log::trace() << "FieldsNV::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

double FieldsNV::dot_product_with(const FieldsNV & other) const {
    double zz = 0.0;
    nv_field_dot_prod_f90(keyFlds_, other.toFortran(), geom_->toFortran(), zz);
    oops::Log::trace() << "FieldsNV::dot_product_with done" << std::endl;
    return zz;
}

// -----------------------------------------------------------------------------

void FieldsNV::schur_product_with(const FieldsNV & other) {
    //nv_field_schur_prod_f90(keyFlds_, other.toFortran());
    oops::Log::trace() << "FieldsNV::shur_product_with not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

//void FieldsNV::interpolate(const LocationsNV & locs, GomsNV & gom) const {
//    //nv_field_interp_f90(keyFlds_, locs.toFortran(), gom.toFortran());
//    oops::Log::trace() << "FieldsNV::interpolate done" << std::endl;
//}

// -----------------------------------------------------------------------------

//void FieldsNV::interpolateTL(const LocationsNV & locs, GomsNV & gom) const {
//    util::Duration elapsed_time(time_ - reftime_);
//    const int current_iteration = elapsed_time.toSeconds()/tstep_.toSeconds(); 
//    oops::Log::debug() << "FieldsNV::interpolateTL current_iteration = "
//                       << current_iteration << std::endl;
//    ASSERT(keyFldstl_ != 0);
//    nv_field_interp_tl_f90(keyFldstl_, geom_->toFortran(), locs.toFortran(), 
//                           gom.toFortran(), current_iteration);
//    oops::Log::trace() << "FieldsNV::interpolateTL done" << std::endl;
//}

// -----------------------------------------------------------------------------

//void FieldsNV::interpolateAD(const LocationsNV & locs, const GomsNV & gom) {
//    util::Duration elapsed_time(time_ - reftime_);
//    const int current_iteration = elapsed_time.toSeconds()/tstep_.toSeconds(); 
//    oops::Log::debug() << "FieldsNV::interpolateAD current_iteration = " 
//                       << current_iteration << std::endl;
//    ASSERT(keyFldstl_ != 0);
//    nv_field_interp_ad_f90(keyFldstl_, geom_->toFortran(), locs.toFortran(), 
//                           gom.toFortran(), current_iteration);
//    oops::Log::trace() << "FieldsNV::interpolateAD done" << std::endl;
//}

// -----------------------------------------------------------------------------

void FieldsNV::toTlFields(const int & skipBalance) {
    nv_field_totlfields_f90(keyFlds_, keyFldstl_, geom_->toFortran(), 
                            vars_->toFortran(), skipBalance);
    ASSERT(keyFldstl_ != 0);
    oops::Log::trace() << "FieldsNV::toTlFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::toIncrFields() {
    ASSERT(keyFldstl_ != 0);
    nv_field_toincrfields_f90(keyFlds_, keyFldstl_, geom_->toFortran(), vars_->toFortran()); 
    oops::Log::trace() << "FieldsNV::toIncrFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::changeResolution(const FieldsNV & rhs) {
    //nv_field_change_resol_f90(keyFlds_, rhs.toFortran());
    oops::Log::trace() << "FieldsNV::changeResolution not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::add(const FieldsNV & rhs) {
    nv_field_add_incr_f90(keyFlds_, rhs.toFortran());
    oops::Log::trace() << "FieldsNV::add done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::diff(const FieldsNV & x1, const FieldsNV & x2) {
    nv_field_diff_incr_f90(keyFlds_, x1.toFortran(), x2.toFortran());
    oops::Log::trace() << "FieldsNV::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

void FieldsNV::read(const eckit::Configuration & config) {
    const eckit::Configuration * conf = &config;
    nv_field_read_file_f90(keyFlds_, geom_->toFortran(), &conf);
    oops::Log::trace() << "FieldsNV::read done" << std::endl;
}

// -----------------------------------------------------------------------------

// serialize (input field - output to array)
void FieldsNV::getarray( int64_t & len, std::vector<double> & array) {
//std::vector<double> & array) {
    std::vector<double> v_inc(len, 0);
    oops::Log::trace() << "FieldsNV::getarray start" << std::endl;
    nv_field_get_f90(keyFlds_, len, v_inc.data());
    array.insert(array.end(), v_inc.begin(), v_inc.end());
    oops::Log::trace() << "FieldsNV::getarray done" << std::endl;
}

// -----------------------------------------------------------------------------

// deserialize (input array - output to field)
void FieldsNV::setarray( int64_t & len, const std::vector<double> & array) {
//std::vector<double> & array) {
    oops::Log::trace() << "FieldsNV::setarray start" << std::endl;
    nv_field_set_f90(keyFlds_, len, array.data());
    oops::Log::trace() << "FieldsNV::setarray done" << std::endl;
}


// -----------------------------------------------------------------------------

void FieldsNV::write(const eckit::Configuration & config) const {
    const eckit::Configuration * conf = &config;
    nv_field_write_file_f90(keyFlds_, geom_->toFortran(), &conf);
    oops::Log::trace() << "FieldsNV::write done" << std::endl;
}

// -----------------------------------------------------------------------------

double FieldsNV::norm() const {
    double zz = 0.0;
    nv_field_norm_f90(keyFlds_, geom_->toFortran(), zz);
    oops::Log::trace() << "FieldsNV::norm done" << std::endl;
    return zz;
}

// -----------------------------------------------------------------------------

void FieldsNV::print(std::ostream & os) const {
    os << "FieldsNV::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace nv
