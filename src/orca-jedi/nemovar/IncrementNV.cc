/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "IncrementNV.h"

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "nemovar/FieldsNV.h"
#include "nemovar/GeometryNV.h"
//include "nemovar/GomsNV.h"
//include "nemovar/LocationsNV.h"
#include "nemovar/StateNV.h"
#include "nemovar/VariablesNV.h"


namespace nv {

// -----------------------------------------------------------------------------
// Constructor, destructor
// -----------------------------------------------------------------------------

IncrementNV::IncrementNV(const GeometryNV & resol, const VariablesNV & vars,
                         const util::DateTime & vt)
    : fields_(resol, vars, vt)
{
    fields_.zero();
    oops::Log::trace() << "IncrementNV::IncrementNV ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

IncrementNV::IncrementNV(const GeometryNV & resol, const IncrementNV & other)
    : fields_(resol, other.fields_.getVariables(), other.refTime())
{
    fields_.validTime() = other.validTime();
    fields_.toFortranTl() = other.fields_.toFortranTl();
    fields_.changeResolution(other.fields_);
    oops::Log::trace() << "IncrementNV::IncrementNV created by interpolation." << std::endl;
}

// -----------------------------------------------------------------------------

IncrementNV::IncrementNV(const IncrementNV & other, const bool copy)
    : fields_(other.fields_, copy)
{
    fields_.toFortranTl() = other.fields_.toFortranTl();
    oops::Log::trace() << "IncrementNV::IncrementNV copy-created." << std::endl;
}

// -----------------------------------------------------------------------------
IncrementNV::~IncrementNV() {
    oops::Log::trace() << "IncrementNV::IncrementNV dtor done" << std::endl;
}
// -----------------------------------------------------------------------------
// Basic operators
// -----------------------------------------------------------------------------

void IncrementNV::diff(const StateNV & x1, const StateNV & x2) {
    ASSERT(this->validTime() == x1.validTime());
    ASSERT(this->validTime() == x2.validTime());
    fields_.diff(x1.getFields(), x2.getFields());
    oops::Log::trace() << "IncrementNV::diff done" << std::endl;
}

// -----------------------------------------------------------------------------

IncrementNV & IncrementNV::operator=(const IncrementNV & rhs) {
    fields_ = rhs.fields_;
    oops::Log::trace() << "IncrementNV::operator= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

IncrementNV & IncrementNV::operator+=(const IncrementNV & dx) {
    ASSERT(this->validTime() == dx.validTime());
    fields_ += dx.fields_;
    oops::Log::trace() << "IncrementNV::operator+= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

IncrementNV & IncrementNV::operator-=(const IncrementNV & dx) {
    ASSERT(this->validTime() == dx.validTime());
    fields_ -= dx.fields_;
    oops::Log::trace() << "IncrementNV::operator-= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

IncrementNV & IncrementNV::operator*=(const double & zz) {
    fields_ *= zz;
    oops::Log::trace() << "IncrementNV::operator*= done" << std::endl;
    return *this;
}

// -----------------------------------------------------------------------------

void IncrementNV::zero() {
    fields_.zero();
    oops::Log::trace() << "IncrementNV::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void IncrementNV::zeroTlFields() {
    fields_.zeroTlFields();
    oops::Log::trace() << "IncrementNV::zeroTlFields done" << std::endl;
}

// -----------------------------------------------------------------------------

void IncrementNV::zero(const util::DateTime & time) {
    fields_.zero(time);
    oops::Log::trace() << "IncrementNV::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

void IncrementNV::axpy(const double & zz, const IncrementNV & dx,
                        const bool check) {
    ASSERT(!check || this->validTime() == dx.validTime());
    fields_.axpy(zz, dx.fields_);
    oops::Log::trace() << "IncrementNV::axpy done" << std::endl;
}

// -----------------------------------------------------------------------------

double IncrementNV::dot_product_with(const IncrementNV & other) const {
    return dot_product(fields_, other.fields_);
    oops::Log::trace() << "IncrementNV::dot_product done" << std::endl;

}

// -----------------------------------------------------------------------------

void IncrementNV::schur_product_with(const IncrementNV & other) {
    fields_.schur_product_with(other.fields_);
    oops::Log::trace() << "IncrementNV::schur_product_with done" << std::endl;
}

// -----------------------------------------------------------------------------

void IncrementNV::accumul(const double & zz, const StateNV & xx) {
    //fields_.axpy(zz, xx.getFields());
    oops::Log::trace() << "IncrementNV::accumul not implemented" << std::endl;
}

// -----------------------------------------------------------------------------
/// Interpolate to observation location
// -----------------------------------------------------------------------------

//void IncrementNV::interpolateTL(const LocationsNV & locs, GomsNV & cols) const {
//    fields_.interpolateTL(locs, cols);
//    oops::Log::trace() << "IncrementNV::interpolateTL done" << std::endl;
//}

// -----------------------------------------------------------------------------

//void IncrementNV::interpolateAD(const LocationsNV & locs, const GomsNV & cols) {
//    fields_.interpolateAD(locs, cols);
//    oops::Log::trace() << "IncrementNV::interpolateAD done" << std::endl;
//}

// -----------------------------------------------------------------------------
// I/O and diagnostics
// -----------------------------------------------------------------------------

void IncrementNV::read(const eckit::Configuration & files) {
    fields_.read(files);
    oops::Log::trace() << "IncrementNV::read done" << std::endl;
}

// -----------------------------------------------------------------------------

void IncrementNV::write(const eckit::Configuration & files) const {
    fields_.write(files);
    oops::Log::trace() << "IncrementNV::write done" << std::endl;
}

// -----------------------------------------------------------------------------

void IncrementNV::print(std::ostream & os) const {
    os << "Valid time: " << this->validTime() << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace nv
