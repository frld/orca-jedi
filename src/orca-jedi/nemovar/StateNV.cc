/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "StateNV.h"

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

#include "nemovar/FieldsNV.h"
#include "nemovar/GeometryNV.h"
//include "nemovar/GomsNV.h"
#include "nemovar/IncrementNV.h"
//include "nemovar/LocationsNV.h"
#include "nemovar/ModelNV.h"
#include "nemovar/VariablesNV.h"

namespace nv {

// -----------------------------------------------------------------------------
// Constructor, destructor
// -----------------------------------------------------------------------------

StateNV::StateNV(const GeometryNV & resol, const ModelNV &, const eckit::Configuration & file)
{
    VariablesNV vars(eckit::LocalConfiguration(file,"Variables"));
    
    oops::Log::trace() << "StateNV::StateNV trying to access date DJL" << std::endl;
    util::DateTime time(file.getString("date"));
    oops::Log::info() << "StateNV::StateNV time: " << time << std::endl;
    fields_.reset(new FieldsNV(resol, vars, time));
    ASSERT(fields_->toFortran() != 0);
    ASSERT(fields_->getVariables().toFortran() != 0);
    fields_->read(file);

    oops::Log::trace() << "StateNV::StateNV created and read in." << std::endl;
}

// -----------------------------------------------------------------------------

StateNV::StateNV(const GeometryNV & resol, const StateNV & other)
{
    fields_.reset(new FieldsNV(resol, other.fields_->getVariables(), other.refTime()));
    fields_->validTime() = other.validTime();
    ASSERT(fields_->toFortran() != 0);
    ASSERT(fields_->getVariables().toFortran() != 0);
    fields_->changeResolution(*other.fields_);
    oops::Log::trace() << "StateNV::StateNV created by interpolation." << std::endl;
}

// -----------------------------------------------------------------------------

StateNV::StateNV(const StateNV & other)
{
    fields_.reset(new FieldsNV(*other.fields_));
    ASSERT(fields_->toFortran() != 0);
    ASSERT(fields_->getVariables().toFortran() != 0);
    oops::Log::trace() << "StateNV::StateNV copy-created." << std::endl;
}

// ----
// zero

StateNV::StateNV(const GeometryNV & resol, const VariablesNV & cvars, const eckit::Configuration & conf)
{
    //VariablesNV vars(eckit::LocalConfiguration(conf,"Variables"));
    oops::Log::trace() << "StateNV zero trying to access date from conf DJL" << std::endl;
    //util::DateTime time(conf.getString("date"));
    util::DateTime time("2021-06-30T12:00:00Z");                 // DJL hardwire
    oops::Log::info() << "StateNV::StateNV time: " << time << std::endl;
    fields_.reset(new FieldsNV(resol, cvars, time));
    oops::Log::info() << "StateNV::StateNV reset finished" << std::endl;
    ASSERT(fields_->toFortran() != 0);
    ASSERT(fields_->getVariables().toFortran() != 0);
    fields_->zero(time);

    oops::Log::trace() << "StateNV::StateNV created zeroed." << std::endl;
}

// -----------------------------------------------------------------------------

StateNV::~StateNV() {
    oops::Log::trace() << "StateNV::StateNV destructed" << std::endl;
}

// -----------------------------------------------------------------------------
// Basic operators
// -----------------------------------------------------------------------------

StateNV & StateNV::operator=(const StateNV & rhs) {
    *fields_ = *rhs.fields_;
    return *this;
}

// -----------------------------------------------------------------------------
// Interpolate to observation location
// -----------------------------------------------------------------------------

//void StateNV::interpolate(const LocationsNV & locs, GomsNV & cols) const {
//    fields_->interpolate(locs, cols);
//}

// -----------------------------------------------------------------------------
// Interactions with Increments
// -----------------------------------------------------------------------------

StateNV & StateNV::operator+=(const IncrementNV & dx) {
    ASSERT(this->validTime() == dx.validTime());
    fields_->add(dx.getFields());
    return *this;
}

// -----------------------------------------------------------------------------
// I/O and diagnostics
// -----------------------------------------------------------------------------

void StateNV::read(const eckit::Configuration & files) {
    fields_->read(files);
}

// -----------------------------------------------------------------------------

void StateNV::write(const eckit::Configuration & files) const {
    fields_->write(files);
}

// -----------------------------------------------------------------------------

void StateNV::print(std::ostream & os) const {
    os << "StateNV Valid time: " << this->validTime() << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace nv
