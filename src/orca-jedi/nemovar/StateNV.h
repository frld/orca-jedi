/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_STATENV_H
#define NV_STATENV_H

#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "nemovar/FieldsNV.h"
#include "nemovar/FortranNV.h"

namespace eckit {
  class Configuration;
}

namespace nv {
    class GeometryNV;
//    class GomsNV;
    class IncrementNV;
//    class LocationsNV;
    class ModelNV;
    class VariablesNV;

// NEMO model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateNV :
        public util::Printable,
        private util::ObjectCounter<StateNV> {
public:
    static const std::string classname() {return "nv::StateNV";}

// Constructor, destructor
    StateNV(const GeometryNV &, const VariablesNV &, const util::DateTime &);
    StateNV(const GeometryNV &, const ModelNV &, const eckit::Configuration &);
    StateNV(const GeometryNV &, const StateNV &);
    StateNV(const StateNV &);
    StateNV(const GeometryNV &, const VariablesNV &, const eckit::Configuration &);
    virtual ~StateNV();
    StateNV & operator=(const StateNV &);

// Interpolate to observation location
//    void interpolate(const LocationsNV &, GomsNV &) const;

// Interpolate full fields
    void forceWith(const StateNV &) {}

// Interactions with Increment
    StateNV & operator+=(const IncrementNV &);

// I/O and diagnostics
    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;
    double norm() const {return fields_->norm();}

// Access to data
    const FieldsNV& getFields() const {return *fields_;}
    std::shared_ptr<FieldsNV> getFieldsPtr() const {return fields_;}
    F90flds& toFortran() {return fields_->toFortran();}
    const F90flds& toFortran() const {return fields_->toFortran();}
    std::shared_ptr<const GeometryNV> geometry() const {
        return fields_->geometry();
    }

// Time
    util::DateTime validTime() {return fields_->validTime();}
    const util::DateTime validTime() const {return fields_->validTime();}
    util::DateTime refTime() {return fields_->refTime();}
    const util::DateTime refTime() const {return fields_->refTime();}
    void updateTime(const util::Duration & dt) {fields_->updateTime(dt);}

// Access to time step
    void setTstep(util::Duration tstep) {fields_->setTstep(tstep);}
    const util::Duration & getTstep() {return fields_->getTstep();}

private:
    void print(std::ostream &) const;

//    std::unique_ptr<FieldsNV> fields_;
    std::shared_ptr<FieldsNV> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace nv

#endif // NV_STATENV_H
