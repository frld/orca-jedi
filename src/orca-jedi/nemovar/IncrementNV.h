/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_INCREMENTNV_H
#define NV_INCREMENTNV_H

#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/base/GeneralizedDepartures.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/dot_product.h"

#include "nemovar/FieldsNV.h"
#include "nemovar/FortranNV.h"

namespace eckit {
    class Configuration;
}

namespace nv {
    class GeometryNV;
//    class GomsNV;
//    class LocationsNV;
    class StateNV;
    class VariablesNV;

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models
 */

// -----------------------------------------------------------------------------

class IncrementNV :
        public util::Printable,
        public oops::GeneralizedDepartures,
        private util::ObjectCounter<IncrementNV> {
public:
    static const std::string classname() {return "nv::IncrementNV";}

// Constructor, destructor
    IncrementNV(const GeometryNV &, const VariablesNV &, const util::DateTime &);
    IncrementNV(const GeometryNV &, const IncrementNV &);
    IncrementNV(const IncrementNV &, const bool copy = true);
    virtual ~IncrementNV();

// Interpolate to observation location
//    void interpolateTL(const LocationsNV &, GomsNV &) const;
//    void interpolateAD(const LocationsNV &, const GomsNV &);

// Interactions with State
    void diff(const StateNV &, const StateNV &);

// Linear algebra operators
    void zero();
    void zeroTlFields();
    void zero(const util::DateTime &);
    IncrementNV & operator =(const IncrementNV &);
    IncrementNV & operator+=(const IncrementNV &);
    IncrementNV & operator-=(const IncrementNV &);
    IncrementNV & operator*=(const double &);
    void axpy(const double &, const IncrementNV &, const bool check = true);
    double dot_product_with(const IncrementNV &) const;
    void schur_product_with(const IncrementNV &);
    void random() {};
    void accumul(const double &, const StateNV &);

// I/O and diagnostics
    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;
    double norm() const {return fields_.norm();}
    double max(const VariablesNV & var) const {return fields_.max(var);}
    double min(const VariablesNV & var) const {return fields_.min(var);}

// Geometry
    std::shared_ptr<const GeometryNV> geometry() const {
        return fields_.geometry();
    }

// Access to data
    const FieldsNV & getFields() const {return fields_;}  
//    std::shared_ptr<FieldsNV> getFieldsPtr() const {return fields_;}
    F90flds& toFortran() {return fields_.toFortran();}
    const F90flds& toFortran() const {return fields_.toFortran();}

    void setarray(int64_t & len, const std::vector<double> & array) {fields_.setarray(len, array);}

// Time
    util::DateTime validTime() {return fields_.validTime();}
    const util::DateTime validTime() const {return fields_.validTime();}
    util::DateTime refTime() {return fields_.refTime();}
    const util::DateTime refTime() const {return fields_.refTime();}
    void updateTime(const util::Duration & dt) {
      fields_.updateTime(dt);}

// Access to time step
    void setTstep(util::Duration tstep) {fields_.setTstep(tstep);}
    const util::Duration & getTstep() {return fields_.getTstep();}

// Balanced/un-balanced
    void toTlFields(const int & skipBalance = 0) {fields_.toTlFields(skipBalance);}
    void toIncrFields() {fields_.toIncrFields();}

private:
    void print(std::ostream &) const;

    FieldsNV fields_;
//    std::shared_ptr<FieldsNV> fields_;
};

// -----------------------------------------------------------------------------

}  // namespace nv

#endif // NV_INCREMENTNV_H
