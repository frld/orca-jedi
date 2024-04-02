/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_FIELDSNV_H
#define NV_FIELDSNV_H

#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "nemovar/FortranNV.h"
#include "nemovar/GeometryNV.h"
#include "nemovar/VariablesNV.h"

// Forward declarations
namespace eckit {
    class Configuration;
}

namespace nv {
//    class GomsNV;
//    class LocationsNV;

// FieldsNV class to handle a set of nemo fields.

class FieldsNV :
        public util::Printable,
        private util::ObjectCounter<FieldsNV> {
public:
    static const std::string classname() {return "nv::FieldsNV";}

// Constructors and basic operators
    FieldsNV(const GeometryNV &, const VariablesNV &, const util::DateTime &);
    FieldsNV(const FieldsNV &, const bool copyFields = true);
    ~FieldsNV();

    FieldsNV & operator=(const FieldsNV &);
    FieldsNV & operator+=(const FieldsNV &);
    FieldsNV & operator-=(const FieldsNV &);
    FieldsNV & operator*=(const double &);
    void zero();
    void zeroTlFields();
    void zero(const util::DateTime &);
    void axpy(const double &, const FieldsNV &);
    double dot_product_with(const FieldsNV &) const;
    void schur_product_with(const FieldsNV &);

// Interpolate to given location
//    void interpolate(const LocationsNV &, GomsNV &) const;
//    void interpolateTL(const LocationsNV &, GomsNV &) const;
//    void interpolateAD(const LocationsNV &, const GomsNV &);

// Interpolate full fields
    void changeResolution(const FieldsNV &);
    void add(const FieldsNV &);
    void diff(const FieldsNV &, const FieldsNV &);

// Utilities
    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &) const;
    void getarray(int64_t &, std::vector<double> &);                      // see lfric serialize and deserialize
    void setarray(int64_t &, const std::vector<double> &);
    double norm() const;
    double max(const VariablesNV &) const {return 0.0;}
    double min(const VariablesNV &) const {return 0.0;}

// Access to data
    const VariablesNV& getVariables() const {return *vars_;}
    F90flds& toFortran() {return keyFlds_;}
    const F90flds& toFortran() const {return keyFlds_;}
    F90fldstl& toFortranTl() {return keyFldstl_;}
    const F90fldstl& toFortranTl() const {return keyFldstl_;}
    std::shared_ptr<const GeometryNV> geometry() const {return geom_;}

// Access to time step
    void setTstep(util::Duration tstep) {tstep_ = tstep;}
    const util::Duration & getTstep() const {return tstep_;}

    util::DateTime validTime() {return time_;}
    const util::DateTime validTime() const {return time_;}
    util::DateTime refTime() {return reftime_;}
    const util::DateTime refTime() const {return reftime_;}
    void updateTime(const util::Duration & dt) {
      ASSERT(std::abs(dt.toSeconds()) == tstep_.toSeconds());
      time_ += dt;
    }

// Balanced/un-balanced
    void toTlFields(const int & skipBalance = 0);
    void toIncrFields();

private:
    void print(std::ostream &) const;
    util::Duration tstep_;
    util::DateTime reftime_;
    util::DateTime time_;

    std::shared_ptr<const GeometryNV> geom_;
    std::shared_ptr<const VariablesNV> vars_;
    F90flds keyFlds_;
    F90fldstl keyFldstl_;
};

std::ostream & operator<< (std::ostream &, const FieldsNV &);

}  // namespace nv

#endif // NV_FIELDSNV_H
