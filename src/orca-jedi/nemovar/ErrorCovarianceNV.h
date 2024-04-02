/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_ErrorCovarianceNV_H
#define NV_ErrorCovarianceNV_H

#include <memory>

#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "orca-jedi/nemovar/GeometryNV.h"
#include "orca-jedi/nemovar/FortranNV.h"
#include "orca-jedi/nemovar/VariablesNV.h"

// Forward declarations
namespace eckit {
    class Configuration;
}

namespace nv {
    class Increment;
//    class State;
    class IncrementNV;
//    class IncrModCtlVecNV;
    class StateNV;

// Background error covariance matrix for NEMO model.

class ErrorCovarianceNV :
        public util::Printable,
        private eckit::NonCopyable,
        private util::ObjectCounter<ErrorCovarianceNV> {
public:
    static const std::string classname() {return "nv::ErrorCovarianceNV";}

// Constructor, destructor
    ErrorCovarianceNV(const GeometryNV &, const VariablesNV &,
                        const eckit::Configuration &, const StateNV &);
    ErrorCovarianceNV(const GeometryNV &, const VariablesNV &,
                        const eckit::Configuration &);
    ErrorCovarianceNV(const VariablesNV &, const eckit::Configuration &);
    ~ErrorCovarianceNV();

// Methods
    void linearize(const StateNV &, const GeometryNV &);
    void multiply(const IncrementNV &, IncrementNV &) const;
    void inverseMultiply(const IncrementNV &, IncrementNV &) const;
//    void multiplySqrt(const IncrModCtlVecNV &, IncrementNV &) const;
//    void multiplySqrtTrans(const IncrementNV &, IncrModCtlVecNV &) const;
    void randomize(IncrementNV &) const;

private:
    void print(std::ostream &) const;

    std::unique_ptr<const GeometryNV> geom_;
    const VariablesNV ctrlvars_;
    util::DateTime time_;
    F90bmat keyBmat_;
};

}  // namespace nv

#endif // NV_ErrorCovarianceNV_H
