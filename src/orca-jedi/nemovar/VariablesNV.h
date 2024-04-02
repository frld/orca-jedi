/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_VARIABLESNV_H
#define NV_VARIABLESNV_H

#include "eckit/exception/Exceptions.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "orca-jedi/nemovar/FortranNV.h"

// Forward declarations
namespace eckit {
    class Configuration;
}

namespace nv {

// VariablesNV class to handle variable lists for NEMO model.

class VariablesNV :
        public util::Printable,
        private util::ObjectCounter<VariablesNV> {
public:
    static const std::string classname() {return "nv::VariablesNV";}

// Constructor, destructor
    explicit VariablesNV(const eckit::Configuration & conf);
    VariablesNV(const int ivar);
    VariablesNV(const VariablesNV & other);
    ~VariablesNV();

// Access to data
    F90vars& toFortran() {return keyVar_;}
    const F90vars& toFortran() const {return keyVar_;}
    
    int variablesnv_int;

private:
    void print(std::ostream &) const;

    F90vars keyVar_;
};

}   // namespace nv

#endif // NV_VARIABLESNV_H
