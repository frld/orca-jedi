/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "VariablesNV.h"

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/util/Logger.h"

using oops::Log;

namespace nv {

// -----------------------------------------------------------------------------

VariablesNV::VariablesNV(const eckit::Configuration & conf)
    : keyVar_(0)
{
    Log::trace() << "VariablesNV: conf is: " << conf << std::endl;
    if(conf.has("Variables")) {
        Log::trace() << "VariablesNV: conf has Variables " << std::endl;
        const eckit::LocalConfiguration varConf(conf, "Variables");
        const eckit::Configuration * pconf = &varConf;
        nv_variables_create_f90(keyVar_, &pconf);
    } else {
        Log::trace() << "VariablesNV: conf does not have Variables " << std::endl;
        const eckit::Configuration * pconf = &conf;
        nv_variables_create_f90(keyVar_, &pconf);
    }
    ASSERT(keyVar_ != 0);
    
    variablesnv_int = 1;
}

// -----------------------------------------------------------------------------

VariablesNV::VariablesNV(const int ivar)
    : keyVar_(0)
{
    Log::trace() << "VariablesNV: ivar = " << ivar << std::endl;
    nv_variables_setup_f90(keyVar_, ivar);
    ASSERT(keyVar_ != 0);
    
    variablesnv_int = 2;
}

// -----------------------------------------------------------------------------

VariablesNV::~VariablesNV() {
    nv_variables_delete_f90(keyVar_);
}

// -----------------------------------------------------------------------------

VariablesNV::VariablesNV(const VariablesNV & other) {
    Log::trace() << "VariablesNV: clone" << std::endl;
    ASSERT(other.keyVar_ != 0);
    nv_variables_clone_f90(keyVar_, other.keyVar_);
    ASSERT(keyVar_ != 0);
    
    variablesnv_int = 3;
}

// -----------------------------------------------------------------------------

void VariablesNV::print(std::ostream & os) const {
    os << "VariablesNV::print not yet implemented." << std::endl;
}

// -----------------------------------------------------------------------------

}   // namespace nv

