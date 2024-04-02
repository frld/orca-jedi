/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_MODELNV_H
#define NV_MODELNV_H

#include <string>

#include "eckit/memory/NonCopyable.h"

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "nemovar/FortranNV.h"
#include "nemovar/GeometryNV.h"

// Forward declarations
namespace eckit {
    class Configuration;
}

namespace nv {
    class ModelBiasNV;
    class StateNV;

/*!
 *  NEMO nonlinear model definition and configuration parameters.
 */
// NEMO model definition.

class ModelNV:
        public util::Printable,
        private eckit::NonCopyable,
        private util::ObjectCounter<ModelNV> {
public:
    static const std::string classname() {return "nv::ModelNV";}

// Constructors, destructors;
    ModelNV(const GeometryNV &, const eckit::Configuration &);
    ~ModelNV();

// Model integration
    void initialize(StateNV &) const;
    void step(StateNV &, const ModelBiasNV &) const;
    void finalize(StateNV &) const;

    F90traj saveTrajectory(StateNV &) const;

// Other utilities
    const util::Duration & timeResolution() const {return tstep_;}

private:
    void print(std::ostream &) const;

    const GeometryNV & geom_;
    util::Duration tstep_;
    F90conf nvConfig_;
};

// -----------------------------------------------------------------------------

}  // namespace nv

#endif // NV_MODELNV_H
