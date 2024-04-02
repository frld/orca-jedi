/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "ModelNV.h"

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "nemovar/FortranNV.h"
#include "nemovar/GeometryNV.h"
#include "nemovar/StateNV.h"


namespace nv {

// -----------------------------------------------------------------------------

ModelNV::ModelNV(const GeometryNV & resol, const eckit::Configuration & config)
    : geom_(resol)
{
    oops::Log::trace() << "ModelNV::ModelNV" << std::endl;
    tstep_ = util::Duration(config.getString("tstep"));
    oops::Log::trace() << "ModelNV::ModelNV tstep is " << tstep_ << std::endl;

    //const eckit::Configuration * configc = &config;
    //nemo_setup_f90(&configc, geom_.toFortran(), nvConfig_);

    oops::Log::trace() << "ModelNV::ModelNV ctor not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

ModelNV::~ModelNV() {
    //nemo_delete_f90(nvConfig_);
    oops::Log::trace() << "ModelNV::~ModelNV dtor not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelNV::initialize(StateNV & xx) const {
    xx.setTstep(tstep_);
    //nemo_prepare_integration_f90(nvConfig_, xx.toFortran());
    oops::Log::trace() << "ModelNV::initialize not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelNV::step(StateNV & xx, const ModelBiasNV &) const {
    //nemo_propagate_f90(nvConfig_, xx.toFortran());
    xx.updateTime(tstep_);
    oops::Log::info() << xx << std::endl;
}

// -----------------------------------------------------------------------------

void ModelNV::finalize(StateNV & xx) const {
    oops::Log::trace() << "ModelNV::finalize not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

F90traj ModelNV::saveTrajectory(StateNV & xx) const {
    F90traj ftraj = 0;
    //nemo_prop_traj_f90(nvConfig_, xx.toFortran(), ftraj);
    //ASSERT(ftraj != 0);
    oops::Log::trace() << "ModelNV::saveTrajectory not implemented" << std::endl;
    return ftraj;
}

// -----------------------------------------------------------------------------

void ModelNV::print(std::ostream & os) const {
  os << "ModelNV::print not yet implemented." << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace nv
