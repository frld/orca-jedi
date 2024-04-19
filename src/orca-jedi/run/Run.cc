/*
 * (C) British Crown Copyright 2017-2018 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/runs/Run.h"

#include "nemovar/FortranNV.h"

#include "Run.h"

namespace orcamodel {
// -----------------------------------------------------------------------------
Run::Run(int argc, char ** argv) : oops::Run(argc, argv) {
  oops::Log::trace() << "Creating Run(orcamodel) DJL" << std::endl;
  const eckit::Configuration * conf = &config();
  // unifiedmodel_init_f90(&conf);

// DJL could do with something to skip if not required
  oops::Log::trace() << "Calling NEMOVAR_INIT" << std::endl;

  nv::nv_init_f90(&conf);
    
  oops::Log::trace() << "Run(orcamodel) created" << std::endl;
}
// -----------------------------------------------------------------------------
Run::~Run() {
  oops::Log::trace() << "Destructing Run(orcamodel) DJL" << std::endl;
//  unifiedmodel_finalize_f90();
  oops::Log::trace() << "Run(orcamodel): MPI finalized" << std::endl;
}
// -----------------------------------------------------------------------------

} // namespace orcamodel
