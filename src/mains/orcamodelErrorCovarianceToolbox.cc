/*
 * (C) British Crown Copyright 2022 MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/Library.h"
#include "oops/generic/instantiateModelFactory.h"

#include "saber/oops/ErrorCovarianceToolbox.h"
//include "oops/runs/Run.h"
#include "orca-jedi/run/Run.h"
#include "saber/oops/instantiateCovarFactory.h"

#include "orca-jedi/utilities/OrcaModelTraits.h"

int main(int argc,  char ** argv) {
//  oops::Run run(argc, argv);
  orcamodel::Run run(argc, argv);
  atlas::Library::instance().initialise();
  saber::instantiateCovarFactory<orcamodel::OrcaModelTraits>();

  saber::ErrorCovarianceToolbox<orcamodel::OrcaModelTraits> var;
  int i = run.execute(var);
  atlas::Library::instance().finalise();
  return i;
}
