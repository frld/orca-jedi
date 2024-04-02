/*
 * (C) British Crown Copyright 2022 MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/Library.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/Run.h"
#include "saber/oops/ErrorCovarianceTraining.h"
#include "saber/oops/instantiateCovarFactory.h"
#include "orca-jedi/utilities/OrcaModelTraits.h"


int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  //  oops::instantiateModelFactory<orcamodel::OrcaModelTraits>();   // required?
  atlas::Library::instance().initialise();
  saber::instantiateCovarFactory<orcamodel::OrcaModelTraits>();
  saber::ErrorCovarianceTraining<orcamodel::OrcaModelTraits> dir;
  int i = run.execute(dir);
  atlas::Library::instance().finalise();
  return i;
}
