/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/base/Variables.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/OptionalParameter.h"

namespace orcamodel {

class OrcaStateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaStateParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> nemoFieldFile{"nemo field file", this};
  oops::OptionalParameter<std::string> varianceFieldFile{
    "variance field file", "", this};
  oops::OptionalParameter<std::string> maskFieldFile{
    "mask field file", "", this};
  oops::OptionalParameter<bool> analyticInit{"analytic initialisation", this};
  oops::OptionalParameter<std::string> outputNemoFieldFile{
    "output nemo field file", "", this};
  oops::RequiredParameter<util::DateTime> date{"date", this};
  oops::RequiredParameter<oops::Variables> stateVariables{"state variables",
    this};
};

}  //  namespace orcamodel
