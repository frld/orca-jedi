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

class NemoFieldParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NemoFieldParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> name {"name", this};
  oops::RequiredParameter<std::string> nemoName {"nemo field name", this};
  oops::RequiredParameter<std::string> modelSpace {"model space", this};
  oops::OptionalParameter<std::string> variableType {"variable type", this};
};

class OrcaGeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaGeometryParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::vector<NemoFieldParameters>> nemoFields
    {"nemo variables", this};
  oops::RequiredParameter<std::string> gridName
    {"grid name", this};
  oops::RequiredParameter<int> nLevels {"number levels", this};
  oops::OptionalParameter<bool> useNemovar {"use nemovar", this};
  oops::OptionalParameter<int> sourceMeshHalo {"source mesh halo", this};
  oops::OptionalParameter<std::string> partitioner {"partitioner", this};
};

}  //  namespace orcamodel
