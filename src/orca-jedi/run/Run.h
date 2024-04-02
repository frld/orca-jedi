/*
 * (C) British Crown Copyright 2017-2018 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


#ifndef SRC_RUN_ORCAMODEL_H_
#define SRC_RUN_ORCAMODEL_H_

#include "oops/runs/Run.h"

namespace orcamodel {

/*!
 *  Run encapsulates one unifiedmodel/OOPS run.
 */

// -----------------------------------------------------------------------------

class Run : public oops::Run {
 public:
  Run(int, char **);
  ~Run();
};

// -----------------------------------------------------------------------------

}  // namespace orcamodel
#endif  // SRC_RUN_ORCAMODEL_H_
