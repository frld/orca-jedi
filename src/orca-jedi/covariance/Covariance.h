/*
 * (C) Copyright 2022  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

//include "orca-jedi/covariance/CovarianceParameters.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/nemovar/GeometryNV.h"
#include "orca-jedi/nemovar/VariablesNV.h"
#include "orca-jedi/nemovar/ErrorCovarianceNV.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace orcamodel {
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for orca model.

class Covariance : public util::Printable,
                   private boost::noncopyable,
                   private util::ObjectCounter<Covariance> {
 public:
//  typedef CovarianceParameters Parameters_;
  static const std::string classname() {return "orcamodel::Covariance";}

  Covariance(const Geometry &, const oops::Variables &,
             const eckit::Configuration &, const State &, const State &);

  void multiply(const Increment &, Increment &) const;
  void inverseMultiply(const Increment &, Increment &) const;
  void randomize(Increment &) const;

 private:
  std::shared_ptr<nv::GeometryNV> nvgeom_;
  std::shared_ptr<nv::VariablesNV> nvvars_;
  std::shared_ptr<nv::ErrorCovarianceNV> nverrorcov_;
  Geometry geom_;
  oops::Variables vars_;
  std::string covtyp_;
//  util::DateTime time_;
  void print(std::ostream & os) const {os << "Covariance";}
// -----------------------------------------------------------------------------

};
// -----------------------------------------------------------------------------

}  // namespace orcamodel
