/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <memory>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace.h"
#include "atlas/mesh.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

#include "eckit/mpi/Comm.h"
#include "eckit/log/Timer.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

#include "orca-jedi/geometry/GeometryParameters.h"

#include "orca-jedi/nemovar/GeometryNV.h"

#include "oops/util/Logger.h"
#include "orca-jedi/utilities/Types.h"

namespace atlas {
  class Field;
  class FieldSet;
  class Mesh;
}

namespace orcamodel {

// -----------------------------------------------------------------------------
/// Geometry handles for ORCA model.

  oops::Variables orcaVariableFactory(const eckit::Configuration & config);

class Geometry : public util::Printable {
 public:
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  ~Geometry();
  Geometry(const Geometry &); 

  std::vector<size_t> variableSizes(const oops::Variables &) const;
  std::vector<std::string> variableNemoSpaces(const oops::Variables & vars)
      const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  const oops::Variables & variables() const;
  void latlon(std::vector<double> & lats, std::vector<double> & lons,
              const bool halo) const;
  const atlas::FunctionSpace & functionSpace() const {return funcSpace_;}
  const atlas::FieldSet & extraFields() const {return extraFields_;}
  const atlas::FieldSet & fields() const {return extraFields_;}
  atlas::FieldSet & extraFields() {return extraFields_;}

  const atlas::Grid & grid() const {return grid_;}
  const atlas::Mesh & mesh() const {return mesh_;}
  const std::string nemo_var_name(const std::string std_name) const;
  const bool variable_in_variable_type(std::string variable_name,
    std::string variable_type) const;
  bool levelsAreTopDown() const {return true;}
  std::string distributionType() const {
      return params_.partitioner.value().value_or("serial");}
  void set_gmask(atlas::Field &) const;
//  void set_hmask(atlas::Field &) const;
////  const nv::GeometryNV & getNVgeometry() const {return nvgeom_;}
  std::shared_ptr<const nv::GeometryNV> getNVgeometryPtr() const {return nvgeom_;}
//  const nv::GeometryNV & getNVgeometry() const {ASSERT(nvgeom_); return *nvgeom_;}
  const nv::GeometryNV & getNVgeometry() const {oops::Log::trace() << "getNVgeometry DJL" << std::endl; return *nvgeom_;}

// Access to data
//    F90geom& toFortran() {return geom_->toFortran();}
//    const F90geom& toFortran() const {return geom_->toFortran();}

//  int nvgeom_avail_;          // DJL
      return params_.partitioner.value();}
  FieldDType fieldPrecision(std::string variable_name) const;
  std::shared_ptr<eckit::Timer> timer() const {return eckit_timer_;}
  void log_status() const;

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  oops::Variables vars_;
  size_t n_levels_;
  OrcaGeometryParameters params_;
  atlas::Grid grid_;
  bool usenemovar_;
  atlas::grid::Partitioner partitioner_;
  atlas::Mesh mesh_;
  atlas::functionspace::NodeColumns funcSpace_;
    atlas::FieldSet extraFields_;
//  nv::GeometryNV nvgeometry_;    // nv::GeometryNV DJL?
    std::shared_ptr<nv::GeometryNV> nvgeom_;

//    std::shared_ptr<GeometryF90> geom_;
  atlas::FieldSet nofields_;
  std::shared_ptr<eckit::Timer> eckit_timer_;
};

std::tuple<int, int> xypt(int jpt);

// -----------------------------------------------------------------------------

}  // namespace orcamodel
