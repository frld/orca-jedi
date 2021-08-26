/*
 * (C) British Crown Copyright 2017-2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>
#include <memory>
#include <math.h>
#include <vector>
#include <functional>
#include <numeric>
#include <iostream>
#include <sstream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/parallel/mpi/mpi.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "eckit/mpi/DataType.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateReadUtils.h"
#include "orca-jedi/model/ModelBias.h"
#include "orca-jedi/model/Model.h"


namespace orcamodel {

// Constructor, destructor
State::State(const Geometry & geom,
             const oops::Variables & vars,
             const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(vars), time_(time) 
{

  stateFields_ = atlas::FieldSet();

  for (size_t i=0; i<vars_.size(); ++i){
    // find the nemo name corresponding to the state variable name
    auto nemo_var_name = geom_->nemo_var_name(vars_[i]);
    // if it isn't already in stateFields, add it
    if (!stateFields_.has_field(nemo_var_name)){
      stateFields_.add( geom_->funcSpace().createField<double> (
            atlas::option::name(nemo_var_name) ) );
    }
  }

  oops::Log::trace() << "State(ORCA)::State created." << std::endl;
}

State::State(const Geometry & geom,
             const eckit::Configuration & conf)
  : geom_(new Geometry(geom))
    , vars_(conf, "state variables")
    , time_(util::DateTime(conf.getString("date")))
    , stateFields_()
{
  std::stringstream config_stream; 
  config_stream << "orcamodel::State:: config " << conf;
  oops::Log::debug() << config_stream.str() << std::endl;
  for (size_t i=0; i<vars_.size(); ++i){
    // find the nemo name corresponding to the state variable name
    auto nemo_var_name = geom_->nemo_var_name(vars_[i]);
    // if it isn't already in stateFields, add it
    if (!stateFields_.has_field(nemo_var_name)){
      stateFields_.add( geom_->funcSpace().createField<double> (
            atlas::option::name(nemo_var_name) ) );
    }
  }
  if (conf.has("analytic_init")) {
    this->analytic_init(conf, *geom_);
  } else {

    if ( !conf.has("nemo field file") ) {
      oops::Log::trace() << "State(ORCA):: nemo field file not in configuration"
                         << std::endl;
      oops::Log::trace() << "State(ORCA):: Zeroing state rather than reading from file"
                         << std::endl;
      this->zero();
    } else {
      readFieldsFromFile(conf, *geom_, stateFields_, "background");
      readFieldsFromFile(conf, *geom_, stateFields_, "background variance");
    }

  }
  oops::Log::trace() << "State(ORCA)::State created." << std::endl;
}

State::State(const Geometry & resol,
             const State & other)
{ 
  std::string err_message =
    "orcamodel::State::State created by interpolation Not implemented ";
  throw eckit::NotImplemented(err_message, Here());

  oops::Log::trace() << "State(ORCA)::State created by interpolation."
                     << "For now we are assuming that there is no resolution"
                     << "change"
                     << std::endl;
}

State::State(const State & other) {
  std::string err_message =
      "orcamodel::State::State created by copy not implemented ";
  throw eckit::NotImplemented(err_message, Here());
  oops::Log::trace() << "State(ORCA)::State copied." << std::endl;
}

State::~State() {
  oops::Log::trace() << "State(ORCA)::State destructed." << std::endl;
}

// Basic operators

State & State::operator=(const State & rhs) {
  std::string err_message =
      "orcamodel::State::State::operator= not implemented";
  throw eckit::NotImplemented(err_message, Here());
  return *this;
}

// Interpolate full fields

void State::changeResolution(const State & other) {
  std::string err_message =
      "orcamodel::State::State::changeResolution not implemented";
  throw eckit::NotImplemented(err_message, Here());
  oops::Log::trace() << "State(ORCA)::change resolution done" << std::endl;
}

// Interactions with Increments

State & State::operator+=(const Increment & dx) {
  std::string err_message =
      "orcamodel::State::State::operator+= not implemented";
  throw eckit::NotImplemented(err_message, Here());
  oops::Log::trace() << "State(ORCA)::add increment starting" << std::endl;
  oops::Log::trace() << "State(ORCA)::add increment done" << std::endl;
  return *this;
}

// I/O and diagnostics

void State::read(const eckit::Configuration & config) {
  oops::Log::trace() << "State(ORCA)::read starting" << std::endl;
  if ( !config.has("nemo field file") ) {
    oops::Log::trace() << "State(ORCA):: nemo field file not in configuration"
                       << std::endl;
    throw eckit::AssertionFailed("State(ORCA):: cannot find field file in configuration", Here());
  } else {
    readFieldsFromFile(config, *geom_, stateFields_, "background");
  }
  oops::Log::trace() << "State(ORCA)::read done" << std::endl;
}

void State::analytic_init(const eckit::Configuration & config,
                          const Geometry & geom) {
  oops::Log::trace() << "State(ORCA)::analytic_init starting" << std::endl;
  //this->zero();
  oops::Log::trace() << "State(ORCA)::analytic_init done" << std::endl;
}

void State::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "State(ORCA)::write starting" << std::endl;
  std::string err_message =
      "orcamodel::State::State::write not implemented";
  throw eckit::NotImplemented(err_message, Here());
}

void State::print(std::ostream & os) const {
  oops::Log::trace() << "State(ORCA)::print starting" << std::endl;

  os << std::endl << " Model state valid at time: " << validTime() << std::endl;
  os << "    " << vars_ <<  std::endl;
  os << "    " << "atlas field norms: ";
  for (atlas::Field field : stateFields_) {
    std::string fieldName = field.name();
    os << fieldName << ": " << norm(fieldName) << " ";
  }
  os <<  std::endl;

  oops::Log::trace() << "State(ORCA)::print done" << std::endl;
}

// For accumulator

void State::zero() {
  oops::Log::trace() << "State(ORCA)::zero starting" << std::endl;

  for (atlas::Field field : stateFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::State::zero:: field name = " << fieldName
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 1>( field );
    for ( atlas::idx_t j = 0; j < field_view.shape( 0 ); ++j ) {
      field_view( j ) = 0;
    }
  }

  oops::Log::trace() << "State(ORCA)::zero done" << std::endl;
}

double State::norm(const std::string & field_name) const {

  double norm = 0;

  auto field_view = atlas::array::make_view<double, 1>( stateFields_[field_name] );
  for ( atlas::idx_t j = 0; j < field_view.shape( 0 ); ++j ) {
    norm += field_view( j )*field_view( j ) ;
  }
  return sqrt(norm)/field_view.shape( 0 );
}

void State::accumul(const double & zz, const State & xx) {
  oops::Log::trace() << "State(ORCA)::accumul starting" << std::endl;

  std::string err_message =
      "orcamodel::State::State::accumul not implemented";
  throw eckit::NotImplemented(err_message, Here());

  oops::Log::trace() << "State(ORCA)::accumul done" << std::endl;
}

}  // namespace orcamodel