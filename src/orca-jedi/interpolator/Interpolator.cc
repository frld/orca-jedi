/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <cstddef>
#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation.h"
#include "atlas/functionspace.h"
#include "atlas/field/MissingValue.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/missingValues.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateIOUtils.h"
#include "orca-jedi/increment/Increment.h"

#include "orca-jedi/interpolator/Interpolator.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
}

namespace orcamodel {
  class State;
  class Geometry;

  atlas::functionspace::PointCloud atlasObsFuncSpaceFactory(
      const std::vector<double>& lats,
      const std::vector<double>& lons) {
      size_t nlocs = lats.size();

      // Setup observation functionspace
      oops::Log::trace() << "orcamodel::Interpolator:: creating "
                         << "atlasObsFuncSpace with nlocs = " << nlocs
                         << std::endl;
      atlas::Field points("lonlat", atlas::array::make_datatype<double>(),
          atlas::array::make_shape(nlocs, 2));
      auto arrv_t = atlas::array::make_view<double, 2>(points);
      for (unsigned int j = 0; j < nlocs; ++j) {
        arrv_t(j, 1) = lats[j];
        arrv_t(j, 0) = lons[j];
      }
      oops::Log::trace() << "orcamodel::Interpolator:: creating "
                         << "atlasObsFuncSpace ... done" << std::endl;

      return atlas::functionspace::PointCloud(std::move(points));
  }

  Interpolator::Interpolator(const eckit::Configuration & conf,
      const Geometry & geom, const std::vector<double>& lats,
      const std::vector<double>& lons) :
      nlocs_(lats.size()),
      atlasObsFuncSpace_(atlasObsFuncSpaceFactory(lats, lons)),
      interpolator_(eckit::LocalConfiguration(conf, "atlas-interpolator"),
                    geom.functionSpace(),
                    atlasObsFuncSpace_),
      comm_(geom.getComm()) {
    params_.validateAndDeserialize(conf);
//    fileCounter_ = 0;
    oops::Log::trace() << "orcamodel::Interpolator:: conf:" << conf
                       << std::endl;
    if (nlocs_ == 0) {
      oops::Log::trace() << "orcamodel::Interpolator:: nlocs == 0" << std::endl;
    }
  }

  void Interpolator::apply(const oops::Variables& vars, const State& state,
      const std::vector<bool> & mask,
      std::vector<double>& result) const {
    oops::Log::trace() << "[" << comm_.rank()
                       << "] orcamodel::Interpolator::apply starting "
                       << std::endl;

    const size_t nvars = vars.size();

    for (size_t j=0; j < nvars; ++j) {
      if (!state.variables().has(vars[j])) {
        std::stringstream err_stream;
        err_stream << "orcamodel::Interpolator::apply varname \" "
                   << "\" " << vars[j]
                   << " not found in the model state." << std::endl;
        err_stream << "    add the variable to the state variables and "
                   << "add a mapping from the geometry to that variable."
                   << std::endl;
        throw eckit::BadParameter(err_stream.str(), Here());
      }
    }

    const std::vector<size_t> varSizes =
      state.geometry()->variableSizes(vars);
    oops::Log::debug() << "orcamodel::Interpolator::apply nvars = " << nvars
                       << " nlocs_ = " << nlocs_;
    size_t nvals = 0;
    for (size_t jvar=0; jvar < nvars; ++jvar) {
      nvals += nlocs_ * varSizes[jvar];
      oops::Log::debug() << " varSizes[" << jvar << "] = " << varSizes[jvar];
    }
    oops::Log::debug() << " nvals = " << nvals << std::endl;
    result.resize(nvals);

    auto res_iter = result.begin();
    for (size_t jvar=0; jvar < nvars; ++jvar) {
      const auto execute = [&](auto typeVal) {
        using T = decltype(typeVal);
        executeInterpolation<T>(vars[jvar], varSizes[jvar], state, res_iter);
      };

      ApplyForFieldType(execute,
                        state.geometry()->fieldPrecision(vars[jvar]),
                        std::string("orcamodel::Interpolator::apply '")
                          + vars[jvar] + "' field type not recognised");
    }

    // input is state output is result

    ASSERT(result.size() == nvals);
    oops::Log::trace() << "orcamodel::Interpolator::apply done "
                       << std::endl;
  }

  void Interpolator::apply(const oops::Variables& vars, const Increment& inc,
             const std::vector<bool> & mask,
             std::vector<double>& result) const {

    // DJL question can the atlas templates help?

    oops::Log::trace() << "orcamodel::Interpolator::apply increment "
                       << std::endl;

    // input is inc output is result

    const size_t nvars = vars.size();

    for (size_t j=0; j < nvars; ++j) {
      if (!inc.variables().has(vars[j])) {
        std::stringstream err_stream;
        err_stream << "orcamodel::Interpolator::apply varname \" "
                   << "\" " << vars[j]
                   << " not found in the model increment." << std::endl;
        err_stream << "    add the variable to the increment variables and "
                   << "add a mapping from the geometry to that variable."
                   << std::endl;
        throw eckit::BadParameter(err_stream.str(), Here());
      }
    }

    const std::vector<size_t> varSizes =
      inc.geometry()->variableSizes(vars);
    size_t nvals = 0;
    for (size_t jvar=0; jvar < nvars; ++jvar) nvals += nlocs_ * varSizes[jvar];
    result.resize(nvals);

    std::size_t out_idx = 0;
    for (size_t jvar=0; jvar < nvars; ++jvar) {
      auto gv_varname = vars[jvar];
      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
          atlas::option::name(gv_varname) |
          atlas::option::levels(varSizes[jvar]));
      interpolator_.execute(inc.incrementFields()[gv_varname], tgt_field);
      auto field_view = atlas::array::make_view<double, 2>(tgt_field);
      atlas::field::MissingValue mv(inc.incrementFields()[gv_varname]);
      bool has_mv = static_cast<bool>(mv);
      for (std::size_t klev=0; klev < varSizes[jvar]; ++klev) {
        for (std::size_t iloc=0; iloc < nlocs_; iloc++) {
          if (has_mv && mv(field_view(iloc, klev))) {
            result[out_idx] = util::missingValue<double>();
          } else {
            result[out_idx] = field_view(iloc, klev);
          }
          ++out_idx;
        }
      }
    }

    oops::Log::trace() << "orcamodel::Interpolator::apply increment done "
                       << std::endl;

  }

  void Interpolator::applyAD(const oops::Variables& vars, Increment& inc,
               const std::vector<bool> & mask,
               const std::vector<double> & resultin) const
  {

    // input is resultin output is inc

    oops::Log::trace() << "orcamodel::Interpolator::applyAD start "
                       << std::endl;

//    throw eckit::NotImplemented("Adjoint interpolation not implemented",
//                                Here());

    oops::Log::debug() << "DJL ** Interpolator::applyAD this needs checking **" << std::endl;

// ** Not sure what I'm doing yet - Just trying to do the opposite of applyAD

    const size_t nvars = vars.size();

    for (size_t j=0; j < nvars; ++j) {
      if (!inc.variables().has(vars[j])) {
        std::stringstream err_stream;
        err_stream << "orcamodel::Interpolator::apply varname \" "
                   << "\" " << vars[j]
                   << " not found in the model increment." << std::endl;
        err_stream << "    add the variable to the increment variables and "
                   << "add a mapping from the geometry to that variable."
                   << std::endl;
        throw eckit::BadParameter(err_stream.str(), Here());
      }
    }

    const std::vector<size_t> varSizes =
      inc.geometry()->variableSizes(vars);
    size_t nvals = 0;

//    boost::uuids::uuid uuid = boost::uuids::random_generator()();    
//    std::shared_ptr<const Geometry> geom = inc.geometry();
//    writeGenFieldsToFile("applyADpre"+ boost::uuids::to_string(uuid) +".nc", *geom, inc.validTime(), inc.incrementFields());

    for (size_t jvar=0; jvar < nvars; ++jvar) nvals += nlocs_ * varSizes[jvar];
//    result.resize(nvals);

    std::size_t out_idx = 0;
    for (size_t jvar=0; jvar < nvars; ++jvar) {
      oops::Log::debug() << "DJL ** jvar " << jvar << " " << nvars
                         << "varSizes " << varSizes[jvar]
                         << std::endl;
      auto gv_varname = vars[jvar];
//      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
//          atlas::option::name(gv_varname) |
//          atlas::option::levels(varSizes[jvar]));
//      atlas::Field incField = inc.incrementFields()[jvar];

      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
        atlas::option::name(gv_varname) |
        atlas::option::levels(varSizes[jvar]));

      auto field_view = atlas::array::make_view<double, 2>(tgt_field);


      for (std::size_t klev=0; klev < varSizes[jvar]; ++klev) {
        for (std::size_t iloc=0; iloc < nlocs_; iloc++) {
//          if (has_mv && mv(field_view(iloc, klev))) {
////            result[out_idx] = util::missingValue<double>();
//          } else {
          oops::Log::debug() << "DJL iloc " << iloc << " klev " << klev << " out_idx "
                            << out_idx << " resultin[out_idx] " << resultin[out_idx] << std::endl;
          field_view(iloc, klev) = resultin[out_idx];
//          }
          ++out_idx;
        }
      }


//      atlas::field::MissingValue mv(inc.incrementFields()[gv_varname]);
//      bool has_mv = static_cast<bool>(mv);
      interpolator_.execute_adjoint(inc.incrementFields()[jvar], tgt_field);

    }    // jvar

    // DJL write to file
    
    boost::uuids::uuid uuid = boost::uuids::random_generator()();    
    std::shared_ptr<const Geometry> geom = inc.geometry();
//    std::ostringstream out;
//    out << std::setfill('0') << std::setw(6) << uuid;
    writeGenFieldsToFile("applyAD"+ boost::uuids::to_string(uuid) +".nc", *geom, inc.validTime(), inc.incrementFields());
//    fileCounter_++;

    oops::Log::trace() << "orcamodel::Interpolator::applyAD done "
                       << std::endl;


  }



  /// \brief Execute the atlas interpolation on a single field.
  /// \param gv_varname The GeoVaLs variable name of the source field.
  /// \param var_size The number of observation locations in the interpolant.
  /// \param state The state object containing the model state.
  /// \param iter Reference to the interator into the output vector.
  template<class T> void Interpolator::executeInterpolation(
      const std::string& gv_varname,
      size_t var_size,
      const State& state,
      std::vector<double>::iterator& iter) const {
    atlas::Field tgt_field = atlasObsFuncSpace_.createField<T>(
        atlas::option::name(gv_varname) |
        atlas::option::levels(var_size));
    interpolator_.execute(state.stateFields()[gv_varname], tgt_field);
    auto field_view = atlas::array::make_view<T, 2>(tgt_field);
    atlas::field::MissingValue mv(state.stateFields()[gv_varname]);
    bool has_mv = static_cast<bool>(mv);
    for (std::size_t klev=0; klev < var_size; ++klev) {
      for (std::size_t iloc=0; iloc < nlocs_; iloc++) {
        if (has_mv && mv(field_view(iloc, klev))) {
          *iter = util::missingValue<double>();
        } else {
          *iter = static_cast<double>(field_view(iloc, klev));
        }
        ++iter;
      }
    }
  }

  template void Interpolator::executeInterpolation<double>(
      const std::string& gv_varname,
      size_t var_size,
      const State& state,
      std::vector<double>::iterator& iter) const;
  template void Interpolator::executeInterpolation<float>(
      const std::string& gv_varname,
      size_t var_size,
      const State& state,
      std::vector<double>::iterator& iter) const;

  void Interpolator::print(std::ostream & os) const {
    os << "orcamodel::Interpolator: " << std::endl;
    os << "  Obs function space " << atlasObsFuncSpace_ << std::endl;
    os << "  Interpolator " << interpolator_ << std::endl;
  }

}  // namespace orcamodel
