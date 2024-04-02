/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/nemovar/GeometryNV.h"

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"

#include "eckit/mpi/Comm.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"


namespace orcamodel {

oops::Variables orcaVariableFactory(const eckit::Configuration & config) {
  OrcaGeometryParameters params;
  params.validateAndDeserialize(config);

  std::vector<int> channels{};
  std::vector<std::string> names{};
  for (const NemoFieldParameters& nemoVariable :
        params.nemoFields.value()) {
    std::string  name = nemoVariable.name.value();
    if (std::find(names.begin(), names.end(), name) == names.end()) {
      names.push_back(name);
    }
  }

  return oops::Variables(names, channels);
}


// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm) :
                      comm_(comm), vars_(orcaVariableFactory(config)),
                      n_levels_(config.getInt("number levels")),
                      grid_(config.getString("grid name")),
                      usenemovar_(config.getBool("use nemovar", false))
{
    oops::Log::debug() << "Setting up orca-jedi geometry DJL" << std::endl;         // DJL

    params_.validateAndDeserialize(config);
    int64_t halo = params_.sourceMeshHalo.value().value_or(0);
    auto meshgen_config = grid_.meshgenerator()
                          | atlas::option::halo(halo);

    atlas::MeshGenerator meshgen(meshgen_config);
    auto partitioner_config = grid_.partitioner();
    partitioner_config.set("type",
        params_.partitioner.value().value_or("serial"));
    partitioner_ = atlas::grid::Partitioner(partitioner_config);
    mesh_ = meshgen.generate(grid_, partitioner_);
    funcSpace_ = atlas::functionspace::NodeColumns(
        mesh_, atlas::option::halo(halo));

    oops::Log::debug() << "Looking at setting up orca-jedi nemovar geometry DJL" << std::endl;         // DJL
    if (usenemovar_) {
       oops::Log::debug() << "Now setting up orca-jedi nemovar geometry DJL" << std::endl;         // DJL
       nvgeom_.reset(new nv::GeometryNV(config));
    }

    // Fill extra geometry fields for BUMP / SABER
    // these are area (needed?), vunit, hmask, gmask
    extraFields_ = atlas::FieldSet();

    // Vertical unit - DJL set to something sensible?
    atlas::Field vunit = funcSpace_.createField<double>(
      atlas::option::name("vunit") | atlas::option::levels(n_levels_));
    auto field_view = atlas::array::make_view<double, 2>(vunit);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
         field_view(j, k) = 1.;
      }
    }
    extraFields_->add(vunit);

    // halo mask / owned
//    atlas::Field hmask = funcSpace_.createField<int>(
// DJL temporary <int to <double * 2
    atlas::Field hmask = funcSpace_.createField<int32_t>(
      atlas::option::name("owned") | atlas::option::levels(n_levels_));
    auto ghost = atlas::array::make_view<int32_t, 1>(mesh_.nodes().ghost());

    auto field_view1 = atlas::array::make_view<int32_t, 2>(hmask);
    for (atlas::idx_t j = 0; j < field_view1.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view1.shape(1); ++k) {
        int x, y;
        std::tie(x, y) = xypt(j);
        // DJL hardwired to orca2 needs generalising
        if (ghost(j) || x < 2 || y >= 146 ) {field_view1(j, k) = 0;
//                      if (k==0) {oops::Log::debug() << "hmask ghost point " << j << std::endl; }
                      }           // 0 mask, 1 ocean 
        else {field_view1(j,k) = 1;}
      }
    }
    // Add field
    oops::Log::debug() << "Geometry adding hmask (set to all ocean except halo)" << std::endl;         // DJL
    extraFields_->add(hmask);

    // geometry mask
//    atlas::Field hmask = funcSpace_.createField<int>(
// DJL temporary <int to <double * 2
    atlas::Field gmask = funcSpace_.createField<int32_t>(
      atlas::option::name("gmask") | atlas::option::levels(n_levels_));
//    auto ghost = atlas::array::make_view<int32_t, 1>(mesh_.nodes().ghost());

    auto field_view2 = atlas::array::make_view<int32_t, 2>(gmask);
    for (atlas::idx_t j = 0; j < field_view2.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view2.shape(1); ++k) {
        int x, y;
        std::tie(x, y) = xypt(j);
        // DJL hardwired to orca2 needs generalising
        if (ghost(j) || x < 2 || y >= 146 ) {field_view2(j, k) = 0;
//                      if (k==0) {oops::Log::debug() << "gmask ghost point " << j << std::endl; }
                      }           // 0 mask, 1 ocean 
        else {field_view2(j,k) = 1;}
      }
    }
    // Add field
    oops::Log::debug() << "Geometry adding gmask (set to all ocean except halo)" << std::endl;         // DJL
    extraFields_->add(gmask);


//    // land mask - use also to blank out unwanted parts of the grid in an attempt to remove duplicated points
//// DJL temporary <int to <double * 2
//    atlas::Field gmask = funcSpace_.createField<int32_t>(
//      atlas::option::name("gmask") | atlas::option::levels(n_levels_));
//    auto field_view2 = atlas::array::make_view<int32_t, 2>(gmask);
////    auto lonlat_view = atlas::array::make_view<double, 2>(funcSpace_.lonlat());
//
//    for (atlas::idx_t j = 0; j < field_view2.shape(0); ++j) {
//      for (atlas::idx_t k = 0; k < field_view2.shape(1); ++k) {
////        int x, y;
////        std::tie(x, y) = xypt(j);
////        if (y >= 84) {field_view2(j, k) = 0;}           // 0 mask, 1 ocean 
////        if (y >= 84 || x > 140 || x < 2 || y < 2 ) {field_view2(j, k) = 0;}           // 0 mask, 1 ocean 
////        else {field_view2(j,k) = 1;}
////      }
////      if (j > 955 && j < 1045) {
////        oops::Log::debug() << lonlat_view(j,0) << "," << lonlat_view(j,1) << " "; 
//      field_view2(j,k) = 1;
//      }
//    }
//    // Add field
//    // DJL

//    oops::Log::debug() << "Geometry adding gmask (set all ocean)" << std::endl;         // DJL
   
//    extraFields_->add(gmask);

    atlas::Field area = funcSpace_.createField<double>(
      atlas::option::name("area") | atlas::option::levels(n_levels_));
    auto field_view3 = atlas::array::make_view<double, 2>(area);
 
    for (atlas::idx_t j = 0; j < field_view3.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view3.shape(1); ++k) {
        field_view3(j,k) = 4e10;           // DJL should change   2 degrees
      }
    }
 
    // Add field
    extraFields_->add(area);

    // create a nemovar geometry object
    //const eckit::Configuration * configc = &config
    //nv::GeometryF90 nvgeometry(config);               // GeometryNV
    // put in GeometryNV DJL
    
//    oops::Log::debug() << "Setting up NEMOVAR geometry DJL" << std::endl;         // DJL
//    nv::GeometryNV _nvgeometry(config);
// set shared pointer
//    nvgeometry_.reset(new nv::GeometryNV(config));
//    nvgeom_avail_=31;
//    
//    oops::Log::debug() << "Trying to grab the nvgeometry back DJL" << std::endl;         // DJL
//    
//    nv::GeometryNV nvgeom2 = this->getNVgeometry();   // Test
//
//    oops::Log::debug() << "Finished trying to grab the nvgeometry back DJL" << std::endl;         // DJL
    
}

// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_), 
    vars_(other.vars_),
    n_levels_(other.n_levels_),
    params_(other.params_),
    grid_(other.grid_),
    partitioner_(other.partitioner_),
    mesh_(other.mesh_),
    funcSpace_(other.funcSpace_),
    extraFields_(other.extraFields_),
    nvgeom_(other.nvgeom_)
{
// should be copying atlas things and the pointer to the Nemovar/F90 geometry

//    geomnv_int = 3;  // DJL
    oops::Log::trace() << "Geometry::Geometry copy ctor done" << std::endl;
}



// -----------------------------------------------------------------------------
Geometry::~Geometry() {

    oops::Log::trace() << "Geometry::Geometry copy dtor done" << std::endl;
}

const std::string Geometry::nemo_var_name(const std::string std_name) const {
  for (const auto & nemoField : params_.nemoFields.value()) {
    if (std_name == nemoField.name.value()) return nemoField.nemoName.value();
  }
  return std_name; // DJL
  std::stringstream err_stream;
  err_stream << "orcamodel::Geometry::nemo_var_name variable name \" ";
  err_stream << "\" " << std_name << " not recognised. " << std::endl;
  throw eckit::BadValue(err_stream.str(), Here());
}

// -----------------------------------------------------------------------------
/// \brief Give the number of levels for each provided level - surface variables
///        have 1 level, volumetric variables have "number levels" levels.
/// \param[in]     vars  variables to check.
/// \return        vector of number of levels in each variable.
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const
{
  std::vector<size_t> varSizes(vars.size());
  std::fill(varSizes.begin(), varSizes.end(), 0);

  auto nemoFields = params_.nemoFields.value();

  for (size_t i=0; i < vars.size(); ++i) {
    for (const auto & nemoField : nemoFields) {
      if (nemoField.name.value() == vars[i]) {
        if (nemoField.modelSpace.value() == "surface") {
          varSizes[i] = 1;
        } else {
          varSizes[i] = n_levels_;
        }
      }
    }
    if (varSizes[i] == 0) {
      std::stringstream err_stream;
      err_stream << "orcamodel::Geometry::variableSizes variable name \" ";
      err_stream << "\" " << vars[i] << " not recognised. " << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }
  }
  return varSizes;
}

void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
    const bool halo) const {
  const auto lonlat = atlas::array::make_view<double, 2>(funcSpace_.lonlat());
  const auto ghosts = atlas::array::make_view<int32_t, 1>(
      mesh_.nodes().ghost());
  const auto haloDistance = atlas::array::make_view<int32_t, 1>(
      mesh_.nodes().halo());
  auto isRequired = [&](const size_t nodeElem) {
    if (halo) {
      return !ghosts(nodeElem) || (haloDistance(nodeElem) > 0);
    }
    return !ghosts(nodeElem);
  };
  const size_t npts = funcSpace_.size();
  for (size_t nodeElem = 0; nodeElem < npts; ++nodeElem) {
    if (isRequired(nodeElem)) {
      lons.emplace_back(lonlat(nodeElem, 0));
      lats.emplace_back(lonlat(nodeElem, 1));
    }
  }
}

// -----------------------------------------------------------------------------
/// \brief Give the space of nemo field for each variable - surface, volume or
///         vertical. at the moment we need this distinction to read 3D depth
///         data from a 1D array
/// \param[in]     vars  variables to check.
/// \return        vector of variable Nemo model spaces.
std::vector<std::string> Geometry::variableNemoSpaces(
    const oops::Variables & vars) const
{
  std::vector<std::string> varNemoSpaces(vars.size(), "");

  auto nemoFields = params_.nemoFields.value();

  for (size_t i=0; i < vars.size(); ++i) {
    for (const auto & nemoField : nemoFields) {
      if (nemoField.name.value() == vars[i]) {
        if (nemoField.modelSpace.value() == "surface" ||
            nemoField.modelSpace.value() == "volume" ||
            nemoField.modelSpace.value() == "vertical" ) {
          varNemoSpaces[i] = nemoField.modelSpace.value();
        } else {
            std::stringstream err_stream;
            err_stream << "orcamodel::Geometry::variableNemoSpaces modelSpace"
                       << " \"" << nemoField.modelSpace.value()
                       << "\" not recognised for field \""
                       << nemoField.name.value() << "\"." << std::endl;
            throw eckit::BadValue(err_stream.str(), Here());
        }
      }
    }
    if (varNemoSpaces[i] == "") {
      std::stringstream err_stream;
      err_stream << "orcamodel::Geometry::variableSizes variable name \"";
      err_stream << vars[i] << "\" not available in the state. ";
      err_stream << "Either add this state variable to the model ";
      err_stream << "configuration or remove the corresponding obs variable";
      err_stream << " from the filters configuration." << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }
  }
  return varNemoSpaces;
}

const oops::Variables & Geometry::variables() const {
  return vars_;
}

/// \brief Check if a variable's data is a member of a type (e.g if it can be
///        sourced from the background file, variance file, or MDT file).
/// \param[in]     variable_name  Name of variable.
/// \param[in]     variable_type  Type of variable.
/// \return        Boolean for membership.
const bool Geometry::variable_in_variable_type(std::string variable_name,
  std::string variable_type) const {
  auto nemoFields = params_.nemoFields.value();
  for (const auto & nemoField : nemoFields) {
    if (nemoField.name.value() == variable_name) {
      std::string type = nemoField.variableType.value().value_or("background");
      return type == variable_type;
    }
  }

  std::stringstream err_stream;
  err_stream << "orcamodel::Geometry::variable_in_variable_type variable name ";
  err_stream << "\"" << variable_name << "\" not recognised. " << std::endl;
  throw eckit::BadValue(err_stream.str(), Here());
}

void Geometry::print(std::ostream & os) const {
  os << "Not Implemented";
}

// Determine x,y location from jpt
// DJL hardwired to orca2 needs generalising
std::tuple<int, int> xypt(int jpt) { int xwid=182; int y = jpt / xwid; int x = jpt - y*xwid; return std::make_tuple(x,y); } 

/*
void Geometry::set_gmask(atlas::Field & maskfield) const {

    oops::Log::debug() << "Geometry setting gmask" << std::endl;         // DJL

//    atlas::Field gmask = extraFields_["gmask"];                          // DJL does this work?
    atlas::Field gmask = extraFields_.field("gmask");                          // DJL does this work?

    auto field_viewin = atlas::array::make_view<double, 2>(maskfield);
    auto field_view1 = atlas::array::make_view<int32_t, 2>(gmask);
//    auto lonlat_view = atlas::array::make_view<double, 2>(funcSpace_.lonlat());

    for (atlas::idx_t j = 0; j < field_view1.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view1.shape(1); ++k) {
//        int x, y;
//        std::tie(x, y) = xypt(j);
//        if (y >= 84) {field_view1(j, k) = 0;}           // 0 mask, 1 ocean 
//        if (y >= 84 || x > 140 || x < 2 || y < 2 ) {field_view1(j, k) = 0;}           // 0 mask, 1 ocean 
//        else {
        if (field_viewin(j,k) == 1){
        field_view1(j,k) = 1;
        } else {
        field_view1(j,k) = 0;
        }
//        }
      }
//      if (j > 955 && j < 1045) {
//        oops::Log::debug() << lonlat_view(j,0) << "," << lonlat_view(j,1) << " "; 
//      }
    }
} 
*/

void Geometry::set_gmask(atlas::Field & maskfield) const {

    oops::Log::debug() << "Geometry setting gmask from maskfield" << std::endl;         // DJL

    atlas::Field gmask = extraFields_.field("gmask");                          // DJL does this work?

    auto field_viewin = atlas::array::make_view<double, 2>(maskfield);
    auto field_view1 = atlas::array::make_view<int32_t, 2>(gmask);
//    auto lonlat_view = atlas::array::make_view<double, 2>(funcSpace_.lonlat());

    for (atlas::idx_t j = 0; j < field_view1.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view1.shape(1); ++k) {

// only change values that are currently unmasked ( 0 mask, 1 ocean )
        if (field_view1(j,k) == 1) {
           if (field_viewin(j,k) == 0) { 
              field_view1(j,k) = 0;
//              if (k == 0) {
//              oops::Log::debug() << "setting gmask " << j << "," << k << " gmask " << field_view1(j,k) << std::endl;
//              } 
              }
        }
      }
//      if (j > 955 && j < 1045) {
//        oops::Log::debug() << lonlat_view(j,0) << "," << lonlat_view(j,1) << " "; 
//      }
    }
} 


}  // namespace orcamodel
