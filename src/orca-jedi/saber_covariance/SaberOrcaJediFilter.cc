/* Copyright Met Office 2023
*/

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/StructuredColumns.h"

#include "orca-jedi/saber_covariance/SaberOrcaJediFilter.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "saber/oops/Utilities.h"

//include "oops/base/IdentityMatrix.h"
//include "oops/assimilation/GMRESR.h"
// DJL (below)
//include <boost/uuid/uuid.hpp>            // uuid class
//include <boost/uuid/uuid_generators.hpp> // generators
//include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
//include "orca-jedi/state/StateIOUtils.h"


namespace saber {
namespace orcajedifilter {

// -----------------------------------------------------------------------------

//static SaberOuterBlockMaker<OrcaJediFilter> makerOrcaJediFilter_(
//        "orca jedi filter");
static SaberCentralBlockMaker<OrcaJediFilter> makerOrcaJediFilter_(
        "orca jedi filter");

// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
/// Background error covariance matrix for orca model.

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
OrcaJediFilter::OrcaJediFilter(const oops::GeometryData & geometryData,
                                               const oops::Variables & activeVars,
                                               const eckit::Configuration & covarConf,
                                               const Parameters_ & params,
                                               const oops::FieldSet3D & xb,
                                               const oops::FieldSet3D & fg)
  : SaberCentralBlockBase(params, xb.validTime()), params_(params),
    geometryData_(geometryData),
    activeVars_(activeVars),
    ctlVecSize_(0)
//    activeVars_(getActiveVars(params, outerVars)),
//    innerGeometryData_(geometryData),
//    innerVars_(outerVars)

{

   oops::Log::trace() << "OrcaJediFilter Covariance setting up" << std::endl;

  // Compute total number of levels
  size_t nlev = 0;
  for (const std::string & var : activeVars.variables()) {
    nlev += activeVars.getLevels(var);
  }

  // Compute control vector size
  ctlVecSize_ = nlev*geometryData_.functionSpace().size();

   oops::Log::trace() << "OrcaJediFilter ctlVecSize " << ctlVecSize_ << std::endl;

   oops::Log::trace() << "OrcaJediFilter Covariance created" << std::endl;
             
}
// -----------------------------------------------------------------------------
void OrcaJediFilter::multiply(oops::FieldSet3D & fieldSet) const { 

//Increment & dxin, Increment & dxout) const {

// convert increments

// where is NEMOVAR geom stored?
// vars_nv from vars?
// initialise NV increments

   oops::Log::trace() << "OrcaJediFilter Covariance multiply" << std::endl;

//   oops::Log::trace() << "OrcaJediFilter Covariance multiply index " << std::endl;
//   dxout = dxin;
   
   atlas::Field field = fieldSet[0];
   auto field_view = atlas::array::make_view<double, 2>(field);
//   atlas::Field fieldin = fieldSet[0];
//   auto field_viewin = atlas::array::make_view<double, 2>(fieldin);

//   oops::Log::trace() << "Covariance multiply 4 len1 " << len1 << std::endl;
//   oops::Log::trace() << "Covariance multiply 4 field_view.shape(0) " << field_view.shape(0) << len1 << std::endl;


// Shapiro ?

//   if (covtyp_ == "Shapiro") {

      atlas::idx_t k = 0;

      auto array = new double[field_view.shape(0)][10];

      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
          array[j][k] = field_view(j,k);
      }

   //   float w1 = 0.5;
   //   float wa = 1 + w1*2;

      float w1 = 0.25;
      float wa = 1 + w1*4;

      for (int it = 0; it < 20; ++it) {
   //   for (atlas::idx_t j = 1; j < field_view.shape(0)-1; ++j) {
   ///      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
   //        field_view(j, k) = (array[j][k] + w1*(array[j-1][k]+array[j+1][k]))/wa;
   //      }
   //    }
   ///    }

         float xwid = 182;
         float ywid = 149;
         for (int xpt = 1; xpt < xwid-1 ; ++xpt){
            for (int ypt = 1; ypt < ywid-1; ++ypt) {
               int j0 = ypt*182 + xpt;    // DJL hardwired to work with orca2  
               int j1 = (ypt-1)*182 + xpt;  
               int j2 = (ypt+1)*182 + xpt;
               int j3 = ypt*182 + xpt-1;
               int j4 = ypt*182 + xpt+1;

               field_view(j0, k) = (array[j0][k] + w1*(array[j1][k]+array[j2][k]+array[j3][k]+array[j4][k]))/wa;      
            }
         }

         for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
             array[j][k] = field_view(j,k);
         }
      }

   // renormalise
      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
   //      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
           field_view(j, k) = field_view(j, k)*31.3688;              // 31.3688 for 20 iterations of the Shapiro
   //      }
       }

/*   } else if (covtyp_ == "Nemovar") {

      std::shared_ptr<const nv::GeometryNV> nvgeom = geom_.getNVgeometryPtr();

      util::DateTime time = dxin.time();

      oops::Log::trace() << "Covariance multiply 1" << std::endl;

      nv::IncrementNV dxin_nv(*nvgeom, *nvvars_, time);

   // Convert from Increment (orca-jedi / atlas) to IncrementNV

      int64_t len1 = 27118;

      std::vector<double> field1array(len1, 0);   

      atlas::Field fieldin = dxin.incrementFields()[0];
      auto field_view1 = atlas::array::make_view<double, 2>(fieldin);

      atlas::idx_t k1 = 0;
      for (atlas::idx_t j = 0; j < field_view1.shape(0); ++j) {
   //      for (atlas::idx_t k = 0; k < field_view1.shape(1); ++k) {
          field1array[j] = field_view1(j, k1); 
   //      }
       }

      dxin_nv.setarray(len1, field1array);

      oops::Log::trace() << "Covariance multiply 2" << std::endl;

      nv::IncrementNV dxout_nv(*nvgeom, *nvvars_, time);

      oops::Log::trace() << "Covariance multiply 3" << std::endl;

      (*nverrorcov_).multiply(dxin_nv, dxout_nv);

   // convert increments back to orcajedi form

   // copy increment
   //   Increment dxout(dxin);
   // paste in data from nv increment

   // convert from IncrementNV out to Increment (orca-jedi / atlas).
      nv::FieldsNV fields2nv = dxout_nv.getFields();

      std::vector<double> field2array;
      fields2nv.getarray(len1, field2array);   

      atlas::Field field = dxout.incrementFields()[0];
      auto field_view = atlas::array::make_view<double, 2>(field);

      oops::Log::trace() << "Covariance multiply 4 len1 " << len1 << std::endl;
      oops::Log::trace() << "Covariance multiply 4 field_view.shape(0) " << field_view.shape(0) << len1 << std::endl;

      atlas::idx_t k = 0;
      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
   //      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
           field_view(j, k) = field2array[j];
   //      }
      }

    } else {  
    
       std::string err_message = "The selected error covariance covtype "+covtyp_+" is not recognised";

       throw eckit::NotImplemented(err_message, Here());
    }
*/
//
//    boost::uuids::uuid uuid = boost::uuids::random_generator()();    
//
//    writeGenFieldsToFile("errorcov_dxin_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxin.validTime(), dxin.incrementFields(), FieldDType::Double);
//    
//    writeGenFieldsToFile("errorcov_dxout_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxout.validTime(), dxout.incrementFields(), FieldDType::Double);
     oops::Log::trace() << "OrcaJediFilter Covariance multiply end" << std::endl;
     }


OrcaJediFilter::~OrcaJediFilter() {
  oops::Log::trace() << classname() << "::~OrcaJediFilter dtor starting" << std::endl;
  oops::Log::trace() << classname() << "::~OrcaJediFilter dtor done" << std::endl;
}

// -----------------------------------------------------------------------------

void OrcaJediFilter::randomize(oops::FieldSet3D & fset) const {
  oops::Log::trace() << classname() << "::randomize starting" << std::endl;
  oops::Log::trace() << classname() << "::randomize warning empty" << std::endl;

  // Consistency check
  for (const auto & var : activeVars_.variables()) {
      ASSERT(fset.has(var));
  }

  // Random initialization
//  fset.randomInit(geometryData_.functionSpace(), activeVars_);

  oops::Log::trace() << classname() << "::randomize done" << std::endl;
}


/* oops::FieldSet3D OrcaJediFilter::generateInnerFieldSet(
  const oops::GeometryData & innerGeometryData,
  const oops::Variables & innerVars) const {
  oops::FieldSet3D fset(this->validTime(), innerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(innerGeometryData.comm(),
                                           innerGeometryData.functionSpace(),
                                           innerVars));

   oops::Log::trace() << "OrcaJediFilter generateInnerFieldSet done" << std::endl;

  return fset;
}

// -----------------------------------------------------------------------------

oops::FieldSet3D OrcaJediFilter::generateOuterFieldSet(
  const oops::GeometryData & outerGeometryData,
  const oops::Variables & outerVars) const {
  oops::FieldSet3D fset(this->validTime(), outerGeometryData.comm());
  fset.deepCopy(util::createSmoothFieldSet(outerGeometryData.comm(),
                                           outerGeometryData.functionSpace(),
                                           outerVars));

   oops::Log::trace() << "OrcaJediFilter generateOuterFieldSet done" << std::endl;

  return fset;
}
*/


// -----------------------------------------------------------------------------
void OrcaJediFilter::print(std::ostream & os) const {
  os << classname();
}

}  // namespace orcajedifilter
}  // namespace saber
