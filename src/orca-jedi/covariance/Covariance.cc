/* Copyright Met Office 2023
*/

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/StructuredColumns.h"

#include "orca-jedi/covariance/Covariance.h"
#include "orca-jedi/nemovar/ErrorCovarianceNV.h"

#include "orca-jedi/state/State.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "nemovar/GeometryNV.h"
#include "nemovar/FortranNV.h"
#include "nemovar/VariablesNV.h"
#include "nemovar/StateNV.h"
#include "nemovar/IncrementNV.h"

#include "oops/util/Logger.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/assimilation/GMRESR.h"
// DJL (below)
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include "orca-jedi/state/StateIOUtils.h"



namespace orcamodel {

// -----------------------------------------------------------------------------
/// Background error covariance matrix for orca model.

// -----------------------------------------------------------------------------
Covariance::Covariance(const Geometry & geom, const oops::Variables & vars,
             const eckit::Configuration & conf, const State & x1, const State & x2)
             : geom_(geom), vars_(vars), covtyp_(conf.getString("covtype"))
{

// convert state

   oops::Log::trace() << "Covariance setting up" << std::endl;

   oops::Log::trace() << "DJL covtype " << covtyp_ << std::endl;

// where is NEMOVAR geom stored?
// vars_nv from vars?
// Need a new constructor to zero state

//   nv::VariablesNV vars_nv(eckit::LocalConfiguration(conf, "Variables"));
//   nv::VariablesNV nvvars_(eckit::LocalConfiguration(conf, "state variables"));          // surely won't work

//   nvvars_.reset(new nv::VariablesNV(0));
   
//   nv::ErrorCovarianceNV errorcovnv(*nvvars_, conf);    // Setup Nemovar covariance
//   nv::ErrorCovarianceNV errorcovnv();    // Setup Nemovar covariance


//   nv::VariablesNV nvvars_(conf);
//   nvvars_.reset(new nv::VariablesNV(conf));

   if (covtyp_ == "Nemovar") {

      nvvars_.reset(new nv::VariablesNV(0));

      oops::Log::trace() << "DJL variablesnv_int " << (*nvvars_).variablesnv_int << std::endl;

   // Alternative to put oops::Variables in VariablesNV

   //   oops::Log::trace() << "DJL nvgeom_avail_ " << geom_.nvgeom_avail_ << std::endl;

   //   nvgeom_.reset(geom_.getNVgeometry());

   //   nv::GeometryNV nvgeom = geom_.getNVgeometry(); // DJL try shared pointer instead
       std::shared_ptr<const nv::GeometryNV> nvgeom = geom_.getNVgeometryPtr();

   //   nv::GeometryNV nvgeom = nv::GeometryNV(conf);

      oops::Log::trace() << "DJL geomnv_int " << (*nvgeom).geomnv_int << std::endl;


   //    = std::make_shared<nv::GeometryNV>(geom.getNVgeometry());
   //   nvgeom_ = geom.getNVgeometryPtr();   // DJL

      oops::Log::trace() << "DJL creating stateNV x1nv " << std::endl;

      nv::StateNV x1nv(*nvgeom, *nvvars_, conf);

      oops::Log::trace() << "DJL calling ErrorCovarianceNV4" << std::endl;

   //   nv::ErrorCovarianceNV nverrorcov4(*nvgeom, *nvvars_, conf, x1nv);  //, x2nv);  // test
      nverrorcov_.reset(new nv::ErrorCovarianceNV(*nvgeom, *nvvars_, conf, x1nv));  //, x2nv);

   }


   oops::Log::trace() << "Covariance created" << std::endl;
             
}
// -----------------------------------------------------------------------------
void Covariance::multiply(const Increment & dxin, Increment & dxout) const {

// convert increments

// where is NEMOVAR geom stored?
// vars_nv from vars?
// initialise NV increments

   oops::Log::trace() << "Covariance multiply" << std::endl;

   oops::Log::trace() << "Covariance multiply index " << std::endl;
//   dxout = dxin;
   
   atlas::Field field = dxout.incrementFields()[0];
   auto field_view = atlas::array::make_view<double, 2>(field);
   atlas::Field fieldin = dxin.incrementFields()[0];
   auto field_viewin = atlas::array::make_view<double, 2>(fieldin);

//   oops::Log::trace() << "Covariance multiply 4 len1 " << len1 << std::endl;
//   oops::Log::trace() << "Covariance multiply 4 field_view.shape(0) " << field_view.shape(0) << len1 << std::endl;


// Shapiro ?

   if (covtyp_ == "Shapiro") {

      oops::Log::trace() << "Covariance multiply shapiro" << std::endl;

      atlas::idx_t k = 0;

      auto array = new double[field_view.shape(0)][10];

      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
          array[j][k] = field_viewin(j,k);
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

   } else if (covtyp_ == "Nemovar") {

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

    boost::uuids::uuid uuid = boost::uuids::random_generator()();    

    writeGenFieldsToFile("errorcov_dxin_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxin.validTime(), dxin.incrementFields(), FieldDType::Double);
    
    writeGenFieldsToFile("errorcov_dxout_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxout.validTime(), dxout.incrementFields(), FieldDType::Double);
   
   oops::Log::trace() << "Covariance multiply end" << std::endl;

}
// -----------------------------------------------------------------------------
void Covariance::inverseMultiply(const Increment & dxi, Increment & dxo) const {

   oops::Log::trace() << "Covariance inverse multiply" << std::endl;
   
   // DJL test
//   multiply(dxin, dxout);

// adapted/copied from saber/ErrorCovariance.h

     // Iterative inverse
   oops::IdentityMatrix<Increment> Id;

   oops::Log::trace() << "Covariance inverse multiply identity matrix done" << std::endl;

   dxo.zero();

   oops::Log::trace() << "Covariance inverse multiply identity zero done" << std::endl;

    boost::uuids::uuid uuid = boost::uuids::random_generator()();    

    writeGenFieldsToFile("cov_invmult_dxi_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxi.validTime(), dxi.incrementFields(), FieldDType::Double);
    
    writeGenFieldsToFile("cov_invmult_dxo1_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxo.validTime(), dxo.incrementFields(), FieldDType::Double);

   oops::GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);

    writeGenFieldsToFile("cov_invmult_dxo2_"+ boost::uuids::to_string(uuid) +".nc", geom_, dxo.validTime(), dxo.incrementFields(), FieldDType::Double);

   oops::Log::trace() << "Covariance inverse multiply done" << std::endl;

}
// -----------------------------------------------------------------------------
void Covariance::randomize(Increment & dx) const {

   oops::Log::trace() << "Covariance randomize" << std::endl;

}
// -----------------------------------------------------------------------------

}  // namespace orcamodel
