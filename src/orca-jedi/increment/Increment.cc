/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <ostream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/StructuredColumns.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/log/CodeLocation.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/GeoVaLs.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateIOUtils.h"
#include "orca-jedi/increment/Increment.h"
// DJL (below)
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include "orca-jedi/state/StateIOUtils.h"

namespace orcamodel {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                     const oops::Variables & vars,
                     const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(vars), time_(time),
    incrementFields_()
{
  incrementFields_ = atlas::FieldSet();

  setupIncrementFields();
  
  this->zero();   // DJL needed?

  oops::Log::debug() << "Increment(ORCA)::Increment created for "<< validTime()
                     << std::endl;

// DJL Debugging output (save extra fields)

//  const atlas::FieldSet & extrafs = geom.extraFields();
//  writeGenFieldsToFile("test.nc", geom, time, extrafs);

}

// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                     const Increment & other)
  : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  std::string err_message =
      "orcamodel::Increment::constructor(geom, other) not implemented";
  throw eckit::NotImplemented(err_message, Here());
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy)
  //: geom_(other.geom_), vars_(other.vars_), time_(other.time_), incrementFields_()  
  //(other.incrementFields_)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{  
  oops::Log::debug() << "Increment(ORCA)::Increment copy " << copy << std::endl;

  incrementFields_ = atlas::FieldSet();

  setupIncrementFields();

  if (copy) {

  for (size_t i=0; i < vars_.size(); ++i) {

    // copy variable from _Fields to new field set
    atlas::Field field = other.incrementFields_[i];
    oops::Log::debug() << "Copying increment field " << field.name() << std::endl;
    incrementFields_->add(field);

  }
  
  }


  // copy from other to incrementFields
/*    auto ghost = atlas::array::make_view<int32_t, 1>(
        geom_->mesh().nodes().ghost());
    for (int i = 0; i< other.incrementFields_.size();i++)
    {
      atlas::Field field = incrementFields_[i];
      atlas::Field field1 = other.incrementFields_[i];
      std::string fieldName = field.name();
      std::string fieldName1 = field1.name();
      oops::Log::debug() << "orcamodel::Increment::copy:: field name = " << fieldName
                         << " field name 1 = " << fieldName1
                         << std::endl;
      auto field_view = atlas::array::make_view<double, 2>(field);
      auto field_view1 = atlas::array::make_view<double, 2>(field1);
      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
  //        if (!ghost(j)) field_view(j, k) += field_view1(j, k);
          field_view(j,k) = field_view1(j,k);
        }
      }

      int jpt = ypt*182 + xpt;  // DJL hardwired to work with orca2  
      oops::Log::debug() << "DJL value " << jpt << " from " << field_view1(jpt, 0) << " to " << field_view(jpt, 0) << std::endl;

    }
*/
  
//  int xpt = 91;  int ypt = 49;  
  
  oops::Log::debug() << "Increment(ORCA)::Increment copied." << std::endl;

  oops::Log::debug() << "increment copy self print";
  print(oops::Log::debug());
  oops::Log::debug() << "increment copy other print";
  other.print(oops::Log::debug());  

}

// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_(other.incrementFields_)
{

//  Increment(other, 1); 

//  std::string err_message =
//      "orcamodel::Increment::copy constructor not implemented";
//  throw eckit::NotImplemented(err_message, Here());
}

/// Basic operators
Increment & Increment::operator=(const Increment & rhs) {
  time_ = rhs.time_;
  incrementFields_ = rhs.incrementFields_;
  vars_ = rhs.vars_;
  geom_.reset();
  geom_ = rhs.geom_;

  oops::Log::debug() << "Increment(ORCA)::= copy ended" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
//  lfric_increment_self_add_f90(keyInc_, dx.keyInc_);

//   incrementFields_ -= dx.incrementFields_;

  oops::Log::debug() << "increment add self print";
  print(oops::Log::debug());
  oops::Log::debug() << "increment add dx print";
  dx.print(oops::Log::debug());  

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size();i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
//  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::add:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
//        if (!ghost(j)) field_view(j, k) += field_view1(j, k);
        field_view(j,k) += field_view1(j,k);
      }
    }
  }

//    boost::uuids::uuid uuid = boost::uuids::random_generator()();    
//
//    writeGenFieldsToFile("increment_add_intern_"+ boost::uuids::to_string(uuid) +".nc", *geom_, validTime(), incrementFields_);
//
//    writeGenFieldsToFile("increment_add_dx_"+ boost::uuids::to_string(uuid) +".nc", *geom_, dx.validTime(), dx.incrementFields());

  oops::Log::debug() << "increment add self print";
  print(oops::Log::debug());
  oops::Log::debug() << "increment add dx print";
  dx.print(oops::Log::debug());  

  oops::Log::debug() << "Increment(ORCA)::+ add ended" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator-=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
//  lfric_increment_self_sub_f90(keyInc_, dx.keyInc_);

//  incrementFields_ -= dx.incrementFields_;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size();i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
//  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::subtract:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
//        if (!ghost(j)) field_view(j, k) -= field_view1(j, k);
        field_view(j, k) -= field_view1(j, k);
      }
    }
  }

  oops::Log::debug() << "Increment(ORCA)::- subtract ended" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
Increment & Increment::operator*=(const double & zz) {
//  lfric_increment_self_mul_f90(keyInc_, zz);

//  incrementFields_ *= zz;

// see state.cc

  oops::Log::debug() << "orcamodel::Increment:multiply start" << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::multiply:: field name = " << fieldName 
                       << " zz " << zz
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
//        if (!ghost(j)) field_view(j, k) *= zz;
        field_view(j, k) *= zz;
      }
    }
  }

  oops::Log::debug() << "Increment(ORCA)::* multiplication ended" << std::endl;
  return *this;
}

void Increment::diff(const State & x1, const State & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size();i++)
  {
    atlas::Field field1 = x1.getField(i);
    atlas::Field field2 = x2.getField(i);
    atlas::Field fieldi = incrementFields_[i];

    std::string fieldName1 = field1.name();
    std::string fieldName2 = field2.name();
    std::string fieldNamei = fieldi.name();
    oops::Log::debug() << "orcamodel::Increment::diff:: field name 1 = " << fieldName1
                       << " field name 2 = " << fieldName2
                       << " field name inc = " << fieldNamei                       
                       << std::endl;
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    auto field_view2 = atlas::array::make_view<double, 2>(field2);
    auto field_viewi = atlas::array::make_view<double, 2>(fieldi);
    for (atlas::idx_t j = 0; j < field_viewi.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_viewi.shape(1); ++k) {
        if (!ghost(j)) { field_viewi(j, k) = field_view1(j, k) - field_view2(j, k); }
        else { field_viewi(j,k) = 0; }
      }
    }
  }

}

void Increment::setval(const double & val) {
  oops::Log::trace() << "Increment(ORCA)::setval starting" << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::setval:: field name = " << fieldName
                       << "value " << val
                       << std::endl;

    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
//          if (!ghost(j)) field_view(j, k) = 0;
        field_view(j, k) = val;
      }
    }
  }

  oops::Log::trace() << "Increment(ORCA)::setval done" << std::endl;
}

void Increment::zero() {

  oops::Log::trace() << "Increment(ORCA)::zero starting" << std::endl;

  this->setval(0);

  oops::Log::trace() << "Increment(ORCA)::zero done" << std::endl;
  
}

void Increment::ones() {

  oops::Log::trace() << "Increment(ORCA)::ones starting" << std::endl;
  
  this->setval(1);

  oops::Log::trace() << "Increment(ORCA)::ones done" << std::endl;
  
}
  
void Increment::zero(const util::DateTime & vt) {
  time_ = vt;

  oops::Log::debug() << "orcamodel::Increment::zero time " << vt << std::endl;

//  Increment::zero();
  this->zero();
}


/*
void Increment::ones() {
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::ones:: field name = " << fieldName
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
//        if (!ghost(j)) field_view(j, k) = 1;
        field_view(j, k) = 1;
      }
    }
  }

}
*/



void Increment::axpy(const double & zz, const Increment & dx, const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size();i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::subtract:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
//        if (!ghost(j)) field_view(j, k) += zz * field_view1(j, k);
        field_view(j, k) += zz * field_view1(j, k);
      }
    }
  }

}


double Increment::dot_product_with(const Increment & dx) const {

  double zz = 0;
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  // How should this deal with multiple fields only want to do this with one field??? DJL
  for (int i = 0; i< incrementFields_.size();i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::dot_product_with:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) zz += field_view(j, k) * field_view1(j, k);
      }
    }
  }

    boost::uuids::uuid uuid = boost::uuids::random_generator()();    

    writeGenFieldsToFile("increment_dotprod_intern_"+ boost::uuids::to_string(uuid) +".nc", *geom_, validTime(), incrementFields_);

    writeGenFieldsToFile("increment_dotprod_dx_"+ boost::uuids::to_string(uuid) +".nc", *geom_, dx.validTime(), dx.incrementFields());

  oops::Log::debug() << "orcamodel::Increment::dot_product_with ended :: zz = " << zz << std::endl;

return zz; }

void Increment::schur_product_with(const Increment & dx) {

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size();i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::schur_product_with:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) *= field_view1(j, k);
      }
    }
  }

}

void Increment::random() {

  oops::Log::debug() << "orcamodel::Increment::random start" << std::endl;
  oops::Log::debug() << "orcamodel::Increment::random seed_ " << seed_ << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::random:: field name = " << fieldName
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    // DJL using the same seed needs fixing...
    util::NormalDistribution<double> xx(field_view.shape(0)*field_view.shape(1), 0.0, 1.0, seed_);
    int idx=0;
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) = xx[idx];
        idx++;
      }
    }
  }


}

void Increment::dirac(const eckit::Configuration & conf) {
// Add a delta function at jpt, kpt
// DJL this is 0,0 in orca2 - assuming a single processor run
//  int xpt = 141;  int ypt = 73;
// DJL this is 30S, 100W in orca2 (-30N, 260E)
//  int xpt = 91;  int ypt = 49; 

// DJL hardwired to fixed point on the orca2 grid
// should be made more flexible in future    
  int jpt = ypt*182 + xpt;
  int kpt = 0; 

  oops::Log::debug() << "orcamodel::Increment::dirac:: delta function at jpt = " << jpt 
                       << " kpt " << kpt << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::dirac:: field name = " << fieldName
                       << std::endl;

    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (j == jpt && k == kpt) { 
            field_view(j, k) = 1;
          }
          else
          {
            field_view(j, k) = 0;
          }
        }
      }
    }
    
  }

}




// ------------------------------------------------------------------------------------------------
void Increment::toFieldSet(atlas::FieldSet & fset) const {
  oops::Log::debug() << "Increment toFieldSet starting" << std::endl;

//  lfricjedi_increment_to_fieldset_f90(keyInc_,
//                                      geom_->toFortran(),
//                                      vars_,
//                                      fset.get());

//  incrementFields_->toFieldSet(fset);

//  fset = incrementFields_;
// DJL

  fset = atlas::FieldSet();

  for (size_t i=0; i < vars_.size(); ++i) {
    // add variable if it isn't already in incrementFields
//    std::vector<size_t> varSizes = geom_->variableSizes(vars_);
//    if (!fset_.has(vars_[i])) {
//      fset_.add(geom_->functionSpace().createField<double>(
//           atlas::option::name(vars_[i]) |
//           atlas::option::levels(varSizes[i])));
//    }
    // copy variable from increments to new field set
    atlas::Field fieldinc = incrementFields_[i];
    std::string fieldName = fieldinc.name();
    oops::Log::debug() << "Copy increment toFieldSet " << fieldName << std::endl;

// debugging output  
    auto field_view = atlas::array::make_view<double, 2>(fieldinc);
//    int xpt = 91;  int ypt = 49;  
//    int jpt = ypt*182 + xpt;  // DJL hardwired to work with orca2  
//    oops::Log::debug() << "DJL value " << jpt << " " << field_view(jpt-1, 0) << " " << field_view(jpt, 0) << " " << field_view(jpt+1, 0) << std::endl;
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        int x, y;                                       // DJL
        std::tie(x, y) = xypt(j);
        if (abs(x-xpt) < 2 && abs(y-ypt) < 2) {
          oops::Log::debug() << "x " << x << " y " << y << " tofieldset field_view " << field_view(j,k) << std::endl; 
        }
      }
    }

    fset->add(fieldinc);
  }

  
  
  oops::Log::debug() << "Increment toFieldSet done" << std::endl;
}
// ------------------------------------------------------------------------------------------------
void Increment::toFieldSetAD(const atlas::FieldSet & fset) {
  oops::Log::debug() << "Increment toFieldSetAD starting" << std::endl;

  std::string err_message =
      "orcamodel::Increment::toFieldSetAD not implemented";
  throw eckit::NotImplemented(err_message, Here());

  oops::Log::debug() << "Increment toFieldSetAD done" << std::endl;
}
// ------------------------------------------------------------------------------------------------
void Increment::fromFieldSet(const atlas::FieldSet & fset) {
  oops::Log::debug() << "Increment fromFieldSet start" << std::endl;

//  lfricjedi_increment_from_fieldset_f90(keyInc_,
//                                        geom_->toFortran(),
//                                        vars_,
//                                        fset.get());


//  incrementFields_ = fset;

  for (int i = 0; i< fset.size();i++) {
    atlas::Field field = fset[i];
    atlas::Field fieldinc = incrementFields_[i];
    oops::Log::debug() << "Increment fromFieldSet field " << i << " " << field.name() << std::endl;
    oops::Log::debug() << "Increment fromFieldSet fieldinc " << i << " " << fieldinc.name() << std::endl;

// copy from field to incrementfields

    auto field_view_to = atlas::array::make_view<double, 2>(fieldinc);
    auto field_view_from = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view_to.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view_to.shape(1); ++k) {
        field_view_to(j, k) = field_view_from(j, k);
        int x, y;                                       // DJL
        std::tie(x, y) = xypt(j);
        if (abs(x-xpt) < 2 && abs(y-ypt) < 2) {
          oops::Log::debug() << "x " << x << " y " << y << " fromfieldset field_view " << field_view_to(j,k) << std::endl; 
        }
      }
    }

//    int xpt = 91;  int ypt = 49;  
    int jpt = ypt*182 + xpt;  // DJL hardwired to work with orca2  
    oops::Log::debug() << "DJL value " << jpt << " from " << field_view_from(jpt, 0) << " to " << field_view_to(jpt, 0) << std::endl;


  }

  oops::Log::debug() << "Increment fromFieldSet done" << std::endl;
}
// ------------------------------------------------------------------------------------------------


void Increment::setupIncrementFields() {
  for (size_t i=0; i < vars_.size(); ++i) {
    // add variable if it isn't already in incrementFields
    std::vector<size_t> varSizes = geom_->variableSizes(vars_);
    if (!incrementFields_.has(vars_[i])) {
      incrementFields_.add(geom_->functionSpace().createField<double>(
           atlas::option::name(vars_[i]) |
           atlas::option::levels(varSizes[i])));
      oops::Log::trace() << "Increment(ORCA)::setupIncrementFields : "
                         << vars_[i] << "has dtype: "
                         << (*(incrementFields_.end()-1)).datatype().str() << std::endl;
      geom_->log_status();
    }
  }
}


/// I/O and diagnostics
void Increment::read(const eckit::Configuration & conf) {
  std::string err_message =
      "orcamodel::Increment::read not implemented";
  throw eckit::NotImplemented(err_message, Here());
  
}

void Increment::write(const eckit::Configuration & conf) const {

  oops::Log::debug() << "orcamodel::increment::write" << std::endl;

  // OrcaStateParameters params 

  writeIncFieldsToFile(conf, *geom_, validTime(), incrementFields_);

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::write:: field name = " << fieldName
                       << std::endl;

    auto field_view = atlas::array::make_view<double, 2>(field);
//    int xpt = 141;  int ypt = 73;  
    int xpt = 91; int ypt = 49;
    int jpt = ypt*182 + xpt;  // DJL hardwired to work with orca2  
    oops::Log::debug() << "DJL value " << jpt << " " << field_view(jpt-1, 0) << " " << field_view(jpt, 0) << " " << field_view(jpt+1, 0) << std::endl;

//    auto field_view = atlas::array::make_view<double, 2>(field);
    
//    oops::Log::debug() << std::endl << "Increment::write DJL ";
//    for (int i=996; i < 1005; ++i) {
//      oops::Log::debug() << i << "," << field_view(i,0) << " ";
//    }
    oops::Log::debug() << std::endl;
  }

  std::cout << "orcamodel::increment::write end DJL" << std::endl;

}

void Increment::print(std::ostream & os) const {

  double norm;
  double min;
  double max;

  oops::Log::trace() << "Increment(ORCA)::print starting" << std::endl;

  os << std::endl << " Increment valid at time: " << validTime() << std::endl;
  os << std::string(4, ' ') << vars_ <<  std::endl;
  os << std::string(4, ' ') << "atlas field:" << std::endl;
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::tie(norm, min, max) = stats(fieldName);
    os << std::string(8, ' ') << fieldName << " norm: " << std::setprecision(5)
       << norm;
    os << " min: " << min << " max: " << max << std::endl;
  }

  oops::Log::trace() << "Increment(ORCA)::print done" << std::endl;
}

std::tuple<double, double, double> Increment::stats(const std::string & fieldName) const {
  double norm = 0;
  int valid_points = 0;
  double min=1e30;
  double max=-1e30;

  auto field_view = atlas::array::make_view<double, 2>(
      incrementFields_[fieldName]);
  oops::Log::trace() << "Increment(ORCA)::norm" << std::endl;
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  atlas::field::MissingValue mv(incrementFields()[fieldName]);
  bool has_mv = static_cast<bool>(mv);
  for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
    for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
      if (!ghost(j)) {
        if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
          if (field_view(j,k) > max) { max=field_view(j,k); }
          if (field_view(j,k) < min) { min=field_view(j,k); }
          norm += field_view(j, k)*field_view(j, k);
          ++valid_points;
          int x, y;                                       // DJL
          std::tie(x, y) = xypt(j);
          if (abs(x-xpt) < 2 && abs(y-ypt) < 2) {
            oops::Log::debug() << "x " << x << " y " << y << " increments field_view " << field_view(j,k) << std::endl; 
          }
        }
      }
    }
  }
  
  return std::make_tuple(sqrt(norm)/valid_points, min, max);
}

double Increment::norm() const {
  double norm = 0;
  double dnorm = 0;
  double min = 0;
  double max = 0;
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::tie(dnorm, min, max) = stats(fieldName);
    norm += dnorm;
  }
  return norm;
}


}  // namespace orcamodel
