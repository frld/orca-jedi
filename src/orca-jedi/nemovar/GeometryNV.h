/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef NV_GEOMETRYNV_H
#define NV_GEOMETRYNV_H

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/memory/NonCopyable.h"


#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "orca-jedi/nemovar/FortranNV.h"

namespace nv {

// -----------------------------------------------------------------------------

// GeometryNV handles geometry for nemo model.

class GeometryF90 :
//        private eckit::NonCopyable,
        private util::ObjectCounter<GeometryF90> {
public:
    static const std::string classname() {return "nv::GeometryF90";}

// Constructor, destructor
//    explicit GeometryF90(const eckit::Configuration & conf);
    GeometryF90(const eckit::Configuration & conf);
    ~GeometryF90();

// Other methods
    void print(std::ostream & os) const {
        os << "GeometryNV::print not yet implemented." << std::endl;
    }

// Access to data
    F90geom& toFortran() {return keyGeom_;}
    const F90geom& toFortran() const {return keyGeom_;}
private:
    F90geom keyGeom_;
};

// -----------------------------------------------------------------------------

class GeometryNV :
        public util::Printable,
        private util::ObjectCounter<GeometryNV> {
public:
    static const std::string classname() {return "nv::GeometryNV";}

// Constructor, destructor
    explicit GeometryNV(const eckit::Configuration & conf);
    GeometryNV(const GeometryNV & other);
    ~GeometryNV();

// Access to data
//    F90geom& toFortran() {return geom_->toFortran();}
//    const F90geom& toFortran() const {return geom_->toFortran();}
    F90geom& toFortran() {return geom_.toFortran();}
    const F90geom& toFortran() const {return geom_.toFortran();}

    int geomnv_int;     // DJL

private:
    void print(std::ostream & os) const;

//    std::shared_ptr<GeometryF90> geom_;
    GeometryF90 geom_;
};

// -----------------------------------------------------------------------------

}  // namespace nv

#endif // NV_GEOMETRYNV_H
