/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "GeometryNV.h"

#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace nv {
// -----------------------------------------------------------------------------

GeometryF90::GeometryF90(const eckit::Configuration & conf) {
    oops::Log::debug() << "GeometryF90 initialise" << std::endl;
    const eckit::Configuration * configc = &conf;
    oops::Log::debug() << "GeometryF90 initialise2" << std::endl;
    nv_geo_setup_f90(keyGeom_, &configc);
    ASSERT(keyGeom_ != 0);
    oops::Log::trace() << "GeometryF90::GeometryF90 ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

GeometryF90::~GeometryF90() {
    oops::Log::trace() << "GeometryF90::GeometryF90 dtor start (dummy)" << std::endl;     // DJL
//    nv_geo_delete_f90(keyGeom_);   // DJL skip
    oops::Log::trace() << "GeometryF90::GeometryF90 dtor done" << std::endl;
}

// -----------------------------------------------------------------------------


GeometryNV::GeometryNV(const eckit::Configuration & conf)
//    : geom_(new GeometryF90(conf))
    : geom_(GeometryF90(conf))
{
    geomnv_int = 2;
    oops::Log::trace() << "GeometryNV::GeometryNV ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

GeometryNV::GeometryNV(const GeometryNV & other)
    : geom_(other.geom_)
{
    geomnv_int = 3;
    oops::Log::trace() << "GeometryNV::GeometryNV copy ctor done" << std::endl;
}

// -----------------------------------------------------------------------------

GeometryNV::~GeometryNV() {
    oops::Log::trace() << "GeometryNV::GeometryNV dtor done" << std::endl;
}

// -----------------------------------------------------------------------------

void GeometryNV::print(std::ostream & os) const {
//    geom_->print(os);
     geom_.print(os);
}

// -----------------------------------------------------------------------------

}  // namespace nv
