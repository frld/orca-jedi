
#include <sstream>

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"
#include "atlas/field.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"

namespace orcamodel {

void readFieldsFromFile(
  const eckit::Configuration & conf,
  const Geometry & geom,
  atlas::FieldSet & fs,
  std::string variable_type) {

    oops::Log::trace() << "orcamodel::readFieldsFromFile:: start "
                       << std::endl;

    // Open Nemo Feedback file
    std::string nemo_file_name = "";

    oops::Log::debug() << "orcamodel::readFieldsFromFile:: configuration "
                       << conf
                       << std::endl;
    if (variable_type == "background" ) {
      nemo_file_name = conf.getString("nemo field file");
    }
    if (variable_type == "background variance" ) {
    oops::Log::debug() << "orcamodel::readFieldsFromFile:: background variance has field file "
                       << conf.has("variance field file") << std::endl;
      nemo_file_name = conf.getString("variance field file");
    }
    oops::Log::debug() << "orcamodel::readFieldsFromFile:: nemo field file " << nemo_file_name
                       << std::endl;

    auto nemo_field_path = eckit::PathName(nemo_file_name);
    oops::Log::debug() << "orcamodel::readFieldsFromFile:: nemo_field_path "
                       << nemo_field_path << std::endl;
    NemoFieldReader nemo_file(nemo_field_path);

    // Read fields from Nemo field file
    // field names in the atlas fieldset are assumed to match their names in
    // the field file, perhaps change this or update "state variables" to use
    // these names, or to be a dictionary between the standard anmes and the
    // NEMO names
    for (atlas::Field field : fs) {
      std::string fieldName = field.name();
      oops::Log::debug() << "orcamodel::readFieldsFromFile:: field name = " << fieldName
                         << std::endl;

      oops::Log::debug() << "orcamodel::readFieldsFromFile:: "
                         << "geom.variable_in_variable_type(\"" << fieldName << "\", \"" << variable_type << "\") "
                         << geom.variable_in_variable_type(fieldName, variable_type) << std::endl;
      if (geom.variable_in_variable_type(fieldName, variable_type)) {

        auto field_view = atlas::array::make_view<double, 1>( field );
        nemo_file.read_surf_var(fieldName, field_view);
      }
    }

    oops::Log::trace() << "orcamodel::readFieldsFromFile:: readFieldsFromFile done "
                       << std::endl;
};

}
