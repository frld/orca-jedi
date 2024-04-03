#
# (C) British Crown Copyright 2024 Met Office
#

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# nemovar
find_package( nemovar REQUIRED )
include_directories( ${NEMOVAR_INCLUDE_DIRS} )

#nemovar/RunNV.h
#nemovar/RunNV.cc


list( APPEND oops_orcamodel_src_files
covariance/Covariance.h
covariance/CovarianceParameters.h
covariance/Covariance.cc
nemovar/ErrorCovarianceNV.h
nemovar/ErrorCovarianceNV.cc
nemovar/ErrorCovarianceNV.interface.F90
nemovar/GeometryNV.cc
nemovar/GeometryNV.h
nemovar/GeometryNV.interface.F90
nemovar/FortranNV.h
nemovar/VariablesNV.cc
nemovar/VariablesNV.h
nemovar/VariablesNV.interface.F90
nemovar/StateNV.cc
nemovar/StateNV.h
nemovar/FieldsNV.interface.F90
nemovar/FieldsTLAD.interface.F90
nemovar/FieldsNV.cc
nemovar/FieldsNV.h
nemovar/ModelNV.h
nemovar/ModelNV.cc
nemovar/IncrementNV.h
nemovar/IncrementNV.cc
nemovar/nv_init_wrapper.F90
nemovar/pm_link_mod.F90
nemovar/pm_linked_list_mod.F90
geometry/GeometryParameters.h
geometry/Geometry.h
geometry/Geometry.cc
interpolator/InterpolatorParameters.h
interpolator/Interpolator.h
interpolator/Interpolator.cc
increment/Increment.cc
increment/Increment.h
model/ModelBiasCovariance.h
model/ModelBias.h
model/ModelBiasIncrement.h
model/Model.h
run/Run.h
run/Run.cc
state/State.h
state/State.cc
state/StateIOUtils.h
state/StateIOUtils.cc
nemo_io/NemoFieldReader.h
nemo_io/NemoFieldReader.cc
nemo_io/NemoFieldWriter.h
nemo_io/NemoFieldWriter.cc
nemo_io/OrcaIndex.h
nemo_io/ReadServer.h
nemo_io/ReadServer.cc
nemo_io/WriteServer.h
nemo_io/WriteServer.cc
utilities/OrcaModelTraits.h
utilities/Types.h
utilities/Types.cc
variablechanges/VariableChangeParameters.h
variablechanges/VariableChange.h
variablechanges/LinearVariableChange.h
variablechanges/LinearVariableChange.cc
variablechanges/LinearVariableChangeParameters.h
)



ecbuild_add_library( TARGET   orcamodel
                     SOURCES  ${oops_orcamodel_src_files}
                     PUBLIC_LIBS  oops saber ufo atlas atlas-orca eckit netcdf nemovar
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
                  #target_link_libraries( orcamodel PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries( orcamodel PUBLIC NetCDF::NetCDF_CXX)
target_include_directories(orcamodel PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
#target_include_directories(orcamodel PUBLIC $<BUILD_INTERFACE:/data/users/frld/jedi/build-mo-bundle/nemovar/module> )
#target_include_directories(orcamodel PUBLIC $<BUILD_INTERFACE:${NEMOVAR_INCLUDE_DIRS}> )
target_include_directories(orcamodel PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/nemovar/module> )

