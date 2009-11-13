# CMake configuration for all code in this dir and below
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

# Don't use aux_source_directory on . because we don't want to compile openMalaria.cpp in to the lib.
set (Model_CPP
  population.cpp InfectionIncidenceModel.cpp
  gzstream.C Simulation.cpp summary.cpp Global.cpp
  intervention.cpp proteome.cpp Human.cpp
  inputData.cpp NeonatalMortality.cpp
  
  ${CMAKE_CURRENT_BINARY_DIR}/scenario.cxx
  
  Clinical/ClinicalModel.cpp
  Clinical/EventScheduler.cpp
  Clinical/ImmediateOutcomes.cpp
  Clinical/OldCaseManagement.cpp
  Clinical/Episode.cpp
  
  Drug/DrugModel.cpp
  Drug/PkPdDrug.cpp
  Drug/drug.cpp
  
  Pathogenesis/PathogenesisModel.cpp
  Pathogenesis/Mueller.cpp
  Pathogenesis/Pyrogen.cpp
  Pathogenesis/Predet.cpp
  
  Transmission/TransmissionModel.cpp
  Transmission/NonVector.cpp
  Transmission/PerHostTransmission.cpp
  Transmission/Vector/VectorTransmission.cpp
  Transmission/Vector/VectorAnopheles.cpp
  Transmission/Vector/HostCategoryAnopheles.cpp
  
  WithinHost/WithinHostModel.cpp
  WithinHost/Infection.cpp
  WithinHost/Descriptive.cpp
  WithinHost/DescriptiveInfection.cpp
  WithinHost/DescriptiveIPT.cpp
  WithinHost/DescriptiveIPTInfection.cpp
  WithinHost/Dummy.cpp
  WithinHost/DummyInfection.cpp
  WithinHost/Empirical.cpp
  WithinHost/EmpiricalInfection.cpp
  
  util/BoincWrapper.cpp
  util/timer.cpp
  util/gsl.cpp
  util/vectors.cpp
  util/WeibullDecayedValue.cpp
)

# Headers - only included so they show up in IDEs:
file (GLOB Model_H "${CMAKE_SOURCE_DIR}/include/*.h")
file (GLOB SModel_H "${CMAKE_SOURCE_DIR}/include/*/*.h")
file (GLOB SSModel_H "${CMAKE_SOURCE_DIR}/include/*/*/*.h")

set (XSD_FILE ${CMAKE_CURRENT_SOURCE_DIR}/scenario.xsd)
add_custom_command (
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/scenario.cxx ${CMAKE_CURRENT_BINARY_DIR}/scenario.ixx ${CMAKE_CURRENT_BINARY_DIR}/scenario.hxx
  COMMAND ${XSD_EXECUTABLE} cxx-tree
      --type-naming ucc --function-naming java
      --namespace-map =scnXml
      --generate-serialization --generate-inline
      --generate-doxygen --generate-intellisense
      ${XSD_FILE}
  MAIN_DEPENDENCY ${XSD_FILE}
  COMMENT "Compiling ${XSD_FILE}"
  VERBATIM
)


include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}
)

add_library (model
  ${Model_CPP}
  ${Model_H}
  ${SModel_H}
  ${SSModel_H}
  ${CMAKE_CURRENT_BINARY_DIR}/scenario.ixx ${CMAKE_CURRENT_BINARY_DIR}/scenario.hxx
)

if (MSVC)
  set_target_properties (model PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)
