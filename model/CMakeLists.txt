# CMake configuration for all code in this dir and below
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

# Don't use aux_source_directory on . because we don't want to compile openMalaria.cpp in to the lib.
set (Model_CPP
  WeibullDecayedValue.cpp population.cpp timer.cpp InfectionIncidenceModel.cpp
  gzstream.C simulation.cpp summary.cpp BoincWrapper.cpp global.cpp
  intervention.cpp proteome.cpp human.cpp
)
aux_source_directory (Clinical Clinical_CPP)
aux_source_directory (Drug Drug_CPP)
aux_source_directory (Pathogenesis Pathogenesis_CPP)
aux_source_directory (Transmission Transmission_CPP)
aux_source_directory (WithinHost WithinHost_CPP)
aux_source_directory (util util_CPP)

# Headers - only included so they show up in the Visual Studio project:
aux_source_directory (${CMAKE_SOURCE_DIR}/include Model_H)
aux_source_directory (${CMAKE_SOURCE_DIR}/include/Clinical Clinical_H)
aux_source_directory (${CMAKE_SOURCE_DIR}/include/Drug Drug_H)
aux_source_directory (${CMAKE_SOURCE_DIR}/include/Pathogenesis Pathogenesis_H)
aux_source_directory (${CMAKE_SOURCE_DIR}/include/Transmission Transmission_H)
aux_source_directory (${CMAKE_SOURCE_DIR}/include/WithinHost WithinHost_H)
aux_source_directory (${CMAKE_SOURCE_DIR}/include/util util_H)


include_directories (SYSTEM
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
)

add_library (model
  ${Model_CPP}
  ${Clinical_CPP}
  ${Drug_CPP}
  ${Pathogenesis_CPP}
  ${Transmission_CPP}
  ${WithinHost_CPP}
  ${util_CPP}
  ${Model_H}
  ${Clinical_H}
  ${Drug_H}
  ${Pathogenesis_H}
  ${Transmission_H}
  ${WithinHost_H}
)

if (MSVC)
  set_target_properties (model PROPERTIES LINK_FLAGS "${OM_LINK_FLAGS}")
endif (MSVC)


#for Visual Studio file groupings:
source_group (model\\Clinical FILES ${Clinical_CPP})
source_group (model\\Drug FILES ${Drug_CPP})
source_group (model\\Pathogenesis FILES ${Pathogenesis_CPP})
source_group (model\\Transmission FILES ${Transmission_CPP})
source_group (model\\WithinHost FILES ${WithinHost_CPP})
source_group (model\\util FILES ${util_CPP})

source_group (include FILES ${Model_H})
source_group (include\\Clinical FILES ${Clinical_H})
source_group (include\\Drug FILES ${Drug_H})
source_group (include\\Pathogenesis FILES ${Pathogenesis_H})
source_group (include\\Transmission FILES ${Transmission_H})
source_group (include\\WithinHost FILES ${WithinHost_H})
source_group (include\\util FILES ${util_H})
