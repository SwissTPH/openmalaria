# CMake configuration for all code in this dir and below
# Copyright © 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

add_subdirectory (../../util)

set (Infection_CPP
  Infection.cpp
  DescriptiveInfection.cpp
  DescriptiveIPTInfection.cpp
  DummyInfection.cpp
  EmpiricalInfection.cpp
  MolineauxInfection.cpp
)

# Headers - only included so they show up in IDEs:
file (GLOB_RECURSE Infection_H "${CMAKE_SOURCE_DIR}/include/WithinHost/Infection/*.h")

include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${SPIRIT_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/model
)

add_library (Infection
  ${Infection_CPP}
  ${Infection_H}
  ${CMAKE_BUILD_DIR}/util/scenario.hxx
)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_BINARY_DIR}/util/scenario.hxx PROPERTIES GENERATED 1)
add_dependencies (Infection util_scenario_xsd)

if (MSVC)
  set_target_properties (Infection PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)
