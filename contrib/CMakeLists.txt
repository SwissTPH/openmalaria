# CMake configuration for all code in this dir and below
# Copyright © 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

set (Contrib_CPP
    gzstream/gzstream.C
)

include_directories (SYSTEM
  ${Z_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_SOURCE_DIR}/contrib/gzstream
  ${CMAKE_SOURCE_DIR}/contrib/cxxtest
  ${CMAKE_SOURCE_DIR}/contrib/FastDelegate_src
  ${CMAKE_SOURCE_DIR}/contrib/floating_point_utilities_v3
)

add_library (contrib
  ${Contrib_CPP}
  ${Contrib_H}
)

if (MSVC)
  set_target_properties (contrib PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)
