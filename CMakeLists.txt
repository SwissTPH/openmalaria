# Root CMake configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
cmake_minimum_required (VERSION 2.4)
project (OpenMalaria CXX)

# The FindXXX.cmake files included are not quite the same as standard ones.
# They abort if not found, and look in a few extra places (e.g. ${CMAKE_SOURCE_DIR}/lib).
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

# Statically link libgcc; isn't going to work when other C++ libraries are dynamically linked
IF( ${CMAKE_COMPILER_IS_GNUCXX} )
  option (STATICALLY_LINK "Attempt to link libraries, including stdc++, statically.")
  if (STATICALLY_LINK)
    link_directories (${CMAKE_SOURCE_DIR}/lib)
    ADD_DEFINITIONS(-static-libgcc)
    SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libgcc")
  endif (STATICALLY_LINK)
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

option (BUILD_SHARED_LIBS "Link xsdcxx, WithinHostModel, etc. libraries dynamically instead of statically (almost certainly not wanted)." OFF)
MARK_AS_ADVANCED (BUILD_SHARED_LIBS)

find_package (XSD)
find_package (XERCESC)

add_subdirectory (xsdcxx)
add_subdirectory (model)

option(OMTEST_ENABLE "Enable 'make test' or Visual Studio build target to run tests" ON)
if (OMTEST_ENABLE)
  enable_testing()
  add_subdirectory (test)
endif (OMTEST_ENABLE)

option(OMTEST_UNIT "Enable lower-level unittests" OFF)
if (OMTEST_UNIT)
  add_subdirectory (unittest)
endif (OMTEST_UNIT)
