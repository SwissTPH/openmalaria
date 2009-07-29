# Root CMake configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
cmake_minimum_required (VERSION 2.4)
project (OpenMalaria CXX)

# The FindXXX.cmake files included are not quite the same as standard ones.
# They abort if not found, and look in a few extra places (e.g. ${CMAKE_SOURCE_DIR}/lib).
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})


# -----  Compile options  -----

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

if (MSVC)
  option (OM_USE_LIBCMT "Use LIBCMT.lib instead of MSVCRDT[D].lib when linking. This removes a dependency on a .net library, but usually requires manually building dependency libraries." OFF)
  # set a temporary variable, since cache variables don't get rewritten
  if (OM_USE_LIBCMT)
    set (OM_LINK_FLAGS "/NODEFAULTLIB:MSVCRDT /NODEFAULTLIB:MSVCRDTD")
  else (OM_USE_LIBCMT)
    set (OM_LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
  endif (OM_USE_LIBCMT)
endif (MSVC)


# -----  Find dependencies  -----

find_package (XSD)
find_package (XERCESC)
find_package (GSL)
find_package (Z)

option (BOINC_INTEGRATION "Set to enable integration of the BOINC client (otherwise it's not used).")
if (BOINC_INTEGRATION)
  find_library (BOINC_LIB boinc PATHS ${CMAKE_SOURCE_DIR}/../boinc/lib)
  find_library (BOINC_API_LIB boinc_api PATHS ${CMAKE_SOURCE_DIR}/../boinc/api)
  include_directories (${CMAKE_SOURCE_DIR}/../boinc ${CMAKE_SOURCE_DIR}/../boinc/api ${CMAKE_SOURCE_DIR}/../boinc/lib)
  if (NOT BOINC_LIB OR NOT BOINC_API_LIB)
    message (FATAL_ERROR "Unable to find a boinc library (boinc or boinc_api)")
  endif (NOT BOINC_LIB OR NOT BOINC_API_LIB)
  set (BOINC_LIBRARIES ${BOINC_API_LIB} ${BOINC_LIB})
  mark_as_advanced (BOINC_API_LIB BOINC_LIB)
else (BOINC_INTEGRATION)
  add_definitions (-DWITHOUT_BOINC)
endif (BOINC_INTEGRATION)


# -----  Compile code  -----

add_subdirectory (xsdcxx)
add_subdirectory (model)


# -----  generate openMalaria  -----

if (UNIX)
  find_library (PTHREAD_LIBRARIES pthread PATHS ${CMAKE_SOURCE_DIR}/lib /usr/lib /usr/local/lib)
  if (NOT PTHREAD_LIBRARIES)
    message (SEND_ERROR "PThread library not found")
  endif (NOT PTHREAD_LIBRARIES)
endif (UNIX)
mark_as_advanced (BOINC_LIBRARIES PTHREAD_LIBRARIES)

include_directories (SYSTEM
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
)

link_directories (${CMAKE_BINARY_DIR}/xsdcxx)

add_executable (openMalaria model/openMalaria.cpp)

target_link_libraries (openMalaria
  model
  xsdcxx
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${BOINC_LIBRARIES}
)

if (MSVC)
  set_target_properties (openMalaria PROPERTIES LINK_FLAGS "${OM_LINK_FLAGS}")
endif (MSVC)

# -----  OMTEST - black-box & unit testing  -----

option(OMTEST_ENABLE "Enable 'make test' or Visual Studio build target to run black-box tests" ON)
if (OMTEST_ENABLE)
  enable_testing()
  add_subdirectory (test)
endif (OMTEST_ENABLE)

option(CXXTEST_ENABLE "Enable lower-level unittests using cxxtest" OFF)
if (CXXTEST_ENABLE)
  enable_testing()
  add_subdirectory (unittest)
endif (CXXTEST_ENABLE)
