# Root CMake configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
cmake_minimum_required (VERSION 2.4)
project (OpenMalaria CXX)

# The FindXXX.cmake files included are not quite the same as standard ones.
# They abort if not found, and look in a few extra places (e.g. ${CMAKE_SOURCE_DIR}/lib).
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

set (OM_LATEST_SCHEMA "scenario_13.xsd")

# -----  Compile options  -----

# Statically link libgcc; isn't going to work when other C++ libraries are dynamically linked
IF( ${CMAKE_COMPILER_IS_GNUCXX} )
  option (OM_STATICALLY_LINK "Attempt to link libraries, including stdc++, statically.")
  if (OM_STATICALLY_LINK)
    link_directories (${CMAKE_SOURCE_DIR}/lib)
    ADD_DEFINITIONS(-static-libgcc)
    SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libgcc")
  endif (OM_STATICALLY_LINK)
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

option (BUILD_SHARED_LIBS "Link xsdcxx, WithinHostModel, etc. libraries dynamically instead of statically (almost certainly not wanted)." OFF)
MARK_AS_ADVANCED (BUILD_SHARED_LIBS)

if (MSVC)
  # This is now the OM_STATICALLY_LINK for windows
  option (OM_USE_LIBCMT "Use LIBCMT.lib instead of MSVCRT[D].lib when linking, and link statically. This removes a dependency on a .net library, but usually requires manually building dependency libraries." OFF)
  # set a temporary variable, since cache variables don't get rewritten
  if (OM_USE_LIBCMT)
    set (OM_LINK_FLAGS "/NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD")
    set (OM_COMPILE_FLAGS "/DXERCES_STATIC_LIBRARY")
    else (OM_USE_LIBCMT)
    set (OM_LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
    set (OM_COMPILE_FLAGS)
  endif (OM_USE_LIBCMT)
endif (MSVC)


option (OM_VECTOR_OUTPUT "Write per-timestep vector statistics to vector.csv" OFF)
if (OM_VECTOR_OUTPUT)
    add_definitions (OM_CSV_REPORTING)
endif (OM_VECTOR_OUTPUT)


# -----  Find dependencies  -----

find_package (XSD)
find_package (XERCESC)
find_package (GSL)
find_package (Z)
# help for windows installer
set (BOOST_ROOT "C:\\Program Files\\boost\\boost_1_40")
find_package( Boost 1.36.0 )
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
  message (FATAL_ERROR "Unable to find boost headers! Please install libboost-dev or boost from boost.org")
endif(Boost_FOUND)

option (OM_BOINC_INTEGRATION "Set to enable integration of the BOINC client (otherwise it's not used).")
if (OM_BOINC_INTEGRATION)
  if (OM_USE_LIBCMT) # pick only the right libs
    set (BOINC_LIB_POSTFIX "_staticcrt")
  else (OM_USE_LIBCMT)
    set (BOINC_LIB_POSTFIX "")
  endif (OM_USE_LIBCMT)
  message (STATUS "Looking for: boinc${BOINC_LIB_POSTFIX}")
  message (STATUS "in: ${CMAKE_SOURCE_DIR}/../boinc/win_build/Build/Win32/Release")
  find_library (BOINC_LIB boinc${BOINC_LIB_POSTFIX} libboinc${BOINC_LIB_POSTFIX}
    PATHS ${CMAKE_SOURCE_DIR}/../boinc/lib
    ${CMAKE_SOURCE_DIR}/../boinc/win_build/Build/Win32/Release
    ${CMAKE_SOURCE_DIR}/../boinc/win_build/Build/Win32/Debug
  )
  find_library (BOINC_API_LIB boinc_api${BOINC_LIB_POSTFIX} libboincapi${BOINC_LIB_POSTFIX}
    PATHS ${CMAKE_SOURCE_DIR}/../boinc/api
    ${CMAKE_SOURCE_DIR}/../boinc/win_build/Build/Win32/Release
    ${CMAKE_SOURCE_DIR}/../boinc/win_build/Build/Win32/Debug
  )
  include_directories (${CMAKE_SOURCE_DIR}/../boinc ${CMAKE_SOURCE_DIR}/../boinc/api ${CMAKE_SOURCE_DIR}/../boinc/lib)
  if (NOT BOINC_LIB OR NOT BOINC_API_LIB)
    message (FATAL_ERROR "Unable to find a boinc library (boinc or boinc_api)")
  endif (NOT BOINC_LIB OR NOT BOINC_API_LIB)
  set (BOINC_LIBRARIES ${BOINC_API_LIB} ${BOINC_LIB})
  mark_as_advanced (BOINC_API_LIB BOINC_LIB)
else (OM_BOINC_INTEGRATION)
  add_definitions (-DWITHOUT_BOINC)
endif (OM_BOINC_INTEGRATION)


# -----  Compile code  -----

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
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/model
)
link_directories (${CMAKE_BINARY_DIR}/xsdcxx)

add_executable (openMalaria model/openMalaria.cpp)

target_link_libraries (openMalaria
  model
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${BOINC_LIBRARIES}
)

if (MSVC)
  set_target_properties (openMalaria PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)

# -----  OM_BOXTEST - black-box & unit testing  -----

option(OM_BOXTEST_ENABLE "Enable black-box testing of openMalaria (use 'make test' or Visual Studio build target)" ON)
if (OM_BOXTEST_ENABLE)
  enable_testing()
  add_subdirectory (test)
endif (OM_BOXTEST_ENABLE)

option(OM_CXXTEST_ENABLE "Enable lower-level unittests using cxx (use 'make test' or Visual Studio build target)test" ON)
if (OM_CXXTEST_ENABLE)
  enable_testing()
  add_subdirectory (unittest)
endif (OM_CXXTEST_ENABLE)
