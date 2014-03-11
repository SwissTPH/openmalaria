# Root CMake configuration for OpenMalaria
# Copyright © 2005-2013 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

cmake_minimum_required (VERSION 2.4)
# We still (try to) support cmake 2.4, but cmake 2.8 policies also work:
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.6)
    cmake_policy(VERSION 2.8)
elseif("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
    cmake_policy(VERSION 2.6)
endif("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.6)

project (OpenMalaria CXX)

# The FindXXX.cmake files included are not quite the same as standard ones.
# They abort if not found, and look in a few extra places (e.g. ${CMAKE_SOURCE_DIR}/lib).
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

# -----  Compile options  -----

# C++11 would be nice, but we have to compile for relatively old platforms
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")

# Statically link libgcc; isn't going to work when other C++ libraries are dynamically linked
IF( ${CMAKE_COMPILER_IS_GNUCXX} )
  # This doesn't show up in the cache, but it is passed to gcc:
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  
  option (OM_STATICALLY_LINK "Attempt to link libraries, including stdc++, statically.")
  if (OM_STATICALLY_LINK)
    link_directories (${CMAKE_SOURCE_DIR}/lib)
    # Note: -static-libstdc++ is only available since GCC 4.5
    # See: http://www.trilithium.com/johan/2005/06/static-libstdc/
    ADD_DEFINITIONS(-static-libgcc -static-libstdc++)
    SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libgcc -static-libstdc++")
  endif (OM_STATICALLY_LINK)
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

option (BUILD_SHARED_LIBS "Link xsdcxx, model, etc. libraries dynamically instead of statically (almost certainly not wanted)." OFF)
MARK_AS_ADVANCED (BUILD_SHARED_LIBS)

if (MSVC)
  # This is now the OM_STATICALLY_LINK for windows
  option (OM_USE_LIBCMT "Use LIBCMT.lib instead of MSVCRT[D].lib when linking, and link statically. This removes a dependency on a .net library, but usually requires manually building dependency libraries." OFF)
  # set a temporary variable, since cache variables don't get rewritten
  if (OM_USE_LIBCMT)
    #piece of code from http://stackoverflow.com/questions/1618927/cmake-microsoft-visual-studio-and-monolithic-runtimes
    #Note: doesn't update cmake cache (but still works)
    #We statically link to reduce dependancies
    foreach(flag_var CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      if(${flag_var} MATCHES "/MDd")
        string(REGEX REPLACE "/MDd" "/MTd" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MDd")
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
    
    set (OM_LINK_FLAGS "/NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD /NODEFAULTLIB:library")
    set (OM_COMPILE_FLAGS "/DXERCES_STATIC_LIBRARY")
  else (OM_USE_LIBCMT)
    set (OM_LINK_FLAGS "/NODEFAULTLIB:LIBCMT")
    set (OM_COMPILE_FLAGS)
  endif (OM_USE_LIBCMT)
  set (OM_COMPILE_FLAGS "${OM_COMPILE_FLAGS} /D_CRT_SECURE_NO_DEPRECATE")
endif (MSVC)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

# -----  Find dependencies  -----

find_package (XSD)
find_package (XERCESC)
find_package (GSL)
find_package (Z)


# Find Boost:

# suggested on windows: install boost such that paths are:
# C:\Program Files\boost\boost_VER\boost\config.hpp
# Add versions greater than 1.38 here:
set (Boost_ADDITIONAL_VERSIONS "1.39" "1.39.0" "1.40" "1.40.0" "1.41" "1.41.0" "1.41" "1.42.0"  "1.43.0" "1.44.0" "1.45.0")
#set (Boost_DEBUG)
find_package( Boost 1.37.0 )
if(Boost_FOUND)
    find_path (SPIRIT_INCLUDE_DIRS boost/spirit/include/version.hpp
	PATHS ${CMAKE_SOURCE_DIR}/../spirit ${Boost_INCLUDE_DIRS}
    )
    if (NOT SPIRIT_INCLUDE_DIRS)
	message (FATAL_ERROR "Not found: boost/spirit/include/version.hpp
Required: Spirit 2.1 or later. (You can either upgrade boost to version
1.41 or later or you can extract spirit headers to ${CMAKE_SOURCE_DIR}/../spirit (containing boost/spirit/include/version.hpp).
See: http://boost-spirit.com/home/articles/spirit2-1-release/using-spirit-v2-1-with-older-versions-of-boost-pre-v1-41/ )"
	)
    endif (NOT SPIRIT_INCLUDE_DIRS)
    
    # Extract SPIRIT_VERSION macro (code taken from FindBoost.cmake)
    # Read the whole file:
    SET(SPIRIT_VERSION 0)
    FILE(READ "${SPIRIT_INCLUDE_DIRS}/boost/spirit/include/version.hpp" _spirit_VERSION_HPP_CONTENTS)
    
    STRING(REGEX REPLACE ".*#define SPIRIT_VERSION 0x([0-9A-Fa-f]+).*" "\\1" SPIRIT_VERSION "${_spirit_VERSION_HPP_CONTENTS}")
    SET(SPIRIT_VERSION ${SPIRIT_VERSION} CACHE INTERNAL "The version number for boost/spirit (base 16)")
    
    # Warning: this is a hex number. The test will probably fail if any digits are A-F.
    # The logic I want is !(SPIRIT_VERSION >= 2010) which is _not_ the same as
    # (SPIRIT_VERSION LESS 2010) (because this is false when SPIRIT_VERSION is not a number).
    if (SPIRIT_VERSION GREATER 2010)
    elseif (SPIRIT_VERSION EQUAL 2010)
    else (SPIRIT_VERSION GREATER 2010)
        message (FATAL_ERROR "Spirit version 0x${SPIRIT_VERSION} found.
Required: Spirit 2.1 or later. (You can either upgrade boost to version
1.41 or later or you can extract spirit headers to ${CMAKE_SOURCE_DIR}/../spirit (containing boost/spirit/include/version.hpp).
See: http://boost-spirit.com/home/articles/spirit2-1-release/using-spirit-v2-1-with-older-versions-of-boost-pre-v1-41/ )"
        )
    endif (SPIRIT_VERSION GREATER 2010)
    message (STATUS "Spirit version: ${SPIRIT_VERSION} (from: ${SPIRIT_INCLUDE_DIRS})")
else(Boost_FOUND)
    message (FATAL_ERROR "Unable to find boost headers! Please install libboost-dev or boost from boost.org (preferably at least version 1.41)")
endif(Boost_FOUND)


# Find BOINC:

option (OM_BOINC_INTEGRATION "Set to enable integration of the BOINC client (otherwise it's not used).")
if (OM_BOINC_INTEGRATION)
  find_package (BOINC)
else (OM_BOINC_INTEGRATION)
  add_definitions (-DWITHOUT_BOINC)
endif (OM_BOINC_INTEGRATION)


# -----  Compile-time optional features  -----
# (must come before add_subdirectory (model))

option (OM_STREAM_VALIDATOR "Compile in StreamValidator (see include/util/StreamValidator.h for usage notes)" OFF)
if (OM_STREAM_VALIDATOR)
  add_definitions (-DOM_STREAM_VALIDATOR)
endif (OM_STREAM_VALIDATOR)


# -----  Compile code  -----

add_subdirectory (contrib)
add_subdirectory (schema)
add_subdirectory (model)
add_subdirectory (util/SchemaTranslator)
add_dependencies (model schema)

# -----  generate openMalaria  -----

if (UNIX)
  find_library (PTHREAD_LIBRARIES pthread PATHS ${CMAKE_SOURCE_DIR}/lib /usr/lib /usr/local/lib)
  if (NOT PTHREAD_LIBRARIES)
    message (SEND_ERROR "PThread library not found")
  endif (NOT PTHREAD_LIBRARIES)
endif (UNIX)
mark_as_advanced (BOINC_LIBRARIES PTHREAD_LIBRARIES)

message (STATUS "BOINC includes: ${BOINC_INCLUDE_DIRS}")
include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${SPIRIT_INCLUDE_DIRS}
  ${BOINC_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}
)

add_executable (openMalaria model/openMalaria.cpp)

target_link_libraries (openMalaria
  model
  schema
  contrib
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${BOINC_LIBRARIES}
  # these need to be mentioned explicitly with Clang:
  m
  stdc++
)

if (MSVC)
  set_target_properties (openMalaria PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)


# -----  OM_BOXTEST - black-box & unit testing  -----

option(OM_CXXTEST_ENABLE "Enable lower-level unittests using cxx (use 'make test' or Visual Studio build target)test" ON)
if (OM_CXXTEST_ENABLE)
  enable_testing()
  add_subdirectory (unittest)
endif (OM_CXXTEST_ENABLE)

option(OM_BOXTEST_ENABLE "Enable black-box testing of openMalaria (use 'make test' or Visual Studio build target)" ON)
if (OM_BOXTEST_ENABLE)
  enable_testing()
  add_subdirectory (test)
endif (OM_BOXTEST_ENABLE)
