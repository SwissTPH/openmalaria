# Finds the BOINC library and includes
# Copyright © 2005-2013 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines:
# BOINC_LIBRARIES
# BOINC_INCLUDE_DIRS
# (plus some cache variables)

# Note: on Unix, this will search for BOINC in the standard paths
# (/usr/include, /usr/local/include for headers, etc.).
# You can install BOINC here by (1) downloading from
# http://boinc.berkeley.edu/trac/wiki/SourceCodeGit
# (2) ./_autosetup
# (3) ./configure --disable-client --disable-manager --disable-server
# (4) make && sudo make install

# On all platforms this will also look for BOINC in ../boinc
# (in this case it doesn't need to be "installed").

unset (BOINC_BASE_PATH CACHE)	# don't use this now
if (MSVC)
else ()
  # unset things not needed elsewhere:
  unset (BOINC_LIB_DBG CACHE)
  unset (BOINC_API_LIB_DBG CACHE)
endif()

# Look for headers in the standard locations:
FIND_PATH(BOINC_INCLUDE_DIR boinc/boinc_api.h api/boinc_api.h
  PATHS ${CMAKE_SOURCE_DIR}/../boinc
  DOC "Either an install path like /usr/local/include which contains boinc/boinc_api.h or a BOINC source directory containing api/boinc_api.h"
)
if (BOINC_INCLUDE_DIR)
  if (EXISTS "${BOINC_INCLUDE_DIR}/api/boinc_api.h")
    # Assume source directory
    # Include all three of these directories:
    set (BOINC_INCLUDE_DIRS ${BOINC_INCLUDE_DIR}
      ${BOINC_INCLUDE_DIR}/api ${BOINC_INCLUDE_DIR}/lib)
  else()
    # Assume install directory; need to include the subdirectory:
    set (BOINC_INCLUDE_DIRS "${BOINC_INCLUDE_DIR}/boinc")
  endif ()
else()
  message (FATAL_ERROR "Not found: boinc/boinc_api.h or api/boinc_api.h (make sure BOINC is installed or has source code available and adjust BOINC_INCLUDE_DIR)")
endif()

# We now co-opt BOINC_INCLUDE_DIR to try to find libraries.
# On Windows it is only useful if it points to the BOINC source directory.

if (OM_STATICALLY_LINK) # unix static
  set (BOINC_LIB_NAME "libboinc.a")
  set (BOINC_API_NAME "libboinc_api.a")
elseif (OM_USE_LIBCMT) # windows static
  set (BOINC_LIB_NAME "libboinc_staticcrt.lib")
  set (BOINC_API_NAME "libboincapi_staticcrt.lib")
else ()
  # generic
  set (BOINC_LIB_NAME boinc boinc_staticcrt)
  set (BOINC_API_NAME boinc_api boincapi boincapi_staticcrt)
endif ()

# Find libboinc and libboinc_api/libboincapi
# Note that we always prefer static libraries when building for BOINC.

find_library (BOINC_LIB_OPT ${BOINC_LIB_NAME}
  PATHS
    ${BOINC_INCLUDE_DIR}/win_build/Build/Win32/Release
    ${BOINC_INCLUDE_DIR}/lib
    ${BOINC_INCLUDE_DIR}/../lib
)
find_library (BOINC_API_LIB_OPT ${BOINC_API_NAME}
  PATHS
    ${BOINC_INCLUDE_DIR}/win_build/Build/Win32/Release
    ${BOINC_INCLUDE_DIR}/api
    ${BOINC_INCLUDE_DIR}/../lib
)

if (MSVC)
  # debug libraries may also be needed on Windows
  find_library (BOINC_LIB_DBG ${BOINC_LIB_NAME}
    PATHS
      ${BOINC_INCLUDE_DIR}/win_build/Build/Win32/Debug
  )
  find_library (BOINC_API_LIB_DBG ${BOINC_API_NAME}
    PATHS
      ${BOINC_INCLUDE_DIR}/win_build/Build/Win32/Debug
  )
endif()

if (BOINC_LIB_OPT)
  if (BOINC_LIB_DBG)
    set (BOINC_LIB optimized ${BOINC_LIB_OPT} debug ${BOINC_LIB_DBG})
  else ()
    set (BOINC_LIB ${BOINC_LIB_OPT})
  endif ()
else ()
  if (BOINC_LIB_DBG)
    set (BOINC_LIB ${BOINC_LIB_DBG})
  else ()
    message (FATAL_ERROR "Unable to find boinc library")
  endif ()
endif ()
if (BOINC_API_LIB_OPT)
  if (BOINC_API_LIB_DBG)
    set (BOINC_API_LIB optimized ${BOINC_API_LIB_OPT} debug ${BOINC_API_LIB_DBG})
  else ()
    set (BOINC_API_LIB ${BOINC_API_LIB_OPT})
  endif ()
else ()
  if (BOINC_API_LIB_DBG)
    set (BOINC_API_LIB ${BOINC_API_LIB_DBG})
  else ()
    message (FATAL_ERROR "Unable to find boinc_api/libboincapi library")
  endif ()
endif ()

set (BOINC_LIBRARIES ${BOINC_API_LIB} ${BOINC_LIB})

mark_as_advanced (BOINC_API_LIB BOINC_LIB
  BOINC_API_LIB_OPT BOINC_API_LIB_DBG
  BOINC_LIB_OPT BOINC_LIB_DBG
  BOINC_INCLUDE_DIR)
