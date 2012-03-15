# Finds the GSL library and includes
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines:
# GSL_LIBRARIES
# GSL_INCLUDE_DIRS

# old cache entry; no longer stored in cache:
unset (GSL_INCLUDE_DIRS CACHE)

# list common paths here:
# The gsl-111 dir is because people have GSL 1.8 installed in C:/Program Files/GnuWin32/include
# which isn't supported and gsl-1.11 extracted to in gsl-111
FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_version.h
 PATHS
  ${CMAKE_SOURCE_DIR}/../gsl/include
  ${CMAKE_SOURCE_DIR}/../gsl
  "C:/Program Files/GnuWin32/gsl-111/Binaries/gsl/include"
  "C:/Program Files/GnuWin32/include"
  "C:/Program Files/GnuWin32/gsl-111/gsl/win32/include"
)

if (NOT GSL_INCLUDE_DIR)
  message (FATAL_ERROR "Unable to find gsl/gsl_version.h")
endif (NOT GSL_INCLUDE_DIR)

FIND_PATH(GSL_INCLUDE_DIR2 gsl_sys.h
 PATHS
  ${GSL_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}/gsl
 NO_DEFAULT_PATH
)
set (GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
if (GSL_INCLUDE_DIR2)
  MARK_AS_ADVANCED(GSL_INCLUDE_DIR2)
  if (NOT ${GSL_INCLUDE_DIR} STREQUAL ${GSL_INCLUDE_DIR2})
    set (GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} ${GSL_INCLUDE_DIR2})
  endif (NOT ${GSL_INCLUDE_DIR} STREQUAL ${GSL_INCLUDE_DIR2})
endif (GSL_INCLUDE_DIR2)
message(STATUS "GSL include path(s): ${GSL_INCLUDE_DIRS}")

# gsl puts libraries in various places, with various suffixes. Copy the ones you want to use to ../gsl/lib and remove suffixes.
set (GSL_LIB_PATHS
  PATHS ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/../gsl/lib
  "C:/Program Files/GnuWin32/gsl-111/gsl/win32/lib"
  "C:/Program Files/GnuWin32/gsl-111/Binaries/gsl/lib"
  "C:/Program Files/GnuWin32/lib"
)
find_library (GSL_LIB_OPT gsl ${GSL_LIB_PATHS})
find_library (GSL_LIB_DBG gsl_d ${GSL_LIB_PATHS})
find_library (GSL_CBLAS_LIB_OPT NAMES gslcblas cblas ${GSL_LIB_PATHS})
find_library (GSL_CBLAS_LIB_DBG NAMES gslcblas_d cblas_d ${GSL_LIB_PATHS})
if (GSL_LIB_OPT AND GSL_LIB_DBG)
  set (GSL_LIB optimized ${GSL_LIB_OPT} debug ${GSL_LIB_DBG})
elseif (GSL_LIB_OPT)
  set (GSL_LIB ${GSL_LIB_OPT})
elseif (GSL_LIB_DBG)
  set (GSL_LIB ${GSL_LIB_DBG})
else (GSL_LIB_OPT AND GSL_LIB_DBG)
  message (FATAL_ERROR "Unable to find gsl library")
endif (GSL_LIB_OPT AND GSL_LIB_DBG)
if (GSL_CBLAS_LIB_OPT AND GSL_CBLAS_LIB_DBG)
  set (GSL_CBLAS_LIB optimized ${GSL_CBLAS_LIB_OPT} debug ${GSL_CBLAS_LIB_DBG})
elseif (GSL_CBLAS_LIB_OPT)
  set (GSL_CBLAS_LIB ${GSL_CBLAS_LIB_OPT})
elseif (GSL_CBLAS_LIB_DBG)
  set (GSL_CBLAS_LIB ${GSL_CBLAS_LIB_DBG})
else (GSL_CBLAS_LIB_OPT AND GSL_CBLAS_LIB_DBG)
  message (FATAL_ERROR "Unable to find gsl cblas library")
endif (GSL_CBLAS_LIB_OPT AND GSL_CBLAS_LIB_DBG)

SET (GSL_LIBRARIES ${GSL_LIB} ${GSL_CBLAS_LIB})

MARK_AS_ADVANCED(
  GSL_INCLUDE_DIR
  GSL_INCLUDE_DIRS
  GSL_LIB
  GSL_CBLAS_LIB
  GSL_LIB_OPT
  GSL_LIB_DBG
  GSL_CBLAS_LIB_OPT
  GSL_CBLAS_LIB_DBG
)
