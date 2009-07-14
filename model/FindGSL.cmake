# list common paths here:
# The gsl-111 dir is because people have GSL 1.8 installed in C:/Program Files/GnuWin32/include
# which isn't supported and gsl-1.11 extracted to in gsl-111
FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_version.h
  PATHS ${CMAKE_SOURCE_DIR}/../gsl/include
  "C:/Program Files/GnuWin32/gsl-111/Binaries/gsl/include"
  "C:/Program Files/GnuWin32/include"
)

if (NOT GSL_INCLUDE_DIR)
  message (SEND_ERROR "Unable to find gsl/gsl_version.h")
endif (NOT GSL_INCLUDE_DIR)
message (STATUS "GSL_INCLUDE_DIR is ${GSL_INCLUDE_DIR}")

find_library (GSL_LIB gsl
  PATHS ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/../gsl/lib
  "C:/Program Files/GnuWin32/gsl-111/Binaries/gsl/lib"
  "C:/Program Files/GnuWin32/lib"
)
find_library (GSL_CBLAS_LIB NAMES gsl_cblas cblas
  PATHS ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/../gsl/lib
  "C:/Program Files/GnuWin32/gsl-111/Binaries/gsl/lib"
  "C:/Program Files/GnuWin32/lib"
)
if (NOT GSL_LIB OR NOT GSL_CBLAS_LIB)
  message (SEND_ERROR "Unable to find gsl library")
endif (NOT GSL_LIB OR NOT GSL_CBLAS_LIB)

MARK_AS_ADVANCED(
  GSL_INCLUDE_DIR
  GSL_LIB
  GSL_CBLAS_LIB
)
