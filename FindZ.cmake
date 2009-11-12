# Finds the zlib library and include files
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines:
# Z_INCLUDE_DIRS
# Z_LIBRARIES

# list common paths here:
FIND_PATH(Z_INCLUDE_DIRS zconf.h
 PATHS
  ${CMAKE_SOURCE_DIR}/../zlib/include
  ${CMAKE_SOURCE_DIR}/../zlib
)

if (NOT Z_INCLUDE_DIRS)
  message (SEND_ERROR "Unable to find zlib includes (zconf.h)")
endif (NOT Z_INCLUDE_DIRS)
message (STATUS "Z_INCLUDE_DIRS is ${Z_INCLUDE_DIRS}")

find_library (Z_LIBRARIES NAMES z zlib zdll
  PATHS ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/../zlib/lib
        ${CMAKE_SOURCE_DIR}/../zlib/projects/visualc6/Win32_LIB_Release
)
if (NOT Z_LIBRARIES)
  message (SEND_ERROR "Unable to find zlib library")
endif (NOT Z_LIBRARIES)

MARK_AS_ADVANCED(
  Z_INCLUDE_DIRS
  Z_LIBRARIES
)
