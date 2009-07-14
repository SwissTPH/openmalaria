# Finds the xerces-c library
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines:
# XERCESC_LIBRARIES

find_library (XERCESC_LIB NAMES xerces-c xerces-c_static xerces-c_static_2 xerces-c_static_3
  PATHS ${CMAKE_SOURCE_DIR}/lib
)
if (NOT XERCESC_LIB)
  message (FATAL_ERROR "Unable to find xerces-c library")
endif (NOT XERCESC_LIB)

SET (XERCESC_LIBRARIES ${XERCESC_LIB})

MARK_AS_ADVANCED(
  XERCESC_LIB
  XERCESC_LIBRARIES
)
