# Finds the xerces-c library
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines:
# XERCESC_LIBRARIES

FIND_PATH(XERCESC_INCLUDE_DIRS xercesc/dom/DOMNode.hpp
  PATHS "[HKEY_CURRENT_USER\\software\\xsd\\include]"
  "[HKEY_CURRENT_USER]\\xsd\\include]"
  $ENV{XSDDIR}/include
  /usr/local/include
  /usr/include
  "C:/Program Files/CodeSynthesis XSD 3.2/include"
  "D:/Program Files/CodeSynthesis XSD 3.2/include"
  "C:/Program Files/CodeSynthesis XSD 3.3/include"
  "C:/Program Files (x86)/CodeSynthesis XSD 3.3/include"
  ${CMAKE_SOURCE_DIR}/../xsd/libxsd
)

set (XERCESC_LIB_PATHS
  PATHS ${CMAKE_SOURCE_DIR}/lib
  "C:/Program Files/CodeSynthesis XSD 3.2/lib/vc-8.0"
  "D:/Program Files/CodeSynthesis XSD 3.2/lib/vc-8.0"
  "C:/Program Files/CodeSynthesis XSD 3.3/lib/vc-8.0"
  "C:/Program Files (x86)/CodeSynthesis XSD 3.3/lib/vc-8.0"
)
find_library (XERCESC_LIB_OPT NAMES xerces-c xerces-c_3 xerces-c_2 xerces-c_static xerces-c_static_3 xerces-c_static_2
  ${XERCESC_LIB_PATHS}
)
find_library (XERCESC_LIB_DBG NAMES xerces-c_D xerces-c_3D xerces-c_2D xerces-c_static_D xerces-c_static_3D xerces-c_static_2D
  ${XERCESC_LIB_PATHS}
)

if (XERCESC_LIB_OPT AND XERCESC_LIB_DBG)
  set (XERCESC_LIB optimized ${XERCESC_LIB_OPT} debug ${XERCESC_LIB_DBG})
elseif (XERCESC_LIB_OPT)
  set (XERCESC_LIB ${XERCESC_LIB_OPT})
elseif (XERCESC_LIB_DBG)
  set (XERCESC_LIB ${XERCESC_LIB_DBG})
else ()
  message (FATAL_ERROR "Unable to find xerces-c library")
endif ()
if (NOT XERCESC_INCLUDE_DIRS)
  message (FATAL_ERROR "Unable to find xsd include files (xsd/cxx/parser/elements.hxx)")
endif (NOT XERCESC_INCLUDE_DIRS)

SET (XERCESC_LIBRARIES ${XERCESC_LIB})

MARK_AS_ADVANCED(
  XERCESC_LIB
  XERCESC_LIBRARIES
  XERCESC_INCLUDE_DIRS
  XERCESC_LIB_DBG
  XERCESC_LIB_OPT
)
