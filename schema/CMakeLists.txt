# CMake configuration for all code in this dir and below
# Copyright © 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)


foreach (XSD_NAME scenario demography)
  set (XSD_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${XSD_NAME}.xsd)
  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.cpp ${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.h
    COMMAND ${XSD_EXECUTABLE} cxx-tree
        --type-naming ucc --function-naming java
        --namespace-map =scnXml
        --generate-serialization
        --generate-doxygen --generate-intellisense
        --hxx-suffix .h --cxx-suffix .cpp
        ${XSD_FILE}
    MAIN_DEPENDENCY ${XSD_FILE}
    COMMENT "Compiling ${XSD_FILE}"
    VERBATIM
  )
endforeach (XSD_NAME)

include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_BINARY_DIR}
)

set (SCHEMA_CPP
  ${CMAKE_CURRENT_BINARY_DIR}/scenario.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/demography.cpp
)
# Headers - only included so they show up in IDEs:
file (GLOB_RECURSE SCHEMA_H "${CMAKE_CURRENT_BINARY_DIR}/*.h")
add_library (schema
  ${SCHEMA_CPP}
  ${SCHEMA_H}
)
