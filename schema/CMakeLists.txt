# CMake configuration for all code in this dir and below
# Copyright © 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

set (SCHEMA_NAMES
  scenario
  demography
  monitoring
  interventions
  healthSystem
  entomology
  pharmacology
  vivax
  util
)
set (SCHEMA_CPP "")
set (SCHEMA_H "")
set (SCHEMA_XSD "")
foreach (XSD_NAME ${SCHEMA_NAMES})
  list (APPEND SCHEMA_CPP "${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.cpp")
  list (APPEND SCHEMA_H "${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.h")
  list (APPEND SCHEMA_XSD "${XSD_NAME}.xsd")
  
  set (XSD_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${XSD_NAME}.xsd)
  add_custom_command (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.cpp ${CMAKE_CURRENT_BINARY_DIR}/${XSD_NAME}.h
    COMMAND ${XSD_EXECUTABLE} cxx-tree
        --std c++11
        --type-naming ucc --function-naming java
        --namespace-map http://openmalaria.org/schema/scenario_45=scnXml
#         --generate-serialization
        --generate-doxygen
        --generate-intellisense
        --hxx-suffix .h --cxx-suffix .cpp
        ${XSD_FILE}
    MAIN_DEPENDENCY ${XSD_FILE}
    COMMENT "Compiling ${XSD_FILE}"
    VERBATIM
  )
endforeach (XSD_NAME)

set( INLINED_XSD ${CMAKE_CURRENT_BINARY_DIR}/scenario_current.xsd )
add_custom_command (OUTPUT ${INLINED_XSD}
  DEPENDS ${SCHEMA_XSD}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/util/inlineSchema.py
  scenario.xsd ${INLINED_XSD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Inlining scenario.xsd and dependencies into ${INLINED_XSD}"
  VERBATIM
)
add_custom_target (inlined_xsd ALL DEPENDS ${INLINED_XSD})

include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_BINARY_DIR}
)

add_library (schema
  ${SCHEMA_CPP}
  ${SCHEMA_H}
)
