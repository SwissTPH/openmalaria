# CMake configuration for XML-reading code
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

set (XSD_FILE ${CMAKE_SOURCE_DIR}/xsdcxx/scenario.xsd)
add_custom_command (OUTPUT ${CMAKE_BINARY_DIR}/xsdcxx/scenario.cxx ${CMAKE_BINARY_DIR}/xsdcxx/scenario.ixx ${CMAKE_BINARY_DIR}/xsdcxx/scenario.hxx
    COMMAND ${XSD_EXECUTABLE} cxx-tree
      --type-naming ucc --function-naming java
      --namespace-map =scnXml
      --generate-serialization --generate-inline
      --generate-doxygen --generate-intellisense
      ${XSD_FILE}
    MAIN_DEPENDENCY ${XSD_FILE}
    COMMENT "Compiling ${XSD_FILE}"
    VERBATIM
)

message (STATUS "Xsd includes: ${XSD_INCLUDE_DIRS}")
include_directories (SYSTEM
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
)
include_directories (
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
)
add_library (xsdcxx
  inputData.cpp
  ${CMAKE_BINARY_DIR}/xsdcxx/scenario.cxx ${CMAKE_BINARY_DIR}/xsdcxx/scenario.ixx ${CMAKE_BINARY_DIR}/xsdcxx/scenario.hxx
)

if (MSVC)
  set_target_properties (xsdcxx PROPERTIES LINK_FLAGS "${OM_LINK_FLAGS}")
endif (MSVC)
