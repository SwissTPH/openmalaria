# CMake configuration for XML-reading code
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

set_source_files_properties (scenario.cxx scenario.ixx scenario.hxx PROPERTIES GENERATED TRUE)

set (XSD_FILE ${CMAKE_SOURCE_DIR}/xsdcxx/scenario.xsd)
add_custom_command (OUTPUT scenario.cxx scenario.ixx scenario.hxx
    COMMAND ${XSD_EXECUTABLE} cxx-tree
      --type-naming ucc --function-naming java
      --namespace-map =scnXml
      --generate-inline --generate-doxygen --generate-intellisense
      ${XSD_FILE}
    MAIN_DEPENDENCY ${XSD_FILE}
    COMMENT "Compiling ${XSD_FILE}"
    VERBATIM
)

include_directories (
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
  ${XSD_INCLUDE_DIR}
)
add_library (xsdcxx inputData.cpp scenario.cxx)
