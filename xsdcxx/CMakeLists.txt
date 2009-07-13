include (FindXsd.cmake)
if (XSD_FOUND)
  add_custom_command (OUTPUT scenario.cxx scenario.ixx scenario.hxx
    COMMAND ${XSD_EXECUTABLE} cxx-tree
      --type-naming ucc --function-naming java
      --namespace-map =scnXml
      --generate-inline --generate-doxygen --generate-intellisense
      scenario.xsd
    MAIN_DEPENDENCY scenario.xsd)
endif (XSD_FOUND)
include_directories (../include .)
add_library (xsdcxx inputData.cpp scenario.cxx)
