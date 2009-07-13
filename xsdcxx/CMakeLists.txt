include (FindXsd.cmake)
if (XSD_FOUND)
  set (XSD_FILE ${CMAKE_SOURCE_DIR}/xsdcxx/scenario.xsd)
  add_custom_command (OUTPUT scenario.cxx scenario.ixx scenario.hxx
    COMMAND ${XSD_EXECUTABLE} cxx-tree
      --type-naming ucc --function-naming java
      --namespace-map =scnXml
      --generate-inline --generate-doxygen --generate-intellisense
      ${XSD_FILE}
    MAIN_DEPENDENCY ${XSD_FILE}
    COMMENT "Compiling ${XSD_FILE}"
    VERBATIM)
else (XSD_FOUND)
  # Note: only checking for .hxx file; we also need .ixx and .cxx files (usually all are generated together)
  find_file (SCENARIO_HXX scenario.hxx)
  if (SCENARIO_HXX_NOT_FOUND)
    message (FATAL_ERROR "Unable to find xsd/xsdcxx executable, or pre-build scenario.hxx file!")
  else (SCENARIO_HXX_NOT_FOUND)
    message (STATUS "Unable to find xsd/xsdcxx executable, but found pre-build scenario.hxx file.")
  endif (SCENARIO_HXX_NOT_FOUND)
endif (XSD_FOUND)
include_directories (${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx)
add_library (xsdcxx inputData.cpp scenario.cxx)
