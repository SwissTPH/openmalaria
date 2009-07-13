add_subdirectory (Clinical)
add_subdirectory (Drug)
add_subdirectory (Pathogenesis)
add_subdirectory (TransmissionModel)
add_subdirectory (WithinHostModel)

#for Visual Studio file groupings:
source_group (Clinical REGULAR_EXPRESSION Clinical/.*)

option (USE_BOINC "Set to enable integration of the BOINC client.")
if (USE_BOINC)
  find_library (BOINC_LIB boinc PATHS ${CMAKE_SOURCE_DIR}/../boinc/lib)
  find_library (BOINC_API_LIB boinc_api PATHS ${CMAKE_SOURCE_DIR}/../boinc/boinc_api)
  include_directories (${CMAKE_SOURCE_DIR}/../boinc ${CMAKE_SOURCE_DIR}/../boinc/api ${CMAKE_SOURCE_DIR}/../boinc/lib)
else (USE_BOINC)
  add_definitions (-DWITHOUT_BOINC)
endif (USE_BOINC)

include_directories (${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx)
link_directories (${CMAKE_BINARY_DIR}/xsdcxx)
aux_source_directory (. SOURCE_CPP)
add_executable (openMalaria ${SOURCE_CPP})
add_dependencies (openMalaria xsdcxx/scenario.hxx)
set_target_properties (openMalaria PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries (openMalaria
  xsdcxx Clinical Drug Pathogenesis TransmissionModel WithinHostModel
  gsl gslcblas xerces-c z pthread
  ${BOINC_API_LIB} ${BOINC_LIB})
