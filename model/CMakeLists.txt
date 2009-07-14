# CMake configuration for all code in this dir and below
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

find_package (GSL)
find_package (Z)


#for Visual Studio file groupings:
source_group (Clinical REGULAR_EXPRESSION Clinical/.*)
source_group (Drug REGULAR_EXPRESSION Drug/.*)
source_group (Pathogenesis REGULAR_EXPRESSION Pathogenesis/.*)
source_group (TransmissionModel REGULAR_EXPRESSION TransmissionModel/.*)
source_group (WithinHostModel REGULAR_EXPRESSION WithinHostModel/.*)

option (BOINC_INTEGRATION "Set to enable integration of the BOINC client (otherwise it's not used).")
if (BOINC_INTEGRATION)
  find_library (BOINC_LIB boinc PATHS ${CMAKE_SOURCE_DIR}/../boinc/lib)
  find_library (BOINC_API_LIB boinc_api PATHS ${CMAKE_SOURCE_DIR}/../boinc/boinc_api)
  include_directories (${CMAKE_SOURCE_DIR}/../boinc ${CMAKE_SOURCE_DIR}/../boinc/api ${CMAKE_SOURCE_DIR}/../boinc/lib)
  if (NOT BOINC_LIB OR NOT BOINC_API_LIB)
    message (FATAL_ERROR "Unable to find a boinc library (boinc or boinc_api)")
  endif (NOT BOINC_LIB OR NOT BOINC_API_LIB)
  set (BOINC_LIBRARIES ${BOINC_API_LIB} ${BOINC_LIB})
  mark_as_advanced (BOINC_API_LIB BOINC_LIB)
else (BOINC_INTEGRATION)
  add_definitions (-DWITHOUT_BOINC)
endif (BOINC_INTEGRATION)

if (UNIX)
  find_library (PTHREAD_LIBRARIES pthread PATHS ${CMAKE_SOURCE_DIR}/lib /usr/lib /usr/local/lib)
  if (NOT PTHREAD_LIBRARIES)
    message (SEND_ERROR "PThread library not found")
  endif (NOT PTHREAD_LIBRARIES)
endif (UNIX)
mark_as_advanced (BOINC_LIBRARIES PTHREAD_LIBRARIES)

include_directories (
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
  ${XSD_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
  ${Z_INCLUDE_DIR}
)

link_directories (${CMAKE_BINARY_DIR}/xsdcxx)

aux_source_directory (. OpenMalaria_CPP)
aux_source_directory (Clinical Clinical_CPP)
aux_source_directory (Drug Drug_CPP)
aux_source_directory (Pathogenesis Pathogenesis_CPP)
aux_source_directory (TransmissionModel TransmissionModel_CPP)
aux_source_directory (WithinHostModel WithinHostModel_CPP)

add_executable (openMalaria
  ${OpenMalaria_CPP}
  ${Clinical_CPP}
  ${Drug_CPP}
  ${Pathogenesis_CPP}
  ${TransmissionModel_CPP}
  ${WithinHostModel_CPP}
)
#Really I want to make openMalaria depend on XSD having been run, but this doesn't work:
#add_dependencies (openMalaria ${CMAKE_BINARY_DIR}/xsdcxx/scenario.hxx)
set_target_properties (openMalaria PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries (openMalaria
  xsdcxx
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${BOINC_LIBRARIES}
)
