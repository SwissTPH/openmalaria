# CMake configuration for all code in this dir and below
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)


#for Visual Studio file groupings:
source_group (Clinical REGULAR_EXPRESSION Clinical/.*)
source_group (Drug REGULAR_EXPRESSION Drug/.*)
source_group (Pathogenesis REGULAR_EXPRESSION Pathogenesis/.*)
source_group (Transmission REGULAR_EXPRESSION Transmission/.*)
source_group (WithinHost REGULAR_EXPRESSION WithinHost/.*)

# Don't use aux_source_directory on . because we don't want to compile openMalaria.cpp in to the lib.
set (Model_CPP
  EntoIntervention.cpp population.cpp timer.cpp InfectionIncidenceModel.cpp
  gzstream.C simulation.cpp summary.cpp BoincWrapper.cpp global.cpp
  intervention.cpp proteome.cpp human.cpp GSLWrapper.cpp
)
aux_source_directory (Clinical Clinical_CPP)
aux_source_directory (Drug Drug_CPP)
aux_source_directory (Pathogenesis Pathogenesis_CPP)
aux_source_directory (Transmission Transmission_CPP)
aux_source_directory (WithinHost WithinHost_CPP)

include_directories (SYSTEM
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
  ${XSD_INCLUDE_DIRS}
  ${XERCESC_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${Z_INCLUDE_DIRS}
)

add_library (model
  ${Model_CPP}
  ${Clinical_CPP}
  ${Drug_CPP}
  ${Pathogenesis_CPP}
  ${Transmission_CPP}
  ${WithinHost_CPP}
)

# May want to remove LIBCMT or MSVCRTD like this, BUT if libraries are compiled with one there'll still be link errors
if (MSVC)
  set_target_properties (model
    PROPERTIES LINK_FLAGS /NODEFAULTLIB:LIBCMT)
endif (MSVC)
