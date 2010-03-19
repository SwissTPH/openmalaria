# CMake configuration for openmalaria's unit tests
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

find_package (Python)
find_program (PERL_EXECUTABLE perl
  DOC "Path to perl. Only one of perl or python is needed.")

if (PERL_EXECUTABLE AND PYTHON_EXECUTABLE)
  option (OM_CXXTEST_USE_PYTHON "Use python instead of perl to generate unittest code." OFF)
  if (OM_CXXTEST_USE_PYTHON)
    set (OM_CXXTEST_INTERPRETER ${PYTHON_EXECUTABLE})
    set (OM_CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.py)
  else (OM_CXXTEST_USE_PYTHON)
    set (OM_CXXTEST_INTERPRETER ${PERL_EXECUTABLE})
    set (OM_CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.pl)
  endif (OM_CXXTEST_USE_PYTHON)
else (PERL_EXECUTABLE AND PYTHON_EXECUTABLE)
  if (PERL_EXECUTABLE)
   set (OM_CXXTEST_INTERPRETER ${PERL_EXECUTABLE})
   set (OM_CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.pl)
  else (PERL_EXECUTABLE)
    if (PYTHON_EXECUTABLE)
      set (OM_CXXTEST_INTERPRETER ${PYTHON_EXECUTABLE})
      set (OM_CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.py)
    else (PYTHON_EXECUTABLE)
      set (OM_CXXTEST_INTERPRETER "NOTFOUND")
    endif (PYTHON_EXECUTABLE)
  endif (PERL_EXECUTABLE)
endif (PERL_EXECUTABLE AND PYTHON_EXECUTABLE)

option (OM_CXXTEST_GUI "Draw a graphical progress bar" OFF)
set (OM_CXXTEST_GUI_ARG "")
if (OM_CXXTEST_GUI)
  # also supports QT, but that has more dependencies.
  if (WIN32)
    set (OM_CXXTEST_GUI_ARG "--gui=Win32Gui")
  elseif (UNIX)
    set (OM_CXXTEST_GUI_ARG "--gui=X11Gui")
    find_library (OM_CXXTEST_GUI_LIB X11)
    if (NOT OM_CXXTEST_GUI_LIB)
      message (SEND_ERROR "Requested use of a gui, but can't find X11 library${OM_CXXTEST_GUI_LIB}!")
    endif (NOT OM_CXXTEST_GUI_LIB)
  endif (WIN32)
endif (OM_CXXTEST_GUI)


# For yaml-cpp: I did use this to read in data for the VectorAnopheles unittest.
# It's another not-the-easiest dependency to install, so maybe don't use it?
#find_library (YAML_CPP_LIB NAMES yaml-cpp)
#if (NOT YAML_CPP_LIB)
#  message (FATAL_ERROR "Unable to find yaml-cpp library. (You can avoid this dependency by disabling OM_CXXTEST_ENABLE)")
#endif (NOT YAML_CPP_LIB)

set (CMAKE_CONFIGURED 1)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/TestPaths.h ${CMAKE_CURRENT_BINARY_DIR}/configured/TestPaths.h @ONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/scenario.xml ${CMAKE_CURRENT_BINARY_DIR}/configured/scenario.xml @ONLY)
include_directories (
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/unittest
  ${CMAKE_BINARY_DIR}/unittest
  ${CMAKE_BINARY_DIR}/xsdcxx
)

set (OM_CXXTEST_HEADERS
  PerHostTransmissionSuite.h
  #VectorAnophelesSuite.h needs updating for vector model changes
  ExtraAsserts.h	# must appear after at least some of the above
  HoshenPkPdSuite.h
  LSTMPkPdSuite.h
  CheckpointSuite.h
  DummyInfectionSuite.h
  EmpiricalInfectionSuite.h
  InfectionImmunitySuite.h
  AgeGroupDataSuite.h
  #DescriptiveInfectionSuite.h
)

if (OM_CXXTEST_INTERPRETER)
  add_custom_command (OUTPUT tests.cpp
    COMMAND ${OM_CXXTEST_INTERPRETER} ${OM_CXXTEST_SCRIPT} ${OM_CXXTEST_GUI_ARG} --runner=ParenPrinter -o ${CMAKE_CURRENT_BINARY_DIR}/tests.cpp
      ${OM_CXXTEST_HEADERS}
    DEPENDS ${OM_CXXTEST_HEADERS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating unittest code with cxxtestgen"
    VERBATIM
  )
else (OM_CXXTEST_INTERPRETER)
  message (SEND_ERROR "Unable to find perl or python executable. Thus I can't build tests.cpp which is required to run the unittests.")
endif (OM_CXXTEST_INTERPRETER)

# This file is needed to run the EmpiricalInfection unittest:
configure_file (${CMAKE_SOURCE_DIR}/test/autoRegressionParameters.csv ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_executable (run tests.cpp)
target_link_libraries (run
  model
  ${OM_CXXTEST_GUI_LIB}
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${BOINC_LIBRARIES}
  #${YAML_CPP_LIB}
)

if (MSVC)
  set_target_properties (run PROPERTIES
    LINK_FLAGS "${OM_LINK_FLAGS}"
    COMPILE_FLAGS "${OM_COMPILE_FLAGS}"
  )
endif (MSVC)

add_test (unittest run -keep)
set_tests_properties (unittest
  PROPERTIES PASS_REGULAR_EXPRESSION "^Running [0-9]+ tests.*OK![ \r\n]*$")

mark_as_advanced (
  OM_CXXTEST_GUI_ARG
  OM_CXXTEST_GUI_LIB
)
