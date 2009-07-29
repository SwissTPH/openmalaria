# CMake configuration for openmalaria's unit tests
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

find_package (Python)
find_program (PERL_EXECUTABLE perl
  DOC "Path to perl. Only one of perl or python is needed.")

if (PERL_EXECUTABLE AND PYTHON_EXECUTABLE)
  option (CXXTEST_USE_PYTHON "Use python instead of perl to generate unittest code." OFF)
  if (CXXTEST_USE_PYTHON)
    set (CXXTEST_INTERPRETER ${PYTHON_EXECUTABLE})
    set (CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.py)
  else (CXXTEST_USE_PYTHON)
    set (CXXTEST_INTERPRETER ${PERL_EXECUTABLE})
    set (CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.pl)
  endif (CXXTEST_USE_PYTHON)
elseif (PERL_EXECUTABLE)
  set (CXXTEST_INTERPRETER ${PERL_EXECUTABLE})
  set (CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.pl)
elseif (PYTHON_EXECUTABLE)
  set (CXXTEST_INTERPRETER ${PYTHON_EXECUTABLE})
  set (CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.py)
else (PERL_EXECUTABLE AND PYTHON_EXECUTABLE)
  set (CXXTEST_INTERPRETER "NOTFOUND")
endif (PERL_EXECUTABLE AND PYTHON_EXECUTABLE)

option (CXXTEST_GUI "Draw a graphical progress bar" OFF)
set (CXXTEST_GUI_ARG "")
if (CXXTEST_GUI)
  # also supports QT, but that has more dependencies.
  if (WIN32)
    set (CXXTEST_GUI_ARG "--gui=Win32Gui")
  elseif (UNIX)
    set (CXXTEST_GUI_ARG "--gui=X11Gui")
    find_library (CXXTEST_GUI_LIB X11)
    if (NOT CXXTEST_GUI_LIB)
      message (SEND_ERROR "Requested use of a gui, but can't find X11 library${CXXTEST_GUI_LIB}!")
    endif (NOT CXXTEST_GUI_LIB)
  endif (WIN32)
endif (CXXTEST_GUI)


include_directories (SYSTEM
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/unittest
  ${CMAKE_BINARY_DIR}/xsdcxx
)

set (CXXTEST_HEADERS PerHostSuite.h VectorEmergenceSuite.h)

if (CXXTEST_INTERPRETER)
  # The second COMMAND, cp, is to update the source copy. We can't use configure_file
  # because that runs at configure time, not build time.
  if (UNIX)
	set (CP cp)
  else (UNIX)
	set (CP copy)
  endif (UNIX)
  add_custom_command (OUTPUT tests.cpp
    COMMAND ${CXXTEST_INTERPRETER} ${CXXTEST_SCRIPT} ${CXXTEST_GUI_ARG} --runner=ParenPrinter -o ${CMAKE_CURRENT_BINARY_DIR}/tests.cpp
      ${CXXTEST_HEADERS}
    COMMAND ${CP} ${CMAKE_CURRENT_BINARY_DIR}/tests.cpp ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp
    DEPENDS ${CXXTEST_HEADERS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating unittest code with cxxtestgen"
    VERBATIM
  )
else (CXXTEST_INTERPRETER)
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp ${CMAKE_CURRENT_BINARY_DIR}/tests.cpp COPYONLY)
    message (STATUS "Unable to find perl or python executable. Thus I won't be able to rebuild tests.cpp if required.")
  else (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp)
    message (SEND_ERROR "Unable to find perl or python executable. Thus I can't build tests.cpp, and I can't find an existing version to use.")
  endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests.cpp)
endif (CXXTEST_INTERPRETER)

# This file is needed. "run" won't search for it, so it must be in the current directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/VectorEmergenceSuite.txt ${CMAKE_CURRENT_BINARY_DIR}/VectorEmergenceSuite.txt COPYONLY)

add_executable (run tests.cpp)
target_link_libraries (run
  model
  xsdcxx
  ${CXXTEST_GUI_LIB}
  ${GSL_LIBRARIES}
  ${XERCESC_LIBRARIES}
  ${Z_LIBRARIES}
  ${PTHREAD_LIBRARIES}
  ${BOINC_LIBRARIES}
)

if (MSVC)
  set_target_properties (model PROPERTIES LINK_FLAGS "${OM_LINK_FLAGS}")
endif (MSVC)

add_test (unittest run -keep)

mark_as_advanced (
  CXXTEST_GUI_ARG
  CXXTEST_GUI_LIB
)
