# CMake configuration for openmalaria's unit tests
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

option (CXXTEST_USE_PYTHON "Use python instead of perl to generate unittest code." OFF)
if (CXXTEST_USE_PYTHON)
  set (CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.py)
else (CXXTEST_USE_PYTHON)
  set (CXXTEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.pl)
endif (CXXTEST_USE_PYTHON)

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

add_custom_command (OUTPUT tests.cpp
  COMMAND ${CXXTEST_SCRIPT} ${CXXTEST_GUI_ARG} --runner=ParenPrinter -o ${CMAKE_CURRENT_BINARY_DIR}/tests.cpp
    PerHostSuite.h
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Generating unittest code with cxxtestgen"
  VERBATIM
)
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

add_test (cxxtest run)

mark_as_advanced (
  CXXTEST_GUI_ARG
  CXXTEST_GUI_LIB
  CXXTEST_SCRIPT
)
