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
set (CXXTEST_GUI_LIB "")
if (CXXTEST_GUI)
  # also QT, but that has more dependencies.
  if (WIN32)
    set (CXXTEST_GUI_ARG "--gui=Win32Gui")
  elseif (UNIX)
    set (CXXTEST_GUI_ARG "--gui=X11Gui")
    set (CXXTEST_GUI_LIB "X11")
  endif (WIN32)
endif (CXXTEST_GUI)


include_directories (SYSTEM
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
)

add_custom_command (OUTPUT tests.cpp
  COMMAND ${CXXTEST_SCRIPT} ${CXXTEST_GUI_ARG} --error-printer -o tests.cpp ${CMAKE_CURRENT_SOURCE_DIR}/myTest.h
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/myTest.h
  COMMENT "Generating unittest code with cxxtestgen"
  VERBATIM
)
add_executable (run tests.cpp)
target_link_libraries (run ${CXXTEST_GUI_LIB})

mark_as_advanced (
  CXXTEST_GUI_ARG
  CXXTEST_GUI_LIB
  CXXTEST_SCRIPT
)
