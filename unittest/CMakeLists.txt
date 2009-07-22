# CMake configuration for openmalaria's unit tests
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)


include_directories (SYSTEM
  ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx
)

foreach (_file myTest)
  add_custom_command (OUTPUT ${_file}.cpp
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/cxxtestgen.pl --error-printer -o ${_file}.cpp ${CMAKE_CURRENT_SOURCE_DIR}/${_file}.h
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${_file}.h
    COMMENT "Running cxxtestgen on ${_file}.h ..."
    VERBATIM
  )
  add_executable (${_file} ${_file}.cpp)
endforeach (_file)
