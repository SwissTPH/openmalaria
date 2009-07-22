# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

foreach (_file densities.csv scenario_6.xsd)
  configure_file(${CMAKE_SOURCE_DIR}/${_file} ${CMAKE_BINARY_DIR}/${_file} COPY_ONLY)
endforeach (_file)
