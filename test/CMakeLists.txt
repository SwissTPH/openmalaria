# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

# Testing via runTests.sh script:
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/run.sh ${CMAKE_CURRENT_BINARY_DIR}/run.sh @ONLY)

option (OMTEST_USE_CMAKE "Use CMake scripts instead of runTests.sh to run tests (not implemented)." OFF)
foreach (_file densities.csv Nv0scenarioVector.txt scenario_6.xsd)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPY_ONLY)
endforeach (_file)
if (OMTEST_USE_CMAKE)
else (OMTEST_USE_CMAKE)
  if (UNIX) # uses a bash script & a python script
    # working tests:
    set (TEST_NAMES 1 2 3 4 5 6 9 11 12 IPT DummyPKPD CevCq)
    foreach (TEST_NAME ${TEST_NAMES})
      add_test (${TEST_NAME} ${CMAKE_CURRENT_BINARY_DIR}/run.sh --checkpoint ${TEST_NAME})
    endforeach (TEST_NAME)
    set_tests_properties (${TEST_NAMES}
      PROPERTIES PASS_REGULAR_EXPRESSION "No differences, ok\\.\\.\\.")
  endif (UNIX)
endif (OMTEST_USE_CMAKE)
