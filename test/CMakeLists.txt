# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

# Testing via runTests.sh script:
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/run.sh ${CMAKE_CURRENT_BINARY_DIR}/run.sh @ONLY)

option (OMTEST_USE_CMAKE "Use CMake scripts instead of runTests.sh to run tests (not implemented)." OFF)
foreach (_file densities.csv Nv0scenarioVector.txt scenario_6.xsd)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPY_ONLY)
endforeach (_file)
set (OMTEST_NAMES 1 2 3 4 5 6 9 11 12 IPT DummyPKPD CevCq)
if (OMTEST_USE_CMAKE)
  message (SEND_ERROR "Don't use OMTEST_USE_CMAKE. CMake doesn't provide a scripting language capable of running our tests it seems.")
  foreach (TEST_NAME ${OMTEST_NAMES})
    # This doesn't work because the run${TEST_NAME} target is not run before the test
    # (unless ALL is specified, which runs it after every build)!
    # In other words a SCRIPTING language is needed, not a build config, to do this.
    file (REMOVE checkpoint output${TEST_NAME}.txt)
    add_custom_target (run${TEST_NAME}
      ${CMAKE_BINARY_DIR}/openMalaria --checkpoint --scenario ${CMAKE_CURRENT_SOURCE_DIR}/scenario${TEST_NAME}.xml
      COMMAND ${CMAKE_BINARY_DIR}/openMalaria --checkpoint --scenario ${CMAKE_CURRENT_SOURCE_DIR}/scenario${TEST_NAME}.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      VERBATIM
    )
    add_test (${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/compareOutputsFloat.py ${CMAKE_CURRENT_SOURCE_DIR}/original${TEST_NAME}.txt ${CMAKE_CURRENT_BINARY_DIR}/output.txt 1)
  endforeach (TEST_NAME)
else (OMTEST_USE_CMAKE)
  if (UNIX) # uses a bash script & a python script
    # working tests:
    foreach (TEST_NAME ${OMTEST_NAMES})
      add_test (${TEST_NAME} ${CMAKE_CURRENT_BINARY_DIR}/run.sh --checkpoint ${TEST_NAME})
    endforeach (TEST_NAME)
  endif (UNIX)
endif (OMTEST_USE_CMAKE)
set_tests_properties (${OMTEST_NAMES}
  PROPERTIES PASS_REGULAR_EXPRESSION "No differences, ok\\.\\.\\.")

mark_as_advanced (OMTEST_NAMES)
