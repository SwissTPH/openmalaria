# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

find_package (Python)

# Testing via run.sh or run.py script:
set (OM_BOXTEST_SCHEMA_NAME ${OM_LATEST_SCHEMA})
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/run.py ${CMAKE_CURRENT_BINARY_DIR}/run.py @ONLY)

# working tests (with checkpointing):
set (OM_BOXTEST_NAMES 1 3 4 5 6 11 12 IPT)
# tests with broken checkpointing:
# some of this checkpointing was broken in r569
set (OM_BOXTEST_NC_NAMES 2 9 VecTest VecConstEIR CevCq DummyPKPD)
# Use python script. There is also a shell-script, but it still uses python to validate results.
if (PYTHON_EXECUTABLE)
  foreach (TEST_NAME ${OM_BOXTEST_NAMES})
    add_test (${TEST_NAME} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py -- --checkpoint ${TEST_NAME})
  endforeach (TEST_NAME)
  foreach (TEST_NAME ${OM_BOXTEST_NC_NAMES})
    add_test (${TEST_NAME} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py -- ${TEST_NAME})
  endforeach (TEST_NAME)
endif (PYTHON_EXECUTABLE)
