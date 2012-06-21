# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

find_package (Python)

# Testing via run.sh or run.py script:
configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/run.py
  ${CMAKE_CURRENT_BINARY_DIR}/run.py
  @ONLY
)

# working tests (with checkpointing):
set (OM_BOXTEST_NAMES
  1 2 3 4 5 6 9 10 11 12
  1DayIntervs
  Cohort
  EffectiveDrug
  Empirical
  ESCMTest
  ESTS
  IPT
  IRS30
  Molineaux
  MSAT
  NamawalaArabiensis
  NoInterv
  Penny
  Rach5IC
  R_0
  VecFullTest
  VecTest
  VecMonthly
  NoMPDLarviciding
  SimpleMPDLarviciding
)
# tests with broken checkpointing:
set (OM_BOXTEST_NC_NAMES)
# Disabled due to "in-progress" work: (none)
# Use python script.
if (PYTHON_EXECUTABLE)
  foreach (TEST_NAME ${OM_BOXTEST_NAMES})
    add_test (${TEST_NAME} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py ${TEST_NAME} -- --checkpoint)
  endforeach (TEST_NAME)
  foreach (TEST_NAME ${OM_BOXTEST_NC_NAMES})
    add_test (${TEST_NAME} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py -- ${TEST_NAME})
  endforeach (TEST_NAME)
else (PYTHON_EXECUTABLE)
  message(WARNING "Tests are disabled (Python is needed to run them)")
endif (PYTHON_EXECUTABLE)
