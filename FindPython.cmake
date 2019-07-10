# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines PYTHON_EXECUTABLE. Doesn't abort if not found.

# We need Python 3.x which is now often the default version.
# On Windows, install paths may help
find_program (PYTHON_EXECUTABLE
  NAMES python3 python3.exe python python.exe
  PATHS
    "C:\\Python37-x64" "C:\\Python37"
    "C:\\Python36-x64" "C:\\Python36"
    "C:\\Python35-x64" "C:\\Python35"
    "C:\\Python34-x64" "C:\\Python34"
    "C:\\Python33-x64" "C:\\Python33"
  DOC "Path to python (required).")
message (STATUS "Found Python interpreter: ${PYTHON_EXECUTABLE}")

mark_as_advanced (PYTHON_EXECUTABLE)
