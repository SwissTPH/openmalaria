# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines PYTHON_EXECUTABLE. Doesn't abort if not found.

find_program (PYTHON_EXECUTABLE python
  DOC "Path to python. Needed for testing.")

mark_as_advanced (PYTHON_EXECUTABLE)
