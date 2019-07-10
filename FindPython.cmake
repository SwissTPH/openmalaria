# CMake tests configuration for OpenMalaria
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)
#
# Defines PYTHON_EXECUTABLE. Doesn't abort if not found.

# We need python3 which is now often the default version
find_program (PYTHON_EXECUTABLE python3 python
  DOC "Path to python (required).")

mark_as_advanced (PYTHON_EXECUTABLE)
