# a script taken from http://www.codesynthesis.com/pipermail/xsd-users/2006-March/000269.html and heavily modified:
# Locate Xsd from code synthesis include paths and binary
# Xsd can be found at http://codesynthesis.com/products/xsd/
# Written by Frederic Heem, frederic.heem _at_ telsey.it
#
# This module defines
# XSD_INCLUDE_DIRS, where to find elements.hxx, etc.
# XSD_EXECUTABLE, where is the xsd compiler

# list common paths here:
FIND_PATH(XSD_INCLUDE_DIRS xsd/cxx/parser/elements.hxx
  PATHS "[HKEY_CURRENT_USER\\software\\xsd\\include]"
  "[HKEY_CURRENT_USER]\\xsd\\include]"
  $ENV{XSDDIR}/include
  /usr/local/include
  /usr/include
  "C:/Program Files/CodeSynthesis XSD 3.2/include"
  "D:/Program Files/CodeSynthesis XSD 3.2/include"
  "C:/Program Files/CodeSynthesis XSD 3.3/include"
  "C:/Program Files/CodeSynthesis XSD 4.0/include"
  "C:/Program Files (x86)/CodeSynthesis XSD 3.3/include"
  "C:/Program Files (x86)/CodeSynthesis XSD 4.0/include"
  "/usr/local/Cellar/xsd/4.0.0/include"
  ${CMAKE_SOURCE_DIR}/../xsd/include
  ${CMAKE_SOURCE_DIR}/../xsd/libxsd
)

FIND_PROGRAM(XSD_EXECUTABLE 
  NAMES xsdcxx xsd
  PATHS "[HKEY_CURRENT_USER\\xsd\\bin]" $ENV{XSDDIR}/bin
  "C:/Program Files/CodeSynthesis XSD 3.2/bin"
  "D:/Program Files/CodeSynthesis XSD 3.2/bin"
  "C:/Program Files/CodeSynthesis XSD 3.3/bin"
  "C:/Program Files/CodeSynthesis XSD 4.0/bin"
  "C:/Program Files (x86)/CodeSynthesis XSD 3.3/bin"
  "C:/Program Files (x86)/CodeSynthesis XSD 4.0/bin"
  "/usr/local/Cellar/xsd/4.0.0/bin"
  ${CMAKE_SOURCE_DIR}/../xsd/bin
)

if (NOT XSD_INCLUDE_DIRS)
  message (FATAL_ERROR "Unable to find xsd include files (xsd/cxx/parser/elements.hxx)")
endif (NOT XSD_INCLUDE_DIRS)
if (NOT XSD_EXECUTABLE)
  message (FATAL_ERROR "Unable to find xsd or xsdcxx executable")
endif (NOT XSD_EXECUTABLE)

MARK_AS_ADVANCED(
  XSD_INCLUDE_DIRS
  XSD_EXECUTABLE
)
