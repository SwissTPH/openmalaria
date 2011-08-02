# CMake configuration for all code in this dir and below
# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

find_program (JAVAC javac DOC "Path to javac (Java compiler)")
find_program (JAR jar DOC "Path to jar (Java archive tool)")

if (JAVAC AND JAR)
  # inputs: sources and manifest
  set (SOURCE_NAMES
    SchemaTranslator
  )
  set (ST_MANIFEST ${CMAKE_CURRENT_SOURCE_DIR}/manifest.mft)
  
  set (ST_CLASSES "")
  foreach (SOURCE_NAME ${SOURCE_NAMES})
    set (ST_CLASS ${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_NAME}.class)
    list(APPEND ST_CLASSES ${SOURCE_NAME}.class)
    set (JAVA_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_NAME}.java)
    add_custom_command (OUTPUT ${ST_CLASS}
      DEPENDS ${JAVA_FILE}
      COMMAND ${JAVAC} -d ${CMAKE_CURRENT_BINARY_DIR} ${JAVA_FILE}
      COMMENT "Compiling ${JAVA_FILE}"
      VERBATIM
    )
  endforeach (SOURCE_NAME)

  set( ST_JAR ${CMAKE_CURRENT_BINARY_DIR}/SchemaTranslator.jar )
  file(GLOB ST_ALL_CLASSES ${CMAKE_CURRENT_BINARY_DIR}/*.class)
  # above gives absolute paths; need relative because jar tool stores full path
  set(ST_REL_CLASSES "")
  foreach(CLASS ${ST_ALL_CLASSES})
    file(RELATIVE_PATH REL_CLASS ${CMAKE_CURRENT_BINARY_DIR} ${CLASS})
    list(APPEND ST_REL_CLASSES ${REL_CLASS})
  endforeach(CLASS)
  add_custom_command (OUTPUT ${ST_JAR}
    DEPENDS ${ST_CLASSES} # not quite true, but good enough
    COMMAND ${JAR} cmf ${ST_MANIFEST} ${ST_JAR} ${ST_REL_CLASSES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Creating jar ${ST_JAR}"
    VERBATIM
  )
  add_custom_target (st_jar ALL DEPENDS ${ST_JAR})
  file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scenarios ${CMAKE_CURRENT_BINARY_DIR}/translatedScenarios )
else()
  message(WARNING "Couldn't find javac or jar (Java compiler and archive tool). I can still build OpenMalaria, but not the SchemaTranslator tool.")
endif(JAVAC AND JAR)
