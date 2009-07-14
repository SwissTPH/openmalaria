# list common paths here:
FIND_PATH(Z_INCLUDE_DIR zconf.h
  PATHS ${CMAKE_SOURCE_DIR}/../zlib/include
)

if (NOT Z_INCLUDE_DIR)
  message (SEND_ERROR "Unable to find zlib includes (zconf.h)")
endif (NOT Z_INCLUDE_DIR)
message (STATUS "Z_INCLUDE_DIR is ${Z_INCLUDE_DIR}")

find_library (Z_LIB NAMES z zlib zdll
  PATHS ${CMAKE_SOURCE_DIR}/lib ${CMAKE_SOURCE_DIR}/../zlib/lib
)
if (NOT Z_LIB)
  message (SEND_ERROR "Unable to find zlib library")
endif (NOT Z_LIB)

MARK_AS_ADVANCED(
  Z_INCLUDE_DIR
)
