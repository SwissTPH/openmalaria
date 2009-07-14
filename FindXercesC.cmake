find_library (XERCES_C_LIB NAMES xerces-c xerces-c_static xerces-c_static_2 xerces-c_static_3
  PATHS ${CMAKE_SOURCE_DIR}/lib
)
if (NOT XERCES_C_LIB)
  message (SEND_ERROR "Unable to find xerces-c library")
endif (NOT XERCES_C_LIB)

MARK_AS_ADVANCED(
  Z_INCLUDE_DIR
)
