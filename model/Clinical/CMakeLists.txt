include_directories (${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/xsdcxx)
aux_source_directory (. SOURCE_CPP)
get_filename_component(DIR_NAME ${CMAKE_CURRENT_BINARY_DIR} NAME_WE)
add_library (${DIR_NAME} ${SOURCE_CPP})
