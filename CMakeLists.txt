cmake_minimum_required (VERSION 2.6)
project (OpenMalaria CXX)

# Statically link libgcc; isn't going to work when other C++ libraries are dynamically linked
IF( ${CMAKE_COMPILER_IS_GNUCXX} )
  option (STATICALLY_LINK "Attempt to link libraries, including stdc++, statically.")
  if (STATICALLY_LINK)
    link_directories (${CMAKE_SOURCE_DIR}/lib)
    ADD_DEFINITIONS(-static-libgcc)
    SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -static-libgcc")
  endif (STATICALLY_LINK)
ENDIF( ${CMAKE_COMPILER_IS_GNUCXX} )

option (BUILD_SHARED_LIBS "Link xsdcxx, WithinHostModel, etc. libraries dynamically instead of statically (almost certainly not wanted)." OFF)
MARK_AS_ADVANCED (BUILD_SHARED_LIBS)

add_subdirectory (model)
add_subdirectory (xsdcxx)
