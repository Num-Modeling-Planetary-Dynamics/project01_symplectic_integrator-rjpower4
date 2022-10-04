# ========================================================================================
# File   : src/eaps/CMakeLists.txt
# Brief  : Coordinating CMake file for the `eaps` library
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

add_subdirectory(include)
add_subdirectory(src)

if (EAPS_BUILD_TESTS)
    add_subdirectory(test)
endif ()