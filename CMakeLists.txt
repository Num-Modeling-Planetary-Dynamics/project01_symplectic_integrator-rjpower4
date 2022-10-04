# ========================================================================================
# File   : CMakeLists.txt
# Brief  : Top-level CMake file for Project 1
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================
cmake_minimum_required(VERSION 3.22.1)

# ----------------------------------------------------------------------------------------
# Project Instantiation
# ----------------------------------------------------------------------------------------
project(
    EapsProject1
    DESCRIPTION "Symplectic integration of solar system bodies"
    HOMEPAGE_URL "https://github.com/Num-Modeling-Planetary-Dynamics/project01_symplectic_integrator-rjpower4"
    LANGUAGES CXX
)

# ----------------------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------------------
option(EAPS_BUILD_TESTS "Build the tests for the project")

# ----------------------------------------------------------------------------------------
# Set up testing
# ----------------------------------------------------------------------------------------
if (EAPS_BUILD_TESTS)
    include(FetchContent)
        FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
    )

    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest) 

    enable_testing()
endif ()

# ----------------------------------------------------------------------------------------
# Subdirectory Descent
# ----------------------------------------------------------------------------------------
add_subdirectory(src)