# ========================================================================================
# File   : 
# Brief  : 
# Author : Rolfe Power <rpower@purdue.edu>
# ========================================================================================

include(GoogleTest)

add_executable(keplerian_elements_test keplerian_elements_test.cpp)
target_link_libraries(keplerian_elements_test PRIVATE GTest::gtest_main EAPS::eaps)
gtest_discover_tests(keplerian_elements_test)

add_executable(conic_test conic_test.cpp)
target_link_libraries(conic_test PRIVATE GTest::gtest_main EAPS::eaps)
gtest_discover_tests(conic_test)