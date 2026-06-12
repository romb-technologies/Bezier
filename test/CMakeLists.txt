find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(unit_tests
  curve_test.cpp
  polycurve_test.cpp
  utils_test.cpp
  coefficients_test.cpp
  approx_test.cpp
  known_bugs_test.cpp
)
target_link_libraries(unit_tests PRIVATE bezier GTest::gtest GTest::gtest_main)

# Per-test ctest registration; TIMEOUT turns infinite-loop bugs into failures
gtest_discover_tests(unit_tests PROPERTIES TIMEOUT 60)
