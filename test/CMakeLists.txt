project(bezier_unit_tests)

find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIR})

add_executable(unit_tests unit_tests.cpp unit_tests.hpp)
target_link_libraries(unit_tests PUBLIC bezier ${GTEST_MAIN_LIBRARY} ${GTEST_LIBRARY} pthread)
target_include_directories(unit_tests PUBLIC ${PROJECT_SOURCE_DIR}/include/Bezier/)

