cmake_minimum_required(VERSION 3.15)
project(bezier_c LANGUAGES CXX)

# Pull in the core C++ library (defines target Bezier::Bezier, links Eigen3 PUBLIC).
# EXCLUDE_FROM_ALL keeps its install rules out of this build; we only need the target.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../.. _bezier_core EXCLUDE_FROM_ALL)

add_library(bezier_c SHARED src/bezier_c.cpp)
target_include_directories(bezier_c PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(bezier_c PRIVATE Bezier::Bezier)  # Eigen comes transitively
target_compile_features(bezier_c PRIVATE cxx_std_17)
