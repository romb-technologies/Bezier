cmake_minimum_required(VERSION 3.15)
project(bezier_python LANGUAGES CXX)

# Pull in the core C++ library (defines target Bezier::Bezier, links Eigen3 PUBLIC).
# EXCLUDE_FROM_ALL keeps its install() rules out of the wheel; we only need the target.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../.. _bezier_core EXCLUDE_FROM_ALL)

find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(bezier NB_STATIC src/bindings.cpp)
target_link_libraries(bezier PRIVATE Bezier::Bezier)  # Eigen comes transitively

install(TARGETS bezier LIBRARY DESTINATION .)
