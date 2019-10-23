cmake_minimum_required(VERSION 2.8)
project(BezierCpp)

add_compile_options(-std=c++11 -mfpmath=sse -fPIC)

find_package(Eigen3 REQUIRED)

include_directories(${BezierCpp_SOURCE_DIR}/include
                    ${EIGEN3_INCLUDE_DIR})

message(STATUS "${Eigen3_DIR}")

add_library(beziercpp STATIC
            src/bezier.cpp
            src/polycurve.cpp
            src/offsetted_point.cpp)

target_compile_definitions(beziercpp PRIVATE LEGENDRE_GAUSS_N=10)
