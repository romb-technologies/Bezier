cmake_minimum_required(VERSION 2.8)
project(BezierCpp)

add_compile_options(-std=c++14 -march=native -mfpmath=sse -fPIC)

find_package(Eigen3 REQUIRED)

include_directories(${BezierCpp_SOURCE_DIR}/BezierCpp
                    ${EIGEN3_INCLUDE_DIR})

add_library(beziercpp STATIC
            BezierCpp/bezier.cpp
            BezierCpp/polycurve.cpp
            BezierCpp/polyline.cpp)
