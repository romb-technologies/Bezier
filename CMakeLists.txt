cmake_minimum_required(VERSION 3.12.4)
project(Bezier
  LANGUAGES CXX
  VERSION 0.3.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-fPIC -O2)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR} include)

set(Bezier_SRC
  ${PROJECT_SOURCE_DIR}/src/bezier.cpp
  ${PROJECT_SOURCE_DIR}/src/polycurve.cpp
)

set(Bezier_INC
  ${PROJECT_SOURCE_DIR}/include/Bezier/declarations.h
  ${PROJECT_SOURCE_DIR}/include/Bezier/bezier.h
  ${PROJECT_SOURCE_DIR}/include/Bezier/polycurve.h
)

option(BUILD_SHARED_LIBS "Build shared library (.so) instead of static one (.a)" OFF)

if(BUILD_SHARED_LIBS)
  add_library(Bezier SHARED ${Bezier_SRC})
else()
  add_library(Bezier STATIC ${Bezier_SRC})
endif()

target_include_directories(Bezier PUBLIC
  $<BUILD_INTERFACE:${Bezier_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(Bezier PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(Bezier PROPERTIES PUBLIC_HEADER "${Bezier_INC}")

install(TARGETS Bezier
  EXPORT Bezier-export DESTINATION "lib"
  PUBLIC_HEADER DESTINATION "include/Bezier")

install(EXPORT Bezier-export NAMESPACE Bezier:: DESTINATION "lib/cmake/Bezier" FILE BezierConfig.cmake)
