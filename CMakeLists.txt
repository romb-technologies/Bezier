cmake_minimum_required(VERSION 3.1)
project(BezierCpp
  LANGUAGES CXX
  VERSION 0.2.0)

add_compile_options(-std=c++11 -march=native -fPIC)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

set(BezierCpp_SRC
  src/bezier.cpp
  src/polycurve.cpp
  )

set(BezierCpp_INC
  include/BezierCpp/declarations.h
  include/BezierCpp/legendre_gauss.h
  include/BezierCpp/bezier.h
  include/BezierCpp/polycurve.h
  )

# Options

if(NOT LEGENDRE_GAUSS_PRECISION)
  message(STATUS "Setting Legendre Gauss precision to default value (10).")
  set(LEGENDRE_GAUSS_PRECISION  10)
endif()

option(BUILD_SHARED_LIBS "Build shared library (.so) instead of static one (/.a)" OFF)

# build rules

if(BUILD_SHARED_LIBS)
  add_library(BezierCpp SHARED ${BezierCpp_SRC})
else()
  add_library(BezierCpp STATIC ${BezierCpp_SRC})
endif()

target_include_directories(BezierCpp PUBLIC
  $<BUILD_INTERFACE:${BezierCpp_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_compile_definitions(BezierCpp PRIVATE LEGENDRE_GAUSS_N=${LEGENDRE_GAUSS_PRECISION})

set_target_properties(BezierCpp PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(BezierCpp PROPERTIES PUBLIC_HEADER "${BezierCpp_INC}")

# install rules
install(TARGETS BezierCpp
  EXPORT BezierCpp-export DESTINATION "lib"
  PUBLIC_HEADER DESTINATION "include/BezierCpp")
install(EXPORT BezierCpp-export DESTINATION "lib/cmake/BezierCpp" FILE BezierCppConfig.cmake)
