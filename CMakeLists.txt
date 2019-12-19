cmake_minimum_required(VERSION 3.1)
project(BezierCpp
  LANGUAGES CXX
  VERSION 0.2.0)

add_compile_options(-std=c++11 -mfpmath=sse -fPIC)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

set(BezierCpp_SRC
  src/bezier.cpp
  src/polycurve.cpp
  src/offsetted_point.cpp
  )

set(BezierCpp_INC
  include/BezierCpp/declarations.h
  include/BezierCpp/legendre_gauss.h
  include/BezierCpp/bezier.h
  include/BezierCpp/polycurve.h
  include/BezierCpp/offsetted_point.h
  )

# Option ot set Legendre_Gauss precision (default is 10)
if(NOT LEGENDRE_GAUSS_PRECISION)
  message(STATUS "Setting Legendre Gauss precision to default value (10).")
  set(LEGENDRE_GAUSS_PRECISION  10)
endif()


# build rules

add_library(BezierCpp SHARED ${BezierCpp_SRC})
target_include_directories(BezierCpp PUBLIC
  $<BUILD_INTERFACE:${BezierCpp_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:lib/BezierCpp-${PROJECT_VERSION}>)

target_compile_definitions(BezierCpp PRIVATE LEGENDRE_GAUSS_N=${LEGENDRE_GAUSS_PRECISION})

set_target_properties(BezierCpp PROPERTIES
  VERSION ${PROJECT_VERSION})


# install rules

install(FILES ${BezierCpp_INC} DESTINATION "include/BezierCpp")
install(TARGETS BezierCpp EXPORT BezierCpp-export DESTINATION "lib")
install(EXPORT BezierCpp-export DESTINATION "lib/cmake/BezierCpp" FILE BezierCpp-config.cmake)

