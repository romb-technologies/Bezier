cmake_minimum_required(VERSION 3.12.4)
project(bezier
  LANGUAGES CXX
  VERSION 0.3.2)

find_package(Eigen3 REQUIRED)

set(Bezier_SRC
  ${PROJECT_SOURCE_DIR}/src/bezier.cpp
  ${PROJECT_SOURCE_DIR}/src/polycurve.cpp
  ${PROJECT_SOURCE_DIR}/src/utils.cpp
  )

set(Bezier_INC
  ${PROJECT_SOURCE_DIR}/include/Bezier/coefficients.h
  ${PROJECT_SOURCE_DIR}/include/Bezier/declarations.h
  ${PROJECT_SOURCE_DIR}/include/Bezier/bezier.h
  ${PROJECT_SOURCE_DIR}/include/Bezier/polycurve.h
  ${PROJECT_SOURCE_DIR}/include/Bezier/utils.h
  )

# Options
option(BUILD_SHARED_LIBS "Build shared library (.so) instead of static one (.a)" OFF)

# build rules

if(BUILD_SHARED_LIBS)
  add_library(bezier SHARED ${Bezier_SRC})
else()
  add_library(bezier STATIC ${Bezier_SRC})
endif()

target_include_directories(bezier PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(bezier PUBLIC Eigen3::Eigen)

set_target_properties(bezier PROPERTIES
  VERSION ${PROJECT_VERSION}
  CXX_EXTENSIONS OFF
  POSITION_INDEPENDENT_CODE ON
  PUBLIC_HEADER "${Bezier_INC}"
)
target_compile_features(bezier PUBLIC cxx_std_17)
target_compile_options(bezier PRIVATE -O3 -DNDEBUG)

# install rules
install(TARGETS bezier
  EXPORT bezier-export DESTINATION "lib"
  PUBLIC_HEADER DESTINATION "include/Bezier")
install(EXPORT bezier-export DESTINATION "lib/cmake/Bezier" FILE BezierConfig.cmake)

add_library(bezier::bezier ALIAS bezier)
