cmake_minimum_required(VERSION 3.15)
project(bezier_wasm LANGUAGES CXX)

# Configure with: emcmake cmake -S bindings/wasm -B bindings/wasm/build

set(CORE ${CMAKE_CURRENT_SOURCE_DIR}/../..)
set(CABI ${CMAKE_CURRENT_SOURCE_DIR}/../c-abi)

# Compile core + shared C-ABI shim directly into one module (no find_package(Eigen3) reliance
# under the Emscripten sysroot).
add_executable(bezier_wasm
  ${CORE}/src/bezier.cpp
  ${CORE}/src/polycurve.cpp
  ${CORE}/src/utils.cpp
  ${CABI}/src/bezier_c.cpp)

target_include_directories(bezier_wasm PRIVATE ${CORE}/include ${CABI}/include)
target_compile_features(bezier_wasm PRIVATE cxx_std_17)

# Eigen is header-only; prefer find_package, else fall back to a host include path.
find_package(Eigen3 QUIET)
if(TARGET Eigen3::Eigen)
  target_link_libraries(bezier_wasm PRIVATE Eigen3::Eigen)
else()
  set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3" CACHE PATH "Path to Eigen headers")
  target_include_directories(bezier_wasm PRIVATE ${EIGEN3_INCLUDE_DIR})
endif()

# Exception support is OFF by default in Emscripten; the shim's try/catch and core throws
# (e.g. length(t) out of range) need it, otherwise the error path aborts the module.
# -faligned-allocation: Emscripten libc++ otherwise rejects std::vector<Eigen::Vector2d>
#   (over-aligned) as non-destructible.
# -include vector: core headers use std::vector via transitive includes that Emscripten's
#   libc++ doesn't provide; force it without editing the core headers.
target_compile_options(bezier_wasm PRIVATE -fexceptions -faligned-allocation -include vector)

# Build the exported-functions list straight from the shim header (every bz_* symbol),
# plus the allocator — so it stays in sync without a hand-maintained ~60-name list.
file(READ ${CABI}/include/bezier_c.h _bezier_c_header)
string(REGEX MATCHALL "bz_[a-z_]+" _bz_names "${_bezier_c_header}")
list(REMOVE_DUPLICATES _bz_names)
set(_exported "'_malloc','_free'")
foreach(_n ${_bz_names})
  string(APPEND _exported ",'_${_n}'")
endforeach()

target_link_options(bezier_wasm PRIVATE
  -fexceptions
  -sMODULARIZE=1
  -sEXPORT_ES6=1
  -sENVIRONMENT=web,node
  -sALLOW_MEMORY_GROWTH=1
  "-sEXPORTED_RUNTIME_METHODS=['UTF8ToString','HEAPF64','HEAP32']"
  "-sEXPORTED_FUNCTIONS=[${_exported}]")

set_target_properties(bezier_wasm PROPERTIES
  SUFFIX ".mjs"
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/dist)
