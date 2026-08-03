cmake_minimum_required(VERSION 3.16)

project(NumericalMethodsLibrary
 VERSION 0.1.0
 DESCRIPTION "A personal C++17 numerical methods library"
 LANGUAGES CXX
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(nml
 src/nml.cc
 include/nml/nml.h
 include/nml/detail/nml_impl.h
)

add_library(NML::nml ALIAS nml)

target_include_directories(nml
 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(nml
 PUBLIC
  Eigen3::Eigen
)

target_compile_features(nml
 PUBLIC
  cxx_std_17
)

set_target_properties(nml PROPERTIES
 CXX_STANDARD 17
 CXX_STANDARD_REQUIRED YES
 CXX_EXTENSIONS OFF
)

if(MSVC)
 target_compile_options(nml PRIVATE /W4)
else()
 target_compile_options(nml PRIVATE -Wall -Wextra -Wpedantic)
endif()
