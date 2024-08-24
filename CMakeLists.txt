cmake_minimum_required(VERSION 3.20)

project(XFDTD_MODEL VERSION 0.0.0 LANGUAGES CXX)

set(XFDTD_MODEL_MAIN_PROJECT OFF)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(XFDTD_MODEL_MAIN_PROJECT ON)
endif()

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(xfdtd_core REQUIRED)

set(ASE_READER_PUBLIC_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/external/ase_reader/include)
set(XFDTD_MODEL_PUBLIC_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include ${ASE_READER_PUBLIC_INCLUDE_DIRS})

add_library(xfdtd_model INTERFACE)
target_link_libraries(xfdtd_model INTERFACE xfdtd::xfdtd_core)
target_include_directories(xfdtd_model INTERFACE 
  $<BUILD_INTERFACE:${XFDTD_MODEL_PUBLIC_INCLUDE_DIRS}>
  $<INSTALL_INTERFACE:include>)

add_subdirectory(test)

set(XFDTD_MODEL_CMAKE_DIR ${CMAKE_CURRENT_BINARY_DIR}/cmake)
set(XFDTD_MODEL_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

install(
  EXPORT xfdtd_model-targets
  FILE xfdtd_model-targets.cmake
  NAMESPACE xfdtd::
  DESTINATION ${XFDTD_MODEL_INSTALL_PREFIX}/lib/cmake/xfdtd_model
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/xfdtd_model-config.cmake.in
  ${XFDTD_MODEL_CMAKE_DIR}/xfdtd_model-config.cmake
  INSTALL_DESTINATION ${XFDTD_MODEL_CMAKE_DIR}
)

export(EXPORT xfdtd_model-targets
  FILE ${XFDTD_MODEL_CMAKE_DIR}/xfdtd_model-targets.cmake
  NAMESPACE xfdtd::
)

install(
  FILES ${XFDTD_MODEL_CMAKE_DIR}/xfdtd_model-config.cmake
  DESTINATION ${XFDTD_MODEL_INSTALL_PREFIX}/lib/cmake/xfdtd_model
)

install(TARGETS xfdtd_model
  EXPORT xfdtd_model-targets
  INCLUDES DESTINATION include
)

install(DIRECTORY ${XFDTD_MODEL_PUBLIC_INCLUDE_DIRS}
  DESTINATION ${XFDTD_MODEL_INSTALL_PREFIX}
)
