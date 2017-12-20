cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(exact NONE)
include(config/core.cmake)

configure_file(${PROJECT_SOURCE_DIR}/config/exact-config.cmake.in
  ${PROJECT_BINARY_DIR}/config/exact-config.cmake @ONLY)

include_directories(${PROJECT_SOURCE_DIR}/integral ${PROJECT_SOURCE_DIR}/lattice ${PROJECT_SOURCE_DIR}/lse)
add_subdirectory(test)
