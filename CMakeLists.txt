cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(exact NONE)
include(config/core.cmake)

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(integral)
add_subdirectory(log_sum_exp)
add_subdirectory(ising)
add_subdirectory(afh)
