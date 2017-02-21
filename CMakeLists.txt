cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(exact NONE)
include(config/core.cmake)

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(test)
