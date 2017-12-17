cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(exact NONE)
include(config/core.cmake)

configure_file(${PROJECT_SOURCE_DIR}/config/exact-config.cmake.in
  ${PROJECT_BINARY_DIR}/config/exact-config.cmake @ONLY)

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(test)

install(FILES ${PROJECT_BINARY_DIR}/config/exact-config.cmake
  ${PROJECT_SOURCE_DIR}/config/use-exact.cmake
  DESTINATION share/exact/)
install(DIRECTORY afh ising DESTINATION include/exact
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(DIRECTORY integral lattice lse DESTINATION include/exact
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
