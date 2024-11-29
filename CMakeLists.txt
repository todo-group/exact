cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(exact CXX)
include(cmake/postfix.cmake)

option(EXACT_INSTALL "Install exact." ON)
option(EXACT_BUILD_TESTS "Build all exact's tests." ON)
option(EXACT_BUILD_EXAMPLES "Build all exact's examples." ON)

include_directories(${PROJECT_SOURCE_DIR})

include(boost)
include(eigen3)
include(googletest)
include(standards)
include(lattice)
FetchContent_MakeAvailable(${FetchContents})
# include_directories(${FetchContent_includes})

configure_file(${PROJECT_SOURCE_DIR}/cmake/exact-config.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/exact-config.cmake @ONLY)

add_subdirectory(ising)
add_subdirectory(afh)
add_subdirectory(tfi)
add_subdirectory(test/ising)
add_subdirectory(gallery)
