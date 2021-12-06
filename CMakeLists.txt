cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(exact CXX)
include(cmake/postfix.cmake)

option(EXACT_INSTALL "Install exact." ON)
option(EXACT_BUILD_TESTS "Build all exact's tests." ON)
option(EXACT_BUILD_EXAMPLES "Build all exact's examples." ON)

include_directories(${PROJECT_SOURCE_DIR})

include(googletest)
include(standards)
include(lattice)
FetchContent_MakeAvailable(googletest standards lattice)
include_directories(${standards_SOURCE_DIR} ${lattice_SOURCE_DIR})

configure_file(${PROJECT_SOURCE_DIR}/cmake/exact-config.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/exact-config.cmake @ONLY)

# Eigen3
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

# Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(ising)
add_subdirectory(afh)
add_subdirectory(gallery)
