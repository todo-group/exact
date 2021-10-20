cmake_minimum_required(VERSION 3.4)
project(exact NONE)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/standards/cmake)
include(core)

include_directories(${PROJECT_SOURCE_DIR})

configure_file(${PROJECT_SOURCE_DIR}/config/exact-config.cmake.in
  ${PROJECT_BINARY_DIR}/config/exact-config.cmake @ONLY)

# Eigen3
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

# Boost
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)

# Google Test
set(GOOGLETEST_VERSION 1.11.0)
add_subdirectory(googletest/googletest)
include_directories(googletest/googletest/include)
enable_testing()

include_directories(standards)
include_directories(lattice)

add_subdirectory(ising)
add_subdirectory(afh)
add_subdirectory(gallery)
