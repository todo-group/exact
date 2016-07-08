cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(exact NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(C CXX)

enable_testing()
include(CTest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(ising)
add_subdirectory(afh)
