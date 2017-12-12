#  Copyright (C) 2012-2017 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# Disable in-source builds
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    message(FATAL_ERROR "In source builds are disabled. Please use a separate build directory")
endif()

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# RPATH fix
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

set(_ALPS_ROOT_ENV $ENV{ALPS_ROOT})
if(ALPS_ROOT_DIR OR _ALPS_ROOT_ENV)
  find_package(ALPS HINTS ${ALPS_ROOT_DIR} ${_ALPS_ROOT_ENV})
  message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
  include(${ALPS_USE_FILE})
else(ALPS_ROOT_DIR OR _ALPS_ROOT_ENV)
  find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  add_definitions(-DALPS_INDEP_SOURCE)
endif(ALPS_ROOT_DIR OR _ALPS_ROOT_ENV)
unset(_ALPS_ROOT_ENV)

# Build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

enable_language(C CXX)

include(add_iotest)
