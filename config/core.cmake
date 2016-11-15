#  Copyright (C) 2012-2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
#
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

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

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

enable_language(C CXX)

include(add_iotest)
