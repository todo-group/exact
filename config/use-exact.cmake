#
# This module is provided as ROKKO_USE_FILE by exact-config.cmake.  It can
# be included in a project to load the needed compiler and linker
# settings to use Exact
#

if(NOT EXACT_USE_FILE_INCLUDED)
  set(EXACT_USE_FILE_INCLUDED 1)
  include_directories(${EXACT_ROOT_DIR}/include)
endif(NOT EXACT_USE_FILE_INCLUDED)
