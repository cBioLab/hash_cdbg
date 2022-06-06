# - Try to find SDSL
# Once done this will define
#  LIBSDSL_FOUND - System has SDSL
#  LIBSDSL_INCLUDE_DIRS - The SDSL include directories
#  LIBSDSL_LIBRARIES - The libraries needed to use SDSL
#  LIBSDSL_DEFINITIONS - Compiler switched required for using SDSL

find_package(PkgConfig)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{HOME})
pkg_check_modules(PC_LIBSDSL QUIET sdsl-lite)
set(LIBSDSL_DEFINITIONS ${LIBSDSL_CFLAGS_OTHER})

# find sdsl include directory
find_path(LIBSDSL_INCLUDE_DIR
  NAMES
    sdsl/int_vector.hpp
  HINTS
    ${PC_LIBSDSL_INCLUDEDIR}
    ${PC_LIBSDSL_INCLUDE_DIRS}
  PATHS
    $ENV{LIBSDSL_ROOT}
    $ENV{LIBSDSL_INCLUDE_DIR}
    ${LIBSDSL_ROOT}
    /usr
    /usr/local
    $ENV{HOME}
    $ENV{HOME}/usr
  PATH_SUFFIXES
    include
)

#find sdsl library
find_library(LIBSDSL_LIBRARY
  NAMES
    sdsl
  HINTS
    ${PC_LIBSDSL_LIBDIR}
    ${PC_LIBSDSL_LIBRARY_DIRS}
  PATHS
    $ENV{LIBSDSL_ROOT}
    $ENV{LIBSDSL_INCLUDE_DIR}
    ${LIBSDSL_ROOT}
    /usr
    /usr/local
    $ENV{HOME}
    $ENV{HOME}/usr
  PATH_SUFFIXES
    lib
)

mark_as_advanced(LIBSDSL_INCLUDE_DIR LIBSDSL_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibSDSL
  DEFAULT_MSG
    LIBSDSL_LIBRARY
    LIBSDSL_INCLUDE_DIR
)

# if (PC_LIBSDSL_FOUND)
#   set(LIBSDSL_LIBRARIES ${LIBSDSL_LIBRARY} ${PC_LIBSDSL_LIBRARIES})
#   set(LIBSDSL_INCLUDE_DIRS ${LIBSDSL_INCLUDE_DIR} ${PC_LIBSDSL_INCLUDE_DIRS})
# else()
  set(LIBSDSL_LIBRARIES ${LIBSDSL_LIBRARY})
  set(LIBSDSL_INCLUDE_DIRS ${LIBSDSL_INCLUDE_DIR})
# endif()
