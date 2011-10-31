# - Try to find PaStiX
# Once done this will define
#
#  PASTIX_FOUND        - system has PaStiX
#  PASTIX_INCLUDE_DIRS - include directories for PaStiX
#  PASTIX_LIBRARIES    - libraries for PaStiX

# Check for header file
find_path(PASTIX_INCLUDE_DIRS pastix.h
  HINTS ${PASTIX_DIR} $ENV{PASTIX_DIR} ${PASTIX_DIR}/include $ENV{PASTIX_DIR}/include
  PATH_SUFFIXES install
  DOC "Directory where the PaStiX header is located"
 )

# Check for PaStiX library
find_library(PASTIX_LIBRARY pastix
  HINTS ${PASTIX_DIR} $ENV{PASTIX_DIR} ${PASTIX_DIR}/lib $ENV{PASTIX_DIR}/lib
  PATH_SUFFIXES install
  DOC "The PaStiX library"
  )

# Check for rt library
find_library(RT_LIBRARY rt
  DOC "The RT library"
  )

# Collect libraries
set(PASTIX_LIBRARIES ${PASTIX_LIBRARY} ${RT_LIBRARY})

mark_as_advanced(
  PASTIX_INCLUDE_DIRS
  PASTIX_LIBRARY
  PASTIX_LIBRARIES
  )

# Try to run a test program that uses PaStiX
if (PASTIX_INCLUDE_DIRS AND PASTIX_LIBRARIES)

  set(CMAKE_REQUIRED_INCLUDES  ${PASTIX_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${PASTIX_LIBRARIES})

  # Build and run test program
  include(CheckCXXSourceRuns)
  check_cxx_source_runs("
/* Test program pastix */

//#include <pastix.h>

int main()
{

  return 0;
}
" PASTIX_TEST_RUNS)

endif()
# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PASTIX
  "PaStiX could not be found. Be sure to set PASTIX_DIR."
  PASTIX_LIBRARIES PASTIX_INCLUDE_DIRS PASTIX_TEST_RUNS)
