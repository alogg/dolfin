# - Try to find ParMETIS
# Once done this wil define
#
#  PARMETIS_FOUND       - system has ParMETIS
#  PARMETIS_INCLUDE_DIR - include directories for ParMETIS
#  PARMETIS_LIBRARIES   - libraries for ParMETIS

find_path(PARMETIS_INCLUDE_DIR parmetis.h
  PATHS  $ENV{PARMETIS_DIR}/include
  DOC "Directory where the ParMETIS header files are located"
)

find_library(PARMETIS_LIBRARY parmetis
  PATHS $ENV{PARMETIS_DIR}/lib
)

find_library(METIS_LIBRARY metis
  PATHS $ENV{PARMETIS_DIR}/lib
)

set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

# Standard package handling
find_package_handle_standard_args(ParMETIS
                                  "ParMETIS could not be found."
                                  PARMETIS_INCLUDE_DIR PARMETIS_LIBRARIES)
