#  Find the SUITESPARSE includes and library
#
# The following variables are set:
#
#
#   SUITESPARSS_INCLUDE_DIR
#   SUITESPARSE_LIBRARIES
#   SUITESPARSE_FOUND
#
#
#=============================================================================
# Copyright 2010, Martin Koehler
# http://www-user.tu-chemnitz.de/~komart/
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# Changelog:
#    - Sep. 11, 2012   Martin Koehler
#    - June 25, 2013   Martin Koehler, add _libdir and _incdir
#    - Aug. 20, 2017    Maximilian Behr

include(FindPackageHandleStandardArgs)

# GET_INC_LIB_DIR MACRO
include(CMakeHelpers)
get_inc_lib_dir(_incdir _libdir)


find_package(AMD      QUIET)
find_package(COLAMD   QUIET)
find_package(UMFPACK  QUIET)
find_package(CHOLMOD  QUIET)
find_package(CXSPARSE  QUIET)

show_variable(AMD_FOUND)
show_variable(COLAMD_FOUND)
show_variable(UMFPACK_FOUND)
show_variable(CHOLMOD_FOUND)
show_variable(CXSPARSE_FOUND)

if(AMD_FOUND AND COLAMD_FOUND AND UMFPACK_FOUND AND CHOLMOD_FOUND AND CXSPARSE_FOUND)

    # search in user given directories
    find_path(SUITESPARSE_INCLUDE_DIR   NAMES SuiteSparse_config.h                  PATHS ${SUITESPARSE} PATH_SUFFIXES include      NO_DEFAULT_PATH)
    find_library(SUITESPARSE_LIBRARIES  NAMES suitesparseconfig suitesparse_config  PATHS ${SUITESPARSE} PATH_SUFFIXES lib lib64    NO_DEFAULT_PATH)

    # search in other directories
    find_path(SUITESPARSE_INCLUDE_DIR   NAMES SuiteSparse_config.h                  PATHS ${_incdir} /usr /opt PATH_SUFFIXES include local/include include/suitesparse local/include/suitesparse)
    find_library(SUITESPARSE_LIBRARIES  NAMES suitesparseconfig suitesparse_config  PATHS ${_libdir} /usr /opt PATH_SUFFIXES lib lib64 local/lib local/lib64)

    if(SUITESPARSE_LIBRARIES AND SUITESPARSE_INCLUDE_DIR)
        set(SUITESPARSE_INCLUDE_DIR
            ${AMD_INCLUDE_DIR}
            ${COLAMD_INCLUDE_DIR}
            ${UMFPACK_INCLUDE_DIR}
            ${CHOLMOD_INCLUDE_DIR}
            ${CXSPARSE_INCLUDE_DIR})
        set(SUITESPARSE_LIBRARIES
            ${AMD_LIBRARIES}
            ${COLAMD_LIBRARIES}
            ${UMFPACK_LIBRARIES}
            ${CHOLMOD_LIBRARIES}
            ${CXSPARSE_LIBRARIES})
    endif()

endif()



find_package_handle_standard_args(SUITESPARSE  DEFAULT_MSG  SUITESPARSE_LIBRARIES SUITESPARSE_INCLUDE_DIR)
mark_as_advanced(SUITESPARSE_INCLUDE_DIR  SUITESPARSE_LIBRARIES)


