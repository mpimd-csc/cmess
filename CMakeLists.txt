#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
#               2009-2018
#



CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(MESS)
SET(MESS_MAJOR_VERSION 1)
SET(MESS_MINOR_VERSION 0)
SET(MESS_PATCH_VERSION 0)
SET(MESS_VERSION  ${MESS_MAJOR_VERSION}.${MESS_MINOR_VERSION}.${MESS_PATCH_VERSION})
ENABLE_LANGUAGE(Fortran C CXX)
STRING(TIMESTAMP MESS_CONFIG_TIME "%Y-%m-%d %H:%M:%S" UTC)


# try to get informations from git
FIND_PROGRAM(GIT_EXECUTABLE git)
IF(GIT_EXECUTABLE)
# Get the current working branch
EXECUTE_PROCESS(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE MESS_GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
EXECUTE_PROCESS(
    COMMAND git log -1 --format=%H
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE MESS_GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
ELSE()
    SET(MESS_GIT_BRANCH "")
    SET(MESS_GIT_COMMIT_HASH "")
ENDIF()


# Set PATH for CmakeModules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/CMakeModules)

SET(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# deactivate mac osx rpath
SET(CMAKE_MACOSX_RPATH 0)

# CmakeHelpers and System Tests
INCLUDE(CMakeHelpers)
INCLUDE(CheckTypeSize)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckAttributeExists)
INCLUDE(ExternalProject)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckBlasZdotcMKL)
INCLUDE(FortranCInterface)


# -----Options for build-----
OPTION(BUILD_SHARED_LIBS    "Build Shared Libraries"            ON)
OPTION(DEBUG                "Debug Mode"                        OFF)
OPTION(OPENMP               "Build with OpenMP support"         OFF)
OPTION(SSE4                 "Build with SSE4 optimization"      OFF)
OPTION(SSE2                 "Build with SSE2 optimization"      OFF)
OPTION(MESS64               "Build with 64-bit integers"        OFF)
OPTION(TESTS                "Build Software Tests"              ON)
OPTION(LARGETESTS           "Build Large Software Tests"        OFF)
OPTION(HUGETESTS            "Build Huge Software Tests"         OFF)
OPTION(COVERAGE             "Enable Code Coverage"              OFF)
OPTION(TUTORIALS            "Build Tutorials"                   ON)
OPTION(DOC                  "Build API documentation"           OFF)
OPTION(DOCPDF               "Build API documentation as PDF"    OFF)
OPTION(MATHJAX              "Enable MATHJAX in Documentation"   OFF)
OPTION(INSTMODULES          "Install Additional CMake Modules"  OFF)

# External packages
OPTION(SUITESPARSE          "Build with SuiteSparse"            ON)
OPTION(ARPACK               "Build with ARPACK"                 OFF)
OPTION(SUPERLU              "Build with SUPERLU"                OFF)
OPTION(X11                  "Build with X11 Plot Support"       OFF)
OPTION(MATIO                "Build Matlab File I/O support"     OFF)

# Interfaces
OPTION(PYTHON               "Build Py-M.E.S.S."                 OFF)
OPTION(MATLAB               "Build MEX-M.E.S.S."                OFF)
OPTION(OCTAVE               "Build OCT-M.E.S.S."                OFF)

# GNU Compatible Install Directories (https://www.gnu.org/prep/standards/html_node/Directory-Variables.html)
INCLUDE(GNUInstallDirs)


# Multplicity of the Leading Dimension
SET_IFNDEF(LD_MULTIPLICITY 1)


# Check for BLAS and LAPACK
IF(((DEFINED BLAS AND NOT DEFINED LAPACK) OR (DEFINED LAPACK AND NOT DEFINED BLAS)))
    MESSAGE(FATAL_ERROR "BLAS and LAPACK have to be defined.")
ENDIF()




# find math library
FOREACH(LIBNAME m)
    FIND_LIBRARY(LIB_${LIBNAME} NAMES ${LIBNAME})
    IF(LIB_${LIBNAME})
        ADD(MESS_LIBS ${LIB_${LIBNAME}})
        UNSET(LIB_${LIBNAME} CACHE)
    ELSE()
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Cannot find Library ${LIBNAME} on your system.")
    ENDIF()
ENDFOREACH()


# add include dirs in src and build
ADD(MESS_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/include)

# Fortran Name Mangling
FORTRANCINTERFACE_HEADER(${CMAKE_BINARY_DIR}/fortran_translate.h MACRO_NAMESPACE F77_)
FORTRANCINTERFACE_VERIFY()

# ------------- System Checks -----------------------------------
ADD(CMAKE_REQUIRED_LIBRARIES ${MESS_LIBS})
ADD(CMAKE_EXTRA_INCLUDE_FILES string.h time.h)

CHECK_TYPE_SIZE("size_t"            MESS_HAVE_SIZE_T        LANGUAGE C)
CHECK_TYPE_SIZE("ptrdiff_t"         MESS_HAVE_PTRDIFF_T     LANGUAGE C)
CHECK_INCLUDE_FILES(sys/sysinfo.h   HAVE_SYSINFO_H)
CHECK_INCLUDE_FILES(stdint.h        HAVE_STDINT_H)

CHECK_FUNCTION_EXISTS(clock_gettime MESS_HAVE_CLOCK_GETTIME)
CHECK_FUNCTION_EXISTS(strndup       MESS_HAVE_STRNDUP)
CHECK_FUNCTION_EXISTS(nearbyint     MESS_HAVE_NEARBYINT)
CHECK_SYMBOL_EXISTS(isfinite        "math.h" MESS_HAVE_ISFINITE)
CHECK_FUNCTION_EXISTS(strncmp       MESS_HAVE_STRNCMP)

CHECK_FUNCTION_EXISTS(creal         HAVE_CREAL)
CHECK_FUNCTION_EXISTS(cimag         HAVE_CIMAG)
CHECK_INCLUDE_FILES(complex.h       HAVE_COMPLEX_H)
CHECK_TYPE_SIZE("double _Complex"   DCPX_T                  LANGUAGE C)


IF(NOT(HAVE_COMPLEX_H OR HAVE_CREAL OR HAVE_CIMAG OR DCPX_T))
    MESSAGE(FATAL_ERROR "Your C compiler does not support complex numbers. C-M.E.S.S. cannot be built without it. Please update your C compiler. ")
ENDIF()



#
# OPTION MESS64, experimental
#
IF(MESS64)
    # add fortran flags
    CHECK_ADD_FORTRAN_FLAG("-fdefault-integer-8" SUCCESS1)
    CHECK_ADD_FORTRAN_FLAG("-i8" SUCCESS2)

    # check success
    IF(NOT(SUCCESS1 OR SUCCESS2))
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Cannot add -i8 or -fdefault-integer-8 flag to the Fortran compiler.")
    ENDIF()

    SET(MESS_USE_INTEGER8 TRUE)
    MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "64 Bit integer support is activated, but experimental.")
ENDIF()



#
# C and Fortran Flags
#
IF(DEBUG)
    SET(FLAGS  "-g3 -ggdb -fPIC" "-Wall" "-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE")
ELSE()
    SET(FLAGS "-g3" "-fPIC" "-Wall"  "-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE")
ENDIF()
FOREACH(FLAG ${FLAGS})
    CHECK_ADD_C_FLAG(${FLAG} SUCCESS)
    IF(NOT SUCCESS)
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Cannot add Flag ${FLAG} to your C Compiler.")
    ENDIF()

    CHECK_ADD_CXX_FLAG(${FLAG} SUCCESS)
    IF(NOT SUCCESS)
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Cannot add Flag ${FLAG} to your CXX Compiler.")
    ENDIF()

    CHECK_ADD_FORTRAN_FLAG(${FLAG} SUCCESS)
    IF(NOT SUCCESS)
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Cannot add Flag ${FLAG} to your Fortran Compiler.")
    ENDIF()
ENDFOREACH()


#
# OPTION SSE2
#
IF(SSE2)
    CHECK_ADD_C_FLAG("-msse2" SUCCESS1)
    CHECK_ADD_FORTRAN_FLAG("-msse2" SUCCESS2)

    IF(NOT(SUCCESS1 AND SUCCESS2))
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Compiler Flags -msse2 are not supported by your compiler")
    ENDIF()

ENDIF()


#
# OPTION SSE4
#
IF(SSE4)
    CHECK_ADD_C_FLAG("-msse4" SUCCESS1)
    CHECK_ADD_FORTRAN_FLAG("-msse4" SUCCESS2)

    IF(NOT(SUCCESS1 AND SUCCESS2))
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Compiler Flags -msse4 are not supported by your compiler")
    ENDIF()

ENDIF()


#
# OPTION DEBUG, remove O3 flags
#
SET(CMAKE_BUILD_TYPE "Release")

SET(MESS_DEBUG FALSE)
IF(DEBUG)
    SET(CMAKE_BUILD_TYPE "Debug")
    SET(MESS_DEBUG TRUE)
ENDIF()
STRING(REGEX REPLACE "-O[0-9]?" "-O2" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
STRING(REGEX REPLACE "-O[0-9]?" "-O2" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
STRING(REGEX REPLACE "-O[0-9]?" "-O2" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})

#
# OPTION COVERAGE
#
IF(COVERAGE)
    # check if we use a gnu compiler
    IF(NOT CMAKE_C_COMPILER_ID  STREQUAL "GNU")
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Code coverage is only possible with GNU Compilers")
    ENDIF()

    # add flags
    CHECK_ADD_C_FLAG("-g -fprofile-arcs -ftest-coverage"            SUCCESS1)
    CHECK_ADD_FORTRAN_FLAG("-g -fprofile-arcs -ftest-coverage"      SUCCESS2)
    SET(CMAKE_EXE_LINKER_FLAGS      "${CMAKE_EXE_LINKER_FLAGS} -g -fprofile-arcs -ftest-coverage")
    SET(CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} -g -fprofile-arcs -ftest-coverage")

    IF(NOT(SUCCESS1 AND SUCCESS2))
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Compiler Flags for coverage report are not supported by your compiler")
    ENDIF()

    # find gcov library, gcov lcov and genhtml program and add new target make coverage
    FIND_PROGRAM(LCOV       lcov)
    FIND_PROGRAM(GENHTML    genhtml)
    FIND_PROGRAM(GCOV       gcov)
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found lcov: ${LCOV}")
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found gcov: ${GCOV}")
    MESSAGE_COLOR(STATUS COLOR_CYAN "Found genhtml: ${GENHTML}")

    IF(LCOV AND GENHTML AND GCOV)
        # add target make coverage
        ADD_CUSTOM_TARGET(
            coverage
            COMMAND ${LCOV} --zerocounters --directory "${CMAKE_BINARY_DIR}/"
            COMMAND make -C "${CMAKE_BINARY_DIR}" test
            COMMAND ${LCOV} -c --gcov-tool ${GCOV} --output-file "${CMAKE_BINARY_DIR}/lcov.info" --directory "${CMAKE_BINARY_DIR}/"
            COMMAND ${GENHTML} --legend --output "${CMAKE_BINARY_DIR}/coverage" "${CMAKE_BINARY_DIR}/lcov.info"
            COMMAND  ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Coverage Report available: ${CMAKE_BINARY_DIR}/coverage"
            )
        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make coverage")
    ELSE()
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Coverage will not work due to missing tools. Need lcov, gcov and genhtml.")
    ENDIF()

ENDIF()


#
# Search for External Packages
#
MESSAGE("")
MESSAGE_COLOR(COLOR_BOLD_GREEN "Search for external libraries.")
MESSAGE("")


## handle MATLAB OPTION ON
IF(MATLAB)

    IF(OPENMP)
        SET(OPENMP OFF)
        MESSAGE_COLOR(STATUS COLOR_MAGENTA "OPTION MATLAB IS ON -> OPTION OPENMP IS TURNED OFF")
    ENDIF()

    IF(SUPERLU)
        SET(SUPERLU OFF)
        MESSAGE_COLOR(STATUS COLOR_MAGENTA "OPTION MATLAB IS ON -> OPTION SUPERLU IS TURNED OFF")
    ENDIF()

    IF(ARPACK)
        SET(ARPACK OFF)
        SET(MESS_HAVE_ARPACK TRUE)
        MESSAGE_COLOR(STATUS COLOR_MAGENTA "OPTION MATLAB IS ON AND OPTION ARPACK IS ON -> TAKE ARPACK FROM MATLAB AND TURN OPTION ARPACK OFF")
    ENDIF()

ENDIF()


# Threading
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for threading library.")
FIND_PACKAGE(Threads REQUIRED)
ADD(MESS_LIBS ${CMAKE_THREAD_LIBS_INIT})

# gfortran
MESSAGE("")
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for gfortran library.")
FIND_PACKAGE(Gfortran REQUIRED)
ADD(MESS_LIBS ${GFORTRAN_LIBRARIES})


# ----------------------- BLAS and LAPACK --------------------------
IF(NOT BLAS)
    # NO BLAS/LAPACK WAS GIVEN VIA CMAKE COMMAND ARGUMENTS

    # BLAS
    IF(MESS64)
        MESSAGE("")
        MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for 64-bit BLAS library.")
        FIND_PACKAGE(BLAS64)
    ELSE()
        MESSAGE("")
        MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for 32-bit BLAS Library.")
        FIND_PACKAGE(BLAS)
    ENDIF()

    # LAPACK
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for LAPACK library.")
    FIND_PACKAGE(LAPACK)

ELSE()
    # BLAS/LAPACK WAS GIVEN VIA CMAKE COMMAND ARGUMENTS
    MESSAGE_COLOR(STATUS COLOR_CYAN "Use given BLAS/LAPACK")
    SET(BLAS_FOUND TRUE)
    SET(BLAS_LIBRARIES ${BLAS})
    SET(LAPACK_FOUND TRUE)
    SET(LAPACK_LIBRARIES ${LAPACK})
ENDIF()

IF(NOT (BLAS_FOUND AND LAPACK_FOUND))
    MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "BLAS and/or LAPACK not found. Please specify both via -DBLAS=<blaslibrary> -DLAPACK=<lapacklibrary>")
ELSE()
    MESSAGE(STATUS "BLAS: ${BLA_VENDOR} - ${BLAS_LIBRARIES}")
    MESSAGE(STATUS "LAPACK: ${LAPACK_LIBRARIES}")
    # store libs
    ADD(MESS_LIBS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

    # Adjust MKL Interfce
    CHECKBLASZDOTCMKL("${BLAS_LIBRARIES}" "${BLAS_LINKER_FLAGS}" ${MESS64})
    IF(ZDOTC_MKL)
        MESSAGE(STATUS "Use BLAS ZDOTC MKL fix.")
        ADD_DEFINITIONS(-DZDOTC_MKL)
    ENDIF()


    # CHECK FOR APPLE BLAS
    STRING(REGEX MATCH "Accelerate" ACC_MATCH ${BLAS_LIBRARIES})
    IF(ACC_MATCH)
        MESSAGE(STATUS "Use Apple BLAS.")
        SET(MESS_USE_APPLE_BLAS TRUE)
    ENDIF()

    #Check If BLAS is OpenBLAS
    ADD(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
    CHECK_FUNCTION_EXISTS(openblas_set_num_threads MESS_USE_OPENBLAS)

ENDIF()
# -------------------------------------------------------------------

# OpenMP support
SET(MESS_HAVE_OPENMP FALSE)
IF(OPENMP)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for OpenMP Library.")
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
        SET(CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        SET(CMAKE_C_FLAGS_DEBUG         "${CMAKE_C_FLAGS_DEBUG} ${OpenMP_C_FLAGS}")
        SET(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        SET(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
        SET(CMAKE_SHARED_LINKER_FLAGS   "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_LD_FLAGS}")
        SET(MESS_HAVE_OPENMP TRUE)
    ENDIF()
ENDIF()


# Check for Suitesparse
SET(MESS_HAVE_AMD           FALSE)
SET(MESS_HAVE_COLAMD        FALSE)
SET(MESS_HAVE_UMFPACK       FALSE)
SET(MESS_HAVE_CHOLMOD       FALSE)
SET(MESS_HAVE_CSPARSE       FALSE)
SET(MESS_USE_SUITESPARSE3   FALSE)

IF(SUITESPARSE)
    #search for SuiteSparse
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for SuiteSparse Library.")
    FIND_PACKAGE(SUITESPARSE)

    IF(SUITESPARSE_FOUND)
        MESSAGE(STATUS "Use new SuiteSparse > 4.0.")
        SET(MESS_HAVE_AMD      TRUE)
        SET(MESS_HAVE_COLAMD   TRUE)
        SET(MESS_HAVE_UMFPACK  TRUE)
        SET(MESS_HAVE_CHOLMOD  TRUE)
        SET(MESS_HAVE_CSPARSE  TRUE)
        ADD(MESS_INCLUDE_DIR    ${SUITESPARSE_INCLUDE_DIR})
        ADD(MESS_LIBS           ${SUITESPARSE_LIBRARIES})
    ELSE()
        MESSAGE(STATUS "Search for old SuiteSparse < 4.0.")

        FIND_PACKAGE(UFCONFIG REQUIRED)
        ADD(MESS_INCLUDE_DIR            ${UFCONFIG_INCLUDE_DIR})
        ADD(SUITESPARSE_INCLUDE_DIR     ${UFCONFIG_INCLUDE_DIR})

        FIND_PACKAGE(AMD REQUIRED)
        ADD(MESS_INCLUDE_DIR            ${AMD_INCLUDE_DIR})
        ADD(MESS_LIBS                   ${AMD_LIBRARIES})
        ADD(SUITESPARSE_INCLUDE_DIR     ${AMD_INCLUDE_DIR})
        ADD(SUITESPARSE_LIBRARIES       ${AMD_LIBRARIES})

        FIND_PACKAGE(COLAMD REQUIRED)
        ADD(MESS_INCLUDE_DIR            ${COLAMD_INCLUDE_DIR})
        ADD(MESS_LIBS                   ${COLAMD_LIBRARIES})
        ADD(SUITESPARSE_INCLUDE_DIR     ${COLAMD_INCLUDE_DIR})
        ADD(SUITESPARSE_LIBRARIES       ${COLAMD_LIBRARIES})

        FIND_PACKAGE(CHOLMOD REQUIRED)
        ADD(MESS_INCLUDE_DIR            ${CHOLMOD_INCLUDE_DIR})
        ADD(MESS_LIBS                   ${CHOLMOD_LIBRARIES})
        ADD(SUITESPARSE_INCLUDE_DIR     ${CHOLAMD_INCLUDE_DIR})
        ADD(SUITESPARSE_LIBRARIES       ${CHOLAMD_LIBRARIES})

        FIND_PACKAGE(UMFPACK REQUIRED)
        ADD(MESS_INCLUDE_DIR            ${UMFPACK_INCLUDE_DIR})
        ADD(MESS_LIBS                   ${UMFPACK_LIBRARIES})
        ADD(SUITESPARSE_INCLUDE_DIR     ${UMFPACK_INCLUDE_DIR})
        ADD(SUITESPARSE_LIBRARIES       ${UMFPACK_LIBRARIES})


        FIND_PACKAGE(CXSPARSE REQUIRED)
        ADD(MESS_INCLUDE_DIR            ${CXSPARSE_INCLUDE_DIR})
        ADD(MESS_LIBS                   ${CXSPARSE_LIBRARIES})
        ADD(SUITESPARSE_INCLUDE_DIR     ${CXSPARSE_INCLUDE_DIR})
        ADD(SUITESPARSE_LIBRARIES       ${CXSPARSE_LIBRARIES})

        SET(MESS_HAVE_AMD           TRUE)
        SET(MESS_HAVE_COLAMD        TRUE)
        SET(MESS_HAVE_UMFPACK       TRUE)
        SET(MESS_HAVE_CHOLMOD       TRUE)
        SET(MESS_HAVE_CSPARSE       TRUE)
        SET(MESS_USE_SUITESPARSE3   TRUE)
    ENDIF()
    #MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "SuiteSparse not found. You can deactivate SuiteSparse using -DSUITESPARSE=OFF.")
ENDIF()


# SUPERLU
SET(MESS_HAVE_SUPERLU FALSE)
SET(MESS_HAVE_SUPERLU_43 FALSE)
SET(MESS_HAVE_SUPERLU_50 FALSE)
SET(MESS_HAVE_SUPERLU_MT_20 FALSE)
SET(MESS_HAVE_SUPERLU_MT_30 FALSE)

IF(SUPERLU)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for SuperLU Library.")

    FIND_PACKAGE(SUPERLU REQUIRED)

    # get variables from SUPERLU Module
    IF(SUPERLU_FOUND)
        SET(MESS_HAVE_SUPERLU TRUE)
        ADD(MESS_INCLUDE_DIR ${SUPERLU_INCLUDE_DIR})
        ADD(MESS_LIBS ${SUPERLU_LIBRARIES})

        SET(MESS_HAVE_SUPERLU_43        ${SUPERLU_MIN_VERSION_43})
        SET(MESS_HAVE_SUPERLU_50        ${SUPERLU_MIN_VERSION_50})
        SET(MESS_HAVE_SUPERLU_MT_20     ${SUPERLU_MT_MIN_VERSION_20})
        SET(MESS_HAVE_SUPERLU_MT_30     ${SUPERLU_MT_MIN_VERSION_30})

        IF((SUPERLU_MT_MIN_VERSION_20 OR SUPERLU_MT_MIN_VERSION_30) AND OPENMP)
            MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "SuperLU MT Version detected but -DOPENMP=OFF. SuperLU MT needs OPENMP.")
        ENDIF()

    ELSE()
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "SuperLU not found. You can deactivate SuperLU using -DSUPERLU=OFF.")
    ENDIF()
ENDIF()


# Find ARPACK
IF(ARPACK)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for ARPACK Library.")
    FIND_PACKAGE(ARPACK REQUIRED)
    IF(ARPACK_FOUND)
        ADD(MESS_LIBS ${ARPACK_LIBRARIES})
        SET(MESS_HAVE_ARPACK TRUE)
    ENDIF()
ELSE()
    SET(MESS_HAVE_ARPACK FALSE)
ENDIF()


# Find PARDISO
SET(MESS_HAVE_PARDISO FALSE)
IF(PARDISO)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for PARDISO Library.")
    MESSAGE(STATUS "Added PARDISO: ${PARDISO}.")
    ADD(MESS_LIBS  ${PARDISO})
    SET(MESS_HAVE_PARDISO TRUE)
ENDIF()


# MATIO
SET(MESS_HAVE_MATIO FALSE)
IF(MATIO)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for MATIO Library.")
    FIND_PACKAGE(MATIO REQUIRED)
    IF(MATIO_FOUND)
        SET(MESS_HAVE_MATIO TRUE)
        ADD(MESS_LIBS ${MATIO_LIBRARIES})
        ADD(MESS_INCLUDE_DIR ${MATIO_INCLUDE_DIR})
        SET(MESS_HAVE_MATIO TRUE)
    ENDIF()
ENDIF()




# X11 support
SET(MESS_HAVE_X11 FALSE)
SET(MESS_HAVE_XPM TRUE)
IF(X11)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for X11.")
    FIND_PACKAGE(X11 REQUIRED)
    IF(X11_FOUND)
        SET(MESS_HAVE_X11 TRUE)
        ADD(MESS_INCLUDE_DIR  ${X11_INCLUDE_DIR})
        ADD(MESS_LIBS ${X11_LIBRARIES})
    ENDIF()

    IF(X11_Xpm_FOUND)
        SET(MESS_HAVE_XPM TRUE)
        ADD(MESS_INCLUDE_DIR ${X11_Xpm_INCLUDE_PATH})
        ADD(MESS_LIBS ${X11_Xpm_LIB})
    ENDIF()
ENDIF()



# Build the csc utils library
MESSAGE("")
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Configure libcscutils.")
SET(CSCUTILS_FEATURES "inifile;io;threading;image;hardware;ds")
ADD(MESS_INCLUDE_DIR   "${CMAKE_SOURCE_DIR}/libcscutils/include")
ADD(MESS_LIBS cscutils)
ADD_SUBDIRECTORY(libcscutils)

# check for ZLIB AND BZIP2 from libcscutils
SET(MESS_HAVE_ZLIB  ${CSC_IO_ZLIB})
SET(MESS_HAVE_BZIP2 ${CSC_IO_BZIP2})


# Build Tests
IF(TESTS OR LARGETESTS OR HUGETESTS)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Configure tests.")
    INCLUDE(CTest)
    ENABLE_TESTING()
    ADD_SUBDIRECTORY(tests)
    MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make test")

    #check if VALGRIND was found by  CTestConfig.cmake
    IF(VALGRIND)
        IF(MESS_HAVE_OPENMP)
            MESSAGE(AUTHOR_WARNING "We recommend to turn OpenMP support of, for memtest!")
        ENDIF()
        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make memtest")
        ADD_CUSTOM_TARGET(memtest COMMAND ctest --timeout 86400 -D ExperimentalMemCheck)
    ENDIF()
ENDIF()


# Build Tutorials
IF(TUTORIALS)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Configure tutorials.")
    ADD_SUBDIRECTORY(tutorials)
ENDIF()


# Documentation
IF(DOCPDF AND NOT DOC)
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "DOCPDF needs DOC to be switched on")
ENDIF()
IF(DOC)
    SET(DOXYFILE_EXTRA_SOURCES "${CMAKE_SOURCE_DIR}/documents/category.dox")
    SET(DOXYFILE_OUTPUT_DIR "${CMAKE_BINARY_DIR}/doc")
    IF(DOCPDF)
        SET(DOXYFILE_LATEX ON)
        SET(DOXYFILE_GENERATE_LATEX ON)
    ELSE()
        SET(DOXYFILE_LATEX OFF)
        SET(DOXYFILE_GENERATE_LATEX OFF)
    ENDIF()

    SET(DOXYFILE_LATEX_DIR  "latex")
    SET(DOXYFILE_HTML_DIR  "html")
    SET(PROJECT_VERSION ${MESS_VERSION})
    IF(MATHJAX)
        SET(DOXYFILE_MATHJAX "YES")
    ELSE()
        SET(DOXYFILE_MATHJAX "NO")
    ENDIF()

    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for Doxygen.")
    INCLUDE(UseDoxygen REQUIRED)
    IF(DOXYGEN_FOUND)
        MESSAGE(STATUS "Doxygen version: ${DOXYGEN_VERSION}")
        IF(${DOXYGEN_VERSION} VERSION_LESS 1.8.11)
            MESSAGE_COLOR(COLOR_BOLD_RED "Doxygen version is below 1.8.11. Errors during documentation generation can happen.")
        ENDIF()

        #find dot
        FIND_PROGRAM(DOT dot)
        IF(DOT)
            MESSAGE(STATUS "Found dot: ${DOT}")
            SET(DOXYFILE_DOT "YES")
            SET(DOXYGEN_DOT_PATH ${DOT})
        ENDIF()
        FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/doc/html/misc)
        FILE(COPY ${CMAKE_SOURCE_DIR}/documents/doc_style.css DESTINATION  ${CMAKE_BINARY_DIR}/doc/html/misc/)
        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make doc")

        # add target doc to install
        INSTALL(DIRECTORY   "${CMAKE_BINARY_DIR}/doc/html"  DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/cmess PATTERN ".git" EXCLUDE)

    ENDIF()
ENDIF()


SET(MESS_TARGET_LINK_LIBRARIES "${MESS_LIBS}")



# Matlab Support
SET(MESS_MATLAB FALSE)
IF(MATLAB)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for Matlab.")
    FIND_PACKAGE(Matlab REQUIRED)
    SET(MESS_MATLAB TRUE)
    ADD_SUBDIRECTORY(matlab)
ENDIF()


#Octave Support
IF(OCTAVE STREQUAL ON)

    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for Octave.")

    CHECK_ADD_CXX_FLAG("-std=c++11" SUCCESS)

    # check success
    IF(NOT(SUCCESS))
        MESSAGE_COLOR(FATAL_ERROR COLOR_BOLD_RED "Cannot add --std=c++11 flag to the CXX compiler.")
    ENDIF()

    FIND_PACKAGE(Octave REQUIRED)

    SET(MESS_OCTAVE TRUE)
    ADD_SUBDIRECTORY(octave)
ELSE()
    SET(MESS_OCTAVE FALSE)
ENDIF()


# Python Support
IF(PYTHON)
    MESSAGE("")
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for Python.")
    ADD_SUBDIRECTORY(python)
ENDIF()



# search for cppcheck
MESSAGE("")
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search for code checker: cppcheck and cppcheck-htmlreport.")
FIND_PACKAGE(CPPCHECK QUIET)
IF(CPPCHECK_FOUND)
    MESSAGE(STATUS "cppcheck: ${CPPCHECK}")
    MESSAGE(STATUS "cppcheck-htmlreport: ${CPPCHECKHTMLREPORT}")
    FOREACH(OS unix32 unix64)
        ADD_CUSTOM_TARGET("cmess-cppcheck-${OS}"
            COMMAND ${CPPCHECK}
            --force
            --enable=all
            --suppressions-list=${CMAKE_SOURCE_DIR}/suppress_cppcheck.txt
            --platform=${OS}
            --library=posix,gnu,std
            --template=gcc
            --xml-version=2
            -I${CMAKE_BINARY_DIR}
            -I${MESS_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR} 2> ${CMAKE_BINARY_DIR}/cmess_cppcheck_${OS}.xml
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/cppcheck_reports/${OS}
            COMMAND ${CPPCHECKHTMLREPORT}
            --title "C-M.E.S.S. ${OS}"
            --file=${CMAKE_BINARY_DIR}/cmess_cppcheck_${OS}.xml
            --report-dir=${CMAKE_BINARY_DIR}/cppcheck_reports/${OS}
            --source-dir=${CMAKE_SOURCE_DIR}
            )
        MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make cmess-cppcheck-${OS}")
    ENDFOREACH()
ENDIF()


# search for cmake-lint
MESSAGE("")
MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Search cmake-lint tool.")
FIND_PROGRAM(CMAKELINT NAMES cmakelint)
IF(CMAKELINT)
    MESSAGE(STATUS "Found cmakelint programm: ${CMAKELINT}")
    FILE(GLOB_RECURSE MESS_CMAKE_FILES CMake*.txt *.cmake)
    ADD_CUSTOM_TARGET("cmess-cmakelint" VERBATIM COMMAND ${CMAKELINT} --spaces=4 --config=${CMAKE_SOURCE_DIR}/cmakelint.rc ${MESS_CMAKE_FILES})
    MESSAGE_COLOR(STATUS COLOR_BOLD_CYAN "Activate New Target: make cmess-cmakelint")
ENDIF()




# add other directories
ADD_SUBDIRECTORY(lib)
ADD_SUBDIRECTORY(CMakeModules)



#Install / Uninstall
########### Add uninstall target ###############
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")


INSTALL(DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}/include/mess/           DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mess PATTERN ".svn" EXCLUDE)
INSTALL(FILES       ${CMAKE_CURRENT_BINARY_DIR}/include/mess/config.h   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mess)



# some postprocessing
LIST(REMOVE_DUPLICATES MESS_INCLUDE_DIR)
LIST(REMOVE_DUPLICATES MESS_LIBS)

FOREACH(VAR CMAKE_C_FLAGS   CMAKE_C_FLAGS_RELEASE CMAKE_C_FLAGS_DEBUG
                            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_DEBUG
                            CMAKE_Fortran_FLAGS CMAKE_Fortran_FLAGS_DEBUG CMAKE_Fortran_FLAGS_RELEASE)
    STRING(STRIP "${${VAR}}" TMP)
    SET(${VAR} ${TMP})
ENDFOREACH()


MESSAGE(STATUS)
MESSAGE(STATUS "--------------------- Configuration Summary ------------------------------------")
MESSAGE(STATUS "Build Configuration Time (UTC): ${MESS_CONFIG_TIME}")
MESSAGE(STATUS "git branch:                     ${MESS_GIT_BRANCH}")
MESSAGE(STATUS "git commit hash:                ${MESS_GIT_COMMIT_HASH}")
MESSAGE(STATUS "----- Enabled Features -----")
MESSAGE(STATUS "OpenMP support:                OPENMP          = ${MESS_HAVE_OPENMP}")
MESSAGE(STATUS "SuiteSparse and Components:    SUITESPARSE     = ${SUITESPARSE}")
MESSAGE(STATUS "SUPERLU support:               SUPERLU         = ${MESS_HAVE_SUPERLU}")
MESSAGE(STATUS "ARPACK support                 ARPACK          = ${MESS_HAVE_ARPACK}")
MESSAGE(STATUS "X11 Plot support               X11             = ${MESS_HAVE_X11}")
MESSAGE(STATUS "MATIO support                  MATIO           = ${MESS_HAVE_MATIO}")
MESSAGE(STATUS "Python support (Py-M.E.S.S.)   PYTHON          = ${PYTHON} ( ${PYCMESS} )")
MESSAGE(STATUS "MATLAB(R)  addon:              MATLAB          = ${MATLAB}")
MESSAGE(STATUS "Octave     addon:              OCTAVE          = ${OCTAVE}")
MESSAGE(STATUS "Documentation:                 DOC             = ${DOC}")
MESSAGE(STATUS "Documentation (PDF):           DOCPDF          = ${DOCPDF}")
MESSAGE(STATUS "Build tests:                   TESTS           = ${TESTS}")
MESSAGE(STATUS "Build large tests:             LARGETESTS      = ${LARGETESTS}")
MESSAGE(STATUS "Build huge tests:              HUGETESTS       = ${HUGETESTS}")
MESSAGE(STATUS "Build tutorials:               TUTORIALS       = ${TUTORIALS}")
MESSAGE(STATUS "Debug:                         DEBUG           = ${MESS_DEBUG}")
MESSAGE(STATUS "Coverage:                      COVERAGE        = ${COVERAGE}")
MESSAGE(STATUS "BLAS vendor:                   BLA_VENDOR      = ${BLA_VENDOR}")
MESSAGE(STATUS "MATHJAX support:               MATHJAX         = ${MATHJAX}")
MESSAGE(STATUS "PARDISO support:               PARDISO         = ${MESS_HAVE_PARDISO}")
MESSAGE(STATUS "LD_MULTIPLICITY:               LD_MULTIPLICITY = ${LD_MULTIPLICITY}")
MESSAGE(STATUS "SSE2 support:                  SSE2            = ${SSE2}")
MESSAGE(STATUS "SSE4 support:                  SSE4            = ${SSE4}")
MESSAGE("")
MESSAGE("--------------------- Found Libraries and Compiler Options --------------------")
MESSAGE("BLAS_LIBRARIES = ${BLAS_LIBRARIES}")
MESSAGE("LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
MESSAGE("OPENMP = ${OpenMP_C_FLAGS}")

MESSAGE("MESS_LIBS =")
FOREACH(L ${MESS_LIBS})
    MESSAGE(" - ${L} ")
ENDFOREACH()
MESSAGE("MESS_INCLUDE_DIR =")
FOREACH(I ${MESS_INCLUDE_DIR})
    MESSAGE("- ${I}")
ENDFOREACH()
MESSAGE("")
IF(DEBUG)
    MESSAGE("CC FLAGS  = ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
    MESSAGE("CXX FLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
    MESSAGE("FC FLAGS  = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
ELSE()
    MESSAGE("CC FLAGS  = ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
    MESSAGE("CXX FLAGS = ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
    MESSAGE("FC FLAGS  = ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
ENDIF()
MESSAGE("")

MESSAGE("-------------------- Directories and Additional Information -------------------")
SHOW_PROJECT_INFO()
MESSAGE("CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
MESSAGE("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
FOREACH(COMP LIBDIR BINDIR INCLUDEDIR DATADIR)
    MESSAGE("CMAKE_INSTALL_${COMP} = ${CMAKE_INSTALL_${COMP}}")
    MESSAGE("CMAKE_INSTALL_FULL_${COMP} = ${CMAKE_INSTALL_FULL_${COMP}}")
ENDFOREACH()

MESSAGE("Documentation (HTML): ${DOXYFILE_OUTPUT_DIR}/${DOXYFILE_HTML_DIR} ")
IF(DOXYGEN_LATEX)
    MESSAGE("Documentation (PDF): ${DOXYFILE_OUTPUT_DIR}/${DOXYFILE_LATEX_DIR}/refman.pdf")
ENDIF()
MESSAGE("-------------------------------------------------------------------------------")


IF(MESS64)
    MESSAGE("=================================================================================")
    MESSAGE("!!               This is an experimental 64 bit integer build                  !!")
    MESSAGE("=================================================================================")
ENDIF()
IF(NOT SUITESPARSE)
    MESSAGE("==================================================================================")
    MESSAGE("!! MESS is built without SuiteSparse support, this will cause slow computations !!")
    MESSAGE("==================================================================================")
ENDIF()


# check if C-M.E.S.S. is build with OpenMP and OpenBLAS, this can cause slow downs
IF("${BLA_VENDOR}" STREQUAL "OpenBLAS" AND "${OPENMP}" STREQUAL ON)
    MESSAGE("")
    MESSAGE_COLOR(COLOR_BOLD_RED "===================================================================")
    MESSAGE_COLOR(COLOR_BOLD_RED "C-M.E.S.S. is built with OpenMP and uses OpenBLAS as BLAS library!!")
    MESSAGE_COLOR(COLOR_BOLD_RED "This may cause slow downs.")
    MESSAGE_COLOR(COLOR_BOLD_RED "Consider turning OpenMP off, via -DOPENMP=OFF.")
    MESSAGE_COLOR(COLOR_BOLD_RED "OpenMP      = ${OPENMP}")
    MESSAGE_COLOR(COLOR_BOLD_RED "BLAS        = ${BLAS_LIBRARIES}")
    MESSAGE_COLOR(COLOR_BOLD_RED "BLAS Vendor = ${BLA_VENDOR}")
    MESSAGE_COLOR(COLOR_BOLD_RED "LAPACK      = ${LAPACK_LIBRARIES}")
    MESSAGE_COLOR(COLOR_BOLD_RED "===================================================================")
    MESSAGE("")
ENDIF()


# Write the config header
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_BINARY_DIR}/include/mess/config.h)


# Set library location
SET(MESS_LIBRARY_FILENAME $<TARGET_FILE_NAME:mess>)
SET(MESS_LIBRARY "${CMAKE_INSTALL_FULL_LIBDIR}/${MESS_LIBRARY_FILENAME}")
SET(MESS_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})

CONFIGURE_FILE(${MESS_SOURCE_DIR}/mess-config.cmake.in ${CMAKE_BINARY_DIR}/mess-config.cmake @ONLY)

CONFIGURE_FILE(${MESS_SOURCE_DIR}/mess-config-version.cmake.in  ${CMAKE_BINARY_DIR}/mess-config-version.cmake @ONLY)

INSTALL(FILES ${CMAKE_BINARY_DIR}/mess-config.cmake ${CMAKE_BINARY_DIR}/mess-config-version.cmake
    DESTINATION ${CMAKE_INSTALL/DATADIR}/cmake
    COMPONENT Development)



