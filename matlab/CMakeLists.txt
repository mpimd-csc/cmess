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


SET(CMAKE_CXX_COMPILER ${MATLAB_MEX})
SET(CMAKE_C_COMPILER ${MATLAB_MEX})

# Suffix and Prefix of the output target file
SET(CMAKE_SHARED_LIBRARY_SUFFIX ${MATLAB_EXT})  # set suffix to .mexa64
SET(CMAKE_SHARED_LIBRARY_PREFIX)                # remove the "lib" prefix

# Variables controlling the build-phrase, reset flags
SET(CMAKE_C_FLAGS)
SET(CMAKE_C_FLAGS_DEBUG)
SET(CMAKE_C_FLAGS_RELEASE)
SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS)
SET(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
SET(CMAKE_SHARED_LINKER_FLAGS)


# if matlab version is higher or equal than 9.5 we use the interleaved api
# https://de.mathworks.com/help/matlab/matlab_external/do-i-need-to-upgrade-my-mex-files-to-use-interleaved-complex.html
IF(MATLAB_MEX_INTERLEAVED_COMPLEX_API)
    MESSAGE_COLOR(STATUS COLOR_BOLD_GREEN "Use the Interleaved Complex API for MEX-M.E.S.S.")
    SET(CMAKE_C_FLAGS "-R2018a CFLAGS='$$CFLAGS -std=gnu89' -DMESS64")
ELSE()
    SET(CMAKE_C_FLAGS "-largeArrayDims CFLAGS='$$CFLAGS -std=gnu89' -DMESS64")
ENDIF()

IF(${CMAKE_VERSION} VERSION_LESS "2.8.9")
    SET(CMAKE_SHARED_LIBRARY_C_FLAGS)       # remove the -fPIC option. mex does not accept the "-fPIC" option
ENDIF()

# optimization flags are not supported
STRING(REGEX REPLACE "-O[0-9]?" "" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})

#
# Some Problems occur with old gcc compiler during linking and the pthread
# library. See also
# http://stackoverflow.com/questions/24768622/imagemagick-pthread-h-multiple-definition
# Thats why we changed to gnu89 standard.
#
IF(${CMAKE_VERSION} VERSION_LESS "2.9.0")
    SET(CMAKE_C_COMPILE_OBJECT
        "<CMAKE_C_COMPILER> <DEFINES> <FLAGS> -outdir <OBJECT_DIR> -c <SOURCE>; mv <OBJECT_DIR>/$$(basename <SOURCE> .c).o <OBJECT>"
        )
ELSE()
    SET(CMAKE_C_COMPILE_OBJECT
        "<CMAKE_C_COMPILER> <DEFINES> <FLAGS> <INCLUDES> -outdir <OBJECT_DIR> -c <SOURCE>; mv <OBJECT_DIR>/$$(basename <SOURCE> .c).o <OBJECT>"
        )
ENDIF()

SET(
    CMAKE_C_CREATE_SHARED_LIBRARY
    "<CMAKE_C_COMPILER> -cxx <LINK_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS> -output <TARGET> <OBJECTS> <LINK_LIBRARIES>"
    )


SET(MATLAB_LIBS "mwcholmod;mwumfpack;mwcolamd;mwamd;mwarpack;mwlapack;mwblas;rt")
SET(MATLAB_EXT_LIBS ${CXSPARSE_LIBRARIES})

# Variables controlling the installation RPATH
# SET( CMAKE_INSTALL_RPATH "\$ORIGIN" )
# CMake will reset RPATH at the installation phase, so we need to specify CMAKE_INSTALL_RPATH



ADD_SUBDIRECTORY(c_interface)
INCLUDE_DIRECTORIES(${MESS_INCLUDE_DIR}  ${MATLAB_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_CURRENT_SOURCE_DIR}/c_interface)
ADD_SUBDIRECTORY(mexmess)


# copy/configure *.m files to build directory
FILE(GLOB_RECURSE MATLABFILES ${CMAKE_SOURCE_DIR}/matlab/*.m)
FOREACH(MATLAB_IN ${MATLABFILES})
    STRING(REPLACE ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} MATLAB_OUT ${MATLAB_IN})
    CONFIGURE_FILE(${MATLAB_IN} ${MATLAB_OUT} @ONLY)
ENDFOREACH()

