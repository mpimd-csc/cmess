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


# set CMESS_SRC and CMESS_FSRC in that scope
ADD_SUBDIRECTORY(contrib)
ADD_SUBDIRECTORY(direct)
ADD_SUBDIRECTORY(dynsys)
ADD_SUBDIRECTORY(easyfrontend)
ADD_SUBDIRECTORY(eigen)
ADD_SUBDIRECTORY(formats)
ADD_SUBDIRECTORY(graph)
ADD_SUBDIRECTORY(io)
ADD_SUBDIRECTORY(itsolver)
ADD_SUBDIRECTORY(lrcf_adi)
ADD_SUBDIRECTORY(matrix)
ADD_SUBDIRECTORY(misc)
ADD_SUBDIRECTORY(plot)
ADD_SUBDIRECTORY(precond)
ADD_SUBDIRECTORY(reorder)
ADD_SUBDIRECTORY(vector)
ADD_SUBDIRECTORY(extblas)



INCLUDE_DIRECTORIES(${MESS_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

ADD_LIBRARY(mess ${CMESS_SRC} ${CMESS_FSRC})

# compile seperate matlab library if MATLAB is on
IF(MATLAB STREQUAL ON)
    IF(${CMAKE_VERSION} VERSION_LESS "2.8.8")
        ADD_LIBRARY(mess_matlab_f STATIC ${CMESS_FSRC})
        ADD_LIBRARY(mess_matlab_c STATIC ${CMESS_SRC})
        IF(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
            SET_TARGET_PROPERTIES(mess_matlab_f PROPERTIES COMPILE_FLAGS "-fdefault-integer-8 -DMESS64")
        ELSEIF(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
            SET_TARGET_PROPERTIES(mess_matlab_f PROPERTIES COMPILE_FLAGS "-i8 -DMESS64")
        ENDIF()
        SET_TARGET_PROPERTIES(mess_matlab_c PROPERTIES COMPILE_FLAGS "-DMESS64 -DMESS_MATLAB")
    ELSE()
        ADD_LIBRARY(mess_matlab_f OBJECT ${CMESS_FSRC})
        ADD_LIBRARY(mess_matlab_c OBJECT ${CMESS_SRC})
        IF(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
            SET_TARGET_PROPERTIES(mess_matlab_f PROPERTIES COMPILE_FLAGS "-fdefault-integer-8 -DMESS64")
        ELSEIF(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
            SET_TARGET_PROPERTIES(mess_matlab_f PROPERTIES COMPILE_FLAGS "-i8 -DMESS64")
        ENDIF()
        SET_TARGET_PROPERTIES(mess_matlab_c PROPERTIES COMPILE_FLAGS "-DMESS64 -DMESS_MATLAB")
        ADD_LIBRARY(mess_matlab STATIC $<TARGET_OBJECTS:mess_matlab_f> $<TARGET_OBJECTS:mess_matlab_c>)
    ENDIF()
ENDIF()



TARGET_LINK_LIBRARIES(mess ${MESS_LIBS})

SET(GENERIC_LIB_VERSION "1.0.0")
SET(GENERIC_LIB_SOVERSION "0")
SET_TARGET_PROPERTIES(mess PROPERTIES VERSION ${GENERIC_LIB_VERSION} SOVERSION ${GENERIC_LIB_SOVERSION})


INSTALL(TARGETS mess
    # IMPORTANT: Add the foo library to the "export-set"
    EXPORT MESSLibraryDepends
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT bin
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT shlib
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT stlib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT dev)









