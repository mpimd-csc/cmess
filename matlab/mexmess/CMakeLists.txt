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


FILE(GLOB_RECURSE SRC *.c)

FOREACH(srcfile ${SRC})
    GET_FILENAME_COMPONENT(srclib ${srcfile} NAME_WE)
    ADD_LIBRARY(${srclib} SHARED ${srcfile})
    ADD_DEFINITIONS(-DMATLAB_MEX_FILE)
    SET_TARGET_PROPERTIES(${srclib}  PROPERTIES LINK_SEARCH_START_STATIC 1)
    SET_TARGET_PROPERTIES(${srclib}  PROPERTIES LINK_SEARCH_END_STATIC 1)
    IF(${CMAKE_VERSION} VERSION_LESS "2.8.8")
        TARGET_LINK_LIBRARIES(${srclib}  matlab_interface mess_matlab_c mess_matlab_f cscutils ${MATLAB_EXT_LIBS} ${MATLAB_LIBS} gomp)
    ELSE()
        TARGET_LINK_LIBRARIES(${srclib}  matlab_interface mess_matlab cscutils ${MATLAB_EXT_LIBS} ${MATLAB_LIBS} gomp)
    ENDIF()

    IF(${CMAKE_VERSION} VERSION_GREATER "2.8.8")
        SET_TARGET_PROPERTIES(${srclib} PROPERTIES POSITION_INDEPENDENT_CODE FALSE)
    ENDIF()
ENDFOREACH()




