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

SET(SRC
    check_cimport
    check_cimport_c2d
    check_convert
    check_mvpcall_matrix
    check_perm
    check_real_imag_part
    )

IF(MESS_HAVE_CHOLMOD)
    SET(SRC ${SRC} check_cholmod)
ENDIF()

FOREACH(srcfile ${SRC})
    GET_FILENAME_COMPONENT(srcname ${srcfile} NAME_WE)
    ADD_EXECUTABLE(${srcname} ${srcfile})
    TARGET_LINK_LIBRARIES(${srcname} ${MESS_LIBS})
ENDFOREACH()


# Add tests
MESS_ADD_TEST(check_cimport)
MESS_ADD_TEST(check_cimport_c2d)
MESS_ADD_TEST(check_convert SUFFIX dense    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
MESS_ADD_TEST(check_convert SUFFIX sparse   "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx")

IF(MESS_HAVE_CHOLMOD)
    ADD_TEST(check_cholmod ${EXECUTABLE_OUTPUT_PATH}/check_cholmod)
ENDIF()

MESS_ADD_TEST(check_perm)


# CHECK MVP CALL
FOREACH(ROWS 1 10 20 25 33 64)
    FOREACH(COLS 1 10 20 25 33 64)
        MESS_ADD_TEST(check_mvpcall_matrix SUFFIX ${ROWS}_${COLS} ${ROWS} ${COLS})
        MESS_ADD_TEST(check_real_imag_part SUFFIX ${ROWS}_${COLS} ${ROWS} ${COLS})
    ENDFOREACH()
ENDFOREACH()







