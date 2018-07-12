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

SET(SRC     tutorial_direct_solvers
            tutorial_backslash
            tutorial_qr
            tutorial_sylvester
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()




ADD_TEST(tutorial_direct_solvers    ${EXECUTABLE_OUTPUT_PATH}/tutorial_direct_solvers)
ADD_TEST(tutorial_backslash         ${EXECUTABLE_OUTPUT_PATH}/tutorial_backslash
    "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"
    )


ADD_TEST(tutorial_qr         ${EXECUTABLE_OUTPUT_PATH}/tutorial_qr
    "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"
    )



ADD_TEST(tutorial_sylvester_sparsedense_standard            ${EXECUTABLE_OUTPUT_PATH}/tutorial_sylvester
    "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"
    )

ADD_TEST(tutorial_sylvester_sparsedense_halfgeneralized     ${EXECUTABLE_OUTPUT_PATH}/tutorial_sylvester
    "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A11.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/E11.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"
    )

ADD_TEST(tutorial_sylvester_sparsedense_generalized     ${EXECUTABLE_OUTPUT_PATH}/tutorial_sylvester
    "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A11.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/E11.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"
    )

