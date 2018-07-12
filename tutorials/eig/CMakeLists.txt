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
    tutorial_eig
    tutorial_eigs
    tutorial_gschur
    tutorial_schur
    tutorial_sign
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()


ADD_TEST(tutorial_eig               ${EXECUTABLE_OUTPUT_PATH}/tutorial_eig                  "${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx")
ADD_TEST(tutorial_eigs              ${EXECUTABLE_OUTPUT_PATH}/tutorial_eigs                 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx")
ADD_TEST(tutorial_sign              ${EXECUTABLE_OUTPUT_PATH}/tutorial_sign                 "${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx")
ADD_TEST(tutorial_gschur_real       ${EXECUTABLE_OUTPUT_PATH}/tutorial_gschur               "${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx" "${CMAKE_SOURCE_DIR}/tests/data/ber_E.mtx" 1)
ADD_TEST(tutorial_gschur_complex    ${EXECUTABLE_OUTPUT_PATH}/tutorial_gschur               "${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx" "${CMAKE_SOURCE_DIR}/tests/data/ber_E.mtx" 2)
ADD_TEST(tutorial_schur_real        ${EXECUTABLE_OUTPUT_PATH}/tutorial_schur                "${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx" 1)
ADD_TEST(tutorial_schur_cpx         ${EXECUTABLE_OUTPUT_PATH}/tutorial_schur                "${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx" 2)



