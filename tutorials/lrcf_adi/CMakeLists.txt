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


SET(SRC    tutorial_lradi
    tutorial_lrnm
    tutorial_lyap
    tutorial_care
    tutorial_lrnm-g
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()





ADD_TEST(tutorial_lradi ${EXECUTABLE_OUTPUT_PATH}/tutorial_lradi
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    )

ADD_TEST(tutorial_lradig ${EXECUTABLE_OUTPUT_PATH}/tutorial_lradi
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    )

ADD_TEST(tutorial_lrnm ${EXECUTABLE_OUTPUT_PATH}/tutorial_lrnm
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    )

ADD_TEST(tutorial_lrnmg ${EXECUTABLE_OUTPUT_PATH}/tutorial_lrnm
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    )

ADD_TEST(tutorial_lyap ${EXECUTABLE_OUTPUT_PATH}/tutorial_lyap
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    )

ADD_TEST(tutorial_lyapg ${EXECUTABLE_OUTPUT_PATH}/tutorial_lyap
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    )


##tutorial_care
ADD_TEST(tutorial_care ${EXECUTABLE_OUTPUT_PATH}/tutorial_care
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    )

ADD_TEST(tutorial_careg ${EXECUTABLE_OUTPUT_PATH}/tutorial_care
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    )


ADD_TEST(tutorial_lrnm-g ${EXECUTABLE_OUTPUT_PATH}/tutorial_lrnm-g
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    )

