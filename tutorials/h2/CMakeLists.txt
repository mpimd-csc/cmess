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
    tutorial_h2norm
    tutorial_irka_biorth
    tutorial_tsia
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()


ADD_TEST(tutorial_h2norm ${EXECUTABLE_OUTPUT_PATH}/tutorial_h2norm
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    )

ADD_TEST(tutorial_irka_biorth ${EXECUTABLE_OUTPUT_PATH}/tutorial_irka_biorth
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    10
    200
    1e-9
    1
    )

ADD_TEST(tutorial_tsia ${EXECUTABLE_OUTPUT_PATH}/tutorial_tsia
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    10
    200
    1e-9
    1
    )

ADD_TEST(tutorial_tsiag ${EXECUTABLE_OUTPUT_PATH}/tutorial_tsia
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C"
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
    10
    200
    1e-9
    1
    )

