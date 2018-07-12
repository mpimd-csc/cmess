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


SET(SRC    tutorial_plot
    tutorial_spy
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()



ADD_TEST(tutorial_plot ${EXECUTABLE_OUTPUT_PATH}/tutorial_plot)

ADD_TEST(tutorial_spy ${EXECUTABLE_OUTPUT_PATH}/tutorial_spy
    "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A"
    1024
    1024
    ${CMAKE_CURRENT_BINARY_DIR}/tutorial_spy_plot.bmp
    )



