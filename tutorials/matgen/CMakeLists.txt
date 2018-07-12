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


SET(SRC    tutorial_gen_fdm_cd
    tutorial_gen_laplace
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()

ADD_TEST(tutorial_gen_fdm_cd ${EXECUTABLE_OUTPUT_PATH}/tutorial_gen_fdm_cd  500 "convection.mtx")
ADD_TEST(tutorial_gen_laplace ${EXECUTABLE_OUTPUT_PATH}/tutorial_gen_laplace  500  "laplace.mtx")

