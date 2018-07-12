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
    GET_FILENAME_COMPONENT(srcname ${srcfile} NAME_WE)
    ADD_EXECUTABLE(${srcname} ${srcfile})
    TARGET_LINK_LIBRARIES(${srcname} ${MESS_LIBS})
ENDFOREACH()



# Add tests
MESS_ADD_TEST(check_vectorio)
MESS_ADD_TEST(check_matrixio)
MESS_ADD_TEST(check_fputs_block)

IF( MESS_HAVE_ZLIB )
    MESS_ADD_TEST(check_vectorio    SUFFIX gzip test.mtx.gz)
    MESS_ADD_TEST(check_fputs_block SUFFIX gzip test.gz)
ENDIF()

IF( MESS_HAVE_BZ2 )
    MESS_ADD_TEST(check_vectorio        SUFFIX bzip2 test.mtx.bz2)
    MESS_ADD_TEST(check_fputs_block     SUFFIX bzip2 test.bz2)
ENDIF()



