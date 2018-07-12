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

INCLUDE_DIRECTORIES(${MESS_INCLUDE_DIR} ${MATLAB_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
FILE(GLOB_RECURSE CINTERFACE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
ADD_LIBRARY(matlab_interface STATIC
    ${CINTERFACE_SRC}
    )


