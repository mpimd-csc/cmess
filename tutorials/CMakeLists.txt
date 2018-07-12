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


INCLUDE_DIRECTORIES(${MESS_INCLUDE_DIR})
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

SET(MESS_LIBS mess ${MESS_LIBS})


ADD_SUBDIRECTORY(bt)
ADD_SUBDIRECTORY(direct)
ADD_SUBDIRECTORY(eig)
ADD_SUBDIRECTORY(h2)
ADD_SUBDIRECTORY(io)
ADD_SUBDIRECTORY(itsolver)
ADD_SUBDIRECTORY(lrcf_adi)
ADD_SUBDIRECTORY(matrix)
ADD_SUBDIRECTORY(matgen)
ADD_SUBDIRECTORY(plot)
