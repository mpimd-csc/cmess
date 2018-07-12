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

# SRC VARIABLE IS SET BY SUBDIRECTORIES
ADD_SUBDIRECTORY(MESS_direct)
ADD_SUBDIRECTORY(MESS_enums)
ADD_SUBDIRECTORY(MESS_equation)
ADD_SUBDIRECTORY(MESS_matrix)
ADD_SUBDIRECTORY(MESS_multidirect)
ADD_SUBDIRECTORY(MESS_options)
ADD_SUBDIRECTORY(MESS_status)
ADD_SUBDIRECTORY(MESS_vector)


INCLUDE_DIRECTORIES(${OCTMESS_INCLUDE_DIRS})

ADD_LIBRARY(octave_interface ${SRC})

TARGET_LINK_LIBRARIES(octave_interface mess ${OCTAVE_LIBRARIES})

