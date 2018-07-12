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

SHOW_VARIABLE(OCTAVE_EXECUTABLE)
SHOW_VARIABLE(OCTAVE_INCLUDE_DIRS)
SHOW_VARIABLE(OCTAVE_LIBRARIES)
SHOW_VARIABLE(OCTAVE_OCTINTERP_LIBRARY)
SHOW_VARIABLE(OCTAVE_OCTAVE_LIBRARY)
SHOW_VARIABLE(OCTAVE_CRUFT_LIBRARY)
SHOW_VARIABLE(OCTAVE_VERSION_STRING)
SHOW_VARIABLE(OCTAVE_MAJOR_VERSION)
SHOW_VARIABLE(OCTAVE_MINOR_VERSION)
SHOW_VARIABLE(OCTAVE_PATCH_VERSION)
SHOW_VARIABLE(OCTAVE_OCT_FILE_DIR)
SHOW_VARIABLE(OCTAVE_OCT_LIB_DIR)
SHOW_VARIABLE(OCTAVE_ROOT_DIR)
SHOW_VARIABLE(OCTAVE_M_SITE_DIR)
SHOW_VARIABLE(OCTAVE_OCT_SITE_DIR)


SET(OCTMESS_INCLUDE_DIRS ${MESS_INCLUDE_DIR} ${OCTAVE_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/cpp_interface/include)



ADD_SUBDIRECTORY(cpp_interface)
ADD_SUBDIRECTORY(octmess)





# copy/configure *.m files to build directory
FILE(GLOB_RECURSE OCTFILES ${CMAKE_SOURCE_DIR}/octave/*.m)
FOREACH(OCT_IN ${OCTFILES})
    STRING(REPLACE ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} OCT_OUT ${OCT_IN})
    CONFIGURE_FILE(${OCT_IN} ${OCT_OUT} @ONLY)
ENDFOREACH()

