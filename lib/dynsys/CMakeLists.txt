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


FILE(GLOB _CMESS_SRC *.c)
FILE(GLOB_RECURSE _CMESS_SRC_BT bt/*.c)
SET(_CMESS_SRC_H2 ${_CMESS_SRC_H2}
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/h2_norm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/h2error.c
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/irka_biorth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/irka_para_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/statusopt.c
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/tsia.c
    ${CMAKE_CURRENT_SOURCE_DIR}/h2/tsiag.c
    )

SET(_CMESS_SRC ${_CMESS_SRC} ${_CMESS_SRC_BT} ${_CMESS_SRC_H2})

SET(CMESS_SRC ${CMESS_SRC} ${_CMESS_SRC} PARENT_SCOPE)





