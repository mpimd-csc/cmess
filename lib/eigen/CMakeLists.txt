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


SET(_CMESS_SRC ${_CMESS_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/arnoldi.c
    ${CMAKE_CURRENT_SOURCE_DIR}/arnoldi_template.c
    ${CMAKE_CURRENT_SOURCE_DIR}/eig.c
    ${CMAKE_CURRENT_SOURCE_DIR}/eig_hessenberg.c
    ${CMAKE_CURRENT_SOURCE_DIR}/lanczos_template_nrm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/schur.c
    ${CMAKE_CURRENT_SOURCE_DIR}/schur_to_evd.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sign.c
    ${CMAKE_CURRENT_SOURCE_DIR}/svd.c
    )

IF(MESS_HAVE_ARPACK)
    SET(_CMESS_SRC ${_CMESS_SRC} ${CMAKE_CURRENT_SOURCE_DIR}/arpack.c)
ENDIF()

SET(CMESS_SRC ${CMESS_SRC} ${_CMESS_SRC} PARENT_SCOPE)

