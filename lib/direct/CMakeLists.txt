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

SET(_CMESS_SRC    ${_CMESS_SRC}
    ${CMAKE_CURRENT_SOURCE_DIR}/backslash.c
    ${CMAKE_CURRENT_SOURCE_DIR}/care.c
    ${CMAKE_CURRENT_SOURCE_DIR}/cholfactor.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dgeqp3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/interface/direct.c
    ${CMAKE_CURRENT_SOURCE_DIR}/interface/multidirect.c
    ${CMAKE_CURRENT_SOURCE_DIR}/multisolver/multilu_adavanced.c
    ${CMAKE_CURRENT_SOURCE_DIR}/pinv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/res2/lyap_res2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/res2/stein_res2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/res2/res2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/res2/sylv_res2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/banded.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/bicgstab.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/cholesky.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/cholmod_chol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/csparse.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/csparse_chol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/glyap3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/hessenberg_lu.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/lapack.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/lapack_qr.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/lu.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/newlu.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/sylvester_dense.c
    ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/sylvester_sparsedense.c
    ${CMAKE_CURRENT_SOURCE_DIR}/splsolve.c
    ${CMAKE_CURRENT_SOURCE_DIR}/trisolve.c
    )

IF(MESS_HAVE_SUPERLU)
    SET(_CMESS_SRC ${_CMESS_SRC} ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/superlu.c)
ENDIF()


IF(MESS_HAVE_UMFPACK)
    SET(_CMESS_SRC   ${_CMESS_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/multisolver/multilu_umfpack.c
        ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/umfpack_control.c
        ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/umfpack.c
        )
ENDIF()

IF(MESS_HAVE_MKLPARDISO)
    SET(_CMESS_SRC ${_CMESS_SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/singlesolver/mklpardiso.c
        )
ENDIF()

SET(CMESS_SRC ${CMESS_SRC} ${_CMESS_SRC} PARENT_SCOPE)


SET(_CMESS_FSRC
    ${CMAKE_CURRENT_SOURCE_DIR}/lu_reuse.f90
    )

SET(CMESS_FSRC ${CMESS_FSRC} ${_CMESS_FSRC} PARENT_SCOPE)




