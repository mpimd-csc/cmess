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

SET(_CMESS_FSRC
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/dgegsgndet.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/dgegsgnnone.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/dgesgndet.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/dgesgnfro.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/dgesgnnone.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/zgegsgndet.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/zgegsgnnone.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/zgesgndet.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/zgesgnfro.f
    ${CMAKE_CURRENT_SOURCE_DIR}/sgn/zgesgnnone.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dgeadd.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dgefdf.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dgefnr.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dgetrp.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzdotu.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzgbtrs.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzgeadd.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzgemm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzgemv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzger.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzgetrs.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dzpotrs.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dztbsv.f
    ${CMAKE_CURRENT_SOURCE_DIR}/dztrsm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/zddotc.f
    ${CMAKE_CURRENT_SOURCE_DIR}/zddotu.f
    ${CMAKE_CURRENT_SOURCE_DIR}/zdgemm.f
    ${CMAKE_CURRENT_SOURCE_DIR}/zgeadd.f
    ${CMAKE_CURRENT_SOURCE_DIR}/zgefdf.f
    ${CMAKE_CURRENT_SOURCE_DIR}/zgefnr.f
    )

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET(_CMESS_FSRC ${_CMESS_FSRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/dtrsyl.f
        ${CMAKE_CURRENT_SOURCE_DIR}/dlasy2.f
        )
    MESSAGE(STATUS "Enable workaround for DTRSYL/DLASY2 on Apple.")
ENDIF()

SET(CMESS_FSRC ${CMESS_FSRC} ${_CMESS_FSRC} PARENT_SCOPE)
