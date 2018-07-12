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

SET(SRC
    check_arnoldi_template
    check_arnoldi_template_nrm
    check_eig
    check_eig2
    check_eig_complex
    check_eig_hessenberg_schur
    check_eig_hessenberg_schur_largest_ev
    check_gsign
    check_schur
    check_sign
    check_svd
    check_svd_complex
    check_svd_econ
    check_svd_econ_complex
    )

IF(MESS_HAVE_ARPACK)
    SET(SRC ${SRC}
        check_arpack
        check_arpack_lanczos
        )
ENDIF()

FOREACH(srcfile ${SRC})
    GET_FILENAME_COMPONENT(srcname ${srcfile} NAME_WE)
    ADD_EXECUTABLE(${srcname} ${srcfile})
    TARGET_LINK_LIBRARIES(${srcname} ${MESS_LIBS})
ENDFOREACH()



# Add tests
MESS_ADD_TEST(check_eig  100)
MESS_ADD_TEST(check_eig SUFFIX 2 100)

MESS_ADD_TEST(check_eig_complex         100)
MESS_ADD_TEST(check_schur   SUFFIX 10   10)
MESS_ADD_TEST(check_schur   SUFFIX 50   50)
MESS_ADD_TEST(check_schur   SUFFIX 99   99)
MESS_ADD_TEST(check_svd                 100)
MESS_ADD_TEST(check_svd_complex         100)
MESS_ADD_TEST(check_svd_econ_complex    100 100)
MESS_ADD_TEST(check_sign    ${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx)
MESS_ADD_TEST(check_gsign   ${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx ${CMAKE_SOURCE_DIR}/tests/data/ber_E.mtx)


IF(MESS_HAVE_ARPACK)
    MESS_ADD_TEST(check_arpack  SUFFIX real         ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A)
    MESS_ADD_TEST(check_arpack  SUFFIX complex      ${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx)
    MESS_ADD_TEST(check_arpack_lanczos              ${CMAKE_SOURCE_DIR}/tests/data/SPD_SYM_625.mtx)
ENDIF()


# check arnoldi stuff
MESS_ADD_TEST(check_arnoldi_template        SUFFIX iss_A            ${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx)
MESS_ADD_TEST(check_arnoldi_template        SUFFIX heat-cont        ${CMAKE_SOURCE_DIR}/tests/data/heat-cont/A.mtx)
MESS_ADD_TEST(check_arnoldi_template        SUFFIX rail371          ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.A)
MESS_ADD_TEST(check_arnoldi_template        SUFFIX rand50           ${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx)
MESS_ADD_TEST(check_arnoldi_template        SUFFIX rand50           ${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx)
MESS_ADD_TEST(check_arnoldi_template        SUFFIX NSE_RE500_A      ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_500_lvl1_A.mtx)

MESS_ADD_TEST(check_arnoldi_template_nrm    SUFFIX iss_A    ${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx            40)
MESS_ADD_TEST(check_arnoldi_template_nrm    SUFFIX rail371  ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.A  40)
MESS_ADD_TEST(check_arnoldi_template_nrm    SUFFIX rand50   ${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx        40)


# check hessenberg stuff
FOREACH(DIMS 10 20 25 33 64)
    FOREACH(CPX 0 1)
        FOREACH(KTRIU -1 0 1)
            MESS_ADD_TEST(check_eig_hessenberg_schur            SUFFIX "dim${DIMS}_cpx${CPX}_ktriu${KTRIU}" ${DIMS} ${CPX} ${KTRIU})
            MESS_ADD_TEST(check_eig_hessenberg_schur_largest_ev SUFFIX "dim${DIMS}_cpx${CPX}_ktriu${KTRIU}" ${DIMS} ${CPX} ${KTRIU})
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

