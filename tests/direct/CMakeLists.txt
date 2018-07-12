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

#############################################linsolve#############################################################
#real matrix, real rhs
MESS_ADD_TEST(linsolve SUFFIX newlu  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)

IF(MESS_HAVE_CSPARSE)
    MESS_ADD_TEST(linsolve SUFFIX csparse  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
ENDIF()

IF( MESS_HAVE_UMFPACK )
    MESS_ADD_TEST(linsolve SUFFIX umfpack  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 2)
ENDIF()


MESS_ADD_TEST(linsolve SUFFIX bicgstab  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 3)
MESS_ADD_TEST(linsolve SUFFIX lapacklu  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 4)
MESS_ADD_TEST(linsolve SUFFIX lapackqr  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 5)
MESS_ADD_TEST(linsolve SUFFIX cholesky  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   6)


IF( MESS_HAVE_CSPARSE )
    ADD_TEST(linsolve_chol_csparse       ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   7)
ENDIF()

IF(MESS_HAVE_CHOLMOD)
    ADD_TEST(linsolve_cholesky_cholmod   ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx" 8)
ENDIF()

ADD_TEST(linsolve_banded           ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 9)

IF(MESS_HAVE_SUPERLU)
    ADD_TEST(linsolve_superlu         ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 10)
ENDIF()

IF(MESS_HAVE_MKLPARDISO)
    ADD_TEST(linsolve_mklpardiso        ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 11)
ENDIF()


#complex matrix, real rhs
ADD_TEST(linsolvec_newlu            ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)
ADD_TEST(linsolvec_bicgstab         ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 3)
ADD_TEST(linsolvec_lapacklu         ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 4)
ADD_TEST(linsolvec_lapackqr         ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 5)
IF(MESS_HAVE_CSPARSE)
    ADD_TEST(linsolvec_csparse      ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)
ENDIF()
IF(MESS_HAVE_UMFPACK)
    ADD_TEST(linsolvec_umfpack      ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 2)
ENDIF()
ADD_TEST(linsolvec_banded           ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 9)
IF(MESS_HAVE_SUPERLU)
    ADD_TEST(linsolvec_superlu      ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 10)
ENDIF()

IF(MESS_HAVE_MKLPARDISO)
    ADD_TEST(linsolvec_mklpardiso       ${EXECUTABLE_OUTPUT_PATH}/linsolve  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 11)
ENDIF()

#############################################trisolve#############################################################
ADD_TEST(check_trisolve   ${EXECUTABLE_OUTPUT_PATH}/check_trisolve)

#############################################backslah#############################################################
ADD_TEST(check_backslash               ${EXECUTABLE_OUTPUT_PATH}/check_backslash "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
#ADD_TEST(check_backslash_rectangular   ${EXECUTABLE_OUTPUT_PATH}/check_backslash "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1 )

#############################################linsolve_multi#############################################################
ADD_TEST(linsolve_multi       ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
ADD_TEST(linsolve_multic       ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)
IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolve_multi_umfpack   ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
    ADD_TEST(linsolve_multi_umfpackc   ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)
ENDIF()

#############################################linsolvem_multi#############################################################
ADD_TEST(linsolvem_multi       ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
ADD_TEST(linsolvem_multic       ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)
IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolvem_multi_umfpack   ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
    ADD_TEST(linsolvem_multi_umfpackc   ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)
ENDIF()

#############################################linsolvem#############################################################
#real matrix, real rhs
ADD_TEST(linsolvem_newlu               ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
ADD_TEST(linsolvem_bicgstab             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 3)
ADD_TEST(linsolvem_lapacklu             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 4)
ADD_TEST(linsolvem_lapackqr             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 5)
ADD_TEST(linsolvem_cholesky             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   6)

IF( MESS_HAVE_CSPARSE )
    ADD_TEST(linsolvem_csparse           ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
    ADD_TEST(linsolvem_chol_csparse       ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   7)
ENDIF()

IF(MESS_HAVE_CHOLMOD)
    ADD_TEST(linsolvem_chol_csparse       ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"    8)
ENDIF()

IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolvem_umfpack           ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 2)
ENDIF()

ADD_TEST(linsolvem_banded               ${EXECUTABLE_OUTPUT_PATH}/linsolvem   "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 9)

IF(MESS_HAVE_SUPERLU)
    ADD_TEST(linsolvem_superlu            ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 10)
ENDIF()

IF(MESS_HAVE_MKLPARDISO)
    ADD_TEST(linsolvem_mklpardiso           ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 11)
ENDIF()


#complex matrix, real rhs
ADD_TEST(linsolvemc_newlu               ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)
ADD_TEST(linsolvemc_bicgstab             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 3)
ADD_TEST(linsolvemc_lapacklu             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 4)
ADD_TEST(linsolvemc_lapackqr             ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 5)

IF( MESS_HAVE_CSPARSE )
    ADD_TEST(linsolvemc_csparse         ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)
ENDIF()

IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolvemc_umfpack          ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 2)
ENDIF()

ADD_TEST(linsolvemc_banded               ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 9)

IF(MESS_HAVE_SUPERLU)
    ADD_TEST(linsolvemc_superlu            ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 10)
ENDIF()

IF(MESS_HAVE_MKLPARDISO)
    ADD_TEST(linsolvemc_mklpardiso          ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 11)
ENDIF()


#real matrix, complex rhs
ADD_TEST(linsolvem_newlu_complexrhs        ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -1)
ADD_TEST(linsolvem_bicgstab_complexrhs     ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -4)
ADD_TEST(linsolvem_lapacklu_complexrhs     ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -5)
ADD_TEST(linsolvem_lapackqr_complexrhs     ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -6)
ADD_TEST(linsolvem_cholesky_complexrhs     ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"    -7)

IF(MESS_HAVE_CSPARSE)
    ADD_TEST(linsolvem_csparse_complexrhs      ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -2)
ENDIF()

IF(MESS_HAVE_CHOLMOD)
    ADD_TEST(linsolvem_chol_csparse_complexrhs  ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   -8)
ENDIF()

IF(MESS_HAVE_UMFPACK)
    ADD_TEST(linsolvem_umfpack_complexrhs      ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -3)
ENDIF()

ADD_TEST(linsolvem_banded_complexrhs       ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -10)

IF(MESS_HAVE_SUPERLU)
    ADD_TEST(linsolvem_superlu_complexrhs    ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -11)
ENDIF()

IF(MESS_HAVE_MKLPARDISO)
    ADD_TEST(linsolvem_mklpardiso_complexrhs      ${EXECUTABLE_OUTPUT_PATH}/linsolvem  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" -12)
ENDIF()


#############################################inverse#############################################################
ADD_TEST(inverse_dense_rand   ${EXECUTABLE_OUTPUT_PATH}/inverse 0)
ADD_TEST(inverse_dense_fdm    ${EXECUTABLE_OUTPUT_PATH}/inverse 0 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx")
ADD_TEST(inverse_dense_iss    ${EXECUTABLE_OUTPUT_PATH}/inverse 0 "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx")
ADD_TEST(inverse_dense_complex  ${EXECUTABLE_OUTPUT_PATH}/inverse 0 "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx")
ADD_TEST(inverse_chol       ${EXECUTABLE_OUTPUT_PATH}/inverse 4 "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx")

IF(MESS_HAVE_CSPARSE)
    ADD_TEST(inverse_csparse_rand     ${EXECUTABLE_OUTPUT_PATH}/inverse 1)
    ADD_TEST(inverse_csparse_fdm      ${EXECUTABLE_OUTPUT_PATH}/inverse 1 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx")
    # ISS matrix is not well suited for csparse
    ADD_TEST(inverse_csparse_r10      ${EXECUTABLE_OUTPUT_PATH}/inverse 1 "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
    ADD_TEST(inverse_csparse_complex   ${EXECUTABLE_OUTPUT_PATH}/inverse 1 "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx")
    ADD_TEST(inverse_csparse_chol     ${EXECUTABLE_OUTPUT_PATH}/inverse 5 "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx")
ENDIF()
IF(MESS_HAVE_CHOLMOD)
    ADD_TEST(inverse_cholmod_chol     ${EXECUTABLE_OUTPUT_PATH}/inverse 6 "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx")
ENDIF()
IF(MESS_HAVE_UMFPACK)
    ADD_TEST(inverse_umfpack_rand     ${EXECUTABLE_OUTPUT_PATH}/inverse 2)
    ADD_TEST(inverse_umfpack_fdm      ${EXECUTABLE_OUTPUT_PATH}/inverse 2 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx")
    ADD_TEST(inverse_umfpack_iss      ${EXECUTABLE_OUTPUT_PATH}/inverse 2 "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx")
    ADD_TEST(inverse_umfpack_complex   ${EXECUTABLE_OUTPUT_PATH}/inverse 2 "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx")
ENDIF()

ADD_TEST(inverse_newlu_rand     ${EXECUTABLE_OUTPUT_PATH}/inverse 3)
ADD_TEST(inverse_newlu_fdm      ${EXECUTABLE_OUTPUT_PATH}/inverse 3 "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx")
ADD_TEST(inverse_newlu_iss      ${EXECUTABLE_OUTPUT_PATH}/inverse 3 "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx")
ADD_TEST(inverse_newlu_complex     ${EXECUTABLE_OUTPUT_PATH}/inverse 3 "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx")


#############################################getLU#############################################################
ADD_TEST(getLU_newlu      ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
IF(MESS_HAVE_UMFPACK)
    ADD_TEST(getLU_umfpack    ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 2)
ENDIF()
#ADD_TEST( getLU_bicgstab    ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 3)
ADD_TEST(getLU_lapacklu    ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 4)
ADD_TEST(getLU_lapackqr    ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 5)
ADD_TEST(getLU_cholesky    ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   6)
IF(MESS_HAVE_CSPARSE)
    ADD_TEST(getLU_csparse    ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
    ADD_TEST(getLU_csparse_chol ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   7)
ENDIF()
IF(MESS_HAVE_CHOLMOD)
    ADD_TEST(getLU_cholmod_chol ${EXECUTABLE_OUTPUT_PATH}/check_getLU  "${CMAKE_SOURCE_DIR}/tests/data/rand_spd.mtx"   8)
ENDIF()

ADD_TEST(linsolve_multi ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolve_multi_umfpack ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
ENDIF()
ADD_TEST(linsolve_multic ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)
IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolve_multi_umfpackc ${EXECUTABLE_OUTPUT_PATH}/linsolve_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)
ENDIF()

ADD_TEST(linsolvem_multi ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolvem_multi_umfpack ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
ENDIF()
ADD_TEST(linsolvem_multic ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)
IF( MESS_HAVE_UMFPACK )
    ADD_TEST(linsolvem_multi_umfpackc ${EXECUTABLE_OUTPUT_PATH}/linsolvem_multi  "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)
ENDIF()


#############################################Dense LyapSolve#############################################################

ADD_TEST(glyap_iss          ${EXECUTABLE_OUTPUT_PATH}/glyap_test    "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx"              "${CMAKE_SOURCE_DIR}/tests/data/iss/E.mtx"              "${CMAKE_SOURCE_DIR}/tests/data/iss/B.mtx"              "${CMAKE_SOURCE_DIR}/tests/data/iss/C.mtx")
ADD_TEST(glyap_rail         ${EXECUTABLE_OUTPUT_PATH}/glyap_test    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.A"    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.E"    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.B"    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.C")

#ADD_TEST(glyap_iss          ${EXECUTABLE_OUTPUT_PATH}/gstein_test    "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx"              "${CMAKE_SOURCE_DIR}/tests/data/iss/E.mtx"              "${CMAKE_SOURCE_DIR}/tests/data/iss/B.mtx"              "${CMAKE_SOURCE_DIR}/tests/data/iss/C.mtx")
ADD_TEST(gstein_rail        ${EXECUTABLE_OUTPUT_PATH}/gstein_test   "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.A"    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.E"    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.B"    "${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.C")


#############################################Sparse-Dense Sylvester solver#############################################################
MESS_ADD_TEST(sylvester_sparsedense SUFFIX standard         "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx"         "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
MESS_ADD_TEST(sylvester_sparsedense SUFFIX halfgeneralized  "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A11.mtx"     "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/E11.mtx"   "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
MESS_ADD_TEST(sylvester_sparsedense SUFFIX generalized      "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A11.mtx"     "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/bips98_606/E11.mtx"     "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")

IF(LARGETESTS STREQUAL ON OR  HUGETESTS STREQUAL ON)

    MESS_ADD_TEST(sylvester_sparsedense SUFFIX standard_filter          "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A" "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
    MESS_ADD_TEST(sylvester_sparsedense SUFFIX halfgeneralized_filter   "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A" "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E"
                                                                        "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
    MESS_ADD_TEST(sylvester_sparsedense SUFFIX generalized_filter       "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A" "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_2.mtx"
                                                                        "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E" "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
ENDIF()

#############################################Dense Sylvester solver#############################################################
MESS_ADD_TEST(sylvester_dense SUFFIX standard 0)
MESS_ADD_TEST(sylvester_dense SUFFIX generalized 1)

#############################################determinant#############################################################
MESS_ADD_TEST(check_det  SUFFIX 1    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"  0.007946720688968)
MESS_ADD_TEST(check_det  SUFFIX 2    "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_2.mtx" -0.002254676768426)
MESS_ADD_TEST(check_det  SUFFIX 3    "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"   -583972.017622754792683)

#############################################pinv####################################################################
MESS_ADD_TEST(check_pinv    SUFFIX real_N   "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"  0)
MESS_ADD_TEST(check_pinv    SUFFIX real_T   "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"  1)
MESS_ADD_TEST(check_pinv    SUFFIX cpx_N    "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 0)
MESS_ADD_TEST(check_pinv    SUFFIX cpx_T    "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" 1)

#############################################check_lapackqr####################################################################
MESS_ADD_TEST(check_lapackqr SUFFIX real_overdetermined     "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"      0)
MESS_ADD_TEST(check_lapackqr SUFFIX real_underdetermined    "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx"      1)
MESS_ADD_TEST(check_lapackqr SUFFIX cpx_overdetermined      "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx"     0)
MESS_ADD_TEST(check_lapackqr SUFFIX cpx_underdetermined     "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx"     1)

#############################################check_dgepf####################################################################
MESS_ADD_TEST(check_dgeqp3  SUFFIX 1  "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")
MESS_ADD_TEST(check_dgeqp3  SUFFIX 2  "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_2.mtx")
MESS_ADD_TEST(check_dgeqp3  SUFFIX 3  "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx")

#############################################check_direct_care####################################################################
MESS_ADD_TEST(check_direct_care     "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/iss/E.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/iss/B.mtx"
    "${CMAKE_SOURCE_DIR}/tests/data/iss/C.mtx")

#############################################check_hessenberg####################################################################
FOREACH(DIM 1 2 4 7 10)
    FOREACH(CPX 0 1)
        MESS_ADD_TEST(check_hessenberg SUFFIX "dim${DIM}_cpx${CPX}" ${DIM} ${CPX})
    ENDFOREACH()
ENDFOREACH()

