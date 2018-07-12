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
MESS_ADD_TEST(check_eps)
MESS_ADD_TEST(check_getelement  "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx")
MESS_ADD_TEST(check_getelement SUFFIX 2 "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx")



FOREACH(ROWS 1 10 20 25 33 64)
    FOREACH(COLS 1 10 20 25 33 64)
        MESS_ADD_TEST(check_setelement SUFFIX "${ROWS}_${COLS}" ${ROWS} ${COLS})
    ENDFOREACH()
ENDFOREACH()

FOREACH(ROWS 1 3 5 7 10)
    FOREACH(COLS 1 3 5 7 10)
        MESS_ADD_TEST(check_tridiag SUFFIX "${ROWS}_${COLS}" ${ROWS} ${COLS})
    ENDFOREACH()
ENDFOREACH()


MESS_ADD_TEST(check_any)
MESS_ADD_TEST(check_set_get_col_row)
MESS_ADD_TEST(check_colscale)
MESS_ADD_TEST(check_trace)
MESS_ADD_TEST(check_fro_inn)
MESS_ADD_TEST(check_colaxpy2)
MESS_ADD_TEST(check_colvecaxpy)
MESS_ADD_TEST(check_coldotE)
MESS_ADD_TEST(check_rowsub)
MESS_ADD_TEST(check_colsub)
MESS_ADD_TEST(check_sub)
MESS_ADD_TEST(check_mgs)
MESS_ADD_TEST(check_mgs_inplace)
MESS_ADD_TEST(check_norm "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx" "${CMAKE_SOURCE_DIR}/tests/data/randc50x50.mtx" "${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx")
MESS_ADD_TEST(check_transpose)
MESS_ADD_TEST(check_proj_sym)
MESS_ADD_TEST(check_mgsadd)
MESS_ADD_TEST(check_coldotcE)
MESS_ADD_TEST(check_orth)
MESS_ADD_TEST(check_gaxpy)
MESS_ADD_TEST(check_condest "${CMAKE_SOURCE_DIR}/tests/data/rand50x50.mtx" 1.962107663241329e+03)
MESS_ADD_TEST(check_addcol2p)
MESS_ADD_TEST(check_matrix_colaxpy)
MESS_ADD_TEST(check_matrix_null)
MESS_ADD_TEST(check_matrix_sum)
MESS_ADD_TEST(check_matrix_gbiorth)
MESS_ADD_TEST(check_matrix_biorth)
MESS_ADD_TEST(check_decomp)


# RANK
MESS_ADD_TEST(check_rank "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx"  10)


# NNORMEST
SET(VALS "10;25;100")
SET(REFS "30;200;3000")
SET(CNT 0)
FOREACH(VAL ${VALS})
    LIST(GET REFS ${CNT} REF)
    MESS_ADD_TEST(check_dynorm2 SUFFIX  ${VAL} ${VAL} ${REF})
    MESS_ADD_TEST(check_dynorm2t SUFFIX ${VAL} ${VAL} ${REF})
    MESS_ADD_TEST(check_dynormf SUFFIX  ${VAL} ${VAL} ${REF})
    MATH(EXPR CNT "${CNT}+1")
ENDFOREACH()


MESS_ADD_TEST(check_dynorm_indefinite2)
MESS_ADD_TEST(check_dynorm_indefinitef)
MESS_ADD_TEST(check_mulnormf2)

# Check mess_matrix_cat
FOREACH(I RANGE 0 15 1)
    MESS_ADD_TEST(check_cat SUFFIX ${I} ${I})
ENDFOREACH()



MESS_ADD_TEST(check_colnorm)
MESS_ADD_TEST(check_coldot)
MESS_ADD_TEST(check_coldotc)
MESS_ADD_TEST(check_colvecdot)
MESS_ADD_TEST(check_colvecdotc)
MESS_ADD_TEST(check_diag)


# Matrix Multiply Test
MESS_ADD_TEST(check_mul)
MESS_ADD_TEST(check_mul_tall)
MESS_ADD_TEST(check_mul_small)

#MVP Tests y = Ax
ADD_TEST(check_mvpA25_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 0)
ADD_TEST(check_mvpB25_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/B25.mtx" 0)
ADD_TEST(check_mvpC25_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/C25.mtx" 0)
ADD_TEST(check_mvpFA_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A" 0)
ADD_TEST(check_mvpFE_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E" 0)
ADD_TEST(check_mvpFB_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B" 0)
ADD_TEST(check_mvpFC_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C" 0)
ADD_TEST(check_mvpR1_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx" 0)
ADD_TEST(check_mvpR2_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_2.mtx" 0)
ADD_TEST(check_mvpR3_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x1668_1.mtx" 0)
ADD_TEST(check_mvpRC_NONE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 0)

#MVP Tests y = A^Tx
ADD_TEST(check_mvpA25_TRANSPOSE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 1)
ADD_TEST(check_mvpB25_TRANSPOSE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/B25.mtx" 1)
ADD_TEST(check_mvpC25_TRANSPOSE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/C25.mtx" 1)
ADD_TEST(check_mvpFA_TRANSPOSE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A" 1)
ADD_TEST(check_mvpFE_TRANSPOSE ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E" 1)
ADD_TEST(check_mvpFB_TRANSPOSE  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B" 1)
ADD_TEST(check_mvpFC_TRANSPOSE  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C" 1)
ADD_TEST(check_mvpR1_TRANSPOSE  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx" 1)
ADD_TEST(check_mvpR2_TRANSPOSE  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_2.mtx" 1)
ADD_TEST(check_mvpR3_TRANSPOSE  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x1668_1.mtx" 1)
ADD_TEST(check_mvpRC_TRANSPOSE  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 1)

#MVP Tests y = A^Hx
ADD_TEST(check_mvpA25_HERMITIAN ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/A25.mtx" 2)
ADD_TEST(check_mvpB25_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/B25.mtx" 2)
ADD_TEST(check_mvpC25_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/FDM-CD/C25.mtx" 2)
ADD_TEST(check_mvpFA_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A" 2)
ADD_TEST(check_mvpFE_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E" 2)
ADD_TEST(check_mvpFB_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B" 2)
ADD_TEST(check_mvpFC_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C" 2)
ADD_TEST(check_mvpR1_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_1.mtx" 2)
ADD_TEST(check_mvpR2_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x10_2.mtx" 2)
ADD_TEST(check_mvpR3_HERMITIAN  ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/rand10x1668_1.mtx" 2)
ADD_TEST(check_mvpRC_HERMITIAN ${EXECUTABLE_OUTPUT_PATH}/check_mvp "${CMAKE_SOURCE_DIR}/tests/data/randc20x20.mtx" 2)

#A+B -> B
MESS_ADD_TEST(check_add)


# Biorth
FOREACH(COLS 10 15 20)
    FOREACH(ADDROWS 0 1 2 3 9 10 100 1000)
        MATH(EXPR ROWS "${COLS} +${ADDROWS}")
        ADD_TEST(check_biorth_${ROWS}_${COLS} ${EXECUTABLE_OUTPUT_PATH}/check_biorth ${ROWS} ${COLS})
    ENDFOREACH()
ENDFOREACH()

# Bandwidth
MESS_ADD_TEST(check_bandwidth)

# ELIMINATE ZEROS
FOREACH(ROWS 1 10 20 25 33 64)
    FOREACH(COLS 1 10 20 25 33 64)
        ADD_TEST(check_eliminate_zeros_${ROWS}_${COLS} ${EXECUTABLE_OUTPUT_PATH}/check_eliminate_zeros ${ROWS} ${COLS})
    ENDFOREACH()
ENDFOREACH()

# CHECK ONE VALUE
FOREACH(ROWS 1 10 20 25 33 64)
    FOREACH(COLS 1 10 20 25 33 64)
        ADD_TEST(check_one_value_${ROWS}_${COLS} ${EXECUTABLE_OUTPUT_PATH}/check_one_value ${ROWS} ${COLS})
        ADD_TEST(check_scalem_${ROWS}_${COLS} ${EXECUTABLE_OUTPUT_PATH}/check_scalem ${ROWS} ${COLS})
    ENDFOREACH()
ENDFOREACH()



