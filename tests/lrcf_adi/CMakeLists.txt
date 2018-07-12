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
    check_lradi_unsym
    check_lradi_glyap
    check_ber_sign_fac
    check_lrnm_unsym
    check_lrnm_g
    check_lrnm_so
    check_lradi_dae1
    check_lradi_dae2
    check_linesearch
    check_dae2_handles_para
    check_dae2_handles_adi
    check_lrnm_dae1
    check_lrnm_dae2
    check_lradi_so
    check_so1_handles
    check_so2_handles
    check_lyap_sign_fac
    check_dense_nm
    check_dense_nm_tridiag
    )

FOREACH(srcfile ${SRC})
    ADD_EXECUTABLE(${srcfile} ${srcfile}.c)
    TARGET_LINK_LIBRARIES(${srcfile} ${MESS_LIBS})
ENDFOREACH()


#######################SETUP PARAMETER LIST FOR SMALL, LARGE AND HUGE TESTS#################
IF(LARGETESTS STREQUAL ON)
    # parameter combinations if LARGETESTS=ON
    SET(DIMS 25)                        # dimension of the fdm generated matrices
    SET(PARATYPES_UNSYM 0 1 3 4)        # possible parameter types for unsymmetric systems
    SET(PARATYPES_SYM 0 1 2 3 4)        # possible parameter type for symmetric systems
    SET(EQNTYPES 0 1)                   # possible equation types (MESS_OP_NONE, MESS_OP_TRANSPOSE)
    SET(NM_GPS 0  5)                    # Galerkin projection frequency in LRNM
    SET(SINGLESHIFTS 0)                 # Single Shifts in NM ON / OFF
    SET(RE_NUMBERS 100 300)             # Reynolds numbers for NSE EXAMPLES
    SET(MEMORY_USAGES 0 1)              # memory usage variant, LOW/MEDIUM/HIGH (only in dae2)
    SET(USE_PLUS 0 1)                   # add instead of subtract quadratic term of Riccati equation (only in pcare)

    # 1<->MESS_LU, 3<->UMFPACK, 4<->SUPERLU, 5<->CSPARSELU, 6<->BANDED, 7<->MKLPARDISO
    SET(DIRECT_SOLVERS 1 6)

    IF(MESS_HAVE_UMFPACK)
        SET(DIRECT_SOLVERS ${DIRECT_SOLVERS} 3)
    ENDIF()

    IF(MESS_HAVE_SUPERLU)
        SET(DIRECT_SOLVERS ${DIRECT_SOLVERS} 4)
    ENDIF()

    IF(MESS_HAVE_CSPARSE)
        SET(DIRECT_SOLVERS ${DIRECT_SOLVERS} 5)
    ENDIF()

    # 0<->MESSLU_MULTI, 1<->UMFPACK_MULTI, 2<->BICGSTAB
    SET(MDIRECT_SOLVERS 0)
    IF(MESS_HAVE_UMFPACK)
        SET(MDIRECT_SOLVERS ${MDIRECT_SOLVERS} 1)
    ENDIF()

    # LINESEARCH
    SET(LINESEARCH 0 1)

    # RESIDUAL METHOD, 0<->INDEFINITE, 1<->SPECTRAL
    SET(RESIDUAL_METHODS 0)

ELSE()
    # parameter combinations if LARGETESTS=OFF
    SET(DIMS 25)                       # dimension of the fdm generated matrices
    SET(PARATYPES_UNSYM 0)              # possible parameter types for unsymmetric systems
    SET(PARATYPES_SYM  2)               # possible parameter types for symmetric systems
    SET(EQNTYPES 0 1)                   # possible equation types (MESS_OP_NONE, MESS_OP_TRANSPOSE)
    SET(NM_GPS  0)                      # Galerkin projection frequency in LRNM
    SET(SINGLESHIFTS 0)                 # Single Shifts in NM ON / OFF
    SET(RE_NUMBERS 100)                 # Reynolds numbers for NSE EXAMPLES
    SET(MEMORY_USAGES 1)                # memory usage variant, LOW/MEDIUM/HIGH (only in dae2)
    SET(USE_PLUS 0 1)                   # add instead of subtract quadratic term of Riccati equation (only in pcare)

    # 1<->MESS_LU, 3<->UMFPACK, 4<->SUPERLU, 5<->CSPARSELU, 6<->BANDED
    IF(MESS_HAVE_UMFPACK AND MESS_HAVE_CSPARSE)
        SET(DIRECT_SOLVERS 1 3)
    ELSE()
        SET(DIRECT_SOLVERS 1)
    ENDIF()

    SET(MDIRECT_SOLVERS 0)

    # LINESEARCH
    SET(LINESEARCH 0 1)

    # RESIDUAL METHOD, 0<->INDEFINITE, 1<->SPECTRAL
    SET(RESIDUAL_METHODS 0 1)

ENDIF()

IF(HUGETESTS STREQUAL ON)
    # parameter combinations if HUGETESTS=ON
    SET(DIMS  25 50 100)                # dimension of the fdm generated matrices
    SET(PARATYPES_UNSYM 0 1 3 4)        # possible parameter types for unsymmetric systems
    SET(PARATYPES_SYM 0 1 2 3 4)        # possible parameter type for symmetric systems
    SET(EQNTYPES 0 1)                   # possible equation types (MESS_OP_NONE, MESS_OP_TRANSPOSE)
    SET(NM_GPS 0 1 5)                   # Galerkin projection frequency in LRNM
    SET(SINGLESHIFTS 0 1)               # Single Shifts in NM ON / OFF
    SET(RE_NUMBERS 100 300 500)         # Reynolds numbers for NSE EXAMPLES
    SET(MEMORY_USAGES 0 1 2)            # memory usage variant, LOW/MEDIUM/HIGH (only in dae2)
    SET(USE_PLUS 0 1)                   # add instead of subtract quadratic term of Riccati equation (only in pcare)

    # 1<->MESS_LU, 3<->UMFPACK, 4<->SUPERLU, 5<->CSPARSELU, 6<->BANDED, 7<->MKLPARDISO
    SET(DIRECT_SOLVERS 1 6)

    IF(MESS_HAVE_UMFPACK)
        SET(DIRECT_SOLVERS ${DIRECT_SOLVERS} 3)
    ENDIF()

    IF(MESS_HAVE_SUPERLU)
        SET(DIRECT_SOLVERS ${DIRECT_SOLVERS} 4)
    ENDIF()

    IF(MESS_HAVE_CSPARSE)
        SET(DIRECT_SOLVERS ${DIRECT_SOLVERS} 5)
    ENDIF()

    IF(MESS_HAVE_MKLPARDISO)
        SET(DIRECT_SOLVERS ${DIRECT_SOVLERS} 7)
    ENDIF()

    # 0<->MESSLU_MULTI, 1<->UMFPACK_MULTI
    SET(MDIRECT_SOLVERS 0)
    IF(MESS_HAVE_UMFPACK)
        SET(MDIRECT_SOLVERS ${MDIRECT_SOLVERS} 1)
    ENDIF()

    # LINESEARCH
    SET(LINESEARCH 0 1)

    # RESIDUAL METHOD, 0<->INDEFINITE, 1<->SPECTRAL
    SET(RESIDUAL_METHODS 0 1)

ENDIF()



########################## LINE SEARCH FUNCTION ###################################
# Check line search functions
ADD_TEST(check_linesearch       ${EXECUTABLE_OUTPUT_PATH}/check_linesearch
    ${CMAKE_SOURCE_DIR}/tests/data/linesearch/W.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/linesearch/W_new.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/linesearch/DeltaK.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/linesearch/DeltaK_new.mtx
    3.79295290418053102493286132812500e+09
    7.27920354517517034786955321123969e-03
    5.56871212201042994856834411621094e+07
    1.42857830908468940833699889481068e+02
    -2.94314807746533453464508056640625e+08
    3.60602436042776943736498651560396e+01
    1.06274076082054391712006236048182e+00
    )


########################### New ADI Code ##########################################
#Unsymmetric
FOREACH(DIM ${DIMS})
    FOREACH(PARATYPE ${PARATYPES_UNSYM})
        FOREACH(DIRECT ${DIRECT_SOLVERS})
            FOREACH(MDIRECT ${MDIRECT_SOLVERS})
                FOREACH(RESIDUAL_METHOD ${RESIDUAL_METHODS})

                    ADD_TEST(lradi-unsym-${DIM}-${PARATYPE}-direct${DIRECT}-mdirect${MDIRECT}-resmeth${RESIDUAL_METHOD}
                        ${EXECUTABLE_OUTPUT_PATH}/check_lradi_unsym
                        ${DIM}
                        ${PARATYPE}
                        ${DIRECT}
                        ${MDIRECT}
                        ${RESIDUAL_METHOD}
                        )
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

#Generalized, filter is symmetric
FOREACH(PARATYPE ${PARATYPES_SYM})
    FOREACH(DIRECT ${DIRECT_SOLVERS})
        FOREACH(MDIRECT ${MDIRECT_SOLVERS})
            FOREACH(RESIDUAL_METHOD ${RESIDUAL_METHODS})

                ADD_TEST(lradi-glyap-filter-${PARATYPE}-direct${DIRECT}-mdirect${MDIRECT}-resmeth${RESIDUAL_METHOD}
                    ${EXECUTABLE_OUTPUT_PATH}/check_lradi_glyap
                    ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A
                    ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E
                    ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B
                    ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C
                    ${PARATYPE}
                    ${DIRECT}
                    ${MDIRECT}
                    ${RESIDUAL_METHOD}
                    )

            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

# Riccati
#lrnm unsym
FOREACH(DIM ${DIMS})
    FOREACH(PARATYPE ${PARATYPES_UNSYM})
        FOREACH(NM_GP ${NM_GPS})
            FOREACH(SINGLESHIFT ${SINGLESHIFTS})
                FOREACH(DIRECT ${DIRECT_SOLVERS})
                    FOREACH(MDIRECT ${MDIRECT_SOLVERS})
                        FOREACH(RESIDUAL_METHOD ${RESIDUAL_METHODS})

                            FOREACH(LINE ${LINESEARCH})
                                ADD_TEST(lrnm-unsym-${DIM}-${PARATYPE}-${SINGLESHIFT}-nmgp${NM_GP}-direct${DIRECT}-mdirect${MDIRECT}-line${LINE}-resmeth${RESIDUAL_METHOD}
                                    ${EXECUTABLE_OUTPUT_PATH}/check_lrnm_unsym
                                    ${DIM}
                                    ${PARATYPE}
                                    ${NM_GP}
                                    ${SINGLESHIFT}
                                    ${DIRECT}
                                    ${MDIRECT}
                                    ${LINE}
                                    ${RESIDUAL_METHOD}
                                    )
                            ENDFOREACH()
                        ENDFOREACH()
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

#lrnm filter
FOREACH(PARATYPE ${PARATYPES_SYM})
    FOREACH(NM_GP ${NM_GPS})
        FOREACH(SINGLESHIFT ${SINGLESHIFTS})
            FOREACH(DIRECT ${DIRECT_SOLVERS})
                FOREACH(MDIRECT ${MDIRECT_SOLVERS})
                    FOREACH(RESIDUAL_METHOD ${RESIDUAL_METHODS})

                        FOREACH(LINE ${LINESEARCH})
                            ADD_TEST(lrnm-filter-${PARATYPE}-${SINGLESHIFT}-nmgp${NM_GP}-direct${DIRECT}-mdirect${MDIRECT}-line${LINE}-resmeth${RESIDUAL_METHOD}
                                ${EXECUTABLE_OUTPUT_PATH}/check_lrnm_g
                                ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.A
                                ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.E
                                ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.B
                                ${CMAKE_SOURCE_DIR}/tests/data/filter2D/filter2D.C
                                ${PARATYPE}
                                ${NM_GP}
                                ${SINGLESHIFT}
                                ${DIRECT}
                                ${MDIRECT}
                                ${LINE}
                                ${RESIDUAL_METHOD}
                                )
                        ENDFOREACH()
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

#############################################second order######################################

ADD_TEST(check_so1_handles  ${EXECUTABLE_OUTPUT_PATH}/check_so1_handles
    ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/M_301.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/D_301.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/K_301.mtx)
ADD_TEST(check_so2_handles  ${EXECUTABLE_OUTPUT_PATH}/check_so2_handles
    ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/M_301.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/D_301.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/K_301.mtx)

#SO1/2 Lyapunov symmetric matrices different shift parameter strategies
FOREACH(VARIANT 1 2)
    FOREACH(EQNTYPE ${EQNTYPES})
        FOREACH(PARATYPE ${PARATYPES_UNSYM})
            IF(EQNTYPE EQUAL 0)
                SET(EQNTYPENAME "N")
            ELSE()
                SET(EQNTYPENAME "T")
            ENDIF()

            FOREACH(DIRECT ${DIRECT_SOLVERS})
                FOREACH(MDIRECT ${MDIRECT_SOLVERS})

                    ADD_TEST(lradi-triple-chain-so${VARIANT}-${PARATYPE}-${EQNTYPENAME}-direct${DIRECT}-mdirect${MDIRECT}
                        ${EXECUTABLE_OUTPUT_PATH}/check_lradi_so
                        ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/M_301.mtx
                        ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/D_301.mtx
                        ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/K_301.mtx
                        ${EQNTYPE}
                        ${VARIANT}
                        ${PARATYPE}
                        ${DIRECT}
                        ${MDIRECT}
                        )
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

# SO Riccati symmetric matrices different shift parameters strategies
FOREACH(VARIANT 1 2)
    FOREACH(EQNTYPE 0 1)
        FOREACH(PARATYPE ${PARATYPES_UNSYM})
            FOREACH(NM_GP ${NM_GPS})
                FOREACH(SINGLESHIFT ${SINGLESHIFTS})
                    IF(EQNTYPE EQUAL 0)
                        SET(EQNTYPENAME "N")
                    ELSE()
                        SET(EQNTYPENAME "T")
                    ENDIF()

                    FOREACH(DIRECT ${DIRECT_SOLVERS})

                        # do not use adaptive shifts for SO2 variant, it does not work very well
                        IF((VARIANT EQUAL 2) AND (PARATYPE EQUAL 4))
                            SET(VARIANT EQUAL 1)
                        ENDIF()



                        FOREACH(LINE ${LINESEARCH})
                            ADD_TEST(lrnm-triple-chain-so${VARIANT}-${PARATYPE}-${SINGLESHIFT}-${EQNTYPENAME}-nmgp${NM_GP}-direct${DIRECT}-line${LINE}
                                ${EXECUTABLE_OUTPUT_PATH}/check_lrnm_so
                                ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/M_301.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/D_301.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/TripleChain/K_301.mtx
                                ${EQNTYPE}
                                ${VARIANT}
                                ${PARATYPE}
                                ${NM_GP}
                                ${SINGLESHIFT}
                                ${DIRECT}
                                ${LINE}
                                )
                        ENDFOREACH()
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()


# Sign Function
ADD_TEST(lyap-sign-fac
    ${EXECUTABLE_OUTPUT_PATH}/check_lyap_sign_fac
    ${CMAKE_SOURCE_DIR}/tests/data/iss/A.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/iss/B.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/iss/E.mtx
    )

# bernoulli solver based on sign function
ADD_TEST(ber-sign-fac
    ${EXECUTABLE_OUTPUT_PATH}/check_ber_sign_fac
    ${CMAKE_SOURCE_DIR}/tests/data/ber_A.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/ber_B.mtx
    ${CMAKE_SOURCE_DIR}/tests/data/ber_E.mtx
    )


# INDEX 2 DAE Functions NSE LRADI
# N   <-> (opt->type==MESS_OP_NONE)
# T   <-> (opt->type==MESS_OP_TRANPOSE)
# RE *  <-> Reynoldsnumber = *
# Feed  <-> with initial feedback
#CHECK DAE 2 HANDLES LOW MEMORY_USAGE / HIGH MEMORY_USAGE
FOREACH(RE ${RE_NUMBERS})
    FOREACH(MEM ${MEMORY_USAGES})
        ADD_TEST(check_dae2_handles_para_mem${MEM}
            ${EXECUTABLE_OUTPUT_PATH}/check_dae2_handles_para
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_M.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_A.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_G.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_B.mtx
            ${MEM}
            )
        ADD_TEST(check_dae2_handles_adi_mem${MEM}
            ${EXECUTABLE_OUTPUT_PATH}/check_dae2_handles_adi
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_M.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_A.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_G.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_B.mtx
            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_MM.mtx
            ${MEM}
            )
    ENDFOREACH()
ENDFOREACH()


######dae2 lradi
FOREACH(RE ${RE_NUMBERS})
    FOREACH(EQNTYPE ${EQNTYPES})
        FOREACH(PARATYPE ${PARATYPES_UNSYM})
            FOREACH(MEMORYUSAGE ${MEMORY_USAGES})
                IF(EQNTYPE EQUAL 0)
                    SET(EQNTYPENAME "N")
                    SET(RHS ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_B.mtx)
                    IF((${RE} GREATER 300) OR (${RE} EQUAL 300))
                        SET(FEED ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_Feed0.mtx)
                    ELSE()
                        SET(FEED "")
                    ENDIF()
                ELSE()
                    SET(EQNTYPENAME "T")
                    SET(RHS ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_C.mtx)
                    IF((${RE} GREATER 300) OR (${RE} EQUAL 300))
                        SET(FEED ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_Feed1.mtx)
                    ELSE()
                        SET(FEED "")
                    ENDIF()
                ENDIF()

                FOREACH(DIRECT ${DIRECT_SOLVERS})
                    FOREACH(MDIRECT ${MDIRECT_SOLVERS})


                        ADD_TEST(lradi_dae2_${EQNTYPENAME}_RE${RE}_mem${MEMORYUSAGE}-${PARATYPE}-direct${DIRECT}-mdirect${MDIRECT}
                            ${EXECUTABLE_OUTPUT_PATH}/check_lradi_dae2
                            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_M.mtx
                            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_A.mtx
                            ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_G.mtx
                            ${RHS}
                            ${MEMORYUSAGE}
                            ${PARATYPE}
                            ${DIRECT}
                            ${MDIRECT}
                            ${FEED}
                            )
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()




#lrnm dae2
FOREACH(RE ${RE_NUMBERS})
    FOREACH(EQNTYPE ${EQNTYPES})
        FOREACH(PARATYPE ${PARATYPES_UNSYM})
            FOREACH(MEMORYUSAGE ${MEMORY_USAGES})
                FOREACH(NM_GP ${NM_GPS})
                    FOREACH(SINGLESHIFT ${SINGLESHIFTS})
                        FOREACH(LINE ${LINESEARCH})
                            IF(EQNTYPE EQUAL 0)
                                SET(EQNTYPENAME "N")
                                IF((${RE} GREATER 300) OR (${RE} EQUAL 300))
                                    SET(FEED ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_Feed1.mtx)
                                ELSE()
                                    SET(FEED "")
                                ENDIF()
                            ELSE()
                                SET(EQNTYPENAME "T")
                                IF((${RE} GREATER 300) OR (${RE} EQUAL 300))
                                    SET(FEED ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_Feed0.mtx)
                                ELSE()
                                    SET(FEED "")
                                ENDIF()
                            ENDIF()

                            #ONLY ADD TEST IF RE IS SMALLER THAN 300 and NEWTON GALERKIN PROJECTION IS OFF, BECAUSE IT IS NOT WORKING
                            IF((${RE} LESS 300) OR ( ${NM_GP} EQUAL 0 ))
                                ADD_TEST(lrnm_dae2_${EQNTYPENAME}_RE${RE}_mem${MEMORYUSAGE}-${PARATYPE}-${SINGLESHIFT}-nmgp${NM_GP}-line${LINE}
                                    ${EXECUTABLE_OUTPUT_PATH}/check_lrnm_dae2
                                    ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_M.mtx
                                    ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_A.mtx
                                    ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_G.mtx
                                    ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_B.mtx
                                    ${CMAKE_SOURCE_DIR}/tests/data/NSE/NSE_RE_${RE}_lvl1_C.mtx
                                    ${EQNTYPE}
                                    ${MEMORYUSAGE}
                                    ${PARATYPE}
                                    ${NM_GP}
                                    ${SINGLESHIFT}
                                    ${LINE}
                                    ${FEED}
                                    )
                            ENDIF()
                        ENDFOREACH()
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()


###INDEX 1 Systems
FOREACH(EQNTYPE ${EQNTYPES})
    FOREACH(PARATYPE ${PARATYPES_UNSYM})
        IF(EQNTYPE EQUAL 0)
            SET(EQNTYPENAME "N")
            SET(RHS ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/B.mtx)
        ELSE()
            SET(EQNTYPENAME "T")
            SET(RHS ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/C.mtx)
            #ELSE()
        ENDIF()

        FOREACH(DIRECT ${DIRECT_SOLVERS})
            FOREACH(MDIRECT ${MDIRECT_SOLVERS})

                # deactivate banded direct solver
                IF(DIRECT EQUAL 6)
                    SET(DIRECT 0)
                ENDIF()


                ADD_TEST(lradi_dae1_${EQNTYPENAME}-${PARATYPE}-direct${DIRECT}-mdirect${MDIRECT}
                    ${EXECUTABLE_OUTPUT_PATH}/check_lradi_dae1
                    ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/E11.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A11.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A12.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A21.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A22.mtx
                    ${RHS}
                    ${PARATYPE}
                    ${DIRECT}
                    ${MDIRECT}
                    )
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()

FOREACH(EQNTYPE ${EQNTYPES})
    FOREACH(PARATYPE ${PARATYPES_UNSYM})
        FOREACH(NM_GP ${NM_GPS})
            FOREACH(SINGLESHIFT ${SINGLESHIFTS})
                FOREACH(LINE ${LINESEARCH})

                    IF(EQNTYPE EQUAL 0)
                        SET(EQNTYPENAME "N")
                    ELSE()
                        SET(EQNTYPENAME "T")
                    ENDIF()
                    FOREACH(DIRECT ${DIRECT_SOLVERS})
                        FOREACH(MDIRECT ${MDIRECT_SOLVERS})

                            # deactivate banded direct solver
                            IF(DIRECT EQUAL 6)
                                SET(DIRECT 0)
                            ENDIF()


                            ADD_TEST(lrnm_dae1_${EQNTYPENAME}-${PARATYPE}-${SINGLESHIFT}-nmgp${NM_GP}-direct${DIRECT}-mdirect${MDIRECT}-line${LINE}

                                ${EXECUTABLE_OUTPUT_PATH}/check_lrnm_dae1
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/E11.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A11.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A12.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A21.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/A22.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/B.mtx
                                ${CMAKE_SOURCE_DIR}/tests/data/bips98_606/C.mtx
                                ${EQNTYPE}
                                ${PARATYPE}
                                ${SINGLESHIFT}
                                ${NM_GP}
                                ${DIRECT}
                                ${MDIRECT}
                                ${LINE}
                                )
                        ENDFOREACH()
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()


### pcare
FOREACH(EQNTYPE ${EQNTYPES})
    FOREACH(PLUS ${USE_PLUS})
        FOREACH(LINE ${LINESEARCH})

            IF(EQNTYPE EQUAL 0)
                SET(EQNTYPENAME "N")
            ELSE()
                SET(EQNTYPENAME "T")
            ENDIF()

            # tridiagonal matrices
            ADD_TEST(dense_pcare_nm_tridiag_${EQNTYPENAME}-line${LINE}

                ${EXECUTABLE_OUTPUT_PATH}/check_dense_nm_tridiag
                100
                0
                ${LINE}
                ${EQNTYPE}
                )

            # small rail matrices
            ADD_TEST(dense_pcare_nm_rail_small_${EQNTYPENAME}-plus${PLUS}-line${LINE}-withE
                ${EXECUTABLE_OUTPUT_PATH}/check_dense_nm
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.A
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.E
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.B
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.C
                ${PLUS}
                ${LINE}
                ${EQNTYPE}
                0
                )

            ADD_TEST(dense_nm_rail_small_${EQNTYPENAME}-plus${PLUS}-line${LINE}-withoutE

                ${EXECUTABLE_OUTPUT_PATH}/check_dense_nm
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.A
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.B
                ${CMAKE_SOURCE_DIR}/tests/data/Rail/rail_371_c60.C
                ${PLUS}
                ${LINE}
                ${EQNTYPE}
                0
                )


            FOREACH(INSTANCE iss build heat-cont)
                ADD_TEST(dense_nm_${INSTANCE}_${EQNTYPENAME}-plus${PLUS}-line${LINE}

                    ${EXECUTABLE_OUTPUT_PATH}/check_dense_nm
                    ${CMAKE_SOURCE_DIR}/tests/data/${INSTANCE}/A.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/${INSTANCE}/E.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/${INSTANCE}/B.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/${INSTANCE}/C.mtx
                    ${PLUS}
                    ${LINE}
                    ${EQNTYPE}
                    0
                    )
            ENDFOREACH()


            IF(LARGETESTS STREQUAL ON OR HUGETESTS STREQUAL ON)

                ADD_TEST(dense_nm_rail_${EQNTYPENAME}-plus${PLUS}-line${LINE}-withE
                    ${EXECUTABLE_OUTPUT_PATH}/check_dense_nm
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/A.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/E.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/B.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/C.mtx
                    ${PLUS}
                    ${LINE}
                    ${EQNTYPE}
                    0
                    )

                ADD_TEST(dense_nm_rail_${EQNTYPENAME}-plus${PLUS}-line${LINE}-withoutE

                    ${EXECUTABLE_OUTPUT_PATH}/check_dense_nm
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/A.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/B.mtx
                    ${CMAKE_SOURCE_DIR}/tests/data/Rail/C.mtx
                    ${PLUS}
                    ${LINE}
                    ${EQNTYPE}
                    0
                    )
            ENDIF()
        ENDFOREACH()
    ENDFOREACH()
ENDFOREACH()
