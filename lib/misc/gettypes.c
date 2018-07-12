//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.
//
// Copyright (C) Peter Benner, Martin Koehler, Jens Saak and others
//               2009-2018
//

/**
 * @file lib/misc/gettypes.c
 * @brief Get human readable string from type identifiers.
 * @author @koehlerm
 * @author @mbehr
 *
 * This file contain function to retrieve human readable strings from type definitions
 * and type enumerations.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"



#define MESS_ENUM_STR_1(RETALT,ARG,ENUM)            if( ARG == ENUM ){ return #ENUM; } return RETALT;
#define MESS_ENUM_STR_2(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_1(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_3(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_2(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_4(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_3(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_5(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_4(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_6(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_5(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_7(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_6(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_8(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_7(RETALT,ARG,__VA_ARGS__) }
#define MESS_ENUM_STR_9(RETALT,ARG,ENUM,...)        if( ARG == ENUM ){ return #ENUM; } else { MESS_ENUM_STR_8(RETALT,ARG,__VA_ARGS__) }



/**
 * @brief Return a human readable name of a data type value.
 * @param[in] data_type input data type value
 * @return constant pointer to a human readable string containing the name of the data type
 *
 * The @ref mess_datatype_t_str function translates a data type constant to the
 * corresponding string and returns a constant pointer to this string.
 *
 * @see mess_datatype_t
 *
 */
const char *mess_datatype_t_str ( const mess_datatype_t data_type ){
    MESS_ENUM_STR_2("Unknown mess_datatype_t", data_type, MESS_REAL, MESS_COMPLEX);
}


/**
 * @brief Return a human readable name of a storage type value.
 * @param[in] store_type input storage type value
 * @return constant pointer to a human readable string containing the name of the storage type
 *
 * The @ref mess_storage_t_str function translates a storage type constant to the
 * corresponding string and returns a constant pointer to this string.
 *
 * @see mess_storage_t
 *
 */
const char *  mess_storage_t_str ( const mess_storage_t  store_type ){
     MESS_ENUM_STR_4("Unknown mess_storage_t", store_type, MESS_CSR, MESS_CSC, MESS_DENSE, MESS_COORD);
}


/**
 * @brief Return a human readable name of a symmetry type value.
 * @param[in] sym input symmetry type value
 * @return constant pointer to a human readable string containing the name of the symmetry type
 *
 * The @ref mess_symmetry_t_str function translates a symmetry type constant to the
 * corresponding string and returns a constant pointer to this string.
 *
 * @see mess_symmetry_t
 *
 */
const char *  mess_symmetry_t_str ( const mess_symmetry_t sym  ){
    MESS_ENUM_STR_4("Unknown mess_symmetry_t", sym, MESS_GENERAL, MESS_SYMMETRIC, MESS_SKEWSYMMETRIC, MESS_HERMITIAN );
}


/**
 * @brief Return the human readable name of an operation type value.
 * @param[in] op input  operation type value
 * @return constant pointer to a human readable string containing the name of the operation type
 *
 * The @ref mess_operation_t_str function translates an operation type constant to the
 * corresponding string and returns a constant pointer to this string.
 *
 * @see mess_operation_t
 *
 */
const char *  mess_operation_t_str ( const mess_operation_t op  ){
    MESS_ENUM_STR_3("Unknown mess_operation_t", op, MESS_OP_NONE, MESS_OP_TRANSPOSE, MESS_OP_HERMITIAN);
}


/**
 * @brief Convert the @ref mess_equation_t to a human readable string.
 * @param[in] eqn_type equation type
 * @return human readable string of @p eqn_type
 *
 * The @ref mess_equation_t_str function converts a @ref mess_equation_t to a human readable string.
 *
 */
const char* mess_equation_t_str(const mess_equation_t eqn_type){
    MESS_ENUM_STR_5("Unknown mess_equation_t", eqn_type, MESS_EQN_NONE, MESS_EQN_LYAP, MESS_EQN_GLYAP, MESS_EQN_RICCATI, MESS_EQN_GRICCATI);
}


/**
 * @brief Convert the @ref mess_parameter_t to a human readable string.
 * @param[in] para_type parameter type
 * @return human readable string of @p para_type
 *
 * The @ref mess_parameter_t_str function converts a @ref mess_parameter_t to a human readable string.
 *
 */
const char* mess_parameter_t_str(const mess_parameter_t para_type){
    MESS_ENUM_STR_5("Unknown mess_parameter_t", para_type, MESS_LRCFADI_PARA_MINMAX, MESS_LRCFADI_PARA_MINMAX_REAL, MESS_LRCFADI_PARA_WACHSPRESS,
                                                           MESS_LRCFADI_PARA_ADAPTIVE_V, MESS_LRCFADI_PARA_ADAPTIVE_Z);
}


/**
 * @brief Convert the @ref mess_memusage_t to a human readable string.
 * @param[in] mem_type memory type
 * @return human readable string of @p mem_type
 *
 * The @ref mess_memusage_t_str function converts a @ref mess_memusage_t to a human readable string.
 *
 */
const char* mess_memusage_t_str(const mess_memusage_t mem_type){
    MESS_ENUM_STR_3("Unknown mess_memusage_t", mem_type, MESS_MEMORY_LOW, MESS_MEMORY_MID, MESS_MEMORY_HIGH);
}


/**
 * @brief Convert the @ref mess_norm_t to a human readable string.
 * @param[in] nrm_type memory type
 * @return human readable string of @p nrm_type
 *
 * The @ref mess_norm_t_str function converts a @ref mess_norm_t to a human readable string.
 *
 */
const char* mess_norm_t_str(const mess_norm_t nrm_type){
    MESS_ENUM_STR_4("Unknown mess_norm_t", nrm_type, MESS_2_NORM, MESS_FROBENIUS_NORM, MESS_1_NORM, MESS_INF_NORM);
}


/**
 * @brief Convert the @ref mess_sign_scale_t to a human readable string.
 * @param[in] sc_type type of sign scaling method
 * @return human readable string of @p sc_type
 *
 * The @ref mess_sign_scale_t_str function converts a @ref mess_sign_scale_t to a human readable string.
 *
 */
const char* mess_sign_scale_t_str(const mess_sign_scale_t sc_type){
    MESS_ENUM_STR_3("Unknown mess_sign_scale_t", sc_type,   MESS_SIGN_SCALE_NONE,
                                                            MESS_SIGN_SCALE_FRO,
                                                            MESS_SIGN_SCALE_DET);
}


/**
 * @brief Convert the @ref mess_direct_lupackage_t to a human readable string.
 * @param[in] lu_type type of direct solver
 * @return human readable string of @p lu_type
 *
 * The @ref mess_direct_lupackage_t_str function converts a @ref mess_direct_lupackage_t to a human readable string.
 *
 */
const char* mess_direct_lupackage_t_str(const mess_direct_lupackage_t lu_type){
    MESS_ENUM_STR_8("Unknown mess_direct_lupackage_t", lu_type, MESS_DIRECT_DEFAULT_LU,
                                                                MESS_DIRECT_SPARSE_LU,
                                                                MESS_DIRECT_LAPACK_LU,
                                                                MESS_DIRECT_UMFPACK_LU,
                                                                MESS_DIRECT_SUPERLU_LU,
                                                                MESS_DIRECT_CSPARSE_LU,
                                                                MESS_DIRECT_BANDED_LU,
                                                                MESS_DIRECT_MKLPARDISO_LU);
}


/**
 * @brief Convert the @ref mess_direct_cholpackage_t to a human readable string.
 * @param[in] chol_type type of direct solver
 * @return human readable string of @p chol_type
 *
 * The @ref mess_direct_cholpackage_t_str function converts a @ref mess_direct_cholpackage_t to a human readable string.
 *
 */
const char* mess_direct_cholpackage_t_str(const mess_direct_cholpackage_t chol_type){
    MESS_ENUM_STR_4("Unknown mess_direct_cholpackage_t", chol_type, MESS_DIRECT_DEFAULT_CHOLESKY,
                                                                    MESS_DIRECT_LAPACK_CHOLESKY,
                                                                    MESS_DIRECT_CSPARSE_CHOLESKY,
                                                                    MESS_DIRECT_CHOLMOD_CHOLESKY);
}


/**
 * @brief Convert the @ref mess_multidirect_t to a human readable string.
 * @param[in] ml_type type of multidirect solver
 * @return human readable string of @p ml_type
 *
 * The @ref mess_multidirect_t_str function converts a @ref mess_multidirect_t to a human readable string.
 *
 */
const char* mess_multidirect_t_str(const mess_multidirect_t ml_type){
    MESS_ENUM_STR_2("Unknown mess_multidirect_t", ml_type,  MESS_MULTIDIRECT_SPARSE_LU, MESS_MULTIDIRECT_UMFPACK_LU);
}


/**
 * @brief Convert the @ref mess_residual_t to a human readable string.
 * @param[in] res_type type of residual method
 * @return human readable string of @p res_type
 *
 * The @ref mess_residual_t_str function converts a @ref mess_residual_t to a human readable string.
 *
 */
const char* mess_residual_t_str(const mess_residual_t res_type){
    MESS_ENUM_STR_2("Unknown mess_residual_t", res_type,  MESS_RESIDUAL_INDEFINITE, MESS_RESIDUAL_SPECTRAL);
}

/**
 * @brief Convert the @ref mess_print_format_t to a human readable string.
 * @param[in] print_type type of print methods
 * @return human readable string of @p print_type
 *
 * The @ref mess_print_format_t_str function converts a @ref mess_print_format_t to a human readable string.
 *
 */
const char* mess_print_format_t_str(const mess_print_format_t print_type){
    MESS_ENUM_STR_2("Unknown mess_print_format_t", print_type,  MESS_PRINT_FORMAT_SHORT, MESS_PRINT_FORMAT_LONG);
}


