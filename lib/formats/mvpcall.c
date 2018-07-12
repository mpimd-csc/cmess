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
 * @file lib/formats/mvpcall.c
 * @brief Handling generic matrix-vector products.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/*-----------------------------------------------------------------------------
 *  Normal MVP
 *-----------------------------------------------------------------------------*/
static int mvp_normal(void *data, mess_operation_t op, mess_vector x, mess_vector y) {
    return mess_matrix_mvp(op, (mess_matrix) data, x, y);
}


/*-----------------------------------------------------------------------------
 *  Tranposed MVP
 *-----------------------------------------------------------------------------*/
static int mvp_transpose(void *data, mess_operation_t op, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    int ret = 0;
    switch(op){
        case MESS_OP_NONE:
            return mess_matrix_mvp(MESS_OP_TRANSPOSE, (mess_matrix) data, x, y);
        case MESS_OP_TRANSPOSE:
            return mess_matrix_mvp(MESS_OP_NONE, (mess_matrix) data, x, y);
        case MESS_OP_HERMITIAN:
            ret = mess_matrix_mvp(MESS_OP_NONE, (mess_matrix) data, x, y);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
            ret = mess_vector_conj(y);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_conj);
            break;
    }
    return 0;

}


/*-----------------------------------------------------------------------------
 *  Hermitian MVP
 *-----------------------------------------------------------------------------*/
static int mvp_hermitian(void *data, mess_operation_t op, mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    int ret = 0;
    switch(op){
        case MESS_OP_NONE:
            return mess_matrix_mvp(MESS_OP_HERMITIAN, (mess_matrix) data, x, y);
        case MESS_OP_TRANSPOSE:
            ret = mess_matrix_mvp(MESS_OP_NONE, (mess_matrix) data, x, y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_mvp );
            ret = mess_vector_conj(y);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_conj );
            break;
        case MESS_OP_HERMITIAN:
            return mess_matrix_mvp(MESS_OP_NONE, (mess_matrix) data, x, y);
            break;
    }
    return 0;
}


/**
 * @brief Generate a mess_mvpcall object for a simple matrix.
 * @param[out]  mvpcall   generated mess_mvpcall object based on the given matrix
 * @param[in]   op   input operation applied to \f$ A \f$
 * @param[in]   A    input matrix of the matrix-vector product
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_mvpcall_matrix function generates a @ref mess_mvpcall object which represents
 * \f[ mvpcall = op(A)  \f]
 * where \f$ A \f$ is a mess_matrix and \f$ op \f$ can be
 * <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 */
int mess_mvpcall_matrix(mess_mvpcall * mvpcall, mess_operation_t op, mess_matrix A){
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_real_or_complex(A);
    //mess_check_square(A);

    /*-----------------------------------------------------------------------------
     *  Setup the structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(*mvpcall, struct mess_mvpcall_st *, sizeof(struct mess_mvpcall_st));

    (*mvpcall)->dim = A->rows;
    (*mvpcall)->data_type = A->data_type;
    (*mvpcall)->data = (void *) A;

    switch(op){
        case MESS_OP_NONE:
            (*mvpcall)->mvp = mvp_normal;
            break;
        case MESS_OP_TRANSPOSE:
            (*mvpcall)->mvp = mvp_transpose;
            break;
        case MESS_OP_HERMITIAN:
            (*mvpcall)->mvp = mvp_hermitian;
            break;
    }
    return 0;
}

/**
 * @brief Generate a mess_mvpcall object for a complex operator.
 * @param[out]  mvpcall generated mess_mvpcall object based on the given function handle
 * @param[in]   dim      input dimension of the represented operator
 * @param[in]   data_type input data type of the represented operator
 * @param[in]   mvp      input function pointer to the matrix-vector product function
 * @param[in]   data     input pointer to additional data for the matrix-vector product
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_mvpcall_operator function generates a @ref mess_mvpcall object which represents
 * \f[ y = A x \f]
 * where \f$ A \f$ is a complex operator represented by a function pointer \f$ mvp \f$.
 *
 */
int mess_mvpcall_operator(mess_mvpcall *mvpcall, mess_int_t dim, mess_datatype_t data_type, int (*mvp)(void *data, mess_operation_t op, mess_vector x, mess_vector y), void *data){
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_positive(dim);
    mess_check_nullpointer(mvp);
    if ( data_type != MESS_REAL && data_type != MESS_COMPLEX) {
        MSG_ERROR("The data type must be real or complex.\n");
        return MESS_ERROR_DATATYPE;
    }

    /*-----------------------------------------------------------------------------
     *  Setup the structure
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(*mvpcall, mess_mvpcall , sizeof(struct mess_mvpcall_st));
    (*mvpcall)->dim = dim;
    (*mvpcall)->mvp = mvp;
    (*mvpcall)->data_type = data_type;
    (*mvpcall)->data = data;
    return 0;
}


/**
 * @brief Apply a generic matrix-vector product.
 * @param[in] mvpcall  input mess_mvpcall object defining the matrix-vector product
 * @param[in] op  input operation applied to the matrix
 * @param[in] x   input right hand side vector
 * @param[out] y vector
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_mvpcall_apply function computes
 * \f[ y = op(A) x \f]
 * where \f$ A \f$ is represented by a @ref mess_mvpcall object \f$ mvpcall \f$.
 */
int mess_mvpcall_apply(mess_mvpcall mvpcall, mess_operation_t op, mess_vector x, mess_vector y){
    return mvpcall->mvp(mvpcall->data, op, x, y);
}


/**
 * @brief Clean a mess_mvpcall object .
 * @param[in] mvpcall    input mess_mvpcall object to clear
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_mvpcall_clear function clears a @ref  mess_mvpcall object and sets it
 * to @c NULL.
 */
int mess_mvpcall_clear(mess_mvpcall *mvpcall){
    if ( mvpcall == NULL) return 0;
    if ( *mvpcall == NULL) return 0;
    mess_free(*mvpcall);
    *mvpcall = NULL;
    return 0;
}

