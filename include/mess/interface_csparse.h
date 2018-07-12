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
 * @file include/mess/interface_csparse.h
 * @brief Interface to @csparse.
 * @ref mess.h and cs.h from @csparse should be included in the program already
 */

#ifndef INTERFACE_CSPARSE_H_
#define INTERFACE_CSPARSE_H_

#include "mess/error_macro.h"
#include "mess/config.h"
#include <string.h>

/** @addtogroup interfaces_csparse
 * @{
 */
/**
 * @brief Copy a real matrix from @csparse to @mess (double-integer variant).
 * @param[in] input input @csparse real matrix with int index values
 * @param[out] output output @ref  mess_matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_matrix_from_csparse_di function copies a real matrix from @csparse to @mess.
 *
 */
__attribute__((unused)) static int mess_matrix_from_csparse_di(const cs_di * input, mess_matrix output)
{
    MSG_FNAME(__func__);
    mess_int_t i = 0;
    int ret = 0;


    mess_check_nullpointer(input);
    mess_check_nullpointer(output);

    ret = mess_matrix_alloc(output, input->m,input->n, input->p[input->n], MESS_CSC,MESS_REAL);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);

    memcpy ( output->values, input->x, output->nnz*sizeof(double));
    for ( i = 0; i < output->nnz; i++) {
        output->rowptr[i] = input->i[i];
    }
    for ( i = 0; i <= input->n;  i++){
        output->colptr[i] = input->p[i];
    }

    return 0;
}

/**
 * @brief Default wrapper for mess_matrix_from_csparse_di.
 * @param[in] input @csparse input matrix
 * @param[out] output mess output matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_matrix_from_csparse is a wrapper to @ref mess_matrix_from_csparse_di.
 */
__attribute__((unused)) static int mess_matrix_from_csparse(const cs * input, mess_matrix output) {
    return mess_matrix_from_csparse_di(input, output);
}


/**
 * @brief Copy a real matrix from @csparse to @mess (double-long variant).
 * @param[in] input input @csparse real matrix with long index values
 * @param[out] output output @ref mess_matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_matrix_from_csparse_dl function copies a real matrix from @csparse to mess.
 *
 */
__attribute__((unused)) static int mess_matrix_from_csparse_dl(const cs_dl * input, mess_matrix output)
{
    MSG_FNAME(__func__);
    mess_int_t i = 0;
    int ret = 0;


    mess_check_nullpointer(input);
    mess_check_nullpointer(output);

    ret = mess_matrix_alloc(output, input->m,input->n, input->p[input->n], MESS_CSC,MESS_REAL);
    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_alloc);

    memcpy ( output->values, input->x, output->nnz*sizeof(double));
    for ( i = 0; i < output->nnz; i++) {
        output->rowptr[i] = input->i[i];
    }
    for ( i = 0; i <= input->n;  i++){
        output->colptr[i] = input->p[i];
    }

    return 0;
}

/**
 * @brief Copy a complex matrix from @csparse to @mess (complex-int variant).
 * @param[in] input input @csparse complex matrix with int index values
 * @param[out] output output @ref mess_matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_matrix_from_csparse_ci function copies a complex matrix from @csparse to @mess.
 *
 */
__attribute__((unused)) static int mess_matrix_from_csparse_ci(const cs_ci * input, mess_matrix output)
{
    MSG_FNAME(__func__);
    mess_int_t i=0;
    int ret=0;

    mess_check_nullpointer(input);
    mess_check_nullpointer(output);

    ret = mess_matrix_alloc(output, input->m,input->n, input->p[input->n], MESS_CSC,MESS_COMPLEX);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    memcpy ( output->values_cpx, input->x, output->nnz*sizeof(mess_double_cpx_t));
    for ( i = 0; i < output->nnz; i++) {
        output->rowptr[i] = input->i[i];
    }
    for ( i = 0; i <= input->n;  i++){
        output->colptr[i] = input->p[i];
    }
    return 0;
}

/**
 * @brief Copy a complex matrix from @csparse to @mess (complex-long variant).
 * @param[in] input input @csparse complex matrix with long index values
 * @param[out] output output @ref mess_matrix
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_matrix_from_csparse_cl function copies a complex matrix from @csparse to @mess.
 *
 */
__attribute__((unused)) static int mess_matrix_from_csparse_cl(const cs_cl * input, mess_matrix output)
{
    MSG_FNAME(__func__);
    mess_int_t i=0;
    int ret=0;

    mess_check_nullpointer(input);
    mess_check_nullpointer(output);

    ret = mess_matrix_alloc(output, input->m,input->n, input->p[input->n], MESS_CSC,MESS_COMPLEX);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    memcpy ( output->values_cpx, input->x, output->nnz*sizeof(mess_double_cpx_t));
    for ( i = 0; i < output->nnz; i++) {
        output->rowptr[i] = input->i[i];
    }
    for ( i = 0; i <= input->n;  i++){
        output->colptr[i] = input->p[i];
    }
    return 0;
}



/**
 * @internal
 * @brief Copy a real matrix from @mess to @csparse (double-int variant).
 * @param[in] input     input @ref mess_matrix
 * @param[out] o    output @csparse matrix with int index values
 *
 * The @ref mess_matrix_to_csparse_di function copies a real matrix from @mess to @csparse.
 *
 * @attention Interal use only.
 */
static int __attribute__((unused)) mess_matrix_to_csparse_di(mess_matrix input, cs_di **o) {
    MSG_FNAME(__func__);
    int conv = -1;
    mess_int_t i;
    mess_matrix wI;
    cs_di *output=NULL;

    mess_check_nullpointer(input);
    mess_check_nullpointer(o);
    mess_check_real(input);

    MESS_MATRIX_CHECKFORMAT(input, wI, conv, MESS_CSC);
    output = cs_di_spalloc(wI->rows, wI->cols, wI->nnz, 1, 0);
    if ( output == NULL ) {
        MSG_ERROR("can not allocated cs_dl matrix for output");
        return MESS_ERROR_MEMORY;
    }


    for ( i = 0; i <=wI->cols; i++){
        output->p[i] = wI -> colptr[i];
    }
    for ( i = 0; i < wI->nnz; i++){
        output->i[i] = wI->rowptr[i];
    }
    memcpy(output->x, wI->values, sizeof(double)*(wI->nnz));

    if (conv == 0) {
        mess_matrix_clear(&wI);
    }
    *o= output;
    return 0;
}

/**
 * @internal
 * @brief Wrapper around mess_matrix_to_csparse_di.
 * @param[in] input     input @ref mess_matrix
 * @param[out] o    output @csparse matrix
 *
 * The @ref mess_matrix_to_csparse function is a wrapper to @ref mess_matrix_to_csparse_di.
 * @attention Interal use only.
 **/
static int __attribute__((unused)) mess_matrix_to_csparse(mess_matrix input, cs **o) {
    return mess_matrix_to_csparse_di(input,(cs_di**) o);
}

/**
 * @internal
 * @brief Copy a real matrix from @mess to @csparse (double-long variant).
 * @param[in] input     input @ref mess_matrix
 * @param[out] o    output @csparse matrix with long index values
 *
 * The @ref mess_matrix_to_csparse_dl function copies a real matrix from @mess to @csparse.
 * @attention Interal use only.
 *
 */
static int __attribute__((unused)) mess_matrix_to_csparse_dl(mess_matrix input, cs_dl **o) {
    MSG_FNAME(__func__);
    int conv = -1;
    mess_int_t i;
    mess_matrix wI;
    cs_dl *output=NULL;

    mess_check_nullpointer(input);
    mess_check_nullpointer(o);
    mess_check_real(input);

    MESS_MATRIX_CHECKFORMAT(input, wI, conv, MESS_CSC);
    output = cs_dl_spalloc(wI->rows, wI->cols, wI->nnz, 1, 0);
    if ( output == NULL ) {
        MSG_ERROR("can not allocated cs_dl matrix for output");
        return MESS_ERROR_MEMORY;
    }


    for ( i = 0; i <=wI->cols; i++){
        output->p[i] = wI -> colptr[i];
    }
    for ( i = 0; i < wI->nnz; i++){
        output->i[i] = wI->rowptr[i];
    }
    memcpy(output->x, wI->values, sizeof(double)*(wI->nnz));

    if (conv == 0) {
        mess_matrix_clear(&wI);
    }
    *o= output;
    return 0;
}

/**
 * @internal
 * @brief Copy a complex matrix from @mess to @csparse (complex-int variant).
 * @param[in] input     input @ref mess_matrix
 * @param[out] o    output CXSparse matrix with int index values
 *
 * The @ref mess_matrix_to_csparse_ci function copies a complex matrix from @mess to @csparse.
 *
 * @attention Interal use only.
 */
static int __attribute__((unused)) mess_matrix_to_csparse_ci(mess_matrix input, cs_ci **o) {
    MSG_FNAME(__func__);
    int conv = -1;
    mess_int_t i;
    mess_matrix wI;
    cs_ci *output=NULL;

    mess_check_nullpointer(input);
    mess_check_nullpointer(o);
    mess_check_complex(input);

    // output = cs_cl_spalloc(input->rows, input->cols, input->nnz, 1, 0);

    MESS_MATRIX_CHECKFORMAT(input, wI, conv, MESS_CSC);
    output = cs_ci_spalloc(wI->rows, wI->cols, wI->nnz, 1, 0);
    if ( output == NULL ) {
        MSG_ERROR("can not allocated cs_dl matrix for output");
        return MESS_ERROR_MEMORY;
    }


    for ( i = 0; i <=wI->cols; i++){
        output->p[i] = wI -> colptr[i];
    }
    for ( i = 0; i < wI->nnz; i++){
        output->i[i] = wI->rowptr[i];
    }
    memcpy(output->x, wI->values_cpx, sizeof(mess_double_cpx_t)*(wI->nnz));

    if (conv == 0) {
        mess_matrix_clear(&wI);
    }
    *o= output;
    return 0;
}

/**
 * @internal
 * @brief Copy a complex matrix from @mess to @csparse (complex-long variant).
 * @param[in] input     input @ref mess_matrix
 * @param[out] o    output @csparse matrix with long index values
 *
 * The @ref mess_matrix_to_csparse_cl function copies a complex matrix from @mess to @csparse.
 *
 * @attention Interal use only.
 */
static int __attribute__((unused)) mess_matrix_to_csparse_cl(mess_matrix input, cs_cl **o) {
    MSG_FNAME(__func__);
    int conv = -1;
    mess_int_t i;
    mess_matrix wI;
    cs_cl *output=NULL;

    mess_check_nullpointer(input);
    mess_check_nullpointer(o);
    mess_check_complex(input);

    // output = cs_cl_spalloc(input->rows, input->cols, input->nnz, 1, 0);

    MESS_MATRIX_CHECKFORMAT(input, wI, conv, MESS_CSC);
    output = cs_cl_spalloc(wI->rows, wI->cols, wI->nnz, 1, 0);
    if ( output == NULL ) {
        MSG_ERROR("can not allocated cs_dl matrix for output");
        return MESS_ERROR_MEMORY;
    }


    for ( i = 0; i <=wI->cols; i++){
        output->p[i] = wI -> colptr[i];
    }
    for ( i = 0; i < wI->nnz; i++){
        output->i[i] = wI->rowptr[i];
    }
    memcpy(output->x, wI->values_cpx, sizeof(mess_double_cpx_t)*(wI->nnz));

    if (conv == 0) {
        mess_matrix_clear(&wI);
    }
    *o= output;
    return 0;
}

/** @}  */
#endif /* INTERFACE_CSPARSE_H_ */
