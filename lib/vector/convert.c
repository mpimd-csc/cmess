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
 * @file lib/vector/convert.c
 * @brief Convert functions of @ref mess_vector structures.
 * @author @koehlerm
 * @author @mbehr
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include "blas_defs.h"
#ifdef _OPENMP_H
#include <omp.h>
#endif

/**
 * On MacOS X 10.6 the zdotc implementation
 * is defect, this is the work around
 */
#ifdef MESS_USE_APPLE_BLAS
#include <Accelerate/Accelerate.h>
#endif





/**
 * @brief Create a @ref mess_vector from a @c Fortran array.
 * @param[out] v        output vector to copy the data to
 * @param[in] dim       dimension of the vector
 * @param[in] vals      input array with real values
 * @param[in] vals_cpx  input array with complex values
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_vector_from_farray function copies a @p vals or @p val_cpx array to
 * a @ref mess_vector.
 * Depending whether @p vals  or @p vals_cpx  is given a real @ref mess_vector or a complex @ref mess_vector is created.
 *
 * @sa mess_matrix_dense_from_carray
 * @sa mess_matrix_dense_from_farray
 *
 */
int mess_vector_from_farray(mess_vector v, mess_int_t dim, double *vals, mess_double_cpx_t *vals_cpx ){
    MSG_FNAME(__func__);
    mess_datatype_t dt;

    /*-----------------------------------------------------------------------------
     * Check Inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(v);
    mess_check_positive(dim);

    if ( (!vals && !vals_cpx) || (vals && vals_cpx ) ) {
        MSG_ERROR("Either vals or vals_cpx must be given.");
        return MESS_ERROR_ARGUMENTS;
    }

    dt = vals? MESS_REAL:MESS_COMPLEX;

    /*-----------------------------------------------------------------------------
     *  reset and set values
     *-----------------------------------------------------------------------------*/
    v->data_type = dt;
    v->dim = dim;
    mess_free(v->values);
    mess_free(v->values_cpx);
    if(dt == MESS_REAL){
        mess_try_alloc(v->values, double*, sizeof(double)*dim);
        memcpy(v->values, vals, sizeof(double)*dim);
    }else{
        mess_try_alloc(v->values_cpx, mess_double_cpx_t*, sizeof(mess_double_cpx_t)*dim);
        memcpy(v->values_cpx, vals_cpx, sizeof(mess_double_cpx_t)*dim);
    }
    return 0;
}       /* -----  end of function mess_vector_from_farray  ----- */




/**
 * @brief Create a @ref mess_vector from two  @c Fortran array.
 * @param[out] v        output vector to copy the data to
 * @param[in] dim       dimension of the vector
 * @param[in] vals_re   input array with real part values, @c NULL if not wanted
 * @param[in] vals_im   input array with imaginary part values, @c NULL if not wanted
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_vector_from_lapack function copies a @p vals_re and @p val_im array to
 * a @ref mess_vector. If @p vals_im points to @c NULL or all entries are zero @p v will be a real vector.
 * The function is mostly used when results from @lapack are splited in real an imaginary part and should be copied to
 * a @ref mess_vector.
 *
 * If  @p vals_re and @p  vals_im point to @c NULL an to @c NULL an error is returned.
 *
 */
int mess_vector_from_lapack(mess_vector v, mess_int_t dim, double *vals_re, double *vals_im){
    MSG_FNAME(__func__);
    mess_datatype_t dt = MESS_REAL;
    mess_int_t i = 0;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     * Check Inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(v);
    mess_check_positive(dim);

    if ( !vals_re ) {
        MSG_ERROR("Either vals_re or vals_im must be given.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  check if complex data are needed
     *-----------------------------------------------------------------------------*/
    if (vals_im){
        for(i=0;i<dim;++i){
            if(!vals_im[i]){
                dt = MESS_COMPLEX;
                break;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     * check if a realloc is necessary
     *-----------------------------------------------------------------------------*/
    if( v->dim != dim || dt != v->data_type){
        ret = mess_vector_reset(v);                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_reset);
        ret = mess_vector_alloc(v, dim, dt);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    }

    /*-----------------------------------------------------------------------------
     *  copying
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(v)){
        for(i=0;i<dim;++i){
            v->values[i] = vals_re[i];
        }
    }else{
        for(i=0;i<dim;++i){
            v->values_cpx[i]=vals_re[i] + I*vals_im[i];
        }
    }

    return 0;
}       /* -----  end of function mess_vector_from_lapack ----- */



/**
 *
 * @brief Convert complex vector to real vector if imaginary part of every entry is small.
 * @param[in,out] v     input/output vector  \f$v\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_convert_if_real converts a complex vector to
 * a real one if the imaginary part of all entries are small.
 * More precise if \f$ \frac{|v_i - \operatorname{Re}(v_i)|}{|\operatorname{Re}(v_i)|}\leq eps\ \forall i \f$ then
 * we perform \f$ v \leftarrow \operatorname{Re}(v) \f$.
 *
 */
int mess_vector_convert_if_real(mess_vector v) {
    MSG_FNAME(__func__);
    int ret = 0;
    double eps = mess_eps();
    mess_int_t i =0;
    mess_double_cpx_t entry;
    double realpart;


    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);
    if(MESS_IS_REAL(v)){
        return 0;
    }


    /*-----------------------------------------------------------------------------
     *  check for complex data
     *-----------------------------------------------------------------------------*/
    for(i=0;i<v->dim;++i){
        entry = v->values_cpx[i];
        realpart = creal(entry);
        if(cabs(entry-realpart)>eps*fabs(realpart))
            break;
    }

    if(i>=v->dim){
        ret = mess_vector_toreal_nowarn(v);    FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_toreal_nowarn);
    }

    return 0;
} /* -----  end of function mess_vector_convert_if_real  ----- */


/**
 * @brief Convert a vector to a matrix.
 * @param[in] v input vector
 * @param[out] mat output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_tomatrix converts a @ref mess_vector structure to a @ref mess_matrix structure with  one column.
 */
int mess_vector_tomatrix(mess_vector v, mess_matrix mat) {
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(v);
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(v);
    mess_check_real_or_complex(mat);

    mess_matrix_reset(mat);
    if ((ret = mess_matrix_alloc(mat, v->dim, 1, v->dim, MESS_DENSE, v->data_type))) {
        MSG_ERROR("mess_matrix_alloc returned an error: %d\n", ret);
        return ret;
    }

    if (MESS_IS_COMPLEX(v)) {
        memcpy(mat->values_cpx, v->values_cpx, v->dim*sizeof(mess_double_cpx_t));
    } else {
        memcpy(mat->values, v->values, v->dim * sizeof(double));
    }
    return 0;
}


/**
 * @brief Convert a single column matrix to a vector.
 * @param[in] mat   input single column input matrix
 * @param[out] v    output vector with the elements from the matrix column
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_vector_frommatrix function converts a single column matrix into a vector.
 * @sa mess_vector_tomatrix
 */
int mess_vector_frommatrix(mess_matrix mat, mess_vector v)
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(mat);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(mat);
    mess_check_real_or_complex(v);

    if ( mat -> cols != 1) {
        MSG_ERROR("Matrix must only have one column.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( MESS_IS_REAL(mat)) {
        ret = mess_vector_toreal_nowarn(v);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        ret = mess_vector_resize(v, mat->rows);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        memcpy(v->values, mat->values, sizeof(double) * v->dim);
    } else if ( MESS_IS_COMPLEX(mat)) {
        ret = mess_vector_tocomplex(v);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        ret = mess_vector_resize(v, mat->rows);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        memcpy(v->values_cpx, mat->values_cpx, sizeof(mess_double_cpx_t) * v->dim);
    }
    return 0;
}




/**
 * @brief Convert a vector to complex one.
 * @param[in,out] v input/output vector to convert
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_tocomplex function converts an arbitrary vector to a complex
 * vector.
 *
 * @see mess_vector_toreal
 * @see mess_vector_toreal_nowarn
 * @see mess_vector_totype
 *
 */
int mess_vector_tocomplex(mess_vector v) {
    MSG_FNAME(__func__);
    mess_int_t i;

    mess_check_nullpointer(v);

    if (MESS_IS_COMPLEX(v))
        return 0;
    mess_try_alloc(v->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t) * v->dim);
    for (i = 0; i < v->dim; i++)
        v->values_cpx[i] = v->values[i];
    mess_free(v->values);
    v->values = NULL;
    v->data_type = MESS_COMPLEX;
    return 0;
}

/**
 * @brief Convert a vector to real one.
 * @param[in,out] v input/output vector to convert
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_toreal function converts an arbitrary vector to a real
 * vector. Imaginary parts are ignored.
 *
 * @see mess_vector_tocomplex
 * @see mess_vector_toreal_nowarn
 * @see mess_vector_totype
 *
 */
int mess_vector_toreal(mess_vector v) {
    MSG_FNAME(__func__);
    mess_int_t i;

    mess_check_nullpointer(v);

    if (MESS_IS_COMPLEX(v)) {
        MSG_WARN("complex vector is converted to real, some information can be lost.\n");
        // abort();
    }
    if (MESS_IS_REAL(v))
        return 0;
    mess_try_alloc(v->values, double *, sizeof(double) * v->dim);
    for (i = 0; i < v->dim; i++)
        v->values[i] = creal(v->values_cpx[i]);
    mess_free(v->values_cpx);
    v->values_cpx = NULL;
    v->data_type = MESS_REAL;
    return 0;
}

/**
 * @brief Convert a vector to real values (without any warnings).
 * @param[in,out] v input/output vector to convert
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_toreal_nowarn function converts an arbitrary vector to a real
 * vector. Imaginary parts are ignored and no warnings are displayed.
 *
 * @see mess_vector_tocomplex
 * @see mess_vector_toreal
 * @see mess_vector_totype
 *
 */
int mess_vector_toreal_nowarn(mess_vector v) {
    MSG_FNAME(__func__);
    mess_int_t i;

    mess_check_nullpointer(v);

    if (MESS_IS_REAL(v))
        return 0;
    mess_try_alloc(v->values, double *, sizeof(double) * v->dim);
    for (i = 0; i < v->dim; i++)
        v->values[i] = creal(v->values_cpx[i]);
    mess_free(v->values_cpx);
    v->values_cpx = NULL;
    v->data_type = MESS_REAL;
    return 0;
}

/**
 * @brief Convert a vector to the given datatype.
 * @param[in,out] v input/output given vector
 * @param[in] dt    input   desired datatype
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_totype function converts a vector to a given datatype.
 *
 * @see mess_vector_tocomplex
 * @see mess_vector_toreal
 * @see mess_vector_toreal_nowarn
 *
 */
int mess_vector_totype(mess_vector v, mess_datatype_t dt) {
    MSG_FNAME(__func__);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);
    if (dt == v->data_type)
        return 0;
    switch (dt) {
        case MESS_REAL:
            return mess_vector_toreal_nowarn(v);
            break;
        case MESS_COMPLEX:
            return mess_vector_tocomplex(v);
            break;
        default:
            MSG_ERROR("The desired data type is not known or supported\n");
            return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_totype  ----- */


