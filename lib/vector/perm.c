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
 * @file lib/vector/perm.c
 * @brief Permutation functions of @ref mess_vector structures.
 * @author @koehlerm
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
 * @internal
 * @brief Compare two complex numbers for @c qsort.
 * @param[in] pc1 input pointer to the first number
 * @param[in] pc2 input pointer to the second number
 * @return -1, 0, 1 if @c cabs(*pc1) is smaller, equal, bigger than  @c cabs(*pc2)
 *
 * The @ref __compare_complex function compares two complex numbers
 * for the @c stdlibc @c qsort function.
 * @attention Internal use only.
 */
static int __compare_complex(const void *pc1, const void *pc2) {
    mess_double_cpx_t
        *c1 = (mess_double_cpx_t *) pc1;
    mess_double_cpx_t
        *c2 = (mess_double_cpx_t *) pc2;

    if (cabs(*c1) == cabs(*c2) && carg(*c1) == carg(*c2))
        return 0;
    if (cabs(*c1) < cabs(*c2))
        return -1;
    if (cabs(*c1) == cabs(*c2) && carg(*c1) < carg(*c2))
        return -1;
    return 1;
}

/**
 * @internal
 * @brief Compare two real numbers for @c qsort.
 * @param[in] pc1 input pointer to the first number
 * @param[in] pc2 input pointer to the second number
 * @return -1, 0, 1 if @c *pc1 is smaller, equal, bigger than @c *pc2
 *
 * The @ref __compare_real function compares two real numbers
 * for the @c stdlibc @c qsort function.
 * @attention Internal use only.
 */
static int __compare_real(const void *pc1, const void *pc2) {
    double *c1 = (double *) pc1;
    double *c2 = (double *) pc2;

    if (*c1 < *c2) {
        return -1;
    } else if (*c1 > *c2) {
        return 1;
    } else {
        return 0;
    }
    return 0;
}

/**
 * @internal
 * @brief Compare real part of two complex numbers for @c  qsort.
 * @param[in] pc1 input pointer to the first number
 * @param[in] pc2 input pointer to the second number
 * @return -1, 0, 1 if real part of @c *pc1 is smaller, equal, bigger than real part of  @c *pc2
 *
 * The @ref __compare_realpart function compares two complex numbers by real part
 * for the @c stdlibc @c qsort function.
 * @attention Internal use only.
 */
static int __compare_realpart(const void *pc1, const void *pc2) {
    mess_double_cpx_t *c1 = (mess_double_cpx_t *) pc1;
    mess_double_cpx_t *c2 = (mess_double_cpx_t *) pc2;

    if (creal(*c1) < creal(*c2)){
        return -1;
    } else if (creal(*c1) > creal(*c2)) {
        return 1;
    } else {
        return 0;
    }
    return 0;
}





/**
 *
 * @internal
 * @brief Compare imaginary part of two complex numbers for qsort.
 * @param[in] pc1 input pointer to the first number
 * @param[in] pc2 input pointer to the second number
 * @return -1, 0, 1 if imaginary part of @c *pc1 is smaller, equal, bigger than imaginary part of @c *pc2
 *
 * The @ref __compare_imagpart function compares two complex numbers by imaginary part
 * for the @c stdlibc @c qsort function.
 * @attention Internal use only.
 */
static int __compare_imagpart(const void *pc1, const void *pc2) {
    mess_double_cpx_t *c1 = (mess_double_cpx_t *) pc1;
    mess_double_cpx_t *c2 = (mess_double_cpx_t *) pc2;

    if (cimag(*c1) < cimag(*c2)){
        return -1;
    } else if (cimag(*c1) > cimag(*c2)) {
        return 1;
    } else {
        return 0;
    }
    return 0;
}






/**
 * @brief Permute of a vector.
 * @param[in] in        input vector
 * @param[in] perm      input   permutation
 * @param[in,out] out   input/output vector
 *
 * The @ref mess_vector_perm function permutes a vector, such that
 * \f[
 * out [ i ] = in [ perm[i] ].
 * \f]
 *
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_combine
 * @see mess_vector_perm_split
 * @see mess_vector_iperm
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_combine
 * @see mess_vector_iperm_split
 */
int mess_vector_perm(mess_vector in, mess_int_t *perm, mess_vector out) {
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    if (out->dim != in->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(out, in->dim); FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }

    if ( MESS_IS_REAL(in) && MESS_IS_REAL(out)) {
        for (i = 0; i < in->dim; i++)
            out->values[i] = in->values[perm ? perm[i] : i];
    } else if ( MESS_IS_COMPLEX(in) && MESS_IS_COMPLEX(out)) {
        for (i = 0; i < in->dim; i++)
            out->values_cpx[i] = in->values_cpx[perm ? perm[i] : i];
    } else if ( MESS_IS_COMPLEX(in) && MESS_IS_REAL(out)) {
        MSG_WARN("permute a complex vector to a real vector. Lost of Information possible.\n");
        for (i = 0; i < in->dim; i++)
            out->values[i] = creal(in->values_cpx[perm ? perm[i] : i]);
    } else if ( MESS_IS_REAL(in) && MESS_IS_COMPLEX(out)) {
        for (i = 0; i < in->dim; i++)
            out->values_cpx[i] = in->values[perm ? perm[i] : i];
    } else {
        MSG_ERROR("unknown/unsupported datatype\n");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
}

/**
 * @brief Permute a vector inplace.
 * @param[in,out] v input/output vector
 * @param[in] perm  input permutation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_perm_inplace function permutes a vector nearly in place
 * using a modified version @p dlapmt from @lapack.
 *  \f[
 * in [ i ] \leftarrow  in [ perm[i] ].
 * \f]
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_combine
 * @see mess_vector_perm_split
 * @see mess_vector_iperm
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_combine
 * @see mess_vector_iperm_split
 *
 */
int mess_vector_perm_inplace(mess_vector v, mess_int_t *perm) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    /*-----------------------------------------------------------------------------
     * shift to lapack indexing
     *-----------------------------------------------------------------------------*/
    mess_int_t i;
    mess_int_t forward=1;
    mess_int_t one=1;

    for(i=0;i<v->dim;++i){ perm[i]+=1;}

    /*-----------------------------------------------------------------------------
     *  permute
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(v)){
        //F77_GLOBAL(idlapmt,IDLAPMT)(&forward, &one, &(v->dim), v->values, &one,perm);
        F77_GLOBAL(dlapmt,DLAPMT)(&forward, &one, &(v->dim), v->values, &one,perm);
    }else{
        //F77_GLOBAL(izlapmt,IZLAPMT)(&forward, &one, &(v->dim), v->values_cpx, &one, perm);
        F77_GLOBAL(zlapmt,ZLAPMT)(&forward, &one, &(v->dim), v->values_cpx, &one, perm);
    }

    /*-----------------------------------------------------------------------------
     *  shift back
     *-----------------------------------------------------------------------------*/
    for(i=0;i<v->dim;++i){ perm[i]-=1;}

    return 0;
} /* -----  end of function mess_vector_perm_inplace  ----- */

/**
 * @brief Inverse permutation of a vector.
 * @param[in] in        input vector
 * @param[in] iperm     input permutation
 * @param[in,out] out   input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_iperm function fullfills an inverse permutation, i.e.
 * \f[
 * out [ perm [i] ] = in [ i ].
 * \f]
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_combine
 * @see mess_vector_perm_split
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_combine
 * @see mess_vector_iperm_split
 */
int mess_vector_iperm(mess_vector in, mess_int_t *iperm, mess_vector out) {
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(in);
    mess_check_nullpointer(out);

    if (out->dim != in->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(out, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }

    if ( MESS_IS_REAL(in) && MESS_IS_REAL(out)) {
        for (i = 0; i < in->dim; i++)
            out->values[iperm ? iperm[i] : i] = in->values[i];
    } else if ( MESS_IS_COMPLEX(in) && MESS_IS_COMPLEX(out)) {
        for (i = 0; i < in->dim; i++)
            out->values_cpx[iperm ? iperm[i] : i] = in->values_cpx[i];
    } else if (MESS_IS_COMPLEX(in) && MESS_IS_REAL(out)) {
        MSG_WARN(
                "permute a complex vector to a real vector. Lost of Information possible.\n");
        for (i = 0; i < in->dim; i++)
            out->values[iperm ? iperm[i] : i] = creal(in->values_cpx[i]);
    } else if (MESS_IS_REAL(in) && MESS_IS_COMPLEX(out)) {
        for (i = 0; i < in->dim; i++)
            out->values_cpx[iperm ? iperm[i] : i] = in->values[i];
    } else {
        MSG_ERROR("unknown/unsupported datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

/**
 * @brief Inverse permutation of a vector and spilts real and imaginary part.
 * @param[in] in            input vector \f$in\f$
 * @param[in] iperm         input permutation
 * @param[out] reout        output real part of vector \f$in\f$
 * @param[out] imout        output imaginary part of vector \f$in\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_iperm_split functions fullfills an inverse permutation to a vector and splits it in
 * real and imaginary part, i.e.
 * \f[
 * \operatorname{Re}(out [ perm [i] ]) + \operatorname{Im} (out[ perm[i])) * I  = in [ i ].
 * \f]
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_combine
 * @see mess_vector_perm_split
 * @see mess_vector_iperm
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_combine
 */
int mess_vector_iperm_split(mess_vector in, mess_int_t *iperm,
        mess_vector reout, mess_vector imout) {
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(in);
    mess_check_nullpointer(reout);
    mess_check_nullpointer(imout);
    mess_check_complex(in);

    if (reout->dim != in->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(reout, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    if (imout->dim != in->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(imout, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    ret = mess_vector_toreal_nowarn(reout);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal_nowarn);
    ret = mess_vector_toreal_nowarn(imout);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal_nowarn);

    for (i = 0; i < in->dim; i++) {
        reout->values[iperm ? iperm[i] : i] = creal(in->values_cpx[i]);
        imout->values[iperm ? iperm[i] : i] = cimagl(in->values_cpx[i]);
    }
    return 0;
}

/**
 * @brief Inverse permutation of a vector and spilts real and imaginary part.
 * @param[in] rein              input real part of vector
 * @param[in] imin              input imaginary part of vector
 * @param[in] iperm             input permutation
 * @param[out] out              output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_iperm_combine function fullfills an inverse permutation to real and imaginary part and combines these
 * to a vector, i.e.
 * \f[
 * out [ perm [i] ]  = \operatorname{Re}( in [ i ] ) + \operatorname{Im}( in [i ] ) *I.
 * \f]
 *
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_combine
 * @see mess_vector_perm_split
 * @see mess_vector_iperm
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_split
 *
 *
 */
int mess_vector_iperm_combine(mess_vector rein, mess_vector imin,
        mess_int_t *iperm, mess_vector out) {
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(rein);
    mess_check_nullpointer(imin);
    mess_check_nullpointer(out);
    mess_check_real(rein);
    mess_check_real(imin);

    if (rein->dim != imin->dim) {
        MSG_ERROR(
                "The two input vectors don't have the same dimension rein->dim = " MESS_PRINTF_INT " \t imin->dim = " MESS_PRINTF_INT "\n",
                rein->dim, imin->dim);
        return MESS_ERROR_DIMENSION;
    }

    if (rein->dim != out->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(out, rein->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    ret = mess_vector_tocomplex(out);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);

    for (i = 0; i < rein->dim; i++) {
        out->values_cpx[iperm ? iperm[i] : i] = (rein->values[i])
            + (imin->values[i] * I);
    }
    return 0;
}




/**
 * @brief Permutate a vector and split in real and imaginary part.
 * @param[in] in            input vector \f$in\f$
 * @param[in] perm          input permutation
 * @param[out] reout        output real part of permuted vector \f$in\f$
 * @param[out] imout        output imaginary part of pemuted vector \f$in\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_perm_split function permutes a vector and splits this vector in real and imaginary part, i.e.
 * \f[
 * \operatorname{Re}(out [ i ]) + \operatorname{Im} (out[i])) * I  = in [ perm[i] ].
 * \f]
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_combine
 * @see mess_vector_iperm
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_combine
 * @see mess_vector_iperm_split
 */
int mess_vector_perm_split(mess_vector in, mess_int_t *perm, mess_vector reout,
        mess_vector imout) {
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(in);
    mess_check_nullpointer(reout);
    mess_check_nullpointer(imout);
    mess_check_complex(in);

    if (reout->dim != in->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(reout, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    if (imout->dim != in->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(imout, in->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    ret = mess_vector_toreal_nowarn(reout);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal_nowarn);
    ret = mess_vector_toreal_nowarn(imout);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_toreal_nowarn);

    for (i = 0; i < in->dim; i++) {
        reout->values[i] = creal(in->values_cpx[perm ? perm[i] : i]);
        imout->values[i] = cimag(in->values_cpx[perm ? perm[i] : i]);
    }
    return 0;
}



/**
 * @brief Inverse permutation of a vector and combination of real and imaginary part.
 * @param[in] rein              input real part of  vector
 * @param[in] imin              input imaginary part of vector
 * @param[in] perm              input permutation
 * @param[out] out              output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_perm_combine function fullfills an inverse permutation of the real part rein and the imaginary
 * part imin of a vector and combines these, i.e.
 * \f[
 * out [ i ]  = \operatorname{Re}(rein [perm[i]]) + \operatorname{Im}(imin [perm[i]]) *I.
 * \f]
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_split
 * @see mess_vector_iperm
 * @see mess_vector_iperm_inplace
 * @see mess_vector_iperm_combine
 * @see mess_vector_iperm_split
 */
int mess_vector_perm_combine(mess_vector rein, mess_vector imin,
        mess_int_t *perm, mess_vector out) {
    MSG_FNAME(__func__);
    mess_int_t i;
    int ret = 0;

    mess_check_nullpointer(rein);
    mess_check_nullpointer(imin);
    mess_check_nullpointer(out);
    mess_check_real(rein);
    mess_check_real(imin);

    if (rein->dim != imin->dim) {
        MSG_ERROR(
                "The two input vectors don't have the same dimension rein->dim = " MESS_PRINTF_INT " \t imin->dim = " MESS_PRINTF_INT "\n",
                rein->dim, imin->dim);
        return MESS_ERROR_DIMENSION;
    }

    if (rein->dim != out->dim) {
        MSG_WARN("out hasn't the right dimension. resize.\n");
        ret = mess_vector_resize(out, rein->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_resize);
    }
    ret = mess_vector_tocomplex(out);
    FUNCTION_FAILURE_HANDLE(ret, (ret != 0), mess_vector_tocomplex);

    for (i = 0; i < rein->dim; i++) {
        out->values_cpx[i] = (rein->values[perm ? perm[i] : i])
            + (imin->values[perm ? perm[i] : i] * I);
    }
    return 0;
}
/**
 * @brief Permute a vector inplace.
 * @param[in,out] v input/output vector
 * @param[in] perm  input permutation
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_iperm_inplace function permutes a vector nearly in place
 * using a modified version of @p dlapmt from @lapack.
 *  \f[
 * in [ perm[i] ] \leftarrow  in [ i ].
 * \f]
 *
 *
 * @see mess_vector_perm
 * @see mess_vector_perm_inplace
 * @see mess_vector_perm_combine
 * @see mess_vector_perm_split
 * @see mess_vector_iperm
 * @see mess_vector_iperm_combine
 * @see mess_vector_iperm_split
 */
int mess_vector_iperm_inplace(mess_vector v, mess_int_t *perm) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    /*-----------------------------------------------------------------------------
     * shift to lapack indexing
     *-----------------------------------------------------------------------------*/
    mess_int_t i;
    mess_int_t forward=0;
    mess_int_t one=1;
    for(i=0;i<v->dim;++i){ perm[i]+=1;}


    /*-----------------------------------------------------------------------------
     *  permute
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(v)){
        //F77_GLOBAL(idlapmt,IDLAPMT)(&forward, &one, &(v->dim), v->values, &one,perm);
        F77_GLOBAL(dlapmt,DLAPMT)(&forward, &one, &(v->dim), v->values, &one,perm);
    }else{
        //F77_GLOBAL(izlapmt,IZLAPMT)(&forward, &one, &(v->dim), v->values_cpx, &one, perm);
        F77_GLOBAL(zlapmt,ZLAPMT)(&forward, &one, &(v->dim), v->values_cpx, &one, perm);
    }

    /*-----------------------------------------------------------------------------
     *  shift back
     *-----------------------------------------------------------------------------*/
    for(i=0;i<v->dim;++i){ perm[i]-=1;}

    return 0;
} /* -----  end of function mess_vector_iperm_inplace  ----- */



/**
 * @brief Sort the components of a vector by real part in ascending order.
 * @param[in,out] vect input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_sort_realpart function sorts the values of a vector by real part in ascending order.
 *
 */
int mess_vector_sort_realpart(mess_vector vect) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(vect);

    if (MESS_IS_REAL(vect)) {
        qsort(vect->values, vect->dim, sizeof(double), __compare_real);
    } else if (MESS_IS_COMPLEX(vect)) {
        qsort(vect->values_cpx, vect->dim, sizeof(mess_double_cpx_t), __compare_realpart);
    } else {
        MSG_ERROR("unknown datatype.\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_sort_realpart  ----- */


/**
 * @brief Sort the components of a vector by imaginary part in ascending order.
 * @param[in,out] vect input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_sort_imagpart function sorts the values of a vector by imaginary part in ascending order.
 *
 */
int mess_vector_sort_imagpart(mess_vector vect) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(vect);

    if (MESS_IS_REAL(vect)) {
        return 0;
    } else if (MESS_IS_COMPLEX(vect)) {
        qsort(vect->values_cpx, vect->dim, sizeof(mess_double_cpx_t), __compare_imagpart);
    } else {
        MSG_ERROR("unknown datatype.\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_sort_imagpart  ----- */


/**
 * @brief Sort the components of a vector in ascending order.
 * @param[in,out] vect  input/output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_sort function sorts the values of a vector in ascending order.
 *
 */
int mess_vector_sort(mess_vector vect) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(vect);

    if (MESS_IS_REAL(vect)) {
        qsort(vect->values, vect->dim, sizeof(double), __compare_real);
    } else if (MESS_IS_COMPLEX(vect)) {
        qsort(vect->values_cpx, vect->dim, sizeof(mess_double_cpx_t), __compare_complex);
    } else {
        MSG_ERROR("unknown datatype.\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_sort  ----- */



