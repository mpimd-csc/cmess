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
 * @file include/mess/vector.h
 * @brief Interface to handle vectors.
 * @author @koehlerm
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include "mess/config.h"
#include "mess/matrix.h"

#ifdef __cplusplus
extern "C" {
#endif
    /** @addtogroup vector_ds
     * @{ */

    /**
     * @brief Representation of a real or complex vector.
     *
     * The @ref mess_vector_st structure holds all information about a vector.
     */
    typedef struct mess_vector_st {
        double *values;                 /**< Real values of a vector*/
        mess_double_cpx_t *values_cpx;  /**< Complex values of a vector */
        mess_int_t dim;                 /**< Dimension of a vector */
        mess_datatype_t data_type;      /**< Datatype of the values in a vector */
    } mess_vector_st;

    /**
     * @brief Type definition of a real or complex vector.
     *
     * The @ref mess_vector type definition defines a @ref mess_vector in order to use it easier.
     */
    typedef mess_vector_st*  mess_vector;
    int mess_vector_init(mess_vector *vect);
    int mess_vector_alloc(mess_vector vect, mess_int_t dim, mess_datatype_t dtype);
    int mess_vector_reset(mess_vector vect);
    int mess_vector_clear(mess_vector *vect);
    int mess_vector_copy(mess_vector x, mess_vector y);
    int mess_vector_copy_tocomplex(mess_vector in, mess_vector out);

    int mess_vector_tocomplex(mess_vector v);
    int mess_vector_toreal(mess_vector v);
    int mess_vector_toreal_nowarn(mess_vector v);
    int mess_vector_totype(mess_vector v, mess_datatype_t dt );
    int mess_vector_from_farray(mess_vector v, mess_int_t dim, double * vals, mess_double_cpx_t * vals_cpx);
    int mess_vector_from_lapack(mess_vector v, mess_int_t dim, double *vals_re, double *vals_im);
    int mess_vector_realpart(mess_vector in, mess_vector out);
    int mess_vector_imagpart(mess_vector in, mess_vector out);
    int mess_vector_complex_from_parts(mess_vector xr, mess_vector xc, mess_vector x);

    /** @}  */

    /** @addtogroup vector_misc
     * @{ */
    int mess_vector_print(mess_vector vect);
    int mess_vector_printshort(mess_vector vect);
    int mess_vector_write(const char *filename, mess_vector vect);
    int mess_vector_read(const char *filename, mess_vector vect);
    int mess_vector_printinfo(mess_vector v );
    int mess_vector_resize(mess_vector x, mess_int_t dim);
    int mess_vector_zeros(mess_vector v);
    int mess_vector_ones(mess_vector v);
    int mess_vector_rand(mess_vector v);
    int mess_vector_get(mess_vector v, mess_int_t i, mess_double_cpx_t* val);
    int mess_vector_rand_init(mess_int_t *seed);

    int mess_vector_iperm(mess_vector in, mess_int_t *iperm, mess_vector out);
    int mess_vector_perm(mess_vector in, mess_int_t *perm, mess_vector out);
    int mess_vector_perm_inplace(mess_vector v, mess_int_t *p);
    int mess_vector_iperm_inplace(mess_vector v, mess_int_t *q);
    int mess_vector_iperm_split(mess_vector in, mess_int_t *iperm, mess_vector reout, mess_vector imout);
    int mess_vector_iperm_combine(mess_vector rein, mess_vector imin , mess_int_t *iperm, mess_vector out);
    int mess_vector_perm_split(mess_vector in, mess_int_t *perm, mess_vector reout, mess_vector imout);
    int mess_vector_perm_combine(mess_vector rein, mess_vector imin , mess_int_t *perm, mess_vector out);

    double *mess_vector_valuereal(mess_vector vec, mess_int_t index);
    mess_double_cpx_t *mess_vector_valuecomplex(mess_vector vec, mess_int_t index);

    int mess_vector_sort(mess_vector vect);
    int mess_vector_sort_realpart(mess_vector vect);
    int mess_vector_sort_imagpart(mess_vector vect);
    int mess_vector_logspace10(mess_vector vect, double a, double b, mess_int_t nsample);
    int mess_vector_logspacee(mess_vector vect, double a, double b, mess_int_t nsample);
    int mess_vector_logspace2(mess_vector vect, double a, double b, mess_int_t nsample);
    int mess_vector_linspace(mess_vector vect, double a, double b, mess_int_t nsample);

    int mess_vector_cat(mess_vector x1, mess_vector x2, mess_vector x);
    int mess_vector_split(mess_vector input, mess_int_t n, mess_vector x1, mess_vector x2 );
    int mess_vector_lift(mess_vector in, mess_int_t n, mess_vector out);
    int mess_vector_kron(mess_vector in1, mess_vector in2,  mess_vector out );

    int mess_vector_filter_stable(mess_vector in);
    int mess_vector_filter(mess_vector in, int (*filter_real)(const double *), int (*filter_complex)(const mess_double_cpx_t *));
    int mess_vector_convert_if_real(mess_vector v);

    int mess_vector_map(mess_vector v, double (*f_real)(double), mess_double_cpx_t (*f_cpx)(mess_double_cpx_t));
    int mess_vector_map_abs(mess_vector);
    int mess_vector_map_acos(mess_vector);
    int mess_vector_map_acosh(mess_vector);
    int mess_vector_map_arg(mess_vector);
    int mess_vector_map_asin(mess_vector);
    int mess_vector_map_asinh(mess_vector);
    int mess_vector_map_atan(mess_vector);
    int mess_vector_map_atanh(mess_vector);
    int mess_vector_map_ceil(mess_vector);
    int mess_vector_map_conj(mess_vector);
    int mess_vector_map_cos(mess_vector);
    int mess_vector_map_cosh(mess_vector);
    int mess_vector_map_exp(mess_vector);
    int mess_vector_map_expm1(mess_vector);
    int mess_vector_map_floor(mess_vector);
    int mess_vector_map_log(mess_vector);
    int mess_vector_map_not(mess_vector);
    int mess_vector_map_round(mess_vector);
    int mess_vector_map_sin(mess_vector);
    int mess_vector_map_sinh(mess_vector);
    int mess_vector_map_sqrt(mess_vector);
    int mess_vector_map_tan(mess_vector);
    int mess_vector_map_tanh(mess_vector);
    int mess_vector_map_isfinite(mess_vector);
    int mess_vector_map_isinf(mess_vector);
    int mess_vector_map_isnan(mess_vector);

    int mess_vector_any(mess_vector mat, mess_int_t (*f_real) (double), mess_int_t (*f_cpx) (mess_double_cpx_t), mess_int_t* anyval);
    mess_int_t mess_vector_memsize(mess_vector v);
    /** @}  */

    /** @addtogroup vector_op
     * @{ */

    int mess_vector_axpy(double a, mess_vector x, mess_vector y);
    int mess_vector_norm2(mess_vector x, double* nrm);
    int mess_vector_norm1(mess_vector x, double* nrm);
    int mess_vector_norminf(mess_vector x, double* nrm);
    int mess_vector_norm(mess_vector x, mess_norm_t nrm_t, double* nrm);

    int mess_vector_scale(double alpha, mess_vector x);
    int mess_vector_dot(mess_vector x, mess_vector y, double* dot);

    int mess_vector_minvalue(mess_vector v, double *min);
    int mess_vector_maxvalue(mess_vector v, double *max);

    int mess_vector_scalec( mess_double_cpx_t alpha, mess_vector x);
    int mess_vector_dotc(mess_vector x, mess_vector y, mess_double_cpx_t* dot);
    int mess_vector_axpyc(mess_double_cpx_t a, mess_vector x, mess_vector y);
    int mess_vector_dotu(mess_vector x, mess_vector y, mess_double_cpx_t * dot);

    int mess_vector_diffnorm(mess_vector x1, mess_vector x2, double *diff);
    int mess_vector_diffnorminf(mess_vector x1, mess_vector x2, double *diff);

    int mess_vector_max(mess_vector v, double * maxval, mess_int_t *maxind);
    int mess_vector_conj(mess_vector vector);

    /*-----------------------------------------------------------------------------
     *  some defines
     *-----------------------------------------------------------------------------*/
    /**
     * @brief Macro to scale a vector.
     *
     * The @ref mess_vector_scalee macro scales a real or complex vector with a value.
     *
     **/
#define mess_vector_scalee(val,vec) if (MESS_IS_REAL(vec)) { mess_vector_scale((val),(vec)); } else { mess_vector_scalec((val), (vec));}

    /** @} */



#ifdef __cplusplus
}
#endif
#endif
