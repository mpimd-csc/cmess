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
 * @file lib/vector/kron.c
 * @brief Kronecker product function of @ref mess_vector structures.
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
 * @brief Compute Kronecker product for two vectors.
 * @param[in]       in1  input vector \f$ in_1\f$
 * @param[in]       in2  input vector \f$ in_2\f$
 * @param[out]      out output vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_kron preforms a Kronecker product for  \f$ in_1  \f$ and  \f$ in_2  \f$, i.e.
 * \f[ \begin{array}{cccc}
 * out &=& \left[ \begin{array}{c} in_{1_1} \\ \vdots \\ in_{1_n} \end{array} \right] \otimes &
 * \left[ \begin{array}{ccc} in_{2_1} & \hdots & in_{2_m} \end{array} \right] \\
 * &=& \left[ \begin{array}{c} in_{1_1} in_2 \\ \vdots \\ in_{1_n} in_2 \end{array} \right] &
 * \end{array}\f]
 * The result is a vector \f$ out  \f$.
 *
 */
int mess_vector_kron(mess_vector in1, mess_vector in2,  mess_vector out ){
    MSG_FNAME(__func__);
    int ret=0;
    double alpha=1.0;
    mess_double_cpx_t calpha = 1.0;
    mess_int_t one=1;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in1);
    mess_check_nullpointer(in2);
    mess_check_real_or_complex(in1);
    mess_check_real_or_complex(in2);
    mess_check_nullpointer(out);

    /*-----------------------------------------------------------------------------
     *  perform in1 "kronecker" in2 -> out
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_COMPLEX(in1)|| MESS_IS_COMPLEX(in2)){
        ret = mess_vector_tocomplex(out);                                                                                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_tocomplex);
    }else{
        ret = mess_vector_toreal_nowarn(out);                                                                                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
    }
    ret = mess_vector_zeros(out);                                                                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_zeros);
    ret = mess_vector_resize(out,in1->dim*in2->dim);                                                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);


    if(MESS_IS_REAL(in1) && MESS_IS_REAL(in2)){
        F77_GLOBAL(dger,DGER)(&(in2->dim),&(in1->dim),&alpha, in2->values, &one, in1->values,&one,out->values, &(in2->dim));
    }else if(MESS_IS_REAL(in1)){
        mess_vector t1;
        ret = mess_vector_init(&t1);                                                                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(t1,in1->dim,MESS_COMPLEX);                                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_copy_tocomplex(in1,t1);                                                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy_tocomplex);
        F77_GLOBAL(zgeru,ZGERU)(&(in2->dim),&(t1->dim),&calpha, in2->values_cpx, &one, t1->values_cpx,&one,out->values_cpx, &(in2->dim));
        ret = mess_vector_clear(&t1);                                                                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_clear);
    }else if (MESS_IS_REAL(in2)){
        mess_vector t2;
        ret = mess_vector_init(&t2);                                                                                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(t2,in2->dim,MESS_COMPLEX);                                                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_copy_tocomplex(in2,t2);                                                                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_copy_tocomplex);
        F77_GLOBAL(zgeru,ZGERU)(&(t2->dim),&(in1->dim),&calpha, t2->values_cpx, &one, in1->values_cpx,&one,out->values_cpx, &(t2->dim));
        ret = mess_vector_clear(&t2);                                                                                           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_clear);
    }else{
        F77_GLOBAL(zgeru,ZGERU)(&(in2->dim),&(in1->dim),&calpha, in2->values_cpx, &one, in1->values_cpx,&one,out->values_cpx, &(in2->dim));
    }

    return 0;
}






