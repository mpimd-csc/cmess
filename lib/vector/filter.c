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
 * @file lib/vector/filter.c
 * @brief Filter functions of @ref mess_vector structures.
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
 * @brief Filter values from vector.
 * @param[in,out]       in  input/output vector
 * @param[in]           filter_real function
 * @param[in]           filter_complex function
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref mess_vector_filter returns a vector constructed from values fullfilling condition
 * given by the filter function pointer @p filter_real and @p filter_complex.
 *
 */
int mess_vector_filter(mess_vector in, int (*filter_real)(const double *), int (*filter_complex)(const mess_double_cpx_t *) ){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t i=0,j=0;

    /*-----------------------------------------------------------------------------
     *  check input data
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in);
    mess_check_real_or_complex(in);

    /*-----------------------------------------------------------------------------
     *  filter instable values
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(in)){
        i=1,j=0;
        for(i=0;i<in->dim;++i){
            if((*filter_real)(in->values+i)){
                in->values[j]=in->values[i];
                ++j;
            }
        }
    }else{
        for(i=0;i<in->dim;++i){
            if((*filter_complex)(in->values_cpx+i)){
                in->values_cpx[j]=in->values_cpx[i];
                ++j;
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  resize vector
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_resize(in,j);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
    return 0;
}



static int __filter_stable_real(const double * val){return (*val)<0;};
static int __filter_stable_complex(const mess_double_cpx_t* val){return creal( *val )<0;};

/**
 * @brief Filter values with positive realpart from vector.
 * @param[in,out]       in  input vector
 * @return zero on success or a non zero error value otherwise
 *
 * The @ref  mess_vector_filter_stable deletes all values with nonnegative realpart from vector and return the resized one.
 *
 */
int mess_vector_filter_stable(mess_vector in){
    return mess_vector_filter(in, __filter_stable_real, __filter_stable_complex);
}




