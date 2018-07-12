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
 * @file lib/vector/any.c
 * @brief Anyfunction of @ref mess_vector structures.
 * @author @mbehr
 *
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
 *
 * @brief Checks a @ref mess_vector  if a given predicate is true.
 * @param[in,out] v     input/output @ref mess_vector \f$v\f$
 * @param[in] f_real    input scalar real function to apply to each entry of \f$v\f$ with @ref mess_int_t return value
 * @param[in] f_cpx     input scalar complex function to apply to each entry of \f$v\f$ with @ref mess_int_t return value
 * @param[in,out] anyval input/output pointer to a  @ref mess_int_t
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_any function tests entries of @p v with @p f_real or @p f_cpx. \n
 * If @p v is @ref MESS_REAL @p f_real is used.   \n
 * If @p v is @ref MESS_COMPLEX @p f_cpx is used. \n
 * @p f_real and @p f_cpx are scalar functions with real/complex argument and @ref mess_int_t return value. \n
 * If @p f_real / @p f_cpx applied to any entry \f$ v_{i}\f$ returns a nonzero value @c 1  @p anyval is set to
 * @p  1 and the @ref mess_vector any terminates.
 * If @p f_real / @p f_cpx returns @c 0 for every entry \f$v_{i}\f$ @ref mess_vector_any set @p anyval to @c 0 and terminates.
 * The @ref mess_vector_any function works in a similiar fashion like the @p any function of @matlab.
 *
 */
int mess_vector_any(mess_vector v, mess_int_t (*f_real) (double), mess_int_t (*f_cpx) (mess_double_cpx_t), mess_int_t* anyval) {
    MSG_FNAME(__func__);
    mess_int_t i=0;

    /*-----------------------------------------------------------------------------
     * check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    if(MESS_IS_REAL(v)){
        if(!f_real){
            MSG_ERROR("Please provide a real scalar function as input argument.");
            return MESS_ERROR_ARGUMENTS;
        }
    }

    if(MESS_IS_COMPLEX(v)){
        if(!f_cpx){
            MSG_ERROR("Please provide a complex scalar function as input argument.");
            return MESS_ERROR_ARGUMENTS;
        }
    }
    mess_check_nullpointer(anyval);

    /*-----------------------------------------------------------------------------
     *  check predicat for each entry
     *-----------------------------------------------------------------------------*/
    *anyval = 0;
    if(MESS_IS_REAL(v)){
        for(i=0;i<v->dim;++i){
            if(f_real(v->values[i])){
                *anyval=1;
                return 0;
            }
        }
    }else{
        for(i=0;i<v->dim;++i){
            if(f_cpx(v->values_cpx[i])){
                *anyval=1;
                return 0;
            }
        }
    }

    return 0;
} /* -----  end of function mess_vector_any  ----- */


