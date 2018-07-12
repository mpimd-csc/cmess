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
 * @file lib/vector/map.c
 * @brief Map functions of @ref mess_vector structures.
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
 *
 * @brief Map a scalar function to a vector.
 * @param[in,out] v     input/output vector \f$v\f$
 * @param[in] f_real    input scalar real function to apply to each entry of \f$v\f$
 * @param[in] f_cpx     input scalar complex function to apply to each entry of \f$v\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_map function is a higher order function, which applies
 * a scalar real @p f_real or a scalar complex @p f_cpx function to a
 * each entry of a vector \f$v\f$. Depending on the data type of \f$v\f$
 * we perform the following operation:
 * <center>
 * |        data type        | operation                                |
 * |:-----------------------:|:----------------------------------------:|
 * | @ref MESS_REAL          | \f$ v(i)\leftarrow f_{real}(v(i)) \f$    |
 * | @ref MESS_COMPLEX       | \f$ v(i)\leftarrow f_{cpx}(v(i))  \f$    |
 * </center>
 *
 */
int mess_vector_map(mess_vector v, double (*f_real)(double), mess_double_cpx_t (*f_cpx)(mess_double_cpx_t)) {
    mess_int_t i=0;
    MSG_FNAME(__func__);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    if(MESS_IS_REAL(v)){
        if(!f_real){
            MSG_ERROR("Please provide a real scalar function as input argument.");
            return MESS_ERROR_ARGUMENTS;
        }
        for(i=0;i<v->dim;++i){
            v->values[i] = f_real(v->values[i]);
        }
    }else{
        if(!f_cpx){
            MSG_ERROR("Please provide a complex scalar function as input argument.");
            return MESS_ERROR_ARGUMENTS;
        }
        for(i=0;i<v->dim;++i){
            v->values_cpx[i] = f_cpx(v->values_cpx[i]);
        }
    }

    return 0;
} /* -----  end of function mess_vector_map  ----- */


#undef MESS_MAP_STATIC_FUNC
#define MESS_MAP_STATIC_FUNC(BASENAME,REALFUNC,CPXFUNC)                                         \
                                                                                                \
    static double               map_##BASENAME##_real(double d){return REALFUNC(d);}            \
    static mess_double_cpx_t    map_##BASENAME##_cpx(mess_double_cpx_t z){return CPXFUNC(z);}   \
    int mess_vector_map_##BASENAME (mess_vector v){                                             \
        return mess_vector_map(v,map_##BASENAME##_real,map_##BASENAME##_cpx);                   \
    }                                                                                           \

#undef MESS_MAP_STATIC_FUNC2
#define MESS_MAP_STATIC_FUNC2(BASENAME,REALFUNC)                                                                            \
    static double               map_##BASENAME##_real(double d){return REALFUNC(d);}                                        \
    static mess_double_cpx_t    map_##BASENAME##_cpx(mess_double_cpx_t z){return REALFUNC(creal(z))+I*REALFUNC(cimag(z));}  \
    int mess_vector_map_##BASENAME (mess_vector v){                                                                         \
        return mess_vector_map(v,map_##BASENAME##_real,map_##BASENAME##_cpx);                                               \
    }                                                                                                                       \


#undef MESS_MAP_STATIC_FUNC3
#define MESS_MAP_STATIC_FUNC3(BASENAME,REALFUNC)                                                                                        \
    static double               map_##BASENAME##_real(double d){return REALFUNC(d);}                                                    \
    static mess_double_cpx_t    map_##BASENAME##_cpx(mess_double_cpx_t z){return abs(REALFUNC(creal(z))) && abs(REALFUNC(cimag(z)));}   \
    int mess_vector_map_##BASENAME (mess_vector v){                                                                                     \
        return mess_vector_map(v,map_##BASENAME##_real,map_##BASENAME##_cpx);                                                           \
    }                                                                                                                                   \



MESS_MAP_STATIC_FUNC(abs,fabs,cabs);
MESS_MAP_STATIC_FUNC(acos,acos,cacos);
MESS_MAP_STATIC_FUNC(acosh,acosh,cacosh);
MESS_MAP_STATIC_FUNC(arg,carg,carg);
MESS_MAP_STATIC_FUNC(asin,asin,casin);
MESS_MAP_STATIC_FUNC(asinh,asinh,casinh);
MESS_MAP_STATIC_FUNC(atan,atan,catan);
MESS_MAP_STATIC_FUNC(atanh,atanh,catanh);
MESS_MAP_STATIC_FUNC2(ceil,ceil);
MESS_MAP_STATIC_FUNC(conj,/* identity */,conj);
MESS_MAP_STATIC_FUNC(cos,cos,ccos);
MESS_MAP_STATIC_FUNC(cosh,cosh,ccosh);
MESS_MAP_STATIC_FUNC(exp,exp,cexp);
MESS_MAP_STATIC_FUNC(expm1,expm1,-1+expm1);
MESS_MAP_STATIC_FUNC2(floor,floor);
MESS_MAP_STATIC_FUNC(log,log,clog);
MESS_MAP_STATIC_FUNC2(round,round);
MESS_MAP_STATIC_FUNC(sin,sin,csin);
MESS_MAP_STATIC_FUNC(sinh,sinh,csinh);
MESS_MAP_STATIC_FUNC(sqrt,sqrt,csqrt);
MESS_MAP_STATIC_FUNC(tan,tan,ctan);
MESS_MAP_STATIC_FUNC(tanh,tanh,ctanh);
MESS_MAP_STATIC_FUNC3(isfinite,isfinite);
MESS_MAP_STATIC_FUNC3(isinf,isinf);
MESS_MAP_STATIC_FUNC3(isnan,isnan);


static double               map_not_real(double d){return d?0:1;}
static mess_double_cpx_t    map_not_cpx(mess_double_cpx_t z){return z?0:1;}
int mess_vector_map_not (mess_vector v){
    return mess_vector_map(v,map_not_real,map_not_cpx);
}









