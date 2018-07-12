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
 * @file lib/matrix/map.c
 * @brief Map a scalar function to a matrix.
 * @author @mbehr
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 *
 * @brief Map a scalar function to a matrix.
 * @param[in,out] mat   input/output matrix \f$mat\f$
 * @param[in] f_real    input scalar real function to apply to each entry of \f$mat\f$
 * @param[in] f_cpx     input scalar complex function to apply to each entry of \f$mat\f$
 * @param[in] full      input scalar if nonzero  apply functions to each  entry of \f$mat\f$, otherwise only to nonzero entries of \f$ mat\f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_map function is a higher order function, which applies
 * a scalar real @p f_real or a scalar complex @p f_cpx function to a
 * each entry of a matrix \f$mat\f$. Depending on the data type of \f$mat\f$
 * we perform the following operation:
 * <center>
 * |        data type        | operation                                        |
 * |:-----------------------:|:------------------------------------------------:|
 * | @ref MESS_REAL          | \f$ mat(i,j)\leftarrow f_{real}(mat(i,j)) \f$    |
 * | @ref MESS_COMPLEX       | \f$ mat(i,j)\leftarrow f_{cpx}(mat(i,j))  \f$    |
 * </center>
 *
 * If @p full is zero then the functions @p f_real and @p f_cpx are only applied to the
 * nonzero entries of @p mat.
 * If @p full is nonzero the the function @p f_real and @p f_cpx are applied to all entries
 * of @p mat. If @p mat is a sparse matrix @ref mess_matrix_map may convert the matrix to a
 * dense one due to the possible fill-in.
 */
int mess_matrix_map(mess_matrix mat, double (*f_real)(double), mess_double_cpx_t (*f_cpx)(mess_double_cpx_t), int full) {
    mess_int_t i,j=0;
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);
    if(MESS_IS_REAL(mat) && !f_real){
        MSG_ERROR("Please provide a real scalar function as input argument.");
        return MESS_ERROR_ARGUMENTS;
    }

    if(MESS_IS_COMPLEX(mat) && !f_cpx){
        MSG_ERROR("Please provide a real scalar function as input argument.");
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  check if we need a conversion to a dense matrix
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_SPARSE(mat) && full && ((MESS_IS_REAL(mat) && f_real(0.0)!=0.0) || (MESS_IS_COMPLEX(mat) && f_real(0.0)!=0.0))){
        int ret = 0;
        mess_matrix mat2;
        ret = mess_matrix_init(&mat2);                          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_convert(mat,mat2,MESS_DENSE);         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_convert);
        ret = mess_matrix_copy(mat2,mat);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_clear(&mat2);                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
    }


    /*-----------------------------------------------------------------------------
     *  apply function to a matrix
     *-----------------------------------------------------------------------------*/
    if(MESS_IS_REAL(mat)){

        if(MESS_IS_DENSE(mat)){
            for(i=0;i<mat->cols;++i){
                for(j=0;j<mat->rows;++j){
                    mat->values[j+i*mat->ld]=f_real(mat->values[j+i*mat->ld]);
                }
            }
        }else{
            for(i=0;i<mat->nnz;++i){
                mat->values[i]=f_real(mat->values[i]);
            }
        }
    }else{
        if(MESS_IS_DENSE(mat)){
            for(i=0;i<mat->cols;++i){
                for(j=0;j<mat->rows;++j){
                    mat->values_cpx[j+i*mat->ld]=f_cpx(mat->values_cpx[j+i*mat->ld]);
                }
            }
        }else{
            for(i=0;i<mat->nnz;++i){
                mat->values_cpx[i]=f_cpx(mat->values_cpx[i]);
            }
        }
    }

    return 0;
} /* -----  end of function mess_matrix_map  ----- */


#define MESS_MAP_STATIC_FUNC(BASENAME,REALFUNC,CPXFUNC)                                         \
    static double               map_##BASENAME##_real(double d){return REALFUNC(d);}            \
    static mess_double_cpx_t    map_##BASENAME##_cpx(mess_double_cpx_t z){return CPXFUNC(z);}   \
    int mess_matrix_map_##BASENAME (mess_matrix mat, int full){                                 \
        return mess_matrix_map(mat,map_##BASENAME##_real,map_##BASENAME##_cpx, full);           \
    }                                                                                           \

#define MESS_MAP_STATIC_FUNC2(BASENAME,REALFUNC)                                                                            \
    static double               map_##BASENAME##_real(double d){return REALFUNC(d);}                                        \
    static mess_double_cpx_t    map_##BASENAME##_cpx(mess_double_cpx_t z){return REALFUNC(creal(z))+I*REALFUNC(cimag(z));}  \
    int mess_matrix_map_##BASENAME (mess_matrix mat, int full){                                                             \
        return mess_matrix_map(mat,map_##BASENAME##_real,map_##BASENAME##_cpx,full);                                        \
    }                                                                                                                       \

#define MESS_MAP_STATIC_FUNC3(BASENAME,REALFUNC)                                                                                        \
    static double               map_##BASENAME##_real(double d){return REALFUNC(d);}                                                    \
    static mess_double_cpx_t    map_##BASENAME##_cpx(mess_double_cpx_t z){return abs(REALFUNC(creal(z))) && abs(REALFUNC(cimag(z)));}   \
    int mess_matrix_map_##BASENAME (mess_matrix mat, int full){                                                                         \
        return mess_matrix_map(mat,map_##BASENAME##_real,map_##BASENAME##_cpx,full);                                                    \
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
int mess_matrix_map_not (mess_matrix mat, int full){
    return mess_matrix_map(mat,map_not_real,map_not_cpx,full);
}




