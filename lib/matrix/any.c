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
 * @file lib/matrix/any.c
 * @brief Any function clone from @matlab.
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
 * @brief Checks a matrix rowwise/columnwise if a given predicate is true.
 * @param[in,out] mat   input/output matrix \f$mat\f$
 * @param[in] f_real    input scalar real function to apply to each entry of \f$mat\f$ with @ref mess_int_t return value
 * @param[in] f_cpx     input scalar complex function to apply to each entry of \f$mat\f$ with @ref mess_int_t return value
 * @param[in] dim       input @p 0 for rowise fashion @p 1 for columnwise fashion.
 * @param[in,out] anyvec input/output @ref mess_vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_any function tests entries of @p mat with @p f_real or @p f_cpx in a rowwise/columnwise fashion. \n
 * If @p dim is equal to @c 0 rowwise fashion is used. \n
 * If @p dim is equal to @c 1 columnwise fashion is used. \n
 * If @p mat is @ref MESS_REAL @p f_real is used.   \n
 * If @p mat is @ref MESS_COMPLEX @p f_cpx is used. \n
 * @p f_real and @p f_cpx are scalar functions with real/complex argument and @ref mess_int_t return value. \n
 * If @p f_real / @p f_cpx applied to the entry \f$ mat_{i,j}\f$ returns a nonzero value @c 1 is written to the @p j -th
 * entry of @p anyvec.  \n
 * Independent of the datatype of @p mat @p anyvec is in every case  a @ref MESS_REAL @ref mess_vector.\n
 * If @p dim is equal to @c 0 the number of entries of @p anyvec is equal to the number of rows of @p mat. \n
 * If @p dim is equal to @c 1 the number of entries of @p anyvec is equal to the number of cols of @p mat. \n
 * The @ref mess_matrix_any function works in a similiar fashion like the @p any function of @matlab.
 *
 */
int mess_matrix_any(mess_matrix mat, mess_int_t (*f_real) (double), mess_int_t (*f_cpx) (mess_double_cpx_t),mess_int_t dim, mess_vector anyvec) {
    MSG_FNAME(__func__);
    mess_int_t i=0,j=0,ret=0;

    /*-----------------------------------------------------------------------------
     * check input arguments
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_real_or_complex(mat);

    if(MESS_IS_REAL(mat)){
        if(!f_real){
            MSG_ERROR("Please provide a real scalar function as input argument.");
            return MESS_ERROR_ARGUMENTS;
        }
    }

    if(MESS_IS_COMPLEX(mat)){
        if(!f_cpx){
            MSG_ERROR("Please provide a complex scalar function as input argument.");
            return MESS_ERROR_ARGUMENTS;
        }
    }
    mess_check_nullpointer(anyvec);

    if(dim!=0 && dim!=1){
        MSG_ERROR("Please specify the fashion for traversing with 0 for rowwise or 1 for columnwise. \n");
    }

    /*-----------------------------------------------------------------------------
     * resize vector and zero the entries
     *-----------------------------------------------------------------------------*/
    if(dim==0){
        //rowwise fashion
        ret = mess_vector_resize(anyvec,mat->rows);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
    }else{
        //columnwise fashion
        ret = mess_vector_resize(anyvec,mat->cols);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_resize);
    }
    ret = mess_vector_toreal_nowarn(anyvec);                    FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_toreal_nowarn);
    ret = mess_vector_zeros(anyvec);                            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_zeros);

    if(MESS_IS_DENSE(mat)){
        /*-----------------------------------------------------------------------------
         *  case MESS_DENSE
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(mat)){
            for(i=0;i<mat->rows;++i){
                for(j=0;j<mat->cols;++j){
                    if(f_real(mat->values[i+j*mat->ld])){
                        anyvec->values[dim?j:i]=1;
                    }
                }
            }
        }else{
            for(i=0;i<mat->rows;++i){
                for(j=0;j<mat->cols;++j){
                    if(f_cpx(mat->values_cpx[i+j*mat->ld])){
                        anyvec->values[dim?j:i]=1;
                    }
                }
            }
        }
    }else if (MESS_IS_CSR(mat)){
        /*-----------------------------------------------------------------------------
         *  case MESS_CSR
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(mat)){
            for(i=0;i<mat->rows;++i){
                for(j=mat->rowptr[i];j<mat->rowptr[i+1];j++){
                    if(f_real(mat->values[j])){
                        anyvec->values[dim?mat->colptr[j]:i]=1;
                    }
                }
            }
        }else{
            for(i=0;i<mat->rows;++i){
                for(j=mat->rowptr[i];j<mat->rowptr[i+1];j++){
                    if(f_cpx(mat->values_cpx[j])){
                        anyvec->values[dim?mat->colptr[j]:i]=1;
                    }
                }
            }
        }
    }else if (MESS_IS_CSC(mat)){
        /*-----------------------------------------------------------------------------
         *  case MESS_CSC
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(mat)){
            for(i=0;i<mat->cols;++i){
                for(j=mat->colptr[i];j<mat->colptr[i+1];++j){
                    if(f_real(mat->values[j])){
                        anyvec->values[dim?i:mat->rowptr[j]]=1;
                    }
                }
            }
        }else{
            for(i=0;i<mat->cols;++i){
                for(j=mat->colptr[i];j<mat->colptr[i+1];++j){
                    if(f_cpx(mat->values_cpx[j])){
                        anyvec->values[dim?i:mat->rowptr[j]]=1;
                    }
                }
            }
        }
    }else if (MESS_IS_COORD(mat)){
        /*-----------------------------------------------------------------------------
         *  case MESS_COORD
         *-----------------------------------------------------------------------------*/
        if(MESS_IS_REAL(mat)){
            for(i=0;i<mat->nnz;++i){
                if(f_real(mat->values[i])){
                    anyvec->values[dim?mat->colptr[i]:mat->rowptr[i]]=1;
                }
            }
        }else{
            for(i=0;i<mat->nnz;++i){
                if(f_cpx(mat->values_cpx[i])){
                    anyvec->values[dim?mat->colptr[i]:mat->rowptr[i]]=1;
                }
            }
        }
    }else{
        /*-----------------------------------------------------------------------------
         *  unknown storage type
         *-----------------------------------------------------------------------------*/
        MSG_ERROR("Unknown storage type.\n");
        return MESS_ERROR_STORAGETYPE;
    }

    return 0;
} /* -----  end of function mess_matrix_any  ----- */




