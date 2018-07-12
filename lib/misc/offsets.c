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
 * @file lib/misc/offsets.c
 * @brief Return offsets of structures defined in @mess.
 * @author @mbehr
 */


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/errors.h"
#include "mess/matrix.h"
#include "mess/vector.h"
#include "mess/error_macro.h"
#include <complex.h>


#define MESS_OFFSET(TYPE,FIELD,COUNTER,OFFSETPTR,FIELDNAMESPTR){                \
    if(OFFSETPTR){                                                              \
        OFFSETPTR[COUNTER] = (int) offsetof(TYPE,FIELD);                        \
    }                                                                           \
    if(FIELDNAMESPTR){                                                          \
        FIELDNAMESPTR[COUNTER] = #FIELD;                                        \
    }                                                                           \
    COUNTER=COUNTER+1;                                                          \
}


/**
 * @brief Return size of @ref mess_int_t.
 *
 * The @ref  mess_int_t_size returns the size of @ref mess_int_t datatype.
 */
int mess_int_t_size(void){
    return sizeof(mess_int_t);
}


/**
 * @brief Information about memory laout of @ref mess_matrix_st
 * @param[in,out] size          input/output contains @c sizeof  value
 * @param[in,out] offsets       input/output contains relative address of fields of @ref mess_matrix_st
 * @param[in,out] fieldnames    input/outout contains correspondig field names of @param offsets
 * @return number of fields of @ref mess_matrix_st
 *
 * The @ref mess_matrix_st_offset function returns informations about the memory layout of @ref mess_matrix_st.
 * We assume that @p offsets and @p fieldnames points to preallocated memory.
 * You can call the function also with @c NULL pointers, in this case the pointer is not used.
 *
 */
int mess_matrix_st_offset(int* size, int* offsets, char** fieldnames){

    int counter = 0;

    if(size){ *size = sizeof(struct mess_matrix_st);}

    MESS_OFFSET(struct mess_matrix_st,rows          ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,cols          ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,ld            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,nnz           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,rowptr        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,colptr        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,values        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,values_cpx    ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,store_type    ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,symmetry      ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_matrix_st,data_type     ,counter,offsets,fieldnames);

    return counter;
}


/**
 * @brief Information about memory laout of @ref mess_vector_st
 * @param[in,out] size          input/output contains @c sizeof  value
 * @param[in,out] offsets       input/output contains relative address of fields of @ref mess_vector_st
 * @param[in,out] fieldnames    input/outout contains correspondig field names of @param offsets
 * @return number of fields of @ref mess_vector_st
 *
 * The @ref mess_vector_st_offset function returns informations about the memory layout of @ref mess_vector_st.
 * We assume that @p offsets and @p fieldnames points to preallocated memory.
 * You can call the function also with @c NULL pointers, in this case the pointer is not used.
 *
 */
int mess_vector_st_offset(int* size, int* offsets, char** fieldnames){

    int counter = 0;

    if(size){ *size = sizeof(struct mess_vector_st);}

    MESS_OFFSET(struct mess_vector_st,values            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_vector_st,values_cpx        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_vector_st,dim               ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_vector_st,data_type         ,counter,offsets,fieldnames);

    return counter;
}


/**
 * @brief Information about memory laout of @ref mess_status_st
 * @param[in,out] size          input/output contains @c sizeof  value
 * @param[in,out] offsets       input/output contains relative address of fields of @ref mess_status_st
 * @param[in,out] fieldnames    input/outout contains correspondig field names of @param offsets
 * @return number of fields of @ref mess_status_st
 *
 * The @ref mess_status_st_offset function returns informations about the memory layout of @ref mess_status_st.
 * We assume that @p offsets and @p fieldnames points to preallocated memory.
 * You can call the function also with @c NULL pointers, in this case the pointer is not used.
 *
 */
int mess_status_st_offset(int* size, int* offsets, char** fieldnames){

    int counter = 0;

    if(size){ *size = sizeof(struct mess_status_st);}

    MESS_OFFSET(struct mess_status_st,res2_norms            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,rel_changes           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,res2_norm             ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,res2_change           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,rel_change            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,res2_0                ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,time_all              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,time_adi              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,it                    ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,n_internal_status     ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,internal_status       ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,stop_res2             ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,stop_res2c            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,stop_rel              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,stop_user             ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_status_st,unstable              ,counter,offsets,fieldnames);

    return counter;
}


/**
 * @brief Information about memory laout of @ref mess_options_st
 * @param[in,out] size          input/output contains @c sizeof  value
 * @param[in,out] offsets       input/output contains relative address of fields of @ref mess_options_st
 * @param[in,out] fieldnames    input/outout contains correspondig field names of @param offsets
 * @return number of fields of @ref mess_options_st
 *
 * The @ref mess_options_st_offset function returns informations about the memory layout of @ref mess_options_st.
 * We assume that @p offsets and @p fieldnames points to preallocated memory.
 * You can call the function also with @c NULL pointers, in this case the pointer is not used.
 *
 */
int mess_options_st_offset(int* size, int* offsets, char** fieldnames){

    int counter = 0;

    if(size){ *size = sizeof(struct mess_options_st);}

    MESS_OFFSET(struct mess_options_st,type                 ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_shifts_p         ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_shifts_paratype  ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_shifts_arp_p     ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_shifts_arp_m     ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_shifts_l0        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_shifts_b0        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,proj_space           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_maxit            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_res2_tol         ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_res2c_tol        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_rel_change_tol   ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,adi_output           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,stepfunction         ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,stepfunction_aux     ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,residual_method      ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_maxit             ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_res2_tol          ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_output            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_linesearch        ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_gpStep            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_stepfunction      ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,nm_stepfunction_aux  ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,memory_usage         ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,K0                   ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_options_st,W                    ,counter,offsets,fieldnames);

    return counter;
}


/**
 * @brief Information about memory laout of @ref mess_equation_st
 * @param[in,out] size          input/output contains @c sizeof  value
 * @param[in,out] offsets       input/output contains relative address of fields of @ref mess_equation_st
 * @param[in,out] fieldnames    input/outout contains correspondig field names of @param offsets
 * @return number of fields of @ref mess_equation_st
 *
 * The @ref mess_equation_st_offset function returns informations about the memory layout of @ref mess_equation_st.
 * We assume that @p offsets and @p fieldnames points to preallocated memory.
 * You can call the function also with @c NULL pointers, in this case the pointer is not used.
 *
 */
int mess_equation_st_offset(int* size, int* offsets, char** fieldnames){

    int counter = 0;

    if(size){ *size = sizeof(struct mess_equation_st);}

    MESS_OFFSET(struct mess_equation_st,child               ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,parent              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,dim                 ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,B                   ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,C                   ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,K                   ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,RHS                 ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,clearRHS            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,clearB              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,clearC              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,AX                  ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,EX                  ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,AINV                ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,EINV                ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,ApEX                ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,ApEINV              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,parameter           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,init_rhs            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,clear               ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_equation_st,aux                 ,counter,offsets,fieldnames);

    return counter;
}


/**
 * @brief Information about memory laout of @ref mess_direct_st
 * @param[in,out] size          input/output contains @c sizeof  value
 * @param[in,out] offsets       input/output contains relative address of fields of @ref mess_direct_st
 * @param[in,out] fieldnames    input/outout contains correspondig field names of @param offsets
 * @return number of fields of @ref mess_direct_st
 *
 * The @ref mess_direct_st_offset function returns informations about the memory layout of @ref mess_direct_st.
 * We assume that @p offsets and @p fieldnames points to preallocated memory.
 * You can call the function also with @c NULL pointers, in this case the pointer is not used.
 *
 */
int mess_direct_st_offset(int* size, int* offsets, char** fieldnames){

    int counter = 0;

    if(size){ *size = sizeof(struct mess_direct_st);}

    MESS_OFFSET(struct mess_direct_st,data_type         ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,rows              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,cols              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,data              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,solve             ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,solvet            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,solvem            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,solvemt           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,solveh            ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,solvemh           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,getL              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,getU              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,getpermp          ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,getpermq          ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,getscalerow       ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,getscalecol       ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,det               ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,detc              ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,inverse           ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,clear             ,counter,offsets,fieldnames);
    MESS_OFFSET(struct mess_direct_st,name              ,counter,offsets,fieldnames);

    return counter;
}

