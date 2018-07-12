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

#include <octave/oct.h>
#include <octave/ops.h>
#include <octave/ov-re-sparse.h>
#include "mess/mess.h"
#include "MESS_matrix.h"
#include "MESS_error.h"

/*-----------------------------------------------------------------------------
 *  unary operations
 *-----------------------------------------------------------------------------*/
octave_value op_not(const MESS_matrix & m){
    int ret = 0;
    mess_matrix nm = NULL;
    ret = mess_matrix_init(&nm);                    OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_copy((mess_matrix)m,nm);      OCTMESS_ERROR(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_map_not(nm,1);                OCTMESS_ERROR(ret, (ret!=0), mess_matrix_map_not);
    return new MESS_matrix(nm);
}

octave_value op_uplus (const MESS_matrix & m){
    return new MESS_matrix(m);
}


octave_value op_uminus (const MESS_matrix & m){
    int ret = 0;
    mess_matrix copy = NULL;
    ret = mess_matrix_init(&copy);                      OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_copy((mess_matrix)m, copy);       OCTMESS_ERROR(ret, (ret!=0), mess_matrix_copy);
    ret = mess_matrix_scale(-1.0, copy);                OCTMESS_ERROR(ret, (ret!=0), mess_matrix_scale);
    return new MESS_matrix(copy);
}


octave_value op_transpose(const MESS_matrix & m){
    int ret = 0;
    mess_matrix hm=NULL;
    ret = mess_matrix_init(&hm);                        OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_transpose((mess_matrix)m,hm);     OCTMESS_ERROR(ret, (ret!=0), mess_matrix_transpose);
    return new MESS_matrix(hm);
}


octave_value op_hermitian(const MESS_matrix & m){
    int ret = 0;
    mess_matrix hm=NULL;
    ret = mess_matrix_init(&hm);                        OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_ctranspose((mess_matrix)m,hm);    OCTMESS_ERROR(ret, (ret!=0), mess_matrix_ctranspose);
    return new MESS_matrix(hm);
}



/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
octave_value op_add (const MESS_matrix & m1, const MESS_matrix & m2) {
    mess_matrix result=NULL;
    int ret = 0;
    if(m1.get_ptr() && m2.get_ptr()){
        ret = mess_matrix_init(&result);                            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_copy((mess_matrix)m1,result);             OCTMESS_ERROR(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_add(1.0,(mess_matrix)m2,1.0,result);      OCTMESS_ERROR(ret,(ret!=0),mess_matrix_add);
    }
    return new MESS_matrix(result);
}


octave_value op_sub (const MESS_matrix & m1, const MESS_matrix & m2) {
    mess_matrix result=NULL;
    int ret = 0;
    if(m1.get_ptr() && m2.get_ptr()){
        ret = mess_matrix_init(&result);                            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_copy((mess_matrix)m1,result);             OCTMESS_ERROR(ret,(ret!=0),mess_matrix_copy);
        ret = mess_matrix_add(-1.0,(mess_matrix)m2,1.0,result);     OCTMESS_ERROR(ret,(ret!=0),mess_matrix_add);
    }
    return new MESS_matrix(result);
}

octave_value op_xmul (const mess_operation_t opA, const MESS_matrix & m1, const mess_operation_t opB, const MESS_matrix & m2){
    mess_matrix result=NULL;
    int ret = 0;
    if(m1.get_ptr() && m2.get_ptr()){
        ret = mess_matrix_init(&result);                                                    OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_multiply(opA,(mess_matrix)m1,opB,(mess_matrix)m2,result);         OCTMESS_ERROR(ret,(ret!=0),mess_matrix_multiply);
    }
    return new MESS_matrix(result);
}


octave_value op_mul (const MESS_matrix & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_NONE,m1,MESS_OP_NONE,m2);}
octave_value op_trans_mul (const MESS_matrix & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_TRANSPOSE,m1,MESS_OP_NONE,m2);}
octave_value op_herm_mul (const MESS_matrix & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_HERMITIAN,m1,MESS_OP_NONE,m2);}
octave_value op_mul_trans (const MESS_matrix & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_NONE,m1,MESS_OP_TRANSPOSE,m2);}
octave_value op_mul_herm (const MESS_matrix & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_NONE,m1,MESS_OP_HERMITIAN,m2);}


octave_value op_xldiv (const MESS_matrix & m1 ,const mess_operation_t opA, const MESS_matrix & m2){
    mess_matrix result=NULL;
    int ret = 0;
    if(m1.get_ptr() && m2.get_ptr()){
        ret = mess_matrix_init(&result);                                                    OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);
        ret = mess_matrix_backslashm(opA,(mess_matrix)m1,(mess_matrix)m2,result);           OCTMESS_ERROR(ret,(ret!=0),mess_matrix_backslashm);
    }
    return new MESS_matrix(result);
}

octave_value op_ldiv (const MESS_matrix & m1, const MESS_matrix & m2){return op_xldiv(m1,MESS_OP_NONE,m2);}
octave_value op_trans_ldiv (const MESS_matrix & m1, const MESS_matrix & m2){return op_xldiv(m1,MESS_OP_TRANSPOSE,m2);}
octave_value op_herm_ldiv (const MESS_matrix & m1, const MESS_matrix & m2){return op_xldiv(m1,MESS_OP_HERMITIAN,m2);}

