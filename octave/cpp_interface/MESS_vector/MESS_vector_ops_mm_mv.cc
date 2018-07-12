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
#include "MESS_vector.h"
#include "MESS_error.h"



/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_xmul (const mess_operation_t opA, const MESS_matrix & m, const MESS_vector & v){
    mess_vector result=NULL;
    int ret = 0;
    if(m && v){
        ret = mess_vector_init(&result);                                        OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);
        ret = mess_matrix_mvp(opA, (mess_matrix)m, (mess_vector)v, result);     OCTMESS_ERROR(ret,(ret!=0),mess_matrix_mvp);
    }
    return new MESS_vector(result);
}


octave_value op_mul (const MESS_matrix & m, const MESS_vector & v){return op_xmul(MESS_OP_NONE,m,v);}
octave_value op_trans_mul (const MESS_matrix & m, const MESS_vector & v){return op_xmul(MESS_OP_TRANSPOSE,m,v);}
octave_value op_herm_mul (const MESS_matrix & m, const MESS_vector & v){return op_xmul(MESS_OP_HERMITIAN,m,v);}


octave_value op_xldiv (const MESS_matrix & m ,const mess_operation_t opA, const MESS_vector & v){
    mess_vector result=NULL;
    int ret = 0;
    if(m && v){
        ret = mess_vector_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);
        ret = mess_matrix_backslash(opA,(mess_matrix)m, (mess_vector)v, result);    OCTMESS_ERROR(ret,(ret!=0),mess_matrix_backslash);
    }
    return new MESS_vector(result);
}

octave_value op_ldiv (const MESS_matrix & m, const MESS_vector & v){return op_xldiv(m, MESS_OP_NONE,v);}
octave_value op_trans_ldiv (const MESS_matrix & m, const MESS_vector & v){return op_xldiv(m, MESS_OP_TRANSPOSE,v);}
octave_value op_herm_ldiv (const MESS_matrix & m, const MESS_vector & v){return op_xldiv(m, MESS_OP_HERMITIAN,v);}
