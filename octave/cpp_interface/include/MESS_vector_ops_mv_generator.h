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

#ifndef MESS_VECTOR_OPS_MV_GENERATOR_OCTAVE_H_
#define MESS_VECTOR_OPS_MV_GENERATOR_OCTAVE_H_




#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(OCTAVE_MATRIX_TYPE, OCTAVE_MATRIX_TYPE_VALUE)                                                                     \
                                                                                                                                            \
    octave_value op_xmul (const mess_operation_t opA, const OCTAVE_MATRIX_TYPE & m,  const MESS_vector & v){                                \
        int ret = 0;                                                                                                                        \
        mess_vector result=NULL;                                                                                                            \
        mess_matrix_st dummy = octmess_matrix_conv(m.OCTAVE_MATRIX_TYPE_VALUE());   /* dummy needs no clear */                              \
        ret = mess_vector_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);           \
        ret = mess_matrix_mvp(opA, &dummy, (mess_vector)v, result);                 OCTMESS_ERROR(ret,(ret!=0),mess_matrix_mvp);            \
        return new MESS_vector(result);                                                                                                     \
    }                                                                                                                                       \
                                                                                                                                            \
                                                                                                                                            \
    octave_value op_mul (const OCTAVE_MATRIX_TYPE & m, const MESS_vector & v){return op_xmul(MESS_OP_NONE, m, v);}                          \
    octave_value op_trans_mul (const OCTAVE_MATRIX_TYPE & m, const MESS_vector & v){return op_xmul(MESS_OP_TRANSPOSE, m, v);}               \
    octave_value op_herm_mul (const OCTAVE_MATRIX_TYPE & m, const MESS_vector & v){return op_xmul(MESS_OP_HERMITIAN, m, v);}                \
                                                                                                                                            \
                                                                                                                                            \
    octave_value op_xldiv (const OCTAVE_MATRIX_TYPE & m ,const mess_operation_t opA, const MESS_vector & v){                                \
        int ret = 0;                                                                                                                        \
        mess_vector result=NULL;                                                                                                            \
        mess_matrix_st dummy = octmess_matrix_conv(m.OCTAVE_MATRIX_TYPE_VALUE());   /* dummy needs no clear */                              \
        ret = mess_vector_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);           \
        ret = mess_matrix_backslash(opA,&dummy,(mess_vector)v,result);              OCTMESS_ERROR(ret,(ret!=0),mess_matrix_backslash);      \
        return new MESS_vector(result);                                                                                                     \
    }                                                                                                                                       \
                                                                                                                                            \
                                                                                                                                            \
    octave_value op_ldiv(const OCTAVE_MATRIX_TYPE & m, const MESS_vector & v){return op_xldiv(m,MESS_OP_NONE,v);}                           \
    octave_value op_trans_ldiv(const OCTAVE_MATRIX_TYPE & m, const MESS_vector & v){return op_xldiv(m,MESS_OP_TRANSPOSE,v);}                \
    octave_value op_herm_ldiv(const OCTAVE_MATRIX_TYPE & m, const MESS_vector & v){return op_xldiv(m,MESS_OP_HERMITIAN,v);}                 \

#endif
