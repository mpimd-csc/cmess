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

#ifndef MESS_MATRIX_OPS_MM_GENERATOR_OCTAVE_H_
#define MESS_MATRIX_OPS_MM_GENERATOR_OCTAVE_H_




#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(OCTAVE_MATRIX_TYPE, OCTAVE_MATRIX_TYPE_VALUE)                                                                             \
    octave_value op_add (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2) {                                                                   \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv((m2.OCTAVE_MATRIX_TYPE_VALUE ()));    /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);                             \
        ret = mess_matrix_copy((mess_matrix)m1,result);                 OCTMESS_ERROR(ret, (ret!=0), mess_matrix_copy);                             \
        ret = mess_matrix_add(1.0,&dummy,1.0,result);                   OCTMESS_ERROR(ret, (ret!=0), mess_matrix_add);                              \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
                                                                                                                                                    \
    octave_value op_add(const OCTAVE_MATRIX_TYPE &m1, const MESS_matrix & m2){                                                                      \
        return op_add(m2,m1);                                                                                                                       \
    }                                                                                                                                               \
                                                                                                                                                    \
                                                                                                                                                    \
    octave_value op_sub (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2) {                                                                   \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv((m2.OCTAVE_MATRIX_TYPE_VALUE ()));    /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);                             \
        ret = mess_matrix_copy((mess_matrix)m1,result);                 OCTMESS_ERROR(ret, (ret!=0), mess_matrix_copy);                             \
        ret = mess_matrix_add(-1.0,&dummy,1.0,result);                  OCTMESS_ERROR(ret, (ret!=0), mess_matrix_add);                              \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
    octave_value op_sub (const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2) {                                                                   \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv((m1.OCTAVE_MATRIX_TYPE_VALUE ()));    /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                OCTMESS_ERROR(ret, (ret!=0), mess_matrix_init);                             \
        ret = mess_matrix_copy((mess_matrix)m2,result);                 OCTMESS_ERROR(ret, (ret!=0), mess_matrix_copy);                             \
        ret = mess_matrix_add(1.0,&dummy,-1.0,result);                  OCTMESS_ERROR(ret, (ret!=0), mess_matrix_add);                              \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
    octave_value op_xmul (const mess_operation_t opA, const MESS_matrix & m1, const mess_operation_t opB, const OCTAVE_MATRIX_TYPE & m2){           \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv(m2.OCTAVE_MATRIX_TYPE_VALUE ());      /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);                   \
        ret = mess_matrix_multiply(opA,(mess_matrix)m1,opB,&dummy,result);          OCTMESS_ERROR(ret,(ret!=0),mess_matrix_multiply);               \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
    octave_value op_xmul (const mess_operation_t opA, const OCTAVE_MATRIX_TYPE & m1, const mess_operation_t opB, const MESS_matrix & m2){           \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv(m1.OCTAVE_MATRIX_TYPE_VALUE ());      /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);                   \
        ret = mess_matrix_multiply(opA,&dummy,opB,(mess_matrix)m2,result);          OCTMESS_ERROR(ret,(ret!=0),mess_matrix_multiply);               \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
                                                                                                                                                    \
    octave_value op_mul (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xmul(MESS_OP_NONE, m1, MESS_OP_NONE, m2);}                \
    octave_value op_mul (const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_NONE, m1, MESS_OP_NONE, m2);}                \
    octave_value op_trans_mul (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xmul(MESS_OP_TRANSPOSE, m1, MESS_OP_NONE, m2);}     \
    octave_value op_trans_mul (const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_TRANSPOSE, m1, MESS_OP_NONE, m2);}     \
    octave_value op_herm_mul (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xmul(MESS_OP_HERMITIAN, m1, MESS_OP_NONE, m2);}      \
    octave_value op_herm_mul (const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_HERMITIAN, m1, MESS_OP_NONE, m2);}      \
    octave_value op_mul_trans (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xmul(MESS_OP_NONE, m1, MESS_OP_TRANSPOSE, m2);}     \
    octave_value op_mul_trans (const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_NONE, m1, MESS_OP_TRANSPOSE, m2);}     \
    octave_value op_mul_herm (const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xmul(MESS_OP_NONE, m1, MESS_OP_HERMITIAN, m2);}      \
    octave_value op_mul_herm (const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xmul(MESS_OP_NONE, m1, MESS_OP_HERMITIAN, m2);}      \
                                                                                                                                                    \
                                                                                                                                                    \
                                                                                                                                                    \
                                                                                                                                                    \
    octave_value op_xldiv (const MESS_matrix & m1 ,const mess_operation_t opA, const OCTAVE_MATRIX_TYPE & m2){                                      \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv(m2.OCTAVE_MATRIX_TYPE_VALUE ());      /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);                   \
        ret = mess_matrix_backslashm(opA,(mess_matrix)m1,&dummy,result);            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_backslashm);             \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
    octave_value op_xldiv (const OCTAVE_MATRIX_TYPE & m1 ,const mess_operation_t opA, const MESS_matrix & m2){                                      \
        int ret = 0;                                                                                                                                \
        mess_matrix result=NULL;                                                                                                                    \
        mess_matrix_st dummy = octmess_matrix_conv(m1.OCTAVE_MATRIX_TYPE_VALUE ());      /* dummy needs no clear */                                 \
        ret = mess_matrix_init(&result);                                            OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);                   \
        ret = mess_matrix_backslashm(opA,&dummy,(mess_matrix)m2, result);           OCTMESS_ERROR(ret,(ret!=0),mess_matrix_backslashm);             \
        return new MESS_matrix(result);                                                                                                             \
    }                                                                                                                                               \
                                                                                                                                                    \
                                                                                                                                                    \
                                                                                                                                                    \
    octave_value op_ldiv(const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xldiv(m1,MESS_OP_NONE,m2);}                               \
    octave_value op_ldiv(const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xldiv(m1,MESS_OP_NONE,m2);}                               \
    octave_value op_trans_ldiv(const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xldiv(m1,MESS_OP_TRANSPOSE,m2);}                    \
    octave_value op_trans_ldiv(const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xldiv(m1,MESS_OP_TRANSPOSE,m2);}                    \
    octave_value op_herm_ldiv(const MESS_matrix & m1, const OCTAVE_MATRIX_TYPE & m2){return op_xldiv(m1,MESS_OP_HERMITIAN,m2);}                     \
    octave_value op_herm_ldiv(const OCTAVE_MATRIX_TYPE & m1, const MESS_matrix & m2){return op_xldiv(m1,MESS_OP_HERMITIAN,m2);}                     \


#endif
