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
 * @addtogroup test_matrix
 * @{
 * @file tests/matrix/check_fro_inn.c
 * @brief Check the frobenius inner product.
 * @author @dykstra
 * @test
 * This function checks the @ref mess_matrix_fro_inn function defined in trace.c, that means it checks if the frobenius inner product
 * \f[ fro = \langle op(A),B \rangle_F = \sum\limits_{i,j} \overline{op(A_{i,j})} B_{i,j} \f]
 * is computed correctly for given matrices \f$ A \f$ and \f$ B \f$.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"

double matrix_data1[16] = {11,0,31,0,12,22,0,42,0,0,33,43,14,24,0,0};
double matrix_data2[16] = {-12,0,1,-7,-5,3,42,2,1,23,2,1,-4,0,1,1337};


int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix A,B,Ac,Bc;
    double fro1,fro2,fro3;
    mess_double_cpx_t fro4,fro5,fro6, im = I;

    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&B));
    CALL(mess_matrix_init(&Ac));
    CALL(mess_matrix_init(&Bc));

    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_dense_from_farray(A,4,4,0,matrix_data1,NULL));
    CALL(mess_matrix_dense_from_farray(B,4,4,0,matrix_data2,NULL));

    CALL(mess_matrix_fro_inn(MESS_OP_NONE,A,B,(void*)&fro1));//sollte 42 sein
    CALL(mess_matrix_fro_inn(MESS_OP_TRANSPOSE,A,B,(void*)&fro2));//sollte 24 sein
    CALL(mess_matrix_fro_inn(MESS_OP_HERMITIAN,A,B,(void*)&fro3));//sollte 24 sein

    if ( fabs(fro1-42)/42 > 10 * eps) err++;
    if ( fabs(fro2-24)/24 > 10 * eps) err++;
    if ( fabs(fro3-24)/24 > 10 * eps) err++;

    printf("fro_inn(A,B):err = %d \nfroNONE=%lg\nfroTRANS=%lg\nfroHERM=%lg\ncorrect result: (42,24,24)\n",err,fro1,fro2,fro3);

    CALL(mess_matrix_copy(A,Ac));
    CALL(mess_matrix_copy(B,Bc));
    mess_matrix_tocomplex(Ac);
    mess_matrix_tocomplex(Bc);

    CALL(mess_matrix_fro_inn(MESS_OP_NONE,Ac,Bc,(void*)&fro4));//sollte 42 sein
    CALL(mess_matrix_fro_inn(MESS_OP_TRANSPOSE,Ac,Bc,(void*)&fro5));//sollte 24 sein
    CALL(mess_matrix_fro_inn(MESS_OP_HERMITIAN,Ac,Bc,(void*)&fro6));//sollte 24 sein

    if ( cabs(fro4-42)/42 > 10 * eps) err++;
    if ( cabs(fro5-24)/24 > 10 * eps) err++;
    if ( cabs(fro6-24)/24 > 10 * eps) err++;

    printf("fro_inn(Ac,Bc):err = %d \nfroNONE=%lg+(%lg*i)\nfroTRANS=%lg+(%lg*i)\nfroHERM=%lg+(%lg*i)\ncorrect result: (42,24,24)\n",err,creal(fro4),cimag(fro4),creal(fro5),cimag(fro5),creal(fro6),cimag(fro6));
    
    CALL(mess_matrix_scalec(im,Ac));
    CALL(mess_matrix_scalec(im,Bc));

    printf("Ac, Bc scaled with i:\n");
    
    CALL(mess_matrix_fro_inn(MESS_OP_NONE,Ac,Bc,(void*)&fro4));//sollte 42 sein
    CALL(mess_matrix_fro_inn(MESS_OP_TRANSPOSE,Ac,Bc,(void*)&fro5));//sollte 24 sein
    CALL(mess_matrix_fro_inn(MESS_OP_HERMITIAN,Ac,Bc,(void*)&fro6));//sollte -24 sein

    if ( cabs(fro4-42)/42 > 10 * eps) err++;
    if ( cabs(fro5-24)/24 > 10 * eps) err++;
    if ( cabs(fro6+24)/24 > 10 * eps) err++;

    printf("fro_inn(Ac,Bc):err = %d \nfroNONE=%lg+(%lg*i)\nfroTRANS=%lg+(%lg*i)\nfroHERM=%lg+(%lg*i)\ncorrect result: (42,24,-24)\n",err,creal(fro4),cimag(fro4),creal(fro5),cimag(fro5),creal(fro6),cimag(fro6));
    
    CALL(mess_matrix_fro_inn(MESS_OP_NONE,Ac,B,(void*)&fro4));//sollte -42i sein
    CALL(mess_matrix_fro_inn(MESS_OP_TRANSPOSE,Ac,B,(void*)&fro5));//sollte -24i sein
    CALL(mess_matrix_fro_inn(MESS_OP_HERMITIAN,Ac,B,(void*)&fro6));//sollte 24i sein

    if ( cabs(fro4+(42*im))/42 > 10 * eps) err++;
    if ( cabs(fro5+(24*im))/24 > 10 * eps) err++;
    if ( cabs(fro6-(24*im))/24 > 10 * eps) err++;

    printf("fro_inn(Ac,B):err = %d \nfroNONE=%lg+(%lg*i)\nfroTRANS=%lg+(%lg*i)\nfroHERM=%lg+(%lg*i)\ncorrect result: (-42i,-24i,24i)\n",err,creal(fro4),cimag(fro4),creal(fro5),cimag(fro5),creal(fro6),cimag(fro6));
    
    CALL(mess_matrix_fro_inn(MESS_OP_NONE,A,Bc,(void*)&fro4));//sollte 42i sein
    CALL(mess_matrix_fro_inn(MESS_OP_TRANSPOSE,A,Bc,(void*)&fro5));//sollte 24i sein
    CALL(mess_matrix_fro_inn(MESS_OP_HERMITIAN,A,Bc,(void*)&fro6));//sollte 24i sein

    if ( cabs(fro4-(42*im))/42 > 10 * eps) err++;
    if ( cabs(fro5-(24*im))/24 > 10 * eps) err++;
    if ( cabs(fro6-(24*im))/24 > 10 * eps) err++;

    printf("fro_inn(A,Bc):err = %d \nfroNONE=%lg+(%lg*i)\nfroTRANS=%lg+(%lg*i)\nfroHERM=%lg+(%lg*i)\ncorrect result: (42i,24i,24i)\n",err,creal(fro4),cimag(fro4),creal(fro5),cimag(fro5),creal(fro6),cimag(fro6));





    mess_matrix_clear(&A);
    mess_matrix_clear(&B);
    mess_matrix_clear(&Ac);
    mess_matrix_clear(&Bc);
    return (err>0)?(1):(0);
}

