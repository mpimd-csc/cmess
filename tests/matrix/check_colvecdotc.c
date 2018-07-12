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
 * @file tests/matrix/check_colvecdotc.c
 * @brief Check the complex column vector product.
 * @author @koehlerm
 * @test
 * This function check the @ref mess_matrix_colvecdotc function defined in colops.c that means it checks if the
 * vector product
 * \f[Q(:,col)^Hv \f]
 * is computed correctly for a given complex matrix \f$ Q \f$, a vector \f$ v \f$ and a column \f$ col \f$.
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


mess_double_cpx_t matrix_datac[16] = {11,0,31,0,12,22,0,42,0,0,33,43+40*I,14,24,0,0};
double matrix_datar[16] = {11,0,31,0,12,22,0,42,0,0,33,43,14,24,0,0};


int main ( int argc, char ** argv) {
    mess_init();
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix Ac,Acsrc, Acscc, Ar,Acsrr,Acscr;
    mess_vector vr,vc;
    mess_double_cpx_t  dot1,dot2,dot3;

    //init real matrices
    CALL(mess_matrix_init(&Ar));
    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acsrr));
    //init complex matrices
    CALL(mess_matrix_init(&Ac));
    CALL(mess_matrix_init(&Acscc));
    CALL(mess_matrix_init(&Acsrc));


    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load real matrices
    CALL(mess_matrix_dense_from_farray(Ar,4,4,0,matrix_datar,NULL));
    CALL(mess_matrix_convert(Ar,Acsrr,MESS_CSR));
    CALL(mess_matrix_convert(Ar,Acscr,MESS_CSC));
    MESS_INIT_VECTORS(&vr);
    CALL(mess_vector_alloc(vr,Ar->rows,MESS_REAL));
    CALL(mess_vector_ones(vr));

    //load complex matrices
    CALL(mess_matrix_dense_from_farray(Ac,4,4,0,NULL,matrix_datac));
    CALL(mess_matrix_convert(Ac,Acsrc,MESS_CSR));
    CALL(mess_matrix_convert(Ac,Acscc,MESS_CSC));
    MESS_INIT_VECTORS(&vc);
    CALL(mess_vector_alloc(vc,Ac->rows,MESS_COMPLEX));
    CALL(mess_vector_ones(vc));
    CALL(mess_vector_scalec(I,vc));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    //complex Matrix complex vector
    CALL(mess_matrix_colvecdotc(Ac,3,vc,&dot1));
    CALL(mess_matrix_colvecdotc(Acsrc,3,vc,&dot2));
    CALL(mess_matrix_colvecdotc(Acscc,3,vc,&dot3));

    if ( cabs(dot1-dot2)/cabs(dot1) > 10 * eps) err++;
    if ( cabs(dot1-dot3)/cabs(dot1) > 10 * eps) err++;

    printf("err = %d \n",err);
    printf("%lg + %lg \n",creal(dot1), cimag(dot1));
    printf("%lg + %lg \n",creal(dot2), cimag(dot2));
    printf("%lg + %lg \n",creal(dot3), cimag(dot3));

    //complex Matrix real vector
    CALL(mess_matrix_colvecdotc(Ac,3,vr,&dot1));
    CALL(mess_matrix_colvecdotc(Acsrc,3,vr,&dot2));
    CALL(mess_matrix_colvecdotc(Acscc,3,vr,&dot3));

    if ( cabs(dot1-dot2)/cabs(dot1) > 10 * eps) err++;
    if ( cabs(dot1-dot3)/cabs(dot1) > 10 * eps) err++;

    printf("err = %d \n",err);
    printf("%lg + %lg \n",creal(dot1), cimag(dot1));
    printf("%lg + %lg \n",creal(dot2), cimag(dot2));
    printf("%lg + %lg \n",creal(dot3), cimag(dot3));

    //real Matrix complex vector
    CALL(mess_matrix_colvecdotc(Ar,3,vc,&dot1));
    CALL(mess_matrix_colvecdotc(Acsrr,3,vc,&dot2));
    CALL(mess_matrix_colvecdotc(Acscr,3,vc,&dot3));

    if ( cabs(dot1-dot2)/cabs(dot1) > 10 * eps) err++;
    if ( cabs(dot1-dot3)/cabs(dot1) > 10 * eps) err++;

    printf("err = %d \n",err);
    printf("%lg + %lg \n",creal(dot1), cimag(dot1));
    printf("%lg + %lg \n",creal(dot2), cimag(dot2));
    printf("%lg + %lg \n",creal(dot3), cimag(dot3));

    //case real matrix real vector should be tested in check_colvecdot

    mess_matrix_clear(&Ar);
    mess_matrix_clear(&Acsrr);
    mess_matrix_clear(&Acscr);
    mess_vector_clear(&vr);
    mess_matrix_clear(&Ac);
    mess_matrix_clear(&Acsrc);
    mess_matrix_clear(&Acscc);
    mess_vector_clear(&vc);
    return (err>0)?(1):(0);
}

