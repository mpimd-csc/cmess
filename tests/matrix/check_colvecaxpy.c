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
 * @file tests/matrix/check_colvecaxpy.c
 * @brief Check the axpy scalar vector product.
 * @author  @mbehr
 * @test
 * This function check the @ref mess_matrix_colvecaxpy function defined in colops.c that means it checks if the axpy scalar
 * vector product
 * \f[v+coeff \cdot Q(:,col) \f]
 * is computed correctly for a given matrix \f$ Q \f$, a vector \f$ v \f$, a column \f$ col \f$ and a coefficient
 * \f$ coeff \f$.
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


#define COMPUTESOL(A,COL,COEFF,V,VSOL,VTEMP){                           \
    CALL(mess_vector_copy(V,VSOL));                                     \
    CALL(mess_matrix_getcol(A,COL,VTEMP));                              \
    CALL(mess_vector_tocomplex(VTEMP));                                 \
    CALL(mess_vector_tocomplex(VSOL));                                  \
    CALL(mess_vector_axpyc(COEFF,VTEMP,VSOL));                          \
}

#define CHECKCOLVCAXPY(A,COL,COEFF,VCHECK,V,VSOL,NRM,ERR,EPS){  \
    CALL(mess_vector_copy(V,VCHECK))                                    \
    CALL(mess_matrix_colvecaxpy(COEFF,COL,A,VCHECK));                   \
    /*mess_vector_print(VCHECK);*/                                      \
    /*mess_vector_print(VSOL);*/                                        \
    CALL(mess_vector_diffnorm(VCHECK,VSOL,&NRM));                       \
    printf("%e\n",NRM);                                                 \
    if(NRM>10*EPS){++ERR;}                                              \
}

int main ( int argc, char ** argv) {
    mess_init();
    int n=4;
    double p = 0.6; //density for random matrix generation
    mess_int_t col=2;
    mess_double_cpx_t coeffc=2+3*I;
    double coeffr=2;
    double nrm=0;
    int ret;
    int err = 0;
    double eps = mess_eps();
    mess_matrix Ac,Acsrc, Acscc, Ar,Acsrr,Acscr;
    mess_vector vr,vc,vcheck,vsol, vtemp;

    //init real matrices
    CALL(mess_matrix_init(&Ar));
    CALL(mess_matrix_init(&Acscr));
    CALL(mess_matrix_init(&Acsrr));

    ///init complex matrices
    CALL(mess_matrix_init(&Ac));
    CALL(mess_matrix_init(&Acscc));
    CALL(mess_matrix_init(&Acsrc));

    /*-----------------------------------------------------------------------------
     *  Load matrices
     *-----------------------------------------------------------------------------*/
    //load real matrices
    CALL(mess_matrix_rand(Ar,n,n,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_convert(Ar,Acsrr,MESS_CSR));
    CALL(mess_matrix_convert(Ar,Acscr,MESS_CSC));
    CALL(mess_vector_init(&vr));
    CALL(mess_vector_alloc(vr,Ar->rows,MESS_REAL));
    CALL(mess_vector_rand(vr));
    CALL(mess_vector_init(&vcheck));
    CALL(mess_vector_alloc(vcheck,Ar->rows,MESS_REAL));

    //load complex matrices
    CALL(mess_matrix_rand(Ac,n,n,MESS_DENSE,MESS_COMPLEX,p));
    CALL(mess_matrix_convert(Ac,Acsrc,MESS_CSR));
    CALL(mess_matrix_convert(Ac,Acscc,MESS_CSC));
    CALL(mess_vector_init(&vc));
    CALL(mess_vector_alloc(vc,Ac->rows,MESS_COMPLEX));
    CALL(mess_vector_rand(vc));


    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&vtemp,&vsol);
    CALL(mess_vector_alloc(vtemp,n,MESS_REAL));
    CALL(mess_vector_alloc(vsol,n,MESS_REAL));



    //test function for case real Q, real v, real coeff
    COMPUTESOL(Ar, col, coeffr, vr, vsol, vtemp);
    CHECKCOLVCAXPY(Ar, col, coeffc, vcheck, vr, vsol, nrm, err, eps);           //imaginary part of coeffc is ignored
    CHECKCOLVCAXPY(Acsrr, col, coeffc, vcheck, vr, vsol, nrm, err, eps);        //imaginary part of coeffc is ignored
    CHECKCOLVCAXPY(Acscr, col, coeffc, vcheck, vr, vsol, nrm, err, eps);        //imaginary part of coeffc is ignored

    //test function for case  real Q, complex v
    COMPUTESOL(Ar, col, coeffc, vc, vsol, vtemp);
    CHECKCOLVCAXPY(Ar, col, coeffc, vcheck, vc, vsol, nrm, err, eps);
    CHECKCOLVCAXPY(Acsrr, col, coeffc, vcheck, vc, vsol, nrm, err, eps);
    CHECKCOLVCAXPY(Acscr, col, coeffc, vcheck, vc, vsol, nrm, err, eps);

    //test function for case  complex Q, complex v
    COMPUTESOL(Ac, col, coeffc, vc, vsol, vtemp);
    CHECKCOLVCAXPY(Ac, col, coeffc, vcheck, vc, vsol, nrm, err, eps);
    CHECKCOLVCAXPY(Acsrc, col, coeffc, vcheck, vc, vsol, nrm, err, eps);
    CHECKCOLVCAXPY(Acscc, col, coeffc, vcheck, vc, vsol, nrm, err, eps);

    //test function for case  complex Q, real v
    COMPUTESOL(Ac, col, coeffc, vr, vsol, vtemp);
    CHECKCOLVCAXPY(Ac, col, coeffc, vcheck, vr, vsol, nrm, err, eps);
    CHECKCOLVCAXPY(Acsrc, col, coeffc, vcheck, vr, vsol, nrm, err, eps);
    CHECKCOLVCAXPY(Acscc, col, coeffc, vcheck, vr, vsol, nrm, err, eps);


    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/


    mess_matrix_clear(&Ar);
    mess_matrix_clear(&Acsrr);
    mess_matrix_clear(&Acscr);
    mess_matrix_clear(&Ac);
    mess_matrix_clear(&Acsrc);
    mess_matrix_clear(&Acscc);
    mess_vector_clear(&vr);
    mess_vector_clear(&vc);
    mess_vector_clear(&vtemp);
    mess_vector_clear(&vsol);
    mess_vector_clear(&vcheck);

    return (err>0)?(1):(0);
}
