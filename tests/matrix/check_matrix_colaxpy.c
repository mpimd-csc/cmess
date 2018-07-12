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
 * @file tests/matrix/check_matrix_colaxpy.c
 * @brief Check the column update of a matrix using an axpy-operation.
 * @author  @mbehr
 * @test
 * This function checks the @ref mess_matrix_colaxpy fucntion defined in colops.c that means it checks if an update of a
 * column of matrix using an axpy-operation
 * \f[Q(:,col)=coeff v + Q(:,col) \f]
 * is computed correctly for a given matrix \f$ Q \f$, a given vector \f$ v \f$, a given column \f$ col \f$ and a given
 * coefficient \f$ coeff \f$.
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

#define CHECKMATRIXOCLAXPY(Q,TEMP,COL,COEFF,TYPE,V1,VTEMP1,VTEMP2,VSOL,ERR,EPS,DIFF){           \
    CALL(mess_matrix_copy(Q,TEMP));                                                         \
    CALL(mess_vector_copy(V1,VTEMP1));                                                      \
    CALL(mess_matrix_getcol(TEMP,COL,VSOL));                                                \
    if(MESS_REAL==TYPE){                                                                    \
        CALL(mess_vector_axpy(COEFF,VTEMP1,VSOL));                                              \
    }else{                                                                                  \
        CALL(mess_vector_tocomplex(VTEMP1));                                                    \
        CALL(mess_vector_tocomplex(VSOL));                                                      \
        CALL(mess_vector_axpyc(COEFF,VTEMP1,VSOL));                                             \
    }                                                                                       \
    CALL(mess_matrix_colaxpy(COEFF,VTEMP1,COL,TEMP));                                       \
    CALL(mess_matrix_getcol(TEMP,COL,VTEMP2));                                              \
    CALL(mess_vector_diffnorm(VTEMP2,VSOL,&DIFF));                                          \
    if(DIFF>EPS){                                                                       \
        printf("Failed diff:%e\n",DIFF);                                                \
        ++ERR;                                                                          \
    }                                                                                   \
}

int main ( int argc, char ** argv) {
    mess_init();
    int dim=10, ret=0, err=0;
    double p=1 ;                    //density for random matrix generation, in dense case p is neglected
    double eps = 1e-12;
    double diff=0;
    double realcoeff = rand();
    mess_double_cpx_t cpxcoeff = (1+I)*realcoeff;

    mess_matrix Qr, Qc,Temp;
    mess_vector v1,vsol,vtemp1,vtemp2;

    int col1=0;

    /*-----------------------------------------------------------------------------
     *  Init matrices and vector
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&Qr));
    CALL(mess_matrix_init(&Qc));
    CALL(mess_matrix_init(&Temp));

    MESS_INIT_VECTORS(&v1,&vsol,&vtemp1,&vtemp2);
    mess_vector_alloc(v1,dim,MESS_REAL);
    mess_vector_alloc(vsol,dim,MESS_REAL);
    mess_vector_alloc(vtemp1,dim,MESS_REAL);
    mess_vector_alloc(vtemp2,dim,MESS_REAL);

    /*-----------------------------------------------------------------------------
     *  Load matrices and vector
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_rand(Qr,dim,dim,MESS_DENSE,MESS_REAL,p));
    //  CALL(mess_vector_rand(v1));
    CALL(mess_vector_ones(v1));
    CALL(mess_matrix_rand(Qc,dim,dim,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_tocomplex(Qc));
    CALL(mess_matrix_scalec(1+I,Qc));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    //real matrix, real coeff, real v1
    CHECKMATRIXOCLAXPY(Qr,Temp,col1,realcoeff,MESS_REAL,v1,vtemp1,vtemp2,vsol,err,eps,diff)
        //real matrix, complex coeff, real v1
        CHECKMATRIXOCLAXPY(Qr,Temp,col1,cpxcoeff,MESS_COMPLEX,v1,vtemp1,vtemp2,vsol,err,eps,diff)
        //complex matrix, real coeff, real v1
        CHECKMATRIXOCLAXPY(Qc,Temp,col1,realcoeff,MESS_REAL,v1,vtemp1,vtemp2,vsol,err,eps,diff)
        //complex matrix, complex coeff, real v1
        CHECKMATRIXOCLAXPY(Qc,Temp,col1,cpxcoeff,MESS_COMPLEX,v1,vtemp1,vtemp2,vsol,err,eps,diff)


        mess_vector_tocomplex(v1);

    //real matrix, real coeff, complex v1
    CHECKMATRIXOCLAXPY(Qr,Temp,col1,realcoeff,MESS_REAL,v1,vtemp1,vtemp2,vsol,err,eps,diff)
        //real matrix, complex coeff, complex v1
        CHECKMATRIXOCLAXPY(Qr,Temp,col1,cpxcoeff,MESS_COMPLEX,v1,vtemp1,vtemp2,vsol,err,eps,diff)
        //complex matrix, real coeff, complex v1
        CHECKMATRIXOCLAXPY(Qc,Temp,col1,realcoeff,MESS_REAL,v1,vtemp1,vtemp2,vsol,err,eps,diff)
        //complex matrix, complex coeff, complex v1
        CHECKMATRIXOCLAXPY(Qc,Temp,col1,cpxcoeff,MESS_COMPLEX,v1,vtemp1,vtemp2,vsol,err,eps,diff)


        /*-----------------------------------------------------------------------------
         *  Clear Memory
         *-----------------------------------------------------------------------------*/
        mess_matrix_clear(&Qr);
    mess_matrix_clear(&Qc);
    mess_matrix_clear(&Temp);
    mess_vector_clear(&v1);
    mess_vector_clear(&vtemp1);
    mess_vector_clear(&vtemp2);
    mess_vector_clear(&vsol);

    return (err>0)?(1):(0);
}
