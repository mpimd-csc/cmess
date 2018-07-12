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
 * @file tests/matrix/check_set_get_col_row.c
 * @brief Check the overwriting of a matrix column.
 * @author  @mbehr
 * @test
 * Checks the @ref mess_matrix_setcol,  @ref mess_matrix_getcol, @ref mess_matrix_setrow, @ref mess_matrix_getrow functions.
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

//some functions to create vectors with difficult sparsity patterns
static double v_map_real(double val){return fabs(val)<0.75?val:0;};
static mess_double_cpx_t v_map_cpx(mess_double_cpx_t val){return cabs(val)<0.75?val:0;};



int main ( int argc, char ** argv) {
    mess_init();

    int ret=0, err=0;
    double diff=.0, eps=mess_eps();
    mess_matrix A, Acopy;
    mess_vector vset,vget,vtmp;

    mess_int_t idtA,idtv;
    mess_datatype_t dts [] = {MESS_REAL, MESS_COMPLEX};

    mess_int_t istA;
    mess_storage_t sts [] = {MESS_DENSE, MESS_CSR, MESS_CSC};

    mess_int_t p;
    double ps [] = {0.01, 0.2, 0.5, 0.8, 1.0};

    mess_int_t rows, rowsmax=600, cols, colsmax=600, colidx, rowidx;


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if(argc!=1){
        printf("Usage: %s\n",argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  init matrices and vectors
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_init(&Acopy));
    MESS_INIT_VECTORS(&vset,&vget,&vtmp);
    CALL(mess_vector_alloc(vset,0,MESS_REAL));
    CALL(mess_vector_alloc(vget,0,MESS_REAL));
    CALL(mess_vector_alloc(vtmp,0,MESS_REAL));


    /*-----------------------------------------------------------------------------
     *  iterate over all datatypes, storage types and test
     *-----------------------------------------------------------------------------*/
    //set random seed
    mess_int_t seed = 1;
    CALL(mess_matrix_rand_init(&seed));
    CALL(mess_vector_rand_init(&seed));


    for(p=0;p<5;++p){
        for(rows=1;rows<rowsmax;rows=rows*2){
            for(cols=1;cols<colsmax;cols=cols*3){
                for(idtA=0;idtA<2;idtA++){
                    for(idtv=0;idtv<2;idtv++){
                        for(istA=0;istA<3;istA++){

                            //create random matrix
                            CALL(mess_matrix_rand(A,rows,cols,sts[istA],dts[idtA],ps[p]));

                            /*-----------------------------------------------------------------------------
                             * setcol and getcol test: set a col and get a col and compare the vectors
                             *-----------------------------------------------------------------------------*/
                            //create a random vector with some "sparsity structure"
                            CALL(mess_vector_resize(vset,rows));
                            CALL(mess_vector_totype(vset,dts[idtv]));
                            CALL(mess_vector_rand(vset));
                            CALL(mess_vector_map(vset, v_map_real, v_map_cpx));
                            for(colidx=0;colidx<cols;colidx=1+colidx*3){

                                CALL(mess_matrix_setcol(A,colidx,vset));
                                CALL(mess_matrix_getcol(A,colidx,vget));

                                //set a complex column to a real matrix, the imaginary part is lost
                                if(MESS_IS_REAL(A) && MESS_IS_COMPLEX(vset)){
                                    CALL(mess_vector_realpart(vset,vtmp));
                                }else{
                                    CALL(mess_vector_copy(vset,vtmp));
                                }

                                //compute difference
                                CALL(mess_vector_diffnorminf(vtmp,vget,&diff));

                                if(diff>eps){
                                    //test failed print some informations
                                    printf("FAILED setcol/getcol failed with diff=%e\n",diff);
                                    printf("Matrix A:\n");
                                    mess_matrix_printinfo(A);
                                    printf("Vector v:\n");
                                    mess_vector_printinfo(vset);
                                    printf("Tried to set column ="MESS_PRINTF_INT"\n",colidx);

                                    err=1;
                                    goto clear;
                                }

                            }

                            /*-----------------------------------------------------------------------------
                             * setcol and getcol test: get a col and set a col and compare the matrices
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_matrix_copy(A,Acopy));
                            for(colidx=0;colidx<cols;colidx=1+colidx*3){

                                CALL(mess_matrix_getcol(A,colidx,vget));
                                CALL(mess_matrix_setcol(A,colidx,vget));

                                //compute difference
                                CALL(mess_matrix_diffnormf(A,Acopy,&diff));

                                if(diff>eps){
                                    //test failed print some informations
                                    printf("FAILED setcol/getcol failed with diff=%e\n",diff);
                                    printf("Matrix A:\n");
                                    mess_matrix_printinfo(A);
                                    mess_matrix_print(A);
                                    printf("Matrix Acopy:\n");
                                    mess_matrix_printinfo(Acopy);
                                    mess_matrix_print(Acopy);
                                    printf("Tried to set column ="MESS_PRINTF_INT"\n",colidx);

                                    err=1;
                                    goto clear;
                                }

                            }




                            /*-----------------------------------------------------------------------------
                             * setrow and getrow test: set a row and get a row and compare the vectors
                             *-----------------------------------------------------------------------------*/
                            //create a random vector with some "sparsity structure"
                            CALL(mess_vector_resize(vset,cols));
                            CALL(mess_vector_totype(vset,dts[idtv]));
                            CALL(mess_vector_rand(vset));
                            CALL(mess_vector_map(vset, v_map_real, v_map_cpx));
                            for(rowidx=0;rowidx<rows;rowidx=1+rowidx*2){

                                CALL(mess_matrix_setrow(A,rowidx,vset));
                                CALL(mess_matrix_getrow(A,rowidx,vget));

                                //set a complex column to a real matrix, the imaginary part is lost
                                if(MESS_IS_REAL(A) && MESS_IS_COMPLEX(vset)){
                                    CALL(mess_vector_realpart(vset,vtmp));
                                }else{
                                    CALL(mess_vector_copy(vset,vtmp));
                                }

                                //compute difference
                                CALL(mess_vector_diffnorminf(vtmp,vget,&diff));

                                if(diff>eps){
                                    //test failed print some informations
                                    printf("FAILED setrow/getrow failed with diff=%e\n",diff);
                                    printf("Matrix A:\n");
                                    mess_matrix_printinfo(A);
                                    printf("Vector v:\n");
                                    mess_vector_printinfo(vset);
                                    printf("Tried to set row ="MESS_PRINTF_INT"\n",rowidx);

                                    err=1;
                                    goto clear;
                                }
                            }

                            /*-----------------------------------------------------------------------------
                             * setrow and getrow test: get a row and set a row and compare the matrices
                             *-----------------------------------------------------------------------------*/
                            CALL(mess_matrix_copy(A,Acopy));
                            for(rowidx=0;rowidx<rows;rowidx=1+rowidx*3){

                                CALL(mess_matrix_getrow(A,rowidx,vget));
                                CALL(mess_matrix_setrow(A,rowidx,vget));

                                //compute difference
                                CALL(mess_matrix_diffnormf(A,Acopy,&diff));

                                if(diff>eps){
                                    //test failed print some informations
                                    printf("FAILED setrow/getrow failed with diff=%e\n",diff);
                                    printf("Matrix A:\n");
                                    mess_matrix_printinfo(A);
                                    mess_matrix_print(A);
                                    printf("Matrix Acopy:\n");
                                    mess_matrix_printinfo(Acopy);
                                    mess_matrix_print(Acopy);
                                    printf("Tried to set row ="MESS_PRINTF_INT"\n",rowidx);

                                    err=1;
                                    goto clear;
                                }

                            }
                        }
                    }
                }
            }
        }
    }

clear:
    /*-----------------------------------------------------------------------------
     *  clear matrices and vectors
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&A,&Acopy);
    MESS_CLEAR_VECTORS(&vset,&vget,&vtmp);

    return err;
}

