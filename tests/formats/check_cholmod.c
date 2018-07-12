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
 * @addtogroup test_format
 * @{
 * @file tests/formats/check_cholmod.c
 * @brief Check cholmod.c functions.
 * @author  @mbehr
 * @test
 *
 * This function checks all function defined in cholmod.c that means it checks if a
 * @c cholmod_dense or @c cholmod_sparse
 * datastructure can be converted in a @ref MESS_DENSE @ref mess_matrix or a
 * @ref MESS_CSC @ref mess_matrix structure, respectively and vice versa. \n
 * It also checks if a @c cholmod_dense datastructure can be converted to a @ref mess_vector.
 *
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/interface_cholmod.h"
#include "../call_macro.h"
#include <complex.h>

#define CHECK_CONVERTCHOLMOD_DENSEMATRIX(A,A_CHOLMOD,TMP,C,ERR,EPS,DIFF){   \
    CALL(mess_matrix_convert_dense_to_cholmod(A,&A_CHOLMOD,&C));        \
    CALL(mess_matrix_convert_cholmod_to_dense(A_CHOLMOD, TMP, &C));     \
    cholmod_l_free_dense(&A_CHOLMOD,&C);                                \
    CALL(mess_matrix_diffnorm(A,TMP,&DIFF));                            \
    if(DIFF>EPS){                                                       \
        printf("CHECK_CONVERTCHOLMOD_DENSEMATRIX failed:\n");           \
        CALL(mess_matrix_printinfo(A));                                 \
        CALL(mess_matrix_printshort(A));                                \
        ERR++;                                                          \
    }                                                                   \
}

#define CHECK_CONVERTCHOLMOD_VECTOR(V,V_CHOLMOD,TMP,C,ERR,EPS,DIFF){        \
    CALL(mess_vector_convert_dense_to_cholmod(V,&V_CHOLMOD,&C));        \
    CALL(mess_vector_convert_cholmod_to_dense(V_CHOLMOD, TMP, &C));     \
    cholmod_l_free_dense(&V_CHOLMOD,&C);                                \
    CALL(mess_vector_diffnorm(V,TMP,&DIFF));                            \
    if(DIFF>EPS){                                                       \
        printf("CHECK_CONVERTCHOLMOD_VECTOR failed:\n");                \
        CALL(mess_vector_printinfo(V));                                 \
        CALL(mess_vector_printshort(V));                                \
        ERR++;                                                          \
    }                                                                   \
}

#define CHECK_CONVERTCHOLMOD_SPARSEMATRIX(A,A_CHOLMOD,TMP,C,ERR,EPS,DIFF){  \
    CALL(mess_matrix_convert_csc_to_cholmod(A,&A_CHOLMOD,&C));          \
    CALL(mess_matrix_convert_cholmod_to_csc(A_CHOLMOD, TMP, &C));       \
    cholmod_l_free_sparse(&A_CHOLMOD,&C);                               \
    CALL(mess_matrix_add(1,A,-1,TMP));                                  \
    CALL(mess_matrix_normf(TMP,&DIFF));                                 \
    if(DIFF>EPS){                                                       \
        printf("CONVERTCHOLMOD_SPARSEMATRIX failed:\n");                \
        CALL(mess_matrix_printinfo(A));                                 \
        CALL(mess_matrix_printshort(A));                                \
        ERR++;                                                          \
    }                                                                   \
}

int cholmod_print_dense_all(cholmod_dense *A) {
    int i, j;
    if(A->xtype==CHOLMOD_REAL){
        double *x = A->x;
        for(i = 0; i < A->nrow; i++) {
            for(j = 0; j < A->ncol; j++) {
                printf("% 10e ", x[i + j * A->d]);
            }
            printf("\n");
        }
    }else if(A->xtype==CHOLMOD_COMPLEX){
        mess_double_cpx_t*x = A->x;
        for(i = 0; i < A->nrow; i++) {
            for(j = 0; j < A->ncol; j++) {
                printf("%10e + %10e*I ", creal(x[i + j * A->d]),cimag(x[i + j * A->d]));
            }
            printf("\n");
        }
    }else{
        //MSG_ERROR("unsupported cholmod datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}

int cholmod_print_sparse_all(cholmod_sparse *A){
    int i,j;
    int     *colptr = A->p; //column pointer
    int     *rowptr = A->i; //row pointer;

    if(A->xtype==CHOLMOD_REAL){
        double  *x = A->x;
        for (i = 0 ; i < A->ncol; i++){
            for (j = colptr[i]; j<colptr[i+1]; j++){
                printf("(%d\t, %d) =\t %10e \n", rowptr[j], i, x[j]);
            }
        }
    }else if(A->xtype==CHOLMOD_COMPLEX){
        mess_double_cpx_t  *x = A->x;
        for (i = 0 ; i < A->ncol; i++){
            for (j = colptr[i]; j<colptr[i+1]; j++){
                printf("(%d\t, %d) =\t %10e + %10e*I \n", rowptr[j], i, creal(x[j]),cimag(x[j]));
            }
        }
    }else{
        //MSG_ERROR("unsupported cholmod datatype\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
}


int main(){

    mess_init();
    int ret=0, err=0;
    mess_int_t rows=3, cols=2;
    double diff, eps=1e-7;

    cholmod_dense *A_cholmod, *v_cholmod;
    cholmod_sparse *A_cholmodsparse;

    /* start CHOLMOD */
    cholmod_common c ;
    cholmod_l_start (&c);
    c.itype=CHOLMOD_LONG;
    c.dtype=CHOLMOD_DOUBLE;
    cholmod_l_print_common("Common",&c);

    mess_vector v_real, v_cpx, v_tmp;
    mess_matrix A_dense_real, A_dense_cpx, A_tmp, A_csc_real, A_csc_cpx;

    /*-----------------------------------------------------------------------------
     *  Init/Load Vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&v_real,&v_cpx,&v_tmp);
    mess_vector_alloc(v_real, rows, MESS_REAL);
    mess_vector_alloc(v_cpx , rows, MESS_COMPLEX);
    mess_vector_alloc(v_tmp , rows, MESS_REAL );
    mess_vector_rand(v_real);
    mess_vector_copy(v_real,v_cpx);
    mess_vector_tocomplex(v_cpx);
    mess_vector_scalec(1+2*I,v_cpx);

    /*-----------------------------------------------------------------------------
     * Init/Load matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A_dense_real, &A_dense_cpx, &A_tmp, &A_csc_real, &A_csc_cpx);

    //dense
    CALL(mess_matrix_rand(A_dense_real, rows, cols, MESS_DENSE, MESS_REAL, 1.0));
    CALL(mess_matrix_rand(A_dense_cpx,  rows, cols, MESS_DENSE, MESS_COMPLEX, 1.0));

    //csc
    CALL(mess_matrix_rand(A_csc_real,   rows, cols, MESS_CSC, MESS_REAL, 1.0));
    CALL(mess_matrix_rand(A_csc_cpx,    rows, cols, MESS_CSC, MESS_COMPLEX, 1.0));

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    //case real matrix -> cholmod -> real matrix
    CHECK_CONVERTCHOLMOD_DENSEMATRIX(A_dense_real,A_cholmod,A_tmp,c,err,eps,diff);

    //case complex matrix -> cholmod -> complex matrix
    CHECK_CONVERTCHOLMOD_DENSEMATRIX(A_dense_cpx,A_cholmod,A_tmp,c,err,eps,diff);

    //case real vector -> cholmod -> real vector
    CHECK_CONVERTCHOLMOD_VECTOR(v_real,v_cholmod,v_tmp,c,err,eps,diff)

        //case complex vector -> cholmod -> complex vector
        CHECK_CONVERTCHOLMOD_VECTOR(v_cpx,v_cholmod,v_tmp,c,err,eps,diff)

        //case real csc matrix -> cholmod_sparse -> real csc matrix
        CHECK_CONVERTCHOLMOD_SPARSEMATRIX(A_csc_real,A_cholmodsparse,A_tmp,c,err,eps,diff);

    //case complex csc matrix -> cholmod_sparse -> complex csc matrix
    CHECK_CONVERTCHOLMOD_SPARSEMATRIX(A_csc_cpx,A_cholmodsparse,A_tmp,c,err,eps,diff);

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    /* finish CHOLMOD */
    cholmod_l_free_dense(&A_cholmod,&c);
    cholmod_l_free_dense(&v_cholmod,&c);
    cholmod_l_free_sparse(&A_cholmodsparse,&c);
    cholmod_l_finish (&c) ;

    CALL(mess_matrix_clear(&A_dense_real));
    CALL(mess_matrix_clear(&A_dense_cpx));
    CALL(mess_matrix_clear(&A_csc_cpx));
    CALL(mess_matrix_clear(&A_csc_real));
    CALL(mess_matrix_clear(&A_tmp));

    CALL(mess_vector_clear(&v_real));
    CALL(mess_vector_clear(&v_cpx));
    CALL(mess_vector_clear(&v_tmp));

    return err ;
}


