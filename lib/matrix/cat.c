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
 * @file lib/matrix/cat.c
 * @brief Concatenate up to four matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"


/**
 * @internal
 * @brief Concatenate two matrices in diagonal.
 * @param[in] A         input upper left matrix
 * @param[in] D         input lower right matrix
 * @param[in] otype     input output storage type
 * @param[out] O        input output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_diag function concatenates two matrices in the following way:
 * @verbatim O = [A, 0; 0, D]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_diag ( mess_matrix A, mess_matrix D, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workA = NULL;
    mess_matrix workD = NULL;
    int convA = -1;
    int convD = -1;
    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(D);
    mess_check_nullpointer(O);


    /*-----------------------------------------------------------------------------
     *  Dense Output
     *-----------------------------------------------------------------------------*/

    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(D)) {
            ret = mess_matrix_alloc(O, A->rows+D->rows, workA->cols+workD->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+D->rows, workA->cols+workD->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        // COMPLEX output
        if ( MESS_IS_COMPLEX(O) ){
            if ( MESS_IS_REAL(A)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j=0;j<A->cols;j++){
                    for(i=0;i<A->rows;i++){
                        O->values_cpx[j*O->ld+i]=workA->values[j*workA->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
            }
            roff =A->rows;
            coff =A->cols;
            if (MESS_IS_REAL(D)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j=0;j<D->cols;j++){
                    for(i=0;i<D->rows;i++){
                        O->values_cpx[(j+coff)*O->ld+(roff+i)] = workD->values[j*workD->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&D->rows, &D->cols, workD->values_cpx, &workD->ld, O->values_cpx+(A->cols*O->ld+A->rows), &O->ld);
            }
        }
        // REAL output
        else {
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&D->rows, &D->cols, workD->values, &workD->ld, O->values+(A->cols*O->ld+A->rows), &O->ld);
        }
        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convD != -1){mess_matrix_clear(&workD); }
    }
    /*-----------------------------------------------------------------------------
     *  Sparse CSR output
     *-----------------------------------------------------------------------------*/
    else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);
        //printf("D=\n");
        //mess_matrix_printdata(workD);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+D->cols,workA->nnz+workD->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+D->cols,workA->nnz+workD->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        roff =A->rows;
        coff =A->cols;

        for (i=0;i<D->rows;i++){
            for(j=workD->rowptr[i];j<workD->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->colptr[onz] = coff+workD->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convD != -1){mess_matrix_clear(&workD); }
    }

    /*-----------------------------------------------------------------------------
     *  Sparse CSC output
     *-----------------------------------------------------------------------------*/
    else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+D->cols,workA->nnz+workD->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+D->cols,workA->nnz+workD->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->colptr[0]=0;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }
        roff =A->rows;
        coff =A->cols;

        for (i=0;i<D->cols;i++){
            for(j=workD->colptr[i];j<workD->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->rowptr[onz] = roff+workD->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convD != -1){mess_matrix_clear(&workD); }
    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }

    return( 0) ;
}       /* -----  end of static function __mess_matrix_cat_diag  ----- */


/**
 * @internal
 * @brief Concatenate two matrices on the anti diagonal.
 * @param[in] B         input upper right matrix
 * @param[in] C         input lower left matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_antidiag function concatenates two matrices in the following way:
 * @verbatim O = [0, B; C, 0]; @endverbatim
 * @attention Internal use only.
 */
static int  __mess_matrix_cat_antidiag ( mess_matrix B, mess_matrix C, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workB = NULL;
    mess_matrix workC = NULL;
    int convB = -1;
    int convC = -1;
    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(O);


    /*-----------------------------------------------------------------------------
     *  Dense Output
     *-----------------------------------------------------------------------------*/

    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        if ( MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)) {
            ret = mess_matrix_alloc(O, B->rows+C->rows, workB->cols+workC->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, B->rows+C->rows, workB->cols+workC->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        // COMPLEX output
        if ( MESS_IS_COMPLEX(O) ){
            if ( MESS_IS_REAL(B)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j=0;j<B->cols;j++){
                    for(i=0;i<B->rows;i++){
                        O->values_cpx[(j+C->cols)*O->ld+i]=workB->values[j*workB->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&B->rows, &B->cols, workB->values_cpx, &workB->ld, O->values_cpx+O->ld*C->cols, &O->ld);
            }
            roff =B->rows;
            if (MESS_IS_REAL(C)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j=0;j<C->cols;j++){
                    for(i=0;i<C->rows;i++){
                        O->values_cpx[j*O->ld+(roff+i)] = workC->values[j*workC->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&C->rows, &C->cols, workC->values_cpx, &workC->ld, O->values_cpx+B->rows, &O->ld);
            }
        }
        // REBL output
        else {
            F77_GLOBAL(dlacpy,DLACPY)("A",&B->rows, &B->cols, workB->values, &workB->ld, O->values+C->cols*O->ld, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&C->rows, &C->cols, workC->values, &workC->ld, O->values+B->rows, &O->ld);
        }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    }
    /*-----------------------------------------------------------------------------
     *  Sparse CSR output
     *-----------------------------------------------------------------------------*/
    else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        // MSG_PRINT("D=\n");
        // mess_matrix_printdata(workD);
        if ( MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, B->rows+C->rows, B->cols+C->cols,workB->nnz+workC->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, B->rows+C->rows, B->cols+C->cols,workB->nnz+workC->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        roff =B->rows;
        coff =C->cols;
        for (i=0;i<B->rows;i++){
            for(j=workB->rowptr[i];j<workB->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->colptr[onz] = coff+workB->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }

        for (i=0;i<C->rows;i++){
            for(j=workC->rowptr[i];j<workC->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->colptr[onz] = workC->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    }

    /*-----------------------------------------------------------------------------
     *  Sparse CSC output
     *-----------------------------------------------------------------------------*/
    else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        if ( MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, B->rows+C->rows, B->cols+C->cols,workB->nnz+workC->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, B->rows+C->rows, B->cols+C->cols,workB->nnz+workC->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->colptr[0]=0;
        roff =B->rows;
        coff =C->cols;


        for (i=0;i<C->cols;i++){
            for(j=workC->colptr[i];j<workC->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->rowptr[onz] = roff+workC->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }
        for (i=0;i<B->cols;i++){
            for(j=workB->colptr[i];j<workB->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->rowptr[onz] = workB->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }



    return 0;
}       /* -----  end of function __mess_matrix_cat_antidiag  ----- */


/**
 * @internal
 * @brief Concatenate two matrices in a block column.
 * @param[in] A         input upper left matrix
 * @param[in] C         input lower left matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_firstcol function concatenates two matrices in the following
 * way: @verbatim O = [A; C]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_firstcol ( mess_matrix A, mess_matrix C, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workA = NULL;
    mess_matrix workC = NULL;
    int convA = -1;
    int convC = -1;
    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(C);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense sparse output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        if ( MESS_IS_REAL(O)) {
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&C->rows, &C->cols, workC->values, &workC->ld, O->values+A->rows,&O->ld);

        } else {
            if (MESS_IS_REAL(A)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j = 0; j<A->cols; j++){
                    for ( i = 0; i < A->rows; i++){
                        O->values_cpx[j*O->ld+i] = workA->values[j*workA->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
            }
            if (MESS_IS_REAL(C)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j = 0; j<C->cols; j++){
                    for ( i = 0; i < C->rows; i++){
                        O->values_cpx[j*O->ld+i+A->rows] = workC->values[j*workC->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&C->rows, &C->cols, workC->values_cpx, &workC->ld, O->values_cpx+A->rows,&O->ld);
            }
        }
        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    }
    /*-----------------------------------------------------------------------------
     *   sparse csr output
     *-----------------------------------------------------------------------------*/
    else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols,workA->nnz+workC->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols,workA->nnz+workC->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX (A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if (MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        roff =A->rows;

        for (i=0;i<C->rows;i++){
            for(j=workC->rowptr[i];j<workC->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX (C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if (MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->colptr[onz] = workC->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    }
    /*-----------------------------------------------------------------------------
     *   sparse csc output
     *-----------------------------------------------------------------------------*/
    else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols,workA->nnz+workC->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols,workA->nnz+workC->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        O->colptr[0]=0;
        roff =A->rows;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX (A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if (MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            for(j=workC->colptr[i];j<workC->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX (C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if (MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->rowptr[onz] = roff+workC->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return( 0 );
}       /* -----  end of static function __mess_matrix_cat_firstcol  ----- */

/**
 * @internal
 * @brief Concatenates two matrices in the first block row.
 * @param[in] A         input upper left matrix
 * @param[in] B         input upper right matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_firstrow function concatenates two matrices in the following way:
 * @verbatim O = [A,B]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_firstrow ( mess_matrix A, mess_matrix B, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workA = NULL;
    mess_matrix workB = NULL;
    int convA = -1;
    int convB = -1;
    int ret = 0;
    mess_int_t i,j;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B)){
            ret = mess_matrix_alloc(O, A->rows, A->cols+B->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows, A->cols+B->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        if ( MESS_IS_REAL(O)){
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&B->rows, &B->cols, workB->values, &workB->ld, O->values+A->cols*O->ld,&O->ld);
        } else {
            if (MESS_IS_REAL(A)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j = 0; j<A->cols; j++){
                    for ( i = 0; i < A->rows; i++){
                        O->values_cpx[j*O->ld+i] = workA->values[j*workA->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
            }
            if (MESS_IS_REAL(B)){
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i,j)
#endif
                for (j = 0; j<B->cols; j++){
                    for ( i = 0; i < B->rows; i++){
                        O->values_cpx[(A->cols+j)*O->ld+i] = workB->values[j*workB->ld+i];
                    }
                }
            } else {
                F77_GLOBAL(zlacpy,ZLACPY)("A",&B->rows, &B->cols, workB->values_cpx, &workB->ld, O->values_cpx+A->cols*O->ld,&O->ld);
            }
        }
        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }

        /*-----------------------------------------------------------------------------
         *  CSC sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B)){
            ret = mess_matrix_alloc(O, A->rows, A->cols+B->cols,workA->nnz+workB->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows, A->cols+B->cols,workA->nnz+workB->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->colptr[0]=0;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }
        coff =A->cols;

        for (i=0;i<B->cols;i++){
            for(j=workB->colptr[i];j<workB->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->rowptr[onz] = workB->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }

        /*-----------------------------------------------------------------------------
         * CSR sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B)){
            ret = mess_matrix_alloc(O, A->rows, A->cols+B->cols,workA->nnz+workB->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows, A->cols+B->cols,workA->nnz+workB->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        O->rowptr[0]=0;
        coff =A->cols;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            for(j=workB->rowptr[i];j<workB->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->colptr[onz] = coff+workB->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return(0);
}

/**
 * @internal
 * @brief Concatenates four matrices.
 * @param[in] A         input upper left matrix
 * @param[in] B         input upper right matrix
 * @param[in] C         input lower left matrix
 * @param[in] D         input lower right matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_full function concatenates 4 matrices in a square:
 * @verbatim O = [A,B; C,D]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_full ( mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix D, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workA = NULL;
    mess_matrix workB = NULL;
    mess_matrix workC = NULL;
    mess_matrix workD = NULL;
    int convA = -1;
    int convB = -1;
    int convC = -1;
    int convD = -1;

    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(D);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)|| MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff=A->rows;
        coff=A->cols;
        if ( MESS_IS_REAL(O)){
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&B->rows, &B->cols, workB->values, &workB->ld, O->values+coff*O->ld, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&C->rows, &C->cols, workC->values, &workC->ld, O->values+roff, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&D->rows, &D->cols, workD->values, &workD->ld, O->values+coff*O->ld+roff, &O->ld);

        } else {
#ifdef _OPENMP
#pragma omp parallel sections default(shared) private(i,j)
#endif
            {
#ifdef _OPENMP
#pragma omp section
#endif
                {
                    if ( MESS_IS_REAL(A)) {
                        for (j = 0; j<A->cols; j++){
                            for ( i = 0; i < A->rows; i++){
                                O->values_cpx[j*O->ld+i] = workA->values[j*workA->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {
                    if ( MESS_IS_REAL(B)) {
                        for (j = 0; j<B->cols; j++){
                            for ( i = 0; i < B->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i] = workB->values[j*workB->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&B->rows, &B->cols, workB->values_cpx, &workB->ld, O->values_cpx+coff*O->ld, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {
                    if ( MESS_IS_REAL(C)) {
                        for (j = 0; j<C->cols; j++){
                            for ( i = 0; i < C->rows; i++){
                                O->values_cpx[j*O->ld+i+roff] = workC->values[j*workC->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&C->rows, &C->cols, workC->values_cpx, &workC->ld, O->values_cpx+roff, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {
                    if ( MESS_IS_REAL(D)) {
                        for (j = 0; j<D->cols; j++){
                            for ( i = 0; i < D->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i+roff] = workD->values[j*workD->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&D->rows, &D->cols, workD->values_cpx, &workD->ld, O->values_cpx+coff*O->ld+roff, &O->ld);
                    }
                }
            }

        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         *  CSC sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)|| MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz+workD->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz+workD->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->colptr[0]=0;
        roff =A->rows;
        coff =A->cols;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            for(j=workC->colptr[i];j<workC->colptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }

                O->rowptr[onz] = roff+workC->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }

        for (i=0;i<B->cols;i++){
            for(j=workB->colptr[i];j<workB->colptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }

                O->rowptr[onz] = workB->rowptr[j];
                onz++;
            }
            for(j=workD->colptr[i];j<workD->colptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->rowptr[onz] = roff+workD->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         * CSR sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)|| MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz+workD->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz+workD->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        roff =A->rows;
        coff =A->cols;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            for(j=workB->rowptr[i];j<workB->rowptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->colptr[onz] = coff+workB->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        for (i=0;i<C->rows;i++){
            for(j=workC->rowptr[i];j<workC->rowptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->colptr[onz] = workC->colptr[j];
                onz++;
            }
            for(j=workD->rowptr[i];j<workD->rowptr[i+1];j++){
                if (MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->colptr[onz] = coff+workD->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }


    } else {
        MSG_ERROR("unsupported storage type: %s\n", mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return(0);
}

/**
 * @internal
 * @brief Concatenate three matrices as upper left triangle.
 * @param[in] A         input upper left matrix
 * @param[in] B         input upper right matrix
 * @param[in] C         input lower left matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref __mess_matrix_cat_LABC function concatenates three  matrices like:
 *  @verbatim O = [A,B;C,0]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_LABC ( mess_matrix A, mess_matrix B, mess_matrix C, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__) ;
    mess_matrix workA = NULL;
    mess_matrix workB = NULL;
    mess_matrix workC = NULL;
    int convA = -1;
    int convB = -1;
    int convC = -1;

    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff=A->rows;
        coff=A->cols;
        if ( MESS_IS_REAL(O)){
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&B->rows, &B->cols, workB->values, &workB->ld, O->values+coff*O->ld, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&C->rows, &C->cols, workC->values, &workC->ld, O->values+roff, &O->ld);
        } else {
#ifdef _OPENMP
#pragma omp parallel sections default(shared) private(i,j)
#endif
            {
#ifdef _OPENMP
#pragma omp section
#endif
                {
                    if ( MESS_IS_REAL(A)) {
                        for (j = 0; j<A->cols; j++){
                            for ( i = 0; i < A->rows; i++){
                                O->values_cpx[j*O->ld+i] = workA->values[j*workA->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(B)) {
                        for (j = 0; j<B->cols; j++){
                            for ( i = 0; i < B->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i] = workB->values[j*workB->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&B->rows, &B->cols, workB->values_cpx, &workB->ld, O->values_cpx+coff*O->ld, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(C)) {
                        for (j = 0; j<C->cols; j++){
                            for ( i = 0; i < C->rows; i++){
                                O->values_cpx[j*O->ld+i+roff] = workC->values[j*workC->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&C->rows, &C->cols, workC->values_cpx, &workC->ld, O->values_cpx+roff, &O->ld);
                    }
                }
            }
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }

        /*-----------------------------------------------------------------------------
         *  CSC sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff =A->rows;
        coff =A->cols;
        O->colptr[0]=0;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            for(j=workC->colptr[i];j<workC->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }

                O->rowptr[onz] = roff+workC->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }

        for (i=0;i<B->cols;i++){
            for(j=workB->colptr[i];j<workB->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }

                O->rowptr[onz] = workB->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }

        /*-----------------------------------------------------------------------------
         * CSR sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+B->cols,workA->nnz+workB->nnz+workC->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        roff =A->rows;
        coff =A->cols;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }

                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            for(j=workB->rowptr[i];j<workB->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }

                O->colptr[onz] = coff+workB->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        for (i=0;i<C->rows;i++){
            for(j=workC->rowptr[i];j<workC->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }

                O->colptr[onz] = workC->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return(0);
}

/**
 * @internal
 * @brief Concatenate three matrices as lower right triangle.
 * @param[in] B         input upper right matrix
 * @param[in] C         input lower left matrix
 * @param[in] D         input lower right matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_RBCD function concatenates three matrices like:
 * @verbatim O = [0,B;C,D]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_RBCD ( mess_matrix B, mess_matrix C, mess_matrix D, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workB = NULL;
    mess_matrix workC = NULL;
    mess_matrix workD = NULL;
    int convB = -1;
    int convC = -1;
    int convD = -1;

    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(D);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, B->rows+C->rows, C->cols+D->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, B->rows+C->rows, C->cols+D->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff=B->rows;
        coff=C->cols;
        if ( MESS_IS_REAL(O)){
            F77_GLOBAL(dlacpy,DLACPY)("A",&B->rows, &B->cols, workB->values, &workB->ld, O->values+coff*O->ld, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&C->rows, &C->cols, workC->values, &workC->ld, O->values+roff, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&D->rows, &D->cols, workD->values, &workD->ld, O->values+coff*O->ld+roff, &O->ld);

        } else {
#ifdef _OPENMP
#pragma omp parallel sections default(shared) private(i,j)
#endif
            {
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(B)) {
                        for (j = 0; j<B->cols; j++){
                            for ( i = 0; i < B->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i] = workB->values[j*workB->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&B->rows, &B->cols, workB->values_cpx, &workB->ld, O->values_cpx+coff*O->ld, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(C)) {
                        for (j = 0; j<C->cols; j++){
                            for ( i = 0; i < C->rows; i++){
                                O->values_cpx[j*O->ld+i+roff] = workC->values[j*workC->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&C->rows, &C->cols, workC->values_cpx, &workC->ld, O->values_cpx+roff, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {
                    if ( MESS_IS_REAL(D)) {
                        for (j = 0; j<D->cols; j++){
                            for ( i = 0; i < D->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i+roff] = workD->values[j*workD->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&D->rows, &D->cols, workD->values_cpx, &workD->ld, O->values_cpx+coff*O->ld+roff, &O->ld);
                    }
                }
            }


        }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         *  CSC sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, B->rows+C->rows, C->cols+B->cols,workB->nnz+workC->nnz+workD->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, B->rows+C->rows, C->cols+B->cols,workB->nnz+workC->nnz+workD->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->colptr[0]=0;
        roff =B->rows;
        coff =C->cols;
        for (i=0;i<C->cols;i++){
            for(j=workC->colptr[i];j<workC->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->rowptr[onz] = roff+workC->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }

        for (i=0;i<B->cols;i++){
            for(j=workB->colptr[i];j<workB->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->rowptr[onz] = workB->rowptr[j];
                onz++;
            }
            for(j=workD->colptr[i];j<workD->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->rowptr[onz] = roff+workD->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         * CSR sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);
        // mess_matrix_printinfo(workB);
        // mess_matrix_printinfo(workC);
        // mess_matrix_printinfo(workD);

        if ( MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(C) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, B->rows+C->rows, C->cols+B->cols,workB->nnz+workC->nnz+workD->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, B->rows+C->rows, C->cols+B->cols,workB->nnz+workC->nnz+workD->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

        // MSG_PRINT("O=\n");
        // mess_matrix_printinfo(O);
        O->rowptr[0]=0;
        roff =B->rows;
        coff =C->cols;
        for (i=0;i<B->rows;i++){
            for(j=workB->rowptr[i];j<workB->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->colptr[onz] = coff+workB->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        for (i=0;i<C->rows;i++){
            for(j=workC->rowptr[i];j<workC->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->colptr[onz] = workC->colptr[j];
                // MSG_PRINT("onz: " MESS_PRINTF_INT "\n", onz);
                onz++;
            }
            for(j=workD->rowptr[i];j<workD->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->colptr[onz] = coff+workD->colptr[j];
                // MSG_PRINT("onz: " MESS_PRINTF_INT "\n", onz);
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return(0);
}

/**
 * @internal
 * @brief Concatenate three matrices as upper right triangle.
 * @param[in] A         input upper left matrix
 * @param[in] B         input upper right matrix
 * @param[in] D         input lower right matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_RABD function concatenates three matrices like:
 * @verbatim O = [A,B;0,D]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_RABD ( mess_matrix A, mess_matrix B, mess_matrix D, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workA = NULL;
    mess_matrix workB = NULL;
    mess_matrix workD = NULL;
    int convA = -1;
    int convB = -1;
    int convD = -1;

    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(D);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+B->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+B->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff=A->rows;
        coff=A->cols;
        if ( MESS_IS_REAL(O)){
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&B->rows, &B->cols, workB->values, &workB->ld, O->values+coff*O->ld, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&D->rows, &D->cols, workD->values, &workD->ld, O->values+coff*O->ld+roff, &O->ld);

        } else {
#ifdef _OPENMP
#pragma omp parallel sections default(shared) private(i,j)
#endif
            {
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(A)) {
                        for (j = 0; j<A->cols; j++){
                            for ( i = 0; i < A->rows; i++){
                                O->values_cpx[j*O->ld+i] = workA->values[j*workA->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(B)) {
                        for (j = 0; j<B->cols; j++){
                            for ( i = 0; i < B->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i] = workB->values[j*workB->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&B->rows, &B->cols, workB->values_cpx, &workB->ld, O->values_cpx+coff*O->ld, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(D)) {
                        for (j = 0; j<D->cols; j++){
                            for ( i = 0; i < D->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i+roff] = workD->values[j*workD->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&D->rows, &D->cols, workD->values_cpx, &workD->ld, O->values_cpx+coff*O->ld+roff, &O->ld);
                    }
                }
            }

        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         *  CSC sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+B->cols,workA->nnz+workB->nnz+workD->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+B->cols,workA->nnz+workB->nnz+workD->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->colptr[0]=0;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }
        roff =A->rows;
        coff =A->cols;

        for (i=0;i<B->cols;i++){
            for(j=workB->colptr[i];j<workB->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->rowptr[onz] = workB->rowptr[j];
                onz++;
            }
            for(j=workD->colptr[i];j<workD->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }

                O->rowptr[onz] = roff+workD->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         * CSR sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(B,workB,convB,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(B) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+B->cols,workA->nnz+workB->nnz+workD->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+D->rows, A->cols+B->cols,workA->nnz+workB->nnz+workD->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        roff =A->rows;
        coff =A->cols;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            for(j=workB->rowptr[i];j<workB->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(B)){
                    O->values_cpx[onz] = workB->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(B)){
                    O->values_cpx[onz] = workB->values[j];
                } else {
                    O->values[onz] = workB->values[j];
                }
                O->colptr[onz] = coff+workB->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        for (i=0;i<D->rows;i++){
            for(j=workD->rowptr[i];j<workD->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->colptr[onz] = coff+workD->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convB != -1){mess_matrix_clear(&workB); }
        if ( convD != -1){mess_matrix_clear(&workD); }


    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return(0);
}

/**
 * @internal
 * @brief Concatenate three matrices as lower left triangle.
 * @param[in] A         input upper left matrix
 * @param[in] C         input lower left matrix
 * @param[in] D         input lower right matrix
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref __mess_matrix_cat_LACD function concatenates three matrices like:
 * @verbatim O = [A,0;C,D]; @endverbatim
 * @attention Internal use only.
 */
static int __mess_matrix_cat_LACD ( mess_matrix A, mess_matrix C, mess_matrix D, unsigned short otype, mess_matrix O )
{
    MSG_FNAME(__func__);
    mess_matrix workA = NULL;
    mess_matrix workC = NULL;
    mess_matrix workD = NULL;
    int convA = -1;
    int convC = -1;
    int convD = -1;

    int ret = 0;
    mess_int_t i,j;
    mess_int_t roff=0;
    mess_int_t coff=0;

    /*-----------------------------------------------------------------------------
     *  recheck some input
     *  rest is done by @ref mess_matrix_cat
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(C);
    mess_check_nullpointer(D);
    mess_check_nullpointer(O);

    /*-----------------------------------------------------------------------------
     *  dense output
     *-----------------------------------------------------------------------------*/
    if (otype == MESS_DENSE){
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(C) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+D->cols,0, MESS_DENSE, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+D->cols,0, MESS_DENSE, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff=A->rows;
        coff=A->cols;
        if ( MESS_IS_REAL(O)){
            F77_GLOBAL(dlacpy,DLACPY)("A",&A->rows, &A->cols, workA->values, &workA->ld, O->values, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&C->rows, &C->cols, workC->values, &workC->ld, O->values+roff, &O->ld);
            F77_GLOBAL(dlacpy,DLACPY)("A",&D->rows, &D->cols, workD->values, &workD->ld, O->values+coff*O->ld+roff, &O->ld);

        } else {
#ifdef _OPENMP
#pragma omp parallel sections default(shared) private(i,j)
#endif
            {
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(A)) {
                        for (j = 0; j<A->cols; j++){
                            for ( i = 0; i < A->rows; i++){
                                O->values_cpx[j*O->ld+i] = workA->values[j*workA->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&A->rows, &A->cols, workA->values_cpx, &workA->ld, O->values_cpx, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(C)) {
                        for (j = 0; j<C->cols; j++){
                            for ( i = 0; i < C->rows; i++){
                                O->values_cpx[j*O->ld+i+roff] = workC->values[j*workC->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&C->rows, &C->cols, workC->values_cpx, &workC->ld, O->values_cpx+roff, &O->ld);
                    }
                }
#ifdef _OPENMP
#pragma omp section
#endif
                {

                    if ( MESS_IS_REAL(D)) {
                        for (j = 0; j<D->cols; j++){
                            for ( i = 0; i < D->rows; i++){
                                O->values_cpx[(coff+j)*O->ld+i+roff] = workD->values[j*workD->ld+i];
                            }
                        }
                    } else {
                        F77_GLOBAL(zlacpy,ZLACPY)("A",&D->rows, &D->cols, workD->values_cpx, &workD->ld, O->values_cpx+coff*O->ld+roff, &O->ld);
                    }
                }
            }

        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         *  CSC sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSC) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(C) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+D->cols,workA->nnz+workC->nnz+workD->nnz, MESS_CSC, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+D->cols,workA->nnz+workC->nnz+workD->nnz, MESS_CSC, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        roff =A->rows;
        coff =A->cols;

        O->colptr[0]=0;
        for (i=0;i<A->cols;i++){
            for(j=workA->colptr[i];j<workA->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->rowptr[onz] = workA->rowptr[j];
                onz++;
            }
            for(j=workC->colptr[i];j<workC->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }
                O->rowptr[onz] = roff+workC->rowptr[j];
                onz++;
            }
            O->colptr[i+1] = onz;
        }

        for (i=0;i<D->cols;i++){
            for(j=workD->colptr[i];j<workD->colptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->rowptr[onz] = roff+workD->rowptr[j];
                onz++;
            }
            O->colptr[coff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }

        /*-----------------------------------------------------------------------------
         * CSR sparse output
         *-----------------------------------------------------------------------------*/
    } else if ( otype == MESS_CSR) {
        mess_int_t onz = 0;
        MESS_MATRIX_CHECKFORMAT(A,workA,convA,otype);
        MESS_MATRIX_CHECKFORMAT(C,workC,convC,otype);
        MESS_MATRIX_CHECKFORMAT(D,workD,convD,otype);

        if ( MESS_IS_COMPLEX(A) || MESS_IS_COMPLEX(C) || MESS_IS_COMPLEX(D)){
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+D->cols,workA->nnz+workC->nnz+workD->nnz, MESS_CSR, MESS_COMPLEX);
        } else {
            ret = mess_matrix_alloc(O, A->rows+C->rows, A->cols+D->cols,workA->nnz+workC->nnz+workD->nnz, MESS_CSR, MESS_REAL);
        }
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        O->rowptr[0]=0;
        roff =A->rows;
        coff =A->cols;
        for (i=0;i<A->rows;i++){
            for(j=workA->rowptr[i];j<workA->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(A)){
                    O->values_cpx[onz] = workA->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(A)){
                    O->values_cpx[onz] = workA->values[j];
                } else {
                    O->values[onz] = workA->values[j];
                }
                O->colptr[onz] = workA->colptr[j];
                onz++;
            }
            O->rowptr[i+1] = onz;
        }
        for (i=0;i<C->rows;i++){
            for(j=workC->rowptr[i];j<workC->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(C)){
                    O->values_cpx[onz] = workC->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(C)){
                    O->values_cpx[onz] = workC->values[j];
                } else {
                    O->values[onz] = workC->values[j];
                }

                O->colptr[onz] = workC->colptr[j];
                onz++;
            }
            for(j=workD->rowptr[i];j<workD->rowptr[i+1];j++){
                if ( MESS_IS_COMPLEX(O) && MESS_IS_COMPLEX(D)){
                    O->values_cpx[onz] = workD->values_cpx[j];
                } else if ( MESS_IS_COMPLEX(O) && MESS_IS_REAL(D)){
                    O->values_cpx[onz] = workD->values[j];
                } else {
                    O->values[onz] = workD->values[j];
                }
                O->colptr[onz] = coff+workD->colptr[j];
                onz++;
            }
            O->rowptr[roff+i+1] = onz;
        }

        if ( convA != -1){mess_matrix_clear(&workA); }
        if ( convC != -1){mess_matrix_clear(&workC); }
        if ( convD != -1){mess_matrix_clear(&workD); }


    } else {
        MSG_ERROR("unsupported storage type: %s\n",mess_storage_t_str(otype));
        return(MESS_ERROR_STORAGETYPE);
    }


    return(0);
}



/**
 * @brief Concatenate up to four matrices.
 * @param[in] A         input matrix \f$A\f$
 * @param[in] B         input matrix \f$B\f$
 * @param[in] C         input matrix \f$C\f$
 * @param[in] D         input matrix \f$D\f$
 * @param[in] otype     input output storage type
 * @param[out] O        output matrix \f$O\f$
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref mess_matrix_cat function concatenates 4 matrices like:
 * \f[  O = \begin{bmatrix} A & B \\ C & D \end{bmatrix}. \f]
 *
 * If a block does not exists it must be set to @c NULL. \n
 * For example:
 * @code{.c}
 * mess_matrix_cat(NULL,B,C,NULL,MESS_CSR,O);
 * @endcode
 * performs
 * \f[0=\begin{bmatrix} 0 & B \\ C & 0 \end{bmatrix} \f]
 * and in addition the output has @ref MESS_CSR storage type.
 * 
 * @attention This function does not yet support the MESS_COORD storage type as output type.
 *
 */
int mess_matrix_cat ( mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix D, mess_storage_t otype, mess_matrix O)
{
    MSG_FNAME(__func__);
    mess_int_t rA=0,rB=0,rC=0,rD=0;
    mess_int_t cA=0,cB=0,cC=0,cD=0;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(O);
    MESS_MATRIX_RESET(O);

    if ( A!=NULL){
        rA=A->rows;
        cA=A->cols;
    }
    if ( B!=NULL){
        rB=B->rows;
        cB=B->cols;
    }
    if ( C!=NULL){
        rC=C->rows;
        cC=C->cols;
    }
    if ( D!=NULL){
        rD=D->rows;
        cD=D->cols;
    }

    /*-----------------------------------------------------------------------------
     *  different cases
     *-----------------------------------------------------------------------------*/
    if ( A!=NULL && B==NULL && C==NULL && D==NULL){
        ret =  mess_matrix_convert(A,O,otype);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    } else if ( A==NULL && B!=NULL && C==NULL && D==NULL){
        ret =  mess_matrix_convert(B,O,otype);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    } else if ( A==NULL && B==NULL && C!=NULL && D==NULL){
        ret =  mess_matrix_convert(C,O,otype);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    } else if ( A==NULL && B==NULL && C==NULL && D!=NULL){
        ret =  mess_matrix_convert(D,O,otype);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
    } else if ( A!=NULL && B==NULL && C==NULL && D!= NULL){
        MSG_INFO("O = [ A 0; 0 D]\n");
        ret = __mess_matrix_cat_diag(A,D,otype,O);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_cat_diag);
    } else if ( A == NULL && D==NULL && B!=NULL && C!=NULL) {
        MSG_INFO("O = [ 0 B; C 0]\n");
        ret = __mess_matrix_cat_antidiag(B,C,otype,O);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_cat_antidiag);
    } else if ( A!=NULL && B==NULL && C!=NULL && D==NULL){
        if ( cA!=cC) {
            MSG_ERROR("A and C must have the same number of columns.\n");
            return( MESS_ERROR_DIMENSION );
        }
        ret = __mess_matrix_cat_firstcol(A,C,otype,O);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_cat_firstcol);
    } else if ( A!=NULL && B!=NULL && C==NULL && D==NULL){
        if ( rA==rB ) {
            ret = __mess_matrix_cat_firstrow(A,B,otype,O);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mess_matrix_cat_firstrow);
        } else if (rA == 0) {
            ret = mess_matrix_convert(B,O,otype);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
        } else {
            MSG_ERROR("A and B must have the same number of rows.\n");
            return( MESS_ERROR_DIMENSION ) ;
        }
    } else if  ( A!=NULL && B!=NULL && C!=NULL && D!=NULL){
        if (rA!=rB || cA != cC || rC != rD || cB != cD){
            MSG_ERROR("Dimension of the matrix have to fit.\n");
            return(MESS_ERROR_DIMENSION);
        }
        ret = __mess_matrix_cat_full(A,B,C,D,otype,O);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __mess_matrix_cat_full);
    } else if  ( A!=NULL && B!=NULL && C!=NULL && D==NULL){
        if (rA!=rB || cA != cC ){
            MSG_ERROR("Dimension of the matrix have to fit.\n");
            return(MESS_ERROR_DIMENSION);
        }
        ret = __mess_matrix_cat_LABC(A,B,C,otype,O);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __mess_matrix_cat_LABC);
    }  else if  ( A==NULL && B!=NULL && C!=NULL && D!=NULL){
        if (cB!=cD || rC != rD ){
            MSG_ERROR("Dimension of the matrix have to fit.\n");
            return(MESS_ERROR_DIMENSION);
        }
        ret = __mess_matrix_cat_RBCD(B,C,D,otype,O);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __mess_matrix_cat_RBCD);
    }  else if  ( A!=NULL && B!=NULL && C==NULL && D!=NULL){
        if (cB!=cD || rA != rB ){
            MSG_ERROR("Dimension of the matrix have to fit.\n");
            return(MESS_ERROR_DIMENSION);
        }
        ret = __mess_matrix_cat_RABD(A,B,D,otype,O);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __mess_matrix_cat_RABD);
    } else if  ( A!=NULL && B==NULL && C!=NULL && D!=NULL){
        if (cA!=cC || rC != rD ){
            MSG_ERROR("Dimension of the matrix have to fit.\n");
            return(MESS_ERROR_DIMENSION);
        }
        ret = __mess_matrix_cat_LACD(A,C,D,otype,O);
        FUNCTION_FAILURE_HANDLE(ret,(ret!=0), __mess_matrix_cat_LACD);
    }
    else{
        MSG_ERROR("to be implemented\n");
        return(MESS_ERROR_NOSUPPORT);
    }




    return(0);
}


