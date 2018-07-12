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
 * @file matlab/c_interface/conv/mess_matrix.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"

/**
 * @brief Copy a matrix from @matlab  to a \ref mess_matrix.
 * @param[in] Amatlab input @matlab  matrix
 * @return @c NULL if an error occured or the matrix is empty or a @ref mess_matrix with the content copied from the @matlab  object
 *
 * The @ref mess_matrix_from_mexmess function copies a matrix from @matlab  to @mess. \n
 * If the matrix is sparse data is stored in the compressed column format.
 */

mess_matrix mess_matrix_from_mexmess (const mxArray *Amatlab)
{
    int ret = 0 ;
    mess_matrix A = NULL;
    mess_int_t i,j;
    mwIndex *colptr = NULL;
    mwIndex *rowptr = NULL;
    double *v = NULL;
    double *vi = NULL;
    mess_int_t rows, cols, nnz;

    if ( Amatlab == NULL ) {
        csc_error_message("Amatlab points to NULL.\n");
        return NULL;
    }

    /* Empty matrix is NULL  */
    if ( mxIsEmpty(Amatlab) ) {
        return NULL;
    }

    if ( ! ( mxIsDouble (Amatlab) || mxIsComplex(Amatlab))){
        csc_error_message("The input has to be double or complex.\n");
        return NULL;
    }

    ret = mess_matrix_init(&A);
    if ( ret != 0 ) {
        csc_warn_message("Cannot initialize A.\n");
        return NULL;
    }


    rows = mxGetM (Amatlab) ;
    cols = mxGetN (Amatlab) ;
    if ( mxIsSparse (Amatlab)) {
        /*-----------------------------------------------------------------------------
         *  Copy a sparse matrix from matlab to @ref mess
         *-----------------------------------------------------------------------------*/
        colptr = mxGetJc(Amatlab);
        rowptr = mxGetIr(Amatlab);

        nnz=colptr[cols];

        if (mxIsComplex(Amatlab)) {
            mess_matrix_alloc(A, rows, cols, nnz, MESS_CSC, MESS_COMPLEX);
        }  else {
            mess_matrix_alloc(A, rows, cols, nnz, MESS_CSC, MESS_REAL);
        }

        for ( i = 0; i < A->cols+1 ;  i++) {
            A->colptr[i] = (mess_int_t) colptr[i];
        }
        for ( i = 0; i < A->nnz; i++ ) {
            A->rowptr[i] = (mess_int_t) rowptr[i];
        }
        if ( !mxIsComplex(Amatlab)) {
#if MX_HAS_INTERLEAVED_COMPLEX
            memcpy(A->values, mxGetDoubles(Amatlab), sizeof(double)*A->nnz);
#else
            v = mxGetPr(Amatlab);
            for ( i = 0; i < A->nnz ; i++ ) {
                A->values[i] = v[i];
            }
#endif
        } else {
#if MX_HAS_INTERLEAVED_COMPLEX
            memcpy(A->values_cpx, mxGetComplexDoubles(Amatlab), sizeof(mess_double_cpx_t)*A->nnz);
#else
            v = mxGetPr(Amatlab);
            vi = mxGetPi(Amatlab);
            for ( i = 0; i < A->nnz; i++) {
                A->values_cpx[i] = v[i] + vi[i]*I;
            }
#endif
        }
    } else {
        /*-----------------------------------------------------------------------------
         * copy a dense matrix to @ref mess
         *-----------------------------------------------------------------------------*/
        if (mxIsComplex(Amatlab)) {
            mess_matrix_alloc(A, rows, cols, rows * cols, MESS_DENSE, MESS_COMPLEX);
        }  else {
            mess_matrix_alloc(A, rows, cols, rows * cols, MESS_DENSE, MESS_REAL);
        }

        if ( !mxIsComplex(Amatlab)){
#if MX_HAS_INTERLEAVED_COMPLEX
            v = mxGetDoubles(Amatlab);
#else
            v = mxGetPr(Amatlab);
#endif
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    A->values[i+j*A->ld] = v[i+j*A->rows];
                }
            }
        }else {
#if MX_HAS_INTERLEAVED_COMPLEX
            mess_double_cpx_t * temp = (mess_double_cpx_t*) mxGetComplexDoubles(Amatlab);
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    A->values_cpx[i+j*A->ld] = temp[i+j*A->rows];
                }
            }
#else
            v = mxGetPr(Amatlab);
            vi = mxGetPi(Amatlab);
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    A->values_cpx[i+j*A->ld] = v[i+j*A->rows] + vi[i+j*A->rows] * I;
                }
            }
#endif
        }
    }

    return A;
}

/**
 * @brief Copy a matrix from @mess to @matlab.
 * @param[in] A input mess_matrix
 * @return NULL in case of an error or a @matlab mxArray containing the input matrix copied into the @matlab  format.
 *
 * The @ref mess_matrix_to_mexmess function copies a matrix from @mess to @matlab .
 *
 */

mxArray * mess_matrix_to_mexmess (mess_matrix A)
{
    mxArray *Amatlab ;
    mess_int_t i,j;

    if ( A == NULL ){
        return mxCreateDoubleMatrix(0,0,mxREAL);
    }

    if ( MESS_IS_DENSE(A) ){
        /*-----------------------------------------------------------------------------
         * copy a dense matrix to matlab
         *----------------------------------------------------------------------------*/
        if ( MESS_IS_COMPLEX(A)) {
            Amatlab = mxCreateDoubleMatrix(A->rows, A->cols, mxCOMPLEX);
#if MX_HAS_INTERLEAVED_COMPLEX
            mxComplexDouble * temp_matlab = mxGetComplexDoubles(Amatlab);
            mxComplexDouble  *temp_mess = (mxComplexDouble*) A->values_cpx;  //mxComplexDouble is struct{double, double} in MATLAB
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    temp_matlab[i+j*A->rows] = temp_mess[i+j*A->ld];
                }
            }
#else
            double *v;
            double *vi;
            v = mxGetPr(Amatlab);
            vi = mxGetPi(Amatlab);
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    v[i+j*A->rows] = creal(A->values_cpx[i+j*A->ld]);
                    vi[i+j*A->rows] = cimag(A->values_cpx[i+j*A->ld]);
                }
            }
#endif
        } else {
            Amatlab = mxCreateDoubleMatrix(A->rows, A->cols, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
            mxDouble * temp = mxGetDoubles(Amatlab);
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    temp[i+j*A->rows] = A->values[i+j*A->ld];
                }
            }
#else
            double *v;
            v = mxGetPr(Amatlab);
            for (j = 0; j < A->cols; j++) {
                for (i = 0; i < A->rows; i++) {
                    v[i+j*A->rows] = A->values[i+j*A->ld];
                }
            }
#endif
        }
    } else {
        /*-----------------------------------------------------------------------------
         * copy a sparse matrix to matlab
         *-----------------------------------------------------------------------------*/
        int conv = -1;
        mess_matrix work = NULL;
        if ( !(MESS_IS_CSC(A)) ) {
            conv = 0 ;
            mess_matrix_init(&work);
            mess_matrix_convert(A, work, MESS_CSC);
        } else {
            work = A;
        }

        if ( MESS_IS_COMPLEX(A)) {
            Amatlab = mxCreateSparse(work->rows, work->cols, work->nnz, mxCOMPLEX);
            mwIndex *colptr, *rowptr;
            colptr = mxGetJc(Amatlab);
            rowptr = mxGetIr(Amatlab);
            for ( i = 0; i < work->cols+1 ; i++ ) {
                colptr[i] = (mwIndex) work->colptr[i];
            }
            for ( i = 0 ; i< work->nnz ; i++ ) {
                rowptr[i] = (mwIndex) work->rowptr[i];
            }
#if MX_HAS_INTERLEAVED_COMPLEX
            mxComplexDouble *temp = mxGetComplexDoubles(Amatlab);
            memcpy(temp, work->values_cpx, sizeof(mess_double_cpx_t)*work->nnz);
#else
            double *v, *vi;
            v = mxGetPr(Amatlab);
            vi = mxGetPi(Amatlab);
            for ( i = 0 ; i< work->nnz ; i++ ) {
                v [i] = creal(work->values_cpx[i]);
                vi[i] = cimag(work->values_cpx[i]);
            }
#endif
        } else {
            Amatlab = mxCreateSparse(work->rows, work->cols, work->nnz, mxREAL);
            mwIndex *colptr, *rowptr;
            colptr = mxGetJc(Amatlab);
            rowptr = mxGetIr(Amatlab);
            for ( i = 0; i < work->cols+1 ; i++ ) {
                colptr[i] = (mwIndex) work->colptr[i];
            }
            for ( i = 0 ; i< work->nnz ; i++ ) {
                rowptr[i] = (mwIndex) work->rowptr[i];
            }
#if MX_HAS_INTERLEAVED_COMPLEX
            mxDouble * temp = mxGetDoubles(Amatlab);
            memcpy(temp,work->values,sizeof(double)*work->nnz);
#else
            double *v;
            v = mxGetPr(Amatlab);
            for ( i = 0 ; i< work->nnz ; i++ ) {
                v [i] = work->values[i];
            }
#endif
        }

        if ( work != A ){
            mess_matrix_clear(&work);
        }

    }
    return Amatlab;
}


