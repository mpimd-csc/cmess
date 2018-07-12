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
 * @file lib/precond/iluk.c
 * @brief Generate an \f$ ILU(k) \f$ preconditioner based on @cite Saa03a.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


#ifdef _OPENMP_H
#include <omp.h>
#endif


#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

typedef struct CCSRFmt {
    mess_datatype_t data_type;
    mess_int_t n;
    mess_int_t  *rowptr;
    mess_int_t  *colptr;
    double *values;
    mess_double_cpx_t *values_cpx;
    mess_int_t nnz;
} CCSRMat, *ccsrptr;

typedef struct __iluk {
    mess_datatype_t data_type;
    mess_int_t n;
    double *D;
    mess_double_cpx_t *Dx;
    ccsrptr L;
    ccsrptr U;
} ILUK, *ilukptr;

static int symbolic_iluk( mess_int_t lofM, mess_matrix mat, ilukptr ilu );
static int __iluk_solve_complex( mess_precond myself, mess_solver_options opt, mess_vector in, mess_vector out);

/*-----------------------------------------------------------------------------
 *  Solve with a real ILU
 *-----------------------------------------------------------------------------*/
static int __iluk_solve( mess_precond myself, mess_solver_options opt, mess_vector in, mess_vector out){
    MSG_FNAME(__func__);
    ilukptr ilu;
    mess_int_t n , i, j;
    int ret = 0;
    double *D, *x,*y;
    ccsrptr L, U;

    /*-----------------------------------------------------------------------------
     *  check
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(myself);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);

    ilu = (ilukptr) myself->data;
    mess_check_nullpointer(ilu);
    if ( MESS_IS_COMPLEX(ilu)) {
        return __iluk_solve_complex(myself, opt, in, out);
    }

    n= ilu->n;
    L = ilu->L;
    U = ilu->U;
    D = ilu->D;

    if ( n != in->dim) {
        MSG_ERROR("input has the wrong number of rows. is = " MESS_PRINTF_INT ",  wanted = " MESS_PRINTF_INT "\n", in->dim, n);
        return MESS_ERROR_DIMENSION;
    }
    if (n != out->dim) {
        MSG_ERROR("output has the wrong number of rows. is = " MESS_PRINTF_INT ", wanted = " MESS_PRINTF_INT "\n", out->dim, n );
        return MESS_ERROR_DIMENSION;
    }
    if ( MESS_IS_REAL(in) ) {
        ret =  mess_vector_toreal_nowarn(out);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        x = out -> values;
        y = in -> values;


        /*-----------------------------------------------------------------------------
         *  solve with L
         *-----------------------------------------------------------------------------*/
        for( i = 0; i < n; i++ ) {
            x[i] = y[i];
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {
                x[i] -= x[L->colptr[j]] * L->values[j];
            }
        }

        /*-----------------------------------------------------------------------------
         *  solve with U
         *-----------------------------------------------------------------------------*/
        for( i = n-1; i >= 0; i-- ) {
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {
                x[i] -= x[U->colptr[j]] * U->values[j];
            }
            x[i] *= D[i];
        }
    } else if ( MESS_IS_COMPLEX(in)) {
        mess_double_cpx_t  *xc, *yc;
        ret = mess_vector_tocomplex(out);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex );
        xc = out -> values_cpx;
        yc = in -> values_cpx;


        /*-----------------------------------------------------------------------------
         *  solve with L
         *-----------------------------------------------------------------------------*/
        for( i = 0; i < n; i++ ) {
            xc[i] = yc[i];
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {
                xc[i] -= xc[L->colptr[j]] * L->values[j];
            }
        }

        /*-----------------------------------------------------------------------------
         *  solve with U
         *-----------------------------------------------------------------------------*/
        for( i = n-1; i >= 0; i-- ) {
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {
                xc[i] -= xc[U->colptr[j]] * U->values[j];
            }
            xc[i] *= D[i];
        }

    }
    return 0;
}

/*-----------------------------------------------------------------------------
 *  Solve with a complex ILU
 *-----------------------------------------------------------------------------*/
static int __iluk_solve_complex( mess_precond myself, mess_solver_options opt, mess_vector in, mess_vector out){
    MSG_FNAME(__func__);
    ilukptr ilu;
    mess_int_t n , i, j;
    mess_double_cpx_t *Dx;
    ccsrptr L, U;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(myself);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);
    mess_check_real_or_complex(out);

    ilu = (ilukptr) myself->data;
    mess_check_nullpointer(ilu);

    if ( MESS_IS_REAL(ilu)) {
        return __iluk_solve(myself, opt, in, out);
    }

    n= ilu->n;
    L = ilu->L;
    U = ilu->U;
    Dx = ilu->Dx;

    if ( n != in->dim) {
        MSG_ERROR("input has the wrong number of rows. is = " MESS_PRINTF_INT ",  wanted = " MESS_PRINTF_INT "\n", in->dim, n);
        return MESS_ERROR_DIMENSION;
    }
    if (n != out->dim) {
        MSG_ERROR("output has the wrong number of rows. is = " MESS_PRINTF_INT ", wanted = " MESS_PRINTF_INT "\n", out->dim, n );
        return MESS_ERROR_DIMENSION;
    }
    if ( MESS_IS_REAL(in) ) {
        ret = mess_vector_tocomplex(out);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        mess_double_cpx_t *x;
        double *y;

        x = out -> values_cpx;
        y = in -> values;


        /*-----------------------------------------------------------------------------
         *  solve with L
         *-----------------------------------------------------------------------------*/
        for( i = 0; i < n; i++ ) {
            x[i] = y[i];
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {
                x[i] -= x[L->colptr[j]] * L->values_cpx[j];
            }
        }

        /*-----------------------------------------------------------------------------
         *  solve with U
         *-----------------------------------------------------------------------------*/
        for( i = n-1; i >= 0; i-- ) {
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {
                x[i] -= x[U->colptr[j]] * U->values_cpx[j];
            }
            x[i] *= Dx[i];
        }
    } else if ( MESS_IS_COMPLEX(in)) {
        mess_double_cpx_t  *xc, *yc;
        ret = mess_vector_tocomplex(out);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        xc = out -> values_cpx;
        yc = in -> values_cpx;


        /*-----------------------------------------------------------------------------
         *  solve with L
         *-----------------------------------------------------------------------------*/
        for( i = 0; i < n; i++ ) {
            xc[i] = yc[i];
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {
                xc[i] -= xc[L->colptr[j]] * L->values_cpx[j];
            }
        }

        /*-----------------------------------------------------------------------------
         *  solve with U
         *-----------------------------------------------------------------------------*/
        for( i = n-1; i >= 0; i-- ) {
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {
                xc[i] -= xc[U->colptr[j]] * U->values_cpx[j];
            }
            xc[i] *= Dx[i];
        }
    }
    return 0;
}

/**
 * @brief Clear a matrix.
 * @param[in,out] mat   matrix to clear
 * @return always zero
 *
 * The @ref clean_ccsrmat functions cleans up a matrix \f$ mat \f$.
 *
 */
int clean_ccsrmat(ccsrptr mat) {
    if (mat == NULL) return 0;
    if (mat->n < 1) return 0;
    if (mat->rowptr)mess_free(mat->rowptr);
    if (mat->colptr)mess_free(mat->colptr);
    if (mat->values)mess_free(mat->values);
    if (mat->values_cpx)mess_free(mat->values_cpx);
    mess_free(mat);
    return 0;
}

static int __iluk_clear(mess_precond myself){
    ilukptr ilu;
    if (myself == NULL) return 0;
    ilu = (ilukptr) myself -> data;
    if (ilu == NULL) return 0;
    if (ilu->D!=NULL)mess_free(ilu->D);
    if (ilu->Dx!=NULL)mess_free(ilu->Dx);

    clean_ccsrmat(ilu->L);
    clean_ccsrmat(ilu->U);
    mess_free(ilu);
    return 0;
}


/*-----------------------------------------------------------------------------
 *  Symbolic Phase
 *-----------------------------------------------------------------------------*/
static int symbolic_iluk( mess_int_t lofM, mess_matrix mat, ilukptr ilu )
{
    MSG_FNAME(__func__);
    mess_int_t  n;
    mess_int_t  *levls = NULL, *jbuf = NULL;
    mess_int_t *iw=NULL;
    mess_int_t  **ulvl;  /*  stores lev-fils for U part of ILU factorization*/
    ccsrptr L,U;
    mess_int_t i, j, k, col, ip, it, jpiv;
    mess_int_t incl, incu, jmin, kmin;



    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_square(mat);

    n = mat -> rows;

    mess_try_alloc(iw, mess_int_t* , sizeof(mess_int_t) * n);
    mess_try_alloc(levls, mess_int_t *, n*sizeof(mess_int_t));
    mess_try_alloc(jbuf, mess_int_t *, n*sizeof(mess_int_t));
    mess_try_alloc(ulvl, mess_int_t**,n*sizeof(mess_int_t*));
    for (i=0; i < n ; i++) ulvl[i]=NULL;
    L = ilu-> L;
    U = ilu-> U;


    /* initilize iw */
    for( j = 0; j < n; j++ ) iw[j] = -1;
    for( i = 0; i < n; i++ ) {
        incl = 0;
        incu = i;
        /*-------------------- assign lof = 0 for matrix elements */
        for( j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++ ) {
            col = mat->colptr[j];
            if( col < i ) {
                /*-------------------- L-part  */
                jbuf[incl] = col;
                levls[incl] = 0;
                iw[col] = incl++;
            }  else if (col > i) {
                /*-------------------- U-part  */
                jbuf[incu] = col;
                levls[incu] = 0;
                iw[col] = incu++;
            }
        }
        /*-------------------- symbolic k,i,j Gaussian elimination  */
        jpiv = -1;
        while (++jpiv < incl) {
            k = jbuf[jpiv] ;
            /*-------------------- select leftmost pivot */
            kmin = k;
            jmin = jpiv;
            for( j = jpiv + 1; j< incl; j++) {
                if( jbuf[j] < kmin ) {
                    kmin = jbuf[j];
                    jmin = j;
                }
            }
            /*-------------------- swap  */
            if( jmin != jpiv ) {
                jbuf[jpiv] = kmin;
                jbuf[jmin] = k;
                iw[kmin] = jpiv;
                iw[k] = jmin;
                j = levls[jpiv] ;
                levls[jpiv] = levls[jmin];
                levls[jmin] = j;
                k = kmin;
            }
            /*-------------------- symbolic linear combinaiton of rows  */
            for( j = U->rowptr[k]; j < U->rowptr[k+1]; j++ ) {
                col = U->colptr[j];
                it = ulvl[k][j-U->rowptr[k]]+levls[jpiv]+1 ;
                if( it > lofM ) continue;
                ip = iw[col];
                if( ip == -1 ) {
                    if( col < i) {
                        jbuf[incl] = col;
                        levls[incl] = it;
                        iw[col] = incl++;
                    } else if( col > i ) {
                        jbuf[incu] = col;
                        levls[incu] = it;
                        iw[col] = incu++;
                    }
                } else
                    levls[ip] = min(levls[ip], it);
            }
        }   /* end - while loop */
        /*-------------------- reset iw */
        for( j = 0; j < incl; j++ ) iw[jbuf[j]] = -1;
        for( j = i; j < incu; j++ ) iw[jbuf[j]] = -1;
        /*-------------------- copy L-part */
        L->rowptr[i+1] = L->rowptr[i]+incl;
        if(incl > 0 ) {
            if (L->nnz < L->rowptr[i+1]) {
                mess_int_t ext = (mat->nnz < mat->rows) ? (mat->rows):(mat->nnz);
                mess_try_realloc (L->colptr, mess_int_t * , sizeof(mess_int_t)*(L->nnz+ext));
                L->nnz += ext;
            }
            memcpy( L->colptr+L->rowptr[i], jbuf, sizeof(mess_int_t)*incl);
        }
        /*-------------------- copy U - part        */
        k = incu-i;
        U->rowptr[i+1] = U->rowptr[i]+k;
        if( k > 0 ) {
            if (U->nnz < U->rowptr[i+1]) {
                mess_int_t ext = (mat->nnz < mat->rows) ? (mat->rows):(mat->nnz);
                mess_try_realloc (U->colptr, mess_int_t * , sizeof(mess_int_t)*(U->nnz+ext));
                U->nnz += ext;
            }
            memcpy(U->colptr+U->rowptr[i], jbuf+i, sizeof(mess_int_t)*k );
            /*-------------------- update matrix of levels */
            mess_try_alloc(ulvl[i], mess_int_t *, sizeof(mess_int_t) *k );
            memcpy( ulvl[i], levls+i, k*sizeof(mess_int_t) );
        }
    }

    mess_free(levls);
    mess_free(jbuf);
    mess_free(iw);
    for(i = 0; i < n-1; i++ ) {
        if (ulvl[i])mess_free(ulvl[i]) ;
    }
    mess_free(ulvl);
    L->nnz = L->rowptr[n];
    U->nnz = U->rowptr[n];
    mess_try_realloc (L->colptr, mess_int_t * , sizeof(mess_int_t)*(L->nnz));
    mess_try_realloc (U->colptr, mess_int_t * , sizeof(mess_int_t)*(U->nnz));
    return 0;
}


/**
 * @brief Setup an \f$ ILU(k) \f$ preconditioner from a matrix.
 * @param[out] pre  generated preconditioner
 * @param[in]  mat   input matrix to compute the preconditioner for
 * @param[in]  lof   input level of fill
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_precond_iluk function computes an \f$ ILU(k) \f$ preconditioner for a
 * real or complex matrix \f$ mat \f$. \n
 * All entries with a level \f$ > lof \f$ are dropped and set to \f$ 0 \f$.
 * If \f$ lof==0 \f$ it results in the classical incomplete LU decomposition.
 *
 * \sa mess_precond
 * \sa mess_solver_bicgstab
 *
 * \remarks
 * \li The quality and the efficiency of the preconditioner can be improved by reordering the matrix before using
 *      the AMD reordering (\ref mess_matrix_reorder_amd).
 * \li All the diagonals of the input matrix must not be zero.
 *
 */
int mess_precond_iluk( mess_precond pre, mess_matrix mat, mess_int_t lof ) {
    MSG_FNAME(__func__);
    int ierr;
    mess_int_t n;
    mess_int_t  *jw, i, j, k, col, jpos, jrow;
    ccsrptr L, U;
    double *D;
    mess_double_cpx_t *Dx;
    ilukptr ilu;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_square(mat);
    mess_check_real_or_complex(mat);

    n = mat->rows;

    mess_try_alloc(ilu, ilukptr, sizeof(ILUK));
    mess_try_alloc(ilu->L, ccsrptr, sizeof(CCSRMat));
    mess_try_alloc(ilu->U, ccsrptr, sizeof(CCSRMat));
    if ( MESS_IS_REAL(mat) ) {
        mess_try_alloc(ilu->D, double*, sizeof(double)*n);
        ilu->Dx = NULL;
        ilu->data_type = MESS_REAL;
    } else {
        mess_try_alloc(ilu->Dx, mess_double_cpx_t *, sizeof(mess_double_cpx_t )*n);
        ilu->D = NULL;
        ilu->data_type = MESS_COMPLEX;
    }

    L = ilu->L;
    U = ilu->U;
    D = ilu->D;
    Dx = ilu->Dx;

    // only alloc the base
    mess_try_alloc(L->rowptr, mess_int_t*, sizeof(mess_int_t)*(n+1));
    mess_try_alloc(L->colptr, mess_int_t*, sizeof(mess_int_t*)*(mat->nnz));
    L->nnz = mat->nnz;
    mess_try_alloc(U->rowptr, mess_int_t*, sizeof(mess_int_t)*(n+1));
    mess_try_alloc(U->colptr, mess_int_t*, sizeof(mess_int_t*)*(mat->nnz));
    U->nnz = mat->nnz;
    for (i=0; i<=n; i++) {
        L->rowptr[i]=U->rowptr[i]=0;
    }

    ilu->n=n;
    L->n=n;
    U->n=n;


    /* symbolic factorization to calculate level of fill index arrays */
    ierr = symbolic_iluk( lof, mat, ilu );
    FUNCTION_FAILURE_HANDLE(ierr, (ierr!=0), symbolic_iluk);

    mess_try_alloc(jw, mess_int_t *, sizeof(mess_int_t) *n );
    if ( MESS_IS_REAL(mat) ) {
        mess_try_alloc(L->values, double *, sizeof(double) * L->nnz);
        mess_try_alloc(U->values, double *, sizeof(double) * U->nnz);
        L->values_cpx = NULL;
        U->values_cpx = NULL;
    } else {
        mess_try_alloc(L->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t ) * L->nnz);
        mess_try_alloc(U->values_cpx, mess_double_cpx_t *, sizeof(mess_double_cpx_t ) * U->nnz);
        L->values = NULL;
        U->values = NULL;
    }

    /* set indicator array jw to -1 */
    for( j = 0; j < n; j++ ) jw[j] = -1;

    /* beginning of main loop */
    if ( MESS_IS_REAL(mat)) {
        for( i = 0; i < n; i++ ) {
            /* set up the i-th row accroding to the nonzero information from
               symbolic factorization */

            /* setup array jw[], and initial i-th row */
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {  /* initialize L part   */
                col = L->colptr[j];
                jw[col] = j;
                L->values[j] = 0;
            }
            jw[i] = i;
            D[i] = 0; /* initialize diagonal */
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {  /* initialize U part   */
                col = U->colptr[j];
                jw[col] = j;
                U->values[j] = 0;
            }

            /* copy row from csmat into lu */
            for( j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++ ) {
                col = mat->colptr[j];
                jpos = jw[col];
                if( col < i )
                    L->values[jpos] = mat->values[j];
                else if( col == i )
                    D[i] = mat->values[j];
                else
                    U->values[jpos] = mat->values[j];
            }

            /* eliminate previous rows */
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {
                jrow = L->colptr[j];
                /* get the multiplier for row to be eliminated (jrow) */
                L->values[j] *= D[jrow];
                /* combine current row and row jrow */
                for( k = U->rowptr[jrow]; k < U->rowptr[jrow+1]; k++ ) {
                    col = U->colptr[k];
                    jpos = jw[col];
                    if( jpos == -1 ) continue;
                    if( col < i )
                        L->values[jpos] -= L->values[j] * U->values[k];
                    else if( col == i )
                        D[i] -= L->values[j] * U->values[k];
                    else
                        U->values[jpos] -= L->values[j] * U->values[k];
                }
            }

            /* reset double-pointer to -1 ( U-part) */
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ){
                col = L->colptr[j];
                jw[col] = -1;
            }
            jw[i] = -1;
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {
                col = U->colptr[j];
                jw[col] = -1;
            }
            if( D[i] == 0 ) {
                /* for( j = i+1; j < n; j++ ) {
                   L->values[j] = NULL;
                   U->values[j] = NULL;
                   } */
                MSG_ERROR("fatal error: Zero diagonal found...\n" );
                return MESS_ERROR_SINGULAR;
            }
            D[i] = 1.0 / D[i];
        }
    } else {
        for( i = 0; i < n; i++ ) {
            /* set up the i-th row accroding to the nonzero information from
               symbolic factorization */

            /* setup array jw[], and initial i-th row */
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {  /* initialize L part   */
                col = L->colptr[j];
                jw[col] = j;
                L->values_cpx[j] = 0;
            }
            jw[i] = i;
            Dx[i] = 0; /* initialize diagonal */
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {  /* initialize U part   */
                col = U->colptr[j];
                jw[col] = j;
                U->values_cpx[j] = 0;
            }

            /* copy row from csmat into lu */
            for( j = mat->rowptr[i]; j < mat->rowptr[i+1]; j++ ) {
                col = mat->colptr[j];
                jpos = jw[col];
                if( col < i )
                    L->values_cpx[jpos] = mat->values_cpx[j];
                else if( col == i )
                    Dx[i] = mat->values_cpx[j];
                else
                    U->values_cpx[jpos] = mat->values_cpx[j];
            }

            /* eliminate previous rows */
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ) {
                jrow = L->colptr[j];
                /* get the multiplier for row to be eliminated (jrow) */
                L->values_cpx[j] *= Dx[jrow];
                /* combine current row and row jrow */
                for( k = U->rowptr[jrow]; k < U->rowptr[jrow+1]; k++ ) {
                    col = U->colptr[k];
                    jpos = jw[col];
                    if( jpos == -1 ) continue;
                    if( col < i )
                        L->values_cpx[jpos] -= L->values_cpx[j] * U->values_cpx[k];
                    else if( col == i )
                        Dx[i] -= L->values_cpx[j] * U->values_cpx[k];
                    else
                        U->values_cpx[jpos] -= L->values_cpx[j] * U->values_cpx[k];
                }
            }

            /* reset double-pointer to -1 ( U-part) */
            for( j = L->rowptr[i]; j < L->rowptr[i+1]; j++ ){
                col = L->colptr[j];
                jw[col] = -1;
            }
            jw[i] = -1;
            for( j = U->rowptr[i]; j < U->rowptr[i+1]; j++ ) {
                col = U->colptr[j];
                jw[col] = -1;
            }
            if( Dx[i] == 0.0 ) {
                MSG_ERROR("fatal error: Zero diagonal found...\n" );
                return MESS_ERROR_SINGULAR;
            }
            Dx[i] = 1.0 / Dx[i];
        }
    }
    mess_free(jw);


    /*-----------------------------------------------------------------------------
     *  finalize precond
     *-----------------------------------------------------------------------------*/
    pre->type = MESS_PRECOND_ILUK;
    pre->data = ( void * ) ilu;
    pre->solve = __iluk_solve;
    pre->clear = __iluk_clear;

    return 0;
}

