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
 * @file lib/eigen/arpack.c
 * @brief Interface to @arpack.
 * @author @koehlerm
 *
 * This file implements an interface to @arpack.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>

/*-----------------------------------------------------------------------------
 *  ARPACK Functions
 *-----------------------------------------------------------------------------*/
void F77_GLOBAL(dnaupd,DNAUPD)( mess_int_t *ido, char * bmat, mess_int_t *n, char * which,
        mess_int_t *nev, double* tol, double *resid, mess_int_t *ncv,
        double *v, mess_int_t *ldv, mess_int_t *iparam, mess_int_t *ipntr,
        double *workd, double *workl, mess_int_t* lworkl, mess_int_t *info );

void F77_GLOBAL(dsaupd,DSAUPD)( mess_int_t *ido, char * bmat, mess_int_t *n, char * which,
        mess_int_t *nev, double* tol, double *resid, mess_int_t *ncv,
        double *v, mess_int_t *ldv, mess_int_t *iparam, mess_int_t *ipntr,
        double *workd, double *workl, mess_int_t* lworkl, mess_int_t *info );

void F77_GLOBAL(dneupd,DNEUPD)(MESS_LAPACK_LOGICAL *rvec , char * howmny, MESS_LAPACK_LOGICAL *select, double *dr , double *di,
        double * z , mess_int_t *ldz , double *sigmar, double *sigmai, double *workev,
        char * bmat, mess_int_t *n, char * which,
        mess_int_t *nev, double* tol, double *resid, mess_int_t *ncv,
        double *v, mess_int_t *ldv, mess_int_t *iparam, mess_int_t *ipntr,
        double *workd, double *workl, mess_int_t* lworkl, mess_int_t *info );

void F77_GLOBAL(dseupd,DSEUPD)(MESS_LAPACK_LOGICAL *rvec , char * howmny, MESS_LAPACK_LOGICAL *select, double *d,
        double * z , mess_int_t *ldz , double *sigma,
        char * bmat, mess_int_t *n, char * which,
        mess_int_t *nev, double* tol, double *resid, mess_int_t *ncv,
        double *v, mess_int_t *ldv, mess_int_t *iparam, mess_int_t *ipntr,
        double *workd, double *workl, mess_int_t* lworkl, mess_int_t *info );

void F77_GLOBAL(znaupd,ZNAUPD)( mess_int_t *ido, char * bmat, mess_int_t *n, char * which,
        mess_int_t *nev, double* tol, mess_double_cpx_t *resid, mess_int_t *ncv,
        mess_double_cpx_t *v, mess_int_t *ldv, mess_int_t *iparam, mess_int_t *ipntr,
        mess_double_cpx_t *workd, mess_double_cpx_t *workl, mess_int_t* lworkl, double *rwork,  mess_int_t *info );

void F77_GLOBAL(zneupd,ZNEUPD)(MESS_LAPACK_LOGICAL *rvec , char * howmny, MESS_LAPACK_LOGICAL *select, mess_double_cpx_t *d,
        mess_double_cpx_t * z , mess_int_t *ldz , mess_double_cpx_t * sigma, mess_double_cpx_t *workev,
        char * bmat, mess_int_t *n, char * which,
        mess_int_t *nev, double* tol, mess_double_cpx_t *resid, mess_int_t *ncv,
        mess_double_cpx_t *v, mess_int_t *ldv, mess_int_t *iparam, mess_int_t *ipntr,
        mess_double_cpx_t *workd, mess_double_cpx_t *workl, mess_int_t* lworkl, double *rwork,  mess_int_t *info );


#define ARPACK_IPARAM_SIZE 11
#define ARPACK_IPNTR_SIZE 14


/*-----------------------------------------------------------------------------
 *  ARPACK Error Codes
 *-----------------------------------------------------------------------------*/
static void arpack_print_error(mess_int_t ierr) {
    switch (ierr){
        case -1:
            fprintf(stderr,"N must be positive.\n");
            break;
        case -2:
            fprintf(stderr,"NEV must be positive.\n");
            break;
        case -3:
            fprintf(stderr,"NCV-NEV >= 2 and less than or equal to N.\n");
            break;
        case -4:
            fprintf(stderr,"The maximum number of Arnoldi update iteration     must be greater than zero.\n");
            break;
        case -5:
            fprintf(stderr,"WHICH has the wrong value.\n");
            break;
        case -6:
            fprintf(stderr,"BMAT must be one of 'I' or 'G'.\n");
            break;
        case -7:
            fprintf(stderr,"Length of private work array is not sufficient.\n");
            break;
        case -8:
            fprintf(stderr,"Error return from LAPACK eigenvalue calculation;\n");
            break;
        case -9:
            fprintf(stderr,"Starting vector is zero or Informational error from LAPACK routine dlahqr.\n");
            break;
        case -10:
            fprintf(stderr,"IPARAM(7) must be 1,2,3,4. (or 5 for xsaupd)\n");
            break;
        case -11:
            fprintf(stderr,"IPARAM(7) = 1 and BMAT = 'G' are incompatable.\n");
            break;
        case -12:
            fprintf(stderr,"IPARAM(1) must be equal to 0 or 1 or HOWMNY = 'S' not yet implemented\n");
            break;
        case -13:
            fprintf(stderr,"HOWMNY must be one of 'A' or 'P' if RVEC = .true. Or NEV and which='BE' are incompatible.\n");
            break;
        case -14:
            fprintf(stderr,"D(N/S)AUPD  did not find any eigenvalues to sufficient accuracy.\n");
            break;
        case -15:
            fprintf(stderr,"DNEUPD got a different count of the number of converged Ritz values than DNAUPD got.  This indicates the user probably made an error in passing data from DNAUPD to DNEUPD or that the data was modified before entering DNEUPD.\n or DNAUPD: HOWMNY must be one of A or S");
            break;
        case -16:
            fprintf(stderr, "DSEUPD: HOWMNY='S' not implemented.\n");
            break;
        case -17:
            fprintf(stderr, "DSEUPD: Data modified between DSAUPD and DSEUPD.\n");
            break;
        case -9999:
            fprintf(stderr,"Could not build an Arnoldi factorization. IPARAM(5) returns the size of the current Arnoldi factorization.\n");
            break;
        default:
            fprintf(stderr,"A unknown error occured.\n");
            break;
    }
}

/*-----------------------------------------------------------------------------
 *  ARPACK IPARAM Print
 *-----------------------------------------------------------------------------*/
/**
 * @brief Print @arpack iparam.
 * @param[in] iparam  input pointer to different shifts
 *
 * The @ref arpack_print_iparam function prints all shifts to stderr.
 */
void arpack_print_iparam(mess_int_t *iparam) {
    fprintf(stderr, "IPARAM(0):  ISHIFT         = %ld\n",(long) iparam[0]);
    fprintf(stderr, "IPARAM(1):  UNUSED         = %ld\n",(long) iparam[1]);
    fprintf(stderr, "IPARAM(2):  MXITER         = %ld\n",(long) iparam[2]);
    fprintf(stderr, "IPARAM(3):  NB             = %ld\n",(long) iparam[3]);
    fprintf(stderr, "IPARAM(4):  NCONV          = %ld\n",(long) iparam[4]);
    fprintf(stderr, "IPARAM(5):  IUPD(UNUSED)   = %ld\n",(long) iparam[5]);
    fprintf(stderr, "IPARAM(6):  MODE           = %ld\n",(long) iparam[6]);
    fprintf(stderr, "IPARAM(7):  NP             = %ld\n",(long) iparam[7]);
    fprintf(stderr, "IPARAM(8):  NUMOP          = %ld\n",(long) iparam[8]);
    fprintf(stderr, "IPARAM(9):  NUMOPB         = %ld\n",(long) iparam[9]);
    fprintf(stderr, "IPARAM(10): NUMREO         = %ld\n",(long) iparam[10]);
}

/*-----------------------------------------------------------------------------
 *  Arnoldi Process for Real operators
 *-----------------------------------------------------------------------------*/
static int arpack_arnoldi_real ( mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V) {
    MSG_FNAME(__func__);
    int ret = 0 ;

    /*-----------------------------------------------------------------------------
     *  Variables
     *-----------------------------------------------------------------------------*/
    // dnaupd
    mess_int_t ido;
    char bmat[]="I";
    mess_int_t n;
    char which[]="LR";
    double tol;
    double *resid;
    mess_int_t ncv;
    double *v;
    mess_int_t ldv;
    mess_int_t iparam[ARPACK_IPARAM_SIZE];
    mess_int_t ipntr[ARPACK_IPNTR_SIZE];
    double *workd;
    double *workl;
    mess_int_t lworkl;
    mess_int_t info;
    mess_int_t i;
    // dneupd
    MESS_LAPACK_LOGICAL rvec;
    char howmny[]="P";
    MESS_LAPACK_LOGICAL *select;
    double *dr;
    double *di;
    double *sigmar;
    double *sigmai;
    double *workev;
    mess_int_t ierr = 0;
    // additional ones.
    mess_vector x, y;
    mess_int_t one = 1, j;

    double *work = NULL;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if ( A->data_type != MESS_REAL) {
        MSG_ERROR("matrix-vector product operator must be real.\n");
        return MESS_ERROR_DATATYPE;
    }
    if (A->dim < 0 ) {
        MSG_ERROR("Operator must at least have dimension 0\n");
        return MESS_ERROR_DIMENSION;
    }
    if (A->dim == 0) {
        return 0;
    }

    // check b0
    if ( opt.b0 != NULL) {
        if ( opt.b0->dim != A->dim ) {
            MSG_ERROR("The given start vector has the wrong dimension. b0->dim = %ld, MVP.dim = %ld\n", (long) opt.b0->dim, (long) A->dim);
            return MESS_ERROR_DIMENSION;
        }
        if ( !MESS_IS_REAL(opt.b0)) {
            MSG_ERROR("The given start vector has the wrong data type. For the real Arnoldiprocess it must be real too.\n");
            return MESS_ERROR_DATATYPE;
        }
    }
    /*-----------------------------------------------------------------------------
     *  Prepare Arnoldi
     *-----------------------------------------------------------------------------*/
    ido = 0;
    bmat[0] = 'I';                  // Type of the eigenvalue problem
    n = A->dim;                     // Dimension of the ep
    switch (opt.which) {            // select the kind of desired eigenvalues
        case MESS_EIGEN_ARPACK_LARGE_MAGNITUDE:
            strncpy(which, "LM",3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_MAGNITUDE:
            strncpy(which, "SM",3);
            break;
        case MESS_EIGEN_ARPACK_LARGE_REALPART:
            strncpy(which, "LR",3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_REALPART:
            strncpy(which, "SR", 3);
            break;
        case MESS_EIGEN_ARPACK_LARGE_IMAGPART:
            strncpy(which, "LI", 3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_IMAGPART:
            strncpy( which, "SI", 3);
            break;
        default:
            MSG_ERROR("opt.which option is not supported by mess_eigen_arpack or dnaupd. The option is may be valid for mess_eigen_arpack_lanczos.\n");
            return MESS_ERROR_ARGUMENTS;
    }
    tol = opt.tol;      // tolerance
    ncv = opt.ncv;      // Size of the Krylov subspace
    if ( ncv >= n ) ncv = n - 1;
    if ( ncv - nev < 2 ) {
        MSG_ERROR("The size of the Krylov subspace must be larger than the number of desired eigenvalues.");
        return MESS_ERROR_ARGUMENTS;
    }

    ldv = n;
    for (i = 0; i < ARPACK_IPARAM_SIZE; i++) {
        iparam[i] = 0;
    }
    iparam[0] = 1;          // Shifts from internal
    iparam[2] = opt.maxit;  // Max IT
    iparam[3] = 1;          // Block size
    iparam[6] = 1;          // Mode


    mess_try_alloc(work, double *, sizeof(double)*(3*n+3*ncv*ncv+6*ncv+n+n*ncv+2*(nev+1)+3*ncv));
    lworkl = 3 * ncv * ncv + 6 * ncv;
    workd = work;  work+= 3*n;
    workl = work; work += lworkl;
    resid = work; work += n;
    v = work; work+=(n*ncv);

    if ( opt.b0 != NULL ) {
        for (i = 0; i < n; i++) {
            resid[i] = opt.b0->values[i];
        }
        info = 1;
    } else {
        info = 0;
    }

    MESS_INIT_VECTORS(&x,&y);
    ret = mess_vector_alloc(x, n, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(y, n, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    /*-----------------------------------------------------------------------------
     *  Start the reverse communication
     *-----------------------------------------------------------------------------*/
    do {
        F77_GLOBAL(dnaupd,DNAUPD)(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam,ipntr,workd, workl, &lworkl, &info);
        if ( ido == -1 || ido == 1) {
            memcpy(x->values, &workd[ipntr[0]-1], n * sizeof(double));
            ret  =  mess_mvpcall_apply(A, MESS_OP_NONE, x, y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
            memcpy(&workd[ipntr[1]-1], y->values, n * sizeof(double));
        }
    } while ( ido == -1 || ido == 1  ) ;
    mess_vector_clear(&x);
    mess_vector_clear(&y);
    if ( info < 0 ) {
        MSG_ERROR("Fatal error in ARPACK DNAUPD. info = %ld\n", (long) info);
        arpack_print_error(info);
        mess_free(workd);
        return MESS_ERROR_MISC;
    }

    /*-----------------------------------------------------------------------------
     *  Post process.
     *-----------------------------------------------------------------------------*/
    if ( V != NULL) {
        rvec = 1;
        strcpy(howmny,"A");
    } else {
        rvec = 0;
        strcpy(howmny,"A");
    }
    mess_try_alloc(select, MESS_LAPACK_LOGICAL *,sizeof(MESS_LAPACK_LOGICAL) * (ncv));
    dr = work; work += (nev+1);
    di = work; work += (nev+1);
    workev = work;
    sigmar = NULL;
    sigmai = NULL;

    F77_GLOBAL(dneupd,DNEUPD)(&rvec, howmny, select, dr, di, v,&ldv, sigmar, sigmai, workev, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &ierr);
    if ( ierr < 0 ) {
        MSG_ERROR("Fatal error in ARPACK DNEUPD. ierr = %ld\n", (long) ierr);
        arpack_print_error(ierr);
        mess_free(select);
        mess_free(workd);
        return MESS_ERROR_MISC;
    }

    /*-----------------------------------------------------------------------------
     *  extract eigenvalues
     *-----------------------------------------------------------------------------*/
    mess_int_t nconv, cpx = 0 ;
    nconv = iparam[4];
    MSG_INFO("Converged Eigenvalues: %ld\n", (long) nconv);
    ret = mess_vector_tocomplex(ev);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);

    ret = mess_vector_resize(ev, nev+1);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    for (i = 0; i < nev; i++) {
        if ( di[i] != 0.0) {
            ev->values_cpx[i] = dr[i] + di[i] *I ;
            ev->values_cpx[i+1] = dr[i] - di[i] *I ;
            i++;
        } else {
            ev->values_cpx[i] = dr[i] + di[i] *I ;
        }
        if ( cpx == 0 && di[i] != 0.0) {
            cpx = 1;
        }
    }
    nev = i;
    ret = mess_vector_resize(ev, nev);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    if (cpx == 0 ) {
        ret = mess_vector_toreal_nowarn(ev);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    }
    /*-----------------------------------------------------------------------------
     *  Copy v to V
     *-----------------------------------------------------------------------------*/
    if ( V != NULL ) {
        ret = mess_matrix_alloc( V, n, nev, n*nev,MESS_DENSE, (cpx == 0) ? MESS_REAL:MESS_COMPLEX);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        if ( cpx == 0 ) {
            for ( i = 0; i < nev; i++) {
                F77_GLOBAL(dcopy,DCOPY)(&n, v+i*n, &one, V->values+i*V->ld,&one);
            }
        } else {
            for (i = 0; i < nev; i++) {
                if ( di[i] == 0.0) {
                    for (j = 0; j < n; j++) {
                        V->values_cpx[j+i*V->ld] = v[i*n+j] + 0.0 * I;
                    }
                } else {
                    for (j = 0; j < n; j++) {
                        V->values_cpx[j+i*V->ld] = v[i*n+j] + v[(i+1)*n+j] * I;
                        V->values_cpx[j+(i+1)*V->ld] = conj(V->values_cpx[j+i*V->ld]) ;
                    }
                    i++;
                }
            }
        }
    }

    if (mess_error_level > 2 ) {
        fprintf(stderr,"====== ARPACK iparam at exit ======\n");
        arpack_print_iparam(iparam);
        fprintf(stderr,"===================================\n");
    }

    mess_free(select);
    mess_free(workd);

    return 0;
}


/*-----------------------------------------------------------------------------
 *  Arnoldi Process for Complex operators.
 *-----------------------------------------------------------------------------*/
static int arpack_arnoldi_complex( mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V) {
    MSG_FNAME(__func__);
    int ret = 0 ;

    /*-----------------------------------------------------------------------------
     *  Variables
     *-----------------------------------------------------------------------------*/
    // dnaupd
    mess_int_t ido;
    char bmat[]="I";
    mess_int_t n;
    char which[]="LR";
    double tol;
    mess_double_cpx_t *resid;
    mess_int_t ncv;
    mess_double_cpx_t *v;
    mess_int_t ldv;
    mess_int_t iparam[ARPACK_IPARAM_SIZE];
    mess_int_t ipntr[ARPACK_IPNTR_SIZE];
    mess_double_cpx_t  *workd;
    mess_double_cpx_t *workl;
    double *rwork;
    mess_int_t lworkl;
    mess_int_t info;
    mess_int_t i;
    // dneupd
    MESS_LAPACK_LOGICAL rvec;
    char howmny[]="P";
    MESS_LAPACK_LOGICAL *select;
    mess_double_cpx_t *d;
    mess_double_cpx_t *sigma;
    mess_double_cpx_t *workev;
    mess_int_t ierr = 0;
    // additional ones.
    mess_vector x, y;

    mess_double_cpx_t *work = NULL;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if ( A->data_type != MESS_COMPLEX) {
        MSG_ERROR("matrix-vector product operator must be real.\n");
        return MESS_ERROR_DATATYPE;
    }
    if (A->dim < 0 ) {
        MSG_ERROR("Operator must at least have dimension 0\n");
        return MESS_ERROR_DIMENSION;
    }
    if (A->dim == 0) {
        return 0;
    }

    // check b0
    if ( opt.b0 != NULL) {
        if ( opt.b0->dim != A->dim ) {
            MSG_ERROR("The given start vector has the wrong dimension. b0->dim = %ld, MVP.dim = %ld\n", (long) opt.b0->dim, (long) A->dim);
            return MESS_ERROR_DIMENSION;
        }
    }
    /*-----------------------------------------------------------------------------
     *  Prepare Arnoldi
     *-----------------------------------------------------------------------------*/
    ido = 0;
    bmat[0] = 'I';          // Type of the eigenvalue problem
    n = A->dim;             // Dimension of the ep
    switch (opt.which) {        // select the kind of desired eigenvalues
        case MESS_EIGEN_ARPACK_LARGE_MAGNITUDE:
            strncpy(which, "LM",3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_MAGNITUDE:
            strncpy(which, "SM",3);
            break;
        case MESS_EIGEN_ARPACK_LARGE_REALPART:
            strncpy(which, "LR",3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_REALPART:
            strncpy(which, "SR", 3);
            break;
        case MESS_EIGEN_ARPACK_LARGE_IMAGPART:
            strncpy(which, "LI", 3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_IMAGPART:
            strncpy( which, "SI", 3);
            break;
        default:
            MSG_ERROR("opt.which option is not supported by mess_eigen_arpack or dnaupd. The option is may be valid for mess_eigen_arpack_lanczos.\n");
            return MESS_ERROR_ARGUMENTS;
    }
    tol = opt.tol;      // tolerance
    ncv = opt.ncv;      // Size of the krylov subspace
    if ( ncv >= n ) ncv = n - 1;
    if ( ncv - nev < 2 ) {
        MSG_ERROR("The size of the Krylov subspace must be larger than the number of desired eigenvalues.");
        return MESS_ERROR_ARGUMENTS;
    }

    ldv = n;
    for (i = 0; i < ARPACK_IPARAM_SIZE; i++) {
        iparam[i] = 0;
    }
    iparam[0] = 1;          // Shifts from internal
    iparam[2] = opt.maxit;      // Max IT
    iparam[3] = 1;          // Block size
    iparam[6] = 1;          // Mode


    mess_try_alloc(work, mess_double_cpx_t *, sizeof(mess_double_cpx_t )*(3*n+3*ncv*ncv+5*ncv+n+n*ncv+(nev+1)+2*ncv));
    lworkl = 3 * ncv * ncv + 5 * ncv;
    workd = work;  work+= 3*n;
    workl = work; work += lworkl;
    resid = work; work += n;
    v = work; work+=(n*ncv);
    mess_try_alloc(rwork, double *, sizeof(double) * n);

    if ( opt.b0 != NULL ) {
        if ( MESS_IS_REAL(opt.b0)) {
            for (i = 0; i < n; i++) {
                resid[i] = opt.b0->values[i];
            }
        } else {
            for (i = 0; i < n; i++) {
                resid[i] = opt.b0->values_cpx[i];
            }
        }
        info = 1;
    } else {
        info = 0;
    }

    MESS_INIT_VECTORS(&x,&y);
    ret = mess_vector_alloc(x, n, MESS_COMPLEX);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(y, n, MESS_COMPLEX);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    /*-----------------------------------------------------------------------------
     *  Start the reverse communication
     *-----------------------------------------------------------------------------*/
    do {
        F77_GLOBAL(znaupd,ZNAUPD)(&ido, bmat, &n, which, &nev, &tol, resid,&ncv,
                v, &ldv, iparam,ipntr,workd, workl, &lworkl,rwork, &info);
        if ( ido == -1 || ido == 1) {
            memcpy(x->values_cpx, &workd[ipntr[0]-1], n * sizeof(mess_double_cpx_t ));
            ret  =  mess_mvpcall_apply(A, MESS_OP_NONE, x, y);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
            memcpy(&workd[ipntr[1]-1], y->values_cpx, n * sizeof(mess_double_cpx_t));
        }
    } while ( ido == -1 || ido == 1  ) ;
    mess_vector_clear(&x);
    mess_vector_clear(&y);
    if ( info < 0 ) {
        MSG_ERROR("Fatal error in ARPACK ZNAUPD. info = %ld\n", (long) info);
        arpack_print_error(info);
        mess_free(workd);
        mess_free(rwork);
        return MESS_ERROR_MISC;
    }

    /*-----------------------------------------------------------------------------
     *  Post process.
     *-----------------------------------------------------------------------------*/
    if ( V != NULL) {
        rvec = 1;
        strcpy(howmny,"A");
    } else {
        rvec = 0;
        strcpy(howmny,"A");
    }
    mess_try_alloc(select, MESS_LAPACK_LOGICAL *,sizeof(MESS_LAPACK_LOGICAL) * (ncv));
    d = work; work += (nev+1);
    workev = work;
    sigma = NULL;

    F77_GLOBAL(zneupd,ZNEUPD)(&rvec, howmny, select, d, v,&ldv, sigma, workev,
            bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

    if ( ierr < 0 ) {
        MSG_ERROR("Fatal error in ARPACK ZNEUPD. ierr = %ld\n", (long) ierr);
        arpack_print_error(ierr);
        mess_free(select);
        mess_free(workd);
        mess_free(rwork);
        return MESS_ERROR_MISC;
    }

    /*-----------------------------------------------------------------------------
     *  extract eigenvalues
     *-----------------------------------------------------------------------------*/
    mess_int_t nconv;
    nconv = iparam[4];
    MSG_INFO("Converged Eigenvalues: %ld\n", (long) nconv);
    nconv = nev ;
    ret = mess_vector_tocomplex(ev);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    ret = mess_vector_resize(ev, nconv);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);

    for (i = 0; i < nconv; i++) {
        ev->values_cpx[i] = d[i];
    }
    /*-----------------------------------------------------------------------------
     *  Copy v to V
     *-----------------------------------------------------------------------------*/
    if ( V != NULL ) {
        ret = mess_matrix_alloc( V, n, nconv, n*nconv,MESS_DENSE, MESS_COMPLEX);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        F77_GLOBAL(zlacpy,ZLACPY)("All", &n,&nconv, v, &ldv, V->values_cpx, &V->ld);
    }

    if (mess_error_level > 2 ) {
        fprintf(stderr,"====== ARPACK iparam at exit ======\n");
        arpack_print_iparam(iparam);
        fprintf(stderr,"===================================\n");
    }

    mess_free(select);
    mess_free(workd);
    mess_free(rwork);
    return 0;
}

/*-----------------------------------------------------------------------------
 *  Lanczos Process for Real operators
 *-----------------------------------------------------------------------------*/
static int arpack_lanczos_real ( mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V) {
    MSG_FNAME(__func__);
    int ret = 0 ;

    /*-----------------------------------------------------------------------------
     *  Variables
     *-----------------------------------------------------------------------------*/
    // dnaupd
    mess_int_t ido;
    char bmat[]="I";
    mess_int_t n;
    char which[]="LR";
    double tol;
    double *resid;
    mess_int_t ncv;
    double *v;
    mess_int_t ldv;
    mess_int_t iparam[ARPACK_IPARAM_SIZE];
    mess_int_t ipntr[ARPACK_IPNTR_SIZE];
    double *workd;
    double *workl;
    mess_int_t lworkl;
    mess_int_t info;
    mess_int_t i;
    // dneupd
    MESS_LAPACK_LOGICAL rvec;
    char howmny[]="P";
    MESS_LAPACK_LOGICAL *select;
    double *d;
    double *sigma;
    mess_int_t ierr = 0;
    // additional ones.
    mess_vector x, y;

    double *work = NULL;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    if ( A->data_type != MESS_REAL) {
        MSG_ERROR("matrix-vector product operator must be real.\n");
        return MESS_ERROR_DATATYPE;
    }
    if (A->dim < 0 ) {
        MSG_ERROR("Operator must at least have dimension 0\n");
        return MESS_ERROR_DIMENSION;
    }
    if (A->dim == 0) {
        return 0;
    }

    // check b0
    if ( opt.b0 != NULL) {
        if ( opt.b0->dim != A->dim ) {
            MSG_ERROR("The given start vector has the wrong dimension. b0->dim = %ld, MVP.dim = %ld\n", (long) opt.b0->dim, (long) A->dim);
            return MESS_ERROR_DIMENSION;
        }
        if ( !MESS_IS_REAL(opt.b0)) {
            MSG_ERROR("The given start vector has the wrong data type. For the real Arnoldiprocess it must be real too.\n");
            return MESS_ERROR_DATATYPE;
        }
    }
    /*-----------------------------------------------------------------------------
     *  Prepare Arnoldi
     *-----------------------------------------------------------------------------*/
    ido = 0;
    bmat[0] = 'I';          // Type of the eigenvalue problem
    n = A->dim;             // Dimension of the ep
    switch (opt.which) {        // select the kind of desired eigenvalues
        case MESS_EIGEN_ARPACK_LARGE_ALGEBRAIC:
            strncpy(which, "LA",3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_ALGEBRAIC:
            strncpy(which, "SA",3);
            break;
        case MESS_EIGEN_ARPACK_LARGE_MAGNITUDE:
            strncpy(which, "LM",3);
            break;
        case MESS_EIGEN_ARPACK_SMALL_MAGNITUDE:
            strncpy(which, "SM", 3);
            break;
        case MESS_EIGEN_ARPACK_BOTH_ENDS:
            strncpy(which, "BE", 3);
            break;
        default:
            MSG_ERROR("opt.which option is not supported by mess_eigen_arpack_lanczos or dsaupd. The option is may be valid for mess_eigen_arpack.\n");
            return MESS_ERROR_ARGUMENTS;
    }
    tol = opt.tol;          // tolerance
    ncv = opt.ncv;          // Size of the krylov subspace
    if ( ncv >= n ) ncv = n;
    if ( ncv - nev < 1 ) {
        MSG_ERROR("The size of the Krylov subspace must be larger than the number of desired eigenvalues.");
        return MESS_ERROR_ARGUMENTS;
    }

    ldv = n;
    for (i = 0; i < ARPACK_IPARAM_SIZE; i++) {
        iparam[i] = 0;
    }
    iparam[0] = 1;              // Shifts from internal
    iparam[2] = opt.maxit;      // Max IT
    iparam[3] = 1;              // Block size
    iparam[6] = 1;              // Mode


    mess_try_alloc(work, double *, sizeof(double)*(n+n*ncv+3*n+ncv*ncv+8*ncv+nev));
    lworkl = ncv*ncv+8*ncv;
    workd = work;  work+= 3*n;
    workl = work; work += lworkl;
    resid = work; work += n;
    v = work; work+=(n*ncv);

    if ( opt.b0 != NULL ) {
        for (i = 0; i < n; i++) {
            resid[i] = opt.b0->values[i];
        }
        info = 1;
    } else {
        info = 0;
    }

    MESS_INIT_VECTORS(&x,&y);
    ret = mess_vector_alloc(x, n, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(y, n, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    /*-----------------------------------------------------------------------------
     *  Start the reverse communication
     *-----------------------------------------------------------------------------*/
    do {
        F77_GLOBAL(dsaupd,DSAUPD)(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam,ipntr,workd, workl, &lworkl, &info);
        if ( ido == -1 || ido == 1) {
            memcpy(x->values, &workd[ipntr[0]-1], n * sizeof(double));
            ret  =  mess_mvpcall_apply(A, MESS_OP_NONE, x, y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
            memcpy(&workd[ipntr[1]-1], y->values, n * sizeof(double));
        }
    } while ( ido == -1 || ido == 1  ) ;
    mess_vector_clear(&x);
    mess_vector_clear(&y);
    if ( info < 0 ) {
        MSG_ERROR("Fatal error in ARPACK DSAUPD. info = %ld\n", (long) info);
        arpack_print_error(info);
        mess_free(workd);
        return MESS_ERROR_MISC;
    }

    /*-----------------------------------------------------------------------------
     *  Post process.
     *-----------------------------------------------------------------------------*/
    if ( V != NULL) {
        rvec = 1;
        strcpy(howmny,"A");
    } else {
        rvec = 0;
        strcpy(howmny,"A");
    }
    mess_try_alloc(select, MESS_LAPACK_LOGICAL *,sizeof(MESS_LAPACK_LOGICAL) * (ncv));
    d = work; work += (nev+1);
    sigma = NULL;

    F77_GLOBAL(dseupd,DSEUPD)(&rvec, howmny, select, d, v,&ldv, sigma,
            bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv,
            iparam, ipntr, workd, workl, &lworkl, &ierr);
    if ( ierr < 0 ) {
        MSG_ERROR("Fatal error in ARPACK DSEUPD. ierr = %ld\n", (long) ierr);
        arpack_print_error(ierr);
        mess_free(select);
        mess_free(workd);
        return MESS_ERROR_MISC;
    }

    /*-----------------------------------------------------------------------------
     *  extract eigenvalues
     *-----------------------------------------------------------------------------*/
    mess_int_t nconv;
    nconv = iparam[4];
    MSG_INFO("Converged Eigenvalues: %ld\n", (long) nconv);
    ret = mess_vector_toreal_nowarn(ev);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_resize(ev, nev);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    for (i = 0; i < nev; i++) {
        ev->values[i]=d[i];
    }
    /*-----------------------------------------------------------------------------
     *  Copy v to V
     *-----------------------------------------------------------------------------*/
    if ( V != NULL ) {
        ret = mess_matrix_alloc( V, n, nev, n*nev,MESS_DENSE, MESS_REAL);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        F77_GLOBAL(dlacpy,DLACPY)("All",&n,&nev, v, &ldv, V->values, &V->ld);
    }

    if (mess_error_level > 2 ) {
        fprintf(stderr,"====== LANCZOS iparam at exit ======\n");
        arpack_print_iparam(iparam);
        fprintf(stderr,"====================================\n");
    }

    mess_free(select);
    mess_free(workd);

    return 0;
}



/**
 * @brief An interface to the Arnoldi process from @arpack.
 * @param[in] A  input matrix to compute its eigenvalues from
 * @param[in] nev  input number of desired eigenvalues
 * @param[in] opt  input options structure to configure the Arnoldi process
 * @param[out] ev  vector containing the computed eigenvalues
 * @param[out] V   matrix containing the computed ritz vectors corresponding to the eigenvalues \n
 * (if not desired set it to \c @c NULL)
 * @return zero on success or a non zero error code
 *
 * The @ref mess_eigen_arpack function calls @arpack \cite LehS96 to solve a large scale eigenvalue problem. \n
 * The internal Arnoldi process is configured using the opt structure. Please adjust the number of the
 * Arnoldi vectors in the \ref mess_eigen_arpack_options_t structure such that it fits to the number of
 * desired eigenvalues. \n
 * The number \c ncv of the Arnoldi vectors and the number of desired eigenvalues must fulfill
 * \f[ ncv - nev \leq 2. \f]
 * If only an operator is given instead of a matrix please use \ref mess_eigen_arpack_template.
 *
 * \see mess_eigen_arpack_template
 * \see mess_eigen_arpack_options_t
 * \see mess_eigen_arpack_which_t
 *
 * \remarks The function does only exists if \c MESS_HAVE_ARPACK is set and @mess is compiled with @arpack support.
 *
 */
int mess_eigen_arpack( mess_matrix A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V)
{
    MSG_FNAME(__func__);
    mess_mvpcall mvpcall;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_real_or_complex(A);
    mess_check_nullpointer(ev);
    mess_check_positive(nev);
    mess_check_positive(opt.maxit);
    mess_check_nonnegative(opt.tol);
    mess_check_positive(opt.ncv);


    /*-----------------------------------------------------------------------------
     *  Setup
     *-----------------------------------------------------------------------------*/
    ret = mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, A); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_matrix);

    if (MESS_IS_REAL(A)){
        ret = arpack_arnoldi_real(mvpcall, nev, opt, ev, V);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), arpack_arnoldi_real );
    } else {
        ret = arpack_arnoldi_complex(mvpcall, nev, opt, ev, V);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), arpack_arnoldi_complex);
    }
    mess_mvpcall_clear(&mvpcall);
    return 0;
}       /* -----  end of function mess_eigen_arpack  ----- */


/**
 * @brief An interface to the Arnoldi process from @arpack (template version).
 * @param[in] A  input matrix-vector product object defining the operator
 * @param[in] nev  input number of desired eigenvalues
 * @param[in] opt  input options structure to configure the Arnoldi process
 * @param[out] ev  vector containing the computed eigenvalues
 * @param[out] V   matrix containing the computed ritz vectors corresponding to the eigenvalues\n
 * (if not desired set it to \c NULL)
 * @return zero on success or a non zero error code
 *
 * The @ref mess_eigen_arpack_template function calls @arpack \cite LehS96 to solve a large scale eigenvalue problem. \n
 * The internal Arnoldi process is configured using the @p opt structure. Please adjust the number of the Arnoldi
 * vectors in the \ref mess_eigen_arpack_options_t structure such that it fits to the number of desired eigenvalues.\n
 * The number \c ncv of the Arnoldi vectors and the number of desired eigenvalues  must fulfill
 * \f[ ncv - nev \leq 2. \f]
 * If a matrix given instead of an operator the \ref mess_eigen_arpack function can be used as well.
 *
 * \see mess_eigen_arpack
 * \see mess_eigen_arpack_options_t
 * \see mess_eigen_arpack_which_t
 *
 * \remarks The function does only exists if \c MESS_HAVE_ARPACK is set and @mess is compiled with @arpack support.
 *
 */
int mess_eigen_arpack_template( mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V)
{
    MSG_FNAME(__func__);
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(ev);
    mess_check_positive(nev);
    mess_check_positive(opt.maxit);
    mess_check_nonnegative(opt.tol);
    mess_check_positive(opt.ncv);
    if ( A->dim < 0 ) {
        MSG_ERROR("The matrix-vector product must be at least of dimension 0. \n");
        return MESS_ERROR_DIMENSION;
    }
    if ( A->dim == 0) {
        ret = mess_vector_zeros(ev);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
        return 0;
    }

    /*-----------------------------------------------------------------------------
     *  Call the arnoldi process
     *-----------------------------------------------------------------------------*/
    if (A->data_type == MESS_REAL){
        ret = arpack_arnoldi_real(A, nev, opt, ev, V);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), arpack_arnoldi_real);
    } else if (A->data_type == MESS_COMPLEX) {
        ret = arpack_arnoldi_complex(A, nev, opt, ev, V);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), arpack_arnoldi_complex);
    }
    else {
        MSG_ERROR("The matrix-vector product must define a real or a complex operator. ");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
}       /* -----  end of function mess_eigen_arpack  ----- */

/**
 * @brief An interface to the Lanczos process from @arpack.
 * @param[in] A  input matrix to compute its eigenvalues from
 * @param[in] nev  input number of desired eigenvalues
 * @param[in] opt  input options structure to configure the Lanczos process
 * @param[out] ev  vector containing the computed eigenvalues
 * @param[out] V   matrix containing the computed ritz vectors corresponding to the eigenvalues \n
 * (if not desired set it to \c NULL)
 * @return zero on success or a non zero error code
 *
 * The @ref mess_eigen_arpack_lanczos function calls @arpack \cite LehS96 to solve a large scale eigenvalue problem.\n
 * The internal Lanczos process is configured using the opt structure. Please adjust the number of the Lanczos
 * vectors in the \ref mess_eigen_arpack_options_t structure such that it fits to the number of desired eigenvalues. \n
 * The number \c ncv of the Lanczos vectors and the number of desired eigenvalues must fulfill
 * \f[ ncv - nev \leq 1. \f]
 * If only an operator is given instead of a matrix please use \ref mess_eigen_arpack_lanczos_template.
 *
 * \attention
 * The input matrix must be real symmetric, i.e. \f$ A=A^T \f$. \n
 * If the matrix is not symmetric please use \ref mess_eigen_arpack instead.
 *
 * \see mess_eigen_arpack_lanczos_template
 * \see mess_eigen_arpack_options_t
 * \see mess_eigen_arpack_which_t
 *
 * \remarks The function does only exists if \c MESS_HAVE_ARPACK is set and @mess is compiled with @arpack support.
 *
 */
int mess_eigen_arpack_lanczos( mess_matrix A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V)
{
    MSG_FNAME(__func__);
    mess_mvpcall mvpcall;
    int ret = 0;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_square(A);
    mess_check_real(A);
    mess_check_nullpointer(ev);
    mess_check_positive(nev);
    mess_check_positive(opt.maxit);
    mess_check_nonnegative(opt.tol);
    mess_check_positive(opt.ncv);


    /*-----------------------------------------------------------------------------
     *  Setup
     *-----------------------------------------------------------------------------*/
    ret = mess_mvpcall_matrix(&mvpcall, MESS_OP_NONE, A);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_matrix);
    ret = arpack_lanczos_real(mvpcall, nev, opt, ev, V);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), arpack_lanczos_real);
    mess_mvpcall_clear(&mvpcall);
    return 0;
}       /* -----  end of function mess_eigen_arpack  ----- */

/**
 * @brief An interface to the Lanczos process from @arpack (template version).
 * @param[in] A  input matrix-vector product object defining the operator
 * @param[in] nev  input number of desired eigenvalues
 * @param[in] opt  input options structure to configure the Arnoldi process
 * @param[out] ev  vector containing the computed eigenvalues
 * @param[out] V   matrix containing the computed ritz vectors corresponding to the eigenvalues \n
 * (if not desired set it to \c NULL)
 * @return zero on success or a non zero error code
 *
 * The @ref mess_eigen_arpack_lanczos_template function calls @arpack \cite LehS96 to solve a large scale eigenvalue problem.\n
 * The internal Lanczos process is configured using the opt structure.
 * Please adjust the number of the Lanczos  vectors in the \ref mess_eigen_arpack_options_t structure such that it
 * fits to the number of desired eigenvalues.\n
 * The number \c ncv of the Lanczos vectors and the number of desired eigenvalues must fulfill
 * \f[ ncv - nev \geq 1 . \f]
 * If a matrix given instead of an operator the \ref mess_eigen_arpack function can be used as well.
 *
 * \attention
 * The input matrix must be real symmetric, i.e. \f$ A=A^T \f$. \n
 * If the matrix is not symmetric please use \ref mess_eigen_arpack instead.
 *
 * \see mess_eigen_arpack_lanczos
 * \see mess_eigen_arpack_options_t
 * \see mess_eigen_arpack_which_t
 *
 * \remarks The function does only exists if \c MESS_HAVE_ARPACK is set and @mess is compiled with @arpack support.
 *
 */
int mess_eigen_arpack_lanczos_template( mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V)
{
    MSG_FNAME(__func__);
    int ret;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(ev);
    mess_check_positive(nev);
    mess_check_positive(opt.maxit);
    mess_check_nonnegative(opt.tol);
    mess_check_positive(opt.ncv);
    if ( A->dim < 0 ) {
        MSG_ERROR("The matrix-vector product must be at least of dimension 0. \n");
        return MESS_ERROR_DIMENSION;
    }
    if ( A->dim == 0) {
        ret = mess_vector_zeros(ev);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_zeros);
        return 0;
    }

    /*-----------------------------------------------------------------------------
     *  Call the arnoldi process
     *-----------------------------------------------------------------------------*/
    if (A->data_type == MESS_REAL){
        ret = arpack_lanczos_real(A, nev, opt, ev, V);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), arpack_lanczos_real );
    } else {
        MSG_ERROR("The matrix-vector product must define a real operator. ");
        return MESS_ERROR_DATATYPE;
    }

    return 0;
}       /* -----  end of function mess_eigen_arpack  ----- */

