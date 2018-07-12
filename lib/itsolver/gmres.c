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
 * @file lib/itsolver/gmres.c
 * @brief Solve a linear equation with GMRES.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


#ifdef _OPENMP_H
#include <omp.h>
#endif

// #define MESS_DEBUG_GMRES


/**
 * @brief Solve a system of linear equations with (preconditioned) generalized minimal residual algorithm.
 * @param[in] matrix    input maxtrix \f$ A \f$ of the equation
 * @param[in] pre       input left preconditioner \n (if not used, set to \c NULL)
 * @param[in] b input        right hand side
 * @param[in,out] x    on input: initial guess, \n on output: approximated solution
 * @param[in] opt input object containing the parameters for the iteration, \n  if it is set to \c NULL defaults are used
 * @param[out] stat output statitics about the convergence \n  (if not used, set to \c NULL)
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_solver_gmres function solves a system of linear equations with the (preconditioned) generalized minimal
 * residual (GMRES) algorithm. \n
 * The GMRES-algorithm is a Krylov subspace method to solve linear systems with general
 * matrices. \n
 * For symmetric or skew-symmetric matrices special matrix-vector products are used. \n
 * An initial solution need to be provided in the vector \f$ x \f$ which is used as \f$ x_0 \f$. The relative
 * tolerance provided by \c tol and the dimension of the Krylov subspace provided by \c maxit are
 * parameters for the algorithm. If the relative tolerance is reached before the complete
 * subspace is spanned the algorihm terminates with a hint that no restart is needed.
 * The function online provides one GMRES run. \n
 * To implement a complete restarted GMRES you can use the restart component of the CSS_Solver_Status structure
 * to determine if you call the function with the solution of the last run as initial solution again. \n
 * A short version could be:
 * @code{.c}
 *  int restart = 0;
 *  do {
 *  mess_solver_gmres(matrix, NULL, b, x, opt, stat);
 *  restart = stat->need_restart;
 *  } while ( restart == 1);
 * @endcode
 * The second parameter \c pre is the preconditioner. If pre is not null, the preconditioner is used.
 */
int mess_solver_gmres(  mess_mvpcall matrix,
        mess_precond pre,
        mess_vector b,
        mess_vector x,
        mess_solver_options opt,
        mess_solver_status stat)
{
    MSG_FNAME(__func__);
    int dim;
    mess_int_t i,j, ii, k;

    double **H;
    double *C;
    double *S;
    double *Gamma;
    double *Alpha;
    mess_vector *V;
    mess_vector tmp;
    mess_vector r;
    mess_vector w;
    double tolb;
    double t;
    double normr;
    double fak;
    double hi, hip;
    double beta;
    double normb;
    double tol;
    mess_int_t maxit;
    mess_int_t mj=0;
    int precond = 0;
    int stepdebug = 0;
    int stats = 1;
    int ret = 0;
    double res = 1;


    /*-----------------------------------------------------------------------------
     *  check input parameter
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_nullpointer(x);
    mess_check_nullpointer(b);
    mess_check_real(matrix);
    mess_check_positive(matrix->dim);

    if ( pre != NULL ){
        if ( pre->solve == NULL) {
            MSG_ERROR("missing precondition solve.\n");
            return MESS_ERROR_NULLPOINTER;
        }
        MSG_INFO("preconditioner used\n");
        precond = 1;
    }
    if ( opt == NULL ) {
        MSG_WARN("no options specified, using defaults\n");
        tol = 1e-6;
        maxit = MESS_MIN(150,matrix->dim);
        stepdebug = 0;
    } else {
        tol = opt->tol;
        maxit = opt->maxit;
        if ( opt->stepdebug != NULL ) stepdebug = 1;
    }
    if ( stat == NULL ) {
        MSG_WARN("no statistics avaialable\n");
        stats = 0;
    }

    if ( x->dim != matrix->dim) {
        MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, matrix->dim);
        ret = mess_vector_resize(x, matrix->dim);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    }
    if ( b->dim != matrix->dim) {
        MSG_ERROR("dimension of b mismatch. b->dim = " MESS_PRINTF_INT " \t matrix->dim = " MESS_PRINTF_INT "\n", b->dim, matrix->dim);
        return MESS_ERROR_DIMENSION;
    }

    ret = mess_vector_toreal(x);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
    ret = mess_vector_toreal(b);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);


    dim = matrix->dim;
    if ( maxit > dim ) {
        maxit = dim;
        MSG_WARN("maxit is bigger than dim. maxit set to " MESS_PRINTF_INT "\n", maxit);
    }

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    mess_try_alloc(H , double **, sizeof(double *) * (maxit+1));
    for ( i = 0; i < maxit+1;i++){
        mess_try_alloc(H[i], double *, sizeof(double) * (maxit+1));
    }
    mess_try_alloc( V, mess_vector *, sizeof(mess_vector) * (maxit+1));
    for ( i = 0; i < maxit+1;i++){
        ret = mess_vector_init(&V[i]);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_alloc(V[i], dim, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret , (ret!=0), mess_vector_alloc);
    }

    mess_try_alloc(C, double *,sizeof(double) * (maxit+1));
    mess_try_alloc(S, double *,sizeof(double) * (maxit+1));
    mess_try_alloc(Gamma , double *, sizeof(double) * (maxit+1));
    mess_try_alloc(Alpha , double *, sizeof(double) * (maxit+1));

    MESS_INIT_VECTORS(&tmp,&r,&w);
    ret = mess_vector_alloc(tmp, dim, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(r, dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc(w, dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);


    /*-----------------------------------------------------------------------------
     *  start checks
     *-----------------------------------------------------------------------------*/
    ret = mess_vector_norm2(b,&normb);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    if ( normb == 0.0) normb = 1.0;
    // MSG_PRINT("tol: %lg \t normb = %lg \t tolb=%lg\n", tol, normb, tol*normb );
    tolb = tol* normb;


    /*-----------------------------------------------------------------------------
     *  apply precondtioner
     *-----------------------------------------------------------------------------*/
    // r =b-Ax
    ret = mess_mvpcall_apply(matrix, MESS_OP_NONE,x, tmp);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply);
    if ( stats) stat->num_mvp = 1;
    ret = mess_vector_copy(b, r);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    fak = -1.0;
    ret = mess_vector_axpy( fak, tmp, r);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy );
    ret = mess_vector_norm2( r,&normr);                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
    if (stats) stat->num_mvp = 1;

    if ( normr/normb < tol ){
        MSG_INFO("return before\n");
        if (stats == 1) {
            stat->it = 0;
            stat->res = normr;
            stat->relres = normr/normb;
            stat->need_restart = ( normr/normb < tol) ? 0 :1;
            stat->num_mvp = 1;
            if( stat->need_restart == 0) stat->converged = 1;
        }
        return 0;
    }
    res = normr;

    if ( precond ) {
        // tmp = r;
        double logdiff = 0;
        ret = mess_vector_copy ( r, tmp);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        // M*r = tmp
        ret =   pre->solve(pre,opt, tmp, r);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), pre->solve);
        ret = mess_vector_norm2( r,&normr);          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);
        logdiff =  (fabs(log(normr)/log(10)-log(res)/log(10)));
        // MSG_PRINT("old tol = %lg\n", tol);
        // MSG_PRINT("normr=%lg \t res = %lg log diff = %lg\n", normr, res, logdiff );
        tol *= pow(10, -logdiff );
        MSG_PRINT("new tol = %lg\n", tol);
    }

    // MSG_INFO("residuum at start: %10e  relres: %lg\n", normr, normr/normb);
    mj = 0;
    // V(:,1) = r / normr;
    mess_vector_copy(r, V[0]);
    mess_vector_scale(1/normr, V[0]);
    // Gamma(1) = normr;
    Gamma[0] = normr;

    for (j = 0; j < maxit; j++){
        // if ( normr < tol*startres) {
        if ( res/normb < tol) {
            break;
        }

        mj++;

        /*-----------------------------------------------------------------------------
         *  matrix-vector product
         *-----------------------------------------------------------------------------*/
        if ( precond ) {
            ret = mess_mvpcall_apply(matrix,MESS_OP_NONE, V[j], tmp);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
            ret  =  pre->solve(pre, opt, tmp, w);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), pre->solve);
        }else{
            ret = mess_mvpcall_apply(matrix,MESS_OP_NONE,V[j], w);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
        }
        if (stats) stat->num_mvp ++;

        for (i = 0; i <= j; i++){
            //H(i,j) = V(:,i)'*w;
            ret = mess_vector_dot(V[i],w,&(H[i][j]));                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_dot);
            // w = w - H(i,j)*V(:,i);
            fak = - H[i][j];
            mess_vector_axpy(fak, V[i], w);
        }

        // H(j+1,j) = norm(w);
        ret = mess_vector_norm2(w,&H[j+1][j]);       //FUNCTION_FAILURE_HANLDE(ret,(ret!=0),mess_vector_norm2);

        for ( i = 0; i < j; i++){
            // hi = c(i)*H(i,j)+s(i)*H(i+1, j);
            hi = C[i] * H[i][j] +S[i]*H[i+1][j];
            // hip = -s(i)*H(i,j)+c(i)*H(i+1,j);
            hip = -S[i]*H[i][j] + C[i]*H[i+1][j];
            H[i][j] = hi;
            H[i+1][j] = hip;
        }
        // beta = sqrt(H(j,j)^2 +H(j+1,j)^2);
        beta = sqrt(H[j][j]*H[j][j]+H[j+1][j]*H[j+1][j]);

        if ( beta !=0.0 ){
            // V(:,j+1) = w/H(j+1,j);
            ret = mess_vector_copy(w, V[j+1]);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
            ret = mess_vector_scale(1.0/H[j+1][j], V[j+1]);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);
            // s(j) = H(j+1,j)/beta;
            S[j] = H[j+1][j]/beta;
            // c(j) = H(j,j)/beta;
            C[j] = H[j][j]/beta;
            // H(j,j)=beta;
            H[j][j]=beta;
            // Gamma(j+1)=-s(j)*Gamma(j);
            Gamma[j+1]=-S[j]*Gamma[j];
            // Gamma(j) = c(j)*Gamma(j);
            Gamma[j] = C[j]*Gamma[j];
        }
        normr = fabs(Gamma[j+1]);
        // res = normr;
        MSG_INFO("it = " MESS_PRINTF_INT " \t res = %lg\r", j, res/normb );
        res = normr;
        // Stepdebug...
        if (stepdebug ){
            opt->stepdebug(normr, normr/normb, j, opt->aux_stepdebug);
        }
    }


    /*-----------------------------------------------------------------------------
     *  postprocessing
     *-----------------------------------------------------------------------------*/
    for (ii = 0; ii < mj; ii++){
        i = mj-ii-1;
        t = 0.0;
        for (k=i+1; k < mj; k++){
            t = t + H[i][k]*Alpha[k];
        }
        Alpha[i] = (Gamma[i]-t)/H[i][i];
    }

    for (i=0; i < mj; i++){
        mess_vector_axpy(Alpha[i], V[i], x);
    }

    // r =b-Ax

    /*-----------------------------------------------------------------------------
     *  final statistics
     *-----------------------------------------------------------------------------*/
    // r =b-Ax
    // ret = mess_matrix_mvp(MESS_OP_NONE,matrix, x, tmp);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_mvpcall_apply(matrix, MESS_OP_NONE, x, tmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_apply );
    ret = mess_vector_copy(b, r);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    fak = -1.0;
    mess_vector_axpy( fak, tmp, r);
    ret = mess_vector_norm2( r,&normr);      // FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_norm2);

    if ( stats == 1) {
        stat->it = mj;
        stat->res = normr;
        stat->relres = normr/normb;
        stat->need_restart = ( normr < tolb) ? 0 :1;
        if( stat->need_restart == 0) stat->converged = 1;
        else stat->converged=0;
    }

    for ( i = 0; i < maxit+1; i++){
        mess_free(H[i]);
        mess_vector_clear(&V[i]);
    }
    mess_free(H);
    mess_free(V);
    mess_free(C);
    mess_free(S);
    mess_free(Gamma);
    mess_free(Alpha);
    mess_vector_clear(&tmp);
    mess_vector_clear(&r);
    mess_vector_clear(&w);

    return 0;
    }


    /**
     * @brief Solve a system of linear equations with restarted GMRES.
     * @param[in] matrix input matrix \f$ A \f$ of the equation
     * @param[in] pre input input   preconditioner
     * @param[in] b      input right hand side
     * @param[out] x    solution
     * @param[in] opt    input options
     * @param[out] stat output status
     * @return zero on success or a non-zero error value otherwise
     *
     * The @ref mess_solver_gmres_restart function solves
     * \f[Ax=b \f]
     * with restarted GMRES.
     *
     */
    int mess_solver_gmres_restart(mess_mvpcall matrix,
            mess_precond pre,
            mess_vector b,
            mess_vector x,
            mess_solver_options opt,
            mess_solver_status stat)
    {
        MSG_FNAME(__func__);
        mess_int_t r = -1;
        int err = 0, ret = 0;
        int restart = 0 ;
        double lastrtol = 0;
        double sttol;
        mess_int_t stcount = 0;
        // mess_precond intpre;

        /*-----------------------------------------------------------------------------
         *  check input
         *-----------------------------------------------------------------------------*/
        mess_check_nullpointer(matrix);
        mess_check_nullpointer(b);
        mess_check_nullpointer(x);
        mess_check_nullpointer(opt);
        mess_check_nullpointer(stat);
        mess_check_real(matrix);

        if ( x->dim != matrix->dim) {
            MSG_WARN("resize x from " MESS_PRINTF_INT " to " MESS_PRINTF_INT "\n", x->dim, matrix->dim);
            ret = mess_vector_resize(x, matrix->dim);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
        }

        if ( b->dim != matrix->dim) {
            MSG_ERROR("dimension of b mismatch. b->dim = " MESS_PRINTF_INT " \t matrix->cols = " MESS_PRINTF_INT "\n", b->dim, matrix->dim);
            return MESS_ERROR_DIMENSION;
        }
        ret = mess_vector_toreal(x);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
        ret = mess_vector_toreal(b);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal);
        // intpre = pre;


        /*-----------------------------------------------------------------------------
         *  main loop
         *-----------------------------------------------------------------------------*/
        sttol = mess_eps()*20.0;

        do {
            err=mess_solver_gmres(matrix, pre, b, x, opt, stat);
            FUNCTION_FAILURE_HANDLE(err, (err!=0), mess_solver_gmres);

            if ( (stat->relres >  opt->tol)  && (r< opt->restarts)){
                restart = 1;
            }
            else {
                restart = 0;
            }

            if ( fabs(lastrtol-stat->relres) < sttol) {
                stcount++;
            } else {
                stcount = 0;
            }

            lastrtol = stat->relres;

            if ( stcount >=5 ){
                MSG_WARN("stagnation of GMRES. cancel.\n");
                restart = 0;
                break;
            }

            if ( stat->need_restart !=0 && r < opt->restarts) {
                restart = 1;
            } else {
                restart = 0;
            }
            MSG_INFO("\nrestart = %d, lastrtol = %lg, stcount = " MESS_PRINTF_INT "\n", restart , lastrtol, stcount);
            r++;
        } while (/* stat->need_restart && r < opt->restarts */ restart);

        stat->restarts = r;
        if (stat->converged == 0){
            MSG_ERROR("not converged\n");
        }

        return 0;
    }

