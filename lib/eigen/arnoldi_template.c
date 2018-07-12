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
 * @file lib/eigen/arnoldi_template.c
 * @brief Arnoldi iteration.
 * @author @koehlerm
 * @author @mbehr
 *
 * This file contains the template for the arnoldi process.
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "mess/config.h"
#include <complex.h>
#include <math.h>



#define ARNOLDI_BREAKDOWN_TOL mess_eps()
#define ARNOLDI_RESIZE_CHUNK  20

/**
 * @internal
 *
 * Internal structure holds data for Arnoldi process.
 * @attention Internal use only.
 */
struct mess_krylov_arnoldi_st{
    mess_int_t j;                       /** Number of Arnoldi steps taken.*/
    mess_int_t dim;                     /** Order of the operator used.*/
    mess_datatype_t data_type;          /** Data type of the operator used.*/
    mess_int_t breakdown;               /** If nonzero (early) breakdown occured during Arnoldi process.*/
    double breakdown_tol;               /** Tolerance for breakdown.*/
    mess_matrix V;                      /** Arnoldi vectors as columns.*/
    mess_matrix H;                      /** Upper Hessnberg Matrix.*/
    mess_vector _w;                     /** Cache.*/
    mess_vector _wtmp;                  /** Cache.*/
};

/**
 * @internal
 *
 * Type definition for Arnoldi process.
 * @attention Internal use only.
 */
typedef struct mess_krylov_arnoldi_st * mess_krylov_arnoldi;


/**
 * @brief Initialize a @ref mess_krylov_arnoldi structure.
 * @param[in,out] arn   pointer to a @ref mess_krylov_arnoldi.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_init function initializes a \ref mess_krylov_arnoldi object.
 *
 * @ref mess_krylov_arnoldi_clear
 * @ref mess_krylov_arnoldi_create
 *
 */
int mess_krylov_arnoldi_init(mess_krylov_arnoldi *arn){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_try_alloc(*arn, mess_krylov_arnoldi, sizeof(struct mess_krylov_arnoldi_st));
    memset(*arn,0,sizeof(struct mess_krylov_arnoldi_st));
    return ret;
}


/**
 * @brief Clear a @ref mess_krylov_arnoldi object.
 * @param[in,out] arn   pointer to a @ref mess_krylov_arnoldi object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_clear function cleans up a @ref mess_krylov_arnoldi object.
 *
 * @ref mess_krylov_arnoldi_init
 * @ref mess_krylov_arnoldi_create
 *
 */
int mess_krylov_arnoldi_clear(mess_krylov_arnoldi *arn){
    MSG_FNAME(__func__);
    int ret = 0;

    /*-----------------------------------------------------------------------------
     * clear fields
     *-----------------------------------------------------------------------------*/
    if(arn==NULL || (*arn)==NULL){
        return 0;
    }

    if((*arn)->V){
        ret = mess_matrix_clear(&((*arn)->V));              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_clear);
    }

    if((*arn)->H){
        ret = mess_matrix_clear(&((*arn)->H));              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_clear);
    }

    if((*arn)->_w){
        ret = mess_vector_clear(&((*arn)->_w));              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_clear);
    }

    if((*arn)->_wtmp){
        ret = mess_vector_clear(&((*arn)->_wtmp));          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_clear);
    }

    mess_free(*arn);
    return 0;
}



/**
 * @brief Create an Arnoldi process from @ref mess_mvpcall.
 * @param[in] arn               input, pointer to a @ref mess_krylov_arnoldi object
 * @param[in] A                 input, represents the operator \f$ A\f$ in Arnoldi process.
 * @param[in] sv                input, nonzero initial vector, must have same datatype as @p mess_mvpcall.
 * @param[in] k_est_max_steps   input, number of arnoldi steps, positive and smaller than the order of the operator @p A.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_create function prepares the Arnoldi process.
 * The argument @p k_est_max_steps is an estimation of the maximal number of iterations needed.
 * If more steps are needed, then memory reallocation occur for the Hessenberg matrix and matrix of Arnoldi vectors.
 *
 * @ref mess_krylov_arnoldi_init
 * @ref mess_krylov_arnoldi_clear
 *
 */
int mess_krylov_arnoldi_create(mess_krylov_arnoldi arn, mess_mvpcall A, mess_vector sv, mess_int_t k_est_max_steps){
    MSG_FNAME(__func__);
    int ret = 0;
    double nrm;

    /*-----------------------------------------------------------------------------
     *  check input args
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(arn);
    mess_check_nullpointer(A);
    mess_check_nullpointer(sv);
    mess_check_same_datatype(A,sv);
    mess_check_same_dim(A,sv);
    mess_check_same_dim(A,sv);
    mess_check_positive(k_est_max_steps);

    if(k_est_max_steps >= A->dim){
        MSG_ERROR("The maximal number of arnoldi steps should be smaller than the order of the operator: k="MESS_PRINTF_INT", order="MESS_PRINTF_INT"\n", k_est_max_steps, A->dim);
        return MESS_ERROR_ARGUMENTS;
    }

    ret = mess_vector_norm(sv, MESS_2_NORM, &nrm);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm);
    if (nrm < sqrt(mess_eps())){
        MSG_ERROR("2-norm of initial vector is too small\n", nrm);
        return MESS_ERROR_ARGUMENTS;
    }

    /*-----------------------------------------------------------------------------
     *  set and prepare arnoldi process
     *-----------------------------------------------------------------------------*/
    mess_int_t kp1 = k_est_max_steps + 1;
    arn->j = 0;
    arn->dim = A->dim;
    arn->data_type = A->data_type;
    arn->breakdown = 0;
    arn->breakdown_tol = ARNOLDI_BREAKDOWN_TOL;

    ret = mess_matrix_init(&(arn->V));                                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_alloc(arn->V, arn->dim, kp1, arn->dim*kp1, MESS_DENSE, arn->data_type);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    ret = mess_matrix_init(&(arn->H));                                                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_alloc(arn->H, kp1, kp1-1, kp1*(kp1-1), MESS_DENSE, arn->data_type);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    ret = mess_vector_init(&(arn->_w));                                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_copy(sv, arn->_w);                                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    ret = mess_vector_scale(1.0/nrm, arn->_w);                                                                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);
    ret = mess_matrix_setcol(arn->V,0,arn->_w);                                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

    ret = mess_vector_init(&(arn->_wtmp));                                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(arn->_wtmp, arn->dim, arn->data_type);                                                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);

    return ret;
}


/**
 * @brief Print fields of @ref mess_krylov_arnoldi structure.
 * @param[in] arn               input, pointer to a @ref mess_krylov_arnoldi object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_print function prints a @ref mess_krylov_arnoldi structure.
 *
 */
int mess_krylov_arnoldi_print(mess_krylov_arnoldi arn){

    /*-----------------------------------------------------------------------------
     *  print fields
     *-----------------------------------------------------------------------------*/
    if(!arn){
        MSG_PRINT("mess_krylov_arnoldi points to NULL\n");
    }else{
        MSG_PRINT("Fields:\n");
        MSG_PRINT("\tNumber of Arnoldi steps taken: "MESS_PRINTF_INT"\n", arn->j);
        MSG_PRINT("\tNumber of Arnoldi Vectors:     "MESS_PRINTF_INT"\n", arn->j+1);
        MSG_PRINT("\tDimension of Operator:         "MESS_PRINTF_INT"\n", arn->dim);
        MSG_PRINT("\tDatatype of Operator:          %s\n", mess_datatype_t_str(arn->data_type));
        MSG_PRINT("\tEarly Breakdown occured:       "MESS_PRINTF_INT"\n", arn->breakdown);
        MSG_PRINT("\tEarly Breakdown tolerance:     %e\n", arn->breakdown_tol);
        if(arn->V){
            MSG_PRINT("\tOrthonormal Basis V:           %s, "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT"\n", mess_datatype_t_str(arn->V->data_type), arn->V->rows, arn->V->cols);
        }else{
            MSG_PRINT("\tOrthonormal Basis V:           NULL\n");
        }
        if(arn->H){
            MSG_PRINT("\tHessenberg Matrix H:           %s, "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT"\n", mess_datatype_t_str(arn->H->data_type), arn->H->rows, arn->H->cols);
        }else{
            MSG_PRINT("\tHessenberg Matrix H:           NULL\n");
        }
    }
    return 0;
}

/**
 * @brief Perform exactly one steps of Arnoldi process.
 * @param[in] arn   input, @ref mess_krylov_arnoldi process.
 * @param[in] A     input, @ref mess_mvpcall A
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_step function performs exactly one step
 * of Arnoldi process.
 * If the number of computed arnoldi vectors is equal to the order of
 * @ref A an @ref MESS_ERROR_CONVERGE is returned.
 * If @p arn has already reached state of early breakdown @ref MESS_ERROR_CONVERGE is
 * also returned.
 *
 */
int mess_krylov_arnoldi_step(mess_krylov_arnoldi arn, mess_mvpcall A){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t reorth, i;
    double gamma, beta;
    mess_double_cpx_t gamma_cpx;

    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(arn);
    mess_check_nullpointer(arn->V);
    mess_check_nullpointer(arn->H);
    mess_check_nullpointer(arn->_w);
    mess_check_nullpointer(arn->_wtmp);
    mess_check_same_colsrows(arn->V, arn->H);
    mess_check_nullpointer(A);

    /*-----------------------------------------------------------------------------
     *  check if additional arnoldi step is possible
     *-----------------------------------------------------------------------------*/
    if(arn->breakdown){
        MSG_WARN("Arnoldi process reached state of early breakdown.\n");
        return MESS_ERROR_CONVERGE;
    }

    if(arn->j+1>=arn->dim){
        MSG_ERROR("The number of Arnoldi vectors is equal or larger the order of the operator, j="MESS_PRINTF_INT", dim="MESS_PRINTF_INT".\n", arn->j+1, arn->dim);
        return MESS_ERROR_CONVERGE;
    }

    /*-----------------------------------------------------------------------------
     *  check if new memory allocation is needed
     *-----------------------------------------------------------------------------*/
    if(arn->V->cols<arn->j+1){
        //number of cols of V and H are the same, we realloc more memory for both
        MSG_WARN("mess_krylov_arnoldi needs to increase the memory for the Arnoldi vectors and Hesenberg matrix.\n");

        MSG_WARN_NOBT("mess_krylov_arnoldi resize V from "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT" to "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT"\n",
                arn->V->rows, arn->V->cols,
                arn->V->rows, arn->V->cols+ARNOLDI_RESIZE_CHUNK);
        ret = mess_matrix_resize(arn->V, arn->V->rows, arn->V->cols+ARNOLDI_RESIZE_CHUNK);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);

        MSG_WARN("mess_krylov_arnoldi resize H from "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT" to "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT"\n",
                arn->H->rows, arn->H->cols,
                arn->H->rows+ARNOLDI_RESIZE_CHUNK, arn->H->cols+ARNOLDI_RESIZE_CHUNK);
        ret = mess_matrix_resize(arn->H, arn->H->rows + ARNOLDI_RESIZE_CHUNK, arn->H->cols + ARNOLDI_RESIZE_CHUNK);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_resize);
    }

    /*-----------------------------------------------------------------------------
     *  perform exactly one arnoldi step
     *-----------------------------------------------------------------------------*/
    //perform A_w -> _wtmp
    ret = mess_matrix_getcol(arn->V, arn->j, arn->_w);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
    ret = mess_mvpcall_apply(A, MESS_OP_NONE, arn->_w, arn->_wtmp);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), A.mvp);

    //perform double reorthogonalization, via modified gram schmidt
    for (reorth = 0; reorth<2; ++reorth){
        for (i = 0; i <= arn->j; i++){
            if (MESS_IS_REAL(arn)){
                ret = mess_matrix_colvecdot(arn->V, i, arn->_wtmp, &gamma);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colvecdot);
                arn->H->values[arn->j*arn->H->ld+i] += gamma;
                ret = mess_matrix_colvecaxpy(-gamma,i, arn->V, arn->_wtmp);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colvecaxpy);
            }else{
                ret = mess_matrix_colvecdotc(arn->V, i, arn->_wtmp, &gamma_cpx);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colvecdotc);
                arn->H->values_cpx[arn->j*arn->H->ld+i] += gamma_cpx;
                ret = mess_matrix_colvecaxpy(-gamma_cpx, i, arn->V, arn->_wtmp);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colvecaxpy);
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  add vector to basis V
     *-----------------------------------------------------------------------------*/
    // normalize new vector
    ret = mess_vector_norm(arn->_wtmp, MESS_2_NORM, &beta);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_norm);

    // update Hessenberg
    ret = mess_matrix_setelement(arn->H, arn->j+1, arn->j, beta);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setelement);

    // check for breakdown
    if(beta<=arn->breakdown_tol){
        arn->breakdown = 1;
        MSG_WARN("%s: early breakdown beta = %e at j = " MESS_PRINTF_INT "\n", __func__, beta, arn->j);

        // resize matrices
        ret = mess_matrix_resize(arn->V, arn->V->rows, arn->j+1);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_resize);
        ret = mess_matrix_resize(arn->H, arn->j+1, arn->j+1);                       FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_resize);

        return 0;
    }

    // add to basis
    ret = mess_vector_scale( 1.0/beta, arn->_wtmp);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);
    ret = mess_matrix_setcol(arn->V, arn->j+1, arn->_wtmp);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);

    // increment counter
    arn->j+=1;

    /*-----------------------------------------------------------------------------
     *  finished one step of arnoldi
     *-----------------------------------------------------------------------------*/
    return 0;
}

/**
 * @brief Return number of computed arnoldi vectors of @ref mess_krylov_arnoldi.
 * @param[in] arn                   input, @ref mess_krylov_arnoldi process.
 * @param[out] no_arnoldi_vectors   output, on exit number of Arnoldi vectors of @p arn
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_number_arnoldi_vectors function returns the number
 * of Arnoldi vectors which @p arn holds.
 *
 */
int mess_krylov_arnoldi_number_arnoldi_vectors(mess_krylov_arnoldi arn, mess_int_t* no_arnoldi_vectors){
    MSG_FNAME(__func__);
    mess_check_nullpointer(arn);
    mess_check_nullpointer(no_arnoldi_vectors);
    *no_arnoldi_vectors=arn->j+1;
    return 0;
}


/**
 * @brief Return matrix which contains the Arnoldi vectors columnwise.
 * @param[in] arn   input, @ref mess_krylov_arnoldi process.
 * @param[out] V    output, matrix of Arnoldi vectors columnwise
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_getV function returns the Arnoldi vectors which @p arn
 * currently holds in a matrix @p V columnwise.
 *
 */
int mess_krylov_arnoldi_getV(mess_krylov_arnoldi arn, mess_matrix V){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(arn);
    mess_check_nullpointer(V);
    mess_check_positive(arn->j);
    ret = mess_matrix_colsub(arn->V, 0, arn->j, V);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_colsub);
    return 0;
}


/**
 * @brief Return Hessenberg matrix generated by Arnoldi process.
 * @param[in] arn   input, @ref mess_krylov_arnoldi process.
 * @param[out] H    output, generated Hessenberg matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_getH function returns the Hessenberg
 * Matrix generated by @p arn,
 */
int mess_krylov_arnoldi_getH(mess_krylov_arnoldi arn, mess_matrix H){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_check_nullpointer(arn);
    mess_check_nullpointer(H);
    mess_check_positive(arn->j);

    if(arn->breakdown){
        ret = mess_matrix_sub(arn->H, 0, arn->j, 0, arn->j, H);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sub);
    }else{
        ret = mess_matrix_sub(arn->H, 0, arn->j, 0, arn->j-1, H);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_sub);
    }
    return 0;
}

/**
 * @brief Return largest eigenvalue of Hessenberg matrix and error estimation.
 * @param[in] arn       input, @ref mess_krylov_arnoldi process.
 * @param[out] lambda   output, largest eigenvalue of Hessenberg matrix in magnitude
 * @param[out] err_est  output, error estimation for eigen residual.
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_krylov_arnoldi_H_largest_eigenvalue returns the largest eigenvalue in magnitude
 * of the Hessenberg matrix generated by the Arnoldi process @p arn.
 * @p err_est estimates the error of the eigenresidual
 * \f$ | Av - \lambda v  |_2 \f$. For more information read
 * Proposition 6.8 of @cite Saa92.
 * Matrix generated by @p arn,
 */
int mess_krylov_arnoldi_H_largest_eigenvalue(mess_krylov_arnoldi arn, mess_double_cpx_t* lambda, double* err_est){
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t no_arnoldi_vectors;
    mess_int_t Hrows, Hcols;
    mess_vector evector;

    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(arn);
    mess_check_nullpointer(lambda);
    mess_check_nullpointer(err_est);

    /*-----------------------------------------------------------------------------
     *  init eigenvector
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&evector);

    /*-----------------------------------------------------------------------------
     *  get number of arnoldi vectors, resize matrix H internally and compute
     *  largest eigenvalue
     *-----------------------------------------------------------------------------*/
    ret = mess_krylov_arnoldi_number_arnoldi_vectors(arn, &no_arnoldi_vectors);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_number_arnoldi_vectors);

    mess_check_positive(no_arnoldi_vectors);
    if(no_arnoldi_vectors==1){
        //this is in general only possible if initial vector is an eigenvector
        *lambda = MESS_IS_REAL(arn->H)? arn->H->values[0]: arn->H->values_cpx[0];
        *err_est = 0;
        return 0;
    }

    Hrows = arn->H->rows;
    Hcols = arn->H->cols;

    //number of arnoldi vectors available
    arn->H->rows = no_arnoldi_vectors-1;
    arn->H->cols = no_arnoldi_vectors-1;

    //compute largest eigenvalue in magnitude of Hessenberg matrix
    ret = mess_eigen_hessenberg_abs_largest(arn->H, evector, lambda);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_hessenberg_abs_largest);

    //resize back
    arn->H->rows = Hrows;
    arn->H->cols = Hcols;

    /*-----------------------------------------------------------------------------
     *  compute error via proposition 6.8, if no breakdown occured
     *-----------------------------------------------------------------------------*/
    if(arn->breakdown){
        *err_est = 0;
    }else{
        //mess_matrix_print(arn->H);
        mess_double_cpx_t hnp1n, yn;
        ret = mess_matrix_getelement(arn->H, no_arnoldi_vectors-1, no_arnoldi_vectors-2, NULL, &hnp1n);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getelement);
        ret = mess_vector_get(evector, evector->dim-1, &yn);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_get);
        //printf("hnp1n               %e + %e I \n", creal(hnp1n), cimag(hnp1n));
        //printf("yn                  %e + %e I \n", creal(yn), cimag(yn));
        //printf("arnoldi vectors     %d \n", no_arnoldi_vectors);
        //mess_vector_print(evector);
        *err_est = cabs(hnp1n)*cabs(yn);
    }

    /*-----------------------------------------------------------------------------
     *  clear and return
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_VECTORS(&evector);

    return 0;
}





/**
 * @brief Compute matrices \f$ H \f$ and \f$ V \f$ with Arnoldis method with respect to \f$ A \f$.
 * @param[in] A  input operator applying the \f$ A \f$ matrix for which the Arnoldi algorithm is supposed to run
 * @param[in] k  input number of Arnoldi steps (usually \f$ k \ll n \f$)
 * @param[in] sv   input start vector
 * @param[out] H  \f$(k+1 \times k) \f$ upper Hessenberg matrix
 * @param[out] V  \f$ (n \times k+1)  \f$ orthogonal matrix
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_arnoldi_template function computes matrices \f$ H  \f$ and  \f$ V  \f$ (if  \f$ V  \ne \f$ @c NULL)
 * fullfilling
 * \f[
 * \begin{array}{ccc}
 *      V(:,1) &=& \frac{sv}{ \Vert sv \Vert }, \\
 *      V^T V  &=& I, \\
 *   AV(:,1:k) &=& VH .
 * \end{array}
 * \f]
 * If an early breakdown occurs then it holds that  \f$ AV =V H\f$.
 * The matrix \f$ A \f$ is given as a function handler which computes \f$ y=Ax \f$.
 *
 */
int mess_eigen_arnoldi_template(mess_mvpcall A, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V){
    MSG_FNAME(__func__);
    int ret;
    mess_int_t j;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sv);
    mess_check_nullpointer(H);
    mess_check_real_or_complex(A);
    mess_check_positive(k);

    mess_check_same_dim(A,sv);
    mess_check_same_datatype(A,sv);


    /*-----------------------------------------------------------------------------
     *  prepare mess_krylov_arnoldi struct
     *-----------------------------------------------------------------------------*/
    mess_krylov_arnoldi arn;

    ret = mess_krylov_arnoldi_init(&arn);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_init);
    ret = mess_krylov_arnoldi_create(arn, A, sv, k);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_create);

    /*-----------------------------------------------------------------------------
     *  run arnoldi process
     *-----------------------------------------------------------------------------*/
    for ( j = 0; j < k; j++){
        // compute next arnoldi vector and check for breakdown
        ret = mess_krylov_arnoldi_step(arn, A);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_step);
        //mess_krylov_arnoldi_print(arn);

        if(arn->breakdown){
            break;
        }
    }

    /*-----------------------------------------------------------------------------
     *  finalize and clear
     *-----------------------------------------------------------------------------*/
    if ( V != 0) {
        ret = mess_krylov_arnoldi_getV(arn, V);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_getV);
    }
    ret = mess_krylov_arnoldi_getH(arn, H);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_getH);

    ret = mess_krylov_arnoldi_clear(&arn);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_clear);

    return(0) ;

}


/**
 * @brief Compute absolut value of the largest eigenvalue using Arnoldi process.
 * @param[in] A         input operator applying the \f$ A \f$ matrix for which the Arnoldi algorithm is supposed to run
 * @param[in] k         input number of Arnoldi steps (usually \f$ k \ll n \f$)
 * @param[in] sv        input start vector, non zero
 * @param[out] abseig   output largest eigenvalue in absolute value
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_eigen_arnoldi_template_nrm computes an approximation of
 * the largest eigenvalue in magnitude and returns the magnitude.
 * The function use an arnoldi process and computes the ritz values.
 *
 * The algorithms stops if the maximum iteration number \f$ k \f$ is reached or
 * if the error estimation is smaller than \f$ n eps() \f$ where \f$ n \f$ is the
 * order of the operator.
 * The matrix \f$ A \f$ is given as a function handler which computes \f$ y = Ax \f$.
 *
 */
int mess_eigen_arnoldi_template_nrm(mess_mvpcall A, mess_int_t k, mess_vector sv, double *abseig){
    MSG_FNAME(__func__);
    int ret;
    mess_int_t j;
    double tol, err_est;
    mess_double_cpx_t largest_ev;

    /*-----------------------------------------------------------------------------
     *  check input parameters
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(sv);
    mess_check_real_or_complex(A);
    mess_check_positive(k);

    mess_check_same_dim(A,sv);
    mess_check_nullpointer(abseig);
    *abseig = 1;

    mess_check_same_datatype(A,sv);

    /*-----------------------------------------------------------------------------
     *  prepare mess_krylov_arnoldi struct, initial vector is added as first
     *  arnoldi vector
     *-----------------------------------------------------------------------------*/
    mess_krylov_arnoldi arn;
    ret = mess_krylov_arnoldi_init(&arn);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_init);
    ret = mess_krylov_arnoldi_create(arn, A, sv, k);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_create);

    /*-----------------------------------------------------------------------------
     *  prepare arnoldi process
     *-----------------------------------------------------------------------------*/
    tol = (A->dim)*mess_eps();

    /*-----------------------------------------------------------------------------
     *  run arnoldi process
     *-----------------------------------------------------------------------------*/
    for ( j = 0; j < k; j++){
        // compute next arnoldi vector and check for breakdown
        ret = mess_krylov_arnoldi_step(arn, A);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_step);
        //mess_krylov_arnoldi_print(arn);

        /*-----------------------------------------------------------------------------
         *  the eigenvalue related part
         *-----------------------------------------------------------------------------*/
        ret = mess_krylov_arnoldi_H_largest_eigenvalue(arn, &largest_ev, &err_est);
        *abseig = cabs(largest_ev);

        /*-----------------------------------------------------------------------------
         *  check for breakdown
         *-----------------------------------------------------------------------------*/
        if(arn->breakdown){
            break;
        }

        /*-----------------------------------------------------------------------------
         *  check for error estimate
         *-----------------------------------------------------------------------------*/
        if(err_est < tol ){
            MSG_INFO("Iteration "MESS_PRINTF_INT",  estimated error = %e, tolerance = %e\n", j, err_est, tol);
            break;
        }

    }

    /*-----------------------------------------------------------------------------
     *  clear
     *-----------------------------------------------------------------------------*/
    ret = mess_krylov_arnoldi_clear(&arn);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_krylov_arnoldi_clear);

    return 0;
}
