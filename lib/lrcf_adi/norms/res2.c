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
 * @file lib/lrcf_adi/norms/res2.c
 * @brief   \f$ 2 \f$ Norm residual computation for the Lyapunov and the Riccati operator.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <assert.h>
#include <complex.h>


#define VECTOR_SCALE(factor,vector) if (MESS_IS_REAL((vector))) { \
    mess_vector_scale((factor), (vector));\
} else { \
    mess_vector_scalec((factor), (vector));\
}
#define VECTOR_AXPY(factor, v1, v2) if (MESS_IS_REAL(v1)) {\
    mess_vector_axpy((factor),(v1),(v2)); \
} else { \
    mess_vector_axpyc((factor),(v1),(v2)); \
}
#define MATRIX_ADD(alpha,A,beta,B) if ( MESS_IS_COMPLEX((B))){ \
    mess_matrix_addc((alpha),(A),(beta),(B));\
} else { \
    mess_matrix_add((alpha),(A),(beta), (B));\
}

#define IT_MAX 20


typedef struct res2_data_st {
    mess_matrix F;
    mess_matrix G;
    mess_matrix Z;
    mess_matrix K;
    mess_matrix B;
    mess_matrix M;
    mess_vector t,t1,t2, t3,t4;
} res2_data;

static int res2_mvp(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->G, x, d->t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->G, d->t, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = AZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,x,d->t1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    // mess_vector_printshort(d->t2);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->F, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    // y += ZZ'A'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->F, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->Z, d->t1, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    // mess_vector_printshort(y);
    // scanf("%d",&ret);
    return 0;
}
static int res2_mvpt(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->G, x, d->t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->G, d->t, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = AZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,x,d->t1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_gaxpy(MESS_OP_HERMITIAN,d->F, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y += ZZ'A'x
    ret = mess_matrix_mvp(MESS_OP_NONE,d->F, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->Z, d->t1, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

/**
 * @brief Computes a matrix-vector product of a Lyapunov Equation and a given vector.
 * @param[in] data input data
 * @param[in] op  input operation
 * @param[in] x input vector
 * @param[out] y vector
 * @return zero on succes or a non zero error value
 *
 * The @ref res2_mvpBK function computes a vector \f$ y \f$ defined by a matrix-vector product of a Lyapunov Equation and
 * a given input vector \f$ x \f$:
 * \f[ y = \left( (F-BK)ZZ^T + ZZ^T(F-BK)^T + GG^T \right). \f]
 */
int res2_mvpBK(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;
    // y = GG'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->G, x, d->t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->G, d->t, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);


    // y += (F-BK)ZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, x, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->F,d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->K,d->t2,d->t3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->B,d->t3,d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_axpy(-1, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y += ZZ'(F'-K'B')x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->F, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->B, x, d->t3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->K, d->t3, d->t4); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_vector_axpy(-1, d->t4, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t4); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t4); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(1, d->t4, y ) ; FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // mess_vector_printshort(y);
    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

static int res2g_mvp(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->G, x, d->t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->G, d->t, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = AZZ'M'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->M,x,d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,d->t2,d->t1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    // mess_vector_printshort(d->t2);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->F, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    // y += MZZ'A'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->F, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->M, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

static int res2gt_mvp(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->G, x, d->t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->G, d->t, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = A'ZZ'Mx
    ret = mess_matrix_mvp(MESS_OP_NONE,d->M,x,d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,d->t2,d->t1);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    // mess_vector_printshort(d->t2);
    ret = mess_matrix_gaxpy(MESS_OP_HERMITIAN,d->F, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y += M'ZZ'Ax
    ret = mess_matrix_mvp(MESS_OP_NONE,d->F, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_gaxpy(MESS_OP_HERMITIAN,d->M, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    // ret = mess_vector_printshort(y);
    return 0;
}

/**
 * @brief Computes a matrix-vector product of a generalized Lyapunov Equation and a given vector.
 * @param[in] data input data
 * @param[in] op  input operation
 * @param[in] x input vector
 * @param[out] y vector
 * @return zero on succes or a non zero error value
 *
 * The @ref res2_mvpBK function computes a vector \f$ y \f$ defined by a matrix-vector product of a generalized
 * Lyapunov Equation and a given input vector \f$ x \f$:
 * \f[ y = \left( (F-BK)ZZ^T M + M ZZ^T(F-BK)^T + GG^T \right). \f]
 */

int res2g_mvpBK(void *data,mess_operation_t op,  mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;
    // y = GG'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->G, x, d->t); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->G, d->t, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);


    // y += (F-BK)ZZ'M'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->M, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_gaxpy(MESS_OP_NONE,d->F,d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->K,d->t2,d->t3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->B,d->t3,d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_axpy(-1, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y += MZZ'(F'-K'B')x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->F, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->B, x, d->t3); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->K, d->t3, d->t4); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_vector_axpy(-1, d->t4, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t4); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t4); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_gaxpy(MESS_OP_NONE,d->M, d->t4, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);
    // ret = mess_vector_axpy(1, d->t4, y ) ; FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // mess_vector_printshort(y);
    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}



/**
 * @brief Compute \f$ 2\f$-norm of the Lyapunov residual via an Arnoldi process.
 * @param[in] F             input system matrix
 * @param[in] G             input right hand side
 * @param[in] Z             input computed \f$ Z \f$
 * @param[out] nrm          output computed \f$ 2 \f$-norm of the residual
 * @return zero on succes or a non zero error value
 *
 * The @ref mess_lrcfadi_res2 function computes the \f$ 2 \f$-norm of the Lyapunov residual of \f$ Z \f$ exploiting the
 * symmetry of the residual operator
 * \f[ R := FZZ^T + ZZ^TF^T + GG^T. \f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploiting the sparsity of \f$ F \f$ and the rectangular
 * structure of \f$ Z \f$ and \f$ G \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm.\n
 * The default maximum number of iterations is \f$ 20 \f$. \n
 * The default tolerance is set to \f$ \sqrt{u}\f$.
 *
 */
int mess_lrcfadi_res2 ( mess_matrix F, mess_matrix G, mess_matrix Z, double *nrm)
{
    MSG_FNAME(__func__);
    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv =NULL;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(F);
    mess_check_nullpointer(G);
    mess_check_nullpointer(Z);
    mess_check_square(F);
    mess_check_real(F);
    mess_check_real(G);

    if (Z->rows != F->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (G->rows != F->rows) {
        MSG_ERROR("G has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    MESS_INIT_VECTORS(&(data.t1),&(data.t2),&(data.t));
    ret = mess_vector_alloc((data.t1), Z->cols, Z->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), F->rows, F->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t), G->cols, G->data_type);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    data.Z=Z;
    data.F=F;
    data.G=G;

    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/
    ret = mess_mvpcall_operator(&mvpcall, F->rows, MESS_REAL, res2_mvp, &data);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
    ret = mess_vector_init(&sv);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);


    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);         FUNCTION_FAILURE_HANDLE(ret,( ret!=0), mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);
    return 0;
}

/**
 * @brief Compute \f$ 2 \f$ norm of the Lyapunov residual (transposed version).
 * @param[in] F             input system matrix
 * @param[in] G             input right hand side
 * @param[in] Z             input computed \f$ Z \f$S
 * @param[out] nrm          output computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2t function computes the \f$ 2 \f$-norm of the Lyapunov residual of \f$ Z \f$ exploiting the
 * symmetry of the residual operator
 * \f[ R := F'ZZ^T + ZZ^TF + GG^T \f]
 * or
 * \f[ R := F'ZZ^T + ZZ^TF + G^TG \f]
 * depending on the dimensions of \f$ G \f$. \n
 * It is automatically transposed internal if \f$ G \f$ has the wrong dimension. \n
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploiting the sparsity of \f$ F \f$ and rectangular
 * structure of \f$ Z \f$ and \f$ G \f$. \n
 * Thus it can be computed in O(n) effort andis therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm. \n
 * The default maximum number of iterations is \f$ 20 \f$. \n
 * The default tolerance is set to \f$ \sqrt{u}\f$.
 *
 */
int mess_lrcfadi_res2t(mess_matrix F, mess_matrix G, mess_matrix Z, double *nrm){

    MSG_FNAME(__func__);
    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv =NULL;
    mess_matrix Gt = NULL;
    int conGt = 0 ;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(F);
    mess_check_nullpointer(G);
    mess_check_nullpointer(Z);
    mess_check_square(F);
    mess_check_real(F);
    mess_check_real(G);

    if (Z->rows != F->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    /*-----------------------------------------------------------------------------
     *  transpose the rhs if it is need.
     *-----------------------------------------------------------------------------*/
    if ( G->rows < G->cols) {
        ret = mess_matrix_init(&Gt);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_ctranspose(G,Gt); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
        conGt = 1;
    }else {
        Gt = G;
    }

    if (Gt->rows != F->rows) {
        MSG_ERROR("G has the wrong number of rows/cols. \n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.t1), &(data.t2), &(data.t))
    ret = mess_vector_alloc((data.t1), Z->cols, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), F->rows, F->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t), Gt->cols, G->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    data.Z=Z;
    data.F=F;
    data.G=Gt;

    ret = mess_mvpcall_operator(&mvpcall, F->rows, MESS_REAL, res2_mvpt, &data);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif



    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    if ( conGt )
        mess_matrix_clear(&Gt);

    return 0;
}

/**
 * @brief Compute the \f$ 2 \f$-norm of the Lyapunov residual (\f$ BK \f$ version).
 * @param[in] F         input system matrix
 * @param[in] G         input right hand side
 * @param[in] B         input martix for \f$ (F - BK) \f$
 * @param[in] K         input matrix for \f$ (F - BK) \f$
 * @param[in] Z         input computed  \f$ Z \f$
 * @param[out] nrm      output computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2BK function computes the \f$ 2 \f$-norm of the Lyapunov residual of \f$ Z \f$ exploiting the
 * symmetry of the residual operator
 * \f[ R := (F-BK)ZZ^T + ZZ^T(F-BK)^T + GG^T .\f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploiting the sparsity of \f$ F \f$ and rectangular
 * structure of \f$ Z \f$ and \f$ G \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm. \n
 * The default maximum number of iterations is \f$ 20 \f$. \n
 * The default tolerance is set to \f$ \sqrt{u}\f$.
 *
 */
int mess_lrcfadi_res2BK (mess_matrix F, mess_matrix G, mess_matrix Z, mess_matrix B, mess_matrix K, double *nrm){
    MSG_FNAME(__func__);
    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv =NULL;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(F);
    mess_check_nullpointer(G);
    mess_check_nullpointer(Z);
    mess_check_square(F);
    mess_check_real(F);
    mess_check_real(G);
    mess_check_real(K);
    mess_check_real(B);

    if (Z->rows != F->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (G->rows != F->rows) {
        MSG_ERROR("G has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.t1),&(data.t2),&(data.t3),&(data.t4),&(data.t));
    ret = mess_vector_alloc((data.t1), Z->cols, Z->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), F->rows, F->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t3), B->cols, K->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t4), F->rows, F->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t), G->cols, G->data_type);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    data.Z=Z;
    data.F=F;
    data.G=G;
    data.K=K;
    data.B=B;

    ret = mess_mvpcall_operator(&mvpcall, F->rows, MESS_REAL, res2_mvpBK, &data);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif


    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t3));
    mess_vector_clear(&(data.t4));
    mess_vector_clear(&(data.t));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    return 0;

}

/**
 * @brief Compute the \f$ 2 \f$-norm of the generalized Lyapunov residual.
 * @param[in] F         input system matrix
 * @param[in] G         input right hand side
 * @param[in] M         input mass matrix
 * @param[in] Z         input computed \f$ Z \f$
 * @param[out] nrm      output computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2g function computes the \f$ 2 \f$-norm of the generalized Lyapunov residual of \f$ Z \f$
 * exploiting the symmetry of the residual operator
 * \f[ R := FZZ^TM^T+ MZZ^TF^T + GG^T. \f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploiting the sparsity of \f$ F \f$ and rectangular
 * structure of \f$ Z \f$ and \f$ G \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm. \n
 * The default maximum number of iterations is \f$ 20 \f$. \n
 * The default
 * tolerance is set to \f$ \sqrt{u}\f$.
 */
int mess_lrcfadi_res2g ( mess_matrix F,mess_matrix M, mess_matrix G, mess_matrix Z, double *nrm){
    MSG_FNAME(__func__);
    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv =NULL;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(F);
    mess_check_nullpointer(G);
    mess_check_nullpointer(M);
    mess_check_nullpointer(Z);
    mess_check_square(F);
    mess_check_real(F);
    mess_check_real(G);
    mess_check_real(M);
    mess_check_same_size(F,M);

    if (Z->rows != F->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (G->rows != F->rows) {
        MSG_ERROR("G has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/

    MESS_INIT_VECTORS(&(data.t1),&(data.t2),&(data.t));
    ret = mess_vector_alloc((data.t1), Z->cols, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), F->rows, F->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t), G->cols, G->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    data.Z=Z;
    data.F=F;
    data.G=G;
    data.M=M;
    ret = mess_mvpcall_operator(&mvpcall, F->rows, MESS_REAL, res2g_mvp, &data);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    MESS_INIT_VECTORS(&sv);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif


    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    return 0;

}

/**
 * @brief Compute the \f$ 2 \f$-norm of the generalized Lyapunov residual (transposed version).
 * @param[in] F         input system matrix
 * @param[in] G         input right hand side
 * @param[in] M         input mass matrix
 * @param[in] Z         input computed \f$ Z \f$
 * @param[out] nrm      output computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2gt function computes the 2 Norm of the generalized Lyapunov residual of Z exploiting the
 * symmetry of the residual operator
 * \f[ R := F^TZZ^TM+ M'ZZ^TF + GG^T. \f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploiting the sparsity of \f$ F \f$ and rectangular
 * structure of \f$ Z \f$ and \f$ G \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm. \n
 * The default maximum number of iterations is \f$ 20 \f$. \n
 * The default tolerance is set to \f$ \sqrt{u}\f$.
 * \n
 */
int mess_lrcfadi_res2gt( mess_matrix F,mess_matrix M, mess_matrix G, mess_matrix Z, double *nrm){
    MSG_FNAME(__func__);
    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv =NULL;
    mess_matrix Gt = NULL;
    int conGt = 0 ;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(F);
    mess_check_nullpointer(G);
    mess_check_nullpointer(M);
    mess_check_nullpointer(Z);
    mess_check_square(F);
    mess_check_real(F);
    mess_check_real(G);
    mess_check_real(M);
    mess_check_same_size(F,M);

    if (Z->rows != F->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    /*-----------------------------------------------------------------------------
     *  transpose the rhs if it is need.
     *-----------------------------------------------------------------------------*/

    if ( G->rows < G->cols) {
        MSG_INFO("G transpose\n");
        ret = mess_matrix_init(&Gt);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_ctranspose(G,Gt); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
        conGt = 1;
    }else {
        Gt = G;
    }

    if (Gt->rows != F->rows) {
        MSG_ERROR("G has the wrong number of rows/cols. \n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.t1),&(data.t2),&(data.t));
    ret = mess_vector_alloc((data.t1), Z->cols, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), F->rows, F->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t), Gt->cols, G->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    data.Z=Z;
    data.F=F;
    data.G=Gt;
    data.M=M;

    ret = mess_mvpcall_operator(&mvpcall, F->rows, MESS_REAL, res2gt_mvp, &data);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    if ( conGt )
        mess_matrix_clear(&Gt);

    return 0;

}

/**
 * @brief Compute the \f$ 2 \f$-norm of the generalized Lyapunov residual (\f$ BK \f$-Version).
 * @param[in] F         input system matrix
 * @param[in] G         input right hand side
 * @param[in] M         input mass matrix
 * @param[in] B         input matrix for \f$(F-BK) \f$
 * @param[in] K         input matrix for \f$(F-BK) \f$
 * @param[in] Z         input computed \f$ Z \f$
 * @param[out] nrm      output computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2gBK function computes the \f$ 2 \f$-norm of the generalized Lyapunov residual of \f$ Z \f$
 * exploiting the symmetry of the residual operator
 * \f[ R := FZZ^TM + ZZ^TF^TM + GG^T. \f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploiting the sparsity of \f$ F \f$ and rectangular
 * structure of \f$ Z \f$ and \f$ G \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm. \n
 * The default maximum number of iterations is \f$ 20 \f$.\n
 * The default
 * tolerance is set to \f$ \sqrt{u}\f$.
 */
int mess_lrcfadi_res2gBK ( mess_matrix F,mess_matrix M, mess_matrix G, mess_matrix Z, mess_matrix B, mess_matrix K, double *nrm){
    MSG_FNAME(__func__);
    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv =NULL;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(F);
    mess_check_nullpointer(G);
    mess_check_nullpointer(M);
    mess_check_nullpointer(B);
    mess_check_nullpointer(K);
    mess_check_nullpointer(Z);
    mess_check_square(F);
    mess_check_real(F);
    mess_check_real(G);
    mess_check_real(M);
    mess_check_real(B);
    mess_check_real(K);
    mess_check_same_size(F,M);

    if (Z->rows != F->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (G->rows != F->rows) {
        MSG_ERROR("G has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.t1),&(data.t2),&(data.t3),&(data.t),&(data.t4));
    ret = mess_vector_alloc((data.t1), Z->cols, Z->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), F->rows, F->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t3), K->rows, F->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t), G->cols, G->data_type);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t4), K->cols, F->data_type);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    data.Z=Z;
    data.F=F;
    data.G=G;
    data.M=M;
    data.K=K;
    data.B=B;
    ret = mess_mvpcall_operator(&mvpcall, F->rows, MESS_REAL, res2g_mvpBK, &data);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t3));
    mess_vector_clear(&(data.t4));
    mess_vector_clear(&(data.t));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    return 0;
}



typedef struct res2nm_data_st {
    mess_matrix A;
    mess_matrix B;
    mess_matrix C;
    mess_matrix Z;
    mess_matrix M;
    mess_vector tc,tz, t1,t2, tb, t3,t4;
} res2nm_data;


static int res2nm_mvp (void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2nm_data *d = ( res2nm_data *) data;
    int ret =0 ;

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y = C'Cx
    ret = mess_matrix_mvp(MESS_OP_NONE,d->C, x, d->tc); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->C, d->tc, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);

    // t1 = ZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, x, d->tz); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y+= A't1
    ret = mess_matrix_gaxpy(MESS_OP_HERMITIAN, d->A, d->t1, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);


    // y+= ZZ'Ax
    ret = mess_matrix_mvp(MESS_OP_NONE,d->A, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(1, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y-= ZZ'BB't1
    // MSG_PRINT("t1 = \n");
    // mess_vector_printshort(d->t1);
    // mess_matrix_printinfo(d->B);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->B, d->t1, d->tb); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->B, d->tb, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(-1, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

static int res2nmg_mvp (void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2nm_data *d = ( res2nm_data *) data;
    int ret =0 ;

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y = C'Cx
    ret = mess_matrix_mvp(MESS_OP_NONE,d->C, x, d->tc); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->C, d->tc, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);

    // t1 = ZZ'Mx
    ret = mess_matrix_mvp(MESS_OP_NONE,d->M, x, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t1, d->tz); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y+= A't1
    ret = mess_matrix_gaxpy(MESS_OP_HERMITIAN, d->A, d->t1, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);


    // y+= M'ZZ'Ax
    ret = mess_matrix_mvp(MESS_OP_NONE,d->A, x, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_gaxpy(MESS_OP_HERMITIAN, d->M,d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);
    // ret = mess_vector_axpy(1, d->t2, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y-= M'ZZ'BB't1
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->B, d->t1, d->tb); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->B, d->tb, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->M, d->t2, d->t1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_vector_axpy(-1, d->t1, y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}
/**
 * @brief Compute the \f$ 2 \f$-norm of the Riccati residual.
 * @param[in] A         input matrix
 * @param[in] B         input matrix
 * @param[in] C         input matrix
 * @param[in] Z         input computed \f$ Z \f$
 * @param[out] nrm      output computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2nm function computes the \f$ 2 \f$-norm of the Riccati residual of \f$ Z \f$ exploiting the
 * symmetry of the residual operator
 * \f[ R := A^T  Z Z^T + Z Z^T A + C^TC - Z Z^T B B^T Z Z^T .\f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploitig the sparsity of \f$ A  \f$ and rectangular
 * structure of \f$ B,C \f$ and \f$ Z \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm.
 *
 */

int mess_lrcfadi_res2nm(mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix Z, double *nrm){
    MSG_FNAME(__func__);
    res2nm_data data;
    mess_mvpcall mvpcall;
    mess_vector sv;
    int ret =0 ;


    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(nrm);
    mess_check_square(A);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_real(C);


    if (Z->rows != A->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (B->rows != A->rows) {
        MSG_ERROR("A has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (C->cols != A->cols) {
        MSG_ERROR("C has the wrong number of columns. \n");
        return MESS_ERROR_DIMENSION;
    }


    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/

    MESS_INIT_VECTORS(&(data.tc),&(data.tz),&(data.t1),&(data.t2),&(data.tb));
    ret = mess_vector_alloc((data.tc), C->rows, C->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.tz), Z->cols, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t1), A->rows, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), A->rows, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.tb), B->cols, B->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

    data.A=A;
    data.B=B;
    data.C=C;
    data.Z=Z;

    ret = mess_mvpcall_operator(&mvpcall, A->rows, MESS_REAL, res2nm_mvp, &data);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.tc));
    mess_vector_clear(&(data.tz));
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.tb));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);


    return 0;
}


/**
 * @brief Compute the \f$ 2 \f$-norm of the generalized Riccati residual.
 * @param[in] A   input matrix
 * @param[in] B  input matrix
 * @param[in] C  input matrix
 * @param[in] M  input matrix
 * @param[in] Z  input computed \f$ Z \f$
 * @param[out] nrm computed \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_res2nmg function computes the \f$ 2 \f$-norm of the generalized Riccati residual of \f$ Z \f$
 * exploiting the symmetry of the residual operator
 * \f[ R := A^T  Z*Z^T M + M^T Z Z^T  A + C^T C - M^T Z Z^T B B^T Z Z^T M .\f]
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploitig the sparsity of \f$ A  \f$ and rectangular
 * structure of \f$ B,C \f$ and \f$ Z \f$. \n
 * Thus it can be computed in O(n) effort and is therefore much cheaper than the computation of the e.g. the
 * Frobenius-norm.
 *
 */

int mess_lrcfadi_res2nmg(mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix M, mess_matrix Z, double * nrm){
    MSG_FNAME(__func__);
    res2nm_data data;
    mess_mvpcall mvpcall;
    mess_vector sv;
    int ret =0 ;

    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_nullpointer(M);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(nrm);
    mess_check_square(A);
    mess_check_square(M);
    mess_check_real(A);
    mess_check_real(B);
    mess_check_real(C);
    mess_check_real(M);
    mess_check_same_size(A,M);

    if (Z->rows != A->rows) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (B->rows != A->rows) {
        MSG_ERROR("A has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }
    if (C->cols != A->cols) {
        MSG_ERROR("C has the wrong number of columns. \n");
        return MESS_ERROR_DIMENSION;
    }


    /*-----------------------------------------------------------------------------
     *  build mvp call operator and init random inital vector
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.tc),&(data.tz),&(data.t1),&(data.t2),&(data.tb));
    ret = mess_vector_alloc((data.tc), C->rows, C->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.tz), Z->cols, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t1), A->rows, A->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), A->rows, A->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.tb), B->cols, B->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

    data.A=A;
    data.B=B;
    data.C=C;
    data.Z=Z;
    data.M=M;

    ret = mess_mvpcall_operator(&mvpcall, A->rows, MESS_REAL, res2nmg_mvp, &data);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
    ret = mess_vector_init(&sv);                                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

    /*-----------------------------------------------------------------------------
     *  use arpack if available otherwise arnoldi
     *-----------------------------------------------------------------------------*/
#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = IT_MAX;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev,2,MESS_REAL);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        if(MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        }else{
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, nrm);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&(data.tc));
    mess_vector_clear(&(data.tz));
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.tb));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    return 0;
}


