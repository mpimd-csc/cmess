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
 * @file lib/lrcf_adi/norms/res2_template.c
 * @brief New LRCF-ADI res2 implementation for the eqn structure.
 *
 * @author @koehlerm
 * @author @mbehr
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#include "mess/mess.h"
#include "mess/error_macro.h"


#define ARNOLDI_ITERATIONS 100
#define ARPACK_NCV 10


static int mess_lrcfadi_residual_lyap(mess_equation eqn, mess_options opt, mess_matrix Z, double *nrm);
static int mess_lrcfadi_residual_nm(mess_equation eqn, mess_options opt, mess_matrix Z, double *nrm);



/**
 * @brief Compute \f$ 2 \f$-norm residual of an equation.
 * @param[in]       eqn     input   object defining an equation
 * @param[in]       opt     input   options for the equation and the solver
 * @param[in]       Z       input   factor \f$ Z \f$ of the solution \f$ X \approx Z Z^T \f$
 * @param[out]      nrm     output          \f$ 2 \f$-norm of the residual
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_residual function computes the \f$ 2 \f$-norm of an equation object.
 * It supports Lyapunov Equations as well  as Riccati Equations.
 *
 * \remarks
 * If @mess is compiled with @arpack support the Arnoldi Process from @arpack is used to compute the largest eigenvalue.
 * Otherwise the \ref mess_eigen_arnoldi_template_nrm function is used.
 *
 */
int mess_lrcfadi_residual(mess_equation eqn, mess_options opt, mess_matrix Z, double *nrm){
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input  */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(nrm);
    mess_check_real(Z);

    switch(eqn->eqn_type) {
        case MESS_EQN_LYAP:
        case MESS_EQN_GLYAP:
            ret = mess_lrcfadi_residual_lyap(eqn,opt,Z,nrm);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_residual_lyap);
            break;
        case MESS_EQN_RICCATI:
        case MESS_EQN_GRICCATI:
            ret = mess_lrcfadi_residual_nm(eqn,opt,Z,nrm);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_lrcfadi_residual_nm);
            break;
        case MESS_EQN_NONE:
        default:
            MSG_ERROR("Equation type is not supported.\n");
            return MESS_ERROR_NOSUPPORT;
    }
    return 0;
}



typedef struct res2_data_st {
    mess_equation eqn;
    mess_matrix Z;
    mess_vector t,t1,t2, t3,t4;
    mess_vector rtc,rtz, rt1,rt2, rtb, rt3,rt4;
} res2_data;

static int res2_mvp(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->RHS, x, d->t);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->RHS, d->t, y);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = AZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,x,d->t1);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_NONE, d->t2, d->t3);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);

    ret = mess_vector_axpy(1,d->t3,y);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y += ZZ'A'x
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_HERMITIAN, x, d->t2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->Z, d->t1, y);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

static int res2t_mvp(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->RHS, x, d->t);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->RHS, d->t, y);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = AZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,x,d->t1);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_HERMITIAN, d->t2, d->t3);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1,d->t3,y);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y += ZZ'A'x
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_NONE, x, d->t2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->Z, d->t1, y);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

static int res2g_mvp(void *data, mess_operation_t op,  mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    // y=B*B'*x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->RHS, x, d->t);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->RHS, d->t, y);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y += AZZ'M'x
    mess_equation_E_apply_vector(d->eqn,MESS_OP_HERMITIAN,x,d->t2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,d->t2,d->t1);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_NONE, d->t2, d->t3);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1,d->t3,y);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y += MZZ'A'x
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_HERMITIAN, x, d->t2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_NONE, d->t2, d->t3);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);

    ret = mess_vector_axpy(1,d->t3,y);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    ret = mess_vector_toreal_nowarn(y);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    return 0;
}

static int res2gt_mvp(void *data, mess_operation_t op,  mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2_data * d = (res2_data *) data;
    int ret = 0;

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->RHS, x, d->t);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->RHS, d->t, y);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);

    // y = A'ZZ'Mx
    ret = mess_equation_E_apply_vector(d->eqn,MESS_OP_NONE,x,d->t2);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z,d->t2,d->t1);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z,d->t1,d->t2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_HERMITIAN, d->t2, d->t3);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1,d->t3,y);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y += M'ZZ'Ax
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_NONE, x, d->t2);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->t1);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->t1, d->t2);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_HERMITIAN, d->t2, d->t3);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);

    ret = mess_vector_axpy(1,d->t3,y);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    ret = mess_vector_toreal_nowarn(y);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    return 0;
}

/* Lyapunov residual   */
static int mess_lrcfadi_residual_lyap(mess_equation eqn, mess_options opt, mess_matrix Z, double *nrm){

    MSG_FNAME(__func__);

    res2_data data;
    mess_mvpcall mvpcall;
    int ret =0;
    mess_vector sv = NULL;
    mess_int_t dim = 0;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(nrm);
    mess_check_real(Z);

    dim = mess_equation_dim(eqn);
    if ( eqn->eqn_type != MESS_EQN_LYAP && eqn->eqn_type != MESS_EQN_GLYAP) {
        MSG_ERROR("This function only supports Lyapunov Equations\n");
        return MESS_ERROR_ARGUMENTS;
    }

    if (Z->rows != dim) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    if ( !mess_equation_has_A(eqn)){
        MSG_ERROR("The equation does not provide a Ax.apply function.\n");
        return MESS_ERROR_ARGUMENTS;
    }
    ret = mess_equation_A_pre(eqn);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_pre);

    /*-----------------------------------------------------------------------------
     *  Prepare the arnoldi/lanczos process
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.t1),&(data.t2),&(data.t3),&(data.t),&sv);

    ret = mess_vector_alloc((data.t1), Z->cols,  MESS_REAL);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc((data.t2), dim, MESS_REAL);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc((data.t3), dim, MESS_REAL);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    ret = mess_vector_alloc((data.t),  eqn->RHS->cols, MESS_REAL);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
    data.Z=Z;
    data.eqn=eqn;

    if ( opt->type  == MESS_OP_NONE ) {
        if ( mess_equation_has_E(eqn) ){
            ret = mess_equation_E_pre(eqn);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_pre);
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2g_mvp, &data);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        } else {
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2_mvp, &data);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

        }
    }  else {
        if ( mess_equation_has_E(eqn) ) {
            ret = mess_equation_E_pre(eqn);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_pre);
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2gt_mvp, &data);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        } else {
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2t_mvp, &data);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        }
    }

    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which=MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        //arpack_opt.tol = 0.0;
        arpack_opt.tol = mess_eps();
        arpack_opt.ncv = ARPACK_NCV;
        arpack_opt.maxit = ARNOLDI_ITERATIONS;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev, 2, MESS_REAL);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        //ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, EV);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        ret = mess_eigen_arpack_lanczos_template(mvpcall, 1, arpack_opt, ev, NULL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);

        if (MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        } else {
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);


    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, ARNOLDI_ITERATIONS, sv, nrm);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif

    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    mess_vector_clear(&sv);
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t3));
    mess_vector_clear(&(data.t));
    mess_mvpcall_clear(&mvpcall);
    return 0;
}


typedef struct res2nm_data_st {
    mess_equation eqn;
    mess_matrix Z;
    mess_vector tc,tz, t1,t2, tb, t3,t4;
} res2nm_data;


/*  AX+XA^T-XC^TCX+BB^T = 0.  */
static int res2nm_mvp_none (void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2nm_data *d = ( res2nm_data *) data;
    int ret =0 ;
    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y = BB^Tx
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->B, x, d->tc);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->B, d->tc, y);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);

    // t1 = ZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, x, d->tz);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t1);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t1);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y+= At1
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_NONE, d->t1, d->t2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1, d->t2, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);



    // y+= ZZ'A'x
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_HERMITIAN, x, d->t2);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(1, d->t2, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y-= ZZ'C'C't1
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->C, d->t1, d->tb);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->C, d->tb, d->t2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(-1, d->t2, y);                                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}


/*  A^TX+X^A-XBB^TX+C^TC = 0 */
static int res2nm_mvp_transposed (void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2nm_data *d = ( res2nm_data *) data;
    int ret =0 ;
    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y = C'Cx
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->C, x, d->tc);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->C, d->tc, y);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);

    // t1 = ZZ'x
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, x, d->tz);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t1);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t1);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y+= A't1
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_HERMITIAN, d->t1, d->t2);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1, d->t2, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y+= ZZ'Ax
    ret = mess_equation_A_apply_vector(d->eqn, MESS_OP_NONE, x, d->t2);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(1, d->t2, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y-= ZZ'BB't1
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->B, d->t1, d->tb);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->B, d->tb, d->t2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_vector_axpy(-1, d->t2, y);                                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

/* AXE^T+EXA^T-EXC^TCXiE^T+BB^T = 0.   */
static int res2nmg_mvp_none(void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2nm_data *d = ( res2nm_data *) data;
    int ret =0 ;

    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y = C'Cx
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->B, x, d->tc);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->B, d->tc, y);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);

    // t1 = ZZ'M'x
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_HERMITIAN, x, d->t1);                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t1, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t1);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t1);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y+= A't1
    ret = mess_equation_A_apply_vector(d->eqn,MESS_OP_NONE,d->t1, d->t2);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1, d->t2, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y+= M'ZZ'Ax
    ret = mess_equation_A_apply_vector(d->eqn,MESS_OP_HERMITIAN ,x ,d->t2);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_NONE, d->t2, d->t3);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);
    ret = mess_vector_axpy(1, d->t3, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y-= M'ZZ'BB't1
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->C, d->t1, d->tb);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->C, d->tb, d->t2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_NONE, d->t2, d->t1);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);
    ret = mess_vector_axpy(-1, d->t1, y);                                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}


/* A^TXE+E^TX^A-E^TXBB^TXE+C^TC = 0  */
static int res2nmg_mvp_transposed (void *data, mess_operation_t op, mess_vector x, mess_vector y){
    MSG_FNAME(__func__);
    res2nm_data *d = ( res2nm_data *) data;
    int ret =0 ;

    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y = C'Cx
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->C, x, d->tc);                                    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->C, d->tc, y);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);

    // t1 = ZZ'Mx
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_NONE, x, d->t1);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t1, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t1);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t1);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);

    // y+= A't1
    ret = mess_equation_A_apply_vector(d->eqn,MESS_OP_HERMITIAN,d->t1, d->t2);                  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_vector_axpy(1, d->t2, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);


    // y+= M'ZZ'Ax
    ret = mess_equation_A_apply_vector(d->eqn,MESS_OP_NONE ,x ,d->t2);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_apply_vector);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_HERMITIAN, d->t2, d->t3);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);
    ret = mess_vector_axpy(1, d->t3, y);                                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);

    // y-= M'ZZ'BB't1
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->eqn->B, d->t1, d->tb);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->eqn->B, d->tb, d->t2);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, d->t2, d->tz);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvpt);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Z, d->tz, d->t2);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_toreal_nowarn(d->t2);                                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_equation_E_apply_vector(d->eqn, MESS_OP_HERMITIAN, d->t2, d->t1);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_apply_vector);
    ret = mess_vector_axpy(-1, d->t1, y);                                                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_axpy);
    ret = mess_vector_toreal_nowarn(y);                                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

/*
 *
 * Computes the 2 Norm of the Riccati residual of Z exploiting the
 * symmetry of the residual operator
 *
 * R := A'*Z*Z' + Z*Z'*A' + C'*C - Z*Z'*B*B'*Z*Z'
 *
 * That means it computes the spectral radius of this operator by a
 * power iteration, exploitig the sparsity of A and rectangular
 * structure of Z and B,C. Thus it can be computed in O(n) effort and
 * is therefore much cheaper than the computation of the e.g. the
 * Fobenius norm.
 * \n
 * The Arnoldi Process is controled using
 * \li mess.arnoldi.maxit
 *
 * from the @ref rtm_config.
 * \remarks
 * If @mess is compiled with @arpack support the Arnoldi Process from @arpack is used to compute the largest eigenvalue.
 * Otherwise the \ref mess_eigen_arnoldi_template_nrm function is used.
 */
static int mess_lrcfadi_residual_nm(mess_equation eqn, mess_options opt, mess_matrix Z, double *nrm){
    MSG_FNAME(__func__);
    res2nm_data data;
    mess_mvpcall mvpcall;
    mess_vector sv;
    int ret =0 ;
    mess_int_t dim;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(opt);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(nrm);
    mess_check_real(Z);
    dim = mess_equation_dim(eqn);
    if (Z->rows != dim) {
        MSG_ERROR("Z has the wrong number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    if (!mess_equation_has_A(eqn)) {
        MSG_ERROR("The equation does not provide a Ax.apply function.\n");
        return MESS_ERROR_ARGUMENTS;
    }

    ret = mess_equation_A_pre(eqn);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_A_pre);

    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&(data.tc),&(data.tb), &(data.tz), &(data.t1), &(data.t2), &(data.t3), &sv);

    if (opt -> type == MESS_OP_NONE ) {
        ret = mess_vector_alloc((data.tc), eqn->B->cols, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc((data.tb), eqn->C->rows, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    } else {
        ret = mess_vector_alloc((data.tc), eqn->C->rows, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc((data.tb), eqn->B->cols, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    }
    ret = mess_vector_alloc((data.tz), Z->cols, MESS_REAL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t1), dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t2), dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((data.t3), dim, MESS_REAL);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vetor_init);
    data.eqn=eqn;
    data.Z=Z;

    // mvpcall.data = &data;
    // mvpcall.data_type = MESS_REAL;
    if ( mess_equation_has_E(eqn) ){
        ret = mess_equation_E_pre(eqn);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_equation_E_pre);
        // mvpcall.mvp = res2nmg_mvp;
        if ( opt->type == MESS_OP_NONE ) {
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2nmg_mvp_none, &data);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        } else {
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2nmg_mvp_transposed, &data);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        }

    } else {
        // mvpcall.mvp = res2nm_mvp;
        if ( opt->type == MESS_OP_NONE ) {
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2nm_mvp_none, &data);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        } else {
            ret = mess_mvpcall_operator(&mvpcall, dim, MESS_REAL, res2nm_mvp_transposed, &data);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);
        }
    }

    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which=MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = mess_eps();
        arpack_opt.ncv = ARPACK_NCV;
        arpack_opt.maxit = ARNOLDI_ITERATIONS;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev, 2, MESS_REAL);                                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        //ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        ret = mess_eigen_arpack_lanczos_template(mvpcall, 1, arpack_opt, ev, NULL);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);

        if (MESS_IS_REAL(ev)){
            *nrm = fabs(ev->values[0]);
        } else {
            *nrm = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    ret = mess_eigen_arnoldi_template_nrm(mvpcall, ARNOLDI_ITERATIONS, sv, nrm);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
#endif



    mess_vector_clear(&(data.tc));
    mess_vector_clear(&(data.tz));
    mess_vector_clear(&(data.t1));
    mess_vector_clear(&(data.t2));
    mess_vector_clear(&(data.t3));
    mess_vector_clear(&(data.tb));
    mess_vector_clear(&sv);
    mess_mvpcall_clear(&mvpcall);

    return 0;
}


