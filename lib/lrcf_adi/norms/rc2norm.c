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
 * @file lib/lrcf_adi/norms/rc2norm.c
 * @brief Compute the \f$ 2 \f$-norm difference between two factored solutions of a matrix equation.
 * @author @koehlerm
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/eigenvalue.h"
#include "mess/error_macro.h"
#include <complex.h>

#define IT_MAX  20

struct rc2data {
    mess_matrix Z;
    mess_matrix Zold;
    mess_vector x1,x2;
};

static int rc2mvp(void *data, mess_operation_t op,  mess_vector x, mess_vector y) {
    MSG_FNAME(__func__);
    struct rc2data * d = (struct rc2data*) data;
    int ret = 0;
    ret = mess_vector_toreal_nowarn(y);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Zold, x, d->x1); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_mvp(MESS_OP_NONE,d->Zold,  d->x1,y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_vector_scale(-1.0,y);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_scale);

    ret = mess_matrix_mvp(MESS_OP_HERMITIAN,d->Z, x, d->x2);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
    ret = mess_matrix_gaxpy(MESS_OP_NONE, d->Z, d->x2, y);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_gaxpy);

    ret = mess_vector_toreal_nowarn(y);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    return 0;
}

/**
 * @brief Compute \f$ 2 \f$-norm difference between two factored solutions of a matrix equation.
 * @param[in] Zold input       first factor
 * @param[in] Z input input             second factor
 * @param[out] chg  \f$ 2 \f$ norm difference
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_fac2diff function computes the \f$ 2 \f$-norm difference of
 * \f$   |\underbrace{Z_{old}Z_{old}^T}_{X_{old}} - \underbrace{ZZ^T}_{X} |_2 \f$
 * without forming \f$ X_{old} \f$ and \f$ X \f$.
 *
 */
int  mess_lrcfadi_fac2diff ( mess_matrix Zold, mess_matrix Z, double *chg )
{
    MSG_FNAME(__func__);
    struct rc2data dat;
    mess_mvpcall mvpcall;
    mess_vector sv = NULL;
    int ret =0;

    mess_check_nullpointer(Zold);
    mess_check_nullpointer(Z);
    mess_check_nullpointer(chg);

    if ( Z->rows != Zold->rows ) {
        MSG_ERROR("Z and Zold must have the same number of rows. \n");
        return MESS_ERROR_DIMENSION;
    }

    MESS_INIT_VECTORS(&(dat.x1),&(dat.x2));
    ret = mess_vector_alloc((dat.x1), Zold->cols, Zold->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc((dat.x2), Z->cols, Z->data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);

    dat.Z=Z;
    dat.Zold=Zold;

    ret = mess_mvpcall_operator(&mvpcall, Z->rows, MESS_REAL, rc2mvp, &dat); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_mvpcall_operator);

    ret = mess_vector_init(&sv);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(sv, Z->rows, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_ones(sv); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_ones);

#ifdef MESS_HAVE_ARPACK
    {
        mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
        mess_vector ev;
        arpack_opt.which=MESS_EIGEN_ARPACK_LARGE_MAGNITUDE;
        arpack_opt.tol = 0.0;
        arpack_opt.ncv = 15;
        arpack_opt.maxit = 60;
        arpack_opt.b0 = sv;
        ret = mess_vector_init(&ev);                                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev, 2, MESS_REAL);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_eigen_arpack_template(mvpcall, 1, arpack_opt, ev, NULL);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arpack_template);
        if (MESS_IS_REAL(ev)){
            *chg = fabs(ev->values[0]);
        } else {
            *chg = cabs(ev->values_cpx[0]);
        }
        mess_vector_clear(&ev);
    }
#else
    {
        ret = mess_eigen_arnoldi_template_nrm(mvpcall, IT_MAX, sv, chg);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_arnoldi_template_nrm);
    }
#endif
    mess_vector_clear(&sv);
    mess_vector_clear(&(dat.x1));
    mess_vector_clear(&(dat.x2));
    mess_mvpcall_clear(&mvpcall);

    return 0;
}       /* -----  end of function mess_lrcfadi_rc2norm  ----- */

