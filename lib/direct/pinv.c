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
 * @file lib/direct/pinv.c
 * @brief Compute Pseudoinverse of a dense matrix.
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/misc.h"
#include "mess/error_macro.h"
#include "blas_defs.h"
#include <complex.h>
#include <math.h>





/**
 * @brief Compute the pseudoinverse of a dense matrix \f$A\f$.
 * @param[in] A   input \f$ A \f$
 * @param[out] Pinv  output \f$ P_{inv} \f$ pseudoinverse of \f$A\f$
 *
 * The @ref mess_matrix_pinv computes the pseudoinverse of a
 * dense matrix using @lapack @c dgels or @c zgels methods.
 *
 * @see mess_direct_inverse
 * @see mess_eigen_svd
 */
int mess_matrix_pinv ( mess_matrix A, mess_matrix Pinv )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_int_t m = A->rows, n = A->cols, mn = m>n?m:n,i,lwork,info;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_real_or_complex(A);
    mess_check_dense(A);
    mess_check_nullpointer(Pinv);


    /*-----------------------------------------------------------------------------
     *  call lapack {dz}gels
     *-----------------------------------------------------------------------------*/
    mess_matrix Acopy;
    mess_matrix_init(&Acopy);
    ret = mess_matrix_copy(A,Acopy);                                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);

    MESS_MATRIX_RESET(Pinv);
    ret = mess_matrix_alloc(Pinv,mn,m,0,MESS_DENSE,A->data_type);   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_alloc);
    if(MESS_IS_REAL(A)){
        //upper block identity, and lower block zero for m>n
        for(i=0;i<m;++i){Pinv->values[i+(Pinv->ld)*i]=1;}

        //workspace query and real work
        double work,*workspace;
        lwork=-1;
        F77_GLOBAL(dgels,DGELS)("N",&(m),&(n),&(Pinv->cols),Acopy->values,&(Acopy->ld), Pinv->values,&(Pinv->ld),&(work),&(lwork),&(info));
        lwork = nearbyint(work+1);
        mess_try_alloc(workspace,double*,sizeof(double)*lwork);
        F77_GLOBAL(dgels,DGELS)("N",&(m),&(n),&(Pinv->cols),Acopy->values,&(Acopy->ld), Pinv->values,&(Pinv->ld),workspace,&(lwork),&(info));
        mess_free(workspace);
    }else{
        //upper block identity, and lower block zero for m>n
        for(i=0;i<m;++i){Pinv->values_cpx[i+(Pinv->ld)*i]=1;}

        //workspace query and real work
        mess_double_cpx_t work_cpx, *workspace_cpx;
        lwork=-1;
        F77_GLOBAL(zgels,ZGELS)("N",&(m),&(n),&(Pinv->cols),Acopy->values_cpx,&(Acopy->ld), Pinv->values_cpx,&(Pinv->ld),&(work_cpx),&(lwork),&(info));
        lwork = nearbyint(creal(work_cpx)+1);
        mess_try_alloc(workspace_cpx,mess_double_cpx_t*,sizeof(mess_double_cpx_t)*lwork);
        F77_GLOBAL(zgels,ZGELS)("N",&(m),&(n),&(Pinv->cols),Acopy->values_cpx,&(Acopy->ld), Pinv->values_cpx,&(Pinv->ld),workspace_cpx,&(lwork),&(info));
        mess_free(workspace_cpx);
    }
    mess_matrix_clear(&Acopy);
    ret = mess_matrix_resize(Pinv,n,m);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_resize);

    return ret;
}       /* -----  end of function mess_matrix_pinv  ----- */

