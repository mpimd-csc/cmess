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
 * @addtogroup test_direct
 * @{
 * @file tests/direct/check_det.c
 * @brief Check the det and detc functions of the solver.
 * @author  @mbehr
 * @test
 * Check the @ref mess_direct_determinant and @ref mess_direct_determinantc functions.
 * @}
 */

#include "../call_macro.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"


int test_det(mess_direct sol, double* det_check, mess_double_cpx_t *detc_check){
    double eps = sqrt(mess_eps());
    double mr, mi, e;
    double abs_err, rel_err;
    double det;
    mess_double_cpx_t detc;
    int ret;

    if(MESS_IS_REAL(sol)){
        CALL(mess_direct_determinant(sol,&mr,&e));
        det = mr*pow(2,e);
        abs_err = fabs(det-(*det_check)); rel_err = abs_err/fabs((*det_check));
        if(rel_err > eps){
            printf("FAILED:%s\n",sol->name);
            printf("det        = %e\n", det);
            printf("det check  = %e\n", (*det_check));
            printf("abs. error = %e\n", abs_err);
            printf("rel. error = %e\n", rel_err);
            return ret+1;
        }else{
            printf("PASSED:%s\n",sol->name);
        }
    }else{
        CALL(mess_direct_determinantc(sol,&mr,&mi,&e));
        detc = mr*pow(2,e)+mi*pow(2,e)*I;
        abs_err = cabs(detc-(*detc_check)); rel_err = abs_err/cabs((*detc_check));
        if(rel_err > eps){
            printf("FAILED:%s\n",sol->name);
            printf("detc       = %e + %e*I\n",creal(detc), cimag(detc));
            printf("detc check = %e + %e*I\n",creal(*detc_check), cimag(*detc_check));
            printf("abs. error = %e\n",abs_err);
            printf("rel. error = %e\n",rel_err);
            return ret+1;
        }else{
            printf("PASSED:%s\n",sol->name);
        }
    }
    return ret;
}

int main(int argc, char ** argv){
    mess_version();
    mess_error_level = 3;
    mess_matrix dense,csc;
    mess_direct sol;
    double det_check;
    mess_double_cpx_t detc_check, scal = cos(rand()) + sin(rand())*I;
    int err=0, ret=0;

    /*-----------------------------------------------------------------------------
     *  get parameters
     *-----------------------------------------------------------------------------*/
    if ( argc != 3) {
        printf("usage: %s matrix.mtx det\n", argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  init and read data
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&dense,&csc);
    CALL(mess_matrix_read_formated(argv[1],dense,MESS_DENSE));
    CALL(mess_matrix_convert(dense, csc, MESS_CSC));
    det_check = atof(argv[2]);
    detc_check = cpow(scal,dense->rows)*det_check;

    /*-----------------------------------------------------------------------------
     *  test real
     *-----------------------------------------------------------------------------*/
    //lapack
    CALL(mess_direct_init(&sol)); CALL(mess_direct_create_lapack_lu(dense,sol));
    err += test_det(sol,&det_check,NULL);
    CALL(mess_direct_clear(&sol));

#ifdef MESS_HAVE_UMFPACK
    //umfpack
    CALL(mess_direct_init(&sol)); CALL(mess_direct_create_umfpack(csc,sol));
    err += test_det(sol,&det_check,NULL);
    CALL(mess_direct_clear(&sol));
#endif

    /*-----------------------------------------------------------------------------
     *  test complex
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_tocomplex(dense)); CALL(mess_matrix_scalec(scal, dense));
    CALL(mess_matrix_tocomplex(csc));   CALL(mess_matrix_scalec(scal, csc));

    //lapack
    CALL(mess_direct_init(&sol));CALL(mess_direct_create_lapack_lu(dense,sol));
    err += test_det(sol,NULL, &detc_check);
    CALL(mess_direct_clear(&sol));

#ifdef MESS_HAVE_UMFPACK
     //umfpack
    CALL(mess_direct_init(&sol)); CALL(mess_direct_create_umfpack(csc,sol));
    err += test_det(sol,NULL, &detc_check);
    CALL(mess_direct_clear(&sol));
#endif

    /*-----------------------------------------------------------------------------
     *  clear data
     *-----------------------------------------------------------------------------*/
    MESS_CLEAR_MATRICES(&dense, &csc);

    return err;

}
