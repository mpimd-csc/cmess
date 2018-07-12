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
 * @file lib/dynsys/stable.c
 * @brief Check if a dynamical system is stable or not.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Check if a dynamical system is stable.
 * @param[in] sys input dynamical system
 * @param[out] stable output which is true if the system is stable, otherwise false
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_dynsys_isstable function checks if a given system is stable or not. \n
 * In case of sparse systems this is done approximately by computing some
 * eigenvalues via an Arnoldi process.
 *
 */
int mess_dynsys_isstable ( mess_dynsys sys, int *stable )
{
    MSG_FNAME(__func__);
    mess_matrix A = NULL;
    mess_matrix E = NULL;
    mess_vector ev = NULL;
    int ret = 0;
    int s = 0;
    mess_int_t i = 0;


    mess_check_nullpointer(sys);
    mess_check_nullpointer(stable);
    *stable = 1 ;

    if (MESS_IS_DYNSYS_LTI(sys)) {
        A = sys->A;
        ret = mess_vector_init(&ev);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev, A->rows, MESS_COMPLEX);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        if ( MESS_IS_DENSE(A)){
            ret = mess_eigen_eig(A, ev, NULL);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
            ret = mess_vector_tocomplex(ev);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = 0; i < ev->dim; i++){
                if ( creal(ev->values_cpx[i])>= 0) s++;
            }
            if ( s > 0) *stable = 0;
            else *stable = 1;
        } else {
            ret = mess_eigen_eigs(A, 50,25, NULL, ev);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eigs);
            ret = mess_vector_tocomplex(ev);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = 0; i < ev->dim; i++){
                if ( creal(ev->values_cpx[i])>= 0) s++;
            }
            if ( s > 0) *stable = 0;
            else *stable = 1;
        }
        mess_vector_clear(&ev);
    } else if (MESS_IS_DYNSYS_GLTI(sys)) {
        A = sys->A;
        E = sys->E;
        ret = mess_vector_init(&ev);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc(ev, A->rows, MESS_COMPLEX);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_allocc);
        if ( MESS_IS_DENSE(A)){
            ret = mess_eigen_eigg(A, E, ev, NULL,NULL,NULL);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eigg);
            ret = mess_vector_tocomplex(ev);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = 0; i < ev->dim; i++){
                if ( creal(ev->values_cpx[i])>= 0) s++;
            }
            if ( s > 0) *stable = 0;
            else *stable = 1;
        } else {
            ret = mess_eigen_eiggs(A,E, 50,25, NULL, ev);           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eiggs);
            ret = mess_vector_tocomplex(ev);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
            for ( i = 0; i < ev->dim; i++){
                if ( creal(ev->values_cpx[i])>= 0) s++;
            }
            if ( s > 0) *stable = 0;
            else *stable = 1;
        }
        mess_vector_clear(&ev);
    } else if (MESS_IS_DYNSYS_2ND(sys)) {
        MSG_ERROR("2nd order systems are not supported yet. \n");
        return MESS_ERROR_DYNSYS;

    } else {
        MSG_ERROR("System type not supported/known.\n");
        return MESS_ERROR_DYNSYS;
    }
    return 0;
}       /* -----  end of function mess_dynsys_isstable  ----- */



/**
 * @brief Check if a matrix is stable.
 * @param[in] A  input matrix
 * @param[out] stable output which is true if \f$ A \f$ is stable, otherwise false
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_matrix_isstable function checks if a matrix is stable. \n
 * If the matrix is sparse this is done approximately by using Arnoldis method.
 *
 */
int  mess_matrix_isstable ( mess_matrix A, int *stable )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_vector ev;
    mess_int_t i = 0 ;

    mess_check_nullpointer(A);
    mess_check_nullpointer(stable);
    mess_check_square(A);

    *stable  = 1;
    ret = mess_vector_init(&ev);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(ev, A->rows, MESS_COMPLEX);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_eigen_eig(A, ev, NULL);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eig);
    ret = mess_vector_tocomplex(ev);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    for ( i =0 ; i < ev->dim; i++){
        if ( creal(ev->values_cpx[i]) >=0 ) {
            *stable = 0;
            break;
        }
    }
    mess_vector_clear(&ev);
    return 0;
}       /* -----  end of function mess_matrix_isstable  ----- */


/**
 * @brief Check if a matrix pair is stable.
 * @param[in] A  input first matrix
 * @param[in] E  input second matrix
 * @param[out] stable output which is true if \f$ (A,E) \f$ is stable, otherwise false
 * @return zero on succes or a non-zero error value
 *
 * The @ref mess_matrix_isstableg function checks if the pair \f$ (A,E) \f$ is stable.\n
 * If it is sparse, the check is done approximately by using Arnoldis method.
 *
 */
int  mess_matrix_isstableg ( mess_matrix A, mess_matrix E, int *stable )
{
    MSG_FNAME(__func__);
    int ret = 0;
    mess_vector ev;
    mess_int_t i = 0 ;

    mess_check_nullpointer(A);
    mess_check_nullpointer(stable);
    mess_check_square(A);

    *stable  = 1;
    ret = mess_vector_init(&ev);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(ev, A->rows, MESS_COMPLEX);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_eigen_eigg(A, E, ev, NULL, NULL, NULL);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_eigen_eigg);
    ret = mess_vector_tocomplex(ev);                            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    for ( i =0 ; i < ev->dim; i++){
        if ( creal(ev->values_cpx[i]) >=0 ) {
            *stable = 0;
            break;
        }
    }
    mess_vector_clear(&ev);
    return 0;

}       /* -----  end of function mess_matrix_isstableg  ----- */

