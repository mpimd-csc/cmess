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
 * @file lib/matrix/rand.c
 * @brief Generate different random matrices.
 * @author @koehlerm
 * @author @mbehr
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>
#include <time.h>
#include "blas_defs.h"



/**
 * @brief Generate a positive random floating point number between -1 and 1.
 * @return A positive random number between -1 and 1
 *
 * The @ref __drand function generates an uniformly distributed random number between
 * -1 and 1.
 *
 *  @attention @ref mess_matrix_rand_coord depends on the scaling to [-1;1].
 */
static double __drand()
{
    double ret = -1.0 + 2.0*( (double) rand() / ( (double) RAND_MAX ) );
    return ret;
}

/**
 * @brief Init random number generator with given seed.
 * @param[in] seed   pointer to @ref mess_int_t for given seed, if @c NULL current calender time is used.
 * @return zero in every case
 *
 * The @ref mess_matrix_rand_init initializes the @c rand  function from @c libc for generating random numbers. \n
 * If @p seed points to @c NULL the current calender time is used for initialization of @c rand function.
 * This function is usefull in combination with:
 * @sa mess_matrix_rand
 * @sa mess_matrix_rand_dense
 * @sa mess_matrix_rand_csc
 * @sa mess_matrix_rand_csr
 * @sa mess_matrix_rand_coord
 *
 */
int mess_matrix_rand_init(mess_int_t *seed){

    if(seed){
        srand( (unsigned) (*seed));
    }else{
        time_t t;
        srand((unsigned) time(&t));
    }
    return 0;
}

#define IDX(i,j,LDA) ((j)*(LDA)+(i))

/**
 * @brief Generate a @ref MESS_REAL @ref MESS_COMPLEX random dense matrix.
 * @param[out] mat   output random matrix
 * @param[in] rows   input number of rows
 * @param[in] cols   input number of cols
 * @param[in] dt     input desired dataype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rand_dense function generates a @ref MESS_REAL or @ref MESS_COMPLEX random dense matrix \f$ mat \f$. \n
 * The matrix is allocated by \ref mess_matrix_alloc and gets the leading dimension from it. \n
 * The matrix is filled with values generated by @c rand  from @c libc.
 *
 * @sa mess_matrix_rand_dense_uniform
 *
 */
int mess_matrix_rand_dense ( mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_datatype_t dt )
{
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int ret = 0 ;

    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);
    mess_check_datatype(dt);

    MESS_MATRIX_RESET(mat);
    //srand(time(NULL));

    ret = mess_matrix_alloc(mat, rows, cols, rows*cols, MESS_DENSE, dt);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    if(dt==MESS_REAL){
        for ( i = 0; i< rows ; i++ ){
            for (j = 0 ; j < cols; j++) {
                mat->values[IDX(i,j,mat->ld)] = __drand();
            }
        }
    }else{
        for ( i = 0; i< rows ; i++ ){
            for (j = 0 ; j < cols; j++) {
                mat->values_cpx[IDX(i,j,mat->ld)] = __drand()+I*(__drand());
            }
        }

    }

    return 0;
}       /* -----  end of function mess_matrix_rand_dense  ----- */

/**
 * @brief Generate a random @ref MESS_REAL dense matrix with a uniform value distribution \f$ (-1,1) \f$.
 * @param[out] mat          output random matrix
 * @param[in,out] seed      input/output initial seed for the random number generator
 * @param[in] rows          input number of rows
 * @param[in] cols          input number of cols
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rand_dense_uniform function generates a random dense matrix \f$ mat \f$ with a uniform
 * value distribution between \f$ (-1,1) \f$ generated by @lapack @c DLANRV. \n
 * The \f$ seed \f$ parameter is an integer array of length \f$ 4 \f$.
 * The value of \f$ seed[3] \f$ need to be an odd value, if not it will
 * be adjusted. The \f$ seed \f$ is updated during the call of @c DLANRV.
 *
 * @sa mess_matrix_rand_dense
 *
 */
int mess_matrix_rand_dense_uniform( mess_matrix mat, mess_int_t *seed , mess_int_t rows, mess_int_t cols )
{
    MSG_FNAME(__func__);
    int ret = 0 ;
    mess_int_t one = 1;
    mess_int_t n2;

    mess_check_nullpointer(mat);
    mess_check_nullpointer(seed);
    mess_check_positive(rows);
    mess_check_positive(cols);
    MESS_MATRIX_RESET(mat);

    ret = mess_matrix_alloc(mat, rows, cols, rows*cols, MESS_DENSE, MESS_REAL);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if ( seed[3] % 2 == 0 ) seed[3]++;
    seed[0] %= 4096;
    seed[1] %= 4096;
    seed[2] %= 4096;
    seed[3] %= 4096;
    n2 = cols*mat->ld;
    F77_GLOBAL(dlarnv,DLARNV)( &one, seed, &n2, mat->values );
    return 0;
}       /* -----  end of function mess_matrix_rand_dense  ----- */


/**
 * @brief Generate a @ref MESS_REAL/@ref MESS_COMPLEX random @ref MESS_CSR matrix.
 * @param[out] mat  output random matrix
 * @param[in] rows  input rows of the matrix
 * @param[in] cols  input cols of the matrix
 * @param[in] p     input percentage of non zero elements
 * @param[in] dt    input the desired datatype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rand_csr function generates a @ref MESS_REAL or @ref MESS_COMPLEX random @ref MESS_CSR
 * matrix with approximately \f$ p\cdot cols \f$ entries per row.
 *
 */
int mess_matrix_rand_csr ( mess_matrix mat, mess_int_t rows, mess_int_t cols , double p, mess_datatype_t dt )
{
    MSG_FNAME(__func__);
    mess_int_t nnz;
    mess_int_t i, j;
    mess_int_t nz = 0;
    mess_int_t nzpl = 0;
    int ret = 0;

    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);
    mess_check_datatype(dt);
    mess_check_true(0<=p && p<=1);

    MESS_MATRIX_RESET(mat);

    nnz = (rows*cols*p);
    ret  = mess_matrix_alloc(mat, rows, cols, nnz, MESS_CSR, dt);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    mat->rowptr[0] = 0;
    for ( i = 0; i < rows; i ++) {
        nzpl = cols * p ;
        if ( nzpl == 0) {
            nzpl = 1;
        }
        if ( nz + nzpl >= mat->nnz) {
            if(dt==MESS_REAL){
                mess_try_realloc(mat->values, double *, (nz+nzpl)*sizeof(double));
            }else{
                mess_try_realloc(mat->values_cpx, mess_double_cpx_t *, (nz+nzpl)*sizeof(mess_double_cpx_t));
            }
            mess_try_realloc(mat->colptr, mess_int_t *, (nz+nzpl)*sizeof(mess_int_t));
            mat->nnz = nz+nzpl;
        }
        if(dt==MESS_REAL){
            for ( j = 0; j < nzpl; j++){
                mat->values[nz] = __drand();
                mat->colptr[nz++] = (rand()+1)%cols;
            }

        }else{
            for ( j = 0; j < nzpl; j++){
                mat->values_cpx[nz] = __drand()+I*(__drand());
                mat->colptr[nz++] = (rand()+1)%cols;
            }

        }
        mat->rowptr[i+1] = nz;
    }


    ret = mess_matrix_sort(mat);
    FUNCTION_FAILURE_HANDLE(ret, ( ret!=0), mess_matrix_sort);
    ret = mess_matrix_dupl(mat);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_dupl);


    return 0;
}       /* -----  end of function mess_matrix_rand_csr  ----- */

/**
 * @brief Generate a @ref MESS_REAL/@ref MESS_COMPLEX random @ref MESS_CSC matrix.
 * @param[out] mat  output random matrix
 * @param[in] rows  input rows of the matrix
 * @param[in] cols  input cols of the matrix
 * @param[in] p     input percentage of non zero elements
 * @param[in] dt    input the desired datatype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rand_csc function generates a @ref MESS_REAL or @ref MESS_COMPLEX random @ref MESS_CSC
 * matrix with approximately \f$ p \cdot rows \f$ entries per column.
 *
 */
int mess_matrix_rand_csc ( mess_matrix mat, mess_int_t rows, mess_int_t cols , double p, mess_datatype_t dt )
{
    MSG_FNAME(__func__);
    mess_int_t nnz;
    mess_int_t i, j;
    mess_int_t nz = 0;
    mess_int_t nzpl = 0;
    int ret = 0;

    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);
    mess_check_datatype(dt);
    mess_check_true(0<=p && p<=1);

    MESS_MATRIX_RESET(mat);

    nnz = (rows*cols*p);
    ret  = mess_matrix_alloc(mat, rows, cols, nnz, MESS_CSC, dt);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    mat->colptr[0] = 0;

    for ( i = 0; i < cols; i ++) {
        nzpl = rows * p ;
        if ( nzpl == 0) {
            nzpl = 1;
        }

        if ( nz + nzpl >= mat->nnz) {
            if(dt==MESS_REAL){
                mess_try_realloc(mat->values, double *, (nz+nzpl)*sizeof(double));
                mess_try_realloc(mat->rowptr, mess_int_t *, (nz+nzpl)*sizeof(mess_int_t));
            }else{
                mess_try_realloc(mat->values_cpx, mess_double_cpx_t *, (nz+nzpl)*sizeof(mess_double_cpx_t));
                mess_try_realloc(mat->rowptr, mess_int_t *, (nz+nzpl)*sizeof(mess_int_t));
            }
            mat->nnz = nz+nzpl;
        }
        if(dt==MESS_REAL){

            for ( j = 0; j < nzpl; j++){
                mat->values[nz] = __drand();
                mat->rowptr[nz++] = (rand()+1)%rows;
            }
        }else{
            for ( j = 0; j < nzpl; j++){
                mat->values_cpx[nz] = __drand()+I*(__drand());
                mat->rowptr[nz++] = (rand()+1)%rows;
            }
        }
        mat->colptr[i+1] = nz;
    }
    ret = mess_matrix_sort(mat);    FUNCTION_FAILURE_HANDLE(ret, ( ret!=0), mess_matrix_sort);
    ret = mess_matrix_dupl(mat);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_dupl);


    return 0;
}       /* -----  end of function mess_matrix_rand_csc  ----- */


/**
 * @brief Generate a @ref MESS_REAL/@ref MESS_COMPLEX random COORD matrix.
 * @param[out] mat  output random matrix
 * @param[in] rows  input rows of the matrix
 * @param[in] cols  input cols of the matrix
 * @param[in] p     input percentage of non zero elements
 * @param[in] dt    input the desired datatype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rand_coord function generates a @ref MESS_REAL or @ref MESS_COMPLEX random Compressed Sparse Row
 * matrix with approximately \f$ p\cdot cols \f$ entries per row.
 *
 * @sa mess_matrix_rand_dense
 * @sa mess_matrix_rand_csr
 * @sa mess_matrix_rand_csc
 *
 */
int mess_matrix_rand_coord ( mess_matrix mat, mess_int_t rows, mess_int_t cols , double p, mess_datatype_t dt )
{
    MSG_FNAME(__func__);
    mess_int_t nnz;
    mess_int_t i, j;
    mess_int_t nz = 0;
    int ret = 0;

    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);
    mess_check_datatype(dt);
    mess_check_true(0<=p && p<=1);

    MESS_MATRIX_RESET(mat);

    nnz = (rows*cols*p);
    ret  = mess_matrix_alloc(mat, rows, cols, nnz, MESS_COORD, dt); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    if (MESS_IS_REAL(mat)){
        for ( i = 0; i < mat->rows; ++i){
            for( j=0; j< mat->cols; ++j){
                //throw a "dice" if entry should be filled, scale __drand() to [0,1]
                if( (__drand() + 1) /2 <= p  && nz< nnz){
                    mat->values[nz] = __drand();
                    mat->rowptr[nz] = i;
                    mat->colptr[nz] = j;
                    ++nz;
                }
            }
        }
    }else{
        for ( i=0; i < mat->rows; ++i){
            for( j=0; j< mat->cols; ++j){
                //throw a "dice" if entry should be filled, scale __drand() to [0,1]
                if( (__drand() + 1) /2 <= p && nz<nnz){
                    mat->values_cpx[nz] = __drand()+ I*(__drand());
                    mat->rowptr[nz] = i;
                    mat->colptr[nz] = j;
                    ++nz;
                }
            }
        }
    }

    mat->nnz = nz;

    return 0;
}       /* -----  end of function mess_matrix_rand_coord  ----- */



/**
 * @brief Wrapper around @ref  mess_matrix_rand_csr, @ref mess_matrix_rand_csc, @ref mess_matrix_rand_dense and @ref mess_matrix_rand_coord.
 * @param[out] mat          output random matrix
 * @param[in] rows          input number of rows
 * @param[in] cols          input number of cols
 * @param[in] storetype     input storage type of the output matrix
 * @param[in] p             input fillrate for the rows or columns of a sparse matrix
 * @param[in] dt            input the desired datatype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_rand function calls depending on the @p storetype and @p dt argument the corresponding
 * random function for sparse or dense, @ref MESS_REAL or @ref MESS_COMPLEX matrices. \n
 * The @p p  parameter is ignored in the case of a dense matrix.
 */
int mess_matrix_rand(mess_matrix mat, mess_int_t rows, mess_int_t cols, mess_storage_t storetype, mess_datatype_t dt, double p)
{
    MSG_FNAME(__func__);


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(mat);
    mess_check_positive(rows);
    mess_check_positive(cols);
    mess_check_datatype(dt);


    switch ( storetype ) {
        case MESS_DENSE:
            return mess_matrix_rand_dense(mat, rows, cols,dt);
            break;
        case MESS_CSR:
            return mess_matrix_rand_csr(mat, rows, cols, p,dt);
            break;
        case MESS_CSC:
            return mess_matrix_rand_csc(mat, rows, cols, p,dt);
            break;
        case MESS_COORD:
            return mess_matrix_rand_coord(mat, rows, cols, p, dt);
            break;
        default:
            MSG_ERROR("no random function for type %u available.\n", storetype);
            return MESS_ERROR_STORAGETYPE;
            break;
    }               /* -----  end switch  ----- */

    return 0;
}
