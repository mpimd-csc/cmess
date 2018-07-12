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
 * @file lib/io/write.c
 * @brief Write matrices to files.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "cscutils/io.h"


#define CHECK_NULL(ret) if ((ret) == 0 ) { MSG_ERROR("Writing to file failed.\n"); csc_io_close(f); return MESS_ERROR_FILEIO; }

/**
 * @brief Write a @ref mess_matrix to a @mm file.
 * @param[in] filename   input filename to write matrix
 * @param[in] matrix     input mess_matrix to save to file
 * @return zero on success or a non zero error code
 *
 * The @ref mess_matrix_write function writes a mess_matrix to a @mm file.\n
 * If the file is sparse, the @mm file is in coordinate storage otherwise
 * it is in dense storage. \n
 * If @mess is configured to use @zlib, @bzip or @lzma you can write compressed files by adding
 * .gz, .bz2 or .xz  as file extension.
 *
 * @sa mess_matrix_read
 *
 */
int mess_matrix_write (const char * filename, mess_matrix matrix){
    MSG_FNAME(__func__);
    csc_io_file_t *f;
    mess_matrix work;
    int usework = 0;
    mess_int_t i, j;
    int ret = 0 ;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(filename);

    /*-----------------------------------------------------------------------------
     *  detect compression
     *-----------------------------------------------------------------------------*/
    f  = csc_io_open ( filename, CSC_IO_FILE_WRITE);
    if ( !f ) {
        MSG_ERROR("error opening %s\n", filename);
        return (MESS_ERROR_FILEIO);
    }

    if ( MESS_IS_DENSE(matrix)) {
        work = matrix;
        ret = csc_io_puts("%%MatrixMarket matrix array ", f); CHECK_NULL(ret);
        if ( MESS_IS_REAL ( work)) {
            ret = csc_io_puts("real ", f ); CHECK_NULL(ret);
        } else if ( MESS_IS_COMPLEX(work)) {
            ret = csc_io_puts("complex ", f); CHECK_NULL(ret);
        } else {
            csc_io_close(f);
            MSG_ERROR("unkown datatype: %s\n", mess_datatype_t_str(matrix->data_type));
            return (MESS_ERROR_DATATYPE);
        }
        if ( MESS_IS_GENERAL(work) || MESS_IS_UNSYMMETRIC ( work)) {
            ret = csc_io_puts("general",f );
        } else if ( MESS_IS_HERMITIAN(work)) {
            ret = csc_io_puts("hermitian", f);
        } else if ( MESS_IS_SKEWSYMMETRIC(work)) {
            ret = csc_io_puts("skew-symmetric",f);
        } else if ( MESS_IS_SYMMETRIC(work)) {
            ret = csc_io_puts("symmetric",f);
        } else {
            csc_io_close(f);
            MSG_ERROR("Unknown symmetry code: %s\n",mess_symmetry_t_str(matrix->symmetry) );
            return (MESS_ERROR_SYMMETRIC);
        }
        CHECK_NULL(ret);
        ret = csc_io_puts("\n", f); CHECK_NULL(ret);
        ret = csc_io_printf(f,MESS_PRINTF_INT " " MESS_PRINTF_INT  "\n", work->rows, work->cols); FUNCTION_FAILURE_HANDLE(ret, (ret<=0), csc_io_printf);

        for ( j = 0; j < work->cols; j++) {
            for ( i = 0; i < work->rows; i++){
                if ( MESS_IS_REAL(work)) {
                    ret= csc_io_printf(f,"%20.17e \n", work->values[i+j*work->ld]);
                } else if(MESS_IS_COMPLEX(work)) {
                    ret = csc_io_printf (f, "%20.17e %20.17e\n", creal(work->values_cpx[i+j*work->ld]), cimag(work->values_cpx[i+j*work->ld]));
                }
                FUNCTION_FAILURE_HANDLE(ret, (ret<=0), csc_io_printf);
            }
        }
        csc_io_close(f);
    } else {
        MESS_MATRIX_CHECKFORMAT(matrix, work, usework, MESS_CSR);
        if (usework > 0 ) {
            MSG_ERROR("error converting %s to CSR\n", mess_storage_t_str(matrix->store_type));
            csc_io_close(f);
            return (MESS_ERROR_CONVERT);
        }
        ret = csc_io_puts("%%MatrixMarket matrix coordinate ", f);  CHECK_NULL(ret);

        if ( MESS_IS_REAL ( work))
            ret = csc_io_puts("real ",f);
        else if ( MESS_IS_COMPLEX(work))
            ret = csc_io_puts( "complex ",f);
        else {
            csc_io_close(f);
            MSG_ERROR("unkown datatype: %s\n", mess_datatype_t_str(matrix->data_type));
            return (MESS_ERROR_DATATYPE);
        }
        CHECK_NULL(ret);

        if ( MESS_IS_GENERAL(work) || MESS_IS_UNSYMMETRIC ( work))
            ret = csc_io_puts("general",f);
        else if ( MESS_IS_HERMITIAN(work))
            ret = csc_io_puts( "hermitian", f);
        else if ( MESS_IS_SKEWSYMMETRIC(work))
            ret = csc_io_puts("skew-symmetric",f);
        else if ( MESS_IS_SYMMETRIC(work))
            ret = csc_io_puts( "symmetric",f);
        else {
            csc_io_close (f);
            return (MESS_ERROR_SYMMETRIC);
        }
        ret = csc_io_puts("\n",f );     CHECK_NULL(ret);
        ret = csc_io_printf(f, MESS_PRINTF_INT " " MESS_PRINTF_INT " " MESS_PRINTF_INT "\n", work->rows, work->cols, work->nnz);
        CHECK_NULL(ret);

        for ( i = 0 ; i < work->rows; i++){
            for (j = work->rowptr[i]; j< work->rowptr[i+1]; j++){
                if (MESS_IS_REAL(work)){
                    ret = csc_io_printf(f,  MESS_PRINTF_INT " " MESS_PRINTF_INT " %.15e\n",i+1, work->colptr[j]+1,work->values[j]);
                }
                if (MESS_IS_COMPLEX(work)){
                    ret = csc_io_printf(f , MESS_PRINTF_INT " " MESS_PRINTF_INT " %.15e %.15e\n",i+1, work->colptr[j]+1,creal(work->values_cpx[j]),cimag(work->values_cpx[j]));
                }
                FUNCTION_FAILURE_HANDLE(ret, (ret<=0), csc_io_printf);
            }
        }
        if (usework == 0) {
            mess_matrix_clear(&work);
        }
        csc_io_close(f);
    }
    return (0);
}




