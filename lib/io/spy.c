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
 * @file lib/io/spy.c
 * @brief Generate spy plots of matrices.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/image.h"
#include <complex.h>

/**
 * @brief Generate a bitmap file from the non zero structure of a matrix.
 * @param[in] matrix   input matrix to visualize
 * @param[in] filename  input filename of output bitmap
 * @param[in] width  input width of output bitmap
 * @param[in] height  input height of output bitmap
 * @return zero on success or a non zero error code.
 *
 * The function @ref mess_matrix_spy produces a bitmap of the non zero structure of a matrix. \n
 * The width and height argument are changed if the dimension of the matrix do not
 * fit good in the bitmap.
 *
 */
int mess_matrix_spy(mess_matrix matrix, const char* filename,int width, int height){
    csc_image_bmp pic;
    MSG_FNAME(__func__);
    mess_matrix work;
    double dx, dy;
    int wx, wy;
    int x, y;
    mess_int_t i, j;
    int conv = 0;

    mess_check_nullpointer(matrix);
    mess_check_nullpointer(filename);

    if ( height <= 0 || width <= 0){
        MSG_ERROR("image size is wrong\n");
        return MESS_ERROR_ARGUMENTS;
    }

    MESS_MATRIX_CHECKFORMAT(matrix, work, conv, MESS_CSR);
    dx = ((double) matrix->cols) / width;
    dy = ((double) matrix->rows)/ height;
    if ( dx > 1 ){
        dx = (int) (dx+1);
        wx = 1;
    } else {
        wx = (int)(1.0/dx);
    }
    if ( dy > 1 ){
        dy = (int) (dy+1);
        wy = 1;
    } else {
        wy = (int)(1.0/dy);
    }
    width = (int) (matrix->cols/dx);
    height = (int) (matrix->rows/dy);

    MSG_INFO("filename = %s\n", filename);
    MSG_INFO("height   = %d\n", height);
    MSG_INFO("width    = %d\n", width);
    MSG_INFO("dx       = %g\n", dx);
    MSG_INFO("dy       = %g\n", dy);
    MSG_INFO("wx       = %d\n", wx);
    MSG_INFO("wy       = %d\n", wy);
    MSG_INFO("cols     = " MESS_PRINTF_INT "\n", matrix->cols);
    MSG_INFO("rows     = " MESS_PRINTF_INT "\n", matrix->rows);
    csc_image_bmp_init(&pic, width, height);

    for ( i = 0; i < matrix->rows; i++){
        for (j =  matrix->rowptr[i]; j < matrix->rowptr[i+1]; j++){
            x = (int)(matrix->colptr[j]/dx);
            y = (int)(i / dy);
            csc_image_bmp_fill_rect ( pic, x, y,wx, wy,  255,0 , 0);
        }
    }

    if ( csc_image_bmp_write  (filename, pic) ){
        MSG_ERROR( "An error while writing the bitmap!\n");
    }
    csc_image_bmp_clear(&pic);
    if (conv == 0){
        mess_matrix_clear(&work);
    }

    return 0;
}


