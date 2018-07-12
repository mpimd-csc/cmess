/*
* CSCUTILS - A collection of various software routines uses in CSC projects
* Copyright (C) 2018 Martin Koehler
*
* This library is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation; either version 2.1 of the License, or
* (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this library; if not, see <http://www.gnu.org/licenses/>.
*
*/

#ifndef CSC_MTX_H

#define CSC_MTX_H

#ifdef __cplusplus
extern "C" {
#endif
#if __GNUC__ >= 4
#define CSC_ATTR_SCANF(pos1, pos2)   __attribute__ ((format (scanf, pos1, pos2)));
#else
#define CSC_ATTR_SCANF(pos1, pos2)
#endif


    #include <stdio.h>
    #include <stdlib.h>
    #include <stdarg.h>
    #include <sys/stat.h>
    /**
     @file libcscutils/include/cscutils/mtx.h
      @defgroup mtx MatrixMarket File Support

       This part of the library contains routines to write matrix market files.

       \attention This part of the library depends on the \ref error_message module and the \ref io module.

      @addtogroup mtx
      @{
    */

    /**
     * @brief Write a Double Precision Matrix to a Matrix Market File.
     * @param[in] fn    Filename
     * @param[in] m     Number of Rows
     * @param[in] n     Number of Columns
     * @param[in] A     Matrix to write
     * @param[in] lda   Leading dimension of the matrix
     * @return zero on success, -1 otherwise.
     *
     * The csc_mtx_write_double function write a double precision matrix to a mtx file.
     *
     */
    int csc_mtx_write_double(char * fn, int m, int n, double *A, int lda);



    /**
     * @}
     */


#ifdef __cplusplus
};
#endif



#endif /* end of include guard: CSC_IO_H */



