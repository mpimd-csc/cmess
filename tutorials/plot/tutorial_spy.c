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
 * @file tutorials/plot/tutorial_spy.c
 * @brief Tutorial for demonstrating the @ref mess_matrix_spy function.
 * @author @mbehr
 * @author @koehlerm
 *
 * # Tutorial Plot Sparsity Pattern of a Matrix
 * In this tutorial we show how to use the @ref mess_matrix_spy function to plot the sparsity pattern of a matrix.
 * We use @ref mess_matrix_read_formated to read a matrix from a @mm file.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/plot/tutorial_spy.c" HEADER
 * \n
 *
 *
 * ## 2. Check Input arguments, Init Matrix and Read Data
 * We check the number of input arguments by using the @c argc variable.
 * We further assume that the first input argument is the path to the @mm file, the second and third one is the
 * width and height of the resulting bitmap plot and the last one is the path to the resulting bitmap file. \n
 * The second step is to initialize the @ref mess_matrix structure using @ref mess_matrix_init and use
 * @ref mess_matrix_read_formated to read directly into @ref MESS_CSR format. We use @c atoi to read the width and height
 * from the command line.
 * @snippet "tutorials/plot/tutorial_spy.c" INITREAD
 * \n
 *
 * ## 3. Generate Spy Plot
 * The calling sequence of @ref mess_matrix_spy is easy.
 * The first input argument is the desired matrix, the second and third one is the width and height of the resulting bitmap image.
 * The last input argument indicates the path to the bitmap image.
 * @snippet "tutorials/plot/tutorial_spy.c" PLOT
 * \n
 *
 *
 * ## 4. Clear Memory
 * We allocated memory by @ref mess_matrix_init and @ref mess_matrix_read_formated, so we have to clear it and use
 * @ref mess_matrix_clear.
 * @snippet "tutorials/plot/tutorial_spy.c" CLEAR
 * \n
 *
 *
 * @sa @ref tutorials
 *
 */

//cond -> code is ignored by doxygen documentation, but snippet references are working


///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc, char **argv){

    /// [DECLARE]
    mess_matrix A;
    mess_int_t w, h;
    /// [DECLARE]

    /// [INITREAD]
    if(argc!=5){
        printf("usage: %s A.mtx width height output.bmp",argv[0]);
        return 1;
    }

    mess_matrix_init(&A);
    mess_matrix_read_formated(argv[1],A,MESS_CSR);
    w = atoi(argv[2]);
    h = atoi(argv[3]);
    /// [INITREAD]

    /// [PLOT]
    mess_matrix_spy(A, argv[4], w, h);
    /// [PLOT]

    /// [CLEAR]
    mess_matrix_clear(&A);
    /// [CLEAR]

    return 0;

}
///@endcond












