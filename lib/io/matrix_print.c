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
 * @file lib/io/matrix_print.c
 * @brief Print matrices and underlying information to the standard output.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>


/**
 * @brief Display a @ref mess_matrix to stdout.
 * @param[in] matrix  input the matrix to display
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_display function prints a complete matrix to the standard output. \n
 * In contrast to @ref mess_matrix_print the output is in a rectangular and better readable form.
 * The function uses the @c csc_io_get_term_width function to retrieve the width of the terminal.
 * In case of a sparse matrix @ref mess_matrix_print is called.
 * We use zero-based counting for indices.
 *
 * @sa mess_matrix_printshort
 * @sa mess_matrix_printinfo
 * @sa mess_matrix_printdata
 *
 */
int mess_matrix_display(mess_matrix matrix){
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int term_width, length_token = 0, used_chars=0;
    const char token_delimiter [] = "    ";
    const int length_token_delimiter = 4; // lenght of token_delimiter

    /*-----------------------------------------------------------------------------
     * check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(matrix);
    mess_check_real_or_complex(matrix);

    /*-----------------------------------------------------------------------------
     *  prepare output
     *-----------------------------------------------------------------------------*/
    term_width = csc_get_term_width(); //get width of terminal

    /*-----------------------------------------------------------------------------
     *  print general information
     *-----------------------------------------------------------------------------*/
    mess_matrix_printinfo(matrix);

    /*-----------------------------------------------------------------------------
     *  print a matrix
     *-----------------------------------------------------------------------------*/
    //nice print only for dense matrices, otherwise call standard mess_matrix_print
    if(!MESS_IS_DENSE(matrix)){
        return mess_matrix_print(matrix);
    }

    for (i = 0; i < matrix->rows; i++){
        for (j = 0; j < matrix->cols; j++){
            if(MESS_IS_REAL(matrix)){
                length_token = mess_print_format_double(matrix->values[i+j*matrix->ld]);
            }else{
                length_token = mess_print_format_double_cpx(matrix->values_cpx[i+j*matrix->ld]);
            }
            used_chars += length_token;

            //break if last col was printed
            if(j==(matrix->cols-1)){break;}

            // check if enough space for next token in terminal line
            if(term_width - (used_chars + 2*length_token + 4*length_token_delimiter) <= 0){
                // not enough place for all columns print ellipsis and last value
                if(j<matrix->cols-2){
                    MSG_PRINT("%s...%s", token_delimiter, token_delimiter);
                }else{
                    MSG_PRINT("%s", token_delimiter);
                }
                if(MESS_IS_REAL(matrix)){
                    length_token = mess_print_format_double(matrix->values[i+(matrix->cols-1)*matrix->ld]);
                }else{
                    length_token = mess_print_format_double_cpx(matrix->values_cpx[i+(matrix->cols-1)*matrix->ld]);
                }
                break;
            }else{
                MSG_PRINT("%s", token_delimiter); // print delimiter
                used_chars += length_token_delimiter;
            }
        }
        MSG_PRINT("\n");
        used_chars = 0;
    }
    return 0;
}



/**
 * @brief Print a  @ref mess_matrix to stdout.
 * @param[in] matrix  input the matrix to print out
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_print function prints a complete matrix to the standard output. \n
 * The output format is almost the same as in @matlab .\n
 * In the case of sparse matrices only the non zero entries are printed.
 * We use zero-based counting for indices.
 *
 * @sa mess_matrix_printshort
 * @sa mess_matrix_printinfo
 * @sa mess_matrix_printdata
 * @sa mess_matrix_display
 *
 */
int mess_matrix_print(mess_matrix matrix){
    MSG_FNAME(__func__);
    mess_int_t i,j;

    mess_check_nullpointer(matrix);

    if (MESS_IS_DENSE(matrix)) {
        MSG_PRINT ("size = ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " ) LD = " MESS_PRINTF_INT " \n", matrix->rows, matrix->cols, matrix->ld);
    } else {
        MSG_PRINT ("size = ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " ) \n", matrix->rows, matrix->cols);
    }
    /*-----------------------------------------------------------------------------
     *  print a user defined matrix
     *-----------------------------------------------------------------------------*/

    if ( MESS_IS_CSR ( matrix )) {
        for ( i = 0 ; i < matrix->rows; i++){
            for (j = matrix->rowptr[i]; j<matrix->rowptr[i+1]; j++){
                if ( MESS_IS_REAL(matrix)) {
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e \n", i, matrix->colptr[j], matrix->values[j]);
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, matrix->colptr[j]);
                    mess_print_format_double(matrix->values[j]); MSG_PRINT("\n");
                } else if(MESS_IS_COMPLEX(matrix)){
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e + %10e i\n", i, matrix->colptr[j], creal(matrix->values_cpx[j]), cimag(matrix->values_cpx[j]));
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, matrix->colptr[j]);
                    mess_print_format_double_cpx(matrix->values_cpx[j]); MSG_PRINT("\n");
                }
            }
        }
    } else if ( MESS_IS_CSC ( matrix )){
        for ( i = 0 ; i < matrix->cols; i++){
            for (j = matrix->colptr[i]; j<matrix->colptr[i+1]; j++){
                if ( MESS_IS_REAL(matrix)) {
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e \n", matrix->rowptr[j], i, matrix->values[j]);
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", matrix->rowptr[j], i);
                    mess_print_format_double(matrix->values[j]); MSG_PRINT("\n");
                } else if(MESS_IS_COMPLEX(matrix)){
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e + %10e i\n", matrix->rowptr[j], i, creal(matrix->values_cpx[j]), cimag(matrix->values_cpx[j]));
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, matrix->rowptr[j]);
                    mess_print_format_double_cpx(matrix->values_cpx[j]); MSG_PRINT("\n");
                }
            }
        }
    } else if (MESS_IS_COORD(matrix)) {
        for ( i = 0; i < matrix->nnz; i++) {
            if ( MESS_IS_REAL(matrix)) {
                //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e \n", matrix->rowptr[i],matrix->colptr[i], matrix->values[i]);
                MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", matrix->rowptr[i],matrix->colptr[i]);
                mess_print_format_double(matrix->values[i]); MSG_PRINT("\n");
            } else if(MESS_IS_COMPLEX(matrix)){
                //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e + %10e i\n", matrix->rowptr[i], matrix->colptr[i], creal(matrix->values_cpx[i]), cimag(matrix->values_cpx[i]));
                MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", matrix->rowptr[i],matrix->colptr[i]);
                mess_print_format_double_cpx(matrix->values_cpx[i]); MSG_PRINT("\n");
            }
        }
    } else if (MESS_IS_DENSE(matrix)){
        for ( j = 0; j < matrix->cols; j++ ) {
            for ( i = 0; i < matrix->rows ; i++ ) {
                if ( MESS_IS_REAL(matrix)) {
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %20.15e \n", i, j, matrix->values[i+j*matrix->ld]);
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, j);
                    mess_print_format_double(matrix->values[i+j*matrix->ld]); MSG_PRINT("\n");
                } else if(MESS_IS_COMPLEX(matrix)){
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %20.15e + %20.15e i\n", i, j, creal(matrix->values_cpx[i+j*matrix->ld]), cimag(matrix->values_cpx[i+j*matrix->ld]));
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, j);
                    mess_print_format_double_cpx(matrix->values_cpx[i+j*matrix->ld]); MSG_PRINT("\n");
                }

            }
        }
    }else {
        MSG_ERROR("The storage type: %s is unknown\n", mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return 0;
}

/**
 * @brief Print the first values of a matrix to stdout.
 * @param[in] matrix  input the matrix to print
 * @return zero on success, a non zero error code otherwise
 *
 * The @ref mess_matrix_printshort function prints the first \f$ 40 \f$ values of a matrix to the standard output. \n
 * The output format is almost the same as in @matlab . \n
 * In case of sparse matrices only the non zero entries are printed to the screen.
 * We use zero-based counting for indices.
 *
 * @sa mess_matrix_print
 * @sa mess_matrix_printinfo
 * @sa mess_matrix_printdata
 * @sa mess_matrix_display
 *
 */
int mess_matrix_printshort(mess_matrix matrix) {
    MSG_FNAME(__func__);
    mess_int_t i,j;
    int brief = 1;

    mess_check_nullpointer(matrix);

    if (MESS_IS_DENSE(matrix)) {
        MSG_PRINT ("size = ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " ) LD = " MESS_PRINTF_INT " \n", matrix->rows, matrix->cols, matrix->ld);
    } else {
        MSG_PRINT ("size = ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " ) \n", matrix->rows, matrix->cols);
    }
    /*-----------------------------------------------------------------------------
     *  print a user defined matrix
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_CSR ( matrix )) {
        for ( i = 0 ; i < matrix->rows; i++){
            for (j = matrix->rowptr[i]; j<matrix->rowptr[i+1]; j++){
                if ( MESS_IS_REAL(matrix)) {
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e \n", i, matrix->colptr[j], matrix->values[j]);
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, matrix->colptr[j]);
                    mess_print_format_double(matrix->values[j]); MSG_PRINT("\n");
                } else if(MESS_IS_COMPLEX(matrix)){
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e + %10e i\n", i, matrix->colptr[j], creal(matrix->values_cpx[j]), cimag(matrix->values_cpx[j]));
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, matrix->colptr[j]);
                    mess_print_format_double_cpx(matrix->values_cpx[j]); MSG_PRINT("\n");
                }

                if ( brief == 1 && j > 40) {
                    MSG_PRINT("...\n");
                    return 0;
                }
            }
        }
    } else if ( MESS_IS_CSC ( matrix )){
        for ( i = 0 ; i < matrix->cols; i++){
            for (j = matrix->colptr[i]; j<matrix->colptr[i+1]; j++){
                if ( MESS_IS_REAL(matrix)) {
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e \n", matrix->rowptr[j], i, matrix->values[j]);
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", matrix->rowptr[j], i);
                    mess_print_format_double(matrix->values[j]); MSG_PRINT("\n");
                } else if(MESS_IS_COMPLEX(matrix)){
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e + %10e i\n", matrix->rowptr[j], i, creal(matrix->values_cpx[j]), cimag(matrix->values_cpx[j]));
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, matrix->rowptr[j]);
                    mess_print_format_double_cpx(matrix->values_cpx[j]); MSG_PRINT("\n");
                }
                if (brief == 1 && j > 40) {
                    MSG_PRINT("...\n");
                    return 0;
                }
            }
        }
    } else if (MESS_IS_COORD(matrix)) {
        for ( i = 0; i < matrix->nnz; i++) {
            if ( MESS_IS_REAL(matrix)) {
                //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e \n", matrix->rowptr[i],matrix->colptr[i], matrix->values[i]);
                MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", matrix->rowptr[i],matrix->colptr[i]);
                mess_print_format_double(matrix->values[i]); MSG_PRINT("\n");
            } else if(MESS_IS_COMPLEX(matrix)){
                //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %10e + %10e i\n", matrix->rowptr[i], matrix->colptr[i], creal(matrix->values_cpx[i]), cimag(matrix->values_cpx[i]));
                MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", matrix->rowptr[i],matrix->colptr[i]);
                mess_print_format_double_cpx(matrix->values_cpx[i]); MSG_PRINT("\n");
            }

            if (brief == 1 && i > 40) {
                MSG_PRINT("...\n");
                return 0;
            }
        }
    } else if (MESS_IS_DENSE(matrix)){
        mess_int_t k = 0;
        for ( j = 0; j < matrix->cols; j++){
            for ( i = 0; i < matrix->rows ; i++){
                k++;
                if ( MESS_IS_REAL(matrix)) {
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %20.15e \n", i, j, matrix->values[i+j*matrix->ld]);
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, j);
                    mess_print_format_double(matrix->values[i+j*matrix->ld]); MSG_PRINT("\n");
                } else if(MESS_IS_COMPLEX(matrix)){
                    //MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %20.15e + %20.15e i\n", i, j, creal(matrix->values_cpx[i+j*matrix->ld]), cimag(matrix->values_cpx[i+j*matrix->ld]));
                    MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t", i, j);
                    mess_print_format_double_cpx(matrix->values_cpx[i+j*matrix->ld]); MSG_PRINT("\n");
                }
                if ( k >= 40) break;
            }
            if (k>=40) {
                MSG_PRINT("...\n");
                break;
            }
        }

    }else {
        MSG_ERROR("The storage type: %s is unknown\n", mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return 0;
}

/**
 * @brief Print underlying data structures of a matrix to stdout.
 * @param[in] mat  input the matrix to print
 * @return zero on success, a non zero error code otherwise
 *
 * The @ref mess_matrix_printdata prints underlying data structures of a matrix to the standard output. \n
 * Depending on the storage type the printed information differs.
 * We use zero-based counting for indices.
 *
 * @sa mess_matrix_printshort
 * @sa mess_matrix_printinfo
 * @sa mess_matrix_print
 * @sa mess_matrix_display
 *
 */
int mess_matrix_printdata(mess_matrix mat)
{
    MSG_FNAME(__func__);
    mess_int_t i = 0;
    mess_check_nullpointer(mat);


    /*-----------------------------------------------------------------------------
     *  CSR
     *-----------------------------------------------------------------------------*/
    if ( MESS_IS_CSR(mat)){
        MSG_PRINT("size   =  ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " )\n", mat->rows, mat->cols);
        MSG_PRINT("nnz    = " MESS_PRINTF_INT "\n", mat->nnz);
        MSG_PRINT("rowptr = \n");
        for ( i = 0; i <= mat->rows; i++){ MSG_PRINT(" " MESS_PRINTF_INT " ", mat->rowptr[i]);}
        MSG_PRINT("\ncolptr / values = \n");
        if ( MESS_IS_COMPLEX(mat)) {
            for ( i = 0; i < mat->nnz; i++){ MSG_PRINT(" " MESS_PRINTF_INT "(%lg+%lgi) ", mat->colptr[i], creal(mat->values_cpx[i]), cimag(mat->values_cpx[i]));}
        } else if ( MESS_IS_REAL(mat)){
            for ( i = 0; i < mat->nnz; i++){ MSG_PRINT(" " MESS_PRINTF_INT "(%lg) ", mat->colptr[i], mat->values[i]);}
        }
    }

    /*-----------------------------------------------------------------------------
     *  CSC
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_CSC(mat)) {
        MSG_PRINT("size   =  ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " )\n", mat->rows, mat->cols);
        MSG_PRINT("nnz    = " MESS_PRINTF_INT "\n", mat->nnz);
        MSG_PRINT("colptr = \n");
        for ( i = 0; i <= mat->rows; i++){ MSG_PRINT(" " MESS_PRINTF_INT " ", mat->colptr[i]);}
        MSG_PRINT("\nrowptr / values = \n");
        if ( MESS_IS_COMPLEX(mat)) {
            for ( i = 0; i < mat->nnz; i++){ MSG_PRINT(" " MESS_PRINTF_INT "(%lg+%lgi) ", mat->rowptr[i], creal(mat->values_cpx[i]), cimag(mat->values_cpx[i]));}
        } else if ( MESS_IS_REAL(mat)){
            for ( i = 0; i < mat->nnz; i++){ MSG_PRINT(" " MESS_PRINTF_INT "(%lg) ", mat->rowptr[i], mat->values[i]);}
        }
    }

    /*-----------------------------------------------------------------------------
     *  DENSE
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_DENSE(mat)){
        MSG_PRINT("size   =  ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " )\n", mat->rows, mat->cols);
        MSG_PRINT("ld     =  " MESS_PRINTF_INT "\n", mat->ld);
        MSG_PRINT("values =\n");
        if ( MESS_IS_REAL(mat)){
            for ( i = 0; i < mat->ld*mat->cols; i++){
                MSG_PRINT("%lg ", mat->values[i]);
            }
        }
        else if ( MESS_IS_COMPLEX(mat)){
            for ( i = 0; i < mat->ld*mat->cols; i++){
                MSG_PRINT("%lg+%lgi ", creal(mat->values_cpx[i]), cimag(mat->values_cpx[i]));
            }
        } else {
            MSG_ERROR("data type not supported\n");
            return MESS_ERROR_DATATYPE;
        }
    }

    /*-----------------------------------------------------------------------------
     *  COORD
     *-----------------------------------------------------------------------------*/
    else if ( MESS_IS_COORD(mat)){
        MSG_PRINT("size =  ( " MESS_PRINTF_INT ", " MESS_PRINTF_INT " )\n", mat->rows, mat->cols);
        MSG_PRINT("nnz  = " MESS_PRINTF_INT "\n", mat->nnz);
        for ( i = 0; i < mat->nnz; i++) {
            if ( MESS_IS_REAL(mat)) {
                MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %lg \n", mat->rowptr[i],mat->colptr[i], mat->values[i]);
            } else if(MESS_IS_COMPLEX(mat)){
                MSG_PRINT("( " MESS_PRINTF_INT " \t, " MESS_PRINTF_INT " \t) =\t %lg + %lg i\n", mat->rowptr[i],mat->colptr[i], creal(mat->values_cpx[i]), cimag(mat->values_cpx[i]));
            }
        }
    } else {
        MSG_ERROR(" storage type not supported\n");
        return MESS_ERROR_STORAGETYPE;
    }

    return 0;
}       /* -----  end of function mess_matrix_printdata  ----- */


/**
 * @brief Print basic information of a matrix to the stdout.
 * @param[in] matrix  input the matrix to print
 * @return zero on success, a non zero error code otherwise
 *
 * The @ref mess_matrix_printinfo function prints basic information of a matrix to
 * the standard output. \n
 * These are
 * <ul>
 * <li> size,
 * <li> number of non zero elements,
 * <li> leading dimension,
 * <li> data type,
 * <li> ... .
 * </ul>
 *
 * @sa mess_matrix_printshort
 * @sa mess_matrix_print
 * @sa mess_matrix_printdata
 * @sa mess_matrix_display
 *
 */
int mess_matrix_printinfo(mess_matrix matrix){
    MSG_FNAME(__func__);
    mess_check_nullpointer(matrix);

    MSG_PRINT("Size:                "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT"\n", matrix->rows, matrix->cols);
    MSG_PRINT("Data Type:           %s\n", mess_datatype_t_str(matrix->data_type));
    MSG_PRINT("Store Type:          %s\n", mess_storage_t_str(matrix->store_type));
    MSG_PRINT("Leading Dimension:   "MESS_PRINTF_INT"\n", matrix->ld);
    MSG_PRINT("Size:                ");
    mess_print_bytes(mess_matrix_memsize(matrix));
    MSG_PRINT("\n");
    MSG_PRINT("Memory Adress:       0x%.12x\n", matrix);
    return 0;
}

