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
 * @file lib/io/vector_print.c
 * @brief Print @ref mess_vector instance and underlying information to the standard output.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <math.h>
#include <complex.h>


/**
 * @brief Print out a @ref mess_vector structure.
 * @param[in] vect      input   vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_print function prints all components of a
 * @ref mess_vector structure on stdout. \n
 * If the given vector points to @c NULL a @ref MESS_ERROR_NULLPOINTER
 * is returned. \n
 * If the data type is unknown @ref MESS_ERROR_DATATYPE is returned.
 * We use zero based indexing.
 *
 */
int mess_vector_print(mess_vector vect) {

    MSG_FNAME(__func__);
    mess_int_t i, dim;

    mess_check_nullpointer(vect);

    dim = vect->dim;

    if(dim<=0){
        MSG_PRINT("Vector is empty.\n");
    }

    mess_vector_printinfo(vect);
    if (MESS_IS_COMPLEX(vect)) {
        for (i = 0; i < dim; i++){
            MSG_PRINT(" [ "MESS_PRINTF_INT"\t ] =\t", i);
            mess_print_format_double_cpx(vect->values_cpx[i]);
            MSG_PRINT("\n");
        }
    } else if (MESS_IS_REAL(vect)){
        for (i = 0; i < dim; i++) {
            MSG_PRINT(" [ "MESS_PRINTF_INT"\t ] =\t", i);
            mess_print_format_double(vect->values[i]);
            MSG_PRINT("\n");
        }
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    return 0;
} /* -----  end of function mess_vector_print  ----- */




/**
 * @brief Print out a @ref mess_vector structure in a short form.
 * @param[in] vect input    vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_printshort function is the same like @ref mess_vector_print
 * except of that only the first \f$ 20 \f$ components are printed.
 * We use zero based indexing.
 *
 */
int mess_vector_printshort(mess_vector vect) {
    MSG_FNAME(__func__);
    mess_int_t i, dim;

    mess_check_nullpointer(vect);

    dim = (vect->dim < 20) ? vect->dim : 20;

    mess_vector_printinfo(vect);
    if (MESS_IS_COMPLEX(vect)) {
        for (i = 0; i < dim; i++){
            MSG_PRINT(" [ "MESS_PRINTF_INT"\t ] =\t", i);
            mess_print_format_double_cpx(vect->values_cpx[i]);
            MSG_PRINT("\n");
        }
    } else if (MESS_IS_REAL(vect)) {
        for (i = 0; i < dim; i++) {
            MSG_PRINT(" [ "MESS_PRINTF_INT"\t ] =\t", i);
            mess_print_format_double(vect->values[i]);
            MSG_PRINT("\n");
        }
    } else {
        MSG_ERROR("unknown/unsupported data type\n");
        return MESS_ERROR_DATATYPE;
    }
    MSG_PRINT(" ...\n");
    return 0;
} /* -----  end of function mess_vector_printshort  ----- */


/**
 * @brief Print dimension and datatype of a @ref mess_vector structure on stdout.
 * @param[in] v input vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_vector_printinfo function prints dimension and datatype of
 * a @ref mess_vector structure on stdout.
 *
 */
int mess_vector_printinfo(mess_vector v) {
    MSG_FNAME(__func__);

    mess_check_nullpointer(v);
    MSG_PRINT("dimension:  " MESS_PRINTF_INT "\n", v->dim);
    MSG_PRINT("datatype:   %s\n", mess_datatype_t_str(v->data_type));

    return 0;
} /* -----  end of function mess_vector_printinfo  ----- */


