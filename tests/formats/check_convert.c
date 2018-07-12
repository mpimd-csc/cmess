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
 * @addtogroup test_format
 * @{
 * @file tests/formats/check_convert.c
 * @brief  Check the conversion subroutine.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_matrix_convert function defined in convert.c that means it checks if a matrix can be
 * converted to another storage scheme. \n
 * Tested conversions are:
 * \li CSR \f$ \to \f$ CSR
 * \li CSR \f$ \to \f$ CSC
 * \li CSR \f$ \to \f$ DENSE
 * \li CSR \f$ \to \f$ COORD
 * \li CSC \f$ \to \f$ CSC
 * \li CSC \f$ \to \f$ CSR
 * \li CSC \f$ \to \f$ DENSE
 * \li CSC \f$ \to \f$ COORD
 * \li COORD \f$ \to \f$ COORD
 * \li COORD \f$ \to \f$ CSR
 * \li COORD \f$ \to \f$ CSC
 * \li COORD \f$ \to \f$ DENSE
 * \li DENSE \f$ \to \f$ DENSE
 * \li DENSE \f$ \to \f$ CSR
 * \li DENSE \f$ \to \f$ CSC
 * \li DENSE \f$ \to \f$ COORD
 *
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mess/mess.h"
#include "../call_macro.h"

/**
 * @brief Check the conversion function.
 * @param[in] A     input matrix
 * @param[in] from  input type
 * @param[in] to     input output type
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref check function converts a matrix \f$ A \f$ to the desired output type \f$ to \f$ using @ref mess_matrix_convert and
 * checks if the computed matrix is in the desired storage format.
 *
 */

int check(mess_matrix A, unsigned short from, unsigned short to ){
    mess_matrix inter1;
    mess_matrix inter2;
    int ret ;

    printf("Check from %s via %s and %s to %s ...", mess_storage_t_str(A->store_type), mess_storage_t_str(from), mess_storage_t_str(to),
            mess_storage_t_str(A->store_type));
    CALL(mess_matrix_init(&inter1));
    CALL(mess_matrix_init(&inter2));
    CALL(mess_matrix_convert(A,inter1,from));
    // printf("orig:\n");
    // mess_matrix_printshort(A);
    // printf("inter1:\n");
    // mess_matrix_printshort(inter1);
    CALL(mess_matrix_convert(inter1, inter2, to));
    // printf("inter2:\n");
    // mess_matrix_printshort(inter2);
    CALL(mess_matrix_convert(inter2, inter1, A->store_type));
    // printf("inter1_final\n");
    // mess_matrix_printshort(inter1);
    CALL(mess_matrix_sort(inter1));

    if ( mess_matrix_equal_verbose(A, inter1)){
        printf("passed.\n");
        CALL(mess_matrix_clear(&inter1));
        CALL(mess_matrix_clear(&inter2));
        return 0;
    } else {
        printf("failed.\n");
        CALL(mess_matrix_clear(&inter1));
        CALL(mess_matrix_clear(&inter2));
        return 1;
    }
    CALL(mess_matrix_clear(&inter1));
    CALL(mess_matrix_clear(&inter2));
    return 1;
}

int main ( int argc, char **argv){
    mess_init();
    mess_matrix A;
    int ret ;
    if ( argc !=2 ) {
        fprintf(stderr,"Need a matrix market file as argument.\n");
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  read data
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_init(&A));
    CALL(mess_matrix_read(argv[1], A));
    CALL(mess_matrix_sort(A));


    /*-----------------------------------------------------------------------------
     *  Tests
     *-----------------------------------------------------------------------------*/
    CALL(check(A,MESS_DENSE, MESS_DENSE));
    CALL(check(A,MESS_DENSE, MESS_CSR));
    CALL(check(A,MESS_DENSE, MESS_CSC));
    CALL(check(A,MESS_DENSE, MESS_COORD));

    CALL(check(A,MESS_CSR, MESS_DENSE));
    CALL(check(A,MESS_CSR, MESS_CSR));
    CALL(check(A,MESS_CSR, MESS_CSC));
    CALL(check(A,MESS_CSR, MESS_COORD));

    CALL(check(A,MESS_CSC, MESS_DENSE));
    CALL(check(A,MESS_CSC, MESS_CSR));
    CALL(check(A,MESS_CSC, MESS_CSC));
    CALL(check(A,MESS_CSC, MESS_COORD));

    CALL(check(A,MESS_COORD, MESS_DENSE));
    CALL(check(A,MESS_COORD, MESS_CSR));
    CALL(check(A,MESS_COORD, MESS_CSC));
    CALL(check(A,MESS_COORD, MESS_COORD));



    /*-----------------------------------------------------------------------------
     *  complex
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_tocomplex(A));

    /*-----------------------------------------------------------------------------
     *  Tests, Part2
     *-----------------------------------------------------------------------------*/
    CALL(check(A,MESS_DENSE, MESS_DENSE));
    CALL(check(A,MESS_DENSE, MESS_CSR));
    CALL(check(A,MESS_DENSE, MESS_CSC));
    CALL(check(A,MESS_DENSE, MESS_COORD));

    CALL(check(A,MESS_CSR, MESS_DENSE));
    CALL(check(A,MESS_CSR, MESS_CSR));
    CALL(check(A,MESS_CSR, MESS_CSC));
    CALL(check(A,MESS_CSR, MESS_COORD));

    CALL(check(A,MESS_CSC, MESS_DENSE));
    CALL(check(A,MESS_CSC, MESS_CSR));
    CALL(check(A,MESS_CSC, MESS_CSC));
    CALL(check(A,MESS_CSC, MESS_COORD));

    CALL(check(A,MESS_COORD, MESS_DENSE));
    CALL(check(A,MESS_COORD, MESS_CSR));
    CALL(check(A,MESS_COORD, MESS_CSC));
    CALL(check(A,MESS_COORD, MESS_COORD));




    /*-----------------------------------------------------------------------------
     *  clean up
     *-----------------------------------------------------------------------------*/
    CALL(mess_matrix_clear(&A));
    return 0;
}



