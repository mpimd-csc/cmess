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
 * @addtogroup test_misc
 * @{
 * @file tests/misc/check_freelist.c
 * @brief Check the mess_freelist_t structure.
 *
 * @author  @mbehr
 * @test
 * This function calls
 * * @ref mess_freelist_add_mess_matrix
 * * @ref mess_freelist_add_mess_vector
 * * @ref mess_freelist_add_ptr
 * * @ref mess_freelist_add_mess_equation
 * * @ref mess_freelist_add_mess_options
 * * @ref mess_freelist_add_mess_status
 * * @ref mess_freelist_print
 * * @ref mess_freelist_clear
 *. You should use this test with valgrind to check if all memory is cleared.
 *
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "mess/mess.h"
#include "../call_macro.h"



int main(int argc, char**argv){
    mess_init();
    int ret=0;
    /*-----------------------------------------------------------------------------
     *  check arguments
     *-----------------------------------------------------------------------------*/
    if(argc!=1){
        printf("Usage: %s\n",argv[0]);
        return 1;
    }

    /*-----------------------------------------------------------------------------
     *  create data for freelist
     *-----------------------------------------------------------------------------*/
    double *ptr1=NULL, *ptr2=NULL, *ptr3=NULL;
    mess_vector v1=NULL, v2=NULL, v3=NULL;
    mess_matrix A1=NULL, A2=NULL, A3=NULL;
    mess_equation eqn1=NULL, eqn2=NULL, eqn3=NULL;
    mess_options opt1=NULL, opt2=NULL, opt3=NULL;
    mess_status stat1=NULL, stat2=NULL, stat3=NULL;


    /*-----------------------------------------------------------------------------
     *   create freelist and add instances to list and print
     *-----------------------------------------------------------------------------*/
    mess_freelist mem;
    CALL(mess_freelist_init(&mem));


    /*-----------------------------------------------------------------------------
     *  add pointers
     *-----------------------------------------------------------------------------*/
    ptr1 = (double*)malloc(10*sizeof(double));
    ptr3 = (double*)malloc(10*sizeof(double));
    CALL(mess_freelist_add_ptr(mem,ptr1));
    CALL(mess_freelist_add_ptr(mem,ptr2));
    CALL(mess_freelist_add_ptr(mem,ptr3));
    CALL(mess_freelist_print(mem));
    printf("\n###############################################################\n\n");


    /*-----------------------------------------------------------------------------
     *  add vectors
     *-----------------------------------------------------------------------------*/
    MESS_INIT_VECTORS(&v1,&v3);
    mess_vector_alloc(v1,1,MESS_REAL);
    mess_vector_alloc(v3,3,MESS_COMPLEX);
    CALL(mess_freelist_add_mess_vector(mem,v1));
    CALL(mess_freelist_add_mess_vector(mem,v2));
    CALL(mess_freelist_add_mess_vector(mem,v3));
    CALL(mess_freelist_print(mem));
    printf("\n###############################################################\n\n");


    /*-----------------------------------------------------------------------------
     *  add matrices
     *-----------------------------------------------------------------------------*/
    MESS_INIT_MATRICES(&A1,&A3);
    CALL(mess_freelist_add_mess_matrix(mem,A1));
    CALL(mess_freelist_add_mess_matrix(mem,A2));
    CALL(mess_freelist_add_mess_matrix(mem,A3));
    CALL(mess_freelist_print(mem));
    printf("\n###############################################################\n\n");


    /*-----------------------------------------------------------------------------
     *  add equations
     *-----------------------------------------------------------------------------*/
    MESS_INIT_EQUATIONS(&eqn1,&eqn3);
    CALL(mess_freelist_add_mess_equation(mem,eqn1));
    CALL(mess_freelist_add_mess_equation(mem,eqn2));
    CALL(mess_freelist_add_mess_equation(mem,eqn3));
    CALL(mess_freelist_print(mem));
    printf("\n###############################################################\n\n");


    /*-----------------------------------------------------------------------------
     *  add options
     *-----------------------------------------------------------------------------*/
    mess_options_init(&opt1);
    mess_options_init(&opt3);
    CALL(mess_freelist_add_mess_options(mem,opt1));
    CALL(mess_freelist_add_mess_options(mem,opt2));
    CALL(mess_freelist_add_mess_options(mem,opt3));
    CALL(mess_freelist_print(mem));
    printf("\n###############################################################\n\n");


    /*-----------------------------------------------------------------------------
     *  add status
     *-----------------------------------------------------------------------------*/
    mess_status_init(&stat1);
    mess_status_init(&stat2);
    CALL(mess_freelist_add_mess_status(mem,stat1));
    CALL(mess_freelist_add_mess_status(mem,stat2));
    CALL(mess_freelist_add_mess_status(mem,stat3));
    CALL(mess_freelist_print(mem));
    printf("\n###############################################################\n\n");


    /*-----------------------------------------------------------------------------
     *  clear tht freelist
     *-----------------------------------------------------------------------------*/
    mess_freelist_clear(&mem);

    return 0;
}


