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
 * @addtogroup test_matrix
 * @{
 * @file tests/matrix/check_trace.c
 * @brief Check the trace computation.
 * @author  @mbehr
 * @test
 * This function checks the mess_matrix_trace and @ref mess_matrix_tracec functions defined in trace.c that means it checks
 * if the returned trace of a given matrix is correct.
 *
 * @}
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/mess.h"

#include "../call_macro.h"




double trreal=0;            //for CHECKTRACE
mess_double_cpx_t trcpx=0;     //for CHECKTRACE
#define CHECKTRACE(A, CHECK, EPS , ERR ){                       \
    trreal=0;                                               \
    trcpx =0;                                               \
    if(MESS_IS_REAL((A))){                                  \
        mess_matrix_trace((A),&trreal);                     \
        if(fabs(trreal-((double)CHECK)) > 10* (EPS)){++(ERR); \
            fprintf(stderr,"err: %20.16e\n", fabs(trreal-((double)CHECK)));\
        }       \
    }else if(MESS_IS_COMPLEX((A))){                         \
        mess_matrix_tracec((A), &trcpx);                    \
        if(cabs(trcpx-((mess_double_cpx_t)CHECK)) > 10* (EPS)){++(ERR);};      \
    }                                                       \
}                                                               \



int main ( int argc, char ** argv) {
    mess_init();

    double eps=mess_eps();
    int ret=0, err=0;

    mess_int_t rows,cols; // size of matrix
    rows=5;
    cols=5;

    //prepare dense matrices
    mess_matrix denser, densec;
    mess_matrix_init(&denser);
    mess_matrix_init(&densec);
    mess_matrix_alloc(denser,rows,cols,rows*cols,MESS_DENSE,MESS_REAL);
    mess_matrix_alloc(densec,rows,cols,rows*cols,MESS_DENSE,MESS_COMPLEX);

    mess_int_t i,j;
    double checktracereal=0;
    mess_double_cpx_t checktracecpx=0;
    for(i=0;i<rows;++i){
        for(j=0;j<cols;++j){
            denser->values[i+j*denser->ld]=(double)rand()/((double) (RAND_MAX));
            densec->values_cpx[i+j*denser->ld]=(double)rand()/((double) (RAND_MAX))+(double)rand()/((double) (RAND_MAX))*I;
            if(i==j){
                checktracereal+=denser->values[i+j*denser->ld];
                checktracecpx+=densec->values_cpx[i+j*denser->ld];
            }

        }
    }

    //prepare csr matrices
    mess_matrix csrr,csrc;
    mess_matrix_init(&csrr);
    mess_matrix_init(&csrc);
    CALL(mess_matrix_convert(denser,csrr,MESS_CSR));
    CALL(mess_matrix_convert(densec,csrc,MESS_CSR));

    //prepare csr matrices
    mess_matrix cscr,cscc;
    mess_matrix_init(&cscr);
    mess_matrix_init(&cscc);
    CALL(mess_matrix_convert(denser,cscr,MESS_CSC));
    CALL(mess_matrix_convert(densec,cscc,MESS_CSC));


    //Check DENSE Matrices
    CHECKTRACE(densec,checktracecpx,eps,err);
    CHECKTRACE(denser,checktracereal,eps,err);

    //Check CSC Matrices
    CHECKTRACE(cscc,checktracecpx,eps,err);
    CHECKTRACE(cscr,checktracereal,eps,err);

    //Check CSR Matrices
    CHECKTRACE(csrc,checktracecpx,eps,err);
    CHECKTRACE(csrr,checktracereal,eps,err);


    mess_matrix_clear(&denser);
    mess_matrix_clear(&densec);
    mess_matrix_clear(&cscr);
    mess_matrix_clear(&cscc);
    mess_matrix_clear(&csrr);
    mess_matrix_clear(&csrc);

    return (err>0)?(1):(0);
}

