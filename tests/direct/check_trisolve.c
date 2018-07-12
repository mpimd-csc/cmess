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
 * @addtogroup test_direct
 * @{
 * @file tests/direct/check_trisolve.c
 * @brief Check the functions in trisolve.c.
 * @author  @mbehr
 * @test
 * This function checks all functions defined in trisolve.c that means it checks if the solution of
 * \f[LU X= B ,\f]
 * \f[op(L) x = b ,\f]
 * \f[op(L) X = B ,\f]
 * \f[op(U) x = b ,\f]
 * \f[op(U) X = B \f]
 * is computed correctly for given matrices \f$ L, U \f$ and \f$ B \f$ and given vector \f$ b \f$.
 * Operations \f$ op (.)\f$ on matrices \f$ L \f$ and \f$ U \f$ can be
 * <ul>
 * <li> \ref MESS_OP_NONE (\f$ op(C)= C \f$),
 * <li> \ref MESS_OP_TRANSPOSE (\f$ op(C) = C^T \f$),
 * <li> \ref MESS_OP_HERMITIAN (\f$ op(C) = C^H \f$).
 * </ul>
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



#define CHECK_TRISOLVE(SOLVER,L,SOL,B,STORETYPE,L2,X,ERR,EPS,DIFF){         \
    CALL(mess_matrix_convert((L),(L2),(STORETYPE)));                        \
    CALL(mess_vector_copy((B),(X)));                                        \
    CALL(SOLVER((L2),(X)));                                                 \
    CALL(mess_vector_diffnorm((X),(SOL),&(DIFF)));                          \
        if(DIFF>EPS){                                                   \
            printf("Failed:\n");                                            \
            printf("Solver: %s\n",#SOLVER);                                 \
            printf("Storetype: %s\n", #STORETYPE);                          \
            printf("L2: %s\n", #L2);                                        \
            CALL(mess_matrix_printinfo(L2));                                \
            printf("RHS: %s\n",#B);                                         \
            CALL(mess_vector_printinfo(B));                                 \
            printf("Diff: %e\n",DIFF);                                      \
            ++ERR;                                                          \
        }                                                                   \
    /*printf("diff=%e\n",diff);*/                                           \
    /*mess_matrix_printinfo(L2);*/                                          \
    /*mess_vector_printinfo(x);*/                                           \
}


#define CHECK_ALL_STORAGETYPE(MESS_DIRECT,SOLVER,OP,L,SOL,B,L2,X,ERR,EPS,DIFF)  \
        CALL(mess_direct_init(&MESS_DIRECT));                                   \
        CALL(mess_direct_create_lapack_lu(L,MESS_DIRECT));                              \
        CALL(mess_direct_solve(OP,MESS_DIRECT,B,SOL));                          \
        CHECK_TRISOLVE(SOLVER,L,SOL,B,MESS_DENSE,L2,X,ERR,EPS,DIFF)             \
        CHECK_TRISOLVE(SOLVER,L,SOL,B,MESS_CSC,L2,X,ERR,EPS,DIFF)               \
        CHECK_TRISOLVE(SOLVER,L,SOL,B,MESS_CSR,L2,X,ERR,EPS,DIFF)               \
        CALL(mess_direct_clear(&MESS_DIRECT));                                  \


/*macros for checking *m trisolve functions where X and B are matrices*/
#define CHECK_TRISOLVEM(SOLVER,L,SOL,B,STORETYPE,L2,X,ERR,EPS,DIFF){        \
    CALL(mess_matrix_convert((L),(L2),(STORETYPE)));                        \
    CALL(mess_matrix_copy((B),(X)));                                        \
    CALL(SOLVER((L2),(X)));                                                 \
    CALL(mess_matrix_diffnorm((X),(SOL),&(DIFF)));                              \
        if(DIFF>EPS){                                                   \
            printf("Failed:\n");                                            \
            printf("Solver: %s\n",#SOLVER);                                 \
            printf("Storetype: %s\n", #STORETYPE);                          \
            printf("L2: %s\n", #L2);                                        \
            CALL(mess_matrix_printinfo(L2));                                \
            printf("RHS: %s\n",#B);                                         \
            CALL(mess_matrix_printinfo(B));                                 \
            printf("Diff: %e\n",DIFF);                                      \
            ++ERR;                                                          \
        }                                                                   \
    /*printf("diff=%e\n",diff);*/                                           \
    /*mess_matrix_printinfo(L2);*/                                          \
    /*mess_vector_printinfo(x);*/                                           \
}


#define CHECK_ALL_STORAGETYPEM(MESS_DIRECT,SOLVER,OP,L,SOL,B,L2,X,ERR,EPS,DIFF) \
        CALL(mess_direct_init(&MESS_DIRECT));                                   \
        CALL(mess_direct_create_lapack_lu(L,MESS_DIRECT));                              \
        CALL(mess_direct_solvem(OP,MESS_DIRECT,B,SOL));                         \
        CHECK_TRISOLVEM(SOLVER,L,SOL,B,MESS_DENSE,L2,X,ERR,EPS,DIFF)            \
        CHECK_TRISOLVEM(SOLVER,L,SOL,B,MESS_CSC,L2,X,ERR,EPS,DIFF)              \
        CHECK_TRISOLVEM(SOLVER,L,SOL,B,MESS_CSR,L2,X,ERR,EPS,DIFF)              \
        CALL(mess_direct_clear(&MESS_DIRECT));                                  \




int main ( int argc, char ** argv) {
    int rows=10, cols=10, i, j, ret=0, err=0;
    double eps = 1e-8, p=0.5, diff=.0;

    mess_vector x, sol, b;
    mess_matrix  L,U, LU2,Brows,Bcols,Xmat,Solmat;
    mess_direct solver ;

    /*-----------------------------------------------------------------------------
     *  Init/Load Vectors
     *-----------------------------------------------------------------------------*/

    MESS_INIT_VECTORS(&x,&sol,&b);
    CALL(mess_vector_alloc(x,cols,MESS_REAL));
    CALL(mess_vector_alloc(sol,cols,MESS_REAL));
    CALL(mess_vector_alloc(b,rows,MESS_REAL ));
    CALL(mess_vector_rand(b));

    /*-----------------------------------------------------------------------------
     * Init/Load matrices
     *-----------------------------------------------------------------------------*/

    CALL(mess_matrix_init(&L));
    CALL(mess_matrix_init(&LU2));
    CALL(mess_matrix_init(&U));
    CALL(mess_matrix_init(&Brows));
    CALL(mess_matrix_init(&Bcols));
    CALL(mess_matrix_init(&Xmat));
    CALL(mess_matrix_init(&Solmat));


    //load real matrices
    CALL(mess_matrix_rand(L,rows,cols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(U,rows,cols,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(Brows,rows,rows,MESS_DENSE,MESS_REAL,p));
    CALL(mess_matrix_rand(Bcols,cols,cols,MESS_DENSE,MESS_REAL,p));
    //delete upper triangular and add 1 to main diagonal to ensure system is solvable
    for(i=0;i<rows;++i){
        for(j=0;j<cols;++j){
            if (i==j){
                L->values[j*U->ld+i]+=1;
            }else if (i<j){
                L->values[j*U->ld+i]=0;
            }
        }
    }

    //delete lower triangular and add 1 to main diagonal to ensure system is solvable
    for(i=0;i<rows;++i){
            for(j=0;j<cols;++j){
                if (i==j){
                    U->values[j*U->ld+i]+=1;
                }else if (i>j){
                    U->values[j*U->ld+i]=0;
                }
            }
        }

    /*-----------------------------------------------------------------------------
     * Init direct solvers
     *-----------------------------------------------------------------------------*/

    //mess_direct_init(&solver);

    /*-----------------------------------------------------------------------------
     *  Test different cases
     *-----------------------------------------------------------------------------*/

    //real L/U and real rhs
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lsolve,MESS_OP_NONE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_usolve,MESS_OP_NONE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_ltsolve,MESS_OP_TRANSPOSE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_utsolve,MESS_OP_TRANSPOSE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lhsolve,MESS_OP_HERMITIAN,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_uhsolve,MESS_OP_HERMITIAN,U,sol,b,LU2,x,err,eps,diff)

    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lsolvem,MESS_OP_NONE,L,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_usolvem,MESS_OP_NONE,U,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_ltsolvem,MESS_OP_TRANSPOSE,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_utsolvem,MESS_OP_TRANSPOSE,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lhsolvem,MESS_OP_HERMITIAN,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_uhsolvem,MESS_OP_HERMITIAN,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)



    //real L/U and complex rhs
    mess_vector_tocomplex(b);
    mess_matrix_tocomplex(Brows);
    mess_matrix_tocomplex(Bcols);
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lsolve,MESS_OP_NONE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_usolve,MESS_OP_NONE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_ltsolve,MESS_OP_TRANSPOSE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_utsolve,MESS_OP_TRANSPOSE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lhsolve,MESS_OP_HERMITIAN,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_uhsolve,MESS_OP_HERMITIAN,U,sol,b,LU2,x,err,eps,diff)

    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lsolvem,MESS_OP_NONE,L,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_usolvem,MESS_OP_NONE,U,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_ltsolvem,MESS_OP_TRANSPOSE,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_utsolvem,MESS_OP_TRANSPOSE,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lhsolvem,MESS_OP_HERMITIAN,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_uhsolvem,MESS_OP_HERMITIAN,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)


    //complex L and complex rhs
    mess_matrix_tocomplex(L);
    mess_matrix_tocomplex(U);
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lsolve,MESS_OP_NONE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_usolve,MESS_OP_NONE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_ltsolve,MESS_OP_TRANSPOSE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_utsolve,MESS_OP_TRANSPOSE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lhsolve,MESS_OP_HERMITIAN,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_uhsolve,MESS_OP_HERMITIAN,U,sol,b,LU2,x,err,eps,diff)

    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lsolvem,MESS_OP_NONE,L,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_usolvem,MESS_OP_NONE,U,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_ltsolvem,MESS_OP_TRANSPOSE,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_utsolvem,MESS_OP_TRANSPOSE,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lhsolvem,MESS_OP_HERMITIAN,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_uhsolvem,MESS_OP_HERMITIAN,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)

    //complex L and complex rhs
    CALL(mess_matrix_scalec(1+I,L));
    CALL(mess_vector_scalec(2-I,b));
    CALL(mess_matrix_scalec(2-I,Brows));
    CALL(mess_matrix_scalec(2-I,Bcols));
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lsolve,MESS_OP_NONE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_usolve,MESS_OP_NONE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_ltsolve,MESS_OP_TRANSPOSE,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_utsolve,MESS_OP_TRANSPOSE,U,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_lhsolve,MESS_OP_HERMITIAN,L,sol,b,LU2,x,err,eps,diff)
    CHECK_ALL_STORAGETYPE(solver,mess_solver_uhsolve,MESS_OP_HERMITIAN,U,sol,b,LU2,x,err,eps,diff)

    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lsolvem,MESS_OP_NONE,L,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_usolvem,MESS_OP_NONE,U,Solmat,Brows,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_ltsolvem,MESS_OP_TRANSPOSE,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_utsolvem,MESS_OP_TRANSPOSE,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_lhsolvem,MESS_OP_HERMITIAN,L,Solmat,Bcols,LU2,Xmat,err,eps,diff)
    CHECK_ALL_STORAGETYPEM(solver,mess_solver_uhsolvem,MESS_OP_HERMITIAN,U,Solmat,Bcols,LU2,Xmat,err,eps,diff)

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/

    CALL(mess_vector_clear(&x));
    CALL(mess_vector_clear(&sol));
    CALL(mess_vector_clear(&b));
    CALL(mess_matrix_clear(&L));
    CALL(mess_matrix_clear(&LU2));
    CALL(mess_matrix_clear(&Brows));
    CALL(mess_matrix_clear(&Bcols));
    CALL(mess_matrix_clear(&U));
    CALL(mess_matrix_clear(&Xmat));
    CALL(mess_matrix_clear(&Solmat));
    CALL(mess_direct_clear(&solver));

    return (err>0)?(1):(0);
}
