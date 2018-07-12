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
 * @file lib/matrix/mult.c
 * @brief Matrix-matrix multipication.
 * @author @koehlerm
 *
 * This file implements \f$ C=op(A)*op(B) \f$ where  \f$ op(X) \f$  can be
 *   <center>
 *  |Operation Type              |   \f$op\f$                 |
 *  |:--------------------------:|:--------------------------:|
 *  |@ref MESS_OP_NONE           |   \f$ op(A)=A\f$           |
 *  |@ref MESS_OP_TRANSPOSE      |   \f$ op(A)=A^T\f$         |
 *  |@ref MESS_OP_HERMITIAN      |   \f$ op(A)=A^H\f$         |
 * </center>
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>

#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "blas_defs.h"


#ifdef MESS_HAVE_CSPARSE
#include <cs.h>
#include "mess/interface_csparse.h"
#endif

/*-----------------------------------------------------------------------------
 *  Dense * Sparse
 *-----------------------------------------------------------------------------*/
void F77_GLOBAL(ddgemm_dense_sparsenn,DDGEMM_DENSE_SPARSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsenn,DZGEMM_DENSE_SPARSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsenn,ZDGEMM_DENSE_SPARSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsenn,ZZGEMM_DENSE_SPARSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparsent,DDGEMM_DENSE_SPARSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsent,DZGEMM_DENSE_SPARSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsent,ZDGEMM_DENSE_SPARSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsent,ZZGEMM_DENSE_SPARSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparsenh,DDGEMM_DENSE_SPARSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsenh,DZGEMM_DENSE_SPARSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsenh,ZDGEMM_DENSE_SPARSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsenh,ZZGEMM_DENSE_SPARSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparsetn,DDGEMM_DENSE_SPARSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsetn,DZGEMM_DENSE_SPARSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsetn,ZDGEMM_DENSE_SPARSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsetn,ZZGEMM_DENSE_SPARSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparsett,DDGEMM_DENSE_SPARSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsett,DZGEMM_DENSE_SPARSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsett,ZDGEMM_DENSE_SPARSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsett,ZZGEMM_DENSE_SPARSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparseth,DDGEMM_DENSE_SPARSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparseth,DZGEMM_DENSE_SPARSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparseth,ZDGEMM_DENSE_SPARSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparseth,ZZGEMM_DENSE_SPARSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparsehn,DDGEMM_DENSE_SPARSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsehn,DZGEMM_DENSE_SPARSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsehn,ZDGEMM_DENSE_SPARSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsehn,ZZGEMM_DENSE_SPARSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparseht,DDGEMM_DENSE_SPARSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparseht,DZGEMM_DENSE_SPARSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparseht,ZDGEMM_DENSE_SPARSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparseht,ZZGEMM_DENSE_SPARSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_dense_sparsehh,DDGEMM_DENSE_SPARSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_dense_sparsehh,DZGEMM_DENSE_SPARSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, double *A        , mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_dense_sparsehh,ZDGEMM_DENSE_SPARSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, double *Bvalues        , mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_dense_sparsehh,ZZGEMM_DENSE_SPARSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_double_cpx_t *A, mess_int_t *LDA, mess_int_t *Browptr, mess_int_t *Bcolptr, mess_double_cpx_t *Bvalues, mess_double_cpx_t *C, mess_int_t *LDC);

/*-----------------------------------------------------------------------------
 *  Sparse * Dense
 *-----------------------------------------------------------------------------*/
void F77_GLOBAL(ddgemm_sparse_densenn,DDGEMM_SPARSE_DENSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densenn,DZGEMM_SPARSE_DENSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densenn,ZDGEMM_SPARSE_DENSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densenn,ZZGEMM_SPARSE_DENSENN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_densent,DDGEMM_SPARSE_DENSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densent,DZGEMM_SPARSE_DENSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densent,ZDGEMM_SPARSE_DENSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densent,ZZGEMM_SPARSE_DENSENT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_densenh,DDGEMM_SPARSE_DENSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densenh,DZGEMM_SPARSE_DENSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densenh,ZDGEMM_SPARSE_DENSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densenh,ZZGEMM_SPARSE_DENSENH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_densetn,DDGEMM_SPARSE_DENSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densetn,DZGEMM_SPARSE_DENSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densetn,ZDGEMM_SPARSE_DENSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densetn,ZZGEMM_SPARSE_DENSETN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_densett,DDGEMM_SPARSE_DENSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densett,DZGEMM_SPARSE_DENSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densett,ZDGEMM_SPARSE_DENSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densett,ZZGEMM_SPARSE_DENSETT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_denseth,DDGEMM_SPARSE_DENSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_denseth,DZGEMM_SPARSE_DENSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_denseth,ZDGEMM_SPARSE_DENSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_denseth,ZZGEMM_SPARSE_DENSETH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_densehn,DDGEMM_SPARSE_DENSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densehn,DZGEMM_SPARSE_DENSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densehn,ZDGEMM_SPARSE_DENSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densehn,ZZGEMM_SPARSE_DENSEHN)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_denseht,DDGEMM_SPARSE_DENSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_denseht,DZGEMM_SPARSE_DENSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_denseht,ZDGEMM_SPARSE_DENSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_denseht,ZZGEMM_SPARSE_DENSEHT)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

void F77_GLOBAL(ddgemm_sparse_densehh,DDGEMM_SPARSE_DENSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , double *B        , mess_int_t *LDB, double *C        , mess_int_t *LDC);
void F77_GLOBAL(dzgemm_sparse_densehh,DZGEMM_SPARSE_DENSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, double *Avalues        , mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zdgemm_sparse_densehh,ZDGEMM_SPARSE_DENSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, double *B        , mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);
void F77_GLOBAL(zzgemm_sparse_densehh,ZZGEMM_SPARSE_DENSEHH)(mess_int_t *M, mess_int_t *N, mess_int_t *K, mess_int_t *Arowptr, mess_int_t *Acolptr, mess_double_cpx_t *Avalues, mess_double_cpx_t *B, mess_int_t *LDB, mess_double_cpx_t *C, mess_int_t *LDC);

/**
 * @brief Computes the matrix-matrix product \f$ C=op(A) op(B) \f$.
 * @param[in] opA    input operation on \f$ A \f$
 * @param[in] A      input left hand side matrix \f$A\f$
 * @param[in] opB    input operation on  \f$ B \f$
 * @param[in] B      input right hand side matrix \f$B\f$
 * @param[out] C    \f$ C=op(A) op(B) \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_multiply function computes
 * \f[C=op(A) op(B) , \f]
 * where \f$ op(X) \f$ can be
 * <center>
 *  |   Operation Type              |   \f$op\f$                    |
 *  |:-----------------------------:|:-----------------------------:|
 *  | @ref MESS_OP_NONE             |   \f$ op(X)=X\f$              |
 *  | @ref MESS_OP_TRANSPOSE        |   \f$ op(X)=X^T\f$            |
 *  | @ref MESS_OP_HERMITIAN        |   \f$ op(X)=X^H\f$            |
 * </center>
 * The function supports Dense and Sparse (@ref MESS_CSR, @ref MESS_CSC), real and complex  matrices with all operations types.
 * All other matrix types are not supported.
 */
int mess_matrix_multiply(mess_operation_t opA, mess_matrix A, mess_operation_t opB, mess_matrix B, mess_matrix C)
{
    MSG_FNAME(__func__);
    mess_int_t m, n, k1, k2, k;
    mess_int_t lda, ldb, ldc;
    int ret = 0;
    char opAi[2]={0,0};
    char opBi[2]={0,0};


    /*-----------------------------------------------------------------------------
     *  check inputs
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(B);
    mess_check_nullpointer(C);
    mess_check_real_or_complex(A);
    mess_check_real_or_complex(B);

    MESS_MATRIX_RESET(C);


    /*-----------------------------------------------------------------------------
     *  prepare
     *-----------------------------------------------------------------------------*/
    k1 = k2 = 0;
    switch(opA){
        case MESS_OP_NONE:
            m=A->rows;
            k1=A->cols;
            opAi[0]='N';
            break;
        case MESS_OP_HERMITIAN:
            m=A->cols;
            k1=A->rows;
            opAi[0]='C';
            break;
        case MESS_OP_TRANSPOSE:
            m=A->cols;
            k1=A->rows;
            opAi[0]='T';
            break;
    }
    switch(opB){
        case MESS_OP_NONE:
            n=B->cols;
            k2=B->rows;
            opBi[0]='N';
            break;
        case MESS_OP_HERMITIAN:
            n=B->rows;
            k2=B->cols;
            ldb = n;
            opBi[0]='C';
            break;
        case MESS_OP_TRANSPOSE:
            n=B->rows;
            k2=B->cols;
            opBi[0]='T';
            break;
    }

    if ( k1 != k2) {
        MSG_ERROR("dimension mismatch k1=" MESS_PRINTF_INT ",\t k2=" MESS_PRINTF_INT "\n",k1,k2);
        return (MESS_ERROR_DIMENSION);
    }

    k = k1;
    lda = A->ld;
    ldb = B->ld;


    /* if ( mess_error_level > 2 ) {
       MSG_INFO("GEMM: A(%c,%s) B(%c,%s)\n", opA[0], mess_storage_t_str(A->store_type), opB[0], mess_storage_t_str(B->store_type));

       }
       MSG_WARN("opA = %c opB = %c A->rows = " MESS_PRINTF_INT " \t  A->cols = " MESS_PRINTF_INT " \t B->rows = " MESS_PRINTF_INT " \t B->cols = " MESS_PRINTF_INT "\n", opA[0], opB[0], m, k1, k2, n); */

    /*-----------------------------------------------------------------------------
     *  dense matrices
     *-----------------------------------------------------------------------------*/

    if ( MESS_IS_DENSE(A) && MESS_IS_DENSE(B)){
#ifdef IDEBUG
        MSG_INFO("Dense * Dense\n");
#endif
        double alpha = 1.0;
        double beta = 0.0;

        if ( MESS_IS_REAL (A) && MESS_IS_REAL(B)) {
            ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            ldc = C->ld;
            F77_GLOBAL(dgemm,DGEMM)(opAi, opBi, &m, &n, &k, &alpha, A->values, &lda, B->values, &ldb, &beta, C->values, &ldc);
        } else {
            mess_double_cpx_t al, bt;
            al = 1;
            bt = 0.0;
            ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B) )
                F77_GLOBAL(dzgemmx,DZGEMMX)(opAi, opBi, &m, &n, &k, &al, A->values, &lda, B->values_cpx, &ldb, &bt, C->values_cpx, &ldc);
            else if ( MESS_IS_COMPLEX(A) && MESS_IS_REAL(B))
                F77_GLOBAL(zdgemmx,ZDGEMMX)(opAi, opBi, &m, &n, &k, &al, A->values_cpx, &lda, B->values, &ldb, &bt, C->values_cpx, &ldc);
            else
                F77_GLOBAL(zgemm,ZGEMM)(opAi, opBi, &m, &n, &k, &al, A->values_cpx, &lda, B->values_cpx, &ldb, &bt, C->values_cpx, &ldc);
        }

    } else if ( MESS_IS_DENSE(A) && !(MESS_IS_DENSE(B))) {
#ifdef IDEBUG
        MSG_INFO("Dense * Sparse \n");
#endif
        /*-----------------------------------------------------------------------------
         *  DENSE * Sparse
         *-----------------------------------------------------------------------------*/
        mess_matrix workB = NULL;
        int cleanB = 0;
        // MSG_WARN("opA = %c opB = %c A->rows = " MESS_PRINTF_INT " \t  A->cols = " MESS_PRINTF_INT " \t B->rows = " MESS_PRINTF_INT " \t B->cols = " MESS_PRINTF_INT "\n", opA[0], opB[0], m, k1, k2, n);

        /*-----------------------------------------------------------------------------
         *  A^T * B -> C
         *-----------------------------------------------------------------------------*/
        if ( opA == MESS_OP_TRANSPOSE && opB==MESS_OP_NONE ){
            if ( !(MESS_IS_CSC(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB=1;
            } else if (MESS_IS_CSC(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsetn,DDGEMM_DENSE_SPARSETN)(&m, &n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsetn,ZZGEMM_DENSE_SPARSETN)(&m, &n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsetn,ZDGEMM_DENSE_SPARSETN)(&m, &n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsetn,DZGEMM_DENSE_SPARSETN)(&m, &n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if ( cleanB )  mess_matrix_clear(&workB);

        }
        /*-----------------------------------------------------------------------------
         *  A^H *B -> C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_HERMITIAN && opB==MESS_OP_NONE ){
            if ( !(MESS_IS_CSC(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB=1;
            } else if (MESS_IS_CSC(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsehn,DDGEMM_DENSE_SPARSEHN)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsehn,ZZGEMM_DENSE_SPARSEHN)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsehn,ZDGEMM_DENSE_SPARSEHN)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsehn,DZGEMM_DENSE_SPARSEHN)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if ( cleanB )  mess_matrix_clear(&workB);
        }
        /*-----------------------------------------------------------------------------
         *  A^T*B^T -> C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_TRANSPOSE && opB==MESS_OP_TRANSPOSE ) {
            if ( !(MESS_IS_CSR(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSR(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsett,DDGEMM_DENSE_SPARSETT)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsett,ZZGEMM_DENSE_SPARSETT)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsett,ZDGEMM_DENSE_SPARSETT)(&m, &n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsett,DZGEMM_DENSE_SPARSETT)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanB == 1 ) mess_matrix_clear(&workB);
        }
        /*-----------------------------------------------------------------------------
         *  A^H*B^T -> C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_HERMITIAN && opB==MESS_OP_TRANSPOSE) {
            if ( !(MESS_IS_CSR(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSR(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparseht,DDGEMM_DENSE_SPARSEHT)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparseht,ZZGEMM_DENSE_SPARSEHT)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparseht,ZDGEMM_DENSE_SPARSEHT)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparseht,DZGEMM_DENSE_SPARSEHT)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanB == 1 ) mess_matrix_clear(&workB);
        }
        /*-----------------------------------------------------------------------------
         *  A^H*B^H -> C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_HERMITIAN  && opB==MESS_OP_HERMITIAN ) {
            if ( !(MESS_IS_CSR(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSR(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsehh,DDGEMM_DENSE_SPARSEHH)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsehh,ZZGEMM_DENSE_SPARSEHH)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsehh,ZDGEMM_DENSE_SPARSEHH)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsehh,DZGEMM_DENSE_SPARSEHH)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanB == 1 ) mess_matrix_clear(&workB);
        }

        /*-----------------------------------------------------------------------------
         *  A^T*B^H -> C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_TRANSPOSE && opB==MESS_OP_HERMITIAN ) {
            if ( !(MESS_IS_CSR(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSR(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparseth,DDGEMM_DENSE_SPARSETH)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparseth,ZZGEMM_DENSE_SPARSETH)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparseth,ZDGEMM_DENSE_SPARSETH)(&m,&n,&k,A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparseth,DZGEMM_DENSE_SPARSETH)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanB == 1 ) mess_matrix_clear(&workB);

        }



        /*-----------------------------------------------------------------------------
         *  A*B^T->C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_NONE && opB==MESS_OP_TRANSPOSE ){
            if ( !(MESS_IS_CSR(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSR(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsent,DDGEMM_DENSE_SPARSENT)(&m, &n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsent,ZZGEMM_DENSE_SPARSENT)(&m, &n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsent,ZDGEMM_DENSE_SPARSENT)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsent,DZGEMM_DENSE_SPARSENT)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanB == 1 ) mess_matrix_clear(&workB);

        }
        /*-----------------------------------------------------------------------------
         *  A*B^H->C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_NONE  && opB==MESS_OP_HERMITIAN){
            if ( !(MESS_IS_CSR(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSR(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsenh,DDGEMM_DENSE_SPARSENH)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsenh,ZZGEMM_DENSE_SPARSENH)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsenh,ZDGEMM_DENSE_SPARSENH)(&m,&n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsenh,DZGEMM_DENSE_SPARSENH)(&m,&n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanB == 1 ) mess_matrix_clear(&workB);
        }


        /*-----------------------------------------------------------------------------
         *  A*B->C
         *-----------------------------------------------------------------------------*/
        else if ( opA == MESS_OP_NONE && opB==MESS_OP_NONE ){
            if ( !(MESS_IS_CSC(B))) {
                // MSG_ERROR("A'*B->C, CSR\n");
                ret = mess_matrix_init(&workB);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanB = 1;
            } else if (MESS_IS_CSC(B)){
                workB= B;
                cleanB = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_dense_sparsenn,DDGEMM_DENSE_SPARSENN)(&m, &n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_dense_sparsenn,ZZGEMM_DENSE_SPARSENN)(&m, &n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_dense_sparsenn,ZDGEMM_DENSE_SPARSENN)(&m, &n,&k, A->values_cpx,&lda, workB->rowptr, workB->colptr, workB->values, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_dense_sparsenn,DZGEMM_DENSE_SPARSENN)(&m, &n,&k, A->values,&lda, workB->rowptr, workB->colptr, workB->values_cpx, C->values_cpx,&ldc);
            } else {
                if (cleanB == 1) mess_matrix_clear(&workB);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if ( cleanB == 1 ) {
                mess_matrix_clear(&workB);
            }
        } else {
            MSG_ERROR("Unkown Operation\n");
            return MESS_ERROR_ARGUMENTS;

        }
    } else if ( !(MESS_IS_DENSE(A)) && MESS_IS_DENSE(B)) {
#ifdef IDEBUG
        MSG_INFO("Sparse * DENSE\n");
#endif
        /* Multiply a Sparse Matrix by a Dense matrix B
         *
         */
        mess_matrix workA;
        int cleanA = 0;

        /*-----------------------------------------------------------------------------
         *  A*B->C
         *-----------------------------------------------------------------------------*/
        if (opA == MESS_OP_NONE && opB == MESS_OP_NONE){
            if ( !(MESS_IS_CSR(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSR(A)){
                workA= A;
                cleanA = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
#ifdef IDEBUG
                MSG_INFO("Real * Real\n");
#endif
                F77_GLOBAL(ddgemm_sparse_densenn,DDGEMM_SPARSE_DENSENN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densenn,ZZGEMM_SPARSE_DENSENN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densenn,ZDGEMM_SPARSE_DENSENN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densenn,DZGEMM_SPARSE_DENSENN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);

        }

        /*-----------------------------------------------------------------------------
         *  A*B^T->C
         *-----------------------------------------------------------------------------*/
        else if (opA == MESS_OP_NONE && opB == MESS_OP_TRANSPOSE){
            if ( !(MESS_IS_CSR(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSR(A)){
                workA= A;
                cleanA = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_densent,DDGEMM_SPARSE_DENSENT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densent,ZZGEMM_SPARSE_DENSENT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densent,ZDGEMM_SPARSE_DENSENT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densent,DZGEMM_SPARSE_DENSENT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);
        }

        /*-----------------------------------------------------------------------------
         *  A*B^H -> C
         *-----------------------------------------------------------------------------*/
        else if (opA == MESS_OP_NONE && opB == MESS_OP_HERMITIAN){
            if ( !(MESS_IS_CSR(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSR); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSR(A)){
                workA= A;
                cleanA = 0;
            }
            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_densenh,DDGEMM_SPARSE_DENSENH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densenh,ZZGEMM_SPARSE_DENSENH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densenh,ZDGEMM_SPARSE_DENSENH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densenh,DZGEMM_SPARSE_DENSENH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);

        }

        /*-----------------------------------------------------------------------------
         *  A^T*B->C
         *-----------------------------------------------------------------------------*/
        else if (opA == MESS_OP_TRANSPOSE && opB == MESS_OP_NONE){
            if ( !(MESS_IS_CSC(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSC(A)){
                workA= A;
                cleanA = 0;
            }

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_densetn,DDGEMM_SPARSE_DENSETN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densetn,ZZGEMM_SPARSE_DENSETN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densetn,ZDGEMM_SPARSE_DENSETN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densetn,DZGEMM_SPARSE_DENSETN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);
        }

        /*-----------------------------------------------------------------------------
         *  A^T*B^T
         *-----------------------------------------------------------------------------*/
        else if (opA == MESS_OP_TRANSPOSE && opB == MESS_OP_TRANSPOSE){
            if ( !(MESS_IS_CSC(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSC(A)){
                workA= A;
                cleanA = 0;
            }

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_densett,DDGEMM_SPARSE_DENSETT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densett,ZZGEMM_SPARSE_DENSETT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densett,ZDGEMM_SPARSE_DENSETT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densett,DZGEMM_SPARSE_DENSETT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);
        }

        /*-----------------------------------------------------------------------------
         *  A^TB^H->C
         *-----------------------------------------------------------------------------*/
        else if (opA == MESS_OP_TRANSPOSE && opB == MESS_OP_HERMITIAN){
            // Use the property that CSC is transposed CSR
            if ( !(MESS_IS_CSC(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSC(A)){
                workA= A;
                cleanA = 0;
            }

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_denseth,DDGEMM_SPARSE_DENSETH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_denseth,ZZGEMM_SPARSE_DENSETH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_denseth,ZDGEMM_SPARSE_DENSETH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_denseth,DZGEMM_SPARSE_DENSETH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);

            /*-----------------------------------------------------------------------------
             *  A^HB->C
             *-----------------------------------------------------------------------------*/
        } else if (opA == MESS_OP_HERMITIAN && opB == MESS_OP_NONE){
            if ( !(MESS_IS_CSC(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSC(A)){
                workA= A;
                cleanA = 0;
            }

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_densehn,DDGEMM_SPARSE_DENSEHN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densehn,ZZGEMM_SPARSE_DENSEHN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densehn,ZDGEMM_SPARSE_DENSEHN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densehn,DZGEMM_SPARSE_DENSEHN)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);
        }
        /*-----------------------------------------------------------------------------
         *  A^HB^T -> C
         *-----------------------------------------------------------------------------*/
        else if (opA == MESS_OP_HERMITIAN && opB == MESS_OP_TRANSPOSE){
            if ( !(MESS_IS_CSC(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSC(A)){
                workA= A;
                cleanA = 0;
            }

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_denseht,DDGEMM_SPARSE_DENSEHT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_denseht,ZZGEMM_SPARSE_DENSEHT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_denseht,ZDGEMM_SPARSE_DENSEHT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_denseht,DZGEMM_SPARSE_DENSEHT)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);

        } else if (opA == MESS_OP_HERMITIAN && opB == MESS_OP_HERMITIAN){
            if ( !(MESS_IS_CSC(A))) {
                ret = mess_matrix_init(&workA);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(A,workA, MESS_CSC); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_matrix_convert);
                cleanA = 1;
            } else if (MESS_IS_CSC(A)){
                workA= A;
                cleanA = 0;
            }

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B) ){
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
            }
            ldc = C->ld;

            if (MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                F77_GLOBAL(ddgemm_sparse_densehh,DDGEMM_SPARSE_DENSEHH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values,&ldb, C->values,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(zzgemm_sparse_densehh,ZZGEMM_SPARSE_DENSEHH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_COMPLEX(A) && MESS_IS_REAL(B)){
                F77_GLOBAL(zdgemm_sparse_densehh,ZDGEMM_SPARSE_DENSEHH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values_cpx, B->values,&ldb, C->values_cpx,&ldc);
            } else if (MESS_IS_REAL(A) && MESS_IS_COMPLEX(B)){
                F77_GLOBAL(dzgemm_sparse_densehh,DZGEMM_SPARSE_DENSEHH)(&m, &n,&k, workA->rowptr, workA->colptr,workA->values, B->values_cpx,&ldb, C->values_cpx,&ldc);
            } else {
                if (cleanA == 1) mess_matrix_clear(&workA);
                MSG_ERROR("unknown data type combination.\n");
                return (MESS_ERROR_DATATYPE);
            }
            if (cleanA == 1 ) mess_matrix_clear(&workA);

        } else {
            MSG_ERROR("Unknown Operation combination\n");
            return (MESS_ERROR_ARGUMENTS);
        }

    }

    /*-----------------------------------------------------------------------------
     *  Both matrices are sparse
     *-----------------------------------------------------------------------------*/
    else if ( !(MESS_IS_DENSE(A)) && !(MESS_IS_DENSE(B))) {
#ifdef MESS_HAVE_CSPARSE
        mess_matrix workA, workB;
        mess_int_t cleanA=0, cleanB=0;

        /*-----------------------------------------------------------------------------
         *  prepare work A
         *-----------------------------------------------------------------------------*/
        if (opA==MESS_OP_NONE){
            workA = A;
        }else if (opA==MESS_OP_TRANSPOSE){
            ret = mess_matrix_init(&workA);                                     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_ctranspose(A,workA);                              FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
            if (MESS_IS_COMPLEX(workA)){
                ret = mess_matrix_conj(workA);
                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
            }
            cleanA = 1;
        }else{
            ret = mess_matrix_init(&workA);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_ctranspose(A,workA);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
            cleanA = 1;
        }

        /*-----------------------------------------------------------------------------
         *  prepare work B
         *-----------------------------------------------------------------------------*/
        if (opB==MESS_OP_NONE){
            workB = B;
        }else if (opB==MESS_OP_TRANSPOSE){
            ret = mess_matrix_init(&workB);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_ctranspose(B,workB);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
            if (MESS_IS_COMPLEX(workB)){
                ret = mess_matrix_conj(workB);                                      FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_conj);
            }
            cleanB = 1;
        }else{
            ret = mess_matrix_init(&workB);                                         FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_ctranspose(B,workB);                                  FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_ctranspose);
            cleanB = 1;
        }


        if(MESS_IS_COMPLEX(workA) || MESS_IS_COMPLEX(workB)){
            ret = mess_matrix_tocomplex(workA);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
            ret = mess_matrix_tocomplex(workB);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_tocomplex);
        }


        /*-----------------------------------------------------------------------------
         *  convert to csparse and multiply
         *-----------------------------------------------------------------------------*/
        if (MESS_IS_REAL(workA) && MESS_IS_REAL(workB)){

            cs_dl *csA=NULL, *csB=NULL, *csC=NULL;

            ret = mess_matrix_to_csparse_dl(workA,&csA);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_to_csparse_dl);
            ret = mess_matrix_to_csparse_dl(workB,&csB);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_to_csparse_dl);
            csC = cs_dl_multiply(csA,csB);
            ret = mess_matrix_from_csparse_dl(csC,C);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_from_csparse_dl);
            cs_dl_spfree(csA); cs_dl_spfree(csB); cs_dl_spfree(csC);
        }else{
            cs_cl *csA=NULL, *csB=NULL, *csC=NULL;

            ret = mess_matrix_to_csparse_cl(workA,&csA);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_to_csparse_cl);
            ret = mess_matrix_to_csparse_cl(workB,&csB);                FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_to_csparse_cl);
            csC = cs_cl_multiply(csA,csB);
            ret = mess_matrix_from_csparse_cl(csC,C);                   FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_from_csparse_cl);
            cs_cl_spfree(csA); cs_cl_spfree(csB); cs_cl_spfree(csC);
        }

        /*-----------------------------------------------------------------------------
         *  clean
         *-----------------------------------------------------------------------------*/
        if(cleanA){
            ret = mess_matrix_clear(&workA);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        }

        if(cleanB){
            ret = mess_matrix_clear(&workB);        FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        }

#else

        int cleanB=0;
        mess_int_t i;
        mess_vector w;
        mess_vector w2;
        mess_matrix workB;
        if ( opB == MESS_OP_NONE ) {
            if ( MESS_IS_CSC(B)){
                cleanB = 0;
                workB = B;
            } else {
                ret = mess_matrix_init(&workB);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB,MESS_CSC);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
                cleanB=1;
            }
            if ( MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
                ret = mess_vector_init(&w);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_init(&w2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w,k,MESS_REAL);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
                ret = mess_vector_alloc(w2, m, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
                ret = mess_vector_init(&w);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_init(&w2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w,k,MESS_COMPLEX);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
                ret = mess_vector_alloc(w2, m, MESS_COMPLEX);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
            }

            for ( i = 0; i  < n; i++) {
                ret = mess_matrix_getcol(workB, i, w);                              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getcol);
                ret = mess_matrix_mvp(opA, A, w, w2);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_setcol(C,i,w2);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            }

            if ( cleanB ){
                mess_matrix_clear(&workB);
            }
            mess_vector_clear(&w);
            mess_vector_clear(&w2);

        } else if (opB == MESS_OP_TRANSPOSE ){
            if ( MESS_IS_CSR(B)){
                cleanB = 0;
                workB = B;
            } else {
                ret = mess_matrix_init(&workB);                                     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB,MESS_CSR);                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
                cleanB=1;
            }
            if ( MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
                ret = mess_vector_init(&w);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_init(&w2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w,k,MESS_REAL);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w2, m, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                for ( i = 0; i  < workB->rows; i++) {
                    ret = mess_matrix_getrow(workB, i, w);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getrow);
                    ret = mess_matrix_mvp(opA, A, w, w2);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_setcol(C,i,w2);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
                }
                mess_vector_clear(&w);
                mess_vector_clear(&w2);

            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
                ret = mess_vector_init(&w);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_init(&w2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w,k,MESS_COMPLEX);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
                ret = mess_vector_alloc(w2, m, MESS_COMPLEX);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
                for ( i = 0; i  < workB->rows; i++) {
                    ret = mess_matrix_getrow(workB, i, w);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getrow);
                    ret = mess_matrix_mvp(opA, A, w, w2);                           FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                    ret = mess_matrix_setcol(C,i,w2);                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
                }
                mess_vector_clear(&w);
                mess_vector_clear(&w2);

            }

            if ( cleanB ){
                mess_matrix_clear(&workB);
            }

        } else if (opB == MESS_OP_HERMITIAN ){
            if ( MESS_IS_CSR(B)){
                cleanB = 0;
                workB = B;
                // ret = mess_matrix_conj(workB);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
            } else {
                ret = mess_matrix_init(&workB);                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                ret = mess_matrix_convert(B,workB,MESS_CSR);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
                // ret = mess_matrix_conj(workB);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
                cleanB=1;
            }
            if ( MESS_IS_REAL(A) && MESS_IS_REAL(B)) {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_REAL);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
                ret = mess_vector_init(&w);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_init(&w2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w,k,MESS_REAL);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w2, m, MESS_REAL);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
            } else {
                ret = mess_matrix_alloc (C, m, n, m*n, MESS_DENSE, MESS_COMPLEX);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
                ret = mess_vector_init(&w);                                         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_init(&w2);                                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
                ret = mess_vector_alloc(w,k,MESS_COMPLEX);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
                ret = mess_vector_alloc(w2, m, MESS_COMPLEX);                       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
            }

            for ( i = 0; i  < workB->rows; i++) {
                ret = mess_matrix_getrow(workB, i, w);              FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_getrow);
                ret = mess_vector_conj(w);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_conj);
                ret = mess_matrix_mvp(opA, A, w, w2);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_mvp);
                ret = mess_matrix_setcol(C,i,w2);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_setcol);
            }

            if ( cleanB ){
                mess_matrix_clear(&workB);
            } else {
                // ret = mess_matrix_conj(workB);                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_conj);
            }
            mess_vector_clear(&w);
            mess_vector_clear(&w2);
        }
#endif
    }
    return (0);
}

