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
 * @file lib/direct/singlesolver/newlu.c
 * @brief Easy left looking sparse LU.
 * @author @koehlerm
 */

#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>
#include    <string.h>
#include    "mess/mess.h"
#include    "mess/error_macro.h"

#define FLIP(x) (-(x)-2)
#define UNFLIP(x) (((x) < 0 )? FLIP(x):(x))
#define MARKED(w,x) ( w[(x)] < 0)
#define MARK(w,x)  {  w[(x)] = FLIP(w[(x)]); }


static mess_int_t *perm_inv (mess_int_t const *p,mess_int_t n)
{
    MSG_FNAME(__func__);
    mess_int_t k, *pinv ;
    if (!p) return (NULL) ;
    mess_try_alloc2(pinv, mess_int_t*, n*sizeof (mess_int_t)) ;
    if (!pinv) return (NULL) ;
    for (k = 0 ; k < n ; k++) pinv [p [k]] = k ;
    return (pinv) ;
}


struct newlu_data {
    mess_matrix L;
    mess_matrix U;
    mess_int_t *pinv;
    mess_int_t *q;
    mess_datatype_t data_type;
    mess_int_t rows;
    mess_int_t cols;
};


static int newlu_getL(void *data, mess_matrix L ){
    struct newlu_data *sol = ( struct newlu_data*) data;
    return mess_matrix_convert(sol->L, L, MESS_CSR);
}

static int newlu_getU( void *data, mess_matrix U ) {
    struct newlu_data *sol = ( struct newlu_data*) data;
    return mess_matrix_convert(sol->U, U, MESS_CSR);
}


static int newlu_getpermp( void *data, mess_int_t * p  ) {
    MSG_FNAME(__func__);
    struct newlu_data *sol = ( struct newlu_data*) data;
    mess_int_t i = 0 ;
    mess_int_t *tmp, *t2;

    mess_try_alloc(tmp, mess_int_t*, sizeof(mess_int_t) * sol->L->rows);
    for (i = 0 ; i < sol->L->rows; i++) {
        tmp [i ] = (sol->pinv == NULL)?i:(sol->pinv[i]);
    }
    t2 = perm_inv(tmp, sol->L->rows);
    for (i = 0 ; i < sol->L->rows; i++) {
        p[i]=t2[i];
    }
    mess_free(tmp);
    mess_free(t2);
    return 0;
}

static int newlu_getpermq( void *data, mess_int_t *q) {
    struct newlu_data *sol = ( struct newlu_data*) data;
    mess_int_t i;
    for ( i=0;i<sol->L->rows ; i++ ) {
        q[i] = (sol->q == NULL)? i: (mess_int_t)( sol->q[i]);
    }

    return 0;
}



static int newlu_clear( void * data) {
    struct newlu_data *sol = (struct newlu_data*) data;
    mess_matrix_clear(&(sol->L));
    mess_matrix_clear(&(sol->U));
    mess_free(sol->pinv);
    mess_free(sol->q);
    mess_free(sol);
    return 0;
}

static int lsolve ( mess_matrix L, double * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = 0;  i < L->cols; i++ ) {
        x[i] /= L->values[L->colptr[i]];
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[L->rowptr[j]] -= L->values[j] * x[i];
        }
    }
    return 0;
}

static int lsolverc ( mess_matrix L, mess_double_cpx_t* x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = 0;  i < L->cols; i++ ) {
        x[i] /= L->values[L->colptr[i]];
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[L->rowptr[j]] -= L->values[j] * x[i];
        }
    }
    return 0;
}


static int lsolvec ( mess_matrix L, mess_double_cpx_t* x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = 0;  i < L->cols; i++ ) {
        x[i] /= L->values_cpx[L->colptr[i]];
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[L->rowptr[j]] -= L->values_cpx[j] * x[i];
        }
    }
    return 0;
}


static int ltsolve(mess_matrix L, double *x ){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = L->cols-1;  i >= 0 ; i-- ) {
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[i] -= L->values[j] * x[L->rowptr[j]];
        }
        x[i] /= L->values[L->colptr[i]];
    }
    return 0;

}

static int ltsolverc(mess_matrix L, mess_double_cpx_t *x ){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = L->cols-1;  i >= 0 ; i-- ) {
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[i] -= L->values[j] * x[L->rowptr[j]];
        }
        x[i] /= L->values[L->colptr[i]];
    }
    return 0;

}


static int ltsolvec(mess_matrix L, mess_double_cpx_t *x ){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = L->cols-1;  i >= 0 ; i-- ) {
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[i] -= L->values_cpx[j] * x[L->rowptr[j]];
        }
        x[i] /= L->values_cpx[L->colptr[i]];
    }
    return 0;

}

static int lhsolvec(mess_matrix L, mess_double_cpx_t *x ){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(L);
    mess_check_nullpointer(x);
    mess_check_square(L);
    mess_check_csc(L);


    for ( i = L->cols-1;  i >= 0 ; i-- ) {
        for ( j = L->colptr[i]+1; j < L->colptr[i+1];  j++){
            x[i] -= conj(L->values_cpx[j]) * x[L->rowptr[j]];
        }
        x[i] /= conj(L->values_cpx[L->colptr[i]]);
    }
    return 0;

}



static int usolve ( mess_matrix U, double * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = U->cols-1;  i>=0 ; i--) {
        x[i] /= U->values[U->colptr[i+1]-1];
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[U->rowptr[j]] -= U->values[j] * x[i];
        }
    }
    return 0;
}

static int usolverc ( mess_matrix U, mess_double_cpx_t  * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = U->cols-1;  i>=0 ; i--) {
        x[i] /= U->values[U->colptr[i+1]-1];
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[U->rowptr[j]] -= U->values[j] * x[i];
        }
    }
    return 0;
}





static int usolvec ( mess_matrix U, mess_double_cpx_t * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = U->cols-1;  i>=0 ; i--) {
        x[i] /= U->values_cpx[U->colptr[i+1]-1];
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[U->rowptr[j]] -= U->values_cpx[j] * x[i];
        }
    }
    return 0;
}


static int utsolve ( mess_matrix U, double * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = 0;  i<U->rows; i++) {
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[i] -= U->values[j] * x[U->rowptr[j]];
        }
        x[i] /= U->values[U->colptr[i+1]-1];
    }
    return 0;
}

static int utsolverc ( mess_matrix U, mess_double_cpx_t * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = 0;  i<U->rows; i++) {
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[i] -= U->values[j] * x[U->rowptr[j]];
        }
        x[i] /= U->values[U->colptr[i+1]-1];
    }
    return 0;
}


static int utsolvec ( mess_matrix U, mess_double_cpx_t * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = 0;  i<U->rows; i++) {
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[i] -= U->values_cpx[j] * x[U->rowptr[j]];
        }
        x[i] /= U->values_cpx[U->colptr[i+1]-1];
    }
    return 0;
}

static int uhsolvec ( mess_matrix U, mess_double_cpx_t * x){
    MSG_FNAME(__func__);
    mess_int_t j,i;
    mess_check_nullpointer(U);
    mess_check_nullpointer(x);
    mess_check_square(U);
    mess_check_csc(U);

    for ( i = 0;  i<U->rows; i++) {
        for ( j = U->colptr[i]; j < U->colptr[i+1]-1;  j++){
            x[i] -= conj(U->values_cpx[j]) * x[U->rowptr[j]];
        }
        x[i] /= conj(U->values_cpx[U->colptr[i+1]-1]);
    }
    return 0;
}




static int newlu_solve(void *data, mess_vector b, mess_vector x){
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_vector treal;
    if ( MESS_IS_REAL(sol) && MESS_IS_REAL(b)) {
        mess_vector_toreal(x);
        mess_vector_init(&treal);
        mess_vector_alloc(treal, sol->L->rows, MESS_REAL);
        mess_vector_iperm(b, sol->pinv, treal);
        lsolve(sol->L, treal->values);
        usolve(sol->U, treal->values);
        mess_vector_iperm(treal, sol->q, x);
        mess_vector_clear(&treal);
    } else if (MESS_IS_REAL(sol) && MESS_IS_COMPLEX(b)) {
        mess_vector_tocomplex(x);
        mess_vector_init(&treal);
        mess_vector_alloc(treal, sol->L->rows, MESS_COMPLEX);
        mess_vector_iperm(b, sol->pinv, treal);
        lsolverc(sol->L, treal->values_cpx);
        usolverc(sol->U, treal->values_cpx);
        mess_vector_iperm(treal, sol->q, x);
        mess_vector_clear(&treal);
    } else {
        mess_vector_tocomplex(x);
        mess_vector_init(&treal);
        mess_vector_alloc(treal, sol->L->rows, MESS_COMPLEX);
        mess_vector_iperm(b, sol->pinv, treal);
        lsolvec(sol->L, treal->values_cpx);
        usolvec(sol->U, treal->values_cpx);
        mess_vector_iperm(treal, sol->q, x);
        mess_vector_clear(&treal);
    }

    return 0;
}


static int newlu_solvem(void *data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_matrix workb, tmp;
    int conv = -1;
    mess_int_t i,j;

    MESS_MATRIX_CHECKFORMAT(b, workb, conv, MESS_DENSE);
    /* MESS_MATRIX_RESET(x);
       mess_matrix_alloc(x, b->rows, b->cols, b->rows*b->cols, MESS_DENSE, b->data_type);  */
    mess_matrix_init(&tmp);
    mess_matrix_alloc(tmp, b->rows, b->cols, b->rows*b->cols, MESS_DENSE, b->data_type);

    if ( MESS_IS_REAL(b)) {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                tmp->values[j*tmp->ld+( sol->pinv ? sol->pinv[i]:i)] = workb->values[j*workb->ld+i];
            }
        }
    } else {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                tmp->values_cpx[j*tmp->ld+( sol->pinv ? sol->pinv[i]:i)] = workb->values_cpx[j*workb->ld+i];
            }
        }

    }

    if ( MESS_IS_COMPLEX(sol)) {
        mess_matrix_tocomplex(tmp);
    }
    if ( MESS_IS_REAL(sol) && MESS_IS_REAL(b)) {
        for (j=0; j<b->cols; j++){
            lsolve(sol->L, tmp->values+j*tmp->ld);
            usolve(sol->U, tmp->values+j*tmp->ld);
        }
    } else if (MESS_IS_REAL(sol) && MESS_IS_COMPLEX(b)) {
        for (j=0; j<b->cols; j++){
            lsolverc(sol->L, tmp->values_cpx+j*tmp->ld);
            usolverc(sol->U, tmp->values_cpx+j*tmp->ld);
        }
    } else {
        for (j=0; j<b->cols; j++){
            lsolvec(sol->L, tmp->values_cpx+j*tmp->ld);
            usolvec(sol->U, tmp->values_cpx+j*tmp->ld);
        }
    }

    if(mess_matrix_need_alloc(x,tmp->rows, tmp->cols, tmp->rows*tmp->cols, MESS_DENSE, tmp->data_type)){
        MESS_MATRIX_RESET(x);
        mess_matrix_alloc(x, tmp->rows, tmp->cols, tmp->rows*tmp->cols, MESS_DENSE, tmp->data_type);
    }

    if ( MESS_IS_REAL(tmp)) {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                x->values[j*x->ld+( sol->q ? sol->q[i]:i)] = tmp->values[j*tmp->ld+i];
            }
        }
    } else {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                x->values_cpx[j*x->ld+( sol->q ? sol->q[i]:i)] = tmp->values_cpx[j*tmp->ld+i];
            }
        }

    }
    if ( conv == 0 ) {
        mess_matrix_clear(&workb);
    }
    mess_matrix_clear(&tmp);
    return 0;
}


static int newlu_solvet(void *data, mess_vector b, mess_vector x){
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_vector treal;
    MESS_INIT_VECTORS(&treal);
    if ( MESS_IS_REAL(sol) && MESS_IS_REAL(b) ){
        mess_vector_alloc(treal, sol->L->rows, MESS_REAL);
        mess_vector_perm(b, sol->q, treal);
        utsolve(sol->U, treal->values);
        ltsolve(sol->L, treal->values);
        mess_vector_perm(treal, sol->pinv, x);
        mess_vector_clear(&treal);
    } else if (MESS_IS_REAL(sol) && MESS_IS_COMPLEX(b)) {
        mess_vector_tocomplex(x);
        mess_vector_alloc(treal, sol->L->rows, MESS_COMPLEX);
        mess_vector_perm(b, sol->q, treal);
        utsolverc(sol->U, treal->values_cpx);
        ltsolverc(sol->L, treal->values_cpx);
        mess_vector_perm(treal, sol->pinv, x);
        mess_vector_clear(&treal);
    } else {
        mess_vector_tocomplex(x);
        mess_vector_alloc(treal, sol->L->rows, MESS_COMPLEX);
        mess_vector_perm(b, sol->q, treal);
        utsolvec(sol->U, treal->values_cpx);
        ltsolvec(sol->L, treal->values_cpx);
        mess_vector_perm(treal, sol->pinv, x);
        mess_vector_clear(&treal);
    }

    return 0;
}


static int newlu_solvemt(void *data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_matrix workb, tmp;
    int conv = -1;
    mess_int_t i,j;

    MESS_MATRIX_CHECKFORMAT(b, workb, conv, MESS_DENSE);
    mess_matrix_init(&tmp);
    mess_matrix_alloc(tmp, b->rows, b->cols, b->rows*b->cols, MESS_DENSE, b->data_type);

    if (MESS_IS_COMPLEX(sol) || MESS_IS_COMPLEX(b)){
        mess_matrix_tocomplex(tmp);
    }
    if ( MESS_IS_REAL(b)) {
        if ( MESS_IS_COMPLEX(tmp) ) {
            for (j=0; j < b->cols; j++){
                for(i=0; i< b->rows; i++){
                    tmp->values_cpx[j*tmp->ld+i] = workb->values[j*workb->ld+(sol->q? sol->q[i] : i )];
                }
            }
        } else {
            for (j=0; j < b->cols; j++){
                for(i=0; i< b->rows; i++){
                    tmp->values[j*tmp->ld+i] = workb->values[j*workb->ld+(sol->q? sol->q[i] : i )];
                }
            }
        }
    } else {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                tmp->values_cpx[j*tmp->ld+i] = workb->values_cpx[j*workb->ld+(sol->q? sol->q[i] : i ) ];
            }
        }

    }

    if ( MESS_IS_REAL(sol) && MESS_IS_REAL(b)) {
        for (j=0; j<b->cols; j++){
            utsolve(sol->U, tmp->values+j*tmp->ld);
            ltsolve(sol->L, tmp->values+j*tmp->ld);
        }
    } else if (MESS_IS_REAL(sol) && MESS_IS_COMPLEX(b)) {
        for (j=0; j<b->cols; j++){
            utsolverc(sol->U, tmp->values_cpx+j*tmp->ld);
            ltsolverc(sol->L, tmp->values_cpx+j*tmp->ld);
        }
    } else {
        for (j=0; j<b->cols; j++){
            utsolvec(sol->U, tmp->values_cpx+j*tmp->ld);
            ltsolvec(sol->L, tmp->values_cpx+j*tmp->ld);
        }
    }

    if(mess_matrix_need_alloc(x,tmp->rows, tmp->cols, tmp->rows*tmp->cols, MESS_DENSE, tmp->data_type)){
        MESS_MATRIX_RESET(x);
        mess_matrix_alloc(x, tmp->rows, tmp->cols, tmp->rows*tmp->cols, MESS_DENSE, tmp->data_type);
    }

    if ( MESS_IS_REAL(tmp)) {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                x->values[j*x->ld+i] = tmp->values[j*tmp->ld +(sol->pinv? sol->pinv[i] : i )];
            }
        }
    } else {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                x->values_cpx[j*x->ld+i] = tmp->values_cpx[j*tmp->ld+(sol->pinv? sol->pinv[i] : i ) ];
            }
        }

    }
    if ( conv == 0 ) {
        mess_matrix_clear(&workb);
    }
    mess_matrix_clear(&tmp);

    return 0;
}



static int newlu_solveh(void *data, mess_vector b, mess_vector x){
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_vector treal;
    if ( MESS_IS_REAL(sol)){
        return newlu_solvet(data, b, x);
    } else {
        mess_vector_tocomplex(x);
        mess_vector_init(&treal);
        mess_vector_alloc(treal, sol->L->rows, MESS_COMPLEX);
        mess_vector_perm(b, sol->q, treal);
        uhsolvec(sol->U, treal->values_cpx);
        lhsolvec(sol->L, treal->values_cpx);
        mess_vector_perm(treal, sol->pinv, x);
        mess_vector_clear(&treal);
    }

    return 0;
}


static int newlu_solvemh(void *data, mess_matrix b, mess_matrix x){
    MSG_FNAME(__func__);
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_matrix workb, tmp;
    int conv = -1;
    mess_int_t i,j;

    if ( MESS_IS_REAL(sol)){
        return newlu_solvemt(data, b, x);
    }

    MESS_MATRIX_CHECKFORMAT(b, workb, conv, MESS_DENSE);
    mess_matrix_init(&tmp);
    mess_matrix_alloc(tmp, b->rows, b->cols, b->rows*b->cols, MESS_DENSE, MESS_COMPLEX);

    if ( MESS_IS_REAL(b)) {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                tmp->values_cpx[j*tmp->ld+i] = workb->values[j*workb->ld+(sol->q? sol->q[i] : i )];
            }
        }
    } else {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                tmp->values_cpx[j*tmp->ld+i] = workb->values_cpx[j*workb->ld+(sol->q? sol->q[i] : i ) ];
            }
        }

    }

    for (j=0; j<b->cols; j++){
        uhsolvec(sol->U, tmp->values_cpx+j*tmp->ld);
        lhsolvec(sol->L, tmp->values_cpx+j*tmp->ld);
    }

    if(mess_matrix_need_alloc(x,tmp->rows, tmp->cols, tmp->rows*tmp->cols, MESS_DENSE, tmp->data_type)){
        MESS_MATRIX_RESET(x);
        mess_matrix_alloc(x, tmp->rows, tmp->cols, tmp->rows*tmp->cols, MESS_DENSE, tmp->data_type);
    }

    if ( MESS_IS_REAL(tmp)) {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                x->values[j*x->ld+i] = tmp->values[j*tmp->ld+(sol->pinv? sol->pinv[i] : i )];
            }
        }
    } else {
        for (j=0; j < b->cols; j++){
            for(i=0; i< b->rows; i++){
                x->values_cpx[j*x->ld+i] = tmp->values_cpx[j*tmp->ld+(sol->pinv? sol->pinv[i] : i ) ];
            }
        }

    }

    if ( conv == 0 ) {
        mess_matrix_clear(&workb);
    }
    mess_matrix_clear(&tmp);

    return 0;
}


static int newlu_inverse ( void *data, mess_matrix inv )
{
    MSG_FNAME(__func__);
    struct newlu_data * sol = (struct newlu_data*) data;
    mess_matrix eye;
    int ret = 0;

    mess_check_nullpointer(data);
    mess_check_nullpointer(inv);
    MESS_MATRIX_RESET(inv);

    ret = mess_matrix_init(&eye);                                   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_eye(eye, sol->rows, sol->cols, MESS_DENSE);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
    ret = newlu_solvem(data, eye, inv);                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), newlu_solvem);
    mess_matrix_clear(&eye);
    return 0;
}       /* -----  end of function csparse_inverse  ----- */



/**
 * @brief Compute a sparse left-looking LU decomposition.
 * @param[in] A      input matrix to factorize
 * @param[out] sol solver object containing the left-looking sparse LU decomposition
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_direct_create_sparse_lu function computes a left looking LU decomposition of a sparse matrix.
 */
int mess_direct_create_sparse_lu(mess_matrix A, mess_direct sol) {
    MSG_FNAME(__func__);
    mess_matrix work;
    int conv = -1;
    struct newlu_data * solver;
    mess_int_t lnz, unz;
    mess_int_t *chi, *iperm;
    double *x;
    mess_double_cpx_t *xc;
    mess_int_t n,i, k,p, top;
    int ret = 0;
    double  pivot, npivot;
    mess_double_cpx_t pivotc = 0 ;
    mess_int_t ipivot, col;
    mess_matrix L, U;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(A);
    mess_check_nullpointer(sol);
    mess_check_square(A);
    mess_check_real_or_complex(A);
    MESS_MATRIX_CHECKFORMAT(A, work, conv, MESS_CSC);

    if ( MESS_IS_REAL(A)) {

        /*-----------------------------------------------------------------------------
         *  real case
         *-----------------------------------------------------------------------------*/
        /*-----------------------------------------------------------------------------
         *  allocate arrays
         *-----------------------------------------------------------------------------*/
        mess_try_alloc(solver, struct newlu_data*, sizeof(struct newlu_data));
        ret = mess_matrix_init(&(solver->L));               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&(solver->U));               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_alloc(solver->L, work->rows, work->cols, work->nnz*2, MESS_CSC, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(solver->U, work->rows, work->cols, work->nnz*2, MESS_CSC, MESS_REAL);       FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        solver->data_type = MESS_REAL;

        L = solver->L;
        U = solver->U;
        n = work->rows;
        mess_try_alloc(chi, mess_int_t*, sizeof(mess_int_t) * 2 *n);
        mess_try_alloc(iperm, mess_int_t *, sizeof(mess_int_t) *n);
        mess_try_alloc(x, double *, sizeof(double) *n);
        mess_try_alloc(solver->q, mess_int_t *, sizeof(mess_int_t)* n);

        for ( i = 0; i < n ; i++ ) {
            x[i] = 0;
            L->colptr [i] = 0;
            iperm[i] = -1;
        }
        L->colptr [n] =0;
        lnz = unz = 0;
#ifdef MESS_HAVE_AMD
        mess_matrix_reorder_amd(work, solver->q);
#else
        for (i=0; i<n;i++) solver->q[i]=i;
#endif
        /*-----------------------------------------------------------------------------
         *  for k=0...n-1
         *-----------------------------------------------------------------------------*/
        for (k=0; k < n; k++){
            // MSG_I<F12>:NFO("k = " MESS_PRINTF_INT "\n", k);
            L->colptr [k] = lnz;
            U->colptr [k] = unz;
            if ( lnz + n > L->nnz) {
                ret = mess_matrix_realloc(L, L->rows, L->cols, 2*L->nnz+n);
                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);
            }
            if ( unz + n > U->nnz) {
                ret = mess_matrix_realloc(U, U->rows, U->cols, 2*U->nnz+n);
                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);
            }
            // to change if we use column permutations too
            col = solver->q[k] ;
            ret = mess_direct_splsolve ( L, work, col, &top, chi, x, iperm);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_splsolve);

            ipivot = -1;
            pivot = -1;

            /*-----------------------------------------------------------------------------
             *  U(1:k-1,k) = x(1:k-1) + pivot search
             *-----------------------------------------------------------------------------*/
            for ( p = top; p < n; p++){
                i=chi[p];
                // MSG_INFO("x[" MESS_PRINTF_INT "] = %lg \t\t iperm[" MESS_PRINTF_INT "] = " MESS_PRINTF_INT "\n",chi[p], x[i], i, iperm[i]);

                if ( iperm[i] < 0) {    // Die zeile wurde noch keiner pivot prüfung unterzogen, ist also definitiv nicht zu U gehörend
                    if ( ( npivot = fabs(x[i])) > pivot) { // neues pivot gefunden
                        pivot = npivot;
                        ipivot = i;
                    }
                } else {    // die zeile wurde schon passend getauscht und gehört zu U
                    // MSG_INFO("\t\t-> in U\n");
                    U->rowptr[unz] = iperm[i];
                    U->values[unz++] = x[i];
                }
            }
            if ( ipivot == -1 || pivot < 0){
                MSG_ERROR("singular matrix\n");
                return MESS_ERROR_SINGULAR;
            }
            if (iperm[col] < 0 && fabs(x[col]) >= pivot) ipivot = col;  // nur pivotisieren wenn es auch sinnmacht, sonst bleibt die diagonale
            // ipivot = col;
            pivot =x[ipivot];   // das pivot ohne betrag holen
            U->values[unz] = pivot;
            U->rowptr[unz++] = k;   // und rein damit ins U
            iperm[ipivot] = k;
            L->values[lnz] = 1;
            L->rowptr[lnz++]= ipivot;
            /*-----------------------------------------------------------------------------
             *  L(k+1:n,k) = x(k+1:n)/pivot
             *-----------------------------------------------------------------------------*/
            for ( p = top; p < n ; p++) {
                i = chi[p];
                if (iperm[i] < 0){
                    L->rowptr[lnz] = i;
                    L->values[lnz++] = x[i]/pivot;
                }
                x[i] = 0;
            }
        }

        L->colptr[n] =lnz;
        U->colptr[n] =unz;
        ret = mess_matrix_realloc(L, L->rows, L->cols, lnz);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);
        ret = mess_matrix_realloc(U, U->rows, U->cols, unz);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);

        for ( p=0; p<lnz;p++ ) { L->rowptr[p] = iperm[L->rowptr[p]];}
        mess_free(x);
        mess_free(chi);
    } else {
        /*-----------------------------------------------------------------------------
         *  complex case
         *-----------------------------------------------------------------------------*/
        /*-----------------------------------------------------------------------------
         *  allocate arrays
         *-----------------------------------------------------------------------------*/
        mess_try_alloc(solver, struct newlu_data*, sizeof(struct newlu_data));
        ret = mess_matrix_init(&(solver->L));               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&(solver->U));               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_alloc(solver->L, work->rows, work->cols, work->nnz*2, MESS_CSC, MESS_COMPLEX);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        ret = mess_matrix_alloc(solver->U, work->rows, work->cols, work->nnz*2, MESS_CSC, MESS_COMPLEX);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        solver->data_type = MESS_COMPLEX;
        L = solver->L;
        U = solver->U;
        n = work->rows;
        mess_try_alloc(chi, mess_int_t*, sizeof(mess_int_t) * 2 *n);
        mess_try_alloc(iperm, mess_int_t *, sizeof(mess_int_t) *n);
        mess_try_alloc(xc, mess_double_cpx_t  *, sizeof(mess_double_cpx_t) *n);
        mess_try_alloc(solver->q, mess_int_t *, sizeof(mess_int_t)* n);

        for ( i = 0; i < n ; i++ ) {
            xc[i] = 0;
            L->colptr [i] = 0;
            iperm[i] = -1;
        }
        L->colptr [n] =0;
        lnz = unz = 0;
#ifdef MESS_HAVE_AMD
        mess_matrix_reorder_amd(work, solver->q);
#else
        for ( i =0 ; i < n ; i++) solver->q[i]=i;
#endif
        /*-----------------------------------------------------------------------------
         *  for k=0...n-1
         *-----------------------------------------------------------------------------*/
        for (k=0; k < n; k++){
            double ap,npivotc;
            // MSG_I<F12>:NFO("k = " MESS_PRINTF_INT "\n", k);
            L->colptr [k] = lnz;
            U->colptr [k] = unz;
            if ( lnz + n > L->nnz) {
                ret = mess_matrix_realloc(L, L->rows, L->cols, 2*L->nnz+n);
                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);
            }
            if ( unz + n > U->nnz) {
                ret = mess_matrix_realloc(U, U->rows, U->cols, 2*U->nnz+n);
                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);
            }
            // to change if we use column permutations too
            col = solver->q[k] ;
            ret = mess_direct_splsolvec ( L, work, col, &top, chi, xc, iperm);
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_direct_splsolvec);

            ipivot = -1;
            ap = 0;

            /*-----------------------------------------------------------------------------
             *  U(1:k-1,k) = x(1:k-1) + pivot search
             *-----------------------------------------------------------------------------*/
            for ( p = top; p < n; p++){
                i=chi[p];
                // MSG_INFO("x[" MESS_PRINTF_INT "] = %lg \t\t iperm[" MESS_PRINTF_INT "] = " MESS_PRINTF_INT "\n",chi[p], x[i], i, iperm[i]);

                if ( iperm[i] < 0) {    // Die zeile wurde noch keiner pivot prüfung unterzogen, ist also definitiv nicht zu U gehörend
                    if ( ( npivotc = cabs(xc[i])) > ap) { // neues pivot gefunden
                        ap = npivotc;
                        ipivot = i;
                    }
                } else {    // die zeile wurde schon passend getauscht und gehört zu U
                    // MSG_INFO("\t\t-> in U\n");
                    U->rowptr[unz] = iperm[i];
                    U->values_cpx[unz++] = xc[i];
                }
            }
            if ( ipivot == -1 || ap == 0  ){
                MSG_ERROR("singular matrix\n");
                return MESS_ERROR_SINGULAR;
            }
            if (iperm[col] < 0 && cabs(xc[col]) >= cabs(pivotc)) ipivot = col;  // nur pivotisieren wenn es auch sinnmacht, sonst bleibt die diagonale
            // ipivot = col;
            pivotc =xc[ipivot]; // das pivot ohne betrag holen
            U->values_cpx[unz] = pivotc;
            U->rowptr[unz++] = k;   // und rein damit ins U
            iperm[ipivot] = k;
            L->values_cpx[lnz] = 1;
            L->rowptr[lnz++]= ipivot;
            /*-----------------------------------------------------------------------------
             *  L(k+1:n,k) = x(k+1:n)/pivot
             *-----------------------------------------------------------------------------*/
            for ( p = top; p < n ; p++) {
                i = chi[p];
                if (iperm[i] < 0){
                    L->rowptr[lnz] = i;
                    L->values_cpx[lnz++] = xc[i]/pivotc;
                }
                xc[i] = 0;
            }
        }

        L->colptr[n] =lnz;
        U->colptr[n] =unz;
        ret = mess_matrix_realloc(L, L->rows, L->cols, lnz);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);
        ret = mess_matrix_realloc(U, U->rows, U->cols, unz);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_realloc);

        for ( p=0; p<lnz;p++ ) { L->rowptr[p] = iperm[L->rowptr[p]];}
        mess_free(xc);
        mess_free(chi);

    }


    solver->pinv = iperm;
    sol->rows = A->rows;
    sol->cols = A->cols;
    sol->clear = newlu_clear;
    sol->data = (void*) solver;
    sol->solve = newlu_solve;
    sol->solvet = newlu_solvet;
    sol->solveh = newlu_solveh;
    sol->solvem = newlu_solvem;
    sol->solvemt = newlu_solvemt;
    sol->solvemh = newlu_solvemh;
    sol->inverse = newlu_inverse;
    sol->getL = newlu_getL;
    sol->getU = newlu_getU;
    sol->getpermp = newlu_getpermp;
    sol->getpermq = newlu_getpermq;
    sol->data_type = solver->data_type;
    sol->rows = solver->rows =  A->rows;
    sol->cols = solver->cols = A->cols;
    SET_SOLVERNAME(sol->name, __func__);


    if ( conv == 0) {
        mess_matrix_clear(&work);
    }
    return 0;
}


