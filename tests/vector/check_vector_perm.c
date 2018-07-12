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
 * @file tests/vector/check_vector_perm.c
 * @brief Check the permutation function for vectors.
 *
 * @author  @mbehr
 * @test
 * This function checks
 * * @ref mess_vector_perm
 * * @ref mess_vector_perm_inplace
 * * @ref mess_vector_perm_combine
 * * @ref mess_vector_perm_split
 * * @ref mess_vector_iperm
 * * @ref mess_vector_iperm_inplace
 * * @ref mess_vector_iperm_combine
 * * @ref mess_vector_iperm_split
 * * @ref
 * for @ref MESS_REAL and  @ref MESS_COMPLEX format.
 * @}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "mess/mess.h"
#include "../call_macro.h"



mess_int_t * permute(mess_int_t dim) {
    srand(time(NULL));
    mess_int_t *p;
    p = (mess_int_t*) malloc(sizeof(mess_int_t)*dim);
    mess_int_t k;
    for (k = 0; k < dim; k++)
        p[k] = k;
    for (k = dim-1; k > 0; k--) {
        mess_int_t j = rand() % (k+1);
        mess_int_t temp = p[j];
        p[j] = p[k];
        p[k] = temp;
    }
    return p;
}

/**
 * @internal
 * @brief Checks correctness of @ref mess_vector_perm and @ref mess_vector_iperm.
 * @attention Internal use only.
 * */
int check_perm(mess_vector v, double eps){
    double diff;
    int err=0;
    mess_int_t *perm = permute(v->dim);
    mess_vector vcopy,vtmp;
    mess_vector_init(&vtmp,v->dim,v->data_type);
    mess_vector_init(&vcopy,v->dim, v->data_type);
    mess_vector_copy(v,vcopy);

    mess_vector_perm(vcopy,perm,vtmp);
    mess_vector_iperm(vtmp,perm,vcopy);
    mess_vector_diffnorm(vcopy,v,&diff);
    if(diff>eps){
        printf("mess_vector_(i)perm failed with:\n");
        printf("eps  = %e\n",eps);
        printf("diff = %e\n",diff);
        mess_vector_printinfo(v);
        err++;
    }

    mess_vector_clear(&vcopy);
    mess_vector_clear(&vtmp);
    free(perm);
    return err;
}

/**
 * @internal
 * @brief Checks correctness of @ref mess_vector_perm_inplace and @ref mess_vector_iperm_inplace.
 * @attention Internal use only.
 * */
int check_perm_inplace(mess_vector v, double eps){
    double diff;
    int err=0;
    mess_int_t *perm = permute(v->dim);
    mess_vector vcopy,vtmp;
    mess_vector_init(&vtmp,v->dim,v->data_type);
    mess_vector_init(&vcopy,v->dim, v->data_type);
    mess_vector_copy(v,vcopy);

    mess_vector_perm_inplace(vcopy,perm);
    mess_vector_iperm_inplace(vcopy,perm);
    mess_vector_diffnorm(vcopy,v,&diff);
    if(diff>eps){
        printf("mess_vector_(i)perm_inplace failed with:\n");
        printf("eps  = %e\n",eps);
        printf("diff = %e\n",diff);
        mess_vector_printinfo(v);
        err++;
    }

    mess_vector_clear(&vcopy);
    mess_vector_clear(&vtmp);
    free(perm);
    return err;
}
/**
 * @internal
 * @brief Checks the correctness of @ref mess_vector_perm_split and @ref mess_vector_iperm_split.
 *
 * @attention Internal use only.
 * */
int check_perm_split(mess_vector v, double eps){
    int err=0;
    double diff;
    mess_vector out1, out2, out3;
    mess_vector_init(&out1,v->dim,MESS_REAL);
    mess_vector_init(&out2,v->dim,MESS_REAL);
    mess_vector_init(&out3,v->dim,MESS_COMPLEX);
    mess_int_t* perm = permute(v->dim);

    mess_vector_perm_split(v,perm,out1,out2);
    mess_vector_axpyc(I,out2,out1);
    mess_vector_iperm_split(out1,perm,out2,out3);
    mess_vector_axpyc(I,out3,out2);
    mess_vector_diffnorm(out2,v,&diff);

    if(diff>eps){
        printf("mess_vector_(i)perm_split failed with:\n");
        printf("eps  = %e\n",eps);
        printf("diff = %e\n",diff);
        mess_vector_printinfo(v);
        err++;
    }

    mess_vector_clear(&out1);
    mess_vector_clear(&out2);
    mess_vector_clear(&out3);
    free(perm);
    return err;
}

/**
 * @internal
 * @brief Checks the correctness of @ref mess_vector_perm_combine and @ref mess_vector_iperm_combine.
 *
 * @attention Internal use only.
 * */
int check_perm_combine(mess_vector v, double eps){
    int err=0;
    double diff;
    mess_vector out1, out2, out3;
    mess_vector_init(&out1,v->dim,MESS_REAL);
    mess_vector_init(&out2,v->dim,MESS_REAL);
    mess_vector_init(&out3,v->dim,MESS_COMPLEX);
    mess_int_t* perm = permute(v->dim);

    mess_vector_perm_split(v,perm,out1,out2);
    mess_vector_iperm_combine(out1,out2,perm,out3);
    mess_vector_diffnorm(out3,v,&diff);

    if(diff>eps){
        printf("mess_vector_perm_combine failed with:\n");
        printf("eps  = %e\n",eps);
        printf("diff = %e\n",diff);
        mess_vector_printinfo(v);
        err++;
    }

    mess_vector_iperm_split(v,perm,out1,out2);
    mess_vector_perm_combine(out1,out2,perm,out3);
    mess_vector_diffnorm(out3,v,&diff);

    if(diff>eps){
        printf("mess_vector_iperm_combine failed with:\n");
        printf("eps  = %e\n",eps);
        printf("diff = %e\n",diff);
        mess_vector_printinfo(v);
        err++;
    }

    mess_vector_clear(&out1);
    mess_vector_clear(&out2);
    mess_vector_clear(&out3);
    free(perm);
    return err;
}

int main(){
    mess_version();
    int ret = 0;
    mess_int_t err = 0, dim = 200;
    mess_int_t* perm;
    double  eps=sqrt(mess_eps());
    mess_vector vreal1, vcpx1;

    /*-----------------------------------------------------------------------------
     *  Init vectors/matrices and create permutation array
     *-----------------------------------------------------------------------------*/
    //init
    CALL(mess_vector_init(&vreal1,dim, MESS_REAL));
    CALL(mess_vector_init(&vcpx1, dim, MESS_COMPLEX));

    //random
    CALL(mess_vector_rand(vreal1));
    CALL(mess_vector_rand(vcpx1));

    CALL(mess_vector_scale(dim*2.0,vreal1));
    CALL(mess_vector_map(vreal1,floor,NULL));

    //permutation
    perm = permute(dim);
    mess_int_t i;
    for(i=0;i<vreal1->dim;++i){
        vreal1->values[i]=i;
    }

    /*-----------------------------------------------------------------------------
     *  check mess_vector_perm
     *-----------------------------------------------------------------------------*/
    err += check_perm(vreal1,eps);
    err += check_perm(vcpx1,eps);

    /*-----------------------------------------------------------------------------
     *  check mess_vector_perm_inplace
     *-----------------------------------------------------------------------------*/
    err += check_perm_inplace(vreal1,eps);
    err += check_perm_inplace(vcpx1,eps);

    /*-----------------------------------------------------------------------------
     *   check mess_vector_perm_split
     *-----------------------------------------------------------------------------*/
    err += check_perm_split(vcpx1,eps);

    /*-----------------------------------------------------------------------------
     *   check mess_vector_perm_combine
     *-----------------------------------------------------------------------------*/
    err += check_perm_combine(vcpx1,eps);

    /*-----------------------------------------------------------------------------
     *  Clear Memory
     *-----------------------------------------------------------------------------*/
    mess_free(perm);
    CALL(mess_vector_clear(&vreal1));
    CALL(mess_vector_clear(&vcpx1));

    return err;
}


