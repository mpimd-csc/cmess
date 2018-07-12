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
 * @file lib/lrcf_adi/options.c
 * @brief LRCF-ADI  options.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Initialize LRCF-ADI options with default values.
 * @param[in,out] opt pointer to a mess_options object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_init function initializes a \ref mess_options object and
 * sets the options to the default values.
 *
 * \ref mess_options
 */
int mess_options_init(mess_options *opt) {
    MSG_FNAME(__func__);
    mess_options o;
    mess_try_alloc( *opt, struct mess_options_st*,sizeof (struct mess_options_st ));
    o = *opt;
    o->adi_shifts_paratype  = MESS_LRCFADI_PARA_ADAPTIVE_V;
    o->adi_maxit            = 500;
    o->adi_output           = 1;
    o->adi_rel_change_tol   = mess_eps()*100.0;
    o->adi_res2_tol         = 1e-10;
    o->adi_res2c_tol        = 1e-11;
    o->adi_shifts_p         = NULL;
    o->memory_usage         = MESS_MEMORY_MID;
    o->nm_gpStep            = 1;
    o->nm_linesearch        = 0;
    o->residual_method      = MESS_RESIDUAL_INDEFINITE;
    o->nm_maxit             = 30;
    o->nm_output            = 1;
    o->nm_res2_tol          = sqrt(mess_eps());
    o->nm_singleshifts      = 0;
    o->nm_stepfunction      = NULL;
    o->nm_stepfunction_aux  = NULL;
    o->stepfunction         = NULL;
    o->stepfunction_aux     = NULL;
    o->type                 = MESS_OP_NONE;
    o->adi_shifts_l0        = 16;
    o->adi_shifts_arp_p     = 50;
    o->adi_shifts_arp_m     = 25;
    o->proj_space           = NULL;
    o->K0                   = NULL;
    o->W                    = NULL;
    o->adi_shifts_b0        = NULL;

    return 0;
}


/**
 * @brief Copy a mess_options object to another one.
 * @param[in]  src   input source object to copy from
 * @param[out] dest     destination object to copy to
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_copy function copies a \ref mess_options object
 * to another preinitialized one.
 *
 * \ref mess_options
 * \ref mess_options_init
 *
 */
int mess_options_copy ( mess_options src, mess_options dest )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(src);
    mess_check_nullpointer(dest);

    memcpy(dest, src, sizeof(struct mess_options_st));
    dest ->adi_shifts_p = NULL;
    if (src -> adi_shifts_p != NULL){
        ret = mess_vector_init(&(dest->adi_shifts_p));                                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_alloc((dest->adi_shifts_p), src->adi_shifts_p->dim, src->adi_shifts_p->data_type);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
        ret = mess_vector_copy(src->adi_shifts_p, dest->adi_shifts_p);                                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    }

    dest->K0 = NULL;
    if ( src->K0 ) {
        ret = mess_matrix_init(&(dest->K0));            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_copy(src->K0, dest->K0);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }

    dest->proj_space = NULL;
    if ( src->proj_space ) {
        ret = mess_matrix_init(&(dest->proj_space));                        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_copy(src->proj_space, dest->proj_space);          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }

    dest->W = NULL;
    if ( src->W ) {
        ret = mess_matrix_init(&(dest->W));             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_copy(src->W, dest->W);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
    }

    return ret;
}       /* -----  end of function mess_options_copy  ----- */


/**
 * @brief Clear a @ref mess_options object.
 * @param[in,out] opt pointer to a @ref mess_options object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_clear function cleans up a @ref mess_options object.
 *
 */
int mess_options_clear(mess_options *opt){
    // MSG_FNAME(__func__);
    if ( opt == NULL) return 0;
    if ( *opt == NULL) return 0;
    if ( (*opt)->adi_shifts_p ) {
        mess_vector_clear(&((*opt)->adi_shifts_p));
    }

    if ( (*opt)->K0 ) {
        mess_matrix_clear(&((*opt)->K0));
    }

    if ( (*opt)->proj_space ) {
        mess_matrix_clear(&((*opt)->proj_space));
    }

    if ( (*opt)->W ) {
        mess_matrix_clear(&((*opt)->W));
    }

    mess_free(*opt);
    *opt = NULL;
    return 0;
}



/**
 * @brief Return the size of a @ref mess_options in bytes.
 * @param[in] opt  input @ref mess_options instance
 * @return approximate size of the @p stat in bytes.
 *
 * The @ref mess_options_memsize function determines the size of @p opt in bytes. \n
 *
 * @sa mess_matrix_memsize
 * @sa mess_vector_memsize
 * @sa mess_status_memsize
 *
 */
mess_int_t mess_options_memsize(mess_options opt){
    MSG_FNAME(__func__);
    mess_int_t size = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(opt);


    /*-----------------------------------------------------------------------------
     *  size of flat fields
     *-----------------------------------------------------------------------------*/
    size += sizeof(mess_options_st);


    /*-----------------------------------------------------------------------------
     *  size of deep fields
     *-----------------------------------------------------------------------------*/
    if(opt->adi_shifts_p){
        size += mess_vector_memsize(opt->adi_shifts_p);
    }

    if(opt->adi_shifts_b0){
        size += mess_vector_memsize(opt->adi_shifts_b0);
    }

    if(opt->proj_space){
        size += mess_matrix_memsize(opt->proj_space);
    }

    if(opt->K0){
        size += mess_matrix_memsize(opt->K0);
    }

    if(opt->W){
        size += mess_matrix_memsize(opt->W);
    }

    return size;
}



/**
 * @brief Display a mess_options object on standard output.
 * @param[in] opt   input object to display
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_print function displays the basis information of
 * a @ref mess_options object on the screen.
 *
 */
int mess_options_print(mess_options opt){
    MSG_FNAME(__func__);
    mess_check_nullpointer(opt);

    MSG_PRINT("type:                \t %s  \n", mess_operation_t_str(opt->type));
    MSG_PRINT("adi_maxit:           \t " MESS_PRINTF_INT " \n", opt->adi_maxit );
    MSG_PRINT("adi_output:          \t " MESS_PRINTF_INT " \n", opt->adi_output);
    MSG_PRINT("res2_tol:            \t %lg \n", opt->adi_res2_tol );
    MSG_PRINT("res2c_tol:           \t %lg \n", opt->adi_res2c_tol );
    MSG_PRINT("rel_change_tol:      \t %lg \n", opt->adi_rel_change_tol  );
    MSG_PRINT("adi_shifts_paratype: \t %s  \n", mess_parameter_t_str(opt->adi_shifts_paratype));
    MSG_PRINT("adi_shifts_arp_p:    \t " MESS_PRINTF_INT " \n", opt->adi_shifts_arp_p );
    MSG_PRINT("adi_shifts_arp_m:    \t " MESS_PRINTF_INT " \n", opt->adi_shifts_arp_m );
    MSG_PRINT("adi_shifts_l0:       \t " MESS_PRINTF_INT " \n", opt->adi_shifts_l0 );
    MSG_PRINT("nm_maxit:            \t " MESS_PRINTF_INT " \n", opt->nm_maxit );
    MSG_PRINT("nm_res2_tol:         \t %lg \n", opt->nm_res2_tol);
    MSG_PRINT("nm_singleshifts:     \t %u  \n", opt->nm_singleshifts);
    MSG_PRINT("nm_gpStep:           \t " MESS_PRINTF_INT " \n", opt->nm_gpStep);
    MSG_PRINT("nm_output:           \t " MESS_PRINTF_INT " \n", opt->nm_output);
    MSG_PRINT("nm_linesearch:       \t " MESS_PRINTF_INT " \n", opt->nm_linesearch );
    MSG_PRINT("memory_usage:        \t %s  \n", mess_memusage_t_str(opt->memory_usage) );
    MSG_PRINT("residual_method:     \t %s  \n", mess_residual_t_str(opt->residual_method));

    if(opt->K0){
        MSG_PRINT("K0:                  \t 0x%.12x,  Size: "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT", %s, %s, %s\n",
        opt->K0, opt->K0->rows, opt->K0->cols, mess_datatype_t_str(opt->K0->data_type), mess_storage_t_str(opt->K0->store_type), mess_symmetry_t_str(opt->K0->symmetry));
    }else{
        MSG_PRINT("K0:                  \t 0x%.12x \n",opt->K0);
    }

    if(opt->W){
        MSG_PRINT("W:                   \t 0x%.12x,  Size: "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT", %s, %s, %s\n",
        opt->W, opt->W->rows, opt->W->cols, mess_datatype_t_str(opt->W->data_type), mess_storage_t_str(opt->W->store_type), mess_symmetry_t_str(opt->W->symmetry));
    }else{
        MSG_PRINT("W:                   \t 0x%.12x \n",opt->W);
    }

    if(opt->proj_space){
        MSG_PRINT("proj_space:          \t 0x%.12x,  Size: "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT", %s, %s, %s\n",
        opt->proj_space, opt->proj_space->rows, opt->proj_space->cols, mess_datatype_t_str(opt->proj_space->data_type), mess_storage_t_str(opt->proj_space->store_type),
        mess_symmetry_t_str(opt->proj_space->symmetry));
    }else{
        MSG_PRINT("proj_space:          \t 0x%.12x \n",opt->proj_space);
    }


    return 0;
}





/**
 * @brief Copy the @ref mess_matrix @p proj_space
 * @param[in] opt           input @ref mess_options instance
 * @param[in] proj_space    input @ref mess_matrix instance
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_set_proj_space set the @ref mess_matrix instance @p proj_space
 * to the @c proj_space field of @p opt. Internally it is done by a copy.
 *
 */
int mess_options_set_proj_space(mess_options opt, mess_matrix proj_space){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(opt);

    if(proj_space != NULL){

        if(opt->proj_space != NULL){
            ret = mess_matrix_reset(opt->proj_space);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
            ret = mess_matrix_copy(proj_space,opt->proj_space);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            return ret;
        }else{
            ret = mess_matrix_init(&(opt->proj_space));             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_copy(proj_space,opt->proj_space);     FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            return ret;
        }

    }else{

        if(opt->proj_space != NULL){
            ret = mess_matrix_clear(&(opt->proj_space));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
            opt->proj_space = NULL;
            return ret;
        }else{
            opt->proj_space = NULL;
            return ret;
        }

    }
    MSG_ERROR("Impossible if-else branching\n");
    return MESS_ERROR_NOSUPPORT;

}

/**
 * @brief Delete the field @p proj_space from @ref mess_options instance.
 * @param[in] opt           input @ref mess_options instance
 * @return zero on success or a non-zero error value otherwise
 *
 *
 */
int mess_options_unset_proj_space(mess_options opt){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(opt);

    if(opt->proj_space != NULL){
        ret = mess_matrix_clear(&(opt->proj_space));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        opt->proj_space = NULL;
    }
    return 0;
}


/**
 * @brief Copy the @ref mess_matrix @p K0
 * @param[in] opt   input @ref mess_options instance
 * @param[in] K0    input @ref mess_matrix instance
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_set_K0 set the @ref mess_matrix instance @p K0
 * to the @c K0 field of @p opt. Internally it is done by a copy.
 *
 */
int mess_options_set_K0(mess_options opt, mess_matrix K0){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(opt);

    if(K0 != NULL){

        if(opt->K0 != NULL){
            ret = mess_matrix_reset(opt->K0);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
            ret = mess_matrix_copy(K0,opt->K0);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            return ret;
        }else{
            ret = mess_matrix_init(&(opt->K0));             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_copy(K0,opt->K0);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            return ret;
        }

    }else{

        if(opt->K0 != NULL){
            ret = mess_matrix_clear(&(opt->K0));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
            opt->K0 = NULL;
            return ret;
        }else{
            opt->K0 = NULL;
            return ret;
        }

    }
    MSG_ERROR("Impossible if-else branching\n");
    return MESS_ERROR_NOSUPPORT;

}

/**
 * @brief Delete the field @p K0 from @ref mess_options instance.
 * @param[in] opt           input @ref mess_options instance
 * @return zero on success or a non-zero error value otherwise
 *
 *
 */
int mess_options_unset_K0(mess_options opt){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(opt);

    if(opt->K0 != NULL){
        ret = mess_matrix_clear(&(opt->K0));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        opt->K0 = NULL;
    }
    return 0;
}


/**
 * @brief Copy the @ref mess_matrix @p W
 * @param[in] opt   input @ref mess_options instance
 * @param[in] W    input @ref mess_matrix instance
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_options_set_W set the @ref mess_matrix instance @p W
 * to the @c W field of @p opt. Internally it is done by a copy.
 *
 */
int mess_options_set_W(mess_options opt, mess_matrix W){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(opt);

    if(W != NULL){

        if(opt->W != NULL){
            ret = mess_matrix_reset(opt->W);               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_reset);
            ret = mess_matrix_copy(W,opt->W);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            return ret;
        }else{
            ret = mess_matrix_init(&(opt->W));             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_init);
            ret = mess_matrix_copy(W,opt->W);             FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_copy);
            return ret;
        }

    }else{

        if(opt->W != NULL){
            ret = mess_matrix_clear(&(opt->W));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
            opt->W = NULL;
            return ret;
        }else{
            opt->W = NULL;
            return ret;
        }

    }
    MSG_ERROR("Impossible if-else branching\n");
    return MESS_ERROR_NOSUPPORT;

}

/**
 * @brief Delete the field @p W from @ref mess_options instance.
 * @param[in] opt           input @ref mess_options instance
 * @return zero on success or a non-zero error value otherwise
 *
 *
 */
int mess_options_unset_W(mess_options opt){
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(opt);

    if(opt->W != NULL){
        ret = mess_matrix_clear(&(opt->W));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        opt->W = NULL;
    }
    return 0;
}
















