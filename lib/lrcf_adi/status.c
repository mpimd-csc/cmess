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
 * @file lib/lrcf_adi/status.c
 * @brief LRCF-ADI status.
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
 * @brief Initialize a mess_status object.
 * @param[in,out] status pointer to the object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_status_init function initializes a mess_status object with empty values. \n
 * In case of erros the following error codes can be returned:
 * \li @ref MESS_ERROR_MEMORY no memory left
 * \li all errors from mess_vector_init
 *
 * \sa mess_status
 * \sa mess_options_init
 */
int mess_status_init( mess_status *status){
    MSG_FNAME(__func__);
    int ret;
    mess_try_alloc( *status , struct mess_status_st*, sizeof (struct mess_status_st ));

    MESS_INIT_VECTORS(&((*status)->res2_norms), &((*status)->rel_changes))
    ret = mess_vector_alloc(((*status)->res2_norms),1, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_alloc(((*status)->rel_changes),1, MESS_REAL); FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_init);

    (*status)->res2_norm = 0;
    (*status)->it = 0;
    (*status)->res2_change = 0;
    (*status)->res2_0 = 0;
    (*status)->rel_change = 0;
    (*status)->stop_res2 = 0;
    (*status)->stop_res2c = 0;
    (*status)->stop_rel = 0;
    (*status)->stop_user = 0;
    (*status)->time_all = 0;
    (*status)->time_adi = 0;
    (*status)->unstable = 0 ;
    (*status)->n_internal_status = 0;
    (*status)->internal_status = NULL;
    return 0;
}

/**
 * @brief Clear a mess_status object.
 * @param[in,out] status pointer to the object
 * @return always zero
 *
 *
 * The @ref mess_status_clear cleans up a @ref mess_status object.
 *
 * \sa mess_status_init
 */
int mess_status_clear( mess_status *status){
    // MSG_FNAME(__func__);
    if ( status == NULL) return 0;
    if (*status == NULL) return 0;
    mess_vector_clear(&((*status)->res2_norms));
    mess_vector_clear(&((*status)->rel_changes));
    if ( (*status)->n_internal_status > 0 ) {
        mess_int_t i = 0;
        for ( i = 0 ; i < (*status)->n_internal_status; i++) {
            mess_status_clear(&((*status)->internal_status[i]));
        }
        mess_free((*status)->internal_status);
    }

    mess_free(*status);
    *status = NULL;
    return 0;
}




/**
 * @brief Copy a @ref mess_status object to another one.
 * @param[in]  src      input source object to copy from
 * @param[out] dest     output  object to copy to
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_status_copy function copies a \ref mess_status object
 * to another preinitialized one.
 *
 * \ref mess_status
 * \ref mess_status_init
 *
 */
int mess_status_copy(mess_status src, mess_status dest){
    MSG_FNAME(__func__);
    int ret = 0;
    int i = 0;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(src);
    mess_check_nullpointer(dest);

    /*-----------------------------------------------------------------------------
     *  make a flat copy
     *-----------------------------------------------------------------------------*/
    memcpy(dest, src, sizeof(mess_status_st));

    /*-----------------------------------------------------------------------------
     *  deep copy
     *-----------------------------------------------------------------------------*/
    if(src->res2_norms){
        dest->res2_norms = NULL;
        ret = mess_vector_init(&(dest->res2_norms));
        ret = mess_vector_alloc(dest->res2_norms, src->res2_norms->dim, src->res2_norms->data_type);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_copy(src->res2_norms,dest->res2_norms);                                               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    }

    if(src->rel_changes){
        dest->rel_changes = NULL;
        ret = mess_vector_init(&(dest->rel_changes));
        ret = mess_vector_alloc(dest->rel_changes, src->rel_changes->dim, src->rel_changes->data_type);         FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_alloc);
        ret = mess_vector_copy(src->rel_changes,dest->rel_changes);                                             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
    }

    if(src->n_internal_status){
        dest->internal_status = NULL;
        mess_try_alloc(dest->internal_status, mess_status*, (sizeof(mess_status_st)*src->n_internal_status));
    }


    /*-----------------------------------------------------------------------------
     *  deeep copy recursively
     *-----------------------------------------------------------------------------*/
    for(i=0;i<src->n_internal_status;++i){
        ret = mess_status_copy(src->internal_status[i],dest->internal_status[i]);                               FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_status_copy);
    }

    return ret;
}/* -----  end of function mess_status_copy  ----- */



/**
 * @brief Return the size of a @ref mess_status in bytes.
 * @param[in] stat  input @ref mess_status instance
 * @return approximate size of the @p stat in bytes.
 *
 * The @ref mess_status_memsize function determines the size of @p stat in bytes. \n
 *
 * @sa mess_matrix_memsize
 * @sa mess_vector_memsize
 *
 */
mess_int_t mess_status_memsize(mess_status stat){
    MSG_FNAME(__func__);
    mess_int_t size = 0;
    int i = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(stat);


    /*-----------------------------------------------------------------------------
     *  size of flat fields
     *-----------------------------------------------------------------------------*/
    size += sizeof(mess_status_st);


    /*-----------------------------------------------------------------------------
     *  size of deep fields
     *-----------------------------------------------------------------------------*/
    if(stat->res2_norms){
        size += mess_vector_memsize(stat->res2_norms);
    }

    if(stat->rel_changes){
        size += mess_vector_memsize(stat->rel_changes);
    }


    /*-----------------------------------------------------------------------------
     *  recursive to other status instances
     *-----------------------------------------------------------------------------*/
    for(i=0;i<stat->n_internal_status;++i){
        size += mess_status_memsize(stat->internal_status[i]);
    }

    return size;
}



/**
 * @brief Print a mess_status object to the standard output.
 * @param[in] stat   input object to print
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_status_print function prints basic infomation (except of the residual and
 * rel_change vectors) on stdout.
 *
 * \sa mess_status
 * \sa mess_status_printfull
 * \sa mess_status_write
 */
int mess_status_print(mess_status stat){
    MSG_FNAME(__func__);
    mess_check_nullpointer(stat);


    /*-----------------------------------------------------------------------------
     *  print fields
     *-----------------------------------------------------------------------------*/
    MSG_PRINT("iterations:      \t " MESS_PRINTF_INT "\n", stat->it);
    MSG_PRINT("res2_norm:       \t %.9e\n", stat->res2_norm);
    MSG_PRINT("res2_change:     \t %.9e\n", stat->res2_change);
    MSG_PRINT("rel_change:      \t %.9e\n", stat->rel_change);
    MSG_PRINT("res2_0:          \t %.9e\n", stat->res2_0);
    MSG_PRINT("stop_res2:       \t %u\n", stat->stop_res2);
    MSG_PRINT("stop_res2c:      \t %u\n", stat->stop_res2c);
    MSG_PRINT("stop_rel:        \t %u\n", stat->stop_rel);
    MSG_PRINT("stop_user:       \t %u\n", stat->stop_user);
    MSG_PRINT("time - all:      \t %lg\n", stat->time_all);
    MSG_PRINT("time - adi:      \t %lg\n", stat->time_adi);
    MSG_PRINT("unstable         \t %u \n", stat->unstable);

    if(stat->res2_norms){
        MSG_PRINT("res2_norms:      \t 0x%.12x,  Dim: "MESS_PRINTF_INT", %s\n",stat->res2_norms, stat->res2_norms->dim,  mess_datatype_t_str(stat->res2_norms->data_type));
    }else{
        MSG_PRINT("res2_norms:      \t 0x%.12x \n", stat->res2_norms);
    }

    if(stat->rel_changes){
        MSG_PRINT("rel_changes:      \t 0x%.12x,  Dim: "MESS_PRINTF_INT", %s\n",stat->rel_changes, stat->rel_changes->dim,  mess_datatype_t_str(stat->rel_changes->data_type));
    }else{
        MSG_PRINT("rel_changes:      \t 0x%.12x \n", stat->rel_changes);
    }


    return 0;
}



/**
 * @brief Print a mess_status object to the standard output including embbed objects.
 * @param[in] stat   input object to print
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_status_printfull function displays all information contained in the given
 * \ref mess_status object to the standard output. This includes embedded status objects.
 *
 * \sa mess_status
 * \sa mess_status_print
 * \sa mess_status_write
 */
int mess_status_printfull(mess_status stat){
    MSG_FNAME(__func__);
    if ( stat == NULL)   {
        MSG_ERROR("stat points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    MSG_PRINT("iterations:    \t " MESS_PRINTF_INT "\n", stat->it);
    MSG_PRINT("res2_norm:     \t %.9e\n", stat->res2_norm);
    MSG_PRINT("res2_change:   \t %.9e\n", stat->res2_change);
    MSG_PRINT("rel_change:    \t %.9e\n", stat->rel_change);
    MSG_PRINT("res2_0:      \t %.9e\n", stat->res2_0);
    MSG_PRINT("stop_res2:   \t %u\n", stat->stop_res2);
    MSG_PRINT("stop_res2c:  \t %u\n", stat->stop_res2c);
    MSG_PRINT("stop_rel:    \t %u\n", stat->stop_rel);
    MSG_PRINT("stop_userl:  \t %u\n", stat->stop_user);
    MSG_PRINT("time - all:    \t %lg\n", stat->time_all);
    MSG_PRINT("time - adi:    \t %lg\n", stat->time_adi);
    MSG_PRINT("unstable       \t %u \n", stat->unstable);

    MSG_PRINT("\n");
    MSG_PRINT("change of res2-norms:\n");
    mess_int_t i = 0;
    for(i=0;i<stat->res2_norms->dim;i++){
        MSG_PRINT("[ " MESS_PRINTF_INT " ] = %.9e\n",i, stat->res2_norms->values[i]);
    }

    MSG_PRINT("\n");
    MSG_PRINT("relative change of the factor Z:\n");
    for(i=0;i<stat->rel_changes->dim;i++){
        MSG_PRINT("[ " MESS_PRINTF_INT " ] = %.9e\n", i, stat->rel_changes->values[i]);
    }
    MSG_PRINT("\n");

    if ( stat->n_internal_status > 0 ) {
        for ( i = 0; i < stat->n_internal_status; i++) {
            MSG_PRINT("\nInternal Iteration " MESS_PRINTF_INT ":\n", i);
            mess_status_print(stat->internal_status[i]);
        }
    }
    return 0;
}

/**
 * @brief Save the information of a mess_status object to a text file.
 * @param[in] stat   input object to save
 * @param[in] filename    input filename of the text file
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_status_write function saves all information contained in a
 * \ref mess_status object to a text file \f$ filename \f$.
 *
 * \sa mess_status
 * \sa mess_status_print
 * \sa mess_status_printfull
 */
int mess_status_write(mess_status stat, const char *filename){
    MSG_FNAME(__func__);
    FILE *fp;
    mess_int_t i;
    int err;
    if ( stat == NULL)   {
        MSG_ERROR("stat points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if (filename == NULL ){
        MSG_ERROR("filename points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }

    if ( (fp = fopen(filename,"w")) == NULL){
        err = errno;
        MSG_ERROR("error open file: %s\n", strerror(err));
        return MESS_ERROR_FILEIO;
    }

    fprintf(fp,"iterations:     \t " MESS_PRINTF_INT "\n", stat->it);
    fprintf(fp,"res2_norm:      \t %.8e\n", stat->res2_norm);
    fprintf(fp,"res2_change:    \t %.8e\n", stat->res2_change);
    fprintf(fp,"rel_change:     \t %.8e\n", stat->rel_change);
    fprintf(fp,"res2_0:         \t %.8e\n", stat->res2_0);
    fprintf(fp,"stop_res2:      \t %hu\n", stat->stop_res2);
    fprintf(fp,"stop_res2c:     \t %hu\n", stat->stop_res2c);
    fprintf(fp,"stop_rel:       \t %hu\n", stat->stop_rel);
    fprintf(fp,"time-all:       \t %lg\n", stat->time_all);
    fprintf(fp,"time - adi:     \t %lg\n", stat->time_adi);
    fprintf(fp,"unstable:       \t %hu\n", stat->unstable);

    fprintf(fp,"\n");
    fprintf(fp,"change of res2-norms:\n");
    for(i=0;i<stat->res2_norms->dim;i++){
        fprintf(fp, "[ " MESS_PRINTF_INT " ] = %.8e\n",i, stat->res2_norms->values[i]);
    }

    fprintf(fp,"\n");
    fprintf(fp,"relative change of the factor Z:\n");
    for(i=0;i<stat->rel_changes->dim;i++){
        fprintf(fp, "[ " MESS_PRINTF_INT " ] = %.8e\n", i, stat->rel_changes->values[i]);
    }
    fprintf(fp,"\n");
    fclose(fp);
    return 0;
}

