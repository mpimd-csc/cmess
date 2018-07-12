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
 * @file lib/misc/freelist.c
 * @brief Memory manager.
 * @author @koehlerm
 * @author @mbehr
 */

#include <string.h>
#include "mess/mess.h"
#include "mess/error_macro.h"

/**
 * @brief Initialized a @ref mess_freelist instance.
 * @param[in,out] list allocated @ref mess_freelist instance.
 *
 *  Initialize a @ref mess_freelist instance.
 */
int mess_freelist_init(mess_freelist *list)
{
    MSG_FNAME(__func__);
    int ret=0;
    mess_try_alloc(*list, mess_freelist, sizeof(struct mess_freelist_st));
    memset(*list,0,sizeof(struct mess_freelist_st));
    return ret;
}

//macro for mess_freelist_add* functions
#define MESS_FREELIST_ADD_GENERIC(FNAME, ADDTYPE, ADDVAR, LISTVARCOUNTER, LISTVARPOINTER)           \
    int FNAME(mess_freelist list, ADDTYPE ADDVAR){                                                  \
        MSG_FNAME(__func__);                                                                        \
        int ret=0;                                                                                  \
        size_t newlen = list->LISTVARCOUNTER+1;                                                     \
        mess_try_realloc(list->LISTVARPOINTER , ADDTYPE*, sizeof(ADDTYPE)* newlen);                 \
        list->LISTVARPOINTER[list-> LISTVARCOUNTER] = ADDVAR;                                       \
        list->LISTVARCOUNTER = newlen;                                                              \
        return ret;                                                                                 \
    }


/**
 * @brief Insert an arbitrary pointer in the memory manager.
 * @param[in,out] list   List containing the allocated objects
 * @param[in]     ptr    pointer to add to the @p list
 *
 * The @ref mess_freelist_add_ptr function inserts an arbitrary pointer in the list
 * of allocated objects.
 */
MESS_FREELIST_ADD_GENERIC(mess_freelist_add_ptr, void*, ptr, n_ptrs, ptrs);


/**
 * @brief Insert a  @ref mess_vector in the memory manager.
 * @param[in,out] list   list containing the allocated objects
 * @param[in]     vec    vector to add to the @p list
 *
 * The @ref mess_freelist_add_mess_vector function inserts @p vec in @p list.
 *
 */
MESS_FREELIST_ADD_GENERIC(mess_freelist_add_mess_vector, mess_vector, vec, n_vectors, vectors);


/**
 * @brief Insert a @ref mess_matrix in the memory manager.
 * @param[in,out] list   list containing the allocated objects
 * @param[in]      mat   matrix to add to @p list
 *
 * The @ref mess_freelist_add_mess_matrix function inserts @p mat in @p list.
 *
 */
MESS_FREELIST_ADD_GENERIC(mess_freelist_add_mess_matrix, mess_matrix, mat, n_matrices, matrices);


/**
 * @brief Insert a @ref mess_equation in the memory manager.
 * @param[in,out] list   list containing the allocated objects
 * @param[in]      eqn   equation to add to @p list
 *
 * The @ref mess_freelist_add_mess_equation function inserts @p eqn in @p list.
 *
 */
MESS_FREELIST_ADD_GENERIC(mess_freelist_add_mess_equation, mess_equation, eqn, n_eqns, eqns);


/**
 * @brief Insert a @ref mess_options structure in the memory manager.
 * @param[in,out] list   List containing the allocated objects
 * @param[in]     opt    options to add to the list
 *
 * The @ref mess_freelist_add_mess_options function inserts a @ref mess_status instance in the list
 * of allocated objects.
 */
MESS_FREELIST_ADD_GENERIC(mess_freelist_add_mess_options, mess_options, opt, n_opts, opts);


/**
 * @brief Insert a @ref mess_status structure in the memory manager.
 * @param[in,out] list   List containing the allocated objects
 * @param[in]     stat   status to add to the @p list
 *
 * The @ref mess_freelist_add_mess_status function inserts @p stat in @p list.
 */
MESS_FREELIST_ADD_GENERIC(mess_freelist_add_mess_status, mess_status, stat, n_stats, stats);




/**
 * @brief Print a @ref mess_freelist instance.
 * @param[in] list   the @ref mess_freelist instance.s
 * @return 0.
 *
 * Print content of @ref mess_freelist instance @p list.
 */
int mess_freelist_print(mess_freelist list){
    MSG_FNAME(__func__);
    size_t j = 0;

    if(list==NULL){
        MSG_WARN("mes_freelist list points to NULL.\n");
        return 0;
    }

    /*-----------------------------------------------------------------------------
     *  print number of contained objects
     *-----------------------------------------------------------------------------*/

    MSG_PRINT("Instances:\n");
    MSG_PRINT("\tpointer:       %d\n",list->n_ptrs);
    MSG_PRINT("\tmess_vector:   %d\n",list->n_vectors);
    MSG_PRINT("\tmess_matrix:   %d\n",list->n_matrices);
    MSG_PRINT("\tmess_equation: %d\n",list->n_eqns);
    MSG_PRINT("\tmess_options:  %d\n",list->n_opts);
    MSG_PRINT("\tmess_status:   %d\n",list->n_stats);

    /*-----------------------------------------------------------------------------
     *  print pointers
     *-----------------------------------------------------------------------------*/
    if(list->n_ptrs){
        MSG_PRINT("\npointers:\n");
        MSG_PRINT("\tIndex |     Memory Address   \n");
        MSG_PRINT("\t------+----------------------\n");
        for(j=0;j<list->n_ptrs;++j){
            MSG_PRINT("\t%5u |     0x%.12x   \n",j,list->ptrs[j]);
        }
    }

    /*-----------------------------------------------------------------------------
     *  print vectors
     *-----------------------------------------------------------------------------*/
    if(list->n_vectors){
        MSG_PRINT("\nmess_vector objects:\n");
        MSG_PRINT("\tIndex |     Memory Address   | Information\n");
        MSG_PRINT("\t------+----------------------+-------------\n");
        for(j=0;j<list->n_vectors;++j){
            if(list->vectors[j]){
                MSG_PRINT("\t%5u |     0x%.12x   | Dim: "MESS_PRINTF_INT", %s\n",
                        j,
                        list->vectors[j],
                        list->vectors[j]->dim,
                        mess_datatype_t_str(list->vectors[j]->data_type)
                        );
            }else{
                MSG_PRINT("\t%5u |     0x%.12x   | ---\n",j,list->vectors[j]);
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  print matrices
     *-----------------------------------------------------------------------------*/
    if(list->n_matrices){
        MSG_PRINT("\nmess_matrix objects:\n");
        MSG_PRINT("\tIndex |     Memory Address   | Information\n");
        MSG_PRINT("\t------+----------------------+-------------\n");
        for(j=0;j<list->n_matrices;++j){
            if(list->matrices[j]){
                MSG_PRINT("\t%5u |     0x%.12x   | Size: "MESS_PRINTF_INT"-x-"MESS_PRINTF_INT", %s, %s, %s\n",
                        j,
                        list->matrices[j],
                        list->matrices[j]->rows,
                        list->matrices[j]->cols,
                        mess_datatype_t_str(list->matrices[j]->data_type),
                        mess_storage_t_str(list->matrices[j]->store_type),
                        mess_symmetry_t_str(list->matrices[j]->symmetry)
                        );
            }else{
                MSG_PRINT("\t%5u |     0x%.12x   | ---\n",j,list->matrices[j]);
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  print equations
     *-----------------------------------------------------------------------------*/
    if(list->n_eqns){
        MSG_PRINT("\nmess_equation objects:\n");
        MSG_PRINT("\tIndex |     Memory Address   | Information\n");
        MSG_PRINT("\t------+----------------------+-------------\n");
        for(j=0;j<list->n_eqns;++j){
            if(list->eqns[j]){
                MSG_PRINT("\t%5u |     0x%.12x   | %s\n",j,list->eqns[j],mess_equation_t_str(list->eqns[j]->eqn_type));
            }else{
                MSG_PRINT("\t%5u |     0x%.12x   | ---\n",j,list->eqns[j]);
            }
        }
    }

    /*-----------------------------------------------------------------------------
     *  print opts
     *-----------------------------------------------------------------------------*/
    if(list->n_opts){
        MSG_PRINT("\nmess_options objects:\n");
        MSG_PRINT("\tIndex |     Memory Address   \n");
        MSG_PRINT("\t------+----------------------\n");
        for(j=0;j<list->n_opts;++j){
            MSG_PRINT("\t%5u |     0x%.12x   \n",j,list->opts[j]);
        }
    }

    /*-----------------------------------------------------------------------------
     *  print status
     *-----------------------------------------------------------------------------*/
    if(list->n_stats){
        MSG_PRINT("\nmess_status objects:\n");
        MSG_PRINT("\tIndex |     Memory Address   \n");
        MSG_PRINT("\t------+----------------------\n");
        for(j=0;j<list->n_stats;++j){
            MSG_PRINT("\t%5u |     0x%.12x   \n",j,list->stats[j]);
        }
    }



    return 0;
}


/**
 * @brief Free all objects stored in the memory manager.
 * @param[in] list  input list of objects to free
 *
 * The @ref mess_freelist_clear function removes all objetcs in the list
 * from memory.
 */
int mess_freelist_clear(mess_freelist *list)
{

    MSG_FNAME(__func__);
    int ret=0;
    size_t i;

    /*-----------------------------------------------------------------------------
     *  clear pointers
     *-----------------------------------------------------------------------------*/

    for (i = 0; i < (*list)->n_ptrs; i++) {
        if ( (*list)->ptrs[i] ) {
            mess_free((*list)->ptrs[i]);
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear vectors
     *-----------------------------------------------------------------------------*/
    for (i = 0; i < (*list)->n_vectors; i++) {
        if ((*list)->vectors[i]) {
            ret = mess_vector_clear(&((*list)->vectors[i]));           FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_vector_clear);
        }
    }


    /*-----------------------------------------------------------------------------
     *  clear matrices
     *-----------------------------------------------------------------------------*/
    for (i = 0; i < (*list)->n_matrices; i++) {
        if ((*list)->matrices[i]){
            ret = mess_matrix_clear(&((*list)->matrices[i]));          FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_matrix_clear);
        }
    }

    /*-----------------------------------------------------------------------------
     *  clear equations
     *-----------------------------------------------------------------------------*/

    for (i = 0; i < (*list)->n_eqns; i++) {
        if ( (*list)->eqns[i] ) {
            ret = mess_equation_clear(&((*list)->eqns[i]));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_equation_clear);
        }
    }


    /*-----------------------------------------------------------------------------
     *  clear options
     *-----------------------------------------------------------------------------*/

    for (i = 0; i < (*list)->n_opts; i++) {
        if ( (*list)->opts[i] ) {
            ret = mess_options_clear(&((*list)->opts[i]));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_options_clear);
        }
    }


    /*-----------------------------------------------------------------------------
     *  clear status
     *-----------------------------------------------------------------------------*/
    for (i = 0; i < (*list)->n_stats; i++) {
        if ( (*list)->stats[i] ) {
            ret = mess_status_clear(&((*list)->stats[i]));            FUNCTION_FAILURE_HANDLE(ret,(ret!=0),mess_status_clear);
        }
    }



    mess_free((*list)->ptrs);       (*list)->ptrs = NULL;       (*list)->n_ptrs=0;
    mess_free((*list)->vectors);    (*list)->vectors = NULL;    (*list)->n_vectors=0;
    mess_free((*list)->matrices);   (*list)->matrices = NULL;   (*list)->n_matrices=0;
    mess_free((*list)->eqns);       (*list)->eqns = NULL;       (*list)->n_eqns=0;
    mess_free((*list)->opts);       (*list)->opts = NULL;       (*list)->n_opts=0;
    mess_free((*list)->stats);      (*list)->stats = NULL;      (*list)->n_stats=0;

    mess_free(*list);
    *list=NULL;

    return ret;
}



