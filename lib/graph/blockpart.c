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
 * @file lib/graph/blockpart.c
 * @brief Detect blocks inside a matrix.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mess/mess.h"
#include "mess/error_macro.h"



/**
 * @brief Depth-First-Search on a lower triangular matrix.
 * @param[in] G  input matrix representing the directed graph
 * @param[in] startnode  input start node for the DFS
 * @param[in,out] stack recursion stack to work on
 * @param[in,out] marked array with the node markers
 * @param[in]     marker  input value to mark nodes
 * @param[out] sT   top of the stack where the visited nodes are
 * @return 0 on success or a non zero error code otherwise
 *
 * The @ref lblock_dfs function performs a Depth-First-Search of a graph beginning
 * at the startnode. \n
 * The result is returned as a part of the \f$ stack \f$ array. The visited
 * nodes are contained in \f$ stack[st,...,G->rows-1] \f$.
 */
static int lblock_dfs(mess_matrix G, mess_int_t startnode, mess_int_t *stack, mess_int_t *marked, mess_int_t marker, mess_int_t *sT){
    MSG_FNAME(__func__);
    mess_int_t i, j, head, done, ii ;

    mess_check_nullpointer(G);
    mess_check_nullpointer(stack);
    mess_check_nullpointer(marked);
    mess_check_nullpointer(sT);
    mess_check_csr(G);
    mess_check_square(G);

    head = 0;
    stack[head] = startnode;
    *sT = G->rows ;
    if ( marked[startnode] == marker) {
        return 0;
    }

    while ( head >= 0 ) {
        j = stack[head];

        done = 1;
        marked[j] = marker;
        for (i = G->rowptr[j]; i < G->rowptr[j+1]; i++){
            ii = G->colptr[i];
            if ( ii == j ) continue;
            if ( marked[ii] == marker) continue;
            head ++;
            stack[head] = ii;
            done = 0;
            break;
        }
        if (done) {
            head--;
            *sT = (*sT) - 1;
            stack[*sT] = j ;
        }
    }
    return 0;
}



/**
 * @brief Inline function to swap to elements in an integer array.
 * @param[in,out]   V input/output array
 * @param[in]   i    input index one
 * @param[in]   j    input index two
 *
 * The @ref swap function swaps \f$ V[i] \leftrightarrow V[j] \f$.
 */
static inline void swap(mess_int_t *V, mess_int_t i, mess_int_t j) {
    mess_int_t x;
    x = V[i];
    V[i]= V[j];
    V[j] = x;
}


/**
 * @brief Compute a lower triangular block partition of a matrix.
 * @param[in] G               input matrix containing graph
 * @param[out] perm       output permutation
 * @param[out] nblocks    number of blocks
 * @param[out] blocks     array containing start indices of the blocks
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_graph_lblockpart function computes a permutation to transform matrix \f$ G \f$ into a lower
 * block triangular matrix. \n
 * It returns a permutation array \f$ perm \f$ which can be used to transform the  matrix \f$ G \f$ to a lower block
 * triangular one. \n
 * The \f$ nblocks \f$ parameter contains the number of blocks and
 * the \f$ blocks \f$ array contains starting nodes, that means the index of the upper left diagonal entry where
 * a block starts.
 *
 * \attention It only work for Compressed Sparse Row matrices.
 *
 */
int  mess_graph_lblockpart ( mess_matrix G, mess_int_t * perm, mess_int_t *nblocks, mess_int_t *blocks )
{
    MSG_FNAME(__func__);
    mess_int_t *marked = NULL;
    mess_int_t ci, i ;
    mess_int_t n;
    mess_int_t *V,sV;
    mess_int_t mini, minv;
    mess_int_t *dfs_stack, sT = 0;
    mess_int_t pperm = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(G);
    mess_check_square ( G);
    mess_check_nullpointer(perm);
    mess_check_nullpointer(nblocks);
    mess_check_nullpointer(blocks);
    mess_check_csr(G);

    if ( ! MESS_IS_GENERAL(G)) {
        MSG_ERROR("G must be a general matrix, not a symmetric stored one. \n");
        return MESS_ERROR_STORAGETYPE;
    }

    n = G->rows;
    mess_try_alloc(V, mess_int_t *, sizeof(mess_int_t) *(n));
    mess_try_alloc(marked, mess_int_t *, sizeof(mess_int_t) *(n));
    mess_try_alloc(dfs_stack, mess_int_t *, sizeof(mess_int_t) *(n));

    /*-----------------------------------------------------------------------------
     *  compute the output degree of all nodes.
     *-----------------------------------------------------------------------------*/
    for ( i = 0; i < G->rows; i++){
        marked[i] = (G->rowptr[i+1]-G->rowptr[i]);
    }


    /*-----------------------------------------------------------------------------
     *  allocate otehr stuff
     *-----------------------------------------------------------------------------*/
    for ( i = 0; i < n  ; i++ ) V[i] = i;
    sV = n;

    ci = 0;
    blocks [0 ] = 0 ;
    *nblocks = 0 ;
    while ( sV > 0 ) {
        ci ++ ;
        /*-----------------------------------------------------------------------------
         *  get the minimal out degree , can be done faster by only using V
         *-----------------------------------------------------------------------------*/
        minv = marked[V[0]]; mini = V[0];
        for ( i = 0 ; i < sV ; i ++) {
            if ( minv > marked[V[i]]) {
                minv = marked[i];
                mini = V[i] ;
            }
        }

        /* minv = marked[0]; mini = 0;

           for ( i = 0 ; i < n ; i ++) {
           if ( minv > marked[i]) {
           minv = marked[i];
           mini = i  ;
           }
           }
           */

        /*-----------------------------------------------------------------------------
         *  DFS on G starting from mini
         *-----------------------------------------------------------------------------*/
        lblock_dfs(G, mini, dfs_stack, marked, n+1, &sT);


        /*-----------------------------------------------------------------------------
         *  debug
         *-----------------------------------------------------------------------------*/
        /* MSG_PRINT("DFS from " MESS_PRINTF_INT " in run " MESS_PRINTF_INT "\n", mini, ci );
           for ( i = sT; i < n; i++) {
           MSG_PRINT(" " MESS_PRINTF_INT " ", dfs_stack[i]);
           }
           MSG_PRINT("\n");
           */

        /*-----------------------------------------------------------------------------
         *  add to the sets
         *-----------------------------------------------------------------------------*/
        // MSG_PRINT("n-ts: &" MESS_PRINTF_INT "  \n", n-sT );
        for ( i = sT; i < n ; i++) {
            perm[pperm++] = dfs_stack[i];
            swap(V,sV-1,dfs_stack[i]);
            sV--;
        }
        blocks[*nblocks+1] = pperm;
        *nblocks+=1;

        // MSG_PRINT("sV = " MESS_PRINTF_INT " \n", sV);
        // MSG_PRINT("V=");
        // for ( i = 0; i < sV ; i ++ ) {
        // MSG_PRINT ( " " MESS_PRINTF_INT " " , V[i]);
        // }
        // MSG_PRINT("\n");
    }
    mess_free(V);
    mess_free(marked);
    mess_free(dfs_stack);
    return 0;
}       /* -----  end of function mess_graph_lblockpart  ----- */



