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
 * @file lib/graph/dfs.c
 * @brief Depth-first-search on a graph and reachability set.
 * @author @koehlerm
 * The functions are based on @cite Dav06.
 *
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


/**
 * @brief Perform a depth-first-search on a directed graph defined by a matrix.
 * @param[in]  G     input matrix defining graph
 * @param[in] pinv   input inverse row permutation of G \n
 *                        ( NULL if it does not exist)
 * @param[in] j      input start node for the DFS
 * @param[out] top  index of the element in the chi stack
 * @param[out] xi   stack \f$ xi [top, \cdots ,n-1] \f$ containing set of reached nodes
 * @param[in] pstack     input work array of dimension \f$ n \f$
 * @return 0 on success or a non zero error code otherwise
 *
 * The @ref mess_graph_dfs function performs a depth-first-search (DFS) on a graph \f$ G \f$ starting at node \f$ j \f$.\n
 * The reached nodes are stored in the \f$ xi \f$ array \f$ xi[top,\cdots ,n-1] \f$.
 *
 * The functions are based on @cite Dav06.
 *
 * \attention It only works for Compressed Sparse Column matrices.
 *
 */
int mess_graph_dfs (mess_matrix G, mess_int_t j, mess_int_t* top, mess_int_t* xi, mess_int_t *pstack, const mess_int_t *pinv)
{
    MSG_FNAME(__func__);
    mess_int_t i, p, p2, jnew, head =0;
    int done = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(G);
    mess_check_nullpointer(top);
    mess_check_nullpointer(xi);
    mess_check_nullpointer(pstack);
    mess_check_csc(G);
    mess_check_square(G);
    if ( j < 0 || j >= G->cols){
        MSG_ERROR("j = " MESS_PRINTF_INT "\n is out of range\n", j);
        return MESS_ERROR_DIMENSION;
    }

    /*-----------------------------------------------------------------------------
     *  perform DFS
     *-----------------------------------------------------------------------------*/
    xi [0] = j ;   /* start on node j */
    while (head >= 0) { /*  until any node is left */
        j = xi [head] ;         /* get j from the top of the recursion stack */
        jnew = pinv ? (pinv [j]) : j ;
        if (!MARKED (G->colptr, j)) {
            MARK (G->colptr, j) ;       /* mark node j as visited */
            pstack [head] = (jnew < 0) ? 0 : UNFLIP (G->colptr [jnew]) ;    /* start of the neighbors */
        }
        done = 1 ;                  /* assume node j done if no unvisited neighbors */
        p2 = (jnew < 0) ? 0 : UNFLIP (G->colptr[jnew+1]) ;
        for (p = pstack [head] ; p < p2 ; p++){  /* examine all neighbors of j */
            i = G->rowptr [p] ;            /* consider neighbor node i */
            if (MARKED (G->colptr, i)) continue ;   /* skip visited node i */
            pstack [head] = p ;     /* pause depth-first search of node j */
            xi [++head] = i ;       /* start dfs at node i */
            done = 0 ;              /* node j is not done */
            break ;                 /* break, to start dfs (i) */
        }
        if (done) {               /* depth-first search at node j is done */
            head-- ;            /* remove j from the recursion stack */
            *top = *top -1;
            xi [*top] = j ;    /* and place in the output stack */
        }
    }
    return 0 ;
}


/**
 * @brief Compute the reachability set from a node in a graph.
 * @param[in] G      input matrix defining graph
 * @param[in] B     input matrix or column with the starting points
 * @param[in] k         input number of the column in \f$ B \f$ which should be used
 * @param[out] top  index of the top element of the stack \f$ xi \f$
 * @param[out] xi   stack with the result \n
 *                  (length \f$ 2*n \f$ !!!!!)
 * @param[in] pinv   input inverse row permutation for \f$ G  \f$
 * @return 0 on success or a non zero error code otherwise
 *
 * The @ref mess_graph_reach function calculates the reachability set of the nodes \f$ B(:,k) \f$ over the
 * graph \f$ G \f$. \n
 * The reached nodes are stored in the \f$ xi \f$ array \f$ xi[ top,\cdots ,n-1] \f$.
 *
 * The functions are based on @cite Dav06.
 *
 * \attention It only works for Compressed Sparse Column matrices.
 *
 */
int mess_graph_reach (mess_matrix G, mess_matrix B, mess_int_t k, mess_int_t *top, mess_int_t *xi, const mess_int_t *pinv)
{
    MSG_FNAME(__func__);
    mess_int_t p, n;
    int ret = 0;


    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(G);
    mess_check_nullpointer(B);
    mess_check_nullpointer(top);
    mess_check_nullpointer(xi);
    mess_check_csc(G);
    mess_check_csc(B);
    mess_check_square(G);
    if(k < 0 || k >= B->cols) {
        MSG_ERROR("k = " MESS_PRINTF_INT " is out of range\n",k);
        return MESS_ERROR_DIMENSION;
    }
    n = G->cols ;
    *top = n ;

    /*-----------------------------------------------------------------------------
     *  perfom DFS on each node
     *-----------------------------------------------------------------------------*/
    for (p = B->colptr [k] ; p < B->colptr [k+1] ; p++) {
        if (!MARKED (G->colptr, B->rowptr [p])){    /* start a dfs at unmarked node i */
            ret = mess_graph_dfs (G,B->rowptr[p],top, xi, &(xi[n]), pinv) ;
            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_graph_dfs);
        }
    }
    for (p = *top ; p < n ; p++) MARK (G->colptr, xi [p]) ;  /* restore G */
    return 0 ;
}
/* -----  end of function mess_graph_reach  ----- */
