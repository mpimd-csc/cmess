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
 * @file lib/reorder/cuthill.c
 * @brief Cuthill MC Kee Reordering.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include <complex.h>

/* return value convention */
#define SUCCESS 0
#define FAIL    1

/* macro for malloc failure returns */
#define TRY_MALLOC(__who__, __somuch__, __type__, __func__) {       \
    mess_try_alloc2(__who__, __type__ *, sizeof(__type__) * (__somuch__)); \
    if ((__who__)==NULL)                        \
    {                                   \
        fprintf(stderr,                         \
                "line-%d in %s, file %s: malloc for %d %s entries to %s failed\n",\
                __LINE__, #__func__, __FILE__, (int)(__somuch__), #__type__, #__who__ );    \
        return FAIL;                            \
    }                                   \
}

/* macro for malloc failure returns */
#define TRY_REALLOC(__who__, __somuch__, __type__, __func__) {      \
    mess_try_realloc2(__who__, __type__ *, (__somuch__) * sizeof(__type__)); \
    if ((__who__)==NULL)                        \
    {                                   \
        fprintf(stderr,                         \
                "line-%d in %s, file %s: realloc for %d %s entries to %s failed\n", \
                __LINE__, #__func__, __FILE__, (int)(__somuch__), #__type__, #__who__ ); \
        return FAIL;                            \
    }                                   \
}

/* macro to handle error returns by subroutines,
   __where__ = function in which this happens
   __who__   = the subroutine called
   __err__   = variable holding the return value of the subroutine */
#undef FUNCTION_FAILURE_HANDLE
#define FUNCTION_FAILURE_HANDLE( __err__, __who__, __where__) {\
    if ((__err__)!=SUCCESS) \
    {\
        fprintf(stderr, \
                "line-%d in %s, file %s: call to %s returned with %d(!=SUCCESS)\n",\
                __LINE__, #__where__, __FILE__, #__who__, __err__ );\
        return FAIL; \
    }\
}

/**
 * @internal
 *
 * Data structure for adjency-, level- and similar lists of integers.
 * @attention Internal use only.
 * */
struct int_list
{
    mess_int_t *ix;           /**< index vector,  j=ix[i] marks the begin of the i-th part of el, k=ix[i+1]-1 marks the end of the i-th part of el, i<ixn<=ixnmax */
    mess_int_t ixn;           /**< Actual length of the index vector \f$ ix \f$ */
    mess_int_t ixnmax;        /**< Maximal length of the index vector \f$  ix \f$ */
    mess_int_t *el;           /**< list of elements, with j and k as described in the ix part, el[j+t] is the t-th element of the i-th part of the list, t<=k-j-1 */
    mess_int_t elnmax;        /**< maximal length of el, therefore   ix[ixn+1]-1<=elnmax has to be asured */
};

/**
 * @brief Allocate memory in the @ref int_list structure for a given size.
 * @param[out] L         output @ref int_list structure
 * @param[in] ixnmax     input maximal number of parts of the list
 * @param[in] elnmax     input maximal number of elements in the list
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref int_list_alloc function allocates memory in the @ref int_list structure for a
 * given size defined by the maximal number
 * of parts \f$ ixnmax \f$ and the maximal number of elements \f$ elnmax \f$.
 *
 */
static int int_list_alloc( struct int_list *L, mess_int_t ixnmax, mess_int_t elnmax  ){
    MSG_FNAME(__func__);
    TRY_MALLOC( L->ix, ixnmax+1, mess_int_t, int_list_alloc);
    TRY_MALLOC( L->el, elnmax, mess_int_t, int_list_alloc);

    L->ixn=0;
    L->ixnmax=ixnmax;
    L->elnmax=elnmax;

    return SUCCESS;
}

/**
 * @brief Free memory in the @ref int_list structure.
 * @param[in,out] L     input/output @ref int_list structure
 * @return Nothing
 * The @ref int_list_free function frees memory in the @ref int_list structure, so it can be freed.
 *
 */
static void int_list_free( struct int_list *L ){
    mess_free(L->ix);
    mess_free(L->el);

    L->ixn=0;
    L->ixnmax=0;
    L->elnmax=0;
}

/**
 * @internal
 * @brief Compare two integers.
 * @param[in] a1        input pointer to first integer \f$a_1\f$
 * @param[in] b1        input pointer to second integer \f$b_1\f$
 * @return \f$ a_1 - b_1 \f$
 *
 * The @ref comp_mess_int_t function is a helper function for qsort and compares two integers. \n
 * This function returns  \f$ a_1 - b_1 \f$.
 *
 * @attention Interal use only.
 */
static int comp_mess_int_t(const void *a1, const void *b1 ){
    mess_int_t *a, *b;

    a=(mess_int_t *) a1;
    b=(mess_int_t *) b1;

    return *a-*b;
}


/* Forward definitions   */
static int bandw_red_perm( mess_int_t nr, mess_int_t *adix, mess_int_t *adel,
        mess_int_t *perm, mess_int_t *iperm,
        mess_int_t *bandw);

static int bandw_mark_connected( mess_int_t label, mess_int_t start, mess_int_t *adix,mess_int_t *adel, mess_int_t n1, mess_int_t *marked);
static int bandw_gipost( mess_int_t n1, mess_int_t n2, struct int_list *ad, mess_int_t *perm, mess_int_t *iperm, mess_int_t *bandw);
static int bandw_level_struct( mess_int_t root, mess_int_t n1, mess_int_t n2, struct int_list *ad, struct int_list *L, mess_int_t *hint);
static int bandw_number_adj_vert(struct int_list *ad, struct int_list *Lv,
        mess_int_t *level, mess_int_t l1, mess_int_t l2, mess_int_t n1, mess_int_t n2,
        mess_int_t *perm, mess_int_t *iperm, mess_int_t *current,
        mess_int_t *sorter);

/**
 * @brief Compute the Reverse Cuthill-McKee reordering of a matrix.
 * @param[in]  A        input matrix \f$A\f$
 * @param[out] perm     output permutation containing the Reverse Cuthill-McKee reordering
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_matrix_reorder_rcm function calls the Gibbs-Poole-Stockmeyer modification of the
 * Cuthill-McKee algorithm to produce a bandwith reducing permutaion \f$ perm \f$
 * for a given adjacency list \f$ adel \f$ of a sparse matrix.
 *
 * @sa mess_matrix_reorder_amd
 * @sa mess_matrix_reorder_colamd
 * @sa mess_matrix_reorder
 *
 */
int mess_matrix_reorder_rcm(mess_matrix A, mess_int_t * perm)
{
    MSG_FNAME(__func__);
    mess_matrix AT;
    mess_matrix eye;
    int ret;
    mess_int_t *perm2;
    mess_int_t i;
    mess_int_t bw;


    mess_check_nullpointer(A);
    mess_check_nullpointer(perm);
    mess_check_square(A);

    if ( MESS_IS_DENSE(A)) {
        for (i = 0; i < A->rows; i++) {
            perm[i] = i;
        }
        return 0;
    }

    ret = mess_matrix_init(&AT);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_matrix_init );
    ret = mess_matrix_init(&eye);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_eye(eye, A->rows, A->cols, A->store_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_eye);
    ret = mess_matrix_ctranspose(A, AT);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
    ret = mess_matrix_add(1, eye, 1, AT);   FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    ret = mess_matrix_add(1, A, 1, AT);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_add);
    mess_matrix_clear(&eye);


    mess_try_alloc( perm2, mess_int_t *, sizeof(mess_int_t) * (A->rows));
    if (perm2 == NULL ) {
        MSG_ERROR("Allocate perm2 failed\n");
        mess_matrix_clear(&AT);
        return MESS_ERROR_MEMORY;
    }

    ret = 0;
    if (MESS_IS_CSR(AT)){
        ret = bandw_red_perm(AT->rows, AT->rowptr, AT->colptr, perm2, perm, &bw);
    } else if ( MESS_IS_CSC(AT)) {
        ret = bandw_red_perm(AT->cols, AT->colptr, AT->rowptr, perm2, perm, &bw);
    } else {
        MSG_ERROR("Storage type %s is not supported\n", mess_storage_t_str(AT->store_type));
        mess_matrix_clear(&AT);
        mess_free(perm2);
        return MESS_ERROR_STORAGETYPE;
    }
    mess_free(perm2);
    mess_matrix_clear(&AT);

    if ( ret != SUCCESS) {
        MSG_ERROR("RCM reordering failed.\n");
        return MESS_ERROR_GENERAL;
    }
    return 0;
}


/**
 * @internal
 * @brief Compute bandwidth and corresponding permutation.
 * @param[in] nr        input number of vertices in graph
 * @param[in] adix      input index of adjacency list
 * @param[in] adel      input element vector of adjacency list
 * @param[out] perm     output bandwidth reducing permuation
 * @param[out] iperm    output inverse permutation of \f$ perm \f$
 * @param[out] bandw    output resulting bandwidth
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref bandw_red_perm function calls the Gibbs-Poole-Stockmeyer modification of the
 * Cuthill-McKee algorithm to produce a bandwith reducing permutaion \f$ perm \f$
 * for a given adjacency list \f$ adel \f$ of a sparse matrix.
 *
 * The element vector of the adjacency list \f$ adel \f$ need to be a full adjacency list, that
 * means it has to contain edges \f$ (i,j) \f$ and \f$ (j,i) \f$. \n
 * The permuation \f$ perm \f$ and its inverse permutation \f$ iperm \f$ are given as integer
 * vectors of length \f$ nr \f$.
 *
 * @attention Internal use only.
 */
static int bandw_red_perm( mess_int_t nr, mess_int_t *adix, mess_int_t *adel,mess_int_t *perm, mess_int_t *iperm,mess_int_t *bandw ) {
    MSG_FNAME(__func__);
    int  err;
    mess_int_t i, j;
    mess_int_t *ihelp, *npart, tc, nc;
    struct int_list parts, pcon;

    /* quick return? */
    if (nr==0) return SUCCESS;

    /* check if the graph is connected, if not find the connected parts
    */
    TRY_MALLOC(ihelp, nr, mess_int_t, bandw_red_perm);
    for (i=0; i<nr; i++) ihelp[i]=0;
    nc=0;
    for (i=0; i<nr; i++)
    {
        /* check if this vertex starts a new part */
        if (ihelp[i]==0)
        {
            /* start a new part, mark all connected nodes recursively */
            nc++;
            err=bandw_mark_connected(nc, i, adix, adel, 0, ihelp);
            FUNCTION_FAILURE_HANDLE(err, bandw_mark_connected,
                    bandw_red_perm);
        }
    }

    /* MSG_PRINT("  test:  nc=%3d, nr=%3d\n", nc, nr); */

    err=int_list_alloc(&parts, nc, nr);
    FUNCTION_FAILURE_HANDLE(err, int_list_alloc, bandw_red_perm);
    /* get the number of vertices in each conected subgraph */
    for (i=0; i<nc+1; i++) parts.ix[i]=0;
    for (i=0; i<nr; i++) parts.ix[ihelp[i]]++;
    /* build a index list from that */
    for (i=0; i<nc; i++) parts.ix[i+1]+=parts.ix[i];
    /* build the element list */
    TRY_MALLOC( npart, nc, mess_int_t, bandw_red_perm);
    for (i=0; i<nc; i++) npart[i]=0;
    for (i=0; i<nr; i++)
    {
        tc=ihelp[i]-1;  /* part numbering starts at 1, index numb. at 0 */
        parts.el[parts.ix[tc]+npart[tc]]=i;
        npart[tc]++;
    }
    mess_free(npart);

#ifdef MESS_DEBUG
    if (nr!=parts.ix[nc])
    {
        MSG_ERROR("bandw_red_perm: nr(whole)!=sum(nr(parts))\n");
        return FAIL;
    }
#endif


    /* now build a full connectivity list as labeled by the permutation
    */
    err=int_list_alloc(&pcon, nr, 2*adix[nr]);
    FUNCTION_FAILURE_HANDLE(err, int_list_alloc, bandw_red_perm);

    /* get the degree of each vertex */
    pcon.ix[0]=0;
    for (i=0; i<nr; i++) pcon.ix[parts.el[i]+1]=adix[i+1]-adix[i];


    /* build a index list from that */
    for (i=0; i<nr; i++) pcon.ix[i+1]+=pcon.ix[i];

#ifdef MESS_DEBUG
    /* if (pcon.ix[nr]!=2*adix[nr]) // old for half adj-list */
    if (pcon.ix[nr]!=adix[nr])
    {
        MSG_ERROR("bandw_red_perm: sum(degree) mismatch\n");
        return FAIL;
    }
#endif

    /* build the element list the fill counter goes into ihelp */
    for (i=0; i<nr; i++) ihelp[i]=0;
    for (i=0; i<nr; i++)
        for (j=adix[i]; j<adix[i+1]; j++)
        {
            mess_int_t a,b;
            /* the two vertices of the edge */
            a=parts.el[i];
            b=parts.el[adel[j]];

            /* MSG_PRINT("i=%3d  j=%3d   a=%3d  b=%3d   ", i, j, a, b);
               MSG_PRINT("pcon.ix[a]=%3d   ihelp[a]=%3d   pcon.ix[b]=%3d   "
               "ihelp[b]=%3d\n",
               pcon.ix[a], ihelp[a], pcon.ix[b], ihelp[b]); */
            /* add the two copies of the edge to the full list */
            pcon.el[pcon.ix[a]+ihelp[a]]=b;
            ihelp[a]++;
            /* old for half adj-list:
               pcon.el[pcon.ix[b]+ihelp[b]]=a;
               ihelp[b]++; */
        }

    /* now call the bandwidth reduction for each of the subgraphs, store
       the permutations of the subgraphs in ihelp, get the max of the
       subgraph bandwidths */
    *bandw=0;
    for (i=0; i<nc; i++)
    {
        mess_int_t sbandw;
        err=bandw_gipost(parts.ix[i], parts.ix[i+1], &pcon,
                ihelp, iperm, &sbandw);
        FUNCTION_FAILURE_HANDLE(err, bandw_gipost, bandw_red_perm);

        /* overall bandwidth is max of subgraph bandwidths */
        if (sbandw>*bandw) *bandw=(mess_int_t) sbandw;
    }

    /* the element list parts.el itself forms a permutation, which
       collects the individual subgraphs one after each other,
       the overall bandwith reducing perutation is the composition of
       parts.el and the permutation now in ihelp */
    for (i=0; i<nr; i++) perm[i]=ihelp[parts.el[i]];

    /* the inverse permutation */
    for (i=0; i<nr; i++) iperm[perm[i]]=i;

    /* done, free local memory and return */
    mess_free(ihelp);
    int_list_free(&parts);
    int_list_free(&pcon);

    return SUCCESS;
}


/**
 * @internal
 * @brief Mark vertices with a label.
 * @param[in] label         input label of start vertex and all connected vertices
 * @param[in] start         input vertex to start labeling
 * @param[in] adix          input vector marking beginning and end of each node in adjacency list in \f$ adel \f$
 * @param[in] adel          input vector defining adjacency of nodes together with \f$ adix \f$
 * @param[in] n1            input offset of lowest node number in considered subgraph
 * @param[in,out] marked    input/output vector containing marked vertices
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref bandw_mark_connected function marks the start vertex \f$ start \f$
 * and all vertices connected to start with \f$ label \f$ and stores them in \f$  marked  \f$.\n
 * The ouput vector \f$ marked \f$ has a length of the number of vertices .
 *
 * @attention Internal use only.
 */
static int bandw_mark_connected( mess_int_t label, mess_int_t start, mess_int_t *adix,
        mess_int_t *adel, mess_int_t n1, mess_int_t *marked )
{
    MSG_FNAME(__func__);
    mess_int_t i, j;

    /* Stack for the problems  */
    mess_int_t *startnodes;
    mess_int_t alen, len;

    alen = 1000;
    TRY_MALLOC(startnodes, alen, mess_int_t, __func__);
    if ( startnodes == NULL ) {
        return FAIL;
    }
    len = 1;
    startnodes[0] = start;

    while ( len > 0 ) {
        start = startnodes[len-1]; len--;

        marked[start-n1]=label;

        for (j=adix[start]; j<adix[start+1]; j++)
        {
            i=adel[j];
            if ((marked[i-n1]>0)&&(marked[i-n1]!=label))
            {
                fprintf(stderr, "bandw_mark_connected: found vertex allready "
                        "marked[%d]=%d  != %d=label\n",
                        (int) i, (int) marked[i-n1], (int) label);
                mess_free(startnodes);
                return FAIL;
            }
            if (marked[i-n1]==0)
            {
                if ( len >= alen ) {
                    alen += 2000;
                    TRY_REALLOC(startnodes, alen, mess_int_t, __func__);
                    if (startnodes == NULL) {
                        return FAIL;
                    }
                }
                // MSG_PRINT("len = %d \n", len );
                startnodes[len] = i;
                len++;
                /* err=bandw_mark_connected(label, i, adix, adel, n1, marked);
                   FUNCTION_FAILURE_HANDLE(err, bandw_mark_connected,
                   bandw_mark_connected); */
            }
        }

    }
    mess_free(startnodes);
    return SUCCESS;
}

/**
 * @internal
 * @brief Compute bandwidth and corresponding permutation of a connected subgraph.
 * @param[in] n1            input index of first node of connected subgraph
 * @param[in] n2            input index of last node \f$+1 \f$ of connected subgraph
 * @param[in] ad            input adjacency list of whole graph
 * @param[in] perm          input bandwidth reducing permuation
 * @param[in] iperm         input reverse of \f$ perm \f$
 * @param[in,out] bandw     input/output bandwith for subgraph
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref bandw_gipost function applies the Gibbs-Poole-Stockmeyer modification of the
 * Cuthill-McKee algorithm to produce a bandwith reducing permutaion
 * for a given full adjency list of a connected subgraph of a sparse matrix.
 *
 * In the adjacency list \f$ ad \f$ only the \f$ n1 \f$ to \f$ n2-1 \f$ parts are read and in the
 * (inverse of the) bandwidth reducing permutation the \f$ n1 \f$ to \f$ n2-1 \f$ entries are written to
 * an integer vector.
 *
 * @attention Internal use only.
 */
static int bandw_gipost( mess_int_t n1, mess_int_t n2, struct int_list *ad,
        mess_int_t *perm, mess_int_t *iperm, mess_int_t *bandw  ){
    MSG_FNAME(__func__);
    int  err;
    mess_int_t i, j, nc, tmp, length, last, u =0, v = 0, tv, current;
    mess_int_t width, widthv, widthu;
    mess_int_t from, to , dir, lvl;

    mess_int_t nLi, nLimax; /* Li is going to hold a set of rooted level structures,
                               nLi will mark the number of defined ones in it,
                               nLimax the number of level structs for which
                               memory is allocated in Li
                               Lu and Lv will eventually be the level
                               structures for the endpoints u and v of a pseudo
                               diameter of the graph,
                               LN will be the final level structure which is
                               then used to define the numbering */
    struct int_list **Li, *Lu, *Lv;

    /* integer help vectors */
    mess_int_t *hint, *hintv, *hintu, *sorter, *Nix, *Lix, *Hix;


    /* test the trivial case, because it breaks stuff */
    if (n2-n1==1)
    {
        /* only one node in this graph, so there is only one permutation
        */

        perm[n1]=n1;
        iperm[n1]=n1;
        (*bandw)=1;

        return SUCCESS;
    }

    /* allocate memory */
    TRY_MALLOC( hint, n2-n1, mess_int_t, bandw_gipost);
    TRY_MALLOC( hintv, n2-n1, mess_int_t, bandw_gipost);
    TRY_MALLOC( hintu, n2-n1, mess_int_t, bandw_gipost);
    nLimax=30;
    TRY_MALLOC( Li, nLimax, struct int_list*, bandw_gipost);
    nLi=0;

    TRY_MALLOC( Lv, 1, struct int_list, bandw_gipost);

    /*************************************************************/
    /*                                                           */
    /* Stage 1: find a pseudo diameter to give starting vertices */
    /*                                                           */
    /*************************************************************/

    /* find a vertex v with min degree */
    tmp=n2-n1; /* init mindegree=max possible degree */
    for (i=n1; i<n2; i++)
        if (ad->ix[i+1]-ad->ix[i]<tmp)
        {
            tmp=ad->ix[i+1]-ad->ix[i];
            v=i;
        }

    /*   MSG_PRINT(" test: min degree=%d\n", tmp); */

    err=bandw_level_struct( v, n1, n2, ad, Lv, hint);
    FUNCTION_FAILURE_HANDLE(err, bandw_level_struct, bandw_gipost);

    last=0;
    while (last==0)
    {
        /* number of possible other endpoints of the diameter = number
           vertices in the highest level of the level struct */
        nLi = Lv->ix[Lv->ixn+1] - Lv->ix[Lv->ixn];
        if (nLi>nLimax)
        {
            nLimax=nLi;
            TRY_REALLOC( Li, nLimax, struct int_list*, bandw_gipost);
        }



        /* length of the pseudo diameter */
        length= Lv->ixn;

        /* sort the other possible endpoints by ascending degree */
        TRY_MALLOC(sorter, 2*nLi, mess_int_t, bandw_gipost);
        for (i=0; i< Lv->ix[Lv->ixn+1] - Lv->ix[Lv->ixn]; i++)
        {
            mess_int_t vcandi=Lv->el[ i + Lv->ix[Lv->ixn]];
            /* in sorter: first degree, then number of the vertex */
            sorter[2*i+0]=ad->ix[vcandi+1]-ad->ix[vcandi];
            sorter[2*i+1]=vcandi;
        }
        qsort( sorter, nLi, 2*sizeof(mess_int_t), comp_mess_int_t);

        /* pretent this is the last try */
        last=1;
        /* for each of these build the level structure */
        for(i=0; (i<nLi)&&(last==1); i++)
        {
            mess_int_t vcandi=sorter[2*i+1];
            TRY_MALLOC( Li[i], 1, struct int_list, bandw_gipost);
            /* MSG_PRINT(" test:  sorter[%2d]=  [ %3d, %3d ],   ", i, sorter[2*i+0],
               sorter[2*i+1]);  */
            err=bandw_level_struct( vcandi, n1, n2, ad, Li[i], hint);
            FUNCTION_FAILURE_HANDLE(err, bandw_level_struct, bandw_gipost);

            /* if the current endpoint gives a larger diameter
               candidate, reject u, set u=this candidate, start again */
            if (Li[i]->ixn > length)
            {
                int_list_free(Lv);
                mess_free(Lv);
                v=vcandi;
                Lv=Li[i];
                /* free all prior attempts */
                for (j=0; j<i; j++)
                {
                    int_list_free(Li[j]);
                    mess_free(Li[j]);
                }
                mess_free(sorter);
                /* next try */
                last=0;
            }
        }
    }
    /* now v and Lv are fixed, define u to be the root of the Li with
       the smallest width */
    widthu=n2-n1; /* init min(width)=max possible width */
    for (i=0; i<nLi; i++)
    {
        width=0;
        if (Li[i]->ixn == length)
        {
            for (j=0; j<=length; j++)
            {
                if (Li[i]->ix[j+1]-Li[i]->ix[j]>width)
                {
                    width=Li[i]->ix[j+1]-Li[i]->ix[j];
                }
            }
            if (width<widthu)
            {
                widthu=width;
                u=i;
            }
        }
    }
    /* MSG_PRINT("  test: minwidth=%d  for  Li[%d]\n", widthu, u);  */
    /* by now u is the number of the Li width the smallest width, so
       define Lu and u correctly */
    Lu=Li[u];
    Li[u]=NULL;
    u=sorter[2*u+1];
    /* now forget the remaining Li */
    for (i=0; i<nLi; i++)
    {
        if (Li[i]!=NULL)
        {
            int_list_free(Li[i]);
            mess_free(Li[i]);
        }
    }
    mess_free(Li);
    nLi=0;
    nLimax=0;
    mess_free(sorter);

    /* calculate width of v */
    widthv=n2-n1; /* init min(width)=max possible width */
    {
        width=0;
        if (Lv->ixn == length)
        {
            for (j=0; j<=length; j++)
            {
                if (Lv->ix[j+1]-Lv->ix[j]>width)
                {
                    width=Lv->ix[j+1]-Lv->ix[j];
                }
            }
            if (width<widthv)
            {
                widthv=width;
            }
        }
    }


    /*************************************************************/
    /*                                                           */
    /* Stage 2: minimizing the level width                       */
    /*                                                           */
    /*************************************************************/

    /* MSG_PRINT("\n Stage 2: test \n\n"); */

    /* rebuild Lu and Lv since we need their hintu and hintv */
    int_list_free(Lu);
    err=bandw_level_struct( u, n1, n2, ad, Lu, hintu);
    FUNCTION_FAILURE_HANDLE(err, bandw_level_struct, bandw_gipost);
    int_list_free(Lv);
    err=bandw_level_struct( v, n1, n2, ad, Lv, hintv);
    FUNCTION_FAILURE_HANDLE(err, bandw_level_struct, bandw_gipost);
    /* hintv and hint now hold the pairs of levels for each vertex */

    /* build N, count the vertices whoes level is clear by the pair of
       levels */
    TRY_MALLOC(Nix, length+1, mess_int_t, bandw_gipost);
    TRY_MALLOC(Hix, length+1, mess_int_t, bandw_gipost);
    TRY_MALLOC(Lix, length+1, mess_int_t, bandw_gipost);
    for(i=0; i<=length; i++) Nix[i]=0;
    for (i=0; i<n2-n1; i++)
    {
        if (hintv[i]==length-hintu[i])
        {
            /* MSG_PRINT("  test:   vertex[%3d]: lvl=%d\n", i, hintv[i]); */
            /* take it out of the game and count it in this level */
            hint[i]=-1;
            Nix[hintv[i]]++;
        }
        else
        {
            /* needs to be assigned to an connected subgraph */
            hint[i]=0;
        }
    }

    /* partition the undecided vertices into connected subgraphs */
    nc=0;
    for (i=0; i<n2-n1; i++)
    {
        /* check if this vertex starts a new part */
        if (hint[i]==0)
        {
            /* start a new part, mark all connected nodes recursively */
            nc++;
            err=bandw_mark_connected(nc, i+n1, ad->ix, ad->el, n1, hint);
            FUNCTION_FAILURE_HANDLE(err, bandw_mark_connected,
                    bandw_gipost);
        }
    }
    /*   MSG_PRINT("  test: %d connected subgraphs need to be assigned to a level\n", */
    /*   nc); */

    TRY_MALLOC(sorter , 2*nc, mess_int_t, bandw_gipost);
    for (i=0; i<nc; i++)
    {
        sorter[2*i  ]=0; /* #vertices in this subgraph */
        sorter[2*i+1]=i;
    }
    for (i=0; i<n2-n1; i++)
        if (hint[i]>0) sorter[(hint[i]-1)*2]++;
    qsort( sorter, nc, 2*sizeof(mess_int_t), comp_mess_int_t);
    /* now subgraphs sorted in ascending order of #vertices, start with
       largest subgraph to decide which of the two alternative level
       sets to use */

    /* for (i=0; i<n2-n1; i++) {
       MSG_PRINT(" test:   vertex[%3d]: subgraph=%3d lvl_v=%3d   lvl_u=%3d\n",
       i, hint[i], hintv[i], hintu[i]);
       }  */

    for (j=nc-1; j>=0; j--)
    {
        mess_int_t Hixmax, Lixmax, choose;
        /* count the vertex numbers per level which would result from
           putting the undecided in the level as per v resp. per u */

        /*       MSG_PRINT(" test: subgraph[%2d]= %3d    #vx=%3d\n", j, sorter[2*j+1], */
        /*       sorter[2*j]); */
        for(i=0; i<=length; i++) { Hix[i]=0; Lix[i]=0; }
        for(i=0; i<n2-n1; i++)
        {
            /* if in this subgraph */
            if (hint[i]-1==sorter[2*j+1])
            {
                Hix[hintv[i]]++;
                Lix[length-hintu[i]]++;
            }
        }
        /*get the maximum of each Lix and Hix (with Nix added to each)*/
        Hixmax=0; Lixmax=0;
        for(i=0; i<=length; i++)
        {
            if (Lix[i]+Nix[i]>Lixmax) Lixmax=Lix[i]+Nix[i];
            if (Hix[i]+Nix[i]>Hixmax) Hixmax=Hix[i]+Nix[i];
            /*    MSG_PRINT(" test: Nix[%3d]=%3d   Lix[]=%3d     Hix[]=%3d\n", */
            /*       i, Nix[i], Lix[i], Hix[i]); */
        }
        if (Hixmax<Lixmax)
        {
            choose=1;
        }
        else if (Hixmax>Lixmax)
        {
            choose=2;
        }
        else
        {
            /* equally good, decide by the width of Lv and Lu */
            if (widthv<=widthu) choose=1;
            else choose=2;
        }
        /*       MSG_PRINT(" test:  choose=%d  Hixmax=%d, Lixmax=%d    " */
        /*       "widthv=%d, widthu=%d\n", */
        /*       choose, Hixmax, Lixmax, widthv, widthu); */
        /* now remember what we have decided by modifying hintv and
           hintu */
        for(i=0; i<n2-n1; i++)
        {
            /* if in this subgraph */
            if (hint[i]-1==sorter[2*j+1])
            {
                if (choose==1) hintu[i]=length-hintv[i];
                else hintv[i]=length-hintu[i];
            }
        }
        /* update Nix */
        for(i=0; i<=length; i++)
        {
            if (choose==1) Nix[i]+=Hix[i];
            else Nix[i]+=Lix[i];
        }
    }
    /* all these guys have now been assigned a level, as can be found in
       hintv[i], the number of vertices in each level can be found in,
       now adjust Lv to give this level structure */
    Lv->ix[0]=0;
    for(i=0; i<=length; i++)
    {
        Lv->ix[i+1]= Lv->ix[i]+Nix[i];
        Nix[i]=0;
    }
    for(i=0; i<n2-n1; i++)
    {
        mess_int_t lvl2=hintv[i];
        /* MSG_PRINT(" test: lvl[%3d]=%3d   Lv->ix=%3d Nix[lvl]=%3d\n",
           i, lvl, Lv->ix[lvl], Nix[lvl]);  */
        Lv->el[Lv->ix[lvl2] + Nix[lvl2]]=i+n1;
        Nix[lvl2]++;
    }

    /* sorter, Lu, Nix, Hix, Lix not needed any more */
    int_list_free(Lu);
    mess_free(Lu);
    mess_free(Nix);
    mess_free(Hix);
    mess_free(Lix);
    mess_free(sorter);
    mess_free(hintu);
    mess_free(hint);


    /*************************************************************/
    /*                                                           */
    /* Stage 3: numbering                                        */
    /*                                                           */
    /*************************************************************/

    /* check the degree of v and u to decide at which end of the
       levelstructure to start */
    if (ad->ix[v+1]-ad->ix[v] <= ad->ix[u+1]-ad->ix[u])
    {
        /* v lower degree, so go upwards in the v levels */
        from=0;
        to=length;
        dir=1;
        tv=v;
    }
    else
    {
        /* v higher degree, so downwards in the v levels */
        from=length;
        to=0;
        dir=-1;
        tv=u;
    }
    /* init the permutation (as undefined) */
    for (i=0; i<n2-n1; i++)
    {
        perm[i+n1]  = -1;
        iperm[i+n1] = -1;
    }
    /* allocate workspace for sorting */
    width=widthu;
    if (widthv>width) width=widthv;
    TRY_MALLOC( sorter, 2*width, mess_int_t, bandw_gipost);

    current=n1;
    /* start with tv, assign all elements of level from with consecutive
       numbers */
    perm[tv]=current;
    iperm[current]=tv;
    current++;
    err=bandw_number_adj_vert( ad, Lv, hintv, from, from, n1, n2,
            perm, iperm, &current, sorter);
    FUNCTION_FAILURE_HANDLE( err, bandw_number_adj_vert, bandw_gipost);

    /* for all levels: assign the adjecent vertices of the lower level
       first with numbers, then the remaining of the higher level */
    /*   MSG_PRINT(" test:  start numbering from=%d  to=%d  dir=%+d   " */
    /*   "n1=%d  current=%d\n", from, to, dir, n1, current); */
    for (lvl=from+dir; dir*lvl<=dir*to; lvl+=dir)
    {
        /* lower -> higher */
        err=bandw_number_adj_vert( ad, Lv, hintv, lvl-dir, lvl, n1, n2,
                perm, iperm, &current, sorter);
        FUNCTION_FAILURE_HANDLE( err, bandw_number_adj_vert, bandw_gipost);
        /* higher -> higher */
        err=bandw_number_adj_vert( ad, Lv, hintv, lvl, lvl, n1, n2,
                perm, iperm, &current, sorter);
        FUNCTION_FAILURE_HANDLE( err, bandw_number_adj_vert, bandw_gipost);
        /*       MSG_PRINT(" test:     lvl=%d   current=%d\n", */
        /*       lvl, current); */
    }



    /*   MSG_PRINT(" test: gipost done,  n2=%d   ?=  current=%d\n", n2, current); */

#ifdef MESS_DEBUG
    if (current!=n2)
    {
        MSG_ERROR("bandw_gipost: failed to assign all vertices \n");
        for (i=n1; i<n2; i++)
            MSG_INFO(" vertex %3d  assigned  %3d, level=%3d\n",
                    (int) i, (int) perm[i], (int) hintv[i]);
        /*if (perm[i]==-1)
          fprintf(stderr," vertex %d  not assigned, level=%d\n",
          i, hintv[i]);*/
        return FAIL;
    }
#endif


    /* clean up */
    mess_free(sorter);
    int_list_free(Lv);
    mess_free(Lv);
    mess_free(hintv);

    /* compute the bandwidth */
    *bandw=0;
    for (i=n1; i<n2; i++)
        for (j=ad->ix[i]; j<ad->ix[i+1]; j++)
        {
#ifdef MESS64
            if (labs(perm[i]-perm[ad->el[j]]) > *bandw)
                *bandw=labs(perm[i]-perm[ad->el[j]]);
#else
            if (abs(perm[i]-perm[ad->el[j]]) > *bandw)
                *bandw=abs(perm[i]-perm[ad->el[j]]);
#endif
        }

    return SUCCESS;
}


/**
 * @internal
 * @brief Define level structure for a given subgraph.
 * @param[in] root      input root of level structure
 * @param[in] n1        input index of first node of connected subgraph
 * @param[in] n2        input index of last node \f$+1 \f$ of connected subgraph
 * @param[in] ad        input adjacency list of whole graph
 * @param[out] L        output level structure
 * @param[in,out] hint  input/output integer help vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref bandw_level_struct function defines the level structure rooted at root for a given subgraph.
 *
 * \f$ L \f$ has to be an empty int_list structure. \f$ L \f$ is initialized and its memory is allocated in
 * this function. \n
 * \f$ hint \f$ has to be provided by the calling routine and has a length of \f$ n2-n1 \f$.
 *
 * @attention Interal use only.
 */
static int bandw_level_struct( mess_int_t root, mess_int_t n1, mess_int_t n2, struct int_list *ad,
        struct int_list *L, mess_int_t *hint)
{
    MSG_FNAME(__func__);
    int  err;
    mess_int_t i, j, nr, found;


    nr=n2-n1;

    err=int_list_alloc(L, nr, nr);
    FUNCTION_FAILURE_HANDLE(err, int_list_alloc, bandw_level_struct);

    for (i=0; i<nr; i++) hint[i]=-1;

    L->ix[0]=0;
    L->el[0]=root;
    hint[root-n1]=0;
    L->ixn=1;
    L->ix[1]=1;
    if (nr==1)
    {
        MSG_ERROR( "bandw_level_struct: "
                "called for subgraph of only one node\n");
        return FAIL;
    }


    found=1; /* to start the search for more nodes */
    /* as @ref mess_int_t as the last level produced vertices in the new level, create
       yet a new level */
    while(found==1)
    {
        mess_int_t lvl;
        found=0;
        lvl=L->ixn;
        L->ix[lvl+1]=L->ix[lvl];
        /* for all nodes of the last level */
        for (i=L->ix[lvl -1]; i<L->ix[lvl]; i++)
        {
            mess_int_t vx, avx;
            /* mark all adjecent vertices which are assined to a level
               as bemess_int_ting to the new level */
            vx=L->el[i];
            for (j=ad->ix[vx]; j<ad->ix[vx+1]; j++)
            {
                avx=ad->el[j];
                if (hint[avx-n1]==-1)
                {
                    hint[avx-n1]=lvl;
                    L->el[L->ix[lvl+1]]=avx;
                    L->ix[lvl+1]++;
                    found=1;
                }
            }
        }
        L->ixn++;
    }
    /* correct the last increment, since nothing was found  */
    L->ixn--;

    /* if the last level is empty, reduce the number of levels again */
    if (L->ix[L->ixn]==L->ix[L->ixn+1]) L->ixn--;

    /* MSG_PRINT("  test: levels=%3d, el in last level=%3d\n", L->ixn,
       L->ix[L->ixn+1]-L->ix[L->ixn]); */

#ifdef MESS_DEBUG
    for (i=0; i<nr; i++)
        if (hint[i]==-1)
        {
            MSG_ERROR("bandw_level_struct: node with no level found\n"
                    "node=%d, n1=%d, n2=%d, root=%d \n",
                    (int) (i+n1), (int) n1, (int) n2, (int) root);

            return FAIL;
        }
#endif


    return SUCCESS;
}


/**
 * @internal
 * @brief Number adjacent vertices.
 * @param[in] ad            input adjacency list of whole graph
 * @param[in] Lv            input level structure
 * @param[in] level         input vector specifying level for each node
 * @param[in] l1            input start level
 * @param[in] l2            input offset of lowest node number in considered subgraph
 * @param[in] n1            input index of first node of connected subgraph
 * @param[in] n2            input index of last node \f$+1 \f$ of connected subgraph
 * @param[in,out] perm      input/output integer vector of permutation
 * @param[in,out] iperm     input/output inverse of \f$ perm \f$
 * @param[in,out] current   input/output current highest assigned number \f$ +1 \f$
 * @param[out] sorter       output integer vector
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref bandw_number_adj_vert function starts with a vertex in level \f$ l1 \f$ of lowest assigned number,
 * assigns all adjecent nodes in level \f$ l2 \f$ (which have not been
 * assigned a number yet) a number in order of increasing degree,
 * moves to the next until no further are left.
 *
 * The permutation vector \f$ perm \f$ has a length of the number of vertices, where for every \f$ i \f$
 * with \f$ perm[i]!=-1 \f$ \f$ perm[i]\f$ is considered already given a number.\n
 * \f$ current \f$ is updated every time a vertex is assigned a number.
 *
 * @attention Internal use only.
 */
static int bandw_number_adj_vert(struct int_list *ad, struct int_list *Lv,
        mess_int_t *level, mess_int_t l1, mess_int_t l2, mess_int_t n1, mess_int_t n2,
        mess_int_t *perm, mess_int_t *iperm, mess_int_t *current,
        mess_int_t *sorter ){
#ifdef MESS_DEBUG
    MSG_FNAME(__func__);
#endif
    mess_int_t i, j, sl;

    mess_int_t lowest_l1, found_new;

    /* find lowest number in level l1 */
    i=*current-1;
    while ((level[iperm[i]-n1]==l1)&&(i>n1)) i--;
    if (level[iperm[i]-n1]!=l1) i++;
    lowest_l1=i;

#ifdef MESS_DEBUG
    if (level[iperm[lowest_l1]-n1]!=l1)
    {
        MSG_ERROR("bandw_number_adj_vert: something wrong, nothing "
                "numbered in lower level!\n");
        return FAIL;
    }
#endif

    /* MSG_PRINT(" test: lowest_l1= %3d   l1=%3d   l2=%3d\n", lowest_l1,
       l1, l2);  */

    /* as @ref mess_int_t as there apear unnumbered nodes in the level... */
    found_new=1;
    while(found_new!=0)
    {
        found_new=0;

        /* check if all adjecent nodes to lowest_l1 are numbered yet: */
        if ((iperm[lowest_l1]>=0)&&(level[iperm[lowest_l1]-n1]==l1))
        {
            i=iperm[lowest_l1];
            /* MSG_PRINT(" test: check connected to %d:", i); */
            /* set sorter length (sl) to 0 = number of found unnumbered
               adjecent vertices */
            sl=0;
            /* for all adjecent nodes */
            for (j=ad->ix[i]; j<ad->ix[i+1]; j++)
            {
                mess_int_t adj_vx;
                adj_vx=ad->el[j];
                /* MSG_PRINT("   test: adj_vx=%3d    level[adj]=%3d   perm[adj]=%3d\n",
                   adj_vx, level[adj_vx-n1], perm[adj_vx]);  */
                /* if in level l2 and not numbered */
                if ((level[adj_vx-n1]==l2)&&(perm[adj_vx]==-1))
                {
                    /* add to the sorter list */
                    sorter[sl*2+1]=adj_vx;
                    /* sorting is by degree of the vertices */
                    sorter[sl*2  ]=ad->ix[adj_vx+1]-ad->ix[adj_vx];

                    found_new++;
                    sl++;
                }
            }
            /* MSG_PRINT("  %3d are connected\n",sl);  */
            /* if vertices are found, sort them by degree, number them
               in increasing order */
            if (sl>0)
            {
                qsort( sorter, sl, 2*sizeof(mess_int_t), comp_mess_int_t);
                for (j=0; j<sl; j++)
                {
                    mess_int_t tv;
                    tv=sorter[j*2+1];
                    /* MSG_PRINT(" test: sorter[%3d]= [ %3d, %3d ]   current=%3d\n",
                       j, sorter[j*2], sorter[j*2+1], *current);  */
                    perm[tv]=*current;
                    iperm[*current]=tv;
                    (*current)++;

#ifdef MESS_DEBUG
                    if (*current>n2)
                    {
                        MSG_ERROR("bandw_number_adj_vert: "
                                "something wrong, try to number over "
                                "limit: current=%d   n2=%d\n",
                                (int) *current, (int) n2);
                        return FAIL;
                    }
#else
                    if (*current==n2) return SUCCESS;
#endif
                }
            }

            /* increase the lowest_l1 counter, see if it is still l1 */
            if (lowest_l1<n2-1)
            {
                lowest_l1++;
                if ((iperm[lowest_l1]>=0)&&(level[iperm[lowest_l1]-n1]==l1))
                {
                    /* still l1 ==> make sure its adjacent vertices are
                       numbered first */
                    if (found_new==0) found_new=1;
                }
            }
        } /* end numbering of adjacent vertices to numbered ones */

        /* if the levels are identical and no adjacent vertices to
           numbered ones are found, check if there are unnumbered ones
           left, if so assign the one of lowest degree with a number */
        if ((l1==l2)&&(found_new==0))
        {
            mess_int_t tv, min_deg, candidate =0;
            min_deg=-1;
            /* check all vertices in this level */
            for (i=Lv->ix[l1]; i<Lv->ix[l1+1]; i++)
            {
                tv=Lv->el[i];
                if (perm[tv]==-1)
                {
                    /* see if it is a condidate for min degree */
                    if (min_deg==-1)
                    {
                        /* the first one -> is candidate */
                        candidate=tv;
                        min_deg=ad->ix[tv+1] - ad->ix[tv];
                    }
                    else if (ad->ix[tv+1] - ad->ix[tv] < min_deg)
                    {
                        /* not the first but lower degree */
                        candidate=tv;
                        min_deg=ad->ix[tv+1] - ad->ix[tv];
                    }
                }
            } /* end all vertices of the level */
            /* if a new node of minimal degree was found, number it */
            if (min_deg!=-1)
            {
                tv=candidate;
                perm[tv]=*current;
                iperm[*current]=tv;
                (*current)++;
                lowest_l1=tv;
                found_new=1;
#ifdef MESS_DEBUG
                if (*current>n2)
                {
                    MSG_ERROR("bandw_number_adj_vert: "
                            "something wrong, try to number over "
                            "limit: current=%d   n2=%d\n",
                            (int) *current, (int) n2);
                    return FAIL;
                }
#else
                if (*current==n2) return SUCCESS;
#endif
            }
        } /* end check for unnumbered if l1==l2 */
    } /* end check for nodes to be numbered */

    /* done */
    return SUCCESS;
}

