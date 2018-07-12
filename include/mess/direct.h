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
 * @file include/mess/direct.h
 * @brief Direct linear systems solvers and single-pattern multi-value LU.
 * @author @koehlerm
 */

#ifndef DIRECT_H_
#define DIRECT_H_

#include "mess/config.h"
#include "mess/matrix.h"
#include "mess/vector.h"
#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif


    /** @addtogroup direct_tri
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Triangular Solvers
     *-----------------------------------------------------------------------------*/
    int mess_solver_usolve(mess_matrix U, mess_vector y);
    int mess_solver_usolvem(mess_matrix U, mess_matrix Y);
    int mess_solver_usolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, double *y);
    int mess_solver_usolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *y);
    int mess_solver_usolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *y);

    int mess_solver_lsolve(mess_matrix L, mess_vector y);
    int mess_solver_lsolvem(mess_matrix L, mess_matrix Y);
    int mess_solver_lsolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr, double *y);
    int mess_solver_lsolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr, mess_double_cpx_t *y);
    int mess_solver_lsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr, mess_double_cpx_t *y);

    int mess_solver_lusolvem_kernelcsr_complex (mess_int_t dim, mess_double_cpx_t *lval, mess_int_t *lcolptr, mess_int_t *lrowptr, mess_double_cpx_t *uval, mess_int_t *ucolptr, mess_int_t *urowptr, mess_int_t *p, mess_int_t *q, mess_matrix b, mess_matrix x);


    // transpose solver
    int mess_solver_ltsolve(mess_matrix L, mess_vector y);
    int mess_solver_ltsolvem(mess_matrix L, mess_matrix Y);
    int mess_solver_lhsolve(mess_matrix L, mess_vector y);
    int mess_solver_lhsolvem(mess_matrix L, mess_matrix Y);
    int mess_solver_ltsolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr , double *y);
    int mess_solver_ltsolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y);
    int mess_solver_utsolve_kernelcsr_real(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, double *x);
    int mess_solver_utsolve_kernelcsr_real_complex(mess_int_t dim, double *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *x);
    int mess_solver_ltsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y);
    int mess_solver_lhsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t* colptr, mess_int_t *rowptr , mess_double_cpx_t *y);
    int mess_solver_utsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *x);
    int mess_solver_uhsolve_kernelcsr_complex(mess_int_t dim, mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr, mess_double_cpx_t *x);

    int mess_solver_utsolve(mess_matrix U, mess_vector y);
    int mess_solver_uhsolve(mess_matrix U, mess_vector y);
    int mess_solver_utsolvem(mess_matrix U, mess_matrix y);
    int mess_solver_uhsolvem(mess_matrix U, mess_matrix y);

    /** @}  */


    /** @addtogroup direct_interface
     * @{ */
    /**
     * @brief Representation of a generic direct solver.
     *
     * The @ref mess_direct_st structure represents a generic direct solver. \n
     * Before it can be used it is necessary to initialize the object using \ref mess_direct_init and to
     * clean it up after usage by calling \ref mess_direct_clear. \n
     * A full list of direct solvers is available in the \ref direct_solvers section.
     * Instead of using this structure directly please use the \ref mess_direct type definition because this
     * hides the pointer to this structure and it is used by all functions. Normally it is not necessary
     * to access components of this structure directly. All operations and actions are performed by functions
     * or macros of this section.
     */
    struct mess_direct_st {
        mess_datatype_t data_type;                                                          /**< Data type of the solver */
        mess_int_t rows;                                                                    /**< Number of rows */
        mess_int_t cols;                                                                    /**< Number of columns */
        void *data;                                                                         /**< Pointer to the data for the solver */
        int (*solve) (void *data, mess_vector b, mess_vector x);                            /**< Function to solve \f$ Ax=b   \f$ */
        int (*solvet) (void *data, mess_vector b, mess_vector x);                           /**< Function to solve \f$ A^Tx=b \f$ */
        int (*solvem) (void *data, mess_matrix b, mess_matrix x);                           /**< Function to solve \f$ AX=B   \f$ */
        int (*solvemt) (void *data, mess_matrix b, mess_matrix x);                          /**< Function to solve \f$ A^TX=B \f$ */
        int (*solveh) (void *data, mess_vector b, mess_vector x);                           /**< Function to solve \f$ A^Hx=b \f$ */
        int (*solvemh) (void *data, mess_matrix b, mess_matrix x);                          /**< Function to solve \f$ A^HX=B \f$ */
        int (*getL)     (void *data, mess_matrix L);                                        /**< Function to get a copy of the \f$ L \f$ factor */
        int (*getU) (void *data, mess_matrix U);                                            /**< Function to get a copy of the \f$ U \f$ factor */
        int (*getpermp) (void *data, mess_int_t *p);                                        /**< Function to get a copy of the row permuatation */
        int (*getpermq) (void *data, mess_int_t *q);                                        /**< Function to get a copy of the column  permutation */
        int (*getscalerow) (void *data, mess_vector r);                                     /**< Function to get a copy of a row scaling numbers */
        int (*getscalecol) (void *data, mess_vector c);                                     /**< Function to get a copy of a col scaling numbers */
        int (*det)  (void *data, double *m, double *e);                                     /**< Function to get the real determinant */
        int (*detc) (void *data, double *mr, double *mi, double *e);                        /**< Function to get the complex determinant */
        int (*inverse)  (void *data, mess_matrix inv);                                      /**< Function to compute the inverse */
        int (*clear)  (void *myself);                                                       /**< Clean up function*/
        char *name;                                                                         /**< Name of the solver */
    };

    /**
     * @brief Type definition to hide the pointer on the \ref mess_direct_st structure.
     *
     * The @ref mess_direct type definition hides one pointer on the \ref mess_direct_st structure. \n
     * This type definition is the way how the direct solver interface should be used.
     **/
    typedef struct mess_direct_st *  mess_direct;

    /*-----------------------------------------------------------------------------
     *  A set of macros for easy handling of @ref mess_direct properties
     *-----------------------------------------------------------------------------*/
    /** Check if the solver provide a solve function.  */
#define MESS_HAVE_SOLVE(solver)     ((solver)->solve != NULL)
    /** Check if the solver provide a solvet function.  */
#define MESS_HAVE_SOLVET(solver)    ((solver)->solvet != NULL)
    /** Check if the solver provide a solveh function.  */
#define MESS_HAVE_SOLVEH(solver)    ((solver)->solveh != NULL)
    /** Check if the solver provide a solvem function.  */
#define MESS_HAVE_SOLVEM(solver)    ((solver)->solvem != NULL)
    /** Check if the solver provide a solvemt function.  */
#define MESS_HAVE_SOLVEMT(solver)   ((solver)->solvemt != NULL)
    /** Check if the solver provide a solvemh function.  */
#define MESS_HAVE_SOLVEMH(solver)   ((solver)->solvemh != NULL)
    /** Check if the solver provide a getL function.  */
#define MESS_HAVE_GETL(solver)      ((solver)->getL != NULL)
    /** Check if the solver provide a getU function.  */
#define MESS_HAVE_GETU(solver)      ((solver)->getU != NULL)
    /** Check if the solver provide a getpermp function.  */
#define MESS_HAVE_GETPERMP(solver)  ((solver)->getpermp != NULL)
    /** Check if the solver provide a getpermq function.  */
#define MESS_HAVE_GETPERMQ(solver)  ((solver)->getpermq != NULL)
    /** Check if the solver provide a getscalerow function.  */
#define MESS_HAVE_GETSCALEROW(solver)   ((solver)->getscalerow != NULL)
    /** Check if the solver provide a getscalecol function.  */
#define MESS_HAVE_GETSCALECOL(solver)   ((solver)->getscalecol != NULL)
    /** Check if the solver provide a det function.  */
#define MESS_HAVE_DET(solver)       ((solver)->det != NULL)
    /** Check if the solver provide a detc function.  */
#define MESS_HAVE_DETC(solver)      ((solver)->det != NULL)
    /** Check if the solver provide a inverse function.  */
#define MESS_HAVE_INVERSE(solver)       ((solver)->inverse != NULL)


    int mess_direct_init(mess_direct *direct);
    int mess_direct_clear(mess_direct *direct);
    int mess_direct_solve(mess_operation_t op, mess_direct solver, mess_vector b, mess_vector x);
    int mess_direct_solvem(mess_operation_t op, mess_direct solver, mess_matrix b, mess_matrix x);
    int mess_direct_determinant(mess_direct solver, double *m, double *e);
    int mess_direct_determinantc(mess_direct solver, double *mr, double *mi, double *e);
    int mess_direct_getL(mess_direct solver, mess_matrix L );
    int mess_direct_getU(mess_direct solver, mess_matrix U );
    int mess_direct_getpermp(mess_direct solver, mess_int_t *p );
    int mess_direct_getpermq(mess_direct solver, mess_int_t *q );
    int mess_direct_getscalerow(mess_direct solver, mess_vector r );
    int mess_direct_getscalecol(mess_direct solver, mess_vector c );
    int mess_direct_inverse(mess_direct solver, mess_matrix inv);
    int mess_matrix_backslash(mess_operation_t op, mess_matrix A, mess_vector b, mess_vector x) ;
    int mess_matrix_backslashm(mess_operation_t op, mess_matrix A, mess_matrix b, mess_matrix x) ;
    int mess_matrix_pinv(mess_matrix A, mess_matrix Pinv );
    /** @} */

    /** @addtogroup direct_solvers
     * @{*/
    /*-----------------------------------------------------------------------------
     * Normal Linear Solvers
     *-----------------------------------------------------------------------------*/
    /**
     * @brief Enumeration type to select the behaviour of \ref mess_direct_lu.
     *
     * The @ref mess_direct_lupackage_t type is an enumeration to select the behavior of
     * \ref mess_direct_lu. \n
     * The value \ref MESS_DIRECT_DEFAULT_LU selects the default
     * behavior. That means depending of the storage type of the @ref mess_matrix object a suitable
     * and available solver will be selected. All other enumeration items fix it to
     * the desired solver.
     */
    typedef enum {
        MESS_DIRECT_DEFAULT_LU      = 0     /**< Select the automatic solver detection. */
            ,MESS_DIRECT_SPARSE_LU      = 1     /**< Select the internal sparse solver. */
            ,MESS_DIRECT_LAPACK_LU      = 2     /**< Select the @lapack dense solver. */
            ,MESS_DIRECT_UMFPACK_LU     = 3     /**< Select the @umfpack sparse solver. */
            ,MESS_DIRECT_SUPERLU_LU     = 4     /**< Select the @superlu sparse solver. */
            ,MESS_DIRECT_CSPARSE_LU     = 5     /**< Select the @csparse sparse solver. */
            ,MESS_DIRECT_BANDED_LU      = 6     /**< Select the banded sparse solver. */
            ,MESS_DIRECT_MKLPARDISO_LU  = 7     /**< Select the @ref mess_direct_create_mklpardiso solver. */
    } mess_direct_lupackage_t;
    int mess_direct_lu(mess_matrix matrix, mess_direct direct);
    int mess_direct_lu_select(mess_direct_lupackage_t lupackage);


    /**
     * @brief Enumeration type to select the behaviour of \ref mess_direct_chol.
     *
     * The @ref mess_direct_cholpackage_t type is an enumeration to select the behavior of
     * \ref mess_direct_chol. \n
     * The value \ref MESS_DIRECT_DEFAULT_CHOLESKY selects the default
     * behavior. That means depending on the @ref mess_storage_t type of the @ref mess_matrix a suitable
     * and available solver will be selected. All other enumeration items fix it to
     * the desired solver.
     */
    typedef enum {
        MESS_DIRECT_DEFAULT_CHOLESKY    = 0     /**< Select the automatic solver detection. */
            , MESS_DIRECT_LAPACK_CHOLESKY   = 1     /**< Select the @lapack dense solver. */
            , MESS_DIRECT_CSPARSE_CHOLESKY  = 2     /**< Select the @csparse sparse solver. */
            , MESS_DIRECT_CHOLMOD_CHOLESKY  = 3     /**< Select the @cholmod sparse solver. */
    } mess_direct_cholpackage_t;
    int mess_direct_chol(mess_matrix matrix, mess_direct direct);
    int mess_direct_chol_select(mess_direct_cholpackage_t cholpackage);


#ifdef MESS_HAVE_UMFPACK
    int mess_direct_create_umfpack(mess_matrix matrix, mess_direct direct);
    int mess_direct_create_umfpack_control(mess_matrix matrix, mess_direct direct, double* Control);
#endif
#ifdef MESS_HAVE_SUPERLU
    int mess_direct_create_superlu(mess_matrix matrix, mess_direct direct);
#endif
#ifdef MESS_HAVE_CSPARSE
    int mess_direct_create_csparse_lu(mess_matrix matrix, mess_direct direct);
    int mess_direct_create_csparse_cholesky(mess_matrix matrix, mess_direct direct);
#endif
#ifdef MESS_HAVE_CHOLMOD
    int mess_direct_create_cholmod_cholesky(mess_matrix matrix, mess_direct direct);
#endif

#ifdef MESS_HAVE_MKLPARDISO
    int mess_direct_create_mklpardiso_control(mess_matrix matrix, mess_direct direct, mess_int_t*iparm, mess_int_t mtype);
    int mess_direct_create_mklpardiso(mess_matrix matrix, mess_direct direct);
#endif

    int mess_direct_create_banded(mess_matrix matrix, mess_direct direct);
    int mess_direct_create_bicgstab(mess_matrix matrix, mess_direct solver);
    int mess_direct_create_lapack_lu(mess_matrix matrix, mess_direct solver);
    int mess_direct_create_cholesky(mess_matrix matrix, mess_direct solver);
    int mess_direct_create_lapack_qr (mess_matrix matrix, mess_direct solver);
    int mess_direct_create_sparse_lu(mess_matrix A, mess_direct sol);
    int mess_direct_create_hessenberg_lu(mess_matrix matrix, mess_direct solver);

    /** @}  */

    /** @addtogroup direct_misc
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Misc
     *-----------------------------------------------------------------------------*/
    int mess_direct_splsolve (mess_matrix G, mess_matrix B, mess_int_t k, mess_int_t *top, mess_int_t *xi, double *x, const mess_int_t *pinv);
    int mess_direct_splsolvec (mess_matrix G, mess_matrix B, mess_int_t k, mess_int_t *top, mess_int_t *xi, mess_double_cpx_t *x, const mess_int_t *pinv);
    int mess_decomp_lureuse_kernelcsr(mess_int_t rows, mess_int_t cols,                                 // basic
            double *values, mess_int_t *colptr, mess_int_t *rowptr,         // input matrix
            mess_int_t *lcolptr, mess_int_t *lrowptr,                       // L pattern
            mess_int_t *ucolptr, mess_int_t *urowptr,                       // U pattern
            mess_int_t *p, mess_int_t *q,                                   // permuatations
            double *lvalues, double *uvalues);
    int mess_decomp_lureuse_kernelcsr_complex(mess_int_t rows, mess_int_t cols,
            mess_double_cpx_t *values, mess_int_t *colptr, mess_int_t *rowptr,         // input matrix
            mess_int_t *lcolptr, mess_int_t *lrowptr,                       // L pattern
            mess_int_t *ucolptr, mess_int_t *urowptr,                       // U pattern
            mess_int_t *p, mess_int_t *q,                                   // permuatations
            mess_double_cpx_t *lvalues, mess_double_cpx_t *uvalues);

    int mess_decomp_lureuse_kernelcsr_flopcount(mess_int_t rows, mess_int_t cols,   // basic
            double *values, mess_int_t *colptr, mess_int_t *rowptr,         // input matrix
            mess_int_t *lcolptr, mess_int_t *lrowptr,                       // L pattern
            mess_int_t *ucolptr, mess_int_t *urowptr,                       // U pattern
            mess_int_t *p, mess_int_t *q,                                   // permuatations
            double *lvalues, double *uvalues, mess_int_t *flops);


    /** @}  */

    /** @addtogroup direct_meq
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Matrix Equations
     *-----------------------------------------------------------------------------*/
    int mess_direct_create_generalized_lyapunovchol(mess_matrix A, mess_matrix E, mess_direct solver);
    int mess_direct_create_generalized_lyapunov(mess_matrix A, mess_matrix E, mess_direct solver);
    int mess_direct_create_generalized_stein(mess_matrix A, mess_matrix E, mess_direct solver);
    int mess_direct_care(mess_matrix A, mess_matrix E, mess_matrix Q, mess_matrix G, mess_matrix X);
    int mess_direct_care_sign(mess_matrix A, mess_matrix Q, mess_matrix G, mess_matrix X);

    int mess_direct_create_sylvester_sparsedense(mess_matrix A, mess_matrix F, mess_matrix E, mess_matrix H, mess_direct solver );
    int mess_direct_create_sylvester_dense(mess_matrix A, mess_matrix F, mess_matrix E, mess_matrix H, mess_direct solver );
    int mess_direct_create_generalized_lyapunov(mess_matrix A, mess_matrix E, mess_direct solver);
    int mess_direct_create_generalized_stein(mess_matrix A, mess_matrix E, mess_direct solver);

    /** @} */

    /** @addtogroup  direct_decompose
     *
     * The decompositions in this section does not rely on the \ref direct_interface because they are not used
     * in the direct solver context. \n
     * This are for example the column pivoted QR decomposition or computing
     * the Cholesky factor of a may be singluar matrix.
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Other Decompositions
     *-----------------------------------------------------------------------------*/
    int mess_direct_dgeqp3(mess_matrix A, mess_matrix Q, mess_matrix R, mess_int_t *perm);
    int mess_direct_cholfactor(mess_matrix ZZ, mess_matrix Z);

    /** @}  */


    /** @addtogroup direct_res
     *
     * This section provides residual computations for many linear matrix equations and ordinary
     * linear systems. \n
     * Please keep in mind that all these operation work in dense arithmetics. That
     * means in the case of large-scale matrix equation they might fail. In this case we refer to the
     * techniques presented in the \ref lrcfadi section.
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Residuals.
     *-----------------------------------------------------------------------------*/
    int mess_direct_res2(mess_operation_t op, mess_matrix A, mess_vector x, mess_vector b, double *resid, double * relresid);
    int mess_direct_res2m(mess_operation_t op, mess_matrix A, mess_matrix x, mess_matrix b, double *resid, double * relresid);

    int mess_direct_sylvester_res2(mess_operation_t op, mess_matrix A, mess_matrix H, mess_matrix M, mess_matrix X, double *res);
    int mess_direct_sylvestersg_res2(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix H, mess_matrix M, mess_matrix X, double *res);
    int mess_direct_generalized_sylvester_res2(mess_operation_t op, mess_matrix A, mess_matrix F, mess_matrix E, mess_matrix H, mess_matrix M, mess_matrix X, double *res);

    int mess_direct_lyapunov_res2(mess_operation_t op, mess_matrix A, mess_matrix M, mess_matrix X, double *res);
    int mess_direct_generalized_lyapunov_res2(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix M , mess_matrix X, double *res2);

    int mess_direct_stein_res2(mess_operation_t op, mess_matrix A, mess_matrix M , mess_matrix X, double *res2);
    int mess_direct_generalized_stein_res2(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix M , mess_matrix X, double *res2);
    int mess_direct_care_res2(mess_matrix A, mess_matrix E, mess_matrix Q,mess_matrix G,  mess_matrix X, double *res2);

    /** @} */


    /** @addtogroup multisolvers_interface
     * @{ */
    /**
     * @brief Representation of a generic direct solver.
     *
     * The concept of @ref mess_multidirect solvers in @mess represents a set of direct solvers of the structure
     * \f[ (t_i A + p_i E) x = b,  \f]
     * where \f$ t_i \f$ and \f$ p_i \f$ are a set of shift parameters. \n
     * The multi-direct solvers try to reuse as much
     * information as possible to squeeze down the memory footprint. \n
     * This structure represents a generic interface for this class of solvers. The structure is used as a
     * pointer and so the \ref mess_multidirect type definition
     * hides one pointer if this structure is used. \n
     * Accessing components of this structure is normally not necessary
     * and all operations are covered by functions and macros of this section.
     *
     * A list of possible @ref mess_multidirect solvers is available in the \ref multisolvers_solvers section.
     */
    struct mess_multidirect_st {
        mess_datatype_t data_type;                                                  /**< Data type of the solver */
        mess_int_t indx;                                                            /**< Number of different values */
        mess_int_t rows;                                                            /**< Number of rows */
        mess_int_t cols;                                                            /**< Number of columns */
        void *data;                                                                 /**< Pointer to the solver data */
        int (*solve)   (void *data,mess_int_t ind, mess_vector b, mess_vector x);   /**< Function to solve \f$ A(i)x=b  \f$  */
        int (*solvet)  (void *data,mess_int_t ind, mess_vector b, mess_vector x);   /**< Function to solve \f$ A(i)^Tx=b \f$ */
        int (*solveh)  (void *data,mess_int_t ind, mess_vector b, mess_vector x);   /**< Function to solve \f$ A(i)^Hx=b \f$ */
        int (*solvem)  (void *data,mess_int_t ind, mess_matrix b, mess_matrix x);   /**< Function to solve \f$ A(i)X=B \f$ */
        int (*solvemt) (void *data,mess_int_t ind, mess_matrix b, mess_matrix x);   /**< Function to solve \f$ A(i)^TX=B \f$ */
        int (*solvemh) (void *data,mess_int_t ind, mess_matrix b, mess_matrix x);   /**< Function to solve \f$ A(i)^HX=B \f$ */
        int (*clear)  (void *myself);                                               /**< Function to clean up */
        size_t (*memsize) (void *data);                                             /**< Function to get the size of the data */
        unsigned short (*getdatatype)(void *myself, mess_int_t ind);                /**< Get the data type of a specified solver */
        int (*getL) (void *data, mess_int_t ind, mess_matrix L);                    /**< Get the L factor of a decomposition */
        int (*getU) (void *data, mess_int_t ind, mess_matrix U);                    /**< Get the U factor of a decomposition */
        int (*getp) (void *data, mess_int_t ind, mess_int_t *p);                    /**< Get the row permutation of a decomposition */
        int (*getq) (void *data, mess_int_t ind, mess_int_t *q);                    /**< Get the column permutation of a decomposition */
        int (*getscalerow) (void *data, mess_int_t ind, mess_vector r);             /**< Function to get a copy of a row scaling numbers */
        int (*getscalecol) (void *data, mess_int_t ind, mess_vector c);             /**< Function to get a copy of a col scaling numbers */
        char *name;                                                                 /**< Name of the solver */
    };

    /**
     * @brief Type definition to hide one pointer of the \ref mess_multidirect_st structure.
     *
     * The @ref mess_multidirect type definition hides one pointer to the \ref mess_multidirect_st structure. \n
     * This is the common way how the @ref mess_multidirect  solver structure should be used. All function of the
     * interface take the type definition instead of the pointer to the structure as input parameter.
     */
    typedef struct mess_multidirect_st *mess_multidirect;

    /**
     * @brief Type definition to select the desired multisolver via \ref mess_multidirect_select.
     *
     * The @ref mess_multidirect_t enumeration contains the possible multi-solvers for
     * \f[ (o_i A + p_i E) x = b.\f]
     * The enumeration can be used to select a multisolver via \ref mess_multidirect_select or inside
     * the ADI processes to select the shifted solver during the iteration. This selection is done in the
     * \ref mess_options structure.
     * */
    typedef enum {
        MESS_MULTIDIRECT_SPARSE_LU      = 0,    /**< Select the internal Single-Pattern Multi-Value multisolver. \see mess_multidirect_create_sparse_lu */
        MESS_MULTIDIRECT_UMFPACK_LU     = 1,    /**< Select the @umfpack based multissolver.  \see mess_multidirect_create_umfpack */
    } mess_multidirect_t;

    int mess_multidirect_init( mess_multidirect *direct);
    int mess_multidirect_clear( mess_multidirect *direct);

    int mess_multidirect_solve(mess_operation_t op, mess_multidirect solver, mess_int_t ind, mess_vector b, mess_vector x);
    int mess_multidirect_solvem(mess_operation_t op, mess_multidirect solver,  mess_int_t ind,mess_matrix b, mess_matrix x);



    /** @} */

    /** @addtogroup multisolvers_solvers
     * @{ */
#ifdef MESS_HAVE_UMFPACK
    int mess_multidirect_create_umfpack(mess_matrix matrix, mess_vector shiftsl, mess_vector shiftsr , mess_multidirect mlu, mess_direct base, mess_matrix shiftmatrix);
#endif
    int mess_multidirect_create_sparse_lu(mess_matrix matrix, mess_vector shiftsl, mess_vector shiftsr , mess_multidirect mlu, mess_direct base, mess_matrix shiftmatrix) ;

    int mess_multidirect_create(mess_matrix matrix, mess_vector shiftsl, mess_vector shiftsr , mess_multidirect mlu, mess_direct base, mess_matrix shiftmatrix);
    int mess_multidirect_select( mess_multidirect_t solver);

    /** @} */


    /** @addtogroup direct_misc
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Parallel Stuff and structural analyze
     *-----------------------------------------------------------------------------*/
    // #warning Kommentieren oder verstecken
    typedef struct mess_direct_levelset_st {
        mess_int_t levels;
        mess_int_t *levelptr;
        mess_int_t *levelind;
    } mess_direct_levelset;
    /** @}  */

    /* Only in this file to avoid a forward reference on @ref mess_direct */
    /** @addtogroup matrix_norm
     * @{
     */
    int mess_matrix_norminvest(mess_direct decomp,double *nrm);
    int mess_matrix_condest(mess_matrix A, double *nrm);
    /** @} */

#ifdef __cplusplus
}
#endif

#endif /* DIRECT_H_ */

/** \}@ */
