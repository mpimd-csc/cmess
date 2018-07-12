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
 * @file include/mess/lrcf_adi.h
 * @brief LRCF-ADI related functions and data structures.
 * @author @koehlerm
 */

#ifndef LRCF_ADI_H_
#define LRCF_ADI_H_

#ifdef  __cplusplus
extern "C" {
#endif
    /** @addtogroup lrcfadi_generic
     *  @{
     */
    /**
     * @brief Enumeration type for the available parameter types.
     *
     * The @ref mess_parameter_t enumeration type is used to select the parameters inside the
     * ADI process. \n
     * At the moment we have implemented three different parameter types. The
     * parameters are computed using \ref mess_lrcfadi_parameter.
     *
     */
    typedef enum {
        MESS_LRCFADI_PARA_MINMAX        = 0,    /**< Compute complex parameters using the heuristic for the min-max problem. */
        MESS_LRCFADI_PARA_MINMAX_REAL   = 1,    /**< Compute real parameters using the heuristic for the min-max problem. */
        MESS_LRCFADI_PARA_WACHSPRESS    = 2,    /**< Compute Wachspress parameters.  */
        MESS_LRCFADI_PARA_ADAPTIVE_V    = 3,    /**< Compute adaptive shifts with \f$V\f$ as projection space. */
        MESS_LRCFADI_PARA_ADAPTIVE_Z    = 4     /**< Compute adaptive shifts with \f$Z\f$ as projection space. */
    } mess_parameter_t;



    /**
     * @brief Enumeration type for the available equation types.
     *
     * The @ref mess_equation_t enumeration type is used to define the type of the equation. Mostly used internally.
     *
     */
    typedef enum {
        MESS_EQN_NONE       = 0,            /**< No equation (only for initialization)*/
        MESS_EQN_LYAP       = 1,            /**< Lyapunov Equation */
        MESS_EQN_GLYAP      = 2,            /**< Generalized Lyapunov Equation*/
        MESS_EQN_RICCATI    = 3,            /**< Riccati Equation */
        MESS_EQN_GRICCATI   = 4             /**< Generalized Riccati Equation*/
    } mess_equation_t;

    /**
     * @brief Enumeration type for the available memory usage option in @ref mess_options.
     *
     * The @ref mess_memusage_t enumeration type is used to define the amount of memory the solvers of equation.
     * @sa mess_equation_std_lyap
     * @sa mess_equation_std_riccati
     */
    typedef enum {
        MESS_MEMORY_LOW     = 0,        /**< Lowest memory consuming variant.*/
        MESS_MEMORY_MID     = 1,        /**< Average memory consuming variant. */
        MESS_MEMORY_HIGH    = 2         /**< Highest memory using variant. */
    } mess_memusage_t;

    /**
     * @brief Enumeration type for the residual computation option in @ref mess_options.
     *
     * The @ref mess_residual_t enumeration type is used to define the algorithm for computation of the residual.
     * @sa mess_equation_std_lyap
     * @sa mess_equation_std_riccati
     */
    typedef enum {
        MESS_RESIDUAL_INDEFINITE    = 0,        /**< Residual via Residual Factor Updates (fastest variant).*/
        MESS_RESIDUAL_SPECTRAL      = 1,        /**< Compute Largest Eigenvalue iteratively.  */
    } mess_residual_t;


    /**
     * @internal
     * @brief Structure to pass information to the step function inside the LRCF iterations.
     *
     * The @ref mess_lrcfadi_step structure is used to pass information to a step function inside the ADI
     * iteration. \n
     * The step function, if set, is executed at the end of every ADI or Newton iteration. \n
     * This structure is passed to the step function to provide information about the current state
     * of the iteration. \n
     * The stop_user component can be used to provide a user supplied cancellation
     * criterion.
     * @attention Interal use only.
     */
    typedef struct __mess_lrcfadi_step {
        double res2;                    /**< Current \f$ 2 \f$-norm residual.  */
        double res2_change;             /**< Current \f$ 2 \f$-norm residual change. */
        double rel_change;              /**< Current Frobenius norm change in the factor.  */
        mess_matrix Z;                  /**< Current solution factor \f$ Z \f$. */
        mess_int_t it;                  /**< Current iteration number. */
        unsigned short stop_user;       /**< Output parameter, flag to cancel the iteration. */
    } mess_lrcfadi_step;


    /* Only a forward reference.  */
    struct mess_equation_st;

    /**
     * @brief Options for the LRCF-ADI algorithm family.
     *
     * The @ref mess_options_st structure is used to configure the LRCF-ADI based solvers. \n
     * It contains information about the cancellation criteria, maximum iteration numbers, the
     * shift parameter computation and similar options. It is used by all LRCF related solvers. \n
     * The initialization is done using \ref mess_options_init and it is freed after usage
     * by \ref mess_options_clear. An overview about the current setting can be printed to the
     * standard output using \ref mess_options_print.
     *
     * \sa mess_lradi
     * \sa mess_lrnm
     *
     **/
    typedef struct mess_options_st {
        /** Type of the equation. */
        mess_operation_t type;

        /* SHIFTS  */
        /** Precomputed shift parameters. */
        mess_vector adi_shifts_p;
        /** Type of shift parameters. */
        mess_parameter_t adi_shifts_paratype;
        /** Number of steps in the Arnoldi process with respect to \f$ A \f$ or \f$ (A-UW) \f$. */
        mess_int_t adi_shifts_arp_p;
        /** Number of steps in the Arnoldi process with respect to \f$ A^{-1} \f$ or \f$ (A-UW)^{-1} \f$. */
        mess_int_t adi_shifts_arp_m;
        /** Number of shifts. */
        mess_int_t adi_shifts_l0;
        /** Start vector for the Arnoldi process. */
        mess_vector adi_shifts_b0;
        /** Projection space given by the columns (only needed by @ref MESS_LRCFADI_PARA_ADAPTIVE_V). */
        mess_matrix proj_space;

        /* ADI PART   */
        /** Maximal number of iteration for the ADI process. */
        mess_int_t adi_maxit;
        /** relative \f$ 2 \f$-norm stopping tolerance in the ADI process. */
        double adi_res2_tol;
        /** relative \f$ 2 \f$-norm change tolerance in the ADI process. */
        double adi_res2c_tol;
        /**  Stopping criteria: Frobenius-norm relative change in the ADI process.  */
        double adi_rel_change_tol;
        /** True if you want some nice output during the iteration. */
        mess_int_t adi_output;
        /** New iteration step function, called at the end of every full ADI step. */
        void (*stepfunction)(void *aux, struct mess_equation_st *eqn, struct mess_options_st *opt, mess_lrcfadi_step *step);
        /** Auxillary data for the stepfunction. */
        void *stepfunction_aux;
        /** Set Method for Residual Computation. */
        mess_residual_t residual_method;

        /* NM Part  */
        /** Maximal number of iteration in the Newton method. */
        mess_int_t nm_maxit;
        /** Stopping criteria: \f$ 2 \f$-norm residual (Newton). */
        double nm_res2_tol;
        /**  Verbosity level for the LRNM. */
        mess_int_t  nm_output;
        /** Turn line search on/off. */
        mess_int_t nm_linesearch;
        /** Set a single shift strategy in LRCF-NM. */
        unsigned short nm_singleshifts;
        /** Galerkin projektion in LRCF-NM. */
        mess_int_t nm_gpStep;
        /** New iteration step function, called at the end of every full Newton step. */
        void (*nm_stepfunction)(void *aux, struct mess_equation_st *eqn, struct mess_equation_st *lyap, struct mess_options_st *opt, mess_lrcfadi_step *step);
        /** Auxillary data for the step function. */
        void *nm_stepfunction_aux;

        /** control the usage of memory, 0 indicates the potential lowest memory usage and slowest performance **/
        mess_memusage_t memory_usage;

        /* Internal   */
        /** Matrix \f$ K0 \f$ initial feedback for LRNM. */
        mess_matrix K0;

        /** Matrix \f$ W \f$  low rank lyapunov residual from LRADI. */
        mess_matrix W;

    } mess_options_st;

    typedef mess_options_st * mess_options;

    /**
     * @brief Structure/object to save the final status of a LRCF algorithm.
     *
     *
     */
    typedef struct mess_status_st {
        mess_vector     res2_norms;                 /**<  Absolute \f$ 2 \f$-norm in every step. */
        mess_vector     rel_changes;                /**<  Relative change of the Iterate in every step in \f$ F\f$-norm for lradi and in \f$ 2 \f$ norm for lrnm.*/
        double          res2_norm;                  /**<  Final absoulute \f$ 2 \f$-norm residual.  */
        double          res2_change;                /**<  Final relative \f$ 2 \f$-norm change of the absolute residual. */
        double          rel_change;                 /**<  Final relative change of \f$ Z \f$. */
        double          res2_0;                     /**<  \f$ 2\f$-norm residual of the right hand side.*/
        double          time_all;                   /**<  Overall time.  */
        double          time_adi;                   /**<  Time for ADI.  */
        mess_int_t      it;                         /**<  Number of iteration the algorithm took.  */
        mess_int_t      n_internal_status;          /**<  Internal status to number of iterations. */
        struct mess_status_st **internal_status;    /**<   Internal status at every iteration. */
        unsigned short  stop_res2;                  /**<  Number of iteration the algorithm stops with respect to relative \f$ 2 \f$-norm residual. */
        unsigned short  stop_res2c;                 /**<  Number of iteration the algorithm stops with respect to \f$ 2 \f$-norm relative change. */
        unsigned short  stop_rel;                   /**<  Number of iteration the algorithm stops with respect to relative change in \f$ Z \f$. */
        unsigned short  stop_user;                  /**<  Flag to cancel the iteration. */
        unsigned short  unstable;                   /**<  True if there occur unstable Ritz values. */
    }  mess_status_st;

    typedef mess_status_st * mess_status;

    /**
     * @brief Type definition for the definition of an equation.
     *
     * The @ref mess_status is a type definition for the definition of an equation.
     *
     *
     */
    typedef struct mess_equation_st *  mess_equation;

    /**
     * @brief Structure defining an equation as a set of function pointers.
     *
     * The @ref mess_equation_st structure contains all function pointers and data that
     * is necessary to define a Lyapunov or Riccati Equation. \n
     * All operations depending on the operator are defined as function pointers. \n
     * The  _generate function is used to initialize data structure once before _apply is called the first
     * time. \n
     * The  _clear function is called when the equation is destroyed.
     * */
    struct mess_equation_st {

        /**Pointer to child equation object.*/
        mess_equation child;

        /**Pointer to parent equation object.*/
        mess_equation parent;

        /** Type of the equation (Lyapunov, Riccati, etc..).  */
        mess_equation_t eqn_type;
        /** State space dimension of the equation.  */
        mess_int_t dim;
        /** Matrix \f$ B \f$ of the operator. */
        mess_matrix B;
        /** Matrix \f$ C \f$ of the operator  */
        mess_matrix C;

        /** Matrix \f$ K \f$ of ther operator \f$ (A-BK)\f$ or \f$ (A-CK)\f$ */
        mess_matrix K;

        /** Right hand side for the current Lyapunov Equation, if @c NULL B is used. */
        mess_matrix RHS;
        /** True if the right hand side needs to be cleared.  */
        mess_int_t clearRHS;
        /** True if \f$ B \f$ matrix needs to be cleared.  */
        mess_int_t clearB;
        /** True if \f$ C \f$ matrix  needs to be cleared.  */
        mess_int_t clearC;

        /**
         * @brief Structure defining \f$ Y = op(A)X \f$.
         *
         * Structure defining \f$ Y = op(A)X \f$ where \f$ X \f$ and \f$ Y \f$ are matrices.
         * */
        struct AX_st{
            /** Post process function if \f$ Y = AX \f$ needs a global initialization.  */
            int (*generate)(struct mess_equation_st * eqn);
            /** Apply the operator to a matrix. Compute \f$ Y=op(A)X \f$. */
            int (*apply)(struct mess_equation_st * eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
            /** Clear \f$ AX \f$ objects separately. */
            int (*clear) (struct mess_equation_st * eqn);
            /** Indicator variable if the \f$ AX \f$ objects needs to be cleared separately. */
            int to_clear;
        } AX;

        /**
         * @brief Structure defining \f$ Y = op(E) X \f$.
         *
         * Structure defining \f$ Y = op(E)X \f$ where \f$ X \f$ and \f$ Y \f$ are matrices.
         */
        struct EX_st{
            /** Post process function if \f$ Y = EX \f$ needs a global initialization.  */
            int (*generate)(struct mess_equation_st * eqn);
            /** Apply the mass matrix to a matrix. Compute \f$ Y=op(E)X \f$. */
            int (*apply)(struct mess_equation_st * eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
            /** Clear \f$ EX \f$ objects separately. */
            int (*clear) (struct mess_equation_st * eqn);
            /** Indicator variable if the \f$ EX \f$ objects needs to be cleared separately. */
            int to_clear;
        } EX;

        /**
         * @brief Structure defining \f$ op(A) Y  = X \f$.
         *
         * Structure defining \f$ op(A)Y = X \f$ where \f$ X \f$ and \f$ Y \f$ are matrices.
         */
        struct AINV_st{
            /** Post process function if \f$ op(A)Y = X \f$ needs a global initialization.  */
            int (*generate)(struct mess_equation_st * eqn);
            /** Apply the inverse operator to a matrix. Compute \f$ op(A)Y=X \f$. */
            int (*apply)(struct mess_equation_st * eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
            /** Clear \f$ AINV \f$ objects separately. */
            int (*clear) (struct mess_equation_st * eqn);
            /** Indicator variable if the \f$ AINV \f$ objects needs to be cleared separately. */
            int to_clear;
        } AINV;

        /**
         * @brief Structure defining \f$ op(E) Y  = X \f$.
         *
         * Structure defining \f$ op(E)Y=X \f$ where \f$ X \f$ and \f$ Y \f$ are matrices.
         */
        struct EINV_st{
            /** Post process function if \f$ op(E)Y = X \f$ needs a global initialization.  */
            int (*generate)(struct mess_equation_st * eqn);
            /** Apply the inverse mass matrix to a matrix. Compute \f$ op(E)Y=X \f$. */
            int (*apply)(struct mess_equation_st * eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
            /** Clear \f$ EINV \f$ objects separately. */
            int (*clear) (struct mess_equation_st * eqn);
            /** Indicator variable if the \f$ EINV \f$ objects needs to be cleared separately. */
            int to_clear;
        } EINV;

        /**
         * @brief Structure defining \f$ Y = op(A+ \alpha  p E)X \f$.
         *
         * Structure defining \f$ Y = op(A+ \alpha p E)X \f$ where \f$ Y \f$ and \f$ X \f$ are matrices.
         */
        struct ApEX_st{
            /** Post process function if \f$ Y = op(A+\alpha p E)X \f$ needs a global initialization.  */
            int (*generate) (struct mess_equation_st * eqn, mess_vector parameters);
            /** Apply \f$ op(A+\alpha*p*E)\f$  to a matrix. Compute \f$ y=op(A+\alpha p E)x \f$. */
            int (*apply) (struct mess_equation_st * eqn, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out);
            /** Clear \f$ ApEX \f$ objects separately. */
            int (*clear) (struct mess_equation_st * eqn);
            /** Indicator variable if the \f$ ApEX \f$ objects needs to be cleared separately. */
            int to_clear;
        }ApEX;

        /**
         * @brief Structure defining \f$ op(A+pE)Y = X \f$.
         *
         * Structure defining \f$ op(A+pE)Y=X \f$ where \f$ X \f$ and \f$ Y \f$ are vectors.
         */
        struct ApEINV_st{
            /** Post process function if \f$ op(A+pE)Y = X \f$ needs a global initialization.  */
            int (*generate) (struct mess_equation_st * eqn, mess_vector parameters);
            /** Apply the inverse of  \f$ op(A+pE) \f$ to a matrix. Compute \f$ op(A+pE)Y=X \f$. */
            int (*apply)    (struct mess_equation_st * eqn, mess_operation_t op, mess_double_cpx_t p, mess_int_t idx_p, mess_matrix in, mess_matrix out);
            /** Clear \f$ ApEINV \f$ objects separately. */
            int (*clear) (struct mess_equation_st * eqn);
            /** Indicator variable if the \f$ ApEINV \f$ objects needs to be cleared separately. */
            int to_clear;
        }ApEINV;

        /** Get shift parameter function. */
        int (*parameter)(struct mess_equation_st *eqn, mess_options opt, mess_status stat);

        /** Initialize the right hand side eqn_parent/eqn structure, if eqn_parent is set to @c NULL eqn is used. */
        int (*init_rhs) (struct mess_equation_st *eqn, mess_options opt);

        /** Final clean up function.  */
        int (*clear) ( struct mess_equation_st *eqn);

        /** Pointer to the internal data structure. */
        void *aux;
    };

    int mess_equation_init(mess_equation * eqn);
    int mess_equation_clean( mess_equation eqn );
    int mess_equation_clear(mess_equation * eqn);
    int mess_equation_set_rhs(mess_equation eqn, mess_options opt, mess_matrix rhs);
    int mess_equation_set_rhs2(mess_equation eqn, mess_matrix rhs);

    int mess_status_init(mess_status *status);
    int mess_status_clear(mess_status *status);
    int mess_status_print(mess_status stat);
    int mess_status_printfull(mess_status stat);
    int mess_status_copy(mess_status src, mess_status dest);
    mess_int_t mess_status_memsize(mess_status status);
    int mess_status_write(mess_status stat, const char *filename);

    int mess_options_init(mess_options *opt);
    int mess_options_clear(mess_options *opt);
    int mess_options_print(mess_options opt);
    int mess_options_set_proj_space(mess_options opt, mess_matrix proj_space);
    int mess_options_unset_proj_space(mess_options opt);
    int mess_options_set_K0(mess_options opt, mess_matrix K0);
    int mess_options_unset_K0(mess_options opt);
    int mess_options_set_W(mess_options opt, mess_matrix W);
    int mess_options_unset_W(mess_options opt);
    mess_int_t mess_options_memsize(mess_options opt);
    int mess_options_copy (mess_options src, mess_options dest);

    /** @} */


    /*-----------------------------------------------------------------------------
     *  ADI Algorithms
     *-----------------------------------------------------------------------------*/
    /** @addtogroup lrcfadi_adi
     * @{ */
    int mess_parameter(mess_equation eqn, mess_options opt, mess_status stat);
    int mess_lrcfadi_parameter(mess_equation eqn, mess_options opt, mess_status stat);
    int mess_lrcfadi_getritzvalues(mess_equation eqn, mess_options opt, mess_vector ritz );
    int mess_lrcfadi_residual(mess_equation eqn, mess_options opt, mess_matrix Z,  double *nrm);

    int mess_lradi(mess_equation eqn, mess_options opt, mess_status stat, mess_matrix Z);
    int mess_lrcfadi_adi(mess_equation eqn, mess_options opt, mess_status stat, mess_matrix Z);

    int mess_lrnm(mess_equation eqn, mess_options opt, mess_status stat, mess_matrix Z);
    int mess_lrcfadi_nm(mess_equation eqn, mess_options opt, mess_status stat, mess_matrix Z);

    int mess_lrcfadi_galerkin ( mess_equation eqn, mess_options opt, mess_equation_t eqn_type, mess_matrix Z );

    // line search functions
    int mess_lrcfadi_nm_linesearch_poly(mess_matrix W, mess_matrix Wnew, mess_matrix deltaK, mess_matrix deltaKnew, double *alpha, double *beta, double *delta, double *gamma, double *epsilon, double *zeta);
    int mess_lrcfadi_nm_linesearch_lambda(double *alpha, double* beta, double *delta, double *gamma, double *epsilon, double *zeta, double *lambda);
    int mess_lrcfadi_nm_linesearch(mess_matrix W, mess_matrix Wnew, mess_matrix deltaK, mess_matrix deltaKnew,double *lambda);

    /** @} */

    /*-----------------------------------------------------------------------------
     *  Misc functions
     *-----------------------------------------------------------------------------*/
    /** @addtogroup lrcfadi_misc
     * @{ */
    int mess_lrcfadi_lp_s(mess_double_cpx_t *p, mess_int_t lp, mess_double_cpx_t *set, mess_int_t lset, double *max_r, mess_int_t *ind);
    int mess_lrcfadi_lp_mnmx(mess_double_cpx_t *rw, mess_int_t lrw,  mess_int_t l0, mess_double_cpx_t *p, mess_int_t *lp);
    int mess_lrcfadi_fac2diff(mess_matrix Zold, mess_matrix Z, double *chg);
    int mess_lrcfadi_facFdiff(mess_matrix Zold, mess_matrix Z, double *chg);
    int mess_lrcfadi_para_wachspress(double a, double b, double alpha, double tol, double *p, mess_int_t *lp);
    int mess_lrcfadi_ccsvd(mess_matrix Z, double ccTol);
    int mess_lrcfadi_colcompress(mess_matrix Z, double ccTol);

    int mess_lrcfadi_res2(mess_matrix F, mess_matrix G, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2BK(mess_matrix F, mess_matrix G, mess_matrix Z, mess_matrix B, mess_matrix K, double *nrm);
    int mess_lrcfadi_res2g(mess_matrix F,mess_matrix M, mess_matrix G, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2gBK(mess_matrix F,mess_matrix M, mess_matrix G, mess_matrix Z, mess_matrix B, mess_matrix K, double *nrm);
    int mess_lrcfadi_res2t(mess_matrix F, mess_matrix G, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2gt(mess_matrix F,mess_matrix M, mess_matrix G, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2g_so(mess_matrix M, mess_matrix G, mess_matrix K, mess_matrix RHS, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2gt_so(mess_matrix M, mess_matrix G, mess_matrix K, mess_matrix RHS, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2nm(mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix Z, double *nrm);
    int mess_lrcfadi_res2nmg(mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix M, mess_matrix Z, double *nrm);

    /** @} */

    /*-----------------------------------------------------------------------------
     *  Additional Matrix Equation Solvers
     *-----------------------------------------------------------------------------*/
    /** @addtogroup lrcfadi_othersolvers
     * @{ */
    int mess_dense_nm_gmare(mess_matrix A, mess_matrix E, mess_matrix Q, mess_matrix G, mess_matrix X);
    int mess_dense_nm_gpare(mess_matrix A, mess_matrix E, mess_matrix Q, mess_matrix G, mess_matrix X);
    int mess_dense_nm_gmpare(mess_matrix X0, mess_matrix A, mess_matrix E, mess_matrix Q, mess_matrix G, mess_int_t plus, mess_int_t linesearch, mess_operation_t trans, mess_int_t maxit, mess_norm_t nrm,
                             double absres_tol, double relres_tol, double *absres, double *relres, mess_matrix RES, mess_matrix X);
    int mess_dense_res_gmpare(mess_matrix A, mess_matrix E,  mess_matrix Q, mess_matrix G, mess_matrix X, mess_int_t plus, mess_operation_t trans, mess_norm_t nrm_t, double *res, double *rel);

    /**
     * @brief Enumeration type for scaling parameters for sign function iteration.
     *
     * The @ref mess_sign_scale_t enumeration type defines scaling type for speedup sign function iteration.
     *
     */
    typedef enum {
        MESS_SIGN_SCALE_NONE  = 0,              /**< no scaling*/
        MESS_SIGN_SCALE_FRO   = 1,              /**< Frobenius  scaling \f$ \sqrt{\frac{||X_k||_F}{||X_k^{-1}||_F}} \f$ or \f$ \sqrt{\frac{||X_k||_F}{||EX_k^{-1}E||_F}} \f$   */
        MESS_SIGN_SCALE_DET   = 2,              /**< Determinante scaling \f$ |\det(X_k)|^{\frac{1}{n}} \f$ or \f$ \left(\frac{|\det(X_k)|}{|\det(E)|}\right)^{\frac{1}{n}} \f$   */
    } mess_sign_scale_t;

    int mess_lrcf_gsignfac(mess_matrix A, mess_matrix E, mess_matrix B, mess_int_t *maxit, double *tol, mess_matrix Z, mess_sign_scale_t scale);
    int mess_lrcf_signfac(mess_matrix A, mess_matrix B, mess_int_t *maxit,double *tol, mess_matrix Z, mess_sign_scale_t scale);
    int mess_lrcf_gsignber(mess_matrix A, mess_matrix E, mess_matrix B, mess_int_t *maxit, double *tol, mess_matrix Z, mess_sign_scale_t scale);
    int mess_lrcf_signber(mess_matrix A, mess_matrix B, mess_int_t *maxit, double *tol, mess_matrix Z, mess_sign_scale_t scale);

    /** @} */

    /*-----------------------------------------------------------------------------
     *  Equation Generators
     *-----------------------------------------------------------------------------*/
    /** @addtogroup lrcfadi_eqn
     * @{ */
    int mess_equation_lyap(mess_equation eqn, mess_options opt, mess_matrix A, mess_matrix E, mess_matrix B);
    int mess_equation_riccati(mess_equation eqn, mess_options opt, mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix C );

    int mess_equation_std_lyap(mess_equation eqn, mess_options opt, mess_matrix A, mess_matrix B);
    int mess_equation_std_riccati(mess_equation eqn, mess_options opt, mess_matrix A, mess_matrix B, mess_matrix C );

    int mess_equation_glyap(mess_equation eqn, mess_options opt, mess_matrix A, mess_matrix E, mess_matrix B);
    int mess_equation_griccati(mess_equation eqn, mess_options opt, mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix C );

    int mess_equation_glyap_so1(mess_equation e, mess_options opt , mess_matrix M, mess_matrix D, mess_matrix K, mess_matrix B,double lowerbound, double upperbound);
    int mess_equation_griccati_so1(mess_equation e, mess_options opt , mess_matrix M, mess_matrix D, mess_matrix K, mess_matrix B, mess_matrix C, double lowerbound, double upperbound);

    int mess_equation_glyap_so2(mess_equation e, mess_options opt , mess_matrix M, mess_matrix D, mess_matrix K, mess_matrix B,double lowerbound, double upperbound);
    int mess_equation_griccati_so2(mess_equation e, mess_options opt , mess_matrix M, mess_matrix D, mess_matrix K, mess_matrix B, mess_matrix C, double lowerbound, double upperbound);

    int mess_equation_glyap_dae2(mess_equation e, mess_options opt, mess_matrix M, mess_matrix A, mess_matrix G, mess_matrix B, double delta);
    int mess_equation_griccati_dae2(mess_equation e, mess_options opt, mess_matrix M, mess_matrix A, mess_matrix G, mess_matrix B, mess_matrix C, double delta);

    int mess_equation_glyap_dae1(mess_equation e, mess_options opt , mess_matrix E11, mess_matrix A11, mess_matrix A12, mess_matrix A21, mess_matrix A22,  mess_matrix rhs);
    int mess_equation_griccati_dae1(mess_equation e, mess_options opt , mess_matrix E11, mess_matrix A11, mess_matrix A12, mess_matrix A21, mess_matrix A22,  mess_matrix B, mess_matrix C);

    int mess_equation_stable(mess_equation e, mess_options opt , mess_equation base , mess_matrix B, mess_matrix K );

    /** @} */

    /*-----------------------------------------------------------------------------
     *  Equation Apply Functions
     *-----------------------------------------------------------------------------*/
    /** @addtogroup lrcfadi_eqn_apply
     * @{ */

    mess_int_t mess_equation_dim(mess_equation eqn);
    mess_int_t mess_equation_has_A(mess_equation eqn);
    mess_int_t mess_equation_has_As(mess_equation eqn);
    mess_int_t mess_equation_has_E(mess_equation eqn);
    mess_int_t mess_equation_has_Es(mess_equation eqn);
    mess_int_t mess_equation_has_ApE(mess_equation eqn);
    mess_int_t mess_equation_has_ApEs(mess_equation eqn);

    int mess_equation_A_pre(mess_equation eqn);
    int mess_equation_A_post(mess_equation eqn);
    int mess_equation_A_apply(mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
    int mess_equation_A_apply_vector(mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out);

    int mess_equation_E_pre(mess_equation eqn);
    int mess_equation_E_post(mess_equation eqn);
    int mess_equation_E_apply(mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
    int mess_equation_E_apply_vector(mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out);

    int mess_equation_As_pre(mess_equation eqn);
    int mess_equation_As_post(mess_equation eqn);
    int mess_equation_As_apply(mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
    int mess_equation_As_apply_vector(mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out);

    int mess_equation_Es_pre(mess_equation eqn);
    int mess_equation_Es_post(mess_equation eqn);
    int mess_equation_Es_apply(mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out);
    int mess_equation_Es_apply_vector(mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out);

    int mess_equation_ApE_pre(mess_equation eqn, mess_vector parameters);
    int mess_equation_ApE_post(mess_equation eqn);
    int mess_equation_ApE_apply(mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_matrix in, mess_matrix out);
    int mess_equation_ApE_apply_vector(mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_vector in, mess_vector out);

    int mess_equation_ApEs_pre(mess_equation eqn, mess_vector parameters);
    int mess_equation_ApEs_post(mess_equation eqn);
    int mess_equation_ApEs_apply(mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_matrix in, mess_matrix out);
    int mess_equation_ApEs_apply_vector(mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_vector in, mess_vector out);
    /** @} */


#ifdef __cplusplus
}
#endif

#endif /* LRCF_ADI_H_ */
/** \}@ */
