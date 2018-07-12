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
 * @file include/mess/solver.h
 * @brief Interface to all iterative linear system solvers.
 * @author @koehlerm
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "mess/matrix.h"
#include "mess/vector.h"
#include "mess/direct.h"

#ifdef __cplusplus
extern "C" {
#endif
    /**
     * @addtogroup itsolver_ds
     * @{ */

    /**
     * @brief Representation of status information of an iterative solvers.
     *
     * The @ref __mess_solver_status structure represents the status information about an iterative solver.
     * */
    struct __mess_solver_status {
        double res;                 /**< Residual norm */
        double relres;              /**< Relative residual */
        int need_restart;           /**< True if a restart is needed */
        int converged;              /**< True if the algorithm converged */
        mess_int_t it;              /**< Number of calculated iterations */
        mess_int_t num_mvp;         /**< Number of computed matrix-vector products */
        mess_int_t restarts;
    };

    /**
     * @brief Type definition for status information of an iterative solvers.
     *
     * The @ref mess_solver_status type definition defines the status information about an iterative solver.
     * */
    typedef struct __mess_solver_status *  mess_solver_status;

    /**
     * @brief Representation of options for iterative solvers.
     *
     * The @ref __mess_solver_options structure represents options for an iterative solver.
     **/
    struct __mess_solver_options {
        mess_int_t maxit;                           /**< Maximum number of iterations */
        mess_int_t restarts;                            /**< Maximum number of restarts, if the value is -1 is runs until the algorithm converged. */
        double tol;                             /**< Relative tolerance */
        // double shift;                /**< solve shifted systems */
        double omega;                               /**< Parameter for SOR and SSOR */
        void (*stepdebug)(double res, double relres, mess_int_t it, void *opt);     /**< Step debug function, evaluated every step */
        void *aux_stepdebug;                            /**< Auxilary data for stepdebug */
    };

    /**
     * @brief Type definition for options of an iterative solvers.
     *
     * The @ref mess_solver_options type definition defines options for an iterative solver.
     * */
    typedef struct __mess_solver_options *  mess_solver_options;


    /**
     * @brief Representation of a preconditioner.
     *
     * The @ref mess_precond_st structure represents a preconditioner.
     **/
    struct mess_precond_st {
        void *data;         /**< Auxiliary data */
        int (*clear)(struct mess_precond_st *myself);
        int (*solve)(struct mess_precond_st *myself, mess_solver_options opt, mess_vector in, mess_vector out); /**< Solve \f$ y = M^{-1} x\f$ */
        int type;           /**< Type or name of the precondtioner */
    };

    /**
     * @brief Type definition for a preconditioner.
     *
     * The @ref mess_precond type definition defines a CSS preconditioner.
     **/
    typedef struct mess_precond_st*  mess_precond;

    /**
     * @brief Definition of a diagonal preconditioner.
     *
     * @ref MESS_PRECOND_DIAG defines a diagonal preconditioner.
     **/
#define MESS_PRECOND_DIAG   1

    /**
     * @brief Definition of a \f$ ILU(0) \f$ preconditioner.
     *
     * @ref MESS_PRECOND_ILU0 defines a \f$ ILU(0) \f$ preconditioner.
     **/
#define MESS_PRECOND_ILU0   2

    /**
     * @brief Definition of a \f$ ILU(k) \f$ preconditioner.
     *
     * @ref MESS_PRECOND_ILUK defines a \f$ ILU(k) \f$ preconditioner.
     **/
#define MESS_PRECOND_ILUK   3

    int mess_solver_options_init(mess_solver_options *opt);
    int mess_solver_options_print ( mess_solver_options opt);
    int mess_solver_options_clear(mess_solver_options *opt);
    int mess_solver_status_init(mess_solver_status *opt);
    int mess_solver_status_print ( mess_solver_status opt);
    int mess_solver_status_clear(mess_solver_status *opt);
    int mess_solver_status_clean ( mess_solver_status stat);
    /** @} */

    /** @addtogroup itsolver_sol
     * @{ */
    int mess_solver_gmres(  mess_mvpcall matrix,
            mess_precond pre,
            mess_vector b,
            mess_vector x,
            mess_solver_options opt,
            mess_solver_status  stat);
    int mess_solver_gmres_restart(  mess_mvpcall matrix,
            mess_precond pre,
            mess_vector b,
            mess_vector x,
            mess_solver_options opt,
            mess_solver_status  stat);
    int mess_solver_bicgstab(   mess_mvpcall matrix,
            mess_precond pre,
            mess_vector b,
            mess_vector x,
            mess_solver_options opt,
            mess_solver_status  stat);
    int mess_solver_cg  (   mess_mvpcall matrix,
            mess_precond pre,
            mess_vector b,
            mess_vector x,
            mess_solver_options opt,
            mess_solver_status stat
            );
    int mess_solver_sor(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt, mess_solver_status stat);
    int mess_solver_ssor(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt,mess_solver_status stat);
    int mess_solver_gaussseidel(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt,mess_solver_status stat);
    int mess_solver_jacobi(mess_matrix matrix, mess_precond pre, mess_vector b, mess_vector x, mess_solver_options opt, mess_solver_status stat);
    int mess_solver_jacobi_convergence(mess_matrix matrix, int *convergence);

    /** @}  */

    /** @addtogroup itsolver_pre
     * @{ */

    int mess_precond_init(mess_precond *precond);
    int mess_precond_clear(mess_precond *precond);

    int mess_precond_diag(mess_precond pre,  mess_matrix mat);
    int mess_precond_ilu0(mess_precond pre, mess_matrix matrix );
    int mess_precond_iluk(mess_precond pre, mess_matrix matrix , mess_int_t leveloffill);
    int mess_precond_sor ( mess_precond pre, mess_matrix matrix, mess_int_t it, double omega );
    int mess_precond_ssor ( mess_precond pre, mess_matrix matrix, mess_int_t it, double omega );

    /** @}  */


#ifdef __cplusplus
}
#endif
#endif /* SOLVER_H_ */
/** \}@ */
