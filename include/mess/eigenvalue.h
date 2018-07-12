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
 * @file include/mess/eigenvalue.h
 * @brief Interface for Eigenvalue computation.
 * @author @koehlerm
*/

#ifndef EIGENVALUE_H_
#define EIGENVALUE_H_

#include "mess/direct.h"
#include "mess/lrcf_adi.h"
#ifdef __cplusplus
extern "C" {
#endif

    /** @addtogroup eigenvalues_sparse
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Sparse matrices
     *-----------------------------------------------------------------------------*/
    // The Arnoldi process
    int mess_eigen_arnoldi(mess_matrix A, mess_matrix B, mess_matrix K, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V);
    int mess_eigen_arnoldig(mess_matrix A, mess_matrix U, mess_matrix W, mess_direct E, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V);
    int mess_eigen_arnoldig_inv(mess_direct A, mess_matrix U, mess_matrix W, mess_matrix E, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V);
    int mess_eigen_arnoldi_inv(mess_direct A, mess_matrix B, mess_matrix K,  mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V);

    int mess_eigen_arnoldi_template(mess_mvpcall A, mess_int_t k, mess_vector sv, mess_matrix H, mess_matrix V);
    int mess_eigen_arnoldi_template_nrm(mess_mvpcall A, mess_int_t k, mess_vector sv, double *abseig);
    int mess_eigen_lanczos_template_nrm(mess_mvpcall A, mess_int_t k, mess_vector sv, mess_vector diag, mess_vector subdiag, mess_matrix V, double *abseig);

    int mess_eigen_eigs(mess_matrix A, mess_int_t kp, mess_int_t km, mess_vector r, mess_vector ev);
    int mess_eigen_eiggs(mess_matrix A, mess_matrix E, mess_int_t kp, mess_int_t km, mess_vector r, mess_vector ev);


#ifdef MESS_HAVE_ARPACK
    // Function only available with @arpack support
    /**
     * @brief Enumeration to select the desired eigenvalues from @arpack.
     *
     * The @ref mess_eigen_arpack_which_t  enumeration selects the desired eigenvalues in the Arnoldi process
     * performed by \ref mess_eigen_arpack.
     * \see mess_eigen_arpack_options_t
     * \see mess_eigen_arpack
     * \see mess_eigen_arpack_template
     * \remarks The enumeration is only available if the preprocessor macro \c MESS_HAVE_APRACK is defined and
     * @mess is compiled with @arpack support.
     * */
    typedef enum { /** @arpack should compute eigenvalues of largest magnitude. */
               MESS_EIGEN_ARPACK_LARGE_MAGNITUDE = 0 ,
               /** @arpack should compute eigenvalues of smallest magnitude.  */
               MESS_EIGEN_ARPACK_SMALL_MAGNITUDE,
               /** @arpack should compute eigenvalues of largest real part.
            *  \attention This option is only valid for \ref mess_eigen_arpack and
            * \ref mess_eigen_arpack_template. */
               MESS_EIGEN_ARPACK_LARGE_REALPART,
               /** @arpack should compute eigenvalues of smallest real part.
            * \attention This option is only valid for \ref mess_eigen_arpack and
            * \ref mess_eigen_arpack_template. */
               MESS_EIGEN_ARPACK_SMALL_REALPART,
               /** @arpack should compute eigenvalues of largest imaginary part.
            * \attention This option is only valid for \ref mess_eigen_arpack and
            * \ref mess_eigen_arpack_template. */
               MESS_EIGEN_ARPACK_LARGE_IMAGPART,
               /** @arpack should compute eigenvalues of smallest imaginary part.
            * \attention This option is only valid for \ref mess_eigen_arpack and
            * \ref mess_eigen_arpack_template.  */
               MESS_EIGEN_ARPACK_SMALL_IMAGPART,
               /** @arpack should compute the largest algebraic eigenvalues.
            * \attention This option is only valid for \ref mess_eigen_arpack_lanczos and
            * \ref mess_eigen_arpack_lanczos_template.*/
               MESS_EIGEN_ARPACK_LARGE_ALGEBRAIC,
               /** @arpack should compute the smallest algebraic eigenvalues.
            * \attention This option is only valid for \ref mess_eigen_arpack_lanczos and
            * \ref mess_eigen_arpack_lanczos_template.*/
               MESS_EIGEN_ARPACK_SMALL_ALGEBRAIC,
               /** @arpack should compute the eigenvalues at both ends of the spectrum.
            * \attention This option is only valid for \ref mess_eigen_arpack_lanczos and
            * \ref mess_eigen_arpack_lanczos_template.*/
               MESS_EIGEN_ARPACK_BOTH_ENDS
        } mess_eigen_arpack_which_t;
    /**
     * @brief Options structure for @arpack .
     *
     * The @ref mess_eigen_arpack_options_t  structure defines the behavior of @arpack. \n
     * Using the \c MESS_EIGEN_ARPACK_DEFAULT marco the structure can be set to @matlab defaults during
     * initialization.  \n
     * The structure should be initialized using
     * \code{.C}
     * @ref mess_eigen_arpack_options_t arpack_opt = MESS_EIGEN_ARPACK_DEFAULT;
     * \endcode
     *
     * \see mess_eigen_arpack
     * \see mess_eigen_arpack_which_t
     *
     * \remarks The structure is only available if the preprocessor macro \c MESS_HAVE_APRACK is defined and
     * @mess is compiled with @arpack support.
     * */
    typedef struct {
        /** Stopping criterion: The relative accuracy of the computed Ritz-values.
         * The default is the machine precision @ref mess_eps ().  */
        double tol;
        /** Starting vector for the Arnoldi process.
         * If the vector is \c NULL a randomly chosen vector is used.  */
        mess_vector b0;
        /** Number of Arnoldi vectors to be generated.
         * The default is \c MIN(MAX(NEV+10,50),N-1), where \c NEV is the number of desired
         * eigenvalues and \c N the size of the eigenvalues problem. */
        mess_int_t ncv;
        /** Maximum number of Arnoldi update iterations allowed. */
        mess_int_t maxit;
        /** Kind of the desired eigenvalues. See \ref mess_eigen_arpack_which_t for details. */
        mess_eigen_arpack_which_t which;
    } mess_eigen_arpack_options_t;

    /**
     * @brief Default initializer for the \ref mess_eigen_arpack_options_t structure.
     *
     * The \c MESS_EIGEN_ARPACK_DEFAULT macro initializes the mess_eigen_arpack_options_t structure to
     * a useful set of default values. \n
     * The tolerance is set to the machine precision, the start vector
     * to \c NULL, the maximum number of iterations to 300 and the number of Arnoldi vectors to 50. */
    #define MESS_EIGEN_ARPACK_DEFAULT { .tol = 0.0, .b0 = NULL, .ncv = 100, .maxit = 300, .which = MESS_EIGEN_ARPACK_LARGE_MAGNITUDE }

    int mess_eigen_arpack(mess_matrix A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V);
    int mess_eigen_arpack_template(mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V);
    int mess_eigen_arpack_lanczos(mess_matrix A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V);
    int mess_eigen_arpack_lanczos_template(mess_mvpcall A, mess_int_t nev, mess_eigen_arpack_options_t opt, mess_vector ev, mess_matrix V);

#endif

    /** @}  */


    /** @addtogroup eigenvalues_dense
     * @{ */
    /*-----------------------------------------------------------------------------
     *  Dense Matrices
     *-----------------------------------------------------------------------------*/
    int mess_eigen_eig(mess_matrix A, mess_vector ew, mess_matrix ev);
    int mess_eigen_eigg(mess_matrix A, mess_matrix B, mess_vector ew, mess_vector alpha, mess_vector beta, mess_matrix ev);

    int mess_eigen_schur_complex( mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV );
    int mess_eigen_schur_complex_zgees(mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV );
    int mess_eigen_schur_complex_post(mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV );

    int mess_eigen_schur(mess_matrix A, mess_matrix T, mess_matrix U, mess_vector EV );
    int mess_eigen_gschur(mess_matrix A, mess_matrix B, mess_matrix S, mess_matrix T, mess_matrix U, mess_matrix V, mess_vector EVA, mess_vector EVB, mess_vector EV);
    int mess_eigen_gschur_complex(mess_matrix A, mess_matrix B, mess_matrix S, mess_matrix T, mess_matrix U, mess_matrix V, mess_vector EVA, mess_vector EVB, mess_vector EV);

    int mess_eigen_svd(mess_matrix A, mess_vector S, mess_matrix U, mess_matrix V);
    int mess_eigen_svd_econ(mess_matrix A, mess_vector S, mess_matrix U, mess_matrix V);

    int mess_eigen_schur_to_evd(mess_matrix T, mess_matrix Q, mess_matrix V);
    int mess_eigen_hessenberg_schur(mess_matrix hess, mess_vector ew, mess_matrix Q, mess_matrix T);
    int mess_eigen_hessenberg_abs_largest(mess_matrix hess, mess_vector evector, mess_double_cpx_t *abs_largest_ev);

    int mess_eigen_sign(mess_matrix A, mess_matrix Z, mess_sign_scale_t scale);
    int mess_eigen_gsign(mess_matrix A, mess_matrix B, mess_matrix Z, mess_sign_scale_t scale);

    /** @}  */




#ifdef __cplusplus
}
#endif

#endif /*EIGENVALUE_H_*/
/** @} */
