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
 * @file include/mess/dynsys.h
 * @brief Basic functions for dynamical systems.
 * @author @koehlerm
 */

#ifndef DYNSYS_H_
#define DYNSYS_H_

#ifdef  __cplusplus
extern "C" {
#endif
    /** @addtogroup dynsys_interface
     * @{ */
    /** Define a LTI system. */
#define MESS_DYNSYS_LTI 0x01
    /** Define a generalized LTI system. */
#define MESS_DYNSYS_GLTI    0x02

    /** Define a second order system. */
#define MESS_DYNSYS_2ND 0x03

    /** Check if a @ref mess_dynsys  object is a LTI system. */
#define MESS_IS_DYNSYS_LTI(x) ((*x).type == MESS_DYNSYS_LTI)
    /** Check if a @ref mess_dynsys  object is a GLTI system. */
#define MESS_IS_DYNSYS_GLTI(x) ((*x).type == MESS_DYNSYS_GLTI)
    /** Check if a @ref mess_dynsys  object is a second order system. */
#define MESS_IS_DYNSYS_2ND(x) ((*x).type == MESS_DYNSYS_2ND)

    /** Check if a @ref mess_dynsys  object is a SISO system.  */
#define MESS_IS_DYNSYS_SISO(x) (((*x).inputs==1) && ((*x).outputs==1))

    /**
     * @brief Representation of a dynamical system.
     *
     * The @ref mess_dynsys_st  structure contains all information to represent a dynamical system.
     */
    struct mess_dynsys_st {
        mess_matrix A;              /**< System matrix \f$ A \f$ of a first order system */
        mess_matrix B;              /**< Input matrix \f$ B \f$ of a first and second order system */
        mess_matrix C;              /**< Output matrix \f$ C \f$ of a first order system */
        mess_matrix E;              /**< Mass matrix \f$ E  \f$ of a first order system */
        mess_matrix M;              /**< Matrix for the second derivative in a second order system */
        mess_matrix G;              /**< Matrix for the first derivative in a second order system */
        mess_matrix K;              /**< Matrix for the state in a second order system */
        mess_matrix Cp;             /**< Postion output matrix \f$ C_p \f$ of a second order system */
        mess_matrix Cv;             /**< Velocity output matrix \f$ C_v \f$ of a second order system */
        mess_int_t dim;             /**< State space dimension of the system */
        mess_int_t inputs;          /**< Number of inputs */
        mess_int_t outputs;         /**< Number of outputs */
        unsigned short type;        /**< Type of the dynamical system */
    };

    /**
     * @brief  Type definition of easy usage of mess_dynsys_st.
     *
     * The mess_dynsys type definition of easy usage of @ref mess_dynsys_st .*/
    typedef struct mess_dynsys_st *  mess_dynsys;

    char *mess_dynsysstr(unsigned short type);
    int mess_dynsys_init(mess_dynsys * sys);
    int mess_dynsys_clear(mess_dynsys * sys);
    int mess_dynsys_lti(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C);
    int mess_dynsys_lti_copy(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C);
    int mess_dynsys_glti(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C,mess_matrix E);
    int mess_dynsys_glti_copy(mess_dynsys lti, mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E);
    int mess_dynsys_2nd(mess_dynsys sys, mess_matrix M, mess_matrix G, mess_matrix K, mess_matrix B, mess_matrix Cp, mess_matrix Cv);
    int mess_dynsys_2nd_copy(mess_dynsys sys, mess_matrix M, mess_matrix G, mess_matrix K,mess_matrix B, mess_matrix Cp, mess_matrix Cv);
    int mess_dynsys_2nd_to_1st( mess_dynsys sys2nd, mess_dynsys sys1st );

    int mess_dynsys_printinfo(mess_dynsys sys);
    int mess_dynsys_isstable(mess_dynsys sys, int *stable);

    int mess_dynsys_evaltransfer(mess_dynsys sys, double aa, double bb, mess_int_t nsample, mess_vector omega, mess_vector G , mess_matrix *Gout);
    int mess_dynsys_transfer_diff(mess_int_t nsample, mess_matrix *Go, mess_matrix *Gr, double *err2, double *rel2, double *errinf, double *relinf, mess_vector diffvec );

    int mess_dynsys_project_to_lti(mess_dynsys sys, mess_matrix V, mess_matrix W, mess_dynsys red );
    int mess_dynsys_project_to_glti(mess_dynsys sys, mess_matrix V, mess_matrix W, mess_dynsys red );
    int mess_dynsys_project_2nd_to_1st(mess_dynsys sys2nd, mess_matrix V, mess_matrix W, mess_dynsys sys1st);
    int mess_dynsys_project_2nd_to_2nd(mess_dynsys sys, mess_matrix V, mess_matrix W, mess_dynsys red );
    int mess_matrix_isstableg(mess_matrix A, mess_matrix E, int *stable);
    int mess_matrix_isstable(mess_matrix A, int *stable);

    /** @} */

#ifdef __cplusplus
}
#endif

#endif /* LRCF_ADI_H_ */
/** \}@ */
