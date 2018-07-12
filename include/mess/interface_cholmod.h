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
 * @file include/mess/interface_cholmod.h
 * @brief Interface of conversion from a @ref mess_matrix or a @ref mess_vector to a @cholmod datastructure (and vice versa).
 * @author @koehlerm
 */

#ifndef MESS_INTERFACE_CHOLMOD_H
#define MESS_INTERFACE_CHOLMOD_H

#ifdef MESS_USE_SUITESPARSE3
#ifdef I
#undef I
#endif
#include <cholmod.h>
#define I _Complex_I
#else
#include <cholmod.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
    /** @addtogroup interface_cholmod
     * @{ */
    /*-----------------------------------------------------------------------------
     *  @cholmod interface functions
     *-----------------------------------------------------------------------------*/
    int mess_matrix_convert_dense_to_cholmod(mess_matrix A, cholmod_dense **A_cholmod, cholmod_common *c);
    int mess_matrix_convert_cholmod_to_dense(cholmod_dense *A_cholmod, mess_matrix A, cholmod_common *c);
    int mess_matrix_convert_csc_to_cholmod(mess_matrix A, cholmod_sparse **A_cholmod, cholmod_common *c);
    int mess_matrix_convert_cholmod_to_csc(cholmod_sparse *A_cholmod, mess_matrix A, cholmod_common *c);
    int mess_vector_convert_dense_to_cholmod(mess_vector v, cholmod_dense **v_cholmod, cholmod_common *c);
    int mess_vector_convert_cholmod_to_dense(cholmod_dense *v_cholmod, mess_vector v, cholmod_common *c);

    /** @}  */

#ifdef __cplusplus
}
#endif



#endif /* end of include guard: INTERFACE_CHOLMOD_H */
