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
 * @file include/mess/freelist.h
 * @brief A structure and manage allocation and free of @mess structures.
 * @author @mbehr
 */

#ifndef FREELIST_H_
#define FREELIST_H_


#ifdef  __cplusplus
extern "C" {
#endif
    /** @addtogroup misc
     * @{ */


    /**
     * @internal
     * @brief Structure containing a minimal memory manage in order to avoid memory leaks in the @python and @mexmess interface.
     *
     * Due to the fact that @c C does not have a gargabe collector and @mess does not implement reference counting
     * we have to track the allocated data in some cases.
     *
     * @attention Interal use only.
     * */
    typedef struct mess_freelist_st {
        void **ptrs;                    /**< Array containing arbitrary pointers.*/
        mess_vector *vectors;           /**< Array containing the @ref mess_vector instances.*/
        mess_matrix *matrices;          /**< Array containing the @ref mess_matrix instances.*/
        mess_equation *eqns;            /**< Array containing the @ref mess_equation instances. */
        mess_options *opts;             /**< Array containing the @ref mess_options instances. */
        mess_status  *stats;            /**< Array containing the @ref mess_status instances. */
        size_t n_ptrs;                  /**< Length of the ptrs array.   */
        size_t n_vectors;               /**< Length of the vector array. */
        size_t n_matrices;              /**< Length of the matrix array. */
        size_t n_eqns;                  /**< Length of the eqns array. */
        size_t n_opts;                  /**< Length of the opts array. */
        size_t n_stats;                 /**< Length of the stats array. */
    } *mess_freelist;

    int mess_freelist_init(mess_freelist *list);
    int mess_freelist_print(mess_freelist list);
    int mess_freelist_add_ptr(mess_freelist list, void *ptr);
    int mess_freelist_add_mess_vector(mess_freelist list, mess_vector vec);
    int mess_freelist_add_mess_matrix(mess_freelist list, mess_matrix mat);
    int mess_freelist_add_mess_equation(mess_freelist list, mess_equation eqn);
    int mess_freelist_add_mess_options(mess_freelist list, mess_options opt);
    int mess_freelist_add_mess_status(mess_freelist list, mess_status stat);
    int mess_freelist_clear(mess_freelist *list);


    /** @} */

#ifdef __cplusplus
}
#endif

#endif
