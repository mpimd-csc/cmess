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
 * @file include/mess/easyfrontend.h
 * @brief Easy frontends for several functions.
 * @author @koehlerm
 */
#ifndef MESS_EASYFRONTEND_H
#define MESS_EASYFRONTEND_H

#include "mess/config.h"
#include "mess/mess.h"


#ifdef __cplusplus
extern "C" {
#endif
    /** @addtogroup easyfrontend
     * @{ */
    int mess_lyap(mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix Z);
    int mess_care(mess_matrix A, mess_matrix E, mess_matrix B, mess_matrix C, mess_matrix Z);

    int mess_glyap(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix Y, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix X);
    int mess_gstein(mess_operation_t op, mess_matrix A, mess_matrix E, mess_matrix Y, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix X);

    int mess_tglyap(mess_operation_t op, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix Y, mess_matrix X);
    int mess_tgstein(mess_operation_t op, mess_matrix Ahat, mess_matrix QA, mess_matrix Ehat, mess_matrix QE, mess_matrix Y,  mess_matrix X);

    /** @}  */


#ifdef __cplusplus
}
#endif


#endif
