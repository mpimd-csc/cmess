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

#include <octave/oct.h>
#include <octave/ops.h>
#include "mess/mess.h"
#include "MESS_matrix.h"
#include "MESS_error.h"


/*-----------------------------------------------------------------------------
 *  implement binary operations MESS_matrix <-> octave_matrix
 *
 *  OCTMESS_GENERATOR macro will declare the necessary binary operations.
 *-----------------------------------------------------------------------------*/
#include "MESS_matrix_ops_mm_generator.h"

OCTMESS_GENERATOR(octave_matrix, matrix_value);


