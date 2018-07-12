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
 * @addtogroup test_matrix
 * @{
 * @file tests/matrix/check_eps.c
 * @brief Check the machine epsilon calculation.
 * @author @koehlerm
 * @test
 * This function checks the @ref mess_eps function definde in eps.c that means it checks if the machine epsilon is computed
 * correctly.
 *
 * @}
 *
 */

#include <stdio.h>
#include "mess/mess.h"

int main (int argc, char **argv) {
    mess_init();
    double eps = mess_eps();
    if ( eps < 3e-16) return 0;
    else return 1;
}

