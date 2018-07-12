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
#include "MESS_equation.h"


/*-----------------------------------------------------------------------------
 *  init the installed variable and define install function
 *-----------------------------------------------------------------------------*/
bool MESS_equation::type_installed = false;


void MESS_equation::install(void){
    if(!MESS_equation::type_installed){
        MESS_equation::register_type();
        octave_stdout << "Install " << MESS_equation::static_type_name() << " at type-id = " << MESS_equation::static_type_id() << "\n";

    }
    MESS_equation::type_installed = true;
}
