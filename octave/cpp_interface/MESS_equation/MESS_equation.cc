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
#include <octave/ov-complex.h>
#include <octave/ov-re-sparse.h>
#include <math.h>
#include "mess/mess.h"
#include "MESS_equation.h"
#include "MESS_error.h"


DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (MESS_equation, "mess_equation", "mess_equation");


/*-----------------------------------------------------------------------------
 *  constructors
 *-----------------------------------------------------------------------------*/
MESS_equation::MESS_equation(void):ptr(NULL){ }
MESS_equation::MESS_equation(mess_equation e):ptr(e){ }
MESS_equation::MESS_equation(const MESS_equation & eqn):ptr(eqn.ptr){ }

/*-----------------------------------------------------------------------------
 *  print functions octave
 *-----------------------------------------------------------------------------*/
void MESS_equation::print(std::ostream & os, bool pr_as_read_syntax) const {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_equation::print(std::ostream & os, bool pr_as_read_syntax) {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_equation::print_raw(std::ostream & os, bool pr_as_read_syntax) const {
    os << MESS_equation::static_type_name();
#ifdef MESS_DEBUG
    os << " ( " << this->ptr << " ) ";
#endif
}


bool MESS_equation::print_name_tag (std::ostream& os, const std::string& name) const {
    os << name << " = ";
    return true;
}


/*-----------------------------------------------------------------------------
 *  getter
 *-----------------------------------------------------------------------------*/
mess_equation MESS_equation::get_ptr(void) const{ return this->ptr; }


/*-----------------------------------------------------------------------------
 *  overrides from octave_base_value
 *-----------------------------------------------------------------------------*/
bool MESS_equation::is_constant(void) const { return true; }
bool MESS_equation::is_defined(void) const { return true; }
bool MESS_equation::print_as_scalar(void) const { return true; }
bool MESS_equation::is_empty(void) const { return this->ptr?false:true; }


/*-----------------------------------------------------------------------------
 *  conversion operator
 *-----------------------------------------------------------------------------*/
MESS_equation::operator mess_equation(void) const { return this->ptr; }
MESS_equation::operator bool(void) const{ return this->ptr != NULL; }


/*-----------------------------------------------------------------------------
 *  destructor
 *-----------------------------------------------------------------------------*/
MESS_equation::~MESS_equation(void){
    if(this->ptr){
        mess_equation_clear(&(this->ptr));
    }
    this->ptr=NULL;
}
