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
#include "MESS_multidirect.h"
#include "MESS_error.h"


DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (MESS_multidirect, "mess_multidirect", "mess_multidirect");


/*-----------------------------------------------------------------------------
 *  constructors
 *-----------------------------------------------------------------------------*/
MESS_multidirect::MESS_multidirect(void):ptr(NULL) { }
MESS_multidirect::MESS_multidirect(mess_multidirect d):ptr(d) { }
MESS_multidirect::MESS_multidirect(const MESS_multidirect & d):ptr(d.ptr){ }


/*-----------------------------------------------------------------------------
 *  getter
 *-----------------------------------------------------------------------------*/
mess_multidirect MESS_multidirect::get_ptr(void) const { return this->ptr; }


/*-----------------------------------------------------------------------------
 *  print functions octave
 *-----------------------------------------------------------------------------*/
void MESS_multidirect::print(std::ostream & os, bool pr_as_read_syntax) const {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_multidirect::print(std::ostream & os, bool pr_as_read_syntax) {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_multidirect::print_raw(std::ostream & os, bool pr_as_read_syntax) const {
    os << MESS_multidirect::static_type_name();
#ifdef MESS_DEBUG
    os << " ( " << this->ptr << " ) ";
#endif
}

bool MESS_multidirect::print_name_tag (std::ostream& os, const std::string& name) const {
    os << name << " = ";
    return true;
}


/*-----------------------------------------------------------------------------
 *  overrides from octave_base_value
 *-----------------------------------------------------------------------------*/
bool MESS_multidirect::is_constant(void) const { return true; }
bool MESS_multidirect::is_defined(void) const { return true; }
bool MESS_multidirect::print_as_scalar(void) const { return true; }
bool MESS_multidirect::is_empty(void) const { return this->ptr?false:true; }


/*-----------------------------------------------------------------------------
 *  conversion operator
 *-----------------------------------------------------------------------------*/
MESS_multidirect::operator mess_multidirect(void) const{ return this->ptr; }
MESS_multidirect::operator bool(void) const{ return this->ptr != NULL; }


/*-----------------------------------------------------------------------------
 *  destructor
 *-----------------------------------------------------------------------------*/
MESS_multidirect::~MESS_multidirect(void){
    if(this->ptr){
        mess_multidirect_clear(&(this->ptr));
    }
    this->ptr=NULL;
}
