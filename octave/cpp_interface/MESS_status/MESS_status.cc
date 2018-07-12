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
#include "MESS_status.h"
#include "MESS_error.h"


DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (MESS_status, "mess_status", "mess_status");


/*-----------------------------------------------------------------------------
 *  constructors
 *-----------------------------------------------------------------------------*/
MESS_status::MESS_status(void):ptr(NULL) { }
MESS_status::MESS_status(mess_status s):ptr(s) { }
MESS_status::MESS_status(const MESS_status & s):ptr(s.ptr){ }


/*-----------------------------------------------------------------------------
 *  getter
 *-----------------------------------------------------------------------------*/
mess_status MESS_status::get_ptr(void) const { return this->ptr; }


/*-----------------------------------------------------------------------------
 *  deep copy
 *-----------------------------------------------------------------------------*/
MESS_status MESS_status::deepcopy(void) const{
    MESS_status stat;
    if(*this){
        mess_status ptr = stat.get_ptr();
        mess_status_init(&ptr);
        mess_status_copy(this->ptr,ptr);
    }
    return stat;
}



/*-----------------------------------------------------------------------------
 *  print functions octave
 *-----------------------------------------------------------------------------*/
void MESS_status::print(std::ostream & os, bool pr_as_read_syntax) const {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_status::print(std::ostream & os, bool pr_as_read_syntax) {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_status::print_raw(std::ostream & os, bool pr_as_read_syntax) const {
    os << MESS_status::static_type_name();
#ifdef MESS_DEBUG
    os << " ( " << this->ptr << " ) ";
#endif
}

bool MESS_status::print_name_tag (std::ostream& os, const std::string& name) const {
    os << name << " = ";
    return true;
}


/*-----------------------------------------------------------------------------
 *  overrides from octave_base_value
 *-----------------------------------------------------------------------------*/
bool MESS_status::is_constant(void) const { return true; }
bool MESS_status::is_defined(void) const { return true; }
bool MESS_status::print_as_scalar(void) const { return true; }
bool MESS_status::is_empty(void) const { return this->ptr?false:true; }
size_t MESS_status::byte_size(void) const{ return this->ptr?mess_status_memsize(this->ptr):0; }


/*-----------------------------------------------------------------------------
 *  conversion operator
 *-----------------------------------------------------------------------------*/
MESS_status::operator mess_status(void) const{ return this->ptr; }
MESS_status::operator bool(void) const{ return this->ptr != NULL; }


/*-----------------------------------------------------------------------------
 *  destructor
 *-----------------------------------------------------------------------------*/
MESS_status::~MESS_status(void){
    if(this->ptr){
        mess_status_clear(&(this->ptr));
    }
    this->ptr=NULL;
}
