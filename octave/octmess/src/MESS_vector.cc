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
#include "mess/mess.h"
#include "MESS_error.h"
#include "MESS_vector.h"



/*-----------------------------------------------------------------------------
 *  constructor
 *-----------------------------------------------------------------------------*/
DEFUN_DLD(mess_vector, args, nargout, "blabla"){
    if(args.length()==1){
        if(args(0).is_matrix_type()){
            if(args(0).is_double_type() && args(0).is_real_matrix()){
                return octave_value(new MESS_vector(args(0).matrix_value()));
            }else if(args(0).is_complex_type() && args(0).is_complex_matrix()){
                return octave_value(new MESS_vector(args(0).complex_matrix_value()));
            }else if(args(0).is_real_scalar()){
                return octave_value(new MESS_vector(args(0).double_value()));
            }else if(args(0).is_complex_scalar()){
                return octave_value(new MESS_vector(args(0).complex_value()));
            }
        }
    }
    error("Cannot create %s instance.",MESS_vector::static_type_name().c_str());
    return octave_value_list();
}

/*-----------------------------------------------------------------------------
 *  print
 *-----------------------------------------------------------------------------*/
#define OCTMESS_MESS_VECTOR_PRINT(PRINTFCT)                                                                                 \
DEFUN_DLD(PRINTFCT, args, nargout, "Interface to " #PRINTFCT){                                                              \
    int ret = 0;                                                                                                            \
    if(nargout || args.length()!=1 || args(0).type_id() != MESS_vector::static_type_id()){                                  \
        print_usage();                                                                                                      \
    }else{                                                                                                                  \
        ret = PRINTFCT((mess_vector)const_cast<MESS_vector &> (dynamic_cast<const MESS_vector &> (args(0).get_rep())));     \
        OCTMESS_ERROR(ret,(ret!=0),PRINTFCT);                                                                               \
    }                                                                                                                       \
    return octave_value_list();                                                                                             \
}

OCTMESS_MESS_VECTOR_PRINT(mess_vector_print)
OCTMESS_MESS_VECTOR_PRINT(mess_vector_printinfo)
OCTMESS_MESS_VECTOR_PRINT(mess_vector_printshort)



