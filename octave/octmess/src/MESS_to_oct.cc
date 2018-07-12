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
#include "MESS_matrix.h"
#include "MESS_vector.h"

/*-----------------------------------------------------------------------------
 *  constructor
 *-----------------------------------------------------------------------------*/
DEFUN_DLD(mess_to_oct, args, nargout, "blabla"){

    if(nargout > 1 || args.length() != 1){
        print_usage();
        return octave_value_list();
    }

    /*-----------------------------------------------------------------------------
     *  mess_matrix to octave
     *-----------------------------------------------------------------------------*/
    octave_value_list ret;
    if(args(0).type_id() == MESS_matrix::static_type_id()){
        if(!args(0).is_sparse_type()){
            if(args(0).is_real_type()){
                ret(0) = args(0).matrix_value();
            }else{
                ret(0) = args(0).complex_matrix_value();
            }
        }else{
            if(args(0).is_real_type()){
                ret(0) = args(0).sparse_matrix_value();
            }else{
                ret(0) = args(0).sparse_complex_matrix_value();
            }
        }
        return ret;
    }else if (args(0).type_id() == MESS_vector::static_type_id()){

    }
    error("Cannot convert object %s with type-id %d to octave",args(0).type_name().c_str(),args(0).type_id());
    return octave_value_list();
}
