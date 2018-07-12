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
#include "MESS_vector.h"
#include "MESS_error.h"


/*-----------------------------------------------------------------------------
 *  binary operations octave_sclar <-> MESS_vector
 *-----------------------------------------------------------------------------*/

octave_value op_mul(const octave_scalar & s, const MESS_vector & v){
    mess_vector result = NULL;
    int ret = 0;
    if(v){
        ret = mess_vector_init(&result);                    OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_copy((mess_vector)v,result);      OCTMESS_ERROR(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_scale(s.double_value(),result);   OCTMESS_ERROR(ret,(ret!=0),mess_vector_scale);
    }
    return new MESS_vector(result);
}

octave_value op_mul(const MESS_vector & m, const octave_scalar &s){return op_mul(s,m);}

