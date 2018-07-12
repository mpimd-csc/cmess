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
#include <octave/ov-re-sparse.h>
#include "mess/mess.h"
#include "MESS_vector.h"
#include "MESS_error.h"

/*-----------------------------------------------------------------------------
 *  unary operations MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_not(const MESS_vector & v){
    int ret = 0;
    mess_vector nv = NULL;
    ret = mess_vector_init(&nv);                    OCTMESS_ERROR(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_copy((mess_vector)v,nv);      OCTMESS_ERROR(ret, (ret!=0), mess_vector_copy);
    ret = mess_vector_map_not(nv);                  OCTMESS_ERROR(ret, (ret!=0), mess_vector_map_not);
    return new MESS_vector(nv);
}

octave_value op_uplus (const MESS_vector & v){
    return new MESS_vector(v);
}


octave_value op_uminus (const MESS_vector & v){
    int ret = 0;
    mess_vector copy = NULL;
    ret = mess_vector_init(&copy);                      OCTMESS_ERROR(ret, (ret!=0), mess_vector_init);
    ret = mess_vector_copy((mess_vector)v, copy);       OCTMESS_ERROR(ret, (ret!=0), mess_vector_copy);
    ret = mess_vector_scale(-1.0, copy);                OCTMESS_ERROR(ret, (ret!=0), mess_vector_scale);
    return new MESS_vector(copy);
}


/*-----------------------------------------------------------------------------
 *  binary operations MESS_vector <-> MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_add (const MESS_vector & v1, const MESS_vector & v2) {
    mess_vector result=NULL;
    int ret = 0;
    if(v1 && v2){
        ret = mess_vector_init(&result);                            OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_copy((mess_vector)v2,result);             OCTMESS_ERROR(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_axpy(1.0,(mess_vector)v1,result);         OCTMESS_ERROR(ret,(ret!=0),mess_vector_axpy);
    }
    return new MESS_vector(result);
}


octave_value op_sub (const MESS_vector & v1, const MESS_vector & v2) {
    mess_vector result=NULL;
    int ret = 0;
    if(v1 && v2){
        ret = mess_vector_init(&result);                            OCTMESS_ERROR(ret,(ret!=0),mess_vector_init);
        ret = mess_vector_copy((mess_vector)v1,result);             OCTMESS_ERROR(ret,(ret!=0),mess_vector_copy);
        ret = mess_vector_axpy(-1.0,(mess_vector)v2,result);        OCTMESS_ERROR(ret,(ret!=0),mess_vector_axpy);
    }
    return new MESS_vector(result);
}


octave_value op_trans_mul (const MESS_vector & v1, const MESS_vector & v2){
    mess_double_cpx_t result;
    double *cpxptr = (double*) &result;
    int ret = 0;
    ret = mess_vector_dotu((mess_vector)v1,(mess_vector)v2,&result);    OCTMESS_ERROR(ret,(ret!=0),mess_vector_dotu);
    if(cpxptr[1] == 0.0){
        return cpxptr[0];
    }
    Complex ov {cpxptr[0], cpxptr[1]};
    return ov;
}


octave_value op_herm_mul (const MESS_vector & v1, const MESS_vector & v2){
    mess_double_cpx_t result;
    double *cpxptr = (double*) &result;
    int ret = 0;
    ret = mess_vector_dotc((mess_vector)v1,(mess_vector)v2,&result);    OCTMESS_ERROR(ret,(ret!=0),mess_vector_dotc);
    if(cpxptr[1] == 0.0){
        return cpxptr[0];
    }
    Complex ov {cpxptr[0], cpxptr[1]};
    return ov;
}
