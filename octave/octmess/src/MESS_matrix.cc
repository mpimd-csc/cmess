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
#include "MESS_matrix.h"



/*-----------------------------------------------------------------------------
 *  constructor
 *-----------------------------------------------------------------------------*/
DEFUN_DLD(mess_matrix, args, nargout, "blabla"){

    if(nargout != 1 || args.length() != 1){
        print_usage();
        return octave_value_list();
    }

    //printf("t-id=%d\n",args(0).type_id());

    if(args.length()==1){
        if(args(0).is_sparse_type()){
            if(args(0).is_double_type() && args(0).is_real_matrix()){
                return octave_value(new MESS_matrix(args(0).sparse_matrix_value()));
            }else if(args(0).is_complex_type() && args(0).is_complex_matrix()){
                return octave_value(new MESS_matrix(args(0).sparse_complex_matrix_value()));
            }
        }else{
            if(args(0).is_double_type() && args(0).is_real_matrix()){
                return octave_value(new MESS_matrix(args(0).matrix_value()));
            }else if(args(0).is_complex_type() && args(0).is_complex_matrix()){
                return octave_value(new MESS_matrix(args(0).complex_matrix_value()));
            }else if(args(0).is_real_scalar()){
                return octave_value(new MESS_matrix(args(0).double_value()));
            }else if(args(0).is_complex_scalar()){
                return octave_value(new MESS_matrix(args(0).complex_value()));
            }
        }
    }

    error("Cannot create %s instance.",MESS_matrix::static_type_name().c_str());
    return octave_value_list();
}

/*
%!shared a, b
%!test
%!  for nrows = 1:5
%!      for ncols = 1:5
%!          for cpx = 0:1
%!              for sp= 0:1
%!                  for p = linspace(0,1,6)
%!                      if sp
%!                          a = sprand(nrows,ncols,p) + cpx*I*sprand(nrows,ncols,p);
%!                      else
%!                          a = rand(nrows,ncols) + cpx*I*rand(nrows,ncols);
%!                      end
%!                      b = mess_matrix(a);
%!                      assert(rows(a)==rows(b) && columns(a)==columns(b) && isreal(a)==isreal(b) && iscomplex(a)==iscomplex(b) && issparse(a)==issparse(b))
%!                      assert(all(size(a)==size(b)) && strcmp(typeinfo(b),"mess_matrix") && 0<sizeof(b) && isnumeric(b))
%!                  end
%!              end
%!          end
%!      end
%!  end
%!
*/

/*-----------------------------------------------------------------------------
 *  print
 *-----------------------------------------------------------------------------*/
#define OCTMESS_MESS_MATRIX_PRINT(PRINTFCT)                                                                                 \
DEFUN_DLD(PRINTFCT, args, nargout, "Interface to " #PRINTFCT){                                                              \
    int ret = 0;                                                                                                            \
    if(nargout || args.length()!=1 || args(0).type_id() != MESS_matrix::static_type_id()){                                  \
        print_usage();                                                                                                      \
    }else{                                                                                                                  \
        ret = PRINTFCT((mess_matrix)const_cast<MESS_matrix &> (dynamic_cast<const MESS_matrix &> (args(0).get_rep())));     \
        OCTMESS_ERROR(ret,(ret!=0),PRINTFCT);                                                                               \
    }                                                                                                                       \
    return octave_value_list();                                                                                             \
}

OCTMESS_MESS_MATRIX_PRINT(mess_matrix_print)
OCTMESS_MESS_MATRIX_PRINT(mess_matrix_printinfo)
OCTMESS_MESS_MATRIX_PRINT(mess_matrix_printshort)
OCTMESS_MESS_MATRIX_PRINT(mess_matrix_printdata)

