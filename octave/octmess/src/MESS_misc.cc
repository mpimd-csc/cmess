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
#include "MESS_direct.h"
#include "MESS_enums.h"
#include "MESS_equation.h"
#include "MESS_matrix.h"
#include "MESS_multidirect.h"
#include "MESS_options.h"
#include "MESS_status.h"
#include "MESS_vector.h"


DEFUN_DLD(mess_init, args, nargout, "Initializes OCT-M.E.S.S."){
    if(args.length()){
        print_usage();
        return octave_value_list();
    }

    /*-----------------------------------------------------------------------------
     *  install basic types for cmess
     *-----------------------------------------------------------------------------*/
    MESS_direct::install();
    MESS_equation::install();
    MESS_matrix::install();
    MESS_multidirect::install();
    MESS_options::install();
    MESS_status::install();
    MESS_vector::install();

    /*-----------------------------------------------------------------------------
     *  install enums
     *-----------------------------------------------------------------------------*/
    MESS_datatype_t::install();
    MESS_direct_cholpackage_t::install();
    MESS_direct_lupackage_t::install();
    MESS_equation_t::install();
    MESS_memusage_t::install();
    MESS_multidirect_t::install();
    MESS_operation_t::install();
    MESS_parameter_t::install();
    MESS_residual_t::install();
    MESS_storage_t::install();

    return octave_value_list();
}
/*
%!test mess_init
%!error a = mess_init
*/



DEFUN_DLD(mess_eps, args, nargout, "Interface to mess_eps."){
    if(args.length()){
        print_usage();
        return octave_value_list();
    }
    return octave_value_list(octave_value(mess_eps()));
}
/*
%!assert(0 < mess_eps)
%!assert(mess_eps < 1e-12)
%!assert(isreal(mess_eps))
%!assert(isscalar(mess_eps))
%!error mess_error(1)
*/


DEFUN_DLD(mess_version, args, nargout, "Interface to mess_version."){
    if(args.length()){
        print_usage();
        return octave_value_list();
    }
    mess_version();
    return octave_value_list();
}
/*
%!test mess_version
%!error a = mess_version
*/






