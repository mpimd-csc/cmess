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
#include "MESS_options.h"



/*-----------------------------------------------------------------------------
 *  constructor
 *-----------------------------------------------------------------------------*/
DEFUN_DLD(mess_options, args, nargout, "blabla"){

    if(args.length()){
        print_usage();
        return octave_value_list();
    }


    /*-----------------------------------------------------------------------------
     *  create options instance
     *-----------------------------------------------------------------------------*/
    int ret = 0;
    mess_options opt;
    ret = mess_options_init(&opt);                  OCTMESS_ERROR(ret,(ret=0),mess_options_init);
    return octave_value(new MESS_options(opt));
}



/*-----------------------------------------------------------------------------
 *  print functions
 *-----------------------------------------------------------------------------*/
DEFUN_DLD(mess_options_print, args, nargout, "blabla"){

    if(nargout || args.length()!=1 || args(0).type_id() != MESS_options::static_type_id()){
        print_usage();
    }

    int ret = 0;
    ret = mess_options_print((dynamic_cast<const MESS_options &> (args(0).get_rep())));    OCTMESS_ERROR(ret,(ret!=0),mess_options_print);
    return octave_value_list();
}
