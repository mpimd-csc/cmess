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

#ifndef MESS_ERROR_OCTAVE_H_
#define MESS_ERROR_OCTAVE_H_

#include <octave/oct.h>
#include "mess/mess.h"


#define OCTMESS_ERROR(RET,COND,FUN){                                                                                        \
    if((COND)){                                                                                                             \
        error("OCT-M.E.S.S.: " __FILE__ ":%d\n\t %s returned Error %d - %s",__LINE__,#FUN,(RET),mess_get_error((RET)));    \
    }                                                                                                                       \
}                                                                                                                           \

#endif
