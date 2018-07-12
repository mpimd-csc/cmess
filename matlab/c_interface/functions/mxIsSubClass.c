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

/**
 * @file matlab/c_interface/functions/mxIsSubClass.c
 * @brief
 * @author @mbehr
 */

#include "interface_matlab.h"

int mxIsSubClass (const mxArray  *instance, char* class){

    mxArray* m_temp;
    mxArray* m_eqn_name = mxCreateString(class);
    mxArray* myinstance = (mxArray*) instance;
    mxArray* call[2]= {myinstance, m_eqn_name};
    mexCallMATLAB(1,&m_temp,2,call,"isa");
    int isa_result = mxGetScalar(m_temp);

    return isa_result;

}


