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
 * @file matlab/c_interface/functions/mess_classname_from_mexmess.c
 * @brief
 * @author @mbehr
 */

#include "interface_matlab.h"

/**
 * @brief Get a class name from @matlab object.
 * @param[in] instance input instance of class
 * @return pointer of charaters with class name or @c NULL if failure occurs.
 *
 * The @ref mess_classname_from_mexmess get the classname from a @matlab instance.
 *
 */
char* mess_classname_from_mexmess (const mxArray  *instance){


    /*-----------------------------------------------------------------------------
     *  check input arguments
     *-----------------------------------------------------------------------------*/
    if(!instance){
        csc_error_message("Input instance points to NULL.\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  call class
     *-----------------------------------------------------------------------------*/
    mxArray* classname = NULL;
    mxArray* myinstance = (mxArray*) instance;
    //mexCallMATLAB(1,&classname,1,&instance,"class");
    mexCallMATLAB(1,&classname,1,&myinstance,"class");
    if(!classname){
        csc_error_message("mexCallMATLAB class failed.\n");
        return NULL;
    }

    /*-----------------------------------------------------------------------------
     *  get class string
     *-----------------------------------------------------------------------------*/
    return mxArrayToString(classname);

}


