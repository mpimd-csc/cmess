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
#include <octave/ov-complex.h>
#include <octave/ov-re-sparse.h>
#include <math.h>
#include "mess/mess.h"
#include "MESS_enums.h"




#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(SHORTNAME,OCTMESSENUM,CMESSENUM)                                                                                                              \
    DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (OCTMESSENUM, #CMESSENUM, #CMESSENUM);                                                                                          \
                                                                                                                                                                        \
    /*-----------------------------------------------------------------------------                                                                                     \
     * implement binary operations                                                                                                                                      \
     *-----------------------------------------------------------------------------*/                                                                                   \
    octave_value op_eq(const OCTMESSENUM & o1, const OCTMESSENUM & o2){return octave_value(o1==o2);}                                                                    \
    octave_value op_ne(const OCTMESSENUM & o1, const OCTMESSENUM & o2){return octave_value(o1!=o2);}                                                                    \
    DEFBINOP(SHORTNAME##_##SHORTNAME##_eq, OCTMESSENUM, OCTMESSENUM){return op_eq(dynamic_cast<const OCTMESSENUM & > (a1), dynamic_cast<const OCTMESSENUM & > (a2));}   \
    DEFBINOP(SHORTNAME##_##SHORTNAME##_ne, OCTMESSENUM, OCTMESSENUM){return op_ne(dynamic_cast<const OCTMESSENUM & > (a1), dynamic_cast<const OCTMESSENUM & > (a2));}   \
                                                                                                                                                                        \
    bool OCTMESSENUM::type_installed = false;                                                                                                                           \
    void OCTMESSENUM::install(void){                                                                                                                                    \
        if(!OCTMESSENUM::type_installed){                                                                                                                               \
            OCTMESSENUM::register_type();                                                                                                                               \
            octave_stdout << "Install " << OCTMESSENUM::static_type_name()                                                                                              \
                          << " at type-id = " <<  OCTMESSENUM::static_type_id() << "\n";                                                                                \
            INSTALL_BINOP(op_eq, OCTMESSENUM, OCTMESSENUM,  SHORTNAME## _ ##SHORTNAME## _eq);                                                                           \
            INSTALL_BINOP(op_ne, OCTMESSENUM, OCTMESSENUM,  SHORTNAME## _ ##SHORTNAME## _ne);                                                                           \
                                                                                                                                                                        \
        }                                                                                                                                                               \
        OCTMESSENUM::type_installed = true;                                                                                                                             \
    }                                                                                                                                                                   \


/*-----------------------------------------------------------------------------
 *  implement enums
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mst,MESS_storage_t,mess_storage_t);
OCTMESS_GENERATOR(mdt,MESS_datatype_t,mess_datatype_t);
OCTMESS_GENERATOR(met,MESS_equation_t,mess_equation_t);
OCTMESS_GENERATOR(mmt,MESS_memusage_t,mess_memusage_t);
OCTMESS_GENERATOR(mot,MESS_operation_t,mess_operation_t);
OCTMESS_GENERATOR(mpt,MESS_parameter_t,mess_parameter_t);
OCTMESS_GENERATOR(mdlt,MESS_direct_lupackage_t,mess_direct_lupackage_t);
OCTMESS_GENERATOR(mdct,MESS_direct_cholpackage_t,mess_direct_cholpackage_t);
OCTMESS_GENERATOR(mmdt,MESS_multidirect_t,mess_multidirect_t);
OCTMESS_GENERATOR(mrt,MESS_residual_t,mess_residual_t);



