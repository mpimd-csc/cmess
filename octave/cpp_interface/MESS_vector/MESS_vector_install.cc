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
#include <octave/ov-scalar.h>
#include <octave/ov-complex.h>
#include <octave/ov-re-sparse.h>
#include "MESS_vector.h"


/*-----------------------------------------------------------------------------
 *  unary operations MESS_vector
 *-----------------------------------------------------------------------------*/
DEFUNOP(mv_not,     MESS_vector){return op_not(      dynamic_cast<const MESS_vector&> (a));}
DEFUNOP(mv_uplus,   MESS_vector){return op_uplus(    dynamic_cast<const MESS_vector&> (a));}
DEFUNOP(mv_uminus,  MESS_vector){return op_uminus(   dynamic_cast<const MESS_vector&> (a));}


/*-----------------------------------------------------------------------------
 *  binary operations octave_scalar <-> MESS_vector
 *-----------------------------------------------------------------------------*/
DEFBINOP(os_mv_mul, octave_scalar, MESS_matrix){return op_mul(dynamic_cast<const octave_scalar&> (a1), dynamic_cast<const MESS_vector&> (a2));}
DEFBINOP(mv_os_mul, MESS_matrix, octave_scalar){return op_mul(dynamic_cast<const MESS_vector&> (a1), dynamic_cast<const octave_scalar&> (a2));}


/*-----------------------------------------------------------------------------
 *  binary operations octave_scalar <-> MESS_vector
 *-----------------------------------------------------------------------------*/
DEFBINOP(osc_mv_mul, octave_complex_scalar, MESS_matrix){return op_mul(dynamic_cast<const octave_complex_scalar&> (a1), dynamic_cast<const MESS_vector&> (a2));}
DEFBINOP(mv_osc_mul, MESS_matrix, octave_complex_scalar){return op_mul(dynamic_cast<const MESS_vector&> (a1), dynamic_cast<const octave_complex_scalar&> (a2));}


/*-----------------------------------------------------------------------------
 *  binary operations MESS_vector <-> MESS_vector
 *-----------------------------------------------------------------------------*/
DEFBINOP(mv_mv_add,         MESS_vector, MESS_vector){return op_add(        dynamic_cast<const MESS_vector&> (a1), dynamic_cast<const MESS_vector&> (a2));}
DEFBINOP(mv_mv_sub,         MESS_vector, MESS_vector){return op_sub(        dynamic_cast<const MESS_vector&> (a1), dynamic_cast<const MESS_vector&> (a2));}
DEFBINOP(mv_mv_trans_mul,   MESS_vector, MESS_vector){return op_trans_mul(  dynamic_cast<const MESS_vector&> (a1), dynamic_cast<const MESS_vector&> (a2));}
DEFBINOP(mv_mv_herm_mul,    MESS_vector, MESS_vector){return op_herm_mul(   dynamic_cast<const MESS_vector&> (a1), dynamic_cast<const MESS_vector&> (a2));}


/*-----------------------------------------------------------------------------
 *  MACRO for generating definitions binary operation TYPEA <-> TYPEB
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(SHORTTYPEA, TYPEA, SHORTTYPEB, TYPEB)                                                                                                     \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _mul,         TYPEA, TYPEB){return op_mul(        dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _trans_mul,   TYPEA, TYPEB){return op_trans_mul(  dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _herm_mul,    TYPEA, TYPEB){return op_herm_mul(   dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _ldiv,        TYPEA, TYPEB){return op_ldiv(       dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _trans_ldiv,  TYPEA, TYPEB){return op_trans_ldiv( dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _herm_ldiv,   TYPEA, TYPEB){return op_herm_ldiv(  dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mm, MESS_matrix, mv, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(om, octave_matrix, mv, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_complex_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(ocm, octave_complex_matrix, mv, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_sparse_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(osm, octave_sparse_matrix, mv, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_sparse_complex_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(oscm, octave_sparse_complex_matrix, mv, MESS_vector);


/*-----------------------------------------------------------------------------
 *  init the installed variable and define install function
 *-----------------------------------------------------------------------------*/
bool MESS_vector::type_installed = false;


/*-----------------------------------------------------------------------------
 *  MACRO for generating installation of binary operation TYPEA <-> TYPEB
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(SHORTTYPEA, TYPEA, SHORTTYPEB, TYPEB)                                 \
    INSTALL_BINOP(op_mul,           TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _mul);         \
    INSTALL_BINOP(op_trans_mul,     TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _trans_mul);   \
    INSTALL_BINOP(op_herm_mul,      TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _herm_mul);    \
    INSTALL_BINOP(op_ldiv,          TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _ldiv);        \
    INSTALL_BINOP(op_trans_ldiv,    TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _trans_ldiv);  \
    INSTALL_BINOP(op_herm_ldiv,     TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _herm_ldiv);   \



void MESS_vector::install(void){
    if(!MESS_vector::type_installed){
        MESS_vector::register_type();
        octave_stdout << "Install " << MESS_vector::static_type_name() << " at type-id = " << MESS_vector::static_type_id() << "\n";

        /*-----------------------------------------------------------------------------
         *  install unary operations
         *-----------------------------------------------------------------------------*/
        INSTALL_UNOP(op_not,        MESS_vector, mv_not);
        INSTALL_UNOP(op_uplus,      MESS_vector, mv_uplus);
        INSTALL_UNOP(op_uminus,     MESS_vector, mv_uminus);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_scalar <-> MESS_vector
         *-----------------------------------------------------------------------------*/
        INSTALL_BINOP(op_mul,           octave_scalar, MESS_vector, os_mv_mul);
        INSTALL_BINOP(op_el_mul,        octave_scalar, MESS_vector, os_mv_mul);

        INSTALL_BINOP(op_mul,           MESS_vector, octave_scalar, mv_os_mul);
        INSTALL_BINOP(op_el_mul,        MESS_vector, octave_scalar, mv_os_mul);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_complex_scalar <-> MESS_vector
         *-----------------------------------------------------------------------------*/
        INSTALL_BINOP(op_mul,           octave_complex_scalar, MESS_vector, osc_mv_mul);
        INSTALL_BINOP(op_el_mul,        octave_complex_scalar, MESS_vector, osc_mv_mul);

        INSTALL_BINOP(op_mul,           MESS_vector, octave_complex_scalar, mv_osc_mul);
        INSTALL_BINOP(op_el_mul,        MESS_vector, octave_complex_scalar, mv_osc_mul);

        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_vector <-> MESS_vector
         *-----------------------------------------------------------------------------*/
        INSTALL_BINOP(op_add,           MESS_vector, MESS_vector, mv_mv_add);
        INSTALL_BINOP(op_sub,           MESS_vector, MESS_vector, mv_mv_sub);
        INSTALL_BINOP(op_trans_mul,     MESS_vector, MESS_vector, mv_mv_trans_mul);
        INSTALL_BINOP(op_herm_mul,      MESS_vector, MESS_vector, mv_mv_herm_mul);

        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_matrix -> MESS_vector
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(mm, MESS_matrix, mv, MESS_vector);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_matrix -> MESS_vector
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(om, octave_matrix, mv, MESS_vector);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_complex_matrix -> MESS_vector
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(ocm, octave_complex_matrix, mv, MESS_vector);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_sparse_matrix -> MESS_vector
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(osm, octave_sparse_matrix, mv, MESS_vector);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_sparse_complex_matrix -> MESS_vector
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(oscm, octave_sparse_complex_matrix, mv, MESS_vector);


    }
    MESS_vector::type_installed = true;
}
