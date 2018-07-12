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
#include "MESS_matrix.h"

/*-----------------------------------------------------------------------------
 *  unary operations MESS_matrix
 *-----------------------------------------------------------------------------*/
DEFUNOP(mm_not,        MESS_matrix){return op_not(         dynamic_cast<const MESS_matrix&> (a));}
DEFUNOP(mm_uplus,      MESS_matrix){return op_uplus(       dynamic_cast<const MESS_matrix&> (a));}
DEFUNOP(mm_uminus,     MESS_matrix){return op_uminus(      dynamic_cast<const MESS_matrix&> (a));}
DEFUNOP(mm_transpose,  MESS_matrix){return op_transpose(   dynamic_cast<const MESS_matrix&> (a));}
DEFUNOP(mm_hermitian,  MESS_matrix){return op_hermitian(   dynamic_cast<const MESS_matrix&> (a));}


/*-----------------------------------------------------------------------------
 *  MACRO for generating definitions binary operation TYPEA <-> TYPEB
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(SHORTTYPEA, TYPEA, SHORTTYPEB, TYPEB)                                                                                                     \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _add,         TYPEA, TYPEB){return op_add(        dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _sub,         TYPEA, TYPEB){return op_sub(        dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _mul,         TYPEA, TYPEB){return op_mul(        dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _mul_trans,   TYPEA, TYPEB){return op_mul_trans(  dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _mul_herm,    TYPEA, TYPEB){return op_mul_herm(   dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _trans_mul,   TYPEA, TYPEB){return op_trans_mul(  dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _herm_mul,    TYPEA, TYPEB){return op_herm_mul(   dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _ldiv,        TYPEA, TYPEB){return op_ldiv(       dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _trans_ldiv,  TYPEA, TYPEB){return op_trans_ldiv( dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \
    DEFBINOP(SHORTTYPEA## _ ##SHORTTYPEB## _herm_ldiv,   TYPEA, TYPEB){return op_herm_ldiv(  dynamic_cast<const TYPEA&> (a1), dynamic_cast<const TYPEB&> (a2));}    \


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mm, MESS_matrix, mm, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_matrix
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mm, MESS_matrix, om, octave_matrix);
OCTMESS_GENERATOR(om, octave_matrix, mm, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_complex_matrix
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mm, MESS_matrix, ocm, octave_complex_matrix);
OCTMESS_GENERATOR(ocm, octave_complex_matrix, mm, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_sparse_matrix
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mm, MESS_matrix, osm, octave_sparse_matrix);
OCTMESS_GENERATOR(osm, octave_sparse_matrix, mm, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_sparse_complex_matrix
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(mm, MESS_matrix, oscm, octave_sparse_complex_matrix);
OCTMESS_GENERATOR(oscm, octave_sparse_complex_matrix, mm, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations octave_scalar <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
DEFBINOP(os_mm_mul, octave_scalar, MESS_matrix){return op_mul(dynamic_cast<const octave_scalar&> (a1), dynamic_cast<const MESS_matrix&> (a2));}
DEFBINOP(mm_os_mul, MESS_matrix, octave_scalar){return op_mul(dynamic_cast<const MESS_matrix&> (a1), dynamic_cast<const octave_scalar&> (a2));}


/*-----------------------------------------------------------------------------
 *  binary operations octave_complex_scalar <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
DEFBINOP(ocs_mm_mul, octave_complex_scalar, MESS_matrix){return op_mul(dynamic_cast<const octave_complex_scalar&> (a1), dynamic_cast<const MESS_matrix&> (a2));}
DEFBINOP(mm_ocs_mul, MESS_matrix, octave_complex_scalar){return op_mul(dynamic_cast<const MESS_matrix&> (a1), dynamic_cast<const octave_scalar&> (a2));}


/*-----------------------------------------------------------------------------
 *  init the installed variable and define install function
 *-----------------------------------------------------------------------------*/
bool MESS_matrix::type_installed = false;


/*-----------------------------------------------------------------------------
 *  MACRO for generating installation of binary operation TYPEA <-> TYPEB
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(SHORTTYPEA, TYPEA, SHORTTYPEB, TYPEB)                                 \
    INSTALL_BINOP(op_add,           TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _add);         \
    INSTALL_BINOP(op_sub,           TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _sub);         \
    INSTALL_BINOP(op_mul,           TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _mul);         \
    INSTALL_BINOP(op_trans_mul,     TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _trans_mul);   \
    INSTALL_BINOP(op_herm_mul,      TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _herm_mul);    \
    INSTALL_BINOP(op_mul_trans,     TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _mul_trans);   \
    INSTALL_BINOP(op_mul_herm,      TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _mul_herm);    \
    INSTALL_BINOP(op_ldiv,          TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _ldiv);        \
    INSTALL_BINOP(op_trans_ldiv,    TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _trans_ldiv);  \
    INSTALL_BINOP(op_herm_ldiv,     TYPEA, TYPEB,  SHORTTYPEA## _ ##SHORTTYPEB## _herm_ldiv);   \


void MESS_matrix::install(void){
    if(!MESS_matrix::type_installed){
        MESS_matrix::register_type();
        octave_stdout << "Install " << MESS_matrix::static_type_name() << " at type-id = " << MESS_matrix::static_type_id() << "\n";

        /*-----------------------------------------------------------------------------
         *  install unary operations
         *-----------------------------------------------------------------------------*/
        INSTALL_UNOP(op_not,        MESS_matrix, mm_not);
        INSTALL_UNOP(op_uplus,      MESS_matrix, mm_uplus);
        INSTALL_UNOP(op_uminus,     MESS_matrix, mm_uminus);
        INSTALL_UNOP(op_transpose,  MESS_matrix, mm_transpose);
        INSTALL_UNOP(op_hermitian,  MESS_matrix, mm_hermitian);


        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_matrix <-> MESS_matrix
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(mm, MESS_matrix, mm, MESS_matrix);

        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_matrix <-> octave_matrix
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(mm, MESS_matrix, om, octave_matrix);
        OCTMESS_GENERATOR(om, octave_matrix, mm, MESS_matrix);

        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_matrix <-> octave_complex_matrix
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(mm, MESS_matrix, ocm, octave_complex_matrix);
        OCTMESS_GENERATOR(ocm, octave_complex_matrix, mm, MESS_matrix);


        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_matrix <-> octave_sparse_matrix
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(mm, MESS_matrix, osm, octave_sparse_matrix);
        OCTMESS_GENERATOR(osm, octave_sparse_matrix, mm, MESS_matrix);


        /*-----------------------------------------------------------------------------
         *  install binary operations MESS_matrix <-> octave_sparse_complex_matrix
         *-----------------------------------------------------------------------------*/
        OCTMESS_GENERATOR(mm, MESS_matrix, oscm, octave_sparse_complex_matrix);
        OCTMESS_GENERATOR(oscm, octave_sparse_complex_matrix, mm, MESS_matrix);


        /*-----------------------------------------------------------------------------
         *  install binary operations octave_scalar <-> MESS_matrix
         *-----------------------------------------------------------------------------*/
        INSTALL_BINOP(op_mul,           octave_scalar, MESS_matrix, os_mm_mul);
        INSTALL_BINOP(op_el_mul,        octave_scalar, MESS_matrix, os_mm_mul);

        INSTALL_BINOP(op_mul,           MESS_matrix, octave_scalar, mm_os_mul);
        INSTALL_BINOP(op_el_mul,        MESS_matrix, octave_scalar, mm_os_mul);

        /*-----------------------------------------------------------------------------
         *  install binary operations octave_complex_scalar <-> MESS_matrix
         *-----------------------------------------------------------------------------*/
        INSTALL_BINOP(op_mul,           octave_complex_scalar, MESS_matrix, ocs_mm_mul);
        INSTALL_BINOP(op_el_mul,        octave_complex_scalar, MESS_matrix, ocs_mm_mul);

        INSTALL_BINOP(op_mul,           MESS_matrix, octave_complex_scalar, mm_ocs_mul);
        INSTALL_BINOP(op_el_mul,        MESS_matrix, octave_complex_scalar, mm_ocs_mul);

    }
    MESS_matrix::type_installed = true;
}



