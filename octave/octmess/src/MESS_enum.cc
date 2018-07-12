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
#include "MESS_enums.h"


/*-----------------------------------------------------------------------------
 *  MACRO for making enumeration types available in
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(ENUMT,ENUMVAL)                                            \
    DEFUN_DLD(ENUMVAL, args, nargout, "Representation of C-M.E.S.S. enum "#ENUMT){  \
        if(nargout > 1 || args.length()){                                           \
            print_usage();                                                          \
            return octave_value_list();                                             \
        }                                                                           \
        return octave_value(new ENUMT(ENUMVAL));                                    \
    }


OCTMESS_GENERATOR(MESS_storage_t, MESS_UNKNOWN);
OCTMESS_GENERATOR(MESS_storage_t, MESS_CSR);
OCTMESS_GENERATOR(MESS_storage_t, MESS_CSC);
OCTMESS_GENERATOR(MESS_storage_t, MESS_DENSE);
OCTMESS_GENERATOR(MESS_storage_t, MESS_COORD);

OCTMESS_GENERATOR(MESS_equation_t, MESS_EQN_NONE);
OCTMESS_GENERATOR(MESS_equation_t, MESS_EQN_LYAP);
OCTMESS_GENERATOR(MESS_equation_t, MESS_EQN_GLYAP);
OCTMESS_GENERATOR(MESS_equation_t, MESS_EQN_RICCATI);
OCTMESS_GENERATOR(MESS_equation_t, MESS_EQN_GRICCATI);

OCTMESS_GENERATOR(MESS_memusage_t, MESS_MEMORY_LOW);
OCTMESS_GENERATOR(MESS_memusage_t, MESS_MEMORY_MID);
OCTMESS_GENERATOR(MESS_memusage_t, MESS_MEMORY_HIGH);

OCTMESS_GENERATOR(MESS_operation_t, MESS_OP_NONE);
OCTMESS_GENERATOR(MESS_operation_t, MESS_OP_TRANSPOSE);
OCTMESS_GENERATOR(MESS_operation_t, MESS_OP_HERMITIAN);

OCTMESS_GENERATOR(MESS_parameter_t, MESS_LRCFADI_PARA_MINMAX);
OCTMESS_GENERATOR(MESS_parameter_t, MESS_LRCFADI_PARA_MINMAX_REAL);
OCTMESS_GENERATOR(MESS_parameter_t, MESS_LRCFADI_PARA_WACHSPRESS);
OCTMESS_GENERATOR(MESS_parameter_t, MESS_LRCFADI_PARA_ADAPTIVE_V);
OCTMESS_GENERATOR(MESS_parameter_t, MESS_LRCFADI_PARA_ADAPTIVE_Z);

OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_DEFAULT_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_SPARSE_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_LAPACK_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_UMFPACK_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_SUPERLU_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_CSPARSE_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_BANDED_LU);
OCTMESS_GENERATOR(MESS_direct_lupackage_t, MESS_DIRECT_MKLPARDISO_LU);

OCTMESS_GENERATOR(MESS_direct_cholpackage_t, MESS_DIRECT_DEFAULT_CHOLESKY);
OCTMESS_GENERATOR(MESS_direct_cholpackage_t, MESS_DIRECT_LAPACK_CHOLESKY);
OCTMESS_GENERATOR(MESS_direct_cholpackage_t, MESS_DIRECT_CSPARSE_CHOLESKY);
OCTMESS_GENERATOR(MESS_direct_cholpackage_t, MESS_DIRECT_CHOLMOD_CHOLESKY);

OCTMESS_GENERATOR(MESS_multidirect_t, MESS_MULTIDIRECT_SPARSE_LU);
OCTMESS_GENERATOR(MESS_multidirect_t, MESS_MULTIDIRECT_UMFPACK_LU);

OCTMESS_GENERATOR(MESS_residual_t, MESS_RESIDUAL_INDEFINITE);
OCTMESS_GENERATOR(MESS_residual_t, MESS_RESIDUAL_SPECTRAL);




