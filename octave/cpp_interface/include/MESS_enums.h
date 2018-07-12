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

#ifndef MESS_ENUMS_OCTAVE_H_
#define MESS_ENUMS_OCTAVE_H_

#include <octave/oct.h>
#include "mess/mess.h"



#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(OCTMESSENUM,MESSENUM,MESSENUMINIT,PRINTFCT)                                                               \
    class OCTMESSENUM : public octave_base_value {                                                                                  \
                                                                                                                                    \
        public:                                                                                                                     \
            /*-----------------------------------------------------------------------------                                         \
            *  constructors                                                                                                         \
            *-----------------------------------------------------------------------------*/                                        \
            OCTMESSENUM(void):m_enum(MESSENUMINIT){ };                                                                              \
            OCTMESSENUM(MESSENUM enu):m_enum(enu){ };                                                                               \
                                                                                                                                    \
            /*-----------------------------------------------------------------------------                                         \
             *  print function for octave                                                                                           \
             *-----------------------------------------------------------------------------*/                                       \
            void print(std::ostream & os, bool pr_as_read_syntax = false) const{print_raw(os, pr_as_read_syntax); newline(os);}     \
            void print(std::ostream & os, bool pr_as_read_syntax = false){print_raw(os, pr_as_read_syntax); newline(os);}           \
            void print_raw(std::ostream & os, bool pr_as_read_syntax) const{os<<#MESSENUM<<": "<<PRINTFCT(this->m_enum);}           \
                                                                                                                                    \
            /*-----------------------------------------------------------------------------                                         \
             * overrides for octave                                                                                                 \
             *-----------------------------------------------------------------------------*/                                       \
            bool is_constant(void) const { return true; };                                                                          \
            bool is_defined(void) const { return true; };                                                                           \
            bool print_as_scalar(void) { return true; };                                                                            \
            bool is_empty(void) const { return true; };                                                                             \
            size_t byte_size(void) const { return sizeof(MESSENUM);};                                                               \
                                                                                                                                    \
            /*-----------------------------------------------------------------------------                                         \
             * override == and != operator                                                                                          \
             *-----------------------------------------------------------------------------*/                                       \
            bool operator ==(const OCTMESSENUM &b) const { return this->m_enum == b.m_enum;}                                        \
            bool operator !=(const OCTMESSENUM &b) const { return this->m_enum != b.m_enum;}                                        \
                                                                                                                                    \
            /*-----------------------------------------------------------------------------                                         \
             *  conversion operator                                                                                                 \
             *-----------------------------------------------------------------------------*/                                       \
            operator MESSENUM() const {return this->m_enum;};                                                                       \
                                                                                                                                    \
            /*-----------------------------------------------------------------------------                                         \
             * installation routines for interpreter                                                                                \
             *-----------------------------------------------------------------------------*/                                       \
            static bool type_installed;                                                                                             \
            static void install(void);                                                                                              \
                                                                                                                                    \
        private:                                                                                                                    \
            MESSENUM m_enum;                                                                                                        \
            DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA;                                                                                   \
};                                                                                                                                  \
                                                                                                                                    \
/*-----------------------------------------------------------------------------                                                     \
 *  binary operations OCTMESSENUM <-> OCTMESSENUM                                                                                   \
 *-----------------------------------------------------------------------------*/                                                   \
octave_value op_eq(const OCTMESSENUM & o1, const OCTMESSENUM & o2);                                                                 \
octave_value op_ne(const OCTMESSENUM & o1, const OCTMESSENUM & o2);                                                                 \


OCTMESS_GENERATOR(MESS_storage_t,               mess_storage_t,             MESS_UNKNOWN,                   mess_storage_t_str);
OCTMESS_GENERATOR(MESS_datatype_t,              mess_datatype_t,            MESS_REAL,                      mess_datatype_t_str);
OCTMESS_GENERATOR(MESS_equation_t,              mess_equation_t,            MESS_EQN_NONE,                  mess_equation_t_str);
OCTMESS_GENERATOR(MESS_memusage_t,              mess_memusage_t,            MESS_MEMORY_LOW,                mess_memusage_t_str);
OCTMESS_GENERATOR(MESS_operation_t,             mess_operation_t,           MESS_OP_NONE,                   mess_operation_t_str);
OCTMESS_GENERATOR(MESS_parameter_t,             mess_parameter_t,           MESS_LRCFADI_PARA_MINMAX,       mess_parameter_t_str);
OCTMESS_GENERATOR(MESS_direct_lupackage_t,      mess_direct_lupackage_t,    MESS_DIRECT_DEFAULT_LU,         mess_direct_lupackage_t_str);
OCTMESS_GENERATOR(MESS_direct_cholpackage_t,    mess_direct_cholpackage_t,  MESS_DIRECT_DEFAULT_CHOLESKY,   mess_direct_cholpackage_t_str);
OCTMESS_GENERATOR(MESS_multidirect_t,           mess_multidirect_t,         MESS_MULTIDIRECT_SPARSE_LU,     mess_multidirect_t_str);
OCTMESS_GENERATOR(MESS_residual_t,              mess_residual_t,            MESS_RESIDUAL_INDEFINITE,       mess_residual_t_str);

#endif

