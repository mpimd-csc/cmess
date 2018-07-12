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

#ifndef MESS_VECTOR_OCTAVE_H_
#define MESS_VECTOR_OCTAVE_H_

#include <octave/oct.h>
#include "mess/mess.h"
#include "MESS_matrix.h"


class MESS_vector : public octave_base_value {

    public:

        /*-----------------------------------------------------------------------------
         *  constructors
         *-----------------------------------------------------------------------------*/
        MESS_vector(void);
        MESS_vector(mess_vector v);
        MESS_vector(Matrix m);
        MESS_vector(ComplexMatrix m);
        MESS_vector(double m);
        MESS_vector(Complex m);
        MESS_vector(const MESS_vector & v);


        /*-----------------------------------------------------------------------------
         *  getter
         *-----------------------------------------------------------------------------*/
        mess_vector get_ptr(void) const;


        /*-----------------------------------------------------------------------------
         *  deep copy
         *-----------------------------------------------------------------------------*/
        MESS_vector deepcopy(void) const;


        /*-----------------------------------------------------------------------------
         *  print function for octave
         *-----------------------------------------------------------------------------*/
        void print(std::ostream & os, bool pr_as_read_syntax = false) const;
        void print(std::ostream & os, bool pr_as_read_syntax = false);
        void print_raw(std::ostream& os, bool pr_as_read_syntax) const;
        bool print_name_tag(std::ostream & os, const std::string & name) const;


        /*-----------------------------------------------------------------------------
         * overrides for octave
         *-----------------------------------------------------------------------------*/
        bool is_constant(void) const;
        bool is_defined(void) const;
        bool print_as_scalar(void) const;
        bool is_matrix_type(void) const;
        bool is_numeric_type(void) const;
        bool is_real_matrix(void) const;
        bool is_real_type(void) const;
        bool is_double_type(void) const;
        bool is_single_type(void) const;
        bool is_complex_matrix(void) const;
        bool is_complex_type(void) const;
        bool is_empty(void) const;
        dim_vector dims(void) const;
        int ndims(void) const;
        size_t length(void) const;
        size_t byte_size(void) const;

        octave_value map (unary_mapper_t umap) const;


        /*-----------------------------------------------------------------------------
         *  conversion operator
         *-----------------------------------------------------------------------------*/
        operator mess_vector() const;
        operator bool() const;

        /*-----------------------------------------------------------------------------
         *  destructor
         *-----------------------------------------------------------------------------*/
        ~MESS_vector(void);


        /*-----------------------------------------------------------------------------
         * installation routines for interpreter
         *-----------------------------------------------------------------------------*/
        static bool type_installed;
        static void install(void);

    private:

        mess_vector ptr;

        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA;

};

/*-----------------------------------------------------------------------------
 *  binary operations octave_scalar <-> MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_mul(const octave_scalar & s, const MESS_vector & v);
octave_value op_mul(const MESS_vector & v, const octave_scalar & s);


/*-----------------------------------------------------------------------------
 *  binary operations octave_complex_scalar <-> MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_mul(const octave_complex_scalar & s, const MESS_vector & v);
octave_value op_mul(const MESS_vector & v, const octave_complex_scalar & s);


/*-----------------------------------------------------------------------------
 *  unary operations MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_not(const MESS_vector & v);
octave_value op_uplus(const MESS_vector & v);
octave_value op_uminus(const MESS_vector & v);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_vector <-> MESS_vector
 *-----------------------------------------------------------------------------*/
octave_value op_add(const MESS_vector & v1, const MESS_vector & v2);
octave_value op_sub(const MESS_vector & v1, const MESS_vector & v2);
octave_value op_trans_mul(const MESS_vector & v1, const MESS_vector & v2);
octave_value op_herm_mul(const MESS_vector & v1, const MESS_vector & v2);


/*-----------------------------------------------------------------------------
 *  MACRO for generating binary operation TYPEA <-> TYPEB
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(TYPEA, TYPEB)                                                                             \
    octave_value op_xmul(const mess_operation_t opA, TYPEA & m1, const TYPEB &m2);                                  \
    octave_value op_mul(const TYPEA & m1, const TYPEB & m2);                                                        \
    octave_value op_trans_mul(const TYPEA & m1, const TYPEB & m2);                                                  \
    octave_value op_herm_mul(const TYPEA & m1, const TYPEB & m2);                                                   \
    octave_value op_xldiv(const TYPEA & m1, const mess_operation_t opA, const TYPEB & m2);                          \
    octave_value op_ldiv(const TYPEA & m1, const TYPEB & m2);                                                       \
    octave_value op_trans_ldiv(const TYPEA & m1, const TYPEB &m2);                                                  \
    octave_value op_herm_ldiv(const TYPEA & m1, const TYPEB &m2);                                                   \


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(MESS_matrix, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(octave_matrix, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_complex_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(octave_complex_matrix, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_sparse_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(octave_sparse_matrix, MESS_vector);


/*-----------------------------------------------------------------------------
 *  binary operations octave_sparse_complex_matrix -> MESS_vector
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(octave_sparse_complex_matrix, MESS_vector);



#endif







