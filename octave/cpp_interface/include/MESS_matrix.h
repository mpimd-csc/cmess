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

#ifndef MESS_MATRIX_OCTAVE_H_
#define MESS_MATRIX_OCTAVE_H_

#include <octave/ov-scalar.h>
#include <octave/oct.h>
#include <octave/ov-complex.h>
#include <octave/ov-re-sparse.h>
#include "mess/mess.h"



class MESS_matrix : public octave_base_value {

    public:

        /*-----------------------------------------------------------------------------
         *  constructors
         *-----------------------------------------------------------------------------*/
        MESS_matrix(void);
        MESS_matrix(mess_matrix);
        MESS_matrix(Matrix m);
        MESS_matrix(ComplexMatrix m);
        MESS_matrix(SparseMatrix m);
        MESS_matrix(SparseComplexMatrix m);
        MESS_matrix(double m);
        MESS_matrix(Complex m);
        MESS_matrix(const MESS_matrix & m);


        /*-----------------------------------------------------------------------------
         *  getter
         *-----------------------------------------------------------------------------*/
        mess_matrix get_ptr(void) const;


        /*-----------------------------------------------------------------------------
         *  deep copy
         *-----------------------------------------------------------------------------*/
        MESS_matrix deepcopy(void) const;

        /*-----------------------------------------------------------------------------
         *  print function for octave
         *-----------------------------------------------------------------------------*/
        void print(std::ostream & os, bool pr_as_read_syntax = false) const;
        void print(std::ostream & os, bool pr_as_read_syntax = false);
        void print_raw(std::ostream & os, bool pr_as_read_syntax) const;
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
        bool is_sparse_type(void) const;
        bool is_empty(void) const;
        dim_vector dims(void) const;
        int ndims(void) const;
        size_t length(void) const;
        size_t byte_size(void) const;
        octave_value full_value (void) const;
        octave_value map (unary_mapper_t umap) const;
        //bool save_ascii (std::ostream& os);

        Matrix matrix_value(bool = false) const;
        ComplexMatrix complex_matrix_value(bool = false) const;
        SparseMatrix sparse_matrix_value(bool = false) const;
        SparseComplexMatrix sparse_complex_matrix_value(bool = false) const;

        /*-----------------------------------------------------------------------------
         *  conversion operator
         *-----------------------------------------------------------------------------*/
        operator mess_matrix() const;
        operator bool() const;

        /*-----------------------------------------------------------------------------
         *  destructor
         *-----------------------------------------------------------------------------*/
        ~MESS_matrix(void);


        /*-----------------------------------------------------------------------------
         * installation routines for interpreter
         *-----------------------------------------------------------------------------*/
        static bool type_installed;
        static void install(void);

    private:

        mess_matrix ptr;

        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA;

};


/*-----------------------------------------------------------------------------
 *  binary operations octave_scalar <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
octave_value op_mul(const octave_scalar & s, const MESS_matrix & m);
octave_value op_mul(const MESS_matrix & m, const octave_scalar & s);


/*-----------------------------------------------------------------------------
 *  binary operations octave_complex_scalar <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
octave_value op_mul(const octave_complex_scalar & s, const MESS_matrix & m);
octave_value op_mul(const MESS_matrix & m, const octave_complex_scalar & s);


/*-----------------------------------------------------------------------------
 *  unary operations MESS_matrix
 *-----------------------------------------------------------------------------*/
octave_value op_not(const MESS_matrix & m);
octave_value op_uplus(const MESS_matrix & m);
octave_value op_uminus(const MESS_matrix & m);
octave_value op_transpose(const MESS_matrix & m);
octave_value op_hermitian(const MESS_matrix & m);

/*-----------------------------------------------------------------------------
 *  MACRO for generating binary operation TYPEA <-> TYPEB
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(TYPEA, TYPEB)                                                                             \
    octave_value op_add(const TYPEA & m1, const TYPEB &m2);                                                         \
    octave_value op_sub(const TYPEA & m1, const TYPEB &m2);                                                         \
    octave_value op_xmul(const mess_operation_t opA, TYPEA & m1, const mess_operation_t opB, const TYPEB &m2);      \
    octave_value op_mul(const TYPEA & m1, const TYPEB & m2);                                                        \
    octave_value op_trans_mul(const TYPEA & m1, const TYPEB & m2);                                                  \
    octave_value op_herm_mul(const TYPEA & m1, const TYPEB & m2);                                                   \
    octave_value op_mul_trans(const TYPEA & m1, const TYPEB & m2);                                                  \
    octave_value op_mul_herm(const TYPEA & m1, const TYPEB & m2);                                                   \
    octave_value op_xldiv(const TYPEA & m1, const mess_operation_t opA, const TYPEB & m2);                          \
    octave_value op_ldiv(const TYPEA & m1, const TYPEB & m2);                                                       \
    octave_value op_trans_ldiv(const TYPEA & m1, const TYPEB &m2);                                                  \
    octave_value op_herm_ldiv(const TYPEA & m1, const TYPEB &m2);                                                   \


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> MESS_matrix
 *-----------------------------------------------------------------------------*/
OCTMESS_GENERATOR(MESS_matrix, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_matrix
 *-----------------------------------------------------------------------------*/
inline mess_matrix_st octmess_matrix_conv(const Matrix & m){
    mess_matrix_st dummy = {
        /*rows*/            m.rows(),
        /*cols*/            m.cols(),
        /*ld*/              m.rows(),
        /*nnz*/             m.rows()*m.cols(),
        /*rowptr*/          NULL,
        /*colptr*/          NULL,
        /*values*/          const_cast<double*>(m.fortran_vec()),
        /*values_cpx*/      NULL,
        /*storage_type*/    MESS_DENSE,
        /*symmetry*/        MESS_GENERAL,
        /*data_type*/       MESS_REAL,
    };
    return dummy;
}

OCTMESS_GENERATOR(MESS_matrix, octave_matrix);
OCTMESS_GENERATOR(octave_matrix, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_complex_matrix
 *-----------------------------------------------------------------------------*/
inline mess_matrix_st octmess_matrix_conv(const ComplexMatrix & m){
    mess_matrix_st dummy = {
        /*rows*/            m.rows(),
        /*cols*/            m.cols(),
        /*ld*/              m.rows(),
        /*nnz*/             m.rows()*m.cols(),
        /*rowptr*/          NULL,
        /*colptr*/          NULL,
        /*values*/          NULL,
        /*values_cpx*/      reinterpret_cast <mess_double_cpx_t*>(const_cast<std::complex<double>*>(m.fortran_vec())),
        /*storage_type*/    MESS_DENSE,
        /*symmetry*/        MESS_GENERAL,
        /*data_type*/       MESS_COMPLEX,
    };
    return dummy;
}

OCTMESS_GENERATOR(MESS_matrix, octave_complex_matrix);
OCTMESS_GENERATOR(octave_complex_matrix, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_sparse_matrix
 *-----------------------------------------------------------------------------*/
inline mess_matrix_st octmess_matrix_conv(const SparseMatrix & m){
    mess_matrix_st dummy = {
        /*rows*/            m.rows(),
        /*cols*/            m.cols(),
        /*ld*/              m.rows(),
        /*nnz*/             m.nnz(),
        /*rowptr*/          m.ridx(),
        /*colptr*/          m.cidx(),
        /*values*/          const_cast<double*>(m.data()),
        /*values_cpx*/      NULL,
        /*storage_type*/    MESS_CSC,
        /*symmetry*/        MESS_GENERAL,
        /*data_type*/       MESS_REAL,
    };
    return dummy;
}

OCTMESS_GENERATOR(MESS_matrix, octave_sparse_matrix);
OCTMESS_GENERATOR(octave_sparse_matrix, MESS_matrix);


/*-----------------------------------------------------------------------------
 *  binary operations MESS_matrix <-> octave_sparse_complex_matrix
 *-----------------------------------------------------------------------------*/
inline mess_matrix_st octmess_matrix_conv(const SparseComplexMatrix & m){
    mess_matrix_st dummy = {
        /*rows*/            m.rows(),
        /*cols*/            m.cols(),
        /*ld*/              m.rows(),
        /*nnz*/             m.nnz(),
        /*rowptr*/          m.ridx(),
        /*colptr*/          m.cidx(),
        /*values*/          NULL,
        /*values_cpx*/      reinterpret_cast <mess_double_cpx_t*>(const_cast<std::complex<double>*>(m.data())),
        /*storage_type*/    MESS_CSC,
        /*symmetry*/        MESS_GENERAL,
        /*data_type*/       MESS_COMPLEX,
    };
    return dummy;
}

OCTMESS_GENERATOR(MESS_matrix, octave_sparse_complex_matrix);
OCTMESS_GENERATOR(octave_sparse_complex_matrix, MESS_matrix);


#endif
