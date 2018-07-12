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
#include "MESS_matrix.h"
#include "MESS_error.h"


DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (MESS_matrix, "mess_matrix", "mess_matrix");


/*-----------------------------------------------------------------------------
 *  constructors
 *-----------------------------------------------------------------------------*/
MESS_matrix::MESS_matrix(void):ptr(NULL){ }
MESS_matrix::MESS_matrix(mess_matrix m):ptr(m) { }

MESS_matrix::MESS_matrix(Matrix m){
    mess_matrix_init(&(this->ptr));
    mess_matrix_dense_from_farray(this->ptr, m.rows(), m.cols(), m.rows(), m.fortran_vec(), NULL);
}


MESS_matrix::MESS_matrix(ComplexMatrix m){
    mess_matrix_init(&(this->ptr));
    mess_matrix_dense_from_farray(this->ptr, m.rows(), m.cols(), m.rows(), NULL, reinterpret_cast <mess_double_cpx_t*>(m.fortran_vec()));
}


MESS_matrix::MESS_matrix(double m){
    mess_matrix_init(&(this->ptr));
    mess_matrix_alloc(this->ptr, 1, 1, 1, MESS_DENSE, MESS_REAL);
    this->ptr->values[0]=m;
}


MESS_matrix::MESS_matrix(Complex m){
    double * ptr;
    mess_matrix_init(&(this->ptr));
    mess_matrix_alloc(this->ptr, 1, 1, 1, MESS_DENSE, MESS_COMPLEX);
    ptr = (double *) this->ptr->values_cpx;
    ptr[0] = m.real();
    ptr[1] = m.imag();
}


MESS_matrix::MESS_matrix(SparseMatrix m){
    mess_matrix_init(&(this->ptr));
    mess_matrix_csc(this->ptr, m.rows(), m.cols(), m.ridx(), m.cidx(), m.data(), NULL);
}


MESS_matrix::MESS_matrix(SparseComplexMatrix m){
    mess_matrix_init(&(this->ptr));
    mess_matrix_csc(this->ptr, m.rows(), m.cols(), m.ridx(), m.cidx(), NULL, reinterpret_cast<mess_double_cpx_t *>(m.data()));
}


MESS_matrix::MESS_matrix(const MESS_matrix & mat): ptr(mat.ptr){ }


/*-----------------------------------------------------------------------------
 *  getter
 *-----------------------------------------------------------------------------*/
mess_matrix MESS_matrix::get_ptr(void) const{ return this->ptr; }


/*-----------------------------------------------------------------------------
 *  deep copy
 *-----------------------------------------------------------------------------*/
MESS_matrix MESS_matrix::deepcopy(void) const{
    MESS_matrix m;
    if(*this){
        mess_matrix ptr = m.get_ptr();
        mess_matrix_init(&ptr);
        mess_matrix_copy(this->ptr,ptr);
    }
    return m;
}


/*-----------------------------------------------------------------------------
 *  print functions octave
 *-----------------------------------------------------------------------------*/
void MESS_matrix::print(std::ostream & os, bool pr_as_read_syntax) const {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_matrix::print(std::ostream & os, bool pr_as_read_syntax) {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_matrix::print_raw(std::ostream & os, bool pr_as_read_syntax) const {
    os << MESS_matrix::static_type_name();
#ifdef MESS_DEBUG
    os << " ( " << this->ptr << " ) ";
#endif
}

bool MESS_matrix::print_name_tag (std::ostream& os, const std::string& name) const {
    os << name << " = ";
    return true;
}


/*-----------------------------------------------------------------------------
 *  overrides from octave_base_value
 *-----------------------------------------------------------------------------*/
bool MESS_matrix::is_constant(void) const { return true; }
bool MESS_matrix::is_defined(void) const { return true; }
bool MESS_matrix::print_as_scalar(void) const { return true; }
bool MESS_matrix::is_matrix_type(void) const { return true; }
bool MESS_matrix::is_numeric_type(void) const { return true; }
bool MESS_matrix::is_real_matrix(void) const { return this->ptr?MESS_IS_REAL(this->ptr):false; }
bool MESS_matrix::is_real_type(void) const { return this->ptr?MESS_IS_REAL(this->ptr):false; }
bool MESS_matrix::is_double_type(void) const { return true; };
bool MESS_matrix::is_single_type(void) const { return false; }
bool MESS_matrix::is_complex_matrix (void) const { return this->ptr?MESS_IS_COMPLEX(this->ptr):false; }
bool MESS_matrix::is_complex_type (void) const { return this->ptr?MESS_IS_COMPLEX(this->ptr):false; }
bool MESS_matrix::is_sparse_type(void) const { return this->ptr?MESS_IS_SPARSE(this->ptr):false; }
bool MESS_matrix::is_empty(void) const { return this->ptr?false:true; }

size_t MESS_matrix::byte_size(void) const {
    if(this->ptr && (MESS_IS_REAL(this->ptr) || MESS_IS_COMPLEX(this->ptr))){
        return mess_matrix_memsize(this->ptr);
    }
    return 0;
}


dim_vector MESS_matrix::dims (void) const { return this->ptr?dim_vector(this->ptr->rows,this->ptr->cols):dim_vector(0,0); }
int MESS_matrix::ndims (void) const { return 2; }
size_t MESS_matrix::length(void) const{ return this->ptr?MESS_MAX(this->ptr->rows,this->ptr->cols):0; }

Matrix MESS_matrix::matrix_value(bool) const{
    if(this->ptr == NULL){
        error("mess_matrix points to NULL.");
        return Matrix();
    }

    if(this->is_real_matrix() && !(this->is_sparse_type())){
        Matrix om(this->dims());
        memcpy(om.fortran_vec(), this->ptr->values, sizeof(double)*this->ptr->nnz);
        return om;
    }

    error("Cannot convert mess_matrix to Matrix");
    return Matrix();
}

ComplexMatrix MESS_matrix::complex_matrix_value(bool) const{
    if(this->ptr == NULL){
        error("mess_matrix points to NULL.");
        return ComplexMatrix();
    }

    if(this->is_complex_matrix() && !(this->is_sparse_type())){
        ComplexMatrix ocm(this->dims());
        memcpy(reinterpret_cast <mess_double_cpx_t*>(ocm.fortran_vec()), this->ptr->values_cpx, sizeof(mess_double_cpx_t)*this->ptr->nnz);
        return ocm;
    }

    error("Cannot convert mess_matrix to ComplexMatrix");
    return ComplexMatrix();
}

SparseMatrix MESS_matrix::sparse_matrix_value(bool) const{
    if(this->ptr == NULL){
        error("mess_matrix points to NULL.");
        return SparseMatrix();
    }

    if(this->is_real_matrix() && this->is_sparse_type()){
        SparseMatrix osm(this->dims(),this->ptr->nnz);
        memcpy(osm.data(), this->ptr->values, sizeof(double)*this->ptr->nnz);
        return osm;
    }

    error("Cannot convert mess_matrix to SparseMatrix");
    return SparseMatrix();
}





/*-----------------------------------------------------------------------------
 *  overrides for built in functions
 *-----------------------------------------------------------------------------*/
octave_value MESS_matrix::full_value(void) const{
    int ret = 0;
    mess_matrix full = NULL;
    if(this->ptr){
        mess_matrix_init(&full);
        ret = mess_matrix_convert(this->ptr,full,MESS_DENSE);     OCTMESS_ERROR(ret,(ret!=0),mess_matrix_convert);
    }
    return new MESS_matrix(full);
}

#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(MAPNAME)                                                                              \
    case umap_##MAPNAME:                                                                                        \
        ret = mess_matrix_copy(this->ptr, result);  OCTMESS_ERROR(ret,(ret!=0),mess_matrix_copy);               \
        ret = mess_matrix_map_##MAPNAME(result, 1); OCTMESS_ERROR(ret,(ret!=0),mess_matrix_map_##MAPNAME);      \
        return new MESS_matrix(result);                                                                         \

octave_value MESS_matrix::map (unary_mapper_t umap) const{

    int ret = 0;
    mess_matrix result = NULL;
    if(!(*this)){
        error("%s instance points to NULL.",MESS_matrix::static_type_name().c_str());
        return octave_value();
    }

    ret = mess_matrix_init(&result);        OCTMESS_ERROR(ret,(ret!=0),mess_matrix_init);
    switch (umap)
    {
        // handle "special" cases first
        case umap_real:
            ret = mess_matrix_copy(this->ptr,result);       OCTMESS_ERROR(ret,(ret!=0),mess_matrix_copy);
            ret = mess_matrix_toreal(result);               OCTMESS_ERROR(ret,(ret!=0),mess_matrix_toreal);
            return new MESS_matrix(result);
        case umap_imag:
            ret = mess_matrix_imagpart(this->ptr,result);   OCTMESS_ERROR(ret,(ret!=0),mess_matrix_imagpart);
            return new MESS_matrix(result);
        case umap_conj:
            ret = mess_matrix_copy(this->ptr,result);       OCTMESS_ERROR(ret,(ret!=0),mess_matrix_copy);
            ret = mess_matrix_conj(result);                 OCTMESS_ERROR(ret,(ret!=0),mess_matrix_conj);
            return new MESS_matrix(result);
        case umap_finite:
            ret = mess_matrix_copy(this->ptr,result);       OCTMESS_ERROR(ret,(ret!=0),mess_matrix_copy);
            ret = mess_matrix_map_isfinite(result,1);       OCTMESS_ERROR(ret,(ret!=0),mess_matrix_map_isfinite);
            return new MESS_matrix(result);

            // use macro for the other cases
            OCTMESS_GENERATOR(abs);
            OCTMESS_GENERATOR(acos);
            OCTMESS_GENERATOR(acosh);
            OCTMESS_GENERATOR(arg);
            OCTMESS_GENERATOR(asin);
            OCTMESS_GENERATOR(asinh);
            OCTMESS_GENERATOR(atan);
            OCTMESS_GENERATOR(atanh);
            OCTMESS_GENERATOR(ceil);
            OCTMESS_GENERATOR(cos);
            OCTMESS_GENERATOR(cosh);
            OCTMESS_GENERATOR(exp);
            OCTMESS_GENERATOR(expm1);
            OCTMESS_GENERATOR(floor);
            OCTMESS_GENERATOR(log);
            OCTMESS_GENERATOR(round);
            OCTMESS_GENERATOR(sin);
            OCTMESS_GENERATOR(sinh);
            OCTMESS_GENERATOR(sqrt);
            OCTMESS_GENERATOR(tan);
            OCTMESS_GENERATOR(tanh);
            OCTMESS_GENERATOR(isinf);
            OCTMESS_GENERATOR(isnan);
        default:
            goto not_defined;
    }

not_defined:
    error ("%s: not defined for %s", get_umap_name (umap), MESS_matrix::static_type_name ().c_str ());
    mess_matrix_clear(&result);
    return octave_value();
}




/*-----------------------------------------------------------------------------
 *  conversion operator
 *-----------------------------------------------------------------------------*/
MESS_matrix::operator mess_matrix(void) const{ return this->ptr; }
MESS_matrix::operator bool(void) const{ return this->ptr != NULL; }

/*-----------------------------------------------------------------------------
 *  destructor
 *-----------------------------------------------------------------------------*/
MESS_matrix::~MESS_matrix(void){
    if(this->ptr){
        mess_matrix_clear(&(this->ptr));
    }
    this->ptr=NULL;
}







