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
#include "MESS_vector.h"
#include "MESS_error.h"


DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (MESS_vector, "mess_vector", "mess_vector");


/*-----------------------------------------------------------------------------
 *  constructors
 *-----------------------------------------------------------------------------*/
MESS_vector::MESS_vector(void):ptr(NULL) { }
MESS_vector::MESS_vector(mess_vector v):ptr(v) { }

MESS_vector::MESS_vector(Matrix m){
    mess_vector_init(&(this->ptr));
    mess_vector_from_farray(this->ptr, m.rows()*m.cols(), m.fortran_vec(), NULL);
}


MESS_vector::MESS_vector(ComplexMatrix m){
    mess_vector_init(&(this->ptr));
    mess_vector_from_farray(this->ptr, m.rows()*m.cols(), NULL, reinterpret_cast <mess_double_cpx_t*>(m.fortran_vec()));
}


MESS_vector::MESS_vector(double m){
    mess_vector_init(&(this->ptr));
    mess_vector_alloc(this->ptr,1,MESS_REAL);
    this->ptr->values[0]=m;
}


MESS_vector::MESS_vector(Complex m){
    double * ptr;
    mess_vector_init(&(this->ptr));
    mess_vector_alloc(this->ptr,1,MESS_COMPLEX);
    ptr = (double *) this->ptr->values_cpx;
    ptr[0] = m.real();
    ptr[1] = m.imag();
}

MESS_vector::MESS_vector(const MESS_vector & v): ptr(v.ptr){ }

/*-----------------------------------------------------------------------------
 *  getter
 *-----------------------------------------------------------------------------*/
mess_vector MESS_vector::get_ptr(void) const{ return this->ptr; }


/*-----------------------------------------------------------------------------
 *  deep copy
 *-----------------------------------------------------------------------------*/
MESS_vector MESS_vector::deepcopy(void) const{
    MESS_vector v;
    if(*this){
        mess_vector ptr = v.get_ptr();
        mess_vector_init(&ptr);
        mess_vector_copy(this->ptr,ptr);
    }
    return v;
}


/*-----------------------------------------------------------------------------
 *  print functions octave
 *-----------------------------------------------------------------------------*/
void MESS_vector::print(std::ostream & os, bool pr_as_read_syntax) const {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_vector::print(std::ostream & os, bool pr_as_read_syntax) {
    print_raw(os, pr_as_read_syntax);
    newline(os);
}


void MESS_vector::print_raw(std::ostream & os, bool pr_as_read_syntax) const {
    os << MESS_vector::static_type_name();
#ifdef MESS_DEBUG
    os << " ( " << this->ptr << " ) ";
#endif
}


bool MESS_vector::print_name_tag (std::ostream& os, const std::string& name) const {
    os << name << " = ";
    return true;
}




/*-----------------------------------------------------------------------------
 *  overrides from octave_base_value
 *-----------------------------------------------------------------------------*/
bool MESS_vector::is_constant(void) const { return true; }
bool MESS_vector::is_defined(void) const { return true; }
bool MESS_vector::print_as_scalar(void) const { return true; }
bool MESS_vector::is_matrix_type(void) const { return true; }
bool MESS_vector::is_numeric_type(void) const { return true; }
bool MESS_vector::is_real_matrix(void) const { return this->ptr?MESS_IS_REAL(this->ptr):false; }
bool MESS_vector::is_real_type(void) const { return this->ptr?MESS_IS_REAL(this->ptr):false; }
bool MESS_vector::is_double_type(void) const { return true; };
bool MESS_vector::is_single_type(void) const { return false; }
bool MESS_vector::is_complex_matrix (void) const { return this->ptr?MESS_IS_COMPLEX(this->ptr):false; }
bool MESS_vector::is_complex_type (void) const { return this->ptr?MESS_IS_COMPLEX(this->ptr):false; }
bool MESS_vector::is_empty(void) const { return this->ptr?false:true; }

size_t MESS_vector::byte_size(void) const{
    if(this->ptr && (MESS_IS_REAL(this->ptr) || MESS_IS_COMPLEX(this->ptr))){
        return mess_vector_memsize(this->ptr);
    }
    return 0;
}


dim_vector MESS_vector::dims (void) const { return this->ptr?dim_vector(this->ptr->dim,1):dim_vector(0,0); }
int MESS_vector::ndims (void) const{ return 1; }
size_t MESS_vector::length(void) const { return this->ptr?this->ptr->dim:0; }


/*-----------------------------------------------------------------------------
 *  overrides for built in functions
 *-----------------------------------------------------------------------------*/
#undef OCTMESS_GENERATOR
#define OCTMESS_GENERATOR(MAPNAME)                                                                          \
    case umap_##MAPNAME:                                                                                    \
        ret = mess_vector_copy(this->ptr, result);  OCTMESS_ERROR(ret,(ret!=0),mess_vector_copy);           \
        ret = mess_vector_map_##MAPNAME(result);    OCTMESS_ERROR(ret,(ret!=0),mess_vector_map_##MAPNAME);  \
        return new MESS_vector(result);                                                                     \

octave_value MESS_vector::map (unary_mapper_t umap) const{

    int ret = 0;
    mess_vector result = NULL;

    ret = mess_vector_init(&result);                        OCTMESS_ERROR(ret, (ret!=0), mess_vector_init);
    switch (umap)
    {
        // handle "special" cases first
        case umap_real:
            ret = mess_vector_realpart(this->ptr,result);   OCTMESS_ERROR(ret,(ret!=0),mess_vector_realpart);
            return new MESS_vector(result);
        case umap_imag:
            ret = mess_vector_imagpart(this->ptr,result);   OCTMESS_ERROR(ret,(ret!=0),mess_vector_imagpart);
            return new MESS_vector(result);
       case umap_finite:
            ret = mess_vector_copy(this->ptr,result);       OCTMESS_ERROR(ret,(ret!=0),mess_vector_copy);
            ret = mess_vector_map_isfinite(result);         OCTMESS_ERROR(ret,(ret!=0),mess_vector_map_isfinite);
            return new MESS_vector(result);

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
            OCTMESS_GENERATOR(conj);
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
    error ("%s: not defined for %s", get_umap_name (umap), MESS_vector::static_type_name ().c_str ());
    mess_vector_clear(&result);
    return octave_value();
}


/*-----------------------------------------------------------------------------
 *  conversion operator
 *-----------------------------------------------------------------------------*/
MESS_vector::operator mess_vector(void) const { return this->ptr; }
MESS_vector::operator bool(void) const { return this->ptr != NULL; }


/*-----------------------------------------------------------------------------
 *  destructor
 *-----------------------------------------------------------------------------*/
MESS_vector::~MESS_vector(void){
    if(this->ptr){
        mess_vector_clear(&(this->ptr));
    }
    this->ptr=NULL;
}
