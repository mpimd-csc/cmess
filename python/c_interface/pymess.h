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
 * @file python/c_interface/pymess.h
 * @brief Interface to @python.
 *
 */


#ifndef PYCMESS_H

#define PYCMESS_H

#undef MESS_DEBUG

#ifdef  __cplusplus
extern "C" {
#endif
    /**
     * @addtogroup interfaces_python
     * @{
     */
    int pymess_is_reimplemented(PyObject * base, PyObject *instance, const char *attr );

    //conversions for equations
    mess_equation eqn_conv_lyap(PyObject *obj, mess_freelist mem);
    mess_equation eqn_conv_riccati(PyObject *obj, mess_freelist mem);

    mess_equation eqn_conv_lyap_so1(PyObject *obj, mess_freelist mem);
    mess_equation eqn_conv_riccati_so1(PyObject *obj, mess_freelist mem);

    mess_equation eqn_conv_lyap_so2(PyObject *obj, mess_freelist mem);
    mess_equation eqn_conv_riccati_so2(PyObject *obj, mess_freelist mem);

    mess_equation eqn_conv_lyap_dae1(PyObject *obj, mess_freelist mem);
    mess_equation eqn_conv_riccati_dae1(PyObject *obj, mess_freelist mem);

    mess_equation eqn_conv_lyap_dae2(PyObject *obj, mess_freelist mem);
    mess_equation eqn_conv_riccati_dae2(PyObject *obj, mess_freelist mem);

    //test functions
    PyObject* pymess_test_matrix(PyObject *self, PyObject* args);
    PyObject* pymess_test_vector(PyObject *self, PyObject* args);
    PyObject* pymess_test_status(PyObject *self, PyObject *args);
    PyObject* pymess_test_equation(PyObject *self, PyObject *args);

    //lradi
    PyObject* pymess_lrnm(PyObject *self, PyObject *args);

    //lrnm
    PyObject* pymess_lradi(PyObject *self, PyObject *args);

    //dense_nm_gmpare
    PyObject* pymess_dense_nm_gmpare(PyObject * self, PyObject * args);

    //lyap
    PyObject* pymess_lyap(PyObject *self, PyObject *args);

    //sylvester_sparsedense
    PyObject* pymess_sylvester_sparsedense(PyObject *self, PyObject *args);

    //care
    PyObject* pymess_care(PyObject *self, PyObject *args);

    //glyap3 stuff: glyap and gstein
    PyObject* pymess_glyap(PyObject *self, PyObject *args);
    PyObject* pymess_gstein(PyObject *self, PyObject *args);

    //eps
    PyObject* pymess_eps(PyObject *self, PyObject* args);

    //mess_version
    PyObject* pymess_mess_version(PyObject *self, PyObject* args);
    PyObject* pymess_mess_version_verbose(PyObject *self, PyObject* args);
    PyObject* pymess_mess_version_major(PyObject *self, PyObject* args);
    PyObject* pymess_mess_version_minor(PyObject *self, PyObject* args);
    PyObject* pymess_mess_version_patch(PyObject *self, PyObject* args);
    PyObject* pymess_mess_is_debug(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_zlib(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_bzip2(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_umfpack(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_amd(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_colamd(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_cholmod(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_csparse(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_superlu(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_mklpardiso(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_arpack(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_matio(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_openmp(PyObject *self, PyObject* args);
    PyObject* pymess_mess_have_mess64(PyObject *self, PyObject* args);

    //direct select
    PyObject* pymess_direct_select(PyObject *self, PyObject* args);

    //direct chol select
    PyObject* pymess_direct_chol_select(PyObject *self, PyObject* args);

    //multidirect select
    PyObject* pymess_multidirect_select(PyObject*self, PyObject *args);

    //conversions matrices
    mess_matrix matrix_to_c(PyObject *data);
    PyObject* matrix_to_python(mess_matrix c_matrix);

    //conversions vectors
    mess_vector vector_to_c(PyObject* data);
    PyObject* vector_to_python(mess_vector c_vector, int copy);

    //conversions options
    mess_options mess_options_from_python(PyObject *options);
    PyObject *mess_status_to_python(mess_status stat);

    //interface to mess_set_error_level
    PyObject* pymess_set_errorlevel(PyObject*self, PyObject *args);

    mess_equation mess_equation_from_python ( PyObject * obj, mess_freelist  mem, mess_equation_t eqntype, int *need_callback );



#define PYCMESS_ERROR(ret, message) { if ( (ret) != 0 ) { PyErr_SetString(PyExc_RuntimeError,message); return NULL; } }

#define  MESS_Py_XDECREF(...) { \
    void * _stopper = (void * ) -1; \
    PyObject*  _M[] = { __VA_ARGS__, _stopper}; \
    int _i; \
    for ( _i = 0; _M[_i] != _stopper ; _i++) {\
        Py_XDECREF(_M[_i]);\
    }\
}

#define PY_SET_LONG(obj, name, value) do { \
    PyObject * _value = PyLong_FromLong((value)); \
    PyObject_SetAttrString((obj), (name), _value); \
    Py_XDECREF(_value); \
} while (0);

#define PY_SET_DOUBLE(obj, name, value) do { \
    PyObject * _value = PyFloat_FromDouble((value)); \
    PyObject_SetAttrString((obj), (name), _value); \
    Py_XDECREF(_value); \
} while (0);

#define PY_SET_VECTOR(obj, name, value) do { \
    PyObject * _value = vector_to_python(value,1); \
    PyObject_SetAttrString((obj), (name), _value); \
    Py_XDECREF(_value); \
} while (0);

#define PY_SET_BOOL(obj, name, value) do { \
    if ( (value) ){ \
        PyObject_SetAttrString((obj), (name), Py_True);\
    } else { \
        PyObject_SetAttrString((obj), (name), Py_False);\
    }\
} while (0);



#ifdef MESS_DEBUG
#define PY_GET_LONG(dest, from, name){                                                  \
    PyObject *pv = PyObject_GetAttrString(from, name);                              \
    if(pv){                                                                         \
        dest = PyLong_AsLong(pv);                                                   \
        MSG_PRINT("call PY_GET_LONG:pv!=NULL\t dest=%d\t name=%s \n",dest, name);       \
        Py_DECREF(pv);                                                              \
    }else{                                                                          \
        MSG_PRINT("call PY_GET_LONG:pv==NULL\t name=%s \n", name );                 \
    }                                                                               \
}
#else
#define PY_GET_LONG(dest, from, name){                                                  \
    PyObject *pv = PyObject_GetAttrString(from, name);                              \
    if (pv){                                                                        \
        dest = PyLong_AsLong(pv);                                                   \
        Py_DECREF(pv);                                                              \
    }                                                                               \
}
#endif



#ifdef MESS_DEBUG
#define PY_GET_DOUBLE(dest, from, name ) {                                                  \
    PyObject *pv = PyObject_GetAttrString(from, name);                                  \
    if(pv){                                                                             \
        dest = PyFloat_AsDouble(pv);                                                    \
        MSG_PRINT("call PY_GET_DOUBLE:pv!=NULL\t dest=%.3e\t name=%s \n",dest, name );  \
        Py_DECREF(pv);                                                                  \
    }else{                                                                              \
        MSG_PRINT("call PY_GET_DOUBLE:pv==NULL\t  name=%s \n", name );                      \
    }                                                                                   \
}
#else
#define PY_GET_DOUBLE(dest, from, name ){                                                   \
    PyObject *pv = PyObject_GetAttrString(from, name);                                  \
    if (pv){                                                                            \
        dest = PyFloat_AsDouble(pv);                                                    \
        Py_DECREF(pv);                                                                  \
    }                                                                                   \
}
#endif


#ifdef MESS_DEBUG
#define PY_GET_VECTOR(dest,from,name) {                                                     \
    PyObject *pv = PyObject_GetAttrString(from,name);                                   \
    if (pv !=NULL && pv != Py_None) {                                                   \
        dest = vector_to_c(pv);                                                         \
        MSG_PRINT("call PY_GET_VECTOR:pv!=NULL\t name=%s \n",name);                     \
        mess_vector_printinfo(dest);                                                    \
        mess_vector_printshort(dest);                                                   \
        Py_DECREF(pv);                                                                  \
    }else{                                                                              \
        Py_XDECREF(pv);                                                                 \
        dest = NULL;                                                                    \
    }                                                                                   \
}
#else
#define PY_GET_VECTOR(dest,from,name) {                                                     \
    PyObject *pv = PyObject_GetAttrString(from,name);                                   \
    if (pv !=NULL && pv != Py_None) {                                                   \
        dest = vector_to_c(pv);                                                         \
        Py_DECREF(pv);                                                                  \
    }else{                                                                              \
        Py_XDECREF(pv);                                                                 \
        dest = NULL;                                                                    \
    }                                                                                   \
}
#endif



#ifdef MESS_DEBUG
#define PY_GET_MATRIX(dest, from, name) {                                                   \
    PyObject *pv = PyObject_GetAttrString(from, name);                                  \
    if (pv && pv!=Py_None){                                                             \
        dest = matrix_to_c(pv);                                                         \
        if(dest) mess_matrix_sort(dest);                                                \
        MSG_PRINT("call PY_GET_MATRIX:pv!=NULL\t name=%s \n", name );                       \
        mess_matrix_printinfo(dest);                                                    \
        Py_DECREF(pv);                                                                  \
    }else{                                                                              \
        dest = NULL;                                                                    \
    }                                                                                   \
}
#else
#define PY_GET_MATRIX(dest, from, name) {                                                   \
    PyObject *pv = PyObject_GetAttrString(from, name);                                  \
    if (pv){                                                                            \
        dest = matrix_to_c(pv);                                                         \
        if(dest) mess_matrix_sort(dest);                                                \
        Py_DECREF(pv);                                                                  \
    }else{                                                                              \
        dest = NULL;                                                                    \
    }                                                                                   \
}
#endif


#define PYMESS_READ_MATRIX(mem,cmatrix,pymatrix,name,errormark)                                                 \
        (cmatrix) = matrix_to_c((pymatrix));                                                                        \
    if(!(cmatrix)){                                                                                             \
        PyErr_SetString(PyExc_RuntimeError,"Cannot transfer matrix "#name" from Python to C-M.E.S.S.");         \
        goto errormark;                                                                                         \
    }                                                                                                           \
    mess_freelist_add_mess_matrix(mem,(cmatrix));




    /** @}  */
#ifdef __cplusplus
    }
#endif

#endif /* end of include guard: PYCMESS_H */
