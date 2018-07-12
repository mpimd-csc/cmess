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
 * @addtogroup interfaces_python
 * @file python/c_interface/conv/matrix.c
 * @brief Convert a @python matrix to a @ref mess_matrix or vice-versa.
 *
 * @author @koehlerm
 * @author @nitin
 *
 *
 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL MESS_VECTOR_MATRIX_PYTHON_C_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "../pymess.h"
#include "cscutils/error_message.h"

#if PY_MAJOR_VERSION < 3
#define ConvStringtoC(X) PyBytes_AsString(X)
#else
#define ConvStringtoC(X) PyBytes_AsString(PyUnicode_AsUTF8String(X))
#endif


#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_OWNDATA NPY_OWNDATA
#define PyArray_ENABLEFLAGS(arr,flag)  ((PyArrayObject*)(arr))->flags |= (flag)
#endif



#define PyArray_CALL_ISSIGNED(x) PyArray_ISSIGNED((PyArrayObject*)(x))
#define PyArray_CALL_ISUNSIGNED(x) PyArray_ISUNSIGNED((PyArrayObject*)(x))
#define PyArray_CALL_DATA(x) PyArray_DATA((PyArrayObject*)(x))


/**
 * @brief Converts @python matrix to @ref mess_matrix.
 * @param[in] data input @python matrix
 * @return converted matrix as @ref mess_matrix
 *
 * The @ref matrix_to_c function supports the numpy ndarray and scipy csc coo and csr matrix and the following datatypes
 *  \li int8
 *  \li int16
 *  \li int32
 *  \li int64
 *  \li float32
 *  \li float64
 *  \li complex64
 *  \li complex128
 *
 *  If the matrix cannot be converted into a MESS internal matrix @c NULL is returned. This either identifies that an
 *  error occured or it is an empty matrix or of unknown type. Please check the return value before your use it.
 *  If a conversion error occured the @python-Execption is set.
 *
 *
 */
mess_matrix matrix_to_c(PyObject *data){
    MSG_FNAME(__func__);
    int ret = 0;
    PyObject *module = NULL;

    if(PyObject_HasAttrString(data, "__module__")) {
        module  = (PyObject*) PyObject_GetAttrString(data, "__module__");
        MSG_INFO("Python modules name: %s\n",ConvStringtoC(module));


        /*-----------------------------------------------------------------------------
         *  CSR Matrix
         *-----------------------------------------------------------------------------*/
        if (strcmp(ConvStringtoC(module),"scipy.sparse.csr") == 0) {
#if MESS_DEBUG
            MSG_PRINT("scipy.sparse.csr\n");
#endif
            PyObject *temp = NULL, *temp2 = NULL, *dtype = NULL ;
            mess_int_t i;
            char *dtype_str = NULL;
            PyObject *value_data = NULL;
            npy_int64 *array = NULL;
            mess_matrix csr_matrix;

            ret = mess_matrix_init(&csr_matrix);
            if ( ret != 0 ) {
                PyErr_SetString(PyExc_RuntimeError,"Cannot initialize internal matrix.");
                Py_DECREF(module);
                return NULL;
            }

            PY_GET_LONG(csr_matrix->nnz, data, "nnz");
            temp = (PyObject*) PyObject_GetAttrString(data, "shape");
            temp2 = (PyObject*) PyTuple_GET_ITEM(temp, 1); // Only a borrowed reference
            csr_matrix->cols = (mess_int_t) PyLong_AsLong(temp2);
            temp2 = (PyObject*) PyTuple_GetItem(temp, 0);
            csr_matrix->rows = (mess_int_t) PyLong_AsLong(temp2);
            Py_DECREF(temp);

            csr_matrix->ld = csr_matrix->rows;
            csr_matrix->store_type = MESS_CSR;
            csr_matrix->symmetry = MESS_GENERAL;

            temp = (PyObject*) PyObject_GetAttrString(data, "dtype");
            dtype  = (PyObject*) PyObject_GetAttrString(temp, "name");
            dtype_str = ConvStringtoC(dtype);
            Py_DECREF(temp);

            value_data = (PyObject *) PyObject_GetAttrString(data,"data");
            if ( strcmp(dtype_str,"float64") == 0 || strcmp(dtype_str,"float32") == 0|| PyArray_CALL_ISSIGNED(value_data)|| PyArray_CALL_ISUNSIGNED(value_data)){
                PyObject *values = NULL;
                csr_matrix->values = malloc(sizeof(double)*csr_matrix->nnz);
                if ( strcmp(dtype_str,"float64") != 0 ) {
                    PyErr_Clear();
                    values = PyArray_FROM_OT(value_data, NPY_DOUBLE);
                } else {
                    values = value_data;
                }
                memcpy(csr_matrix->values,(double*) PyArray_CALL_DATA(values),sizeof(double)*csr_matrix->nnz);
                csr_matrix->values_cpx = NULL;
                csr_matrix->data_type = MESS_REAL;
                if ( values != value_data) { Py_DECREF(values); }
            }else if (strcmp(dtype_str,"complex128") == 0 || strcmp(dtype_str,"complex64") == 0)  {
                PyObject *values = NULL;
                csr_matrix->values_cpx = malloc(sizeof(mess_double_cpx_t)*csr_matrix->nnz);
                if ( strcmp(dtype_str,"complex64") == 0 ) {
                    PyErr_Clear();
                    values = PyArray_FROM_OT(value_data, NPY_COMPLEX128);
                } else {
                    values = value_data;
                }
                memcpy(csr_matrix->values_cpx,(mess_double_cpx_t *) PyArray_CALL_DATA(values),sizeof(mess_double_cpx_t)*csr_matrix->nnz);
                csr_matrix->values = NULL;
                csr_matrix->data_type = MESS_COMPLEX;
                if ( values != value_data) { Py_DECREF(values); }
            } else {
                mess_matrix_clear(&csr_matrix);
                Py_DECREF(dtype);
                Py_DECREF(value_data);
                Py_DECREF(module);
                PyErr_SetString(PyExc_TypeError,"Argument Matrices should have either integer, float32/64 or complex64/128 entries");
                return NULL;
            }
            Py_DECREF(dtype);
            Py_DECREF(value_data);

            // Get row pointer
            temp = (PyObject*) PyObject_GetAttrString(data, "indptr");
            PyErr_Clear();
            PyObject * indptr = PyArray_FROM_OT(temp, NPY_INT64);
            array = (npy_int64 *) PyArray_CALL_DATA(indptr);
            csr_matrix->rowptr = malloc(sizeof(mess_int_t)*(csr_matrix->rows+1));
            for (i = 0; i < csr_matrix->rows+1; i++) {
                csr_matrix->rowptr[i] = (mess_int_t) array[i];
            }
            Py_DECREF(temp);
            Py_DECREF(indptr);

            // Get col pointer
            temp = (PyObject*) PyObject_GetAttrString(data, "indices");
            PyErr_Clear();
            PyObject * indices = PyArray_FROM_OT(temp, NPY_INT64);
            array = (npy_int64 *) PyArray_CALL_DATA(indices);
            csr_matrix->colptr = malloc(sizeof(mess_int_t)*csr_matrix->nnz);
            for (i = 0; i < csr_matrix->nnz; i++) {
                csr_matrix->colptr[i]=(mess_int_t) array[i];
            }
            Py_DECREF(temp);
            Py_DECREF(indices);

            Py_DECREF(module);
            return csr_matrix;
        }


        /*-----------------------------------------------------------------------------
         *  CSC Matrix
         *-----------------------------------------------------------------------------*/

        if (strcmp(ConvStringtoC(module),"scipy.sparse.csc") == 0){
            PyObject *temp = NULL, *temp2 = NULL, *dtype = NULL ;
            mess_int_t i;
            char *dtype_str = NULL;
            PyObject *value_data = NULL;
            npy_int64 *array = NULL;
            mess_matrix csc_matrix;

            ret = mess_matrix_init(&csc_matrix);
            if ( ret != 0 ) {
                PyErr_SetString(PyExc_RuntimeError,"Cannot initialize internal matrix.");
                Py_DECREF(module);
                return NULL;
            }

            PY_GET_LONG(csc_matrix->nnz, data, "nnz");
            temp = (PyObject*) PyObject_GetAttrString(data, "shape");
            temp2 = (PyObject*) PyTuple_GetItem(temp, 1);
            csc_matrix->cols = (mess_int_t) PyLong_AsLong(temp2);
            temp2 = (PyObject*) PyTuple_GetItem(temp, 0);
            csc_matrix->rows = (mess_int_t) PyLong_AsLong(temp2);
            Py_DECREF(temp);

            csc_matrix->ld = csc_matrix->rows;
            csc_matrix->store_type = MESS_CSC;
            csc_matrix->symmetry = MESS_GENERAL;

            /*-----------------------------------------------------------------------------
             *  get values
             *-----------------------------------------------------------------------------*/
            temp = (PyObject*) PyObject_GetAttrString(data, "dtype");
            dtype = (PyObject*) PyObject_GetAttrString(temp, "name");
            dtype_str = ConvStringtoC(dtype);
            Py_DECREF(temp);

            value_data = (PyObject *) PyObject_GetAttrString(data,"data");
            if ( strcmp(dtype_str,"float64") == 0 || strcmp(dtype_str,"float32") == 0 || PyArray_CALL_ISSIGNED(value_data) || PyArray_CALL_ISUNSIGNED(value_data)){
                PyObject *values = NULL;
                csc_matrix->values = malloc(sizeof(double)*csc_matrix->nnz);
                if ( strcmp(dtype_str,"float64") != 0 ) {
                    PyErr_Clear();
                    values = PyArray_FROM_OT(value_data, NPY_DOUBLE);
                } else {
                    values = value_data;
                }
                memcpy(csc_matrix->values,(double*) PyArray_CALL_DATA(values),sizeof(double)*csc_matrix->nnz);
                csc_matrix->values_cpx = NULL;
                csc_matrix->data_type = MESS_REAL;
                if ( values != value_data) { Py_DECREF(values); }
            }else if (strcmp(dtype_str,"complex128") == 0 || strcmp(dtype_str,"complex64") == 0)  {
                PyObject *values = NULL;
                csc_matrix->values_cpx = malloc(sizeof(mess_double_cpx_t)*csc_matrix->nnz);
                if ( strcmp(dtype_str,"complex64") == 0 ) {
                    PyErr_Clear();
                    values = PyArray_FROM_OT(value_data, NPY_COMPLEX128);
                } else {
                    values = value_data;
                }
                memcpy(csc_matrix->values_cpx,(mess_double_cpx_t *) PyArray_CALL_DATA(values),sizeof(mess_double_cpx_t)*csc_matrix->nnz);
                csc_matrix->values = NULL;
                csc_matrix->data_type = MESS_COMPLEX;
                if ( values != value_data) { Py_DECREF(values); }
            } else {
                mess_matrix_clear(&csc_matrix);
                Py_DECREF(dtype);
                Py_DECREF(value_data);
                Py_DECREF(module);
                PyErr_SetString(PyExc_TypeError,"Argument Matrices should have either integer, float32/64 or complex64/128 entries");
                return NULL;
            }
            Py_DECREF(dtype);
            Py_DECREF(value_data);

            // get col pointer
            temp = (PyObject*) PyObject_GetAttrString(data, "indptr");
            PyErr_Clear();
            PyObject * indptr = PyArray_FROM_OT(temp, NPY_INT64);
            array= (npy_int64 *) PyArray_CALL_DATA(indptr);
            csc_matrix->colptr = malloc(sizeof(mess_int_t)*(csc_matrix->cols+1));
            for (i = 0; i < csc_matrix->cols+1; i++) {
                csc_matrix->colptr[i]=(mess_int_t)array[i];
            }
            Py_DECREF(temp);
            Py_DECREF(indptr);

            // get row ptr
            temp = (PyObject*) PyObject_GetAttrString(data, "indices");
            PyErr_Clear();
            PyObject *indices = PyArray_FROM_OT(temp, NPY_INT64);
            array= (npy_int64*) PyArray_CALL_DATA(indices);
            csc_matrix->rowptr = malloc(sizeof(mess_int_t)*csc_matrix->nnz);
            for (i = 0; i < csc_matrix->nnz; i++) {
                csc_matrix->rowptr[i]=(mess_int_t) array[i];
            }

            Py_DECREF(temp);
            Py_DECREF(indices);
            Py_DECREF(module);
            return csc_matrix;
        }

        /*-----------------------------------------------------------------------------
         *  Coordinate Matrix
         *-----------------------------------------------------------------------------*/
        if (strcmp(ConvStringtoC(module),"scipy.sparse.coo") == 0){
            PyObject *temp = NULL, *temp2 = NULL, *dtype = NULL ;
            mess_int_t i;
            char *dtype_str = NULL;
            PyObject *value_data = NULL;
            npy_int64 *array = NULL;
            mess_matrix coo_matrix;

            ret = mess_matrix_init(&coo_matrix);
            if ( ret != 0 ) {
                PyErr_SetString(PyExc_RuntimeError,"Cannot initialize internal matrix.");
                Py_DECREF(module);
                return NULL;
            }

            PY_GET_LONG(coo_matrix->nnz, data, "nnz");
            temp = (PyObject*) PyObject_GetAttrString(data, "shape");
            temp2 = (PyObject*) PyTuple_GetItem(temp, 1);
            coo_matrix->cols = (mess_int_t) PyLong_AsLong(temp2);
            temp2 = (PyObject*) PyTuple_GetItem(temp, 0);
            coo_matrix->rows = (mess_int_t) PyLong_AsLong(temp2);
            Py_DECREF(temp);

            coo_matrix->ld = coo_matrix->rows;
            coo_matrix->store_type = MESS_COORD;
            coo_matrix->symmetry = MESS_GENERAL;

            // get col pointer
            temp = (PyObject*) PyObject_GetAttrString(data, "col");

            // MSG_PRINT("temp: %lu ndim: %d dim: %d\n", temp, PyArray_NDIM(temp), PyArray_DIM(temp,0));
            // MSG_PRINT("temp-dtype: %s\n", ConvStringtoC(PyObject_GetAttrString(PyObject_GetAttrString(temp,"dtype"),"name")));
            // fflush(stdout);
            PyErr_Clear();
            PyObject *colptr = PyArray_FROM_OT(temp, NPY_INT64);
            // MSG_PRINT("colptr: %lu\n", colptr);
            // MSG_PRINT("colptr-dtype: %s\n", ConvStringtoC(PyObject_GetAttrString(PyObject_GetAttrString(colptr,"dtype"),"name")));
            array= (npy_int64 *) PyArray_CALL_DATA(colptr);

            coo_matrix->colptr = malloc(sizeof(mess_int_t)*(coo_matrix->nnz));
            for (i = 0; i < coo_matrix->nnz; i++) {
                coo_matrix->colptr[i]=(mess_int_t)array[i];
            }
            Py_DECREF(colptr);
            Py_DECREF(temp);
            // get row ptr
            temp = (PyObject*) PyObject_GetAttrString(data, "row");
            PyErr_Clear();
            PyObject *rowptr = PyArray_FROM_OT(temp, NPY_INT64);
            array= (npy_int64*) PyArray_CALL_DATA(rowptr);
            coo_matrix->rowptr = malloc(sizeof(mess_int_t)*coo_matrix->nnz);
            for (i = 0; i < coo_matrix->nnz; i++) {
                coo_matrix->rowptr[i]=(mess_int_t) array[i];
            }
            Py_DECREF(rowptr);
            Py_DECREF(temp);

            /*-----------------------------------------------------------------------------
             *  get values
             *-----------------------------------------------------------------------------*/
            temp = (PyObject*) PyObject_GetAttrString(data, "dtype");
            dtype = (PyObject*) PyObject_GetAttrString(temp, "name");
            dtype_str = ConvStringtoC(dtype);
            Py_DECREF(temp);

            value_data = (PyObject *) PyObject_GetAttrString(data,"data");
            if ( strcmp(dtype_str,"float64") == 0 || strcmp(dtype_str,"float32") == 0|| PyArray_CALL_ISSIGNED(value_data)|| PyArray_CALL_ISUNSIGNED(value_data)){
                PyObject *values = NULL;
                coo_matrix->values = malloc(sizeof(double)*coo_matrix->nnz);
                if ( strcmp(dtype_str,"float64") != 0 ) {
                    PyErr_Clear();
                    values = PyArray_FROM_OT(value_data, NPY_DOUBLE);
                } else {
                    values = value_data;
                }
                memcpy(coo_matrix->values,(double*) PyArray_CALL_DATA(values),sizeof(double)*coo_matrix->nnz);
                coo_matrix->values_cpx = NULL;
                coo_matrix->data_type = MESS_REAL;
                if ( values != value_data) { Py_DECREF(values); }
            }else if (strcmp(dtype_str,"complex128") == 0 || strcmp(dtype_str,"complex64") == 0)  {
                PyObject *values = NULL;
                coo_matrix->values_cpx = malloc(sizeof(mess_double_cpx_t)*coo_matrix->nnz);
                if ( strcmp(dtype_str,"complex64") == 0 ) {
                    PyErr_Clear();
                    values = PyArray_FROM_OT(value_data, NPY_COMPLEX128);
                } else {
                    values = value_data;
                }
                memcpy(coo_matrix->values_cpx,(mess_double_cpx_t *) PyArray_CALL_DATA(values),sizeof(mess_double_cpx_t)*coo_matrix->nnz);
                coo_matrix->values = NULL;
                coo_matrix->data_type = MESS_COMPLEX;
                if ( values != value_data) { Py_DECREF(values); }
            } else {
                mess_matrix_clear(&coo_matrix);
                Py_DECREF(dtype);
                Py_DECREF(value_data);
                Py_DECREF(module);
                PyErr_SetString(PyExc_TypeError,"Argument Matrices should have either integer, float32/64 or complex64/128 entries");
                return NULL;
            }
            Py_DECREF(dtype);
            Py_DECREF(value_data);

            Py_DECREF(module);
            return coo_matrix;
        }


        Py_DECREF(module);

    }

    /*-----------------------------------------------------------------------------
     *  Dense Matrix
     *-----------------------------------------------------------------------------*/
    if (PyArray_Check(data))  {
        PyObject *temp=NULL, *dtype = NULL;
        char *dtype_str = NULL;
        mess_int_t i = 0,j;
        size_t type_size = 1;
        npy_intp stride_rows, stride_cols, rows, cols;
        mess_datatype_t data_type = MESS_REAL;
        int decref_data = 0;
        mess_int_t ndim = 0;

        PY_GET_LONG(ndim, data, "ndim");
        if (ndim < 1 || ndim>2){
            MSG_ERROR("Argument type should be a dense matrix, 2d-array or a scipy.sparse matrix (coo, csr or csc)\n");
            PyErr_SetString(PyExc_TypeError,
                    "Argument type should be a dense matrix, 2d-array or a scipy.sparse matrix (coo, csr or csc)");
            return NULL;
        }
        if ( PyArray_CALL_ISUNSIGNED(data) || PyArray_CALL_ISSIGNED(data)) {
            PyErr_Clear();
            data = PyArray_FROM_OT(data, NPY_DOUBLE);
            decref_data = 1;
        }
        mess_matrix dense_matrix;
        mess_matrix_init(&dense_matrix);
        dense_matrix->rowptr = NULL;
        dense_matrix->colptr = NULL;

        // Get size
        if ( ndim == 1) {
            rows = (mess_int_t) PyArray_DIM((PyArrayObject*) data,0);
            cols = 1;
            stride_rows = PyArray_STRIDE((PyArrayObject*) data,0);
            stride_cols = 1;
        } else {
            rows = (mess_int_t) PyArray_DIM((PyArrayObject*) data,0);
            cols = (mess_int_t) PyArray_DIM((PyArrayObject*) data,1);
            stride_rows = PyArray_STRIDE((PyArrayObject*) data,0);
            stride_cols = PyArray_STRIDE((PyArrayObject*) data,1);
        }
        /* MSG_PRINT("rows( dim 0 )  = %d \t stride = %d\n", rows, stride_rows);
           MSG_PRINT("cols( dim 1 )  = %d \t stride = %d\n", cols, stride_cols);  */

        temp = (PyObject*) PyObject_GetAttrString(data, "dtype");
        dtype = (PyObject*) PyObject_GetAttrString(temp, "name");
        dtype_str = strdup(ConvStringtoC(dtype));
        Py_DECREF(temp);
        Py_DECREF(dtype);
        // MSG_PRINT("dtype_str: %s\n", dtype_str);
        if ( strcmp(dtype_str, "float64") == 0  ) {
            type_size = sizeof(double);
            data_type = MESS_REAL;
        } else if ( strcmp(dtype_str, "float32") == 0 ) {
            type_size = sizeof(float);
            data_type = MESS_REAL;
        } else if ( strcmp(dtype_str, "complex128") == 0 ) {
            type_size = sizeof(mess_double_cpx_t);
            data_type = MESS_COMPLEX;
        } else if ( strcmp(dtype_str, "complex64") == 0 ){
            type_size = sizeof(mess_float_cpx_t);
            data_type = MESS_COMPLEX;
        } else {
            mess_matrix_clear(&dense_matrix);
            free(dtype_str);
            dtype_str = NULL;
            PyErr_SetString(PyExc_TypeError, "Argument Matrices should have either float32/64 or complex64/128 entries");
            if (decref_data) {Py_DECREF(data); }
            return NULL;
        }
        stride_rows/= type_size;
        stride_cols/= type_size;
        /* MSG_PRINT("rows( dim 0 )  = %d \t stride = %d\n", rows, stride_rows);
           MSG_PRINT("cols( dim 1 )  = %d \t stride = %d\n", cols, stride_cols);  */

        ret = mess_matrix_alloc(dense_matrix, rows, cols, rows*cols,MESS_DENSE, data_type);
        if ( ret != 0 ) {
            free(dtype_str);
            dtype_str = NULL;
            PyErr_SetString(PyExc_RuntimeError,"Cannot allocate matrix in C");
            if (decref_data) {Py_DECREF(data); }
            return NULL;
        }

        if ( strcmp(dtype_str, "float64") == 0  ) {
            double *values = PyArray_CALL_DATA(data);
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    dense_matrix->values[i+j*dense_matrix->ld] = values[i*stride_rows+j*stride_cols];
                }
            }
        } else if ( strcmp(dtype_str, "float32") == 0 ) {
            float *values = PyArray_CALL_DATA(data);
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    dense_matrix->values[i+j*dense_matrix->ld] = (double) values[i*stride_rows+j*stride_cols];
                }
            }
        } else if ( strcmp(dtype_str, "complex128") == 0 ) {
            mess_double_cpx_t *values = PyArray_CALL_DATA(data);
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    dense_matrix->values_cpx[i+j*dense_matrix->ld] = values[i*stride_rows+j*stride_cols];
                }
            }
        } else if ( strcmp(dtype_str, "complex64") == 0 ){
            mess_float_cpx_t *values = PyArray_CALL_DATA(data);
            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    dense_matrix->values_cpx[i+j*dense_matrix->ld] = (mess_double_cpx_t) values[i*stride_rows+j*stride_cols];
                }
            }
        }
        free(dtype_str); dtype_str = NULL;
        if (decref_data) {Py_DECREF(data); }
        return dense_matrix;
    }
    /* This error needs to be commented out since one can use the NULL return value for identifying that a matrix does not exists or is
     * empty. */
    /* PyErr_SetString(PyExc_TypeError,"Argument type should be a dense matrix, 2d-array or a scipy.sparse matrix (csr/csc/coo)");
       MSG_ERROR("HIER\n");  */
    return NULL;
}

/**
 * @brief Converts mess_matrix to @python matrix.
 * @param[in] c_matrix  input @ref mess_matrix
 * @return @python matrix
 * Attention the matrix is moved from C to @python. The matrix does not longer
 * exist correctly in C.
 */
PyObject* matrix_to_python(mess_matrix c_matrix) {
    npy_intp i = 0, foo = 0 ;
    npy_int * array = NULL, *array2 = NULL ;

    /*-----------------------------------------------------------------------------
     *  CSR Matrix
     *-----------------------------------------------------------------------------*/
    if (c_matrix->store_type == MESS_CSR) {
        PyObject* ReturnDataArray    = NULL;
        PyObject* ReturnIndicesArray = NULL;
        PyObject* ReturnIndptrArray  = NULL;

        if(c_matrix->data_type == MESS_REAL){
            npy_intp nnz = (npy_intp) c_matrix->nnz;
            ReturnDataArray = PyArray_SimpleNewFromData(1, &nnz, NPY_FLOAT64, (void *) c_matrix->values);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnDataArray, NPY_ARRAY_OWNDATA);
        }
        else if(c_matrix->data_type == MESS_COMPLEX) {
            npy_intp nnz = (npy_intp) c_matrix->nnz;
            ReturnDataArray = PyArray_SimpleNewFromData(1, &nnz, NPY_COMPLEX128, (void *) c_matrix->values_cpx);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnDataArray, NPY_ARRAY_OWNDATA);
        }
        else {
            PyErr_SetString(PyExc_TypeError, "C-MESS Matrices can only be Real or Complex");
            return NULL;
        }

        // rowptr
        foo = (npy_intp) c_matrix->rows +1 ;
        ReturnIndptrArray = PyArray_SimpleNew(1, &foo, NPY_INT);
        array = (npy_int *) PyArray_CALL_DATA(ReturnIndptrArray);
        for ( i = 0; i < foo; i++) {
            array[i] = (npy_int ) c_matrix->rowptr[i];
        }

        foo = (npy_intp) c_matrix->nnz;
        ReturnIndicesArray = PyArray_SimpleNew(1, &foo, NPY_INT);
        array = (npy_int *) PyArray_CALL_DATA(ReturnIndicesArray);
        for (i = 0; i < foo; i++) {
            array[i] = (npy_int) c_matrix->colptr[i];
        }
        //Import Scipy.sparse
        PyObject *ScipySparse = PyImport_ImportModule((char*)"scipy.sparse");
        PyObject *Create_CSRmatrix = PyObject_GetAttrString(ScipySparse,(char*)"csr_matrix");
        PyObject* ArrayReturn = Py_BuildValue("(OOO)", ReturnDataArray, ReturnIndicesArray, ReturnIndptrArray);

        //Create tuple of cols and rows
        PyObject* Dim_Return = Py_BuildValue("(ii)", c_matrix->rows, c_matrix->cols);
        //Create a python csr_matrix
        PyObject *t =   Py_BuildValue("(OO)", ArrayReturn, Dim_Return);
        PyObject *Return_CSRMatrix =  PyObject_CallObject(Create_CSRmatrix, t);
        Py_DECREF(t);

        //clean up
        Py_DECREF(ScipySparse);
        Py_DECREF(Create_CSRmatrix);
        Py_DECREF(ReturnDataArray);
        Py_DECREF(ReturnIndicesArray);
        Py_DECREF(ReturnIndptrArray);
        Py_DECREF(ArrayReturn);
        Py_DECREF(Dim_Return);
        c_matrix->values = NULL;
        c_matrix->values_cpx = NULL;
        MESS_MATRIX_RESET(c_matrix);
        return Return_CSRMatrix;
    }


    /*-----------------------------------------------------------------------------
     *  CSC Matrix
     *-----------------------------------------------------------------------------*/
    else if (c_matrix->store_type == MESS_CSC){
        PyObject* ReturnDataArray    = NULL;
        PyObject* ReturnIndicesArray = NULL;
        PyObject* ReturnIndptrArray  = NULL;

        if(c_matrix->data_type == MESS_REAL){
            npy_intp nnz = (npy_intp) c_matrix->nnz;
            ReturnDataArray = PyArray_SimpleNewFromData(1, &nnz, NPY_FLOAT64, c_matrix->values);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnDataArray, NPY_ARRAY_OWNDATA);
        }
        else if(c_matrix->data_type == MESS_COMPLEX) {
            npy_intp nnz = (npy_intp) c_matrix->nnz;
            ReturnDataArray = PyArray_SimpleNewFromData(1, &nnz, NPY_COMPLEX128, c_matrix->values_cpx);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnDataArray, NPY_ARRAY_OWNDATA);
        }
        else {
            PyErr_SetString(PyExc_TypeError, "C-MESS Matrices can only be Real or Complex");
            return NULL;
        }

        foo = (npy_intp) c_matrix->cols+1;
        ReturnIndptrArray = PyArray_SimpleNew(1, &foo, NPY_INT);
        array = (npy_int *)  PyArray_CALL_DATA(ReturnIndptrArray);
        for (i = 0; i < foo; i++) {
            array[i]=(npy_int) c_matrix->colptr[i];
        }

        foo = (npy_intp) c_matrix->nnz;
        ReturnIndicesArray = PyArray_SimpleNew(1, &foo, NPY_INT);
        array = (npy_int *) PyArray_CALL_DATA(ReturnIndicesArray);
        for (i = 0; i < foo; i++) {
            array[i]=(npy_int)c_matrix->rowptr[i];
        }

        //Import Scipy.sparse
        PyObject *ScipySparse = PyImport_ImportModule((char*)"scipy.sparse");
        PyObject *Create_CSCmatrix = PyObject_GetAttrString(ScipySparse,(char*)"csc_matrix");

        PyObject* ArrayReturn = Py_BuildValue("(OOO)", ReturnDataArray, ReturnIndicesArray, ReturnIndptrArray);

        //Create tuple of cols and rows
        PyObject* Dim_Return = Py_BuildValue("(ii)", c_matrix->rows, c_matrix->cols);
        //Create a python csr_matrix

        PyObject *t = Py_BuildValue("(OO)", ArrayReturn, Dim_Return);
        PyObject *Return_CSCMatrix =  PyObject_CallObject(Create_CSCmatrix, t);
        Py_DECREF(t);

        //clean up
        Py_DECREF(ScipySparse);
        Py_DECREF(Create_CSCmatrix);
        Py_DECREF(ReturnDataArray);
        Py_DECREF(ReturnIndicesArray);
        Py_DECREF(ReturnIndptrArray);
        Py_DECREF(ArrayReturn);
        Py_DECREF(Dim_Return);
        c_matrix->values = NULL;
        c_matrix->values_cpx = NULL;
        MESS_MATRIX_RESET(c_matrix);

        return Return_CSCMatrix;
    }

    /*-----------------------------------------------------------------------------
     *  Coordinate Matrix
     *-----------------------------------------------------------------------------*/
    else if (c_matrix->store_type == MESS_COORD){
        PyObject* ReturnDataArray    = NULL;
        PyObject* ReturnRowArray = NULL;
        PyObject* ReturnColArray  = NULL;
        npy_intp nnz = (npy_intp) c_matrix->nnz;

        if(c_matrix->data_type == MESS_REAL){
            ReturnDataArray = PyArray_SimpleNewFromData(1, &nnz, NPY_FLOAT64, c_matrix->values);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnDataArray, NPY_ARRAY_OWNDATA);
        }
        else if(c_matrix->data_type == MESS_COMPLEX) {
            ReturnDataArray = PyArray_SimpleNewFromData(1, &nnz, NPY_COMPLEX128, c_matrix->values_cpx);
            PyArray_ENABLEFLAGS((PyArrayObject *)ReturnDataArray, NPY_ARRAY_OWNDATA);
        }
        else {
            PyErr_SetString(PyExc_TypeError, "C-MESS Matrices can only be Real or Complex");
            return NULL;
        }

        ReturnColArray = PyArray_SimpleNew(1, &nnz, NPY_INT);
        ReturnRowArray = PyArray_SimpleNew(1, &nnz, NPY_INT);
        array = (npy_int *)  PyArray_CALL_DATA(ReturnColArray);
        array2 = (npy_int *)  PyArray_CALL_DATA(ReturnRowArray);
        for (i = 0; i < nnz; i++) {
            array[i]=(npy_int) c_matrix->colptr[i];
            array2[i]=(npy_int) c_matrix->rowptr[i];
        }
        //Import Scipy.sparse
        PyObject *ScipySparse = PyImport_ImportModule((char*)"scipy.sparse");
        PyObject *Create_CSCmatrix = PyObject_GetAttrString(ScipySparse,(char*)"coo_matrix");

        PyObject* ArrayReturn = Py_BuildValue("(O(OO))", ReturnDataArray, ReturnRowArray, ReturnColArray);

        //Create tuple of cols and rows
        PyObject* Dim_Return = Py_BuildValue("(ii)", c_matrix->rows, c_matrix->cols);
        //Create a python csr_matrix
        PyObject *t = Py_BuildValue("(OO)", ArrayReturn, Dim_Return);
        PyObject *Return_COOMatrix =  PyObject_CallObject(Create_CSCmatrix, t);
        Py_DECREF(t);

        //clean up
        Py_DECREF(ScipySparse);
        Py_DECREF(Create_CSCmatrix);
        Py_DECREF(ReturnDataArray);
        Py_DECREF(ReturnRowArray);
        Py_DECREF(ReturnColArray);
        Py_DECREF(ArrayReturn);
        Py_DECREF(Dim_Return);
        c_matrix->values = NULL;
        c_matrix->values_cpx = NULL;
        MESS_MATRIX_RESET(c_matrix);

        return Return_COOMatrix;
    }

    /*-----------------------------------------------------------------------------
     *  DENSE Matrix
     *-----------------------------------------------------------------------------*/
    else if (c_matrix->store_type == MESS_DENSE) {
        PyObject* ReturnArray = NULL;
        npy_intp  dim[2] = {c_matrix->rows, c_matrix->cols};
        npy_intp stride_rows, stride_cols, rows, cols;
        mess_int_t j;


        if ( MESS_IS_REAL(c_matrix)) {
            ReturnArray = PyArray_SimpleNew(2, dim , NPY_FLOAT64);
        } else if (MESS_IS_COMPLEX(c_matrix)) {
            ReturnArray = PyArray_SimpleNew(2, dim , NPY_COMPLEX128);
        } else {
            PyErr_SetString(PyExc_TypeError, "C-MESS Matrices can only be Real or Complex");
            return NULL;
        }
        if ( ReturnArray == NULL ) {
            PyErr_SetString(PyExc_RuntimeError, "Cannot allocate Python matrix");
            return NULL;
        }
        // MSG_PRINT("\t ReturnArray  %d  ref: %d\n",__LINE__, ReturnArray->ob_refcnt);
        rows = (mess_int_t) PyArray_DIM((PyArrayObject*) ReturnArray,0);
        cols = (mess_int_t) PyArray_DIM((PyArrayObject*) ReturnArray,1);
        stride_rows = PyArray_STRIDE((PyArrayObject*) ReturnArray,0);
        stride_cols = PyArray_STRIDE((PyArrayObject*) ReturnArray,1);
        // MSG_PRINT("\t ReturnArray  %d  ref: %d\n",__LINE__, ReturnArray->ob_refcnt);

        if(c_matrix->data_type == MESS_REAL){
            double *values = PyArray_CALL_DATA(ReturnArray);
            stride_cols /= sizeof(double);
            stride_rows /= sizeof(double);

            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    values[i*stride_rows+j*stride_cols] = c_matrix->values[i+j*c_matrix->ld];
                }
            }

        }
        else if(c_matrix->data_type == MESS_COMPLEX){
            mess_double_cpx_t *values = PyArray_CALL_DATA(ReturnArray);
            stride_cols /= sizeof(mess_double_cpx_t );
            stride_rows /= sizeof(mess_double_cpx_t );

            for (i = 0; i < rows; i++) {
                for (j = 0; j < cols; j++) {
                    values[i*stride_rows+j*stride_cols] = c_matrix->values_cpx[i+j*c_matrix->ld];
                }
            }
        }
        // MSG_PRINT("\t ReturnArray  %d  ref: %d\n",__LINE__, ReturnArray->ob_refcnt);

        /*-----------------------------------------------------------------------------
         *  Call Numpy matrix constructor
         *-----------------------------------------------------------------------------*/
        PyObject *Numpy = PyImport_ImportModule((char*)"numpy");
        PyObject *Create_dense_matrix = PyObject_GetAttrString(Numpy,(char*)"matrix");
        PyObject *t = Py_BuildValue("(O)", ReturnArray);
        PyObject *ReturnMatrix = PyObject_CallObject(Create_dense_matrix, t );
        Py_DECREF(t);

        // MSG_PRINT("\t ReturnArray  %d  ref: %d\n",__LINE__, ReturnArray->ob_refcnt);
        // MSG_PRINT("\t t            %d  ref: %d\n",__LINE__, t->ob_refcnt);
        //clean up
        Py_DECREF(Numpy);
        Py_DECREF(Create_dense_matrix);
        Py_DECREF(ReturnArray);
        // MESS_MATRIX_RESET(c_matrix);
        // MSG_PRINT("\t ReturnArray  %d  ref: %d\n",__LINE__, ReturnArray->ob_refcnt);
        return ReturnMatrix;
    }
    else {
        PyErr_SetString(PyExc_TypeError,"C-MESS Matrices can have only CSR, CSC, COORD, or DENSE structures");
        return NULL;
    }
}



