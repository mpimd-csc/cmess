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
 * @file lib/io/vector_read.c
 * @brief Read a vector from a matrix market file.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/error_macro.h"
#include "mess/mess.h"
#include "cscutils/io.h"


#ifdef MESS_HAVE_MATIO
#include <matio.h>
#define _MATIO_VERSION ( (MATIO_MAJOR_VERSION) *10000 + (MATIO_MINOR_VERSION) * 100 + (MATIO_RELEASE_LEVEL))
#define _MATIO_15  10500
#endif

/* Sets the last entry of an array to value */
#define     SET_LAST(array, len, value)     (array)[(len)-1] = (value)

/* Definitions regarding the file format */
#define LINE_LENGTH 1025
#define TOKEN_LENGTH 30
#define MM_BANNER       "%%MatrixMarket"
#define MM_STR_VECTOR       "vector"
#define MM_STR_REAL     "real"
#define MM_STR_INT      "integer"
#define MM_STR_CPX      "complex"

/**
 * @brief Convert a string to lower case.
 * @param[in,out] str  String to be converted to lower case
 *
 * Convert a string to lowercase.
 */
static void _lowercase (char *str)
{
    char *p;
    if ( str == NULL) return;
    for (p=str; *p!='\0'; *p=tolower(*p),p++);
}


/**
 * @brief Reads the header of a Matrix Market vector file.  (internal version)
 * @param[in] file  input filepointer
 * @param[out] vect  mess_vector to store the read information
 * @return zero on success or otherwise a non zero error code.
 *
 * The @ref __mm_read_vector_info function reads the header of a matrix market vector
 * file and writes the  basic information to matrix_vect.
 */
static int __mm_read_vector_info ( csc_io_file_t * file, mess_vector vect )
{
    MSG_FNAME(__func__);
    char line[LINE_LENGTH];
    char mmbanner[TOKEN_LENGTH];
    char vbanner[TOKEN_LENGTH];
    char datatype[TOKEN_LENGTH];
    char *cret = NULL;
    int ret = 0;
    mess_int_t dim;

    mess_check_nullpointer(file);
    mess_check_nullpointer(vect);

    if ( csc_io_scanf(file, "%s %s %s", mmbanner, vbanner, datatype) != 3 ) {
        MSG_ERROR("Can not read the complete header. \n");
        return (MESS_ERROR_WRONG_HEADER);
    }

    SET_LAST(mmbanner, TOKEN_LENGTH, '\0');
    SET_LAST(vbanner, TOKEN_LENGTH, '\0');
    SET_LAST(datatype, TOKEN_LENGTH, '\0');
    _lowercase(vbanner);
    _lowercase(datatype);

    if ( strncmp(mmbanner, MM_BANNER, strlen(MM_BANNER)) !=0 ) {
        MSG_ERROR("Wrong header information: %s\n",mmbanner);
        return (MESS_ERROR_WRONG_HEADER);
    }
    if (strncmp(vbanner, MM_STR_VECTOR, strlen (MM_STR_VECTOR)) !=0) {
        MSG_ERROR("Wrong header information: %s\n",vbanner);
        return (MESS_ERROR_WRONG_HEADER);
    }

    if (strncmp(datatype, MM_STR_REAL, strlen(MM_STR_REAL)) == 0 ) {
        ret = mess_vector_toreal_nowarn(vect);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        vect->data_type = MESS_REAL;
    } else if ( strncmp(datatype, MM_STR_INT, strlen(MM_STR_INT)) ==0 ) {
        ret = mess_vector_toreal_nowarn(vect);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
        vect->data_type = MESS_REAL;
    } else if ( strncmp(datatype, MM_STR_CPX, strlen(MM_STR_CPX)) == 0 ){
        ret = mess_vector_tocomplex(vect);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
        vect->data_type = MESS_COMPLEX;
    } else {
        MSG_ERROR("Unknown data type: %s\n",datatype);
        return (MESS_ERROR_DATATYPE);
    }

    do {
        cret = csc_io_gets(line, LINE_LENGTH, file);
        if ( cret == NULL ) {
            MSG_ERROR("error reading from file: %s \n",line);
            return (MESS_ERROR_FILEIO);
        }
    } while ( *line == '%' || strlen(line) <= 1 ) ;

    if ( sscanf(line, "" MESS_PRINTF_INT "", &dim) != 1 ){
        MSG_ERROR("error reading dimension from file: %s.\n", line);
        return (MESS_ERROR_FILEIO);
    }
    ret = mess_vector_resize(vect, dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize);
    vect->dim=dim;

    return 0 ;

}       /* -----  end of function mm_read_vector_info  ----- */

#ifdef MESS_HAVE_MATIO
static int __matlab_read(const char *filename, mess_vector vector ) {
    MSG_FNAME(__func__);
    int ret = 0;

    char *fn;
    char *var;

    mat_t *fp;
#if _MATIO_VERSION >= _MATIO_15
    enum mat_ft ft_ver;
#endif
    matvar_t *v;

    mess_int_t rows, cols, dim ;
    mess_int_t i;

    /*  Exract var and file name  */
    mess_try_alloc(fn, char*, (strlen(filename)+1) * sizeof(char));
    strncpy(fn, filename, strlen(filename)+1);
    var = strrchr(fn,':');
    *var = '\0';
    var ++;

    MSG_INFO("Read \"%s\" from %s\n",var, fn);

    /* Open matlab file.  */
    fp = Mat_Open(fn,MAT_ACC_RDONLY);
    if ( fp == NULL ) {
        MSG_ERROR("Error opening MAT file \"%s\"!\n",fn);
        mess_free(fn);
        return MESS_ERROR_FILEIO;
    }
#if _MATIO_VERSION >= _MATIO_15
    ft_ver = Mat_GetVersion(fp);
    switch(ft_ver){
        case MAT_FT_MAT73:
            MSG_INFO("Opened MATLAB v7.3 file\n");
            break;
        case MAT_FT_MAT5:
            MSG_INFO("Opened MATLAB v5 file\n");
            break;
        case MAT_FT_MAT4:
            MSG_INFO("Opened MATLAB v4 file\n");
            break;
        default:
            MSG_INFO("Opened Unknown MATLAB file.\n");
    }
#endif
    /*  Search for the variable  */
    v = Mat_VarRead(fp, var);
    if ( v == NULL ) {
        MSG_ERROR("Cannot find \"%s\" inside the file %s\n", var, fn );
        Mat_Close(fp);
        mess_free(fn);
        return MESS_ERROR_ARGUMENTS;
    }
    if (v->rank != 2 ) {
        MSG_ERROR("Rank != 2 in file %s for variable %s\n", fn, var);
        return MESS_ERROR_DATA;
    }
    if (v->dims[0] > 1 && v->dims[1]>1) {
        MSG_ERROR("Variable %s is not a vector.\n", var);
        Mat_VarFree(v);
        Mat_Close(fp);
        mess_free(fn);

    }
    rows = v->dims[0];
    cols = v->dims[1];
    if ( rows == 1) dim = cols;
    else dim = rows;

    if ( v->isComplex ) {
        ret = mess_vector_tocomplex(vector); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tocomplex);
    } else {
        ret = mess_vector_toreal_nowarn(vector); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_toreal_nowarn);
    }
    ret = mess_vector_resize(vector, dim); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_resize );

    switch(v->class_type) {
        case MAT_C_DOUBLE:
            if ( MESS_IS_REAL(vector)) {
                double * realv = v ->data;
                for (i = 0; i < dim; i++) {
                    vector->values[i] = realv[i];
                }
            } else if ( MESS_IS_COMPLEX(vector)) {
#if _MATIO_VERSION < _MATIO_15
                struct ComplexSplit *cpxdata = v->data;
#else
                mat_complex_split_t *cpxdata = v->data;
#endif
                double *realv = cpxdata->Re;
                double *iv    = cpxdata->Im;
                for (i = 0; i < dim; i++) {
                    vector->values_cpx[i] = realv[i] + iv[i]*I;
                }
            }
            break;
        case MAT_C_SINGLE:
            if ( MESS_IS_REAL(vector )) {
                float * realv = v ->data;
                for (i = 0; i < dim; i++) {
                    vector->values[i] = realv[i];
                }
            } else if ( MESS_IS_COMPLEX(vector)) {
#if _MATIO_VERSION < _MATIO_15
                struct ComplexSplit *cpxdata = v->data;
#else
                mat_complex_split_t *cpxdata = v->data;
#endif
                float  *realv = cpxdata->Re;
                float  *iv    = cpxdata->Im;
                for (i = 0; i < dim; i++) {
                    vector->values_cpx[i] = realv[i] + iv[i]*I;
                }
            }
            break;

        default:
            MSG_ERROR("The class of the variable %s is not supported.\n", var);
            return MESS_ERROR_DATATYPE;
    }

    Mat_VarFree(v);
    Mat_Close(fp);
    mess_free(fn);
    return 0;
}
#endif

/**
 * @brief Read a vector from a matrix market vector file.
 * @param[in]  filename  input filename of matrix market vector file
 * @param[in,out]  vect     vector storing read data
 * @return zero on success or a non zero error code
 *
 * The @ref mess_vector_read function reads a vector from a file. \n
 * The file need to be in @mm vector format, that means the header
 * need to be
 * <center>"%%MatrikMarket vector..." </center>
 * and the dimension is only one value.
 *
 */
int mess_vector_read ( const char * filename, mess_vector vect )
{
    MSG_FNAME(__func__);
    csc_io_file_t *f;
    int ret;
    int err = 0 ;
    mess_int_t i;

    mess_check_nullpointer(filename);
    mess_check_nullpointer(vect);
#ifdef __WIN32__
    if (( strrchr(filename, ':') != NULL) && (strrchr(filename, ':') !=  strchr(filename, ':')))
#else
        if ( strrchr(filename, ':') != NULL )
#endif
        {
#ifdef MESS_HAVE_MATIO
            ret = __matlab_read(filename, vect);
            if ( ret == 0 ) {
                return 0;
            }
            MSG_WARN("File does not seem to be a Matlab file. Try MatrixMarket/HB instead.\n");
#else
            MSG_WARN("Reading variables from Matlab files not support. Please enable MATIO support in MESS. Try MatrixMarket/HB instead.\n");
#endif
        }


    f = csc_io_open(filename, CSC_IO_FILE_READ);
    if ( !f ) {
        MSG_ERROR("error opening %s\n", filename);
        return (MESS_ERROR_FILEIO);
    }
    ret = __mm_read_vector_info(f, vect);               FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __mm_read_vector_info);
    if ( MESS_IS_REAL(vect)) {
        for ( i =0 ; i < vect->dim; i++){
            ret = csc_io_scanf(f, "%lg", &(vect->values[i]));
            if ( ret != 1 ){
                MSG_ERROR("cannot read two values for pattern format in line " MESS_PRINTF_INT "\n", i);
                err =  MESS_ERROR_DATA;
                break;
            }
        }
    } else if ( MESS_IS_COMPLEX(vect)) {
        double re, im;
        for ( i =0 ; i < vect->dim; i++){
            ret = csc_io_scanf(f, "%lg %lg", &(re), &im);
            if ( ret != 2 ){
                MSG_ERROR("cannot read two values for pattern format in line " MESS_PRINTF_INT "\n", i);
                err =  MESS_ERROR_DATA;
                break;
            }
            vect->values_cpx[i] = re + im*I;
        }
    }
    FUNCTION_FAILURE_HANDLE(err, (err!=0), RedingData);
    csc_io_close(f);


    return (0);
}       /* -----  end of function mess_read_vector  ----- */

/**
 * @brief Write a vector to a matrix market vector file.
 * @param[in] filename  input filename where to store the matrix
 * @param[in,out] vect     mess_vector to write
 * @return zero on success or a non zero error code
 *
 * The @ref mess_vector_write function writes a vector to a file. \n
 * The output is similar to a matrix market file but contains only
 * a one dimensional array. \n
 * As consequence  the header is
 * <center>"%%MatrixMarket vector ..." </center>
 * and the stored dimension is only one value.
 */
int mess_vector_write ( const char *filename, mess_vector vect)
{
    MSG_FNAME(__func__);
    csc_io_file_t *f;
    char line[80];
    mess_int_t i;
    int ret = 0;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(filename);
    mess_check_nullpointer(vect);

    f = csc_io_open ( filename, CSC_IO_FILE_WRITE);
    if ( !f ) {
        MSG_ERROR("error opening %s\n", filename);
        return (MESS_ERROR_FILEIO);
    }
    if (MESS_IS_REAL(vect)) snprintf(line,79,"%%%%MatrixMarket vector real\n");
    else if (MESS_IS_COMPLEX(vect))  snprintf(line,79,"%%%%MatrixMarket vector complex\n");
    else {
        MSG_ERROR("unknown data type: %s\n", mess_datatype_t_str(vect->data_type));
        csc_io_close(f);
        return (MESS_ERROR_DATATYPE);
    }

    ret = csc_io_puts(line,f);  FUNCTION_FAILURE_HANDLE(ret, (ret<=0), csc_io_puts);
    ret = csc_io_printf(f,"" MESS_PRINTF_INT "\n", vect->dim); FUNCTION_FAILURE_HANDLE(ret, (ret<=0), csc_io_printf);

    if (MESS_IS_REAL(vect)){
        for ( i=0; i < vect->dim ; i++ ) {
            csc_io_printf(f,"%.18e\n", vect->values[i]);
        }
    } else {
        for ( i=0; i < vect->dim ; i++ ) {
            csc_io_printf(f,"%.18e %.18e\n", creal(vect->values_cpx[i]),cimag(vect->values_cpx[i]));
        }
    }
    csc_io_close(f);
    return (0);
}       /* -----  end of function mess_vector_write  ----- */


/**
 * @brief Write a vector into a @matlab -mfile.
 * @param[in] filename   input file to write
 * @param[in] name   input name of variable in mfile
 * @param[in] v              input vector to write to the mfile
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_vector_mwrite function writes a vector as a @matlab  column
 * vector in to a file. \n
 * It could be loaded into @matlab  by calling the script.
 *
 */
int  mess_vector_mwrite ( const char *filename, const char * name, mess_vector v  )
{
    MSG_FNAME(__func__);
    csc_io_file_t *f;
    mess_int_t i;
    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(filename);
    mess_check_nullpointer(name);
    mess_check_nullpointer(v);
    mess_check_real_or_complex(v);

    /*-----------------------------------------------------------------------------
     * Write vector
     *-----------------------------------------------------------------------------*/
    f = csc_io_open(filename,CSC_IO_FILE_WRITE);
    if ( !f) {
        MSG_ERROR("error opening %s\n", filename);
        return (MESS_ERROR_FILEIO);
    }
    csc_io_printf(f,"%s = [ ", name);
    if (MESS_IS_REAL(v)){
        for (i=0; i < v->dim; i++){
            csc_io_printf(f,"%20.17e ", v->values[i]);
        }
    } else {
        for (i=0; i < v->dim; i++){
            csc_io_printf(f,"%20.17e+%20.17ei ", creal(v->values_cpx[i]), cimag(v->values_cpx[i]));
        }
    }
    csc_io_printf(f," ]';\n");
    csc_io_close(f);

    return (0);
}       /* -----  end of function mess_vector_mwrite  ----- */


