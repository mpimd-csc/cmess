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
 * @file lib/io/read.c
 * @brief Read matrices from files.
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

/* Compressed IO */
#include "cscutils/io.h"

#ifdef MESS_HAVE_MATIO
#include <matio.h>
#define _MATIO_VERSION ( (MATIO_MAJOR_VERSION) *10000 + (MATIO_MINOR_VERSION) * 100 + (MATIO_RELEASE_LEVEL))
#define _MATIO_15  10500
#endif

#define LINE_LENGTH 1025
#define TOKEN_LENGTH 30

/* Sets the last entry of an array to value */
#define     SET_LAST(array, len, value)     (array)[(len)-1] = (value)
/* Definitions for the Matrix Market header */
#define MM_BANNER           "%%MatrixMarket"
#define MM_STR_MATRIX       "matrix"
#define MM_STR_COORD        "coordinate"
#define MM_STR_ARR          "array"
#define MM_STR_REAL         "real"
#define MM_STR_INT          "integer"
#define MM_STR_CPX          "complex"
#define MM_STR_PAT          "pattern"
#define MM_STR_GE           "general"
#define MM_STR_SYM          "symmetric"
#define MM_STR_SKSYM        "skew-symmetric"
#define MM_STR_HER          "hermitian"

/**
 * @internal
 * @brief Convert a string to lower case.
 * @param[in,out] str  input/output string to be converted to lower case
 * @return Nothing.
 *
 * Convert a string to lowercase. This is done inplace.
 *
 * @attention Internal use only.
 */
static void _lowercase (char *str)
{
    char *p;
    if ( str == NULL) return;
    for (p=str; *p!='\0'; *p=tolower(*p),p++);
}

/**
 * @internal
 * @brief Convert a string to upper case.
 * @param[in,out] str  inout/output String to be converted to upper case
 * @return Nothing.
 *
 * Convert a string to upper case. This is done inplace.
 *
 * @attention Internal use only.
 */
static void _uppercase(char* str)
{
    char *p;
    if ( str == NULL) return;
    for (p=str; *p!='\0'; *p=toupper(*p),p++);
}


/**
 * @internal
 * @brief Copy a sub-string out of a given one.
 * @param[in] S      input Input string
 * @param[in] pos    input Start position of the requested sub-string
 * @param[in] len    input Length of the requested sub-string
 * @return A copy of the sub-string or @c NULL in case of an error.
 *
 *
 * The @ref _substr function copies a sub-string beginning at \p pos with \p len characters from \p S to a
 * newly allocated string. The returned string must be freed after it is not longer used.
 *
 * @attention Internal use only.
 */
static char* _substr(const char* S, const int pos, const int len)
{
    MSG_FNAME(__func__);
    int i;
    char *SubS;
    if ( pos+len <= strlen(S)) {
        mess_try_alloc2(SubS, char *, (len+1) * sizeof(char));
        if ( ! SubS) return NULL;
        for (i=0;i<len;i++){
            SubS[i] = S[pos+i];
        }
        SubS[len] = '\0';
    } else {
        SubS = NULL;
    }
    return SubS;
}

/*************************************************/
/*  Parse an *integer* format field to determine */
/*  width and number of elements per line.       */
/*************************************************/
static int _ParseIfmt(char* fmt, int* perline, int* width)
{
    char *tmp;
    int t;
    if (fmt == NULL ) {
        *perline = 0; *width = 0; return 0;
    }
    _uppercase(fmt);
    tmp = strchr(fmt,'(');
    tmp = _substr(fmt,tmp - fmt + 1, strchr(fmt,'I') - tmp - 1);
    *perline = atoi(tmp);
    mess_free(tmp);
    tmp = strchr(fmt,'I');
    tmp = _substr(fmt,tmp - fmt + 1, strchr(fmt,')') - tmp - 1);
    *width = atoi(tmp);
    t = *width;
    mess_free(tmp);
    return t;

}

/*************************************************/
/*  Parse a *real* format field to determine     */
/*  width and number of elements per line.       */
/*  Also sets flag indicating 'E' 'F' 'P' or 'D' */
/*  format.                                      */
/*************************************************/
static int _ParseRfmt(char* fmt, int* perline, int* width, int* prec, int* flag)
{
    MSG_FNAME(__func__);
    char* tmp;
    char* tmp2;
    char* tmp3;
    int len;
    int t;

    if (fmt == NULL ) {
        *perline = 0;
        *width = 0;
        *flag = 0;
        return 0;
    }

    _uppercase(fmt);
    if (strchr(fmt,'(') != NULL)  fmt = strchr(fmt,'(');
    if (strchr(fmt,')') != NULL)  {
        tmp2 = strchr(fmt,')');
        while ( strchr(tmp2+1,')') != NULL ) {
            tmp2 = strchr(tmp2+1,')');
        }
        tmp2[1] = '\0';
    }
    if (strchr(fmt,'P') != NULL)  /* Remove any scaling factor, which */
    {                             /* affects output only, not input */
        if (strchr(fmt,'(') != NULL) {
            tmp = strchr(fmt,'P');
            if ( *(++tmp) == ',' ) tmp++;
            tmp3 = strchr(fmt,'(')+1;
            len = tmp-tmp3;
            tmp2 = tmp3;
            while ( *(tmp2+len) != '\0' ) {
                *tmp2=*(tmp2+len);
                tmp2++;
            }
            *(strchr(fmt,')')+1) = '\0';
        }
    }
    if (strchr(fmt,'E') != NULL) {
        *flag = 'E';
    } else if (strchr(fmt,'D') != NULL) {
        *flag = 'D';
    } else if (strchr(fmt,'F') != NULL) {
        *flag = 'F';
    } else {
        MSG_ERROR("Real format %s in H/B file not supported.\n",fmt);
        return 0;
    }
    tmp = strchr(fmt,'(');
    tmp = _substr(fmt,tmp - fmt + 1, strchr(fmt,*flag) - tmp - 1);
    *perline = atoi(tmp);
    mess_free(tmp);
    tmp = strchr(fmt,*flag);
    if ( strchr(fmt,'.') ) {
        tmp2 = _substr( fmt, strchr(fmt,'.') - fmt + 1, strchr(fmt,')') - strchr(fmt,'.')-1);
        *prec = atoi(tmp2 );
        mess_free(tmp2);
        tmp = _substr(fmt,tmp - fmt + 1, strchr(fmt,'.') - tmp - 1);
    } else {
        tmp = _substr(fmt,tmp - fmt + 1, strchr(fmt,')') - tmp - 1);
    }
    t = atoi(tmp);
    mess_free(tmp);
    return *width = t;
}


/**
 * @internal
 * @brief Read the header of a Matrix Market file.  (internal version)
 * @param[in] file  input filepointer
 * @param[out] matrix_info mess_matrix where the information should be saved
 * @return zero on success, otherwise a non zero error code.
 *
 * The @ref __mm_read_info function reads the header of a @mm file and writes the
 * basic information to matrix_info.
 *
 * @attention Internal use only.
 */
static int __mm_read_info(csc_io_file_t *file, mess_matrix matrix_info)
{
    MSG_FNAME(__func__);
    char line[LINE_LENGTH];
    char mmbanner[TOKEN_LENGTH];
    char mtx[TOKEN_LENGTH];
    char type[TOKEN_LENGTH];
    char datatype[TOKEN_LENGTH];
    char shape[TOKEN_LENGTH];
    int num_items_read;
    char  * cret =NULL;

    mess_int_t N,M,NZ;

    mess_check_nullpointer(file);
    mess_check_nullpointer(matrix_info);


    // %%MatrixMarket matrix type datatype shape
    if (csc_io_scanf(file, "%s %s %s %s %s", mmbanner, mtx, type, datatype, shape) != 5){
        // MSG_ERROR("wrong information in header: %s\n", line);
        return MESS_ERROR_WRONG_HEADER;
    }

    // prevent memory overflow
    SET_LAST(mmbanner, TOKEN_LENGTH, '\0');
    SET_LAST(mtx, TOKEN_LENGTH, '\0');
    SET_LAST(type, TOKEN_LENGTH, '\0');
    SET_LAST(datatype, TOKEN_LENGTH, '\0');
    SET_LAST(shape, TOKEN_LENGTH, '\0');

    _lowercase(mtx);
    _lowercase(type);
    _lowercase(datatype);
    _lowercase(shape);

    /* check if the matrix market banner exists */
    if (strncmp(mmbanner, MM_BANNER, strlen(MM_BANNER)) != 0){
        // MSG_ERROR("wrong header information: %s\n", mmbanner);
        return MESS_ERROR_WRONG_HEADER;
    }

    /* parse the rest of the 1st line*/
    if (strncmp(mtx, MM_STR_MATRIX, strlen(MM_STR_MATRIX)) != 0){
        // MSG_ERROR("wrong header information: %s\n", mtx);
        return MESS_ERROR_WRONG_HEADER;
    }

    if (strncmp(type, MM_STR_ARR,strlen(MM_STR_ARR)) == 0){
        MSG_INFO("read a dense matrix\n");
        matrix_info->store_type = MESS_DENSE;
    } else if (strncmp(type, MM_STR_COORD, strlen(MM_STR_COORD)) == 0){
        MSG_INFO("read a sparse matrix\n");
        matrix_info->store_type = MESS_COORD;
    } else {
        MSG_ERROR("wrong storage format: %s\n", type );
        return MESS_ERROR_WRONG_HEADER;
    }

    if (strncmp(datatype, MM_STR_REAL,strlen(MM_STR_REAL)) == 0) {
        MSG_INFO(" - real entries\n");
        matrix_info->data_type = MESS_REAL;
    } else if (strncmp(datatype, MM_STR_INT, strlen(MM_STR_INT)) == 0){
        // Integermatrices are treated as real-value matrices
        matrix_info->data_type = MESS_REAL;
        MSG_INFO(" - real(integer) entries\n");
    } else if (strncmp(datatype, MM_STR_CPX, strlen(MM_STR_CPX)) == 0){
        matrix_info->data_type = MESS_COMPLEX;
        MSG_INFO(" - complex entries\n");
    }
    else if (strncmp(datatype, MM_STR_PAT, strlen(MM_STR_PAT)) == 0) {
        MSG_INFO(" - only pattern\n");
        matrix_info->data_type = MESS_REAL;
    } else {
        MSG_ERROR("wrong datatype: %s\n", datatype);
        return MESS_ERROR_WRONG_HEADER;
    }

    if (strncmp(shape, MM_STR_GE,strlen(MM_STR_GE)) == 0){
        MSG_INFO(" - General\n");
        matrix_info->symmetry = MESS_GENERAL;
    } else if (strncmp(shape, MM_STR_SYM, strlen(MM_STR_SYM)) == 0){
        MSG_INFO(" - Symmetric\n");
        matrix_info->symmetry = MESS_SYMMETRIC;
    } else if (strncmp(shape, MM_STR_SKSYM, strlen(MM_STR_SKSYM)) == 0) {
        MSG_INFO(" - Skew-Symmetric\n");
        matrix_info->symmetry = MESS_SKEWSYMMETRIC;
    } else if (strncmp(shape, MM_STR_HER, strlen(MM_STR_HER)) == 0) {
        MSG_INFO(" - Hermitian\n");
        matrix_info->symmetry = MESS_HERMITIAN;
    }  else {
        MSG_ERROR("wrong symmetric property: %s\n", shape);
        return MESS_ERROR_WRONG_HEADER;
    }

    // read until matrix size is read.
    do{
        cret = csc_io_gets(line, LINE_LENGTH, file );
        if ( cret == NULL  ) {
            MSG_ERROR("error while reading a line from file");
            return MESS_ERROR_FILEIO;
        }
    } while (line[0] == '%');


    if (MESS_IS_COORD(matrix_info)){
        do {
            num_items_read = sscanf(line, "" MESS_PRINTF_INT " " MESS_PRINTF_INT " " MESS_PRINTF_INT "", &N, &M, &NZ);
            if (num_items_read != 3) {
                cret = csc_io_gets(line, LINE_LENGTH, file);
                if ( cret == NULL ) {
                    MSG_ERROR("error while reading a line from file\n");
                    return MESS_ERROR_FILEIO;
                }
            }
        } while (num_items_read != 3);
        MSG_INFO("read sparse matrix - rows = " MESS_PRINTF_INT ", cols = " MESS_PRINTF_INT ", non-zeros = " MESS_PRINTF_INT "\n", N, M, NZ);
    }else{
        NZ = 0;
        do {
            num_items_read = sscanf(line,"" MESS_PRINTF_INT " " MESS_PRINTF_INT "", &N, &M);
            if (num_items_read != 2) {
                cret = csc_io_gets(line, LINE_LENGTH, file);
                if ( cret == NULL   ) {
                    MSG_ERROR("error while reading a line from file\n");
                    return MESS_ERROR_FILEIO;
                }
            }
        } while (num_items_read != 2);
        NZ = N * M;
        MSG_INFO("read dense matrix - rows = " MESS_PRINTF_INT ", cols = " MESS_PRINTF_INT ", non-zeros = " MESS_PRINTF_INT "\n", N, M, NZ);
    }

    matrix_info->rows = N;
    matrix_info->cols = M;
    matrix_info->nnz = NZ;
    return 0;
}

/**
 * @internal
 * @brief Read the data part of a @mm matrix.
 * @param[in]  f     input File to read from.
 * @param[in,out] matrix    Matrix object to read.
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __mm_read_data function reads the data part out of a @mm file.
 * The \p matrix object has to be filled up with the values read by \ref __mm_read_info.
 *
 * @attention Internal use only.
 */
static int __mm_read_data(csc_io_file_t *f, mess_matrix matrix) {
    MSG_FNAME(__func__);
    int ret;
    mess_int_t i, j = 0;
    mess_symmetry_t sym  = MESS_GENERAL;

    /*  Read matrix market  */
    if (MESS_IS_COORD(matrix)) {
        sym = matrix->symmetry;
        ret = mess_matrix_alloc ( matrix, matrix->rows, matrix->cols, matrix->nnz, MESS_COORD, matrix->data_type);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        matrix->symmetry = sym;
        if (MESS_IS_REAL( matrix ) ) {
            for (i=0; i < matrix->nnz; i++){
                ret = csc_io_scanf(f, "" MESS_PRINTF_INT " " MESS_PRINTF_INT " %lg\n", &(matrix->rowptr[i]), &(matrix->colptr[i]), &(matrix->values[i]));
                if ( ret != 3 ){
                    MSG_ERROR("cannot read value for real format in line " MESS_PRINTF_INT "\n", i);
                    return MESS_ERROR_DATA;
                }
                /* adjust to C-style indexing */
                matrix->rowptr[i]--;
                matrix->colptr[i]--;
            }
        }
        else if (MESS_IS_COMPLEX(matrix)) {
            for (i=0; i < matrix->nnz; i++){
                double re, im;
                ret = csc_io_scanf(f, "" MESS_PRINTF_INT " " MESS_PRINTF_INT " %lg %lg\n", &(matrix->rowptr[i]), &(matrix->colptr[i]), &(re), &(im));
                if ( ret != 4 ){
                    MSG_ERROR("cannot read two values for complex format in line " MESS_PRINTF_INT "\n", i);
                    return MESS_ERROR_DATA;
                }
                matrix->values_cpx[i] = re + im*I;
                /* adjust to C-style indexing */
                matrix->rowptr[i]--;
                matrix->colptr[i]--;
            }
        } else {
            MSG_ERROR("wrong data type: %s\n", mess_datatype_t_str(matrix->data_type));
            return MESS_ERROR_DATATYPE;
        }
        if ( MESS_IS_SYMMETRIC(matrix) ) {
            MSG_INFO ( "read symmetric matrix, convert to a generic one\n");
            ret = mess_matrix_symfillup(matrix);
            FUNCTION_FAILURE_HANDLE(ret , (ret!=0), mess_matrix_symfillup);
        }
    } else if ( MESS_IS_DENSE(matrix) && MESS_IS_GENERAL(matrix)){
        ret = mess_matrix_alloc ( matrix, matrix->rows, matrix->cols, matrix->nnz, MESS_DENSE, matrix->data_type);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        if ( MESS_IS_REAL(matrix)) {
            for ( j = 0; j < matrix->cols ; j++) {
                for ( i = 0; i < matrix->rows ; i++){
                    ret = csc_io_scanf(f, "%lg\n", &(matrix->values[i+j*matrix->ld]));
                    if ( ret != 1){
                        MSG_ERROR("cannot read two values for real format in line " MESS_PRINTF_INT "\n", i);
                        return MESS_ERROR_DATA;
                    }
                }
            }
        } else if ( MESS_IS_COMPLEX(matrix)){
            for ( j = 0 ; j < matrix->cols; j++){
                for (i = 0; i < matrix->rows; i++){
                    double re, im;
                    ret = csc_io_scanf(f, "%lg %lg\n", &(re), &(im));
                    if ( ret != 2 ){
                        MSG_ERROR("cannot read two values for complex format in line " MESS_PRINTF_INT "\n", i);
                        return   MESS_ERROR_DATA;
                    }
                    matrix->values_cpx[i+j*matrix->ld] = re + im*I;
                }
            }
        } else {
            MSG_ERROR("wrong data type\n");
            return MESS_ERROR_DATATYPE;
        }
    } else if ( MESS_IS_DENSE(matrix) && (MESS_IS_SYMMETRIC(matrix) || MESS_IS_HERMITIAN(matrix) || MESS_IS_SKEWSYMMETRIC(matrix))){
        mess_int_t colstart = 0;
        sym = matrix->symmetry;
        ret = mess_matrix_alloc ( matrix, matrix->rows, matrix->cols, matrix->nnz, MESS_DENSE, matrix->data_type);
        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);
        matrix->symmetry = sym;
        if ( MESS_IS_SKEWSYMMETRIC(matrix)) {
            ret = mess_matrix_zeros(matrix);    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_zeros);
            colstart = 1;
        }
        if ( MESS_IS_REAL(matrix)) {
            for ( j = colstart; j < matrix->cols ; j++) {
                for ( i = j; i < matrix->rows ; i++){
                    ret = csc_io_scanf(f, "%lg\n", &(matrix->values[i+j*matrix->ld]));
                    if ( ret != 1){
                        MSG_ERROR("cannot read two values for real format in line " MESS_PRINTF_INT "\n", i);
                        return MESS_ERROR_DATA;
                    }
                    if ( MESS_IS_SKEWSYMMETRIC(matrix) ) {
                        matrix->values[j+i*matrix->ld] =  -matrix->values[i+j*matrix->ld];
                    } else {
                        matrix->values[j+i*matrix->ld] =  matrix->values[i+j*matrix->ld];
                    }
                }
            }
        } else if ( MESS_IS_COMPLEX(matrix)){
            for ( j = colstart ; j < matrix->cols; j++){
                for (i = j; i < matrix->rows; i++){
                    double re, im;
                    ret = csc_io_scanf(f, "%lg %lg\n", &(re), &(im));
                    if ( ret != 2 ){
                        MSG_ERROR("cannot read two values for complex format in line " MESS_PRINTF_INT "\n", i);
                        return   MESS_ERROR_DATA;
                    }
                    if ( MESS_IS_SKEWSYMMETRIC(matrix)) {
                        matrix->values_cpx[i+j*matrix->ld] = re + im*I;
                        matrix->values_cpx[j+i*matrix->ld] = -conj(re + im*I);
                    } else {
                        matrix->values_cpx[i+j*matrix->ld] = re + im*I;
                        matrix->values_cpx[j+i*matrix->ld] = conj(re + im*I);
                    }
                }
            }
        } else {
            MSG_ERROR("wrong data type\n");
            return MESS_ERROR_DATATYPE;
        }

    }else {
        MSG_ERROR ("unkown data storage format: %s\n", mess_storage_t_str(matrix->store_type));
        return MESS_ERROR_STORAGETYPE;
    }
    return 0;
}


/**
 * @internal
 * @brief Read the header information out of a @hb matrix.
 * @param[in] in_file   input file to read from
 * @param Type
 * @param Nrow
 * @param Ncol
 * @param Nnzero
 * @param Nrhs
 * @param Ptrfmt
 * @param Indfmt
 * @param Valfmt
 * @param Rhsfmt
 * @param Ptrcrd
 * @param Indcrd
 * @param Valcrd
 * @param Rhscrd
 * @param Rhstype
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref __hb_read_info function reads the header information out of a @hb
 * matrix.
 *
 * @attention Internal use only.
 */

static int __hb_read_info(csc_io_file_t * in_file, char* Type,
        int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
        int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd,
        char *Rhstype)
{
    /*************************************************************************/
    /*  Read header information from the named H/B file...                   */
    /*************************************************************************/
    MSG_FNAME(__func__);
    int Totcrd,Neltvl,Nrhsix;
    char *ret = 0 ;
    char line[LINE_LENGTH];
    char Title[80];
    char Key[10];



    /*  First line:   */
    if ( csc_io_scanf(in_file, "%72c%8[^\n]", Title, Key) != 2 ) {
        MSG_ERROR("First line did not match Title and Key.\n");
        return MESS_ERROR_DATA;
    }
    Key[8] = '\0';
    Title[72] ='\0';
    MSG_INFO("Read Matrix Title = \"%s\", Key = \"%s\"\n", Title, Key);

    /*  Second line:  */
    ret = csc_io_gets(line, LINE_LENGTH, in_file);
    if ( ret == NULL ) {
        MSG_ERROR("Cannot read line from file\n");
        return MESS_ERROR_FILEIO;
    }
    if ( sscanf(line,"%i",&Totcrd) != 1) Totcrd = 0;
    if ( sscanf(line,"%*i%i",Ptrcrd) != 1) *Ptrcrd = 0;
    if ( sscanf(line,"%*i%*i%i",Indcrd) != 1) *Indcrd = 0;
    if ( sscanf(line,"%*i%*i%*i%i",Valcrd) != 1) *Valcrd = 0;
    if ( sscanf(line,"%*i%*i%*i%*i%i",Rhscrd) != 1) *Rhscrd = 0;
    MSG_INFO("Totcrd = %d, Ptrcrd = %d, Indcrd = %d, Valcrd = %d, Rhscrd = %d\n", Totcrd, *Ptrcrd, *Indcrd, *Valcrd, *Rhscrd);

    /*  Third line:   */
    ret = csc_io_gets( line, LINE_LENGTH, in_file);
    if ( ret == NULL ) {
        MSG_ERROR("Cannot read line from file\n");
        return MESS_ERROR_FILEIO;
    }

    if ( sscanf(line, "%3c", Type) != 1) {
        MSG_ERROR("Cannot read Type info from HB file.\n");
        return MESS_ERROR_DATA;
    }
    Type[3] = '\0';
    _uppercase(Type);
    if ( sscanf(line,"%*3c%i",Nrow) != 1) *Nrow = 0 ;
    if ( sscanf(line,"%*3c%*i%i",Ncol) != 1) *Ncol = 0 ;
    if ( sscanf(line,"%*3c%*i%*i%i",Nnzero) != 1) *Nnzero = 0 ;
    if ( sscanf(line,"%*3c%*i%*i%*i%i",&Neltvl) != 1) Neltvl = 0 ;
    MSG_INFO("Type: %s, Nrow = %d, Ncol = %d, Nnzero = %d, Neltvl = %d\n", Type, *Nrow, *Ncol, *Nnzero, Neltvl);


    /*  Fourth line:  */
    ret = csc_io_gets(line, LINE_LENGTH,in_file);
    if ( ret == NULL ) {
        MSG_ERROR("Cannot read line from file\n");
        return MESS_ERROR_FILEIO;
    }

    if ( sscanf(line, "%16c",Ptrfmt) != 1) {
        MSG_ERROR("Cannot read Pointer format.\n");
        return MESS_ERROR_DATA;
    }
    if ( sscanf(line, "%*16c%16c",Indfmt) != 1){
        MSG_ERROR("Cannot read Indices format.\n");
        return MESS_ERROR_DATA;
    }
    if ( sscanf(line, "%*16c%*16c%20c",Valfmt) != 1) {
        MSG_ERROR("Cannot read Value format.\n");
        return MESS_ERROR_DATA;
    }
    sscanf(line, "%*16c%*16c%*20c%20c",Rhsfmt);
    Ptrfmt[16] = '\0';
    Indfmt[16] = '\0';
    Valfmt[20] = '\0';
    Rhsfmt[20] = '\0';
    MSG_INFO("Ptrfmt: %s, Indfmt: %s, Valfmt: %s, Rhsfmt: %s\n",Ptrfmt, Indfmt, Valfmt, Rhsfmt);

    /*  (Optional) Fifth line: */
    if (*Rhscrd != 0 )
    {
        ret = csc_io_gets(line, LINE_LENGTH, in_file);
        if ( ret == NULL ) {
            MSG_ERROR("Cannot read line from file\n");
            return MESS_ERROR_FILEIO;
        }

        if ( sscanf(line, "%3c", Rhstype) != 1) {
            MSG_ERROR("Cannot read the RHS type information. \n");
            return MESS_ERROR_DATA;
        }
        Rhstype[3]='\0';
        if ( sscanf(line, "%*3c%i", Nrhs) != 1) *Nrhs = 0;
        if ( sscanf(line, "%*3c%*i%i", &Nrhsix) != 1) Nrhsix = 0;
        MSG_INFO("RHS Type: %s, Nrhs = %d, Nrhsix = %d\n", Rhstype, *Nrhs, Nrhsix);
    }
    return 0;
}

static int __hb_read_data(csc_io_file_t *in_file, char *Type, mess_int_t Nrow, mess_int_t Ncol,
        mess_int_t Nnzero, mess_int_t Nrhs,
        char *Ptrfmt, char *Indfmt, char *Valfmt, char *Rhsfmt,
        mess_int_t Ptrcrd, mess_int_t Indcrd, mess_int_t Valcrd, mess_int_t Rhscrd,
        char *Rhstype,
        mess_int_t *colptr,  mess_int_t *rowind,  double *val)
{
    MSG_FNAME(__func__);
    mess_int_t i,j,ind,col,offset,count,last;
    mess_int_t Nentries;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Valperline=0, Valwidth=0, Valprec=0;
    int Valflag='0';           /* Indicates 'E','D', or 'F' float format */
    char* ThisElement;
    char line[LINE_LENGTH];
    char *ret = NULL;
    char *bu_ThisElement;

    /*  Parse the array input formats from Line 3 of HB file  */
    _ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    _ParseIfmt(Indfmt,&Indperline,&Indwidth);
    _ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
    if (Valflag == 'D') {
        *strchr(Valfmt,'D') = 'E';
    }

    offset = 1;  /* if base 0 storage is declared (via macro definition), */
    /* then storage entries are offset by 1                  */
    mess_try_calloc( ThisElement, char *, Ptrwidth+1, sizeof(char));
    count=0;
    for (i=0;i<Ptrcrd;i++)
    {
        ret = csc_io_gets(line, LINE_LENGTH, in_file);
        if ( ret == NULL ) {
            MSG_ERROR("Cannot read line from file\n");
            return MESS_ERROR_FILEIO;
        }

        col =  0;
        for (ind = 0;ind<Ptrperline;ind++)
        {
            if (count > Ncol) break;
            strncpy(ThisElement,line+col,Ptrwidth);
            /*ThisElement = substr(line,col,Ptrwidth);*/
            colptr[count] = atoi(ThisElement)-offset;
            count++; col += Ptrwidth;
        }
    }
    mess_free(ThisElement);

    /*  Read row index array:  */
    mess_try_calloc( ThisElement, char *, Indwidth+1, sizeof(char));
    count = 0;
    for (i=0;i<Indcrd;i++)
    {
        ret = csc_io_gets(line, LINE_LENGTH, in_file);
        if ( ret == NULL ) {
            MSG_ERROR("Cannot read line from file\n");
            return MESS_ERROR_FILEIO;
        }

        col =  0;
        for (ind = 0;ind<Indperline;ind++)
        {
            if (count == Nnzero) break;
            strncpy(ThisElement,line+col,Indwidth);
            /*ThisElement = substr(line,col,Indwidth);*/
            rowind[count] = atoi(ThisElement)-offset;
            count++; col += Indwidth;
        }
    }
    mess_free(ThisElement);

    /*  Read array of values:  AS CHARACTERS*/
    if ( Type[0] == 'C' ) Nentries = 2*Nnzero;
    else Nentries = Nnzero;

    mess_try_calloc( ThisElement, char *, Valwidth+1, sizeof(char));
    bu_ThisElement = ThisElement;
    count = 0;
    for (i=0;i<Valcrd;i++)
    {
        ret = csc_io_gets(line, LINE_LENGTH, in_file);
        if ( ret == NULL ) {
            MSG_ERROR("Cannot read line from file\n");
            return MESS_ERROR_FILEIO;
        }


        if (Valflag == 'D') {
            while( strchr(line,'D') ) *strchr(line,'D') = 'E';
        }
        col =  0;
        for (ind = 0;ind<Valperline;ind++)
        {
            if (count == Nentries) break;
            strncpy(ThisElement,line+col,Valwidth);
            /*ThisElement = substr(line,col,Valwidth);*/
            if ( Valflag != 'F' && strchr(ThisElement,'E') == NULL ) {
                /* insert a char prefix for exp */
                last = strlen(ThisElement);
                for (j=last+1;j>=0;j--) {
                    ThisElement[j] = ThisElement[j-1];
                    if ( ThisElement[j] == '+' || ThisElement[j] == '-' ) {
                        ThisElement[j-1] = Valflag;
                        break;
                    }
                }
            }
            val[count] = atof(ThisElement);
            count++; col += Valwidth;
        }
    }
    mess_free(bu_ThisElement);
    return 0;
}



static int __hb_read(csc_io_file_t *in_file, mess_matrix matrix) {
    MSG_FNAME(__func__);
    char  Type[4];
    int Nrow, Ncol, Nnzero, Nrhs;
    char Ptrfmt[17], Indfmt[17];
    char Valfmt[21], Rhsfmt[21];
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    char Rhstype[4];
    int ret = 0;

    mess_datatype_t data_type;
    mess_symmetry_t sym;

    /*-----------------------------------------------------------------------------
     *  Check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(in_file);
    mess_check_nullpointer(matrix);

    ret =  __hb_read_info( in_file, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
            Ptrfmt, Indfmt, Valfmt, Rhsfmt,
            &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd,
            Rhstype);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __hb_read_info);
    if ( Nrow <= 0) {
        MSG_ERROR("Read a matrix with less than 1 rows.\n");
        return MESS_ERROR_DIMENSION;
    }
    if ( Ncol <= 0) {
        MSG_ERROR("Read a matrix with less than 1 columns.\n");
        return MESS_ERROR_DIMENSION;
    }

    if ( Type[0] == 'R' ){
        data_type = MESS_REAL;
    } else if ( Type[0] == 'C' ) {
        data_type = MESS_COMPLEX;
    } else if ( Type[0] == 'P' ){
        MSG_ERROR("Pattern matrices are not supported by MESS\n");
        return MESS_ERROR_NOPATTERN;
    } else {
        MSG_ERROR("Unknown data type.\n");
        return MESS_ERROR_DATATYPE;
    }
    if ( Type[1] == 'U' || Type[1] == 'R' ){
        sym = MESS_GENERAL;
    } else if ( Type[1] == 'S' ){
        sym = MESS_SYMMETRIC;
    } else if ( Type[1] == 'Z' ) {
        sym = MESS_SKEWSYMMETRIC;
    } else if ( Type[1] == 'H' ){
        sym = MESS_HERMITIAN;
    } else {
        MSG_ERROR("Symmetry data field is invalid.");
        return MESS_ERROR_SYMMETRIC;
    }
    if ( Type[2] != 'A' ){
        MSG_ERROR("Only assembled HB matrices are supported.\n");
        return MESS_ERROR_NOSUPPORT;
    }

    /* Allocate matrix   */
    ret = mess_matrix_alloc( matrix, (mess_int_t) Nrow, (mess_int_t) Ncol, (mess_int_t) Nnzero, MESS_CSC, data_type);
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);

    ret = __hb_read_data(in_file, Type, Nrow, Ncol, Nnzero, Nrhs,
            Ptrfmt, Indfmt, Valfmt, Rhsfmt,
            Ptrcrd, Indcrd, Valcrd, Rhscrd,
            Rhstype,
            matrix->colptr, matrix->rowptr, (matrix->values==NULL)?((double *)matrix->values_cpx):(matrix->values));
    FUNCTION_FAILURE_HANDLE(ret, (ret!=0), __hb_read_data);

    /* Fill up symmetry  */
    matrix->symmetry = sym;

    if ( MESS_IS_SYMMETRIC(matrix) ) {
        MSG_WARN ( "read symmetric matrix, convert to a generic one\n");
        ret = mess_matrix_symfillup(matrix);
        FUNCTION_FAILURE_HANDLE(ret , (ret!=0), mess_matrix_symfillup);
    }
    return 0;
}

#ifdef MESS_HAVE_MATIO
static int __matlab_read(const char *filename, mess_matrix matrix ) {
    MSG_FNAME(__func__);
    int ret = 0;

    char *fn;
    char *var;

    mat_t *fp;
    matvar_t *v;

    mess_int_t rows, cols, nzeros;
    mess_storage_t store_type;
    mess_datatype_t data_type;
    mess_int_t r,c;
#if _MATIO_VERSION < _MATIO_15
    struct sparse_t *sparse = NULL;
#else
    enum mat_ft ft_ver;
    mat_sparse_t *sparse = NULL;
#endif

    /*  Exract var and file name  */
    mess_try_alloc(fn , char *, (strlen(filename)+1) * sizeof(char));
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
    /* if ( Mat_VarReadDataAll(fp, v)) {
       MSG_ERROR("Cannot read the data for  \"%s\" from %s\n", var, fn );
       mess_free(fn);
       return MESS_ERROR_FILEIO;
       } */

    /* Extract basic information   */
    /* MSG_PRINT("Ranks: %d\n", v->rank);
       int i;
       for (i = 0 ; i < v->rank; i++){
       MSG_PRINT("%d - %d\n", i ,(int) v->dims[i]);
       } */
    if (v->rank != 2 ) {
        MSG_ERROR("Rank != 2 in file %s for variable %s\n", fn, var);
        Mat_VarFree(v);
        Mat_Close(fp);
        mess_free(fn);
        return MESS_ERROR_DATA;
    }
    rows = v->dims[0];
    cols = v->dims[1];
    if ( v->isComplex ) {
        data_type = MESS_COMPLEX;
    } else {
        data_type = MESS_REAL;
    }
    if ( v->class_type == MAT_C_SPARSE) {
        sparse =  v->data;
        store_type = MESS_CSC;
        nzeros = sparse->ndata;
    } else {
        store_type = MESS_DENSE;
        nzeros = rows *cols;
    }

    ret = mess_matrix_alloc(matrix, rows, cols, nzeros, store_type, data_type); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_alloc);


    switch(v->class_type) {
        case MAT_C_SPARSE:
            sparse = v->data;
            for ( c = 0; c<cols+1; c++) {
                matrix->colptr[c] = (mess_int_t )sparse->jc[c];
            }
            for (r = 0; r < nzeros; r++) {
                matrix->rowptr[r] = (mess_int_t )sparse->ir[r];
            }
            if (MESS_IS_REAL(matrix)) {
                double *realv = sparse->data;
                for (r = 0; r < nzeros; r++) {
                    matrix->values[r] = realv[r];
                }
            } else if (MESS_IS_COMPLEX(matrix)){
#if _MATIO_VERSION < _MATIO_15
                struct ComplexSplit *cpxdata = sparse->data;
#else
                mat_complex_split_t *cpxdata = sparse->data;
#endif
                double *realv = cpxdata->Re;
                double *iv    = cpxdata->Im;
                for (r = 0; r < nzeros; r++) {
                    matrix->values_cpx[r] = realv[r]+iv[r]*I;
                }
            }
            break;
        case MAT_C_DOUBLE:
            if ( MESS_IS_REAL(matrix)) {
                double * realv = v ->data;
                for (c = 0; c < cols; c++) {
                    for (r = 0; r < rows; r++) {
                        matrix->values[r+c*matrix->ld] = realv[r+c*rows];
                    }
                }
            } else if ( MESS_IS_COMPLEX(matrix)) {
#if _MATIO_VERSION < _MATIO_15
                struct ComplexSplit *cpxdata = v->data;
#else
                mat_complex_split_t *cpxdata = v->data;
#endif
                double *realv = cpxdata->Re;
                double *iv    = cpxdata->Im;
                for (c = 0; c < cols; c++) {
                    for (r = 0; r < rows; r++) {
                        matrix->values_cpx[r+c*matrix->ld] = realv[r+c*rows] + iv[r+c*rows]*I;
                    }
                }
            }
            break;
        case MAT_C_SINGLE:
            if ( MESS_IS_REAL(matrix)) {
                float * realv = v ->data;
                for (c = 0; c < cols; c++) {
                    for (r = 0; r < rows; r++) {
                        matrix->values[r+c*matrix->ld] = (double) realv[r+c*rows];
                    }
                }
            } else if ( MESS_IS_COMPLEX(matrix)) {
#if _MATIO_VERSION < _MATIO_15
                struct ComplexSplit *cpxdata = v->data;
#else
                mat_complex_split_t *cpxdata = v->data;
#endif
                float *realv = cpxdata->Re;
                float  *iv    = cpxdata->Im;
                for (c = 0; c < cols; c++) {
                    for (r = 0; r < rows; r++) {
                        matrix->values_cpx[r+c*matrix->ld] = ((double)realv[r+c*rows]) + ((double)iv[r+c*rows])*I;
                    }
                }
            }
            break;
        default:
            MSG_ERROR("The class of the variable %s is not supported.\n", var);
            Mat_VarFree(v);
            Mat_Close(fp);
            mess_free(fn);
            return MESS_ERROR_DATATYPE;
    }

    Mat_VarFree(v);
    Mat_Close(fp);
    mess_free(fn);
    return 0;
}
#endif


/**
 * @brief Read a matrix from a file into a mess_matrix structure.
 * @param[in] filename  input file to be read
 * @param[out] matrix structure storing all information
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref mess_matrix_read function reads a matrix from a file. \n
 * The supported file formats are @mm, assembled Harwell-Boeing and @matlab (if @mess is compiled with @matio
 * support). \n
 * If @mess supports @zlib and @bzip, files can be compressed as well and it is automatically detected
 * if it is a compressed file or not.\n
 * In the case of a @matlab  file, the filename need to be formed like
 * <center> FILENAME:VARIABLENAME </center>
 * in order to select the variable to read. \n
 * If the matrix is sparse, the output will be a coordinate matrix.\n
 * Dense matrices are read directly into a dense one. \n
 * If a matrix has a symmetry property the missing parts are filled up and the matrix
 * becomes a general full stored one. \n
 * If a matrix should be read to a given storage format use \ref mess_matrix_read_formated instead.
 *
 *
 * @sa mess_matrix_write
 * @sa mess_matrix_symfillup
 */
int mess_matrix_read(const char *filename, mess_matrix matrix)
{
    csc_io_file_t *f;
    MSG_FNAME(__func__);
    int ret;


    mess_check_nullpointer(matrix);
    mess_check_nullpointer(filename);
    MESS_MATRIX_RESET(matrix);

#ifdef ___WIN32__
    if (( strrchr(filename, ':') != NULL) && (strrchr(filename, ':') !=  strchr(filename, ':')))
#else
        if ( strrchr(filename, ':') != NULL )
#endif
        {
#ifdef MESS_HAVE_MATIO
            ret = __matlab_read(filename, matrix);
            if ( ret == 0 ) {
                return 0;
            }
            MSG_WARN("File does not seem to be a Matlab file. Try MatrixMarket/HB instead.\n");
#else
            MSG_WARN("Reading variables from Matlab files not support. Please enable MATIO support in MESS. Try MatrixMarket/HB instead.\n");
            return MESS_ERROR_NOSUPPORT;
#endif

        }



    /*-----------------------------------------------------------------------------
     *  Try Matrix Market
     *-----------------------------------------------------------------------------*/
    f = csc_io_open ( filename, CSC_IO_FILE_READ);
    if ( !f ) {
        MSG_ERROR("error opening %s\n", filename);
        return MESS_ERROR_FILEIO;
    }
    ret = __mm_read_info(f, matrix);

    if ( ret == 0 ) {
        ret = __mm_read_data(f,matrix);
        MSG_INFO("Read MM file successfully\n");

        if ( ret > 0 ) {
            MSG_ERROR("Error %d - %s occurred\n", (int) ret, mess_get_error((int)ret));
            if (matrix->rowptr)mess_free(matrix->rowptr);
            if (matrix->colptr)mess_free(matrix->colptr);
            if (MESS_IS_REAL(matrix) && matrix->values )mess_free(matrix->values);
            if (MESS_IS_COMPLEX(matrix) && matrix->values_cpx)mess_free(matrix->values_cpx);
            matrix->rowptr= NULL;
            matrix->colptr= NULL;
            matrix->values= NULL;
            matrix->values_cpx = NULL;
        }
        csc_io_close(f);
        return 0;

    }
    /* Try HB   */
    csc_io_close(f);
    f = csc_io_open ( filename, CSC_IO_FILE_READ);
    if ( !f ) {
        MSG_ERROR("error opening %s\n", filename);
        return MESS_ERROR_FILEIO;
    }

    ret = __hb_read(f, matrix);
    csc_io_close(f);
    if ( ret == 0 ) {
        MSG_INFO("Read HB file successfully\n");
        return 0;
    }

    return MESS_ERROR_FILEIO;
}


/**
 * @brief Read a matrix from a file and convert it to the desired format.
 * @param[in]  filename  input file to read
 * @param[out] matrix   matrix to save the read data
 * @param[in]  format    input desired format
 * @return zero on success or a non-zero error value otherwise
 *
 *
 * The @ref mess_matrix_read_formated function reads a matrix from a file and converts it
 * to the desired format. \n
 * It works like \ref mess_matrix_read and \ref mess_matrix_convert together.
 *
 * @sa mess_matrix_convert
 * @sa mess_matrix_read
 *
 *
 */
int mess_matrix_read_formated ( const char *filename, mess_matrix matrix, mess_storage_t format )
{
    MSG_FNAME(__func__);
    int ret =0 ;
    mess_matrix internal;

    mess_check_nullpointer(filename);
    mess_check_nullpointer(matrix);

    if ( !(format==MESS_DENSE || format==MESS_COORD || format==MESS_CSC || format==MESS_CSR)) {
        MSG_ERROR("Unknown storage format: %s\n", mess_storage_t_str(format));
        return MESS_ERROR_STORAGETYPE;
    }

    MESS_MATRIX_RESET(matrix);

    ret = mess_matrix_init(&internal);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
    ret = mess_matrix_read(filename, internal);     FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_read);
    ret = mess_matrix_convert(internal, matrix,format); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_convert);
    mess_matrix_clear(&internal);

    return 0;
}       /* -----  end of function mess_matrix_read_formated  ----- */

