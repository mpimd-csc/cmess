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
 * @file matlab/c_interface/interface_matlab.c
 * @brief
 * @author @mbehr
 *
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/error_message.h"
#include "interface_matlab.h"
#include <complex.h>

void my_mexPrintf(const char * str) {
    mexPrintf("%s", str);
    mexEvalString("drawnow;");
}

void my_mexError(const char * str) {
    mexErrMsgIdAndTxt("MEX:MESS:internalerror", str);
    mexEvalString("drawnow;");
}

void my_mexWarn(const char * str) {
    mexWarnMsgTxt(str);
    mexEvalString("drawnow;");
}


extern csc_error_print_t error_handle;
extern csc_error_print_t warn_handle;
extern csc_error_print_t info_handle;
extern csc_error_print_t print_handle;

#ifdef MATLAB_LINUX_32
void error_message2(const char *fmt, ...) {
    char str[2048];
    va_list ap;

    /* Build the @ref message  */
    va_start(ap, fmt);
    /*  str = make_message(fmt, ap); */
    vsnprintf(str, 2048, fmt, ap);
    va_end(ap);

    /* Print the @ref message  */
    error_handle(str);
    /*  /mexErrMsgTxt(str); */

    /*  error_free(str); */
}
#endif

static void * matlab_malloc(size_t len)
{
    return mxMalloc((mwSize) len);
}
static void * matlab_realloc(void *ptr, size_t len)
{
    return mxRealloc(ptr, len);
}
static void * matlab_calloc(size_t nmeb, size_t smeb)
{
    return mxCalloc(nmeb,smeb);
}
static void matlab_free(void *ptr)
{
    return mxFree(ptr);
}

/**
 * @brief Initialize the @mess to be used in a MeX function.
 * @return zero on success or a non zero error code otherwise
 *
 * The @ref init_mexmess function initializes the mex helper functions. It needs to called once at the begin of a mex file.
 */
void init_mexmess() {
    /* Set the  Message handlers */
    csc_error_message_handle(my_mexError);
    csc_warn_message_handle(my_mexWarn);
    csc_info_message_handle(my_mexPrintf);
    csc_print_message_handle(my_mexPrintf);
    csc_error_message_memory((void *) mxMalloc, (void*) mxRealloc, (void *) mxFree);
#ifdef MATLAB_LINUX_32
    error_handle = error_message2;
#endif
    /* Set the output level to error_only  */
    mess_error_level = 1;

    /* Set the Malloc handler  */
    mess_set_malloc( (void *) matlab_malloc, (void*) matlab_realloc, (void *) matlab_calloc, (void*) matlab_free);

}


/**
 * @brief Get a real valued scalar from @matlab.
 * @param[in] Amatlab input pointer to a mxArray
 * @return first real double value from \f$ Amatlab \f$
 *
 * The @ref mess_double_from_mexmess function get a real double value from @matlab .
 *
 */

double mess_double_from_mexmess (const mxArray *Amatlab)
{
    double ret = 0;
    MSG_FNAME(__func__);
    if ( Amatlab == NULL){
        MSG_ERROR("Amatlab points to NULL\n");
        return 0;
    }
#if MX_HAS_INTERLEAVED_COMPLEX
    ret = (*mxGetDoubles(Amatlab));
#else
    ret = (*mxGetPr(Amatlab));
#endif
    return ret;
}


/**
 * @brief Get a complex valued scalar from @matlab.
 * @param[in] Amatlab input pointer to a mxArray
 *@return first complex double value from \f$ Amatlab \f$
 *
 * The @ref mess_complex_from_mexmess function get a complex double value from @matlab .
 *
 */

mess_double_cpx_t mess_complex_from_mexmess (const mxArray *Amatlab)
{
    mess_double_cpx_t ret = 0;
    MSG_FNAME(__func__);
    if ( Amatlab == NULL){
        MSG_ERROR("Amatlab points to NULL\n");
        return 0;
    }

#if MX_HAS_INTERLEAVED_COMPLEX
    if (!mxIsComplex(Amatlab)) {
        ret = (*mxGetDoubles(Amatlab));
    } else {
        ret = mxGetComplexDoubles(Amatlab)->real +  mxGetComplexDoubles(Amatlab)->imag*I;;
    }
#else
    if (!mxIsComplex(Amatlab)) {
        ret = (*mxGetPr(Amatlab));
    } else {
        ret = (*mxGetPr(Amatlab)) + (*mxGetPi(Amatlab))*I;
    }
#endif
    return ret;
}



int mess_int_array_to_matlab (mess_int_t * p, mess_int_t len , mxArray **Amat)
{
    MSG_FNAME(__func__);
    mxArray *Amatlab ;
    mess_int_t i;

    if ( p == NULL) {
        MSG_ERROR("A points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( Amat == NULL ){
        MSG_ERROR("Amat points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }

    Amatlab = mxCreateDoubleMatrix(len, 1, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *v;
    v = mxGetDoubles(Amatlab);
    for ( i = 0; i < len ; i++ ) {
        v [ i ] = (double) p [i];
    }
#else
    double *v;
    v = mxGetPr(Amatlab);
    for ( i = 0; i < len ; i++ ) {
        v [ i ] = (double) p [i];
    }
#endif
    *Amat = Amatlab;
    return 0;
}


mxArray * mxCreateComplexScalar(mess_double_cpx_t s )
{
    mxArray * retv = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    if ( retv == NULL ) return NULL;

#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble * temp = mxGetComplexDoubles(retv);
    temp->real = creal(s);
    temp->imag = cimag(s);
#else
    double *rp, *ip ;
    rp = mxGetPr(retv);
    ip = mxGetPi(retv);
    *rp = creal ( s );
    *ip = cimag ( s );
#endif
    return retv;
}


