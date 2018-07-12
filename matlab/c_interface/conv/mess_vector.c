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
 * @file matlab/c_interface/conv/mess_vector.c
 * @brief
 * @author @mbehr
 */


#include "interface_matlab.h"

/**
 * @brief Copy a vector/matrix from @matlab  to a column vector in @mess.
 * @param[in] Amatlab input @matlab  object
 * @return @c NULL if an error occured or the vector  is empty or a @ref mess_vector with the content copied from the @matlab  object
 *
 * The @ref mess_vector_from_mexmess function copies a vector or a matrix from @matlab  to @mess. \n
 * Data is stored in a column-wise way if the input is a matrix.
 *
 */
mess_vector mess_vector_from_mexmess (const mxArray *Amatlab)
{
    mess_int_t i;
    mess_int_t nz;
    mess_vector A;
    int ret = 0;

    if ( Amatlab == NULL ) {
        csc_warn_message("Input mxArray points to NULL\n");
        return NULL;
    }
    if ( mxIsSparse(Amatlab)) {
        csc_warn_message("mess_vector_from_matlab is not available for sparse data\n");
        return NULL;
    }


    nz = mxGetM(Amatlab) * mxGetN(Amatlab);
    if(nz==0){
        return NULL;
    }

    if ( mxIsComplex (Amatlab)){
        ret = mess_vector_init(&A);
        ret = mess_vector_alloc(A, nz, MESS_COMPLEX);
        if ( ret ) {
            csc_warn_message("Failed to allocate output vector\n");
            return NULL;
        }
#if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble *temp = mxGetComplexDoubles(Amatlab);
        memcpy(A->values_cpx, temp, sizeof(mess_double_cpx_t)*nz);
#else
        double *v, *vi;
        v = mxGetPr(Amatlab);
        vi = mxGetPi(Amatlab);
        for ( i = 0; i < nz ; i++ ) {
            A->values_cpx[i] = v[i] + vi[i]*I;
        }
#endif
    } else {
        ret = mess_vector_init(&A);
        ret = mess_vector_alloc(A, nz, MESS_REAL);
        if ( ret ) {
            csc_warn_message("Failed to allocate output vector\n");
            return NULL;
        }
#if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble * temp = mxGetDoubles(Amatlab);
        memcpy(A->values,temp,sizeof(double)*nz);
#else
        double *v;
        v = mxGetPr(Amatlab);
        for ( i = 0; i < nz ; i++ ) {
            A->values[i] = v[i];
        }
#endif
    }
    return A;
}

/**
 * @brief Copy a vector from @mess to @matlab.
 * @param[in] A input MESS vector
 * @return the @matlab vector representaion of the input vector.
 *
 * The @ref mess_vector_to_mexmess function copies a vector from @mess to a column vector in @matlab.
 *
 */

mxArray * mess_vector_to_mexmess (mess_vector A)
{
    mxArray *Amatlab ;
    mess_int_t i;

    if ( A == NULL) {
        return mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    if ( MESS_IS_COMPLEX(A)) {

        Amatlab = mxCreateDoubleMatrix(A->dim, 1, mxCOMPLEX);
#if MX_HAS_INTERLEAVED_COMPLEX
        mxComplexDouble * temp = mxGetComplexDoubles(Amatlab);
        memcpy(temp,A->values_cpx,sizeof(mess_double_cpx_t)*(A->dim));
#else
        double *v, *vi;
        v = mxGetPr(Amatlab);
        vi = mxGetPi (Amatlab);
        for ( i = 0; i < A->dim ; i++ ) {
            v [ i ] = creal(A->values_cpx[i]);
            vi [ i ] = cimag (A->values_cpx[i]);
        }
#endif
    } else if ( MESS_IS_REAL(A)) {
        Amatlab = mxCreateDoubleMatrix(A->dim, 1, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *temp = mxGetDoubles(Amatlab);
        memcpy(temp, A->values, sizeof(double)*(A->dim));
#else
        double *v;
        v = mxGetPr(Amatlab);
        for ( i = 0; i < A->dim ; i++ ) {
            v [ i ] = A->values[i];
        }
#endif
    } else {
        return mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    return Amatlab;
}


