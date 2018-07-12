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
 * @file lib/matrix/shellsort.c
 * @brief Compact shellsort implementation for matrix data structures.
 * @author @koehlerm
 */

/**
 *  @brief Shellsort algorithm for matrix structures.
 *  @param[in,out] ind index array
 *  @param[in,out] values real value array
 *  @param[in,out] values_cpx complex value array
 *  @param[in] l  input first position to sort
 *  @param[in] r  input last postion to sort
 *  @param[in] type  input selects the data type
 *
 *  Perform a shell soft of parts if the ind and values/values_cpx array.
 *  This file does not provide a user callable routine. It is only included
 *  by other files which needs this algorithm and where using the qsort from
 *  the standard libc is too expensive because of the restructuring of the input
 *  data.
 *
 */
static void  __attribute__((unused))__shellsort( mess_int_t *ind,
        double *values,
        mess_double_cpx_t *values_cpx,
        mess_int_t l,
        mess_int_t r,
        unsigned short type)
{
    if ( r <= l ) return ;
    mess_int_t h;
    mess_int_t i, j;
    mess_int_t vi;
    double vv=0;
    mess_double_cpx_t vim=0;
    for (h=1; h <= (r-l)/9; h=3*h+1);
    for ( ; h>0; h=h/3){
        for (i = l+h; i<=r; i++){
            j = i;
            vi=ind[i];
            if ( type == MESS_REAL) vv = values[i];
            if ( type == MESS_COMPLEX) vim = values_cpx[i];

            while (j >= l+h && vi < ind[j-h]){
                ind[j] = ind [j-h];
                if ( type == MESS_REAL) values[j] = values[j-h];
                if ( type == MESS_COMPLEX) values_cpx[j] = values_cpx[j-h];
                j=j-h;
            }
            ind[j] = vi;
            if ( type == MESS_REAL) values[j] =vv;
            if ( type == MESS_COMPLEX) values_cpx[i] = vim;
        }

    }
}

