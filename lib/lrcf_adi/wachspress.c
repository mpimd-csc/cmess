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
 * @file lib/lrcf_adi/wachspress.c
 * @brief Auxiliary routines for Wachspress parameter computation.
 * @author @koehlerm
 * @author @saak
 *
 */



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mess/mess.h"
#include "mess/error_macro.h"



static double min(double x,double y){/* gives the minimum of x and y*/
    if (x<y) return x;
    else return y;
}

static double dn(double u,double k,double TOL){
    double phi1,phi0;
    double a[100],b[100],c[100];
    int i,j;
    a[0]=1.0;
    c[0]=k;
    b[0]=min(1-mess_eps(),sqrt(1-k*k));
    i=0;
    while (fabs(c[i])>TOL){
        i++;
        a[i]=(a[i-1]+b[i-1])/2;
        b[i]=sqrt(a[i-1]*b[i-1]);
        c[i]=(a[i-1]-b[i-1])/2;
    }
    phi1=pow(2.0,i)*a[i]*u;
    phi0=0.0;
    j=i;
    while (j>0){
        if (j<i) phi1=phi0;
        phi0=(phi1+asin(c[j]*sin(phi1)/a[j]))/2;
        j--;
    }
    if ((1.0-k)<1e-5)
        return cos(phi0)/cos(phi1-phi0);
    else
        return sqrt(1.0-k*k*sin(phi0)*sin(phi0));
}

static double ellipK(double k,double TOL){
    double a0,a1,b0,b1,c1;
    int i1;

    a1=a0=1;
    b0=sqrt(1.0-k*k);
    b0=min(1-mess_eps(),b0);
    i1=0;
    c1=1.0;
    while(c1>TOL){
        a1=(a0+b0)/2;
        b1=sqrt(a0*b0);
        c1=(a0-b0)/2;
        i1++;
        a0=a1;
        b0=b1;
    }
    return M_PI/(2*a1);
}

static double F(double phi,double k,double TOL){
    int i;
    double a0,a1,b0,b1,c1,phi_temp;

    a1=a0=1;
    b0=min(1-mess_eps(),sqrt(1.0-k*k));
    i=0;
    c1=1.0;
    while(c1>TOL){
        a1=(a0+b0)/2;
        b1=sqrt(a0*b0);
        c1=(a0-b0)/2;
        phi_temp=phi+atan((b0/a0)*tan(phi));
        phi = phi_temp + M_PI*floor(phi_temp/M_PI+.5);
        i++;
        a0=a1;
        b0=b1;
    }
    return phi/(pow(2,i)*a1);
}

/**
 * @brief Compute Wachspress parameters.
 * @param[in] a  input lower spectral bound
 * @param[in] b  input upper spectral bound
 * @param[in] alpha  input sector angle
 * @param[in] tol  input tolerance
 * @param[out] p   vector containing Wachspress parameters
 * @param[out] lp length of vector \f$ p \f$
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_lrcfadi_para_wachspress function computes optimal ADI-shift parameters according to
 * the formulas given in \cite Wac95 for matrices with real spectrum. \n
 * The spectral bounds \f$ a ,b \f$ and the sector angle \f$ alpha \f$ are given as
 * \f[\begin{array}{cc} \displaystyle
 *     a &= \min\limits_{i=1, \cdots , n} \Big( \operatorname{Re}(\lambda_i) \Big) , \\
 *     b &= \max\limits_{i=1, \cdots , n} \Big( \operatorname{Re}(\lambda_i) \Big) , \\
 *    \alpha &= \tan^{-1} \max\limits_{i=1, \cdots , n} \Bigg\vert
 *            \displaystyle \frac{\operatorname{Im}(\lambda_i)}{\operatorname{Re}(\lambda_i)} \Bigg\vert,
 * \end{array}
 * \f]
 * where \f$ \lambda_1 , \cdots, \lambda_n \f$ are eigenvalues of the negative matrix. \n
 * Elliptic functions and integrals are computed to a tolerance \f$ tol \f$.
 *
 * More details can be found in \cite Wac95.
 *
 */
int mess_lrcfadi_para_wachspress(double a, double b, double alpha, double tol, double *p, mess_int_t *lp){
    MSG_FNAME(__func__);
    double k,kstrich,cK,phi,v,c2,m;
    const double localtol=1e-15;
    mess_int_t i,J=0;
    int computeJ = 0;

    if ( tol < 1.0) {
        computeJ = 1;
    }

    if (alpha==0)
        kstrich = a/b;
    else {
        c2=2.0/(1.0+(a/b+b/a)/2);
        m = 2.0*cos(alpha)*cos(alpha)/c2 - 1.0;
        if (m<1.0) {
            MSG_ERROR("Parameters would be complex.\n Please check your input data!\n");
            return MESS_ERROR_ARGUMENTS;
        }
        kstrich = 1.0 / (m+sqrt(m*m-1));
    }
    k = sqrt(1.0-kstrich*kstrich);
    k = min(1-mess_eps(),k);
    cK = ellipK(k,localtol);
    if (alpha==0)
        v = ellipK(kstrich,localtol);
    else {
        phi = asin(sqrt(a/(b*kstrich)));
        v = F(phi,kstrich,localtol);
    }

    if ( computeJ ){
        J = ceil(cK/(2*v*M_PI)*log(4.0/tol));
    } else {
        J = (mess_int_t) tol;
    }
    if ( J > 50 ) {
        MSG_WARN("set J from " MESS_PRINTF_INT " to 50\n", J);
        J = 50;
    }

    MSG_INFO(" J: " MESS_PRINTF_INT "\n v: %1.30e\n k: %1.30e\n kstrich: %1.30e\n cK: %1.30e\n\n",J,v,k,kstrich,cK);
    for(i=0;i<J;i++){
        p[i]= -sqrt(a*b/kstrich)*dn(((i+1)-0.5)*cK/J,k,localtol);
        // MSG_INFO("p(" MESS_PRINTF_INT ") = %1.8e\n",i,p[i]);
    }
    *lp = J;
    return 0;
}

