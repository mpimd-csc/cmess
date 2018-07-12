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
 * @file tutorials/eig/tutorial_eig.c
 * @brief Tutorial for demonstrating dense eigenvalue and svd decomposition using @ref mess_eigen_eig and @ref mess_eigen_svd.
 * @author @mbehr
 * @author @koehlerm
 *
 * # Tutorial Eigenvalue Computation and Singular Value Decomposition
 * In this tutorial we show how to solve the eigenvalue problem
 * \f[ Ax = \lambda x, \f]
 * where \f$ A\f$ is a matrix and \f$ \lambda \f$ is an eigenvalue with \f$ x \neq 0 \f$ as its corresponding
 * eigenvector.
 * Furthermore we show how to compute the singular value decomposition \f$ A= U \Sigma V^H\f$.
 * We also plot the eigenvalues as well as the singular values.
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/eig/tutorial_eig.c" HEADER
 * \n
 *
 *
 * ## 2. Check Input arguments, Init Matrix and Read Data
 * We check the number of input arguments by using the @c argc variable.
 * We further assume that both input argument are paths to the @mm file that is needed for
 * reading a matrix.
 * @snippet "tutorials/eig/tutorial_eig.c" INITREAD
 * \n
 *
 * ## 3. Compute the Eigenvalues and Singular Values
 * The function @ref mess_eigen_eig is able to compute the eigenvalues and
 *  @ref mess_eigen_svd computes the singular value decomposition of \f$A\f$.
 * We do not want to compute the eigenvectors and singularvectors, therefore
 * we supress the computation by inserting a @c NULL pointer to the corresponding argument.
 * @snippet "tutorials/eig/tutorial_eig.c" EIGSVD
 * \n
 *
 * ## 4. Plotting the Eigenvalues and Singular Values
 * Now we want to plot the eigenvalues and the singular values on the complex plane.
 * Therefore we use a @ref mess_plotExport instance.
 * The function @ref mess_plotExport_init set the name of the plot the scaling of
 * the axes, the labels for each axes and the number of datasets.
 * We label every dataset with help of @ref mess_plotExport_setLabel.
 * Since we want to plot discrete points we set the type to
 * @c "points". \n
 *
 *
 * ## 5. Add data to the plot
 * We add the singular as well as the eigenvalues to the plot.
 * In general the eigenvalues are complex. Therefore we have
 * to check the type of the @ref mess_vector which contains
 * the eigenvalues. Adding singular values works in a similar fashion and
 * it is easier because singular values are real. \n
 *
 * ## 6. Create a gnuplot script
 * In order to make the results viewable we export the plot to
 * @gnuplot with @ref mess_plotExport_createGnuScript. \n
 * After execution of the programm a file @c Eigenvalues_and_Singular_Values_script.gnu
 * is created. We can check our results by using @gnuplot.
 * @snippet "tutorials/eig/tutorial_eig.c" PLOT
 * \n
 *
 * ## 7. Clear everything
 * We can clear the plot by @ref mess_plotExport_clear.
 * Another usefull way to clear @ref mess_matrix and @ref mess_vector instances
 * is to use the macros @ref MESS_CLEAR_MATRICES and @ref MESS_CLEAR_VECTORS.
 * The advantage is that these macros are able to handle several instances in one call.
 * However using @ref mess_matrix_clear and @ref mess_vector_clear is also possible.
 * @snippet "tutorials/eig/tutorial_eig.c" CLEAR
 *
 * @sa @ref tutorials
 *
*/

//cond -> code is ignored by doxygen documentation, but snippet references are working


///@cond
///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include "mess/mess.h"
///[HEADER]

int main ( int argc, char **argv){

    /// [DECLARE]
    mess_matrix A;
    mess_vector evals,svals;
    /// [DECLARE]

    /// [INITREAD]
    if(argc!=2){
        printf("usage: %s A.mtx",argv[0]);
        return 1;
    }

    MESS_INIT_MATRICES(&A);
    mess_matrix_read_formated(argv[1],A,MESS_DENSE);
    if(A->rows!=A->cols){
        mess_matrix_printinfo(A);
        MESS_CLEAR_MATRICES(&A);
        printf("First matrix has wrong size.\n");
    }
    /// [INITREAD]

    /// [EIGSVD]
    mess_vector_init(&evals);
    mess_vector_alloc(evals,A->rows,MESS_COMPLEX);
    mess_eigen_eig(A,evals,NULL);
    mess_vector_init(&svals);
    mess_vector_alloc(svals,A->rows,MESS_REAL);
    mess_eigen_svd(A,svals,NULL,NULL);
    /// [EIGSVD]

    /// [PLOT]
    mess_plotExport plot;
    mess_plotExport_init(&plot, "Eigenvalues_and_Singular_Values", MESS_PLOT_LIN, MESS_PLOT_LIN, "Real", "Imag", 2);
    mess_plotExport_setLegendPos(plot, "belowCenter");
    mess_plotExport_setLabel(plot, 0, "Eigenvalues");
    mess_plotExport_setLabel(plot, 1, "Singular Values");
    mess_plotExport_setType(plot, 0, "points");
    mess_plotExport_setType(plot, 1, "points");
    mess_plotExport_setColor(plot, 0, "red");
    mess_plotExport_setColor(plot, 1, "blue");

    //add eigenvalues
    mess_int_t i;
    if(MESS_IS_REAL(evals)){
        for(i=0;i<evals->dim;++i){ mess_plotExport_addData(plot,0,evals->values[i],0.0); }
    }else{
        for(i=0;i<evals->dim;++i){ mess_plotExport_addData(plot,0,creal(evals->values_cpx[i]),cimag(evals->values_cpx[i])); }
    }

    //add singular values, singular values are in every case real
    for(i=0;i<svals->dim;++i){ mess_plotExport_addData(plot,1,svals->values[i],0.0); }

    // save plot
    mess_plotExport_createGnuScript(plot);
    /// [PLOT]

    /// [CLEAR]
    mess_plotExport_clear(&plot);
    MESS_CLEAR_MATRICES(&A);
    MESS_CLEAR_VECTORS(&evals,&svals);
    /// [CLEAR]

    return 0;

}
///@endcond





