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
 *
 * @file tutorials/bt/tutorial_bt.c
 * @author mbehr
 * @author koehlerm
 *
 * # Tutorial: Balanced Truncation  Model Order Reduction of LTI System
 * In this Tutorial we show how @mess is used for Model Order Reduction using Balanced Truncation.
 *
 * Therefore, we consider a large scale generalized LTI system:
 * \f[
 *  \begin{array}{ccccc}
 *    E\dot{x}(t)&=& A x(t) &+ B u(t) & , \\
 *    y(t) &=& C x(t) & & ,
 *  \end{array}
 * \f]
 * where \f$ E,\,A \in \mathbb{R}^{n \times n},\; B \in \mathbb{R}^{n\times p} \f$ and \f$C \in \mathbb{R}^{q\times n}\f$.
 * Furthermore, the matrix \f$ E \f$ is nonsingular and the pencil  \f$ (A,\,E)\f$ is stable.  For this LTI system we compute
 * factors of the Controllability Gramian \f$ Z_CZ_C^T \approx X_C \f$ and the Observability Gramian \f$ Z_OZ_O^T\approx X_O \f$,
 * which are the solutions of the generalized Lyapunov Equations
 * \f[
 *  \begin{array}{ccccc}
 *    AX_C E^T + & EX_C A^T &=& -BB^T &, \\
 *    A^T \hat{X_O} E + &E^T \hat{X_O} A &=& -CC^T &, \\
 *    & X_O &=& E^T \hat{X_O} E &.
 *  \end{array}
 *  \f]
 * Using this factors we compute the SVD \cite morMoo81  of
 * \f[
 *    U_O \Sigma U_C^T = Z_O^TZ_C
 * \f]
 * to get the balancing tranformation matrices
 * \f[
 *   S_C = Z_CU_C(:,1:k)\Sigma(1:k,1:k)^{-\frac{1}{2}}\quad\text{and}\quad S_O = Z_OU_O(:,1:k)\Sigma(1:k,1:k)^{-\frac{1}{2}},
 * \f]
 * where \f$ k \f$ is selected such that the first \f$ k \f$ singluar values in \f$ \Sigma \f$ fulfill a given threshold.
 * Using the transformation matrices \f$ S_C \f$ and \f$ S_O \f$ we compute the reduced order model
 * \f[
 *  E := I,\; A_r := S_O^TE^{-1}AS_C, \; B_r := S_O^TB, \;\text{and}\; C_r = C S_C.
 * \f]
 * The threshold for the selection of \f$k\f$ using the following approximation for the \f$ \mathcal{H}\infty \f$-error:
 * \f[ \Vert H-H_r \Vert_{\infty} < 2 \sum\limits_{i= k+1}^{n} \sigma_i < \tau,  \f]
 * where \f$\tau\f$ is a given upper bound for the error and \f$\sigma_i\f$ are the singular values contained in \f$\Sigma\f$.
 *
 *
 * ## 1. Include Header Files
 * We include standard @c C header files and in addtion the @c mess.h header file
 * to make the @mess library available.
 * @snippet "tutorials/bt/tutorial_bt.c" HEADER
 * \n
 *
 * ## 2. Declare Matices, Dynamical Systems and Options Structures
 * We use the \ref mess_matrix structure to declare our matrices \f$ A, B, C, E\f$.
 * Furthermore, we need the matrices \f$S_C\f$ and \f$ S_O \f$ containing the left and the right projection
 * matrix of the Balanced Truncation process. The balanced truncation process needs options for solving the Lyapunov
 * equation and for the truncation process. These are represendted by \ref mess_options, \ref mess_bt_options
 * and \ref mess_bt_status. Because the four matrices \f$ A, B, C, E \f$ represent an LTI system we need two intances
 * of \ref mess_dynsys for the original and the reduced order model.
 * @snippet "tutorials/bt/tutorial_bt.c" DECLARE
 * By using the \c argc variable we can check the number of input arguments and set the error level to 1 to removing
 * debug and warnings from the output.
 * @snippet "tutorials/bt/tutorial_bt.c" INPUT
 * \n
 *
 * ## 3. Init Matrices and Read Data
 *
 * Before we are able to use the matrices we have to init the fields of the \ref mess_matrix structure variables
 * using \ref mess_matrix_init and read the matrices from @mm file format. The remaining command line parameters
 * are used for the tolerance \f$\tau\f$ and the maximum order of the reduced order model.
 * @snippet "tutorials/bt/tutorial_bt.c" READ
 * \n
*
* ## 4. Converting the Input Matrices into a (Generalized) LTI System
* Because the matrices represent a generalized LTI system we have to transfer them into a \ref mess_dynsys object.
* This either be done using a flat copy by \ref mess_dynsys_glti or as a deep copy by \ref mess_dynsys_glti_copy. Deep copy
* means that the matrices a copied into the object and not only referenced. Otherwise if one changes the matrix A it will
* be change in the \ref mess_dynsys object as well.
* @snippet "tutorials/bt/tutorial_bt.c" DYNSYS
* \n
*
* ## 5. Setup the Balanced Truncation and Compute \f$S_C\f$ and \f$S_O\f$
* We set the options to the Balanced Truncation process and compute the transformation matrices. The \ref mess_bt_gslrsrm
* computes transformation matrices such that the input is a generalized LTI system and the reduced order model will be
* a standard state space system. The function internally solves the two Lyapunov equations described above with a low-rank
* solver after wards it computes the SVD of the product of both low rank factors and selects the order \f$k\f$ of the
* reduced order model by the singular values of the SVD as described above. Furthermore, a maximum oder is defined in the
* Balanced Truncation options. If the matrix \f$ E \f$ is the identity the \ref mess_bt_lrsrm function can be used.
* If the reduced order model should have a matrix which is not the identity one has to use \ref mess_bt_gglrsrm.
* @snippet "tutorials/bt/tutorial_bt.c" BT
* \n
*
* ## 6. Compute the Reduced Order Model using \f$S_C\f$ and \f$S_O\f$
* After computing the transformation matrices \f$S_C\f$ and \f$S_O\f$ the reduced order model
* \f[
    *  \begin{array}{ccccc}
    *    \dot{x_r}(t)&=& A_r x_r(t) &+ B_r u(t) & , \\
        *    y_r(t) &=& C_r x_r(t) & & ,
    *  \end{array}
    * \f]
    * is computed by
    * \f[
        *  A_r := S_O^TE^{-1}AS_C, \; B_r := S_O^TB, \;\text{and}\; C_r = C S_C.
        * \f]
        * @snippet "tutorials/bt/tutorial_bt.c" REDUCE
        * \n
        *
        * ## 7. Compute the \f$\mathcal{H}_2\f$
        * Although the \f$ \mathcal{H}_2\f$ norm does not correspond to the Balanced Truncation process it is a good identifier
        * how good the reduce order model approximates the original model. Computing the distance or the error between two LTI
        * systems is done using the \ref mess_h2_error function.
        * @snippet "tutorials/bt/tutorial_bt.c" H2ERROR
        * \n
        *
        * ## 8. Evaluate the Transferfunction of the Original and the Reduced Order Model
        * Because the transferfunction is a nice tool to visualize the input-output behaviour of an LTI in the frequency domain,
        * we want to evaluate it for the original and the reduced order model in order to check if the computed model still shows
        * the same input-output behaviour. Therefore, we sample the transferfunction of both models over a interval
        *  \f$\omega \in \left[ 10^{-1}, 10^6 \right]\f$ and compute \f$ max(svd(H(i\cdot \omega)))\f$. The whole procedure is done by
        *  \ref mess_dynsys_evaltransfer.
        * @snippet "tutorials/bt/tutorial_bt.c" TRANSFER
        * \n
        *
        * ## 9. Clear memory
        * Print the final status and clear all matrices and data structures.
        * @snippet "tutorials/bt/tutorial_bt.c" STATUS_AND_CLEAR
        * \n
        *
        *
        * @sa @ref tutorials
        */


        ///@cond
        ///[HEADER]
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "mess/mess.h"
        ///[HEADER]

        int main ( int argc, char **argv){

            ///[DECLARE]
            mess_matrix A, B, C, E,SO, SC;
            mess_options adiopt;
            mess_bt_options btopt;
            mess_bt_status  btstat;
            mess_dynsys lti, red;
            double tol;
            mess_int_t maxr;
            ///[DECLARE]

            ///[INPUT]
            if ( argc != 7 ) {
                printf("usage: %s E.mtx A.mtx B.mtx C.mtx tol maxr\n", argv[0]);
                return 0;
            }
            mess_set_errorlevel(1);
            ///[INPUT]

            ///[READ]
            mess_matrix_init(&A);
            mess_matrix_init(&B);
            mess_matrix_init(&C);
            mess_matrix_init(&E);

            mess_matrix_read_formated(argv[1],E,MESS_CSR);
            mess_matrix_read_formated(argv[2],A,MESS_CSR);
            mess_matrix_read_formated(argv[3],B,MESS_DENSE);
            mess_matrix_read_formated(argv[4],C,MESS_DENSE);

            tol  = atof(argv[5]);
            maxr = atoi(argv[6]);
            ///[READ]

            ///[DYNSYS]
            mess_dynsys_init(&lti);
            mess_dynsys_glti_copy(lti, A, B, C, E);

            printf("Dimension: " MESS_PRINTF_INT "\n", A->rows);
            printf("Inputs:    " MESS_PRINTF_INT "\n", B->cols);
            printf("Outputs:   " MESS_PRINTF_INT "\n", C->rows);
            ///[DYNSYS]

            ///[BT]
            mess_bt_status_init(&btstat);
            mess_bt_options_init(&btopt);
            btopt->rdim =  maxr;
            btopt->tol  = tol;

            printf("\nBalanced Truncation Options:\n");
            mess_bt_options_print(btopt);
            printf("\n");

            mess_options_init(&adiopt);
            adiopt->adi_output = 0;
            adiopt->adi_shifts_paratype = MESS_LRCFADI_PARA_MINMAX;

            mess_matrix_init(&SO);
            mess_matrix_init(&SC);

            mess_bt_gslrsrm(lti,btopt, adiopt, SO,SC, btstat);

            printf("Balanced Truncation Status:\n");
            mess_bt_status_print(btstat);
            printf("\n");
            ///[BT]

            ///[REDUCE]
            double time_pro = mess_wtime();
            mess_dynsys_init(&red);
            mess_dynsys_project_to_lti(lti, SO,SC,red);
            time_pro = mess_wtime()-time_pro;
            ///[REDUCE]

            ///[H2ERROR]
            double err = 0;
            mess_h2_error(lti, red, &err);
            printf("H2 error = %lg\n", err);
            ///[H2ERROR]

            ///[TRANSFER]
            double time_tf = mess_wtime();
            mess_vector omega, G1, G2;
            MESS_INIT_VECTORS(&omega,&G1,&G2);
            mess_vector_alloc(omega,100,MESS_REAL);
            mess_vector_alloc(G1, 100, MESS_REAL);
            mess_vector_alloc(G2, 100, MESS_REAL);

            mess_dynsys_evaltransfer(lti,-1,6,100,omega,G1,NULL);
            mess_dynsys_evaltransfer(red,-1,6,100,omega,G2,NULL);

            mess_vector_diffnorm(G1,G2,&err);
            printf("error discrete ||G1-G2||_2 = %lg\n", err);
            mess_vector_diffnorminf(G1,G2,&err);
            printf("error discrete ||G1-G2||_infty = %lg\n", err);
            time_tf = mess_wtime()-time_tf;
            ///[TRANSFER]

            ///[STATUS_AND_CLEAR]
            printf("\nTime BT: %lg\nTime Projection: %lg\nTime Transfer Function: %lg\n", btstat->time, time_pro, time_tf);

            mess_vector_clear(&omega);
            mess_vector_clear(&G1);
            mess_vector_clear(&G2);

            mess_matrix_clear(&A);
            mess_matrix_clear(&B);
            mess_matrix_clear(&C);
            mess_matrix_clear(&E);
            mess_matrix_clear(&SO);
            mess_matrix_clear(&SC);

            mess_dynsys_clear(&lti);
            mess_dynsys_clear(&red);

            mess_options_clear(&adiopt);
            mess_bt_options_clear(&btopt);
            mess_bt_status_clear(&btstat);
            return 0;
            ///[STATUS_AND_CLEAR]
        }
///@endcond
