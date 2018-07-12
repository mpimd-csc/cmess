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
 * @file lib/lrcf_adi/equation_apply.c
 * @brief Wrapper for @ref mess_equation  function pointers.
 * @author @mbehr
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Initialization function of the linear mapping \f$x=\mathcal{A}y\f$.
 * @param[in] eqn input  Equation object \f$\mathcal{A}\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_A_pre function initializes the \f$x=\mathcal{A}y\f$ operation of a given
 * equation object. It is called normally once before the operator is uses the first time.
 * If the equation object does not provide an initialization functions the @ref mess_equation_A_pre
 * function returns immediately.
 *
 * @sa mess_equation_A_apply
 * @sa mess_equation_A_post
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_A_pre ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(eqn);

    if ( eqn->AX.generate == NULL) {
        return 0;
    }

    if ( eqn->AX.to_clear > 0 ) {
        MSG_INFO("The equation object's x = Ay function is already initialized.\n");
    }

    ret = eqn->AX.generate(eqn);
    eqn->AX.to_clear = 1;
    return ret;
}       /* -----  end of function mess_equation_A_pre  ----- */


/**
 * @brief Finalizes/Clears the linear mapping \f$x=\mathcal{A}y\f$.
 * @param[in]   eqn input     Equation object \f$x=\mathcal{A}y\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_A_post function finalizes the \f$x=\mathcal{A}y\f$ operator. That
 * means it calls the clear function of the equation object if it exits otherwise it will
 * return immediately.
 *
 * @sa mess_equation_A_pre
 * @sa mess_equation_A_apply
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_A_post ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(eqn);

    if ( eqn->AX.to_clear == 0 ) {
        return 0;
    }

    ret = eqn->AX.clear(eqn);
    eqn->AX.to_clear = 0;
    return ret ;
}       /* -----  end of function mess_equation_A_post  ----- */



/**
 * @brief Compute \f$x = \mathcal{A}y\f$.
 * @param[in] eqn input Equation object \f$x = \mathcal{A}y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input matrix y
 * @param[in] out input Ouput matrix x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_A_apply function computes \f$ x = op(\mathcal{A}) x\f$, where
 * x and y are matrices. For vectors use \ref mess_equation_A_apply_vector instead. If
 * the object does not provide a AX.apply function, \f$\mathcal{A} = I \f$ is assumed.
 *
 * @sa mess_equation_A_apply_vector
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_A_apply ( mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->AX.apply == NULL ) {
        /* Assume Identity */
        ret = mess_matrix_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        return 0;
    } else {
        ret = eqn->AX.apply(eqn, op, in, out);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_A_apply  ----- */

/**
 * @brief Compute \f$x = \mathcal{A}y\f$ vector valued.
 * @param[in] eqn input Equation object \f$x = \mathcal{A}y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input vector y
 * @param[in] out input Ouput vector x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_A_apply_vector function computes \f$ x = op(\mathcal{A}) x\f$, where
 * x and y are vectors. For matrices use \ref mess_equation_A_apply. If
 * the object does not provide a AX.apply function, \f$\mathcal{A} = I \f$ is assumed.
 *
 * @sa mess_equation_A_apply
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_A_apply_vector ( mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->AX.apply == NULL ) {
        /* Assume Identity */
        ret = mess_vector_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        return 0;
    } else {
        mess_matrix input, output;
        ret = mess_matrix_init(&input);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&output);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

        ret = mess_vector_tomatrix(in, input);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tomatrix);
        ret = eqn->AX.apply(eqn, op, input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
        ret = mess_vector_frommatrix(output, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_frommatrix);
        mess_matrix_clear(&input);
        mess_matrix_clear(&output);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_A_apply  ----- */

/**
 * @brief Initialization function of the linear mapping \f$x=\mathcal{E}y\f$.
 * @param[in] eqn input  Equation object \f$\mathcal{E}\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_E_pre function initializes the \f$x=\mathcal{E}y\f$ operation of a given
 * equation object. It is called normally once before the operator is uses the first time.
 * If the equation object does not provide an initialization functions the @ref mess_equation_E_pre
 * function returns immediately.
 *
 * @sa mess_equation_E_apply
 * @sa mess_equation_E_post
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_E_pre ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(eqn);

    if ( eqn->EX.generate == NULL) {
        return 0;
    }

    if ( eqn->EX.to_clear > 0 ) {
        MSG_INFO("The equation object's x = Ey function is already initialized.\n");
    }

    ret = eqn->EX.generate(eqn);
    eqn->EX.to_clear = 1;
    return ret;
}       /* -----  end of function mess_equation_E_pre  ----- */


/**
 * @brief Finalizes/Clears the linear mapping \f$x=\mathcal{E}y\f$.
 * @param[in]   eqn input     Equation object \f$x=\mathcal{E}y\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_E_post function finalizes the \f$x=\mathcal{E}y\f$ operator. That
 * means it calls the clear function of the equation object if it exits otherwise it will
 * return immediately.
 *
 * @sa mess_equation_E_pre
 * @sa mess_equation_E_apply
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_E_post ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(eqn);

    if ( eqn->EX.to_clear == 0 ) {
        return 0;
    }

    ret = eqn->EX.clear(eqn);
    eqn->EX.to_clear = 0;
    return ret ;
}       /* -----  end of function mess_equation_E_post  ----- */



/**
 * @brief Compute \f$x = \mathcal{E}y\f$.
 * @param[in] eqn input Equation object \f$x = \mathcal{E}y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input matrix y
 * @param[in] out input Ouput matrix x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_E_apply function computes \f$ x = op(\mathcal{E}) x\f$, where
 * x and y are matrices. For vectors use \ref mess_equation_E_apply_vector instead. If
 * the object does not provide a EX.apply function, \f$\mathcal{E} = I \f$ is assumed.
 *
 * @sa mess_equation_E_apply_vector
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_E_apply ( mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->EX.apply == NULL ) {
        /* Essume Identity */
#ifdef DEBUG
        MSG_WARN("EX apply is empty using E = I\n");
#endif
        ret = mess_matrix_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        return 0;
    } else {
        ret = eqn->EX.apply(eqn, op, in, out);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->EX.apply);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_E_apply  ----- */

/**
 * @brief Compute \f$x = \mathcal{E}y\f$ vector valued.
 * @param[in] eqn input Equation object \f$x = \mathcal{E}y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input vector y
 * @param[in] out input Ouput vector x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_E_apply_vector function computes \f$ x = op(\mathcal{E}) x\f$, where
 * x and y are vectors. For matrices use \ref mess_equation_E_apply. If
 * the object does not provide a EX.apply function, \f$\mathcal{E} = I \f$ is assumed.
 *
 * @sa mess_equation_E_apply
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_E_apply_vector ( mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->EX.apply == NULL ) {
        /* Essume Identity */
#ifdef DEBUG
        MSG_WARN("EX apply is empty using E = I\n");
#endif
        ret = mess_vector_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        return 0;
    } else {
        mess_matrix input, output;
        ret = mess_matrix_init(&input);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&output);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

        ret = mess_vector_tomatrix(in, input);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tomatrix);
        ret = eqn->EX.apply(eqn, op, input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->EX.apply);
        ret = mess_vector_frommatrix(output, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_frommatrix);
        mess_matrix_clear(&input);
        mess_matrix_clear(&output);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_E_apply  ----- */

/**
 * @brief Initialization function of the linear mapping \f$\mathcal{A}x=y\f$.
 * @param[in] eqn input  Equation object \f$\mathcal{A}\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_As_pre function initializes the \f$\mathcal{A}x=y\f$ operation of a given
 * equation object. It is called normally once before the operator is uses the first time.
 * If the equation object does not provide an initialization functions the @ref mess_equation_As_pre
 * function returns immediately.
 *
 * @sa mess_equation_As_apply
 * @sa mess_equation_As_post
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_As_pre ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(eqn);

    if ( eqn->AINV.generate == NULL) {
        return 0;
    }

    if ( eqn->AINV.to_clear > 0 ) {
        MSG_INFO("The equation object's  Ax=y function is already initialized.\n");
    }

    ret = eqn->AINV.generate(eqn);
    eqn->AINV.to_clear = 1;
    return ret;
}       /* -----  end of function mess_equation_As_pre  ----- */


/**
 * @brief Finalizes/Clears the linear mapping \f$\mathcal{A}x=y\f$.
 * @param[in]   eqn input     Equation object \f$\mathcal{A}x=y\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_As_post function finalizes the \f$\mathcal{A}x= y\f$ operator. That
 * means it calls the clear function of the equation object if it exits otherwise it will
 * return immediately.
 *
 * @sa mess_equation_As_pre
 * @sa mess_equation_As_apply
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_As_post ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(eqn);

    if ( eqn->AINV.to_clear == 0 ) {
        return 0;
    }

    ret = eqn->AINV.clear(eqn);
    eqn->AINV.to_clear = 0;
    return ret ;
}       /* -----  end of function mess_equation_As_post  ----- */



/**
 * @brief Compute \f$\mathcal{A}x=y\f$.
 * @param[in] eqn input Equation object \f$\mathcal{A}x=y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input matrix y
 * @param[in] out input Ouput matrix x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_As_apply function computes \f$op(\mathcal{A})x = y \f$, where
 * x and y are matrices. For vectors use \ref mess_equation_As_apply_vector instead. If
 * the object does not provide a AINV.apply function, \f$\mathcal{A} = I \f$ is assumed.
 *
 * @sa mess_equation_As_apply_vector
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_As_apply ( mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->AINV.apply == NULL ) {
        /* Assume Identity */
        ret = mess_matrix_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        return 0;
    } else {
        ret = eqn->AINV.apply(eqn, op, in, out);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AINV.apply);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_As_apply  ----- */

/**
 * @brief Compute \f$\mathcal{A}x = y\f$ vector valued.
 * @param[in] eqn input Equation object \f$\mathcal{A}x=y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input vector y
 * @param[in] out input Ouput vector x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_As_apply_vector function computes \f$ op(\mathcal{A})x = y \f$, where
 * x and y are vectors. For matrices use \ref mess_equation_As_apply. If
 * the object does not provide a AX.apply function, \f$\mathcal{A} = I \f$ is assumed.
 *
 * @sa mess_equation_As_apply
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_As_apply_vector ( mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->AINV.apply == NULL ) {
        /* Assume Identity */
        ret = mess_vector_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        return 0;
    } else {
        mess_matrix input, output;
        ret = mess_matrix_init(&input);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&output);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

        ret = mess_vector_tomatrix(in, input);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tomatrix);
        ret = eqn->AINV.apply(eqn, op, input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
        ret = mess_vector_frommatrix(output, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_frommatrix);
        mess_matrix_clear(&input);
        mess_matrix_clear(&output);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_As_apply  ----- */


/**
 * @brief Initialization function of the linear mapping \f$\mathcal{E}x=y\f$.
 * @param[in] eqn input  Equation object \f$\mathcal{E}\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_Es_pre function initializes the \f$\mathcal{E}x=y\f$ operation of a given
 * equation object. It is called normally once before the operator is uses the first time.
 * If the equation object does not provide an initialization functions the @ref mess_equation_Es_pre
 * function returns immediately.
 *
 * @sa mess_equation_Es_apply
 * @sa mess_equation_Es_post
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_Es_pre ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(eqn);

    if ( eqn->EINV.generate == NULL) {
        return 0;
    }

    if ( eqn->EINV.to_clear > 0 ) {
        MSG_INFO("The equation object's  Ax=y function is already initialized.\n");
    }

    ret = eqn->EINV.generate(eqn);
    eqn->EINV.to_clear = 1;
    return ret;
}       /* -----  end of function mess_equation_Es_pre  ----- */


/**
 * @brief Finalizes/Clears the linear mapping \f$\mathcal{E}x=y\f$.
 * @param[in]   eqn input     Equation object \f$\mathcal{E}x=y\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_Es_post function finalizes the \f$\mathcal{E}x= y\f$ operator. That
 * means it calls the clear function of the equation object if it exits otherwise it will
 * return immediately.
 *
 * @sa mess_equation_Es_pre
 * @sa mess_equation_Es_apply
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_Es_post ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(eqn);

    if ( eqn->EINV.to_clear == 0 ) {
        return 0;
    }

    ret = eqn->EINV.clear(eqn);
    eqn->EINV.to_clear = 0;
    return ret ;
}       /* -----  end of function mess_equation_Es_post  ----- */



/**
 * @brief Compute \f$\mathcal{E}x=y\f$.
 * @param[in] eqn input Equation object \f$\mathcal{E}x=y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input matrix y
 * @param[in] out input Ouput matrix x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_Es_apply function computes \f$op(\mathcal{E})x = y \f$, where
 * x and y are matrices. For vectors use \ref mess_equation_Es_apply_vector instead. If
 * the object does not provide a EINV.apply function, \f$\mathcal{E} = I \f$ is assumed.
 *
 * @sa mess_equation_Es_apply_vector
 * @sa mess_equation_As_apply
 *
 */
int mess_equation_Es_apply ( mess_equation eqn, mess_operation_t op, mess_matrix in, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->EINV.apply == NULL ) {
        /* Essume Identity */
#ifdef DEBUG
        MSG_WARN("EINV apply is empty using E = I\n");
#endif


        ret = mess_matrix_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_copy);
        return 0;
    } else {
        ret = eqn->EINV.apply(eqn, op, in, out);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->EINV.apply);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_Es_apply  ----- */

/**
 * @brief Compute \f$\mathcal{Ex} = y\f$ vector valued.
 * @param[in] eqn input Equation object \f$\mathcal{E}x=y\f$ blongs to
 * @param[in] op input  Opertion to apply
 * @param[in] in input  Input vector y
 * @param[in] out input Ouput vector x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_Es_apply_vector function computes \f$ op(\mathcal{E})x = y \f$, where
 * x and y are vectors. For matrices use \ref mess_equation_Es_apply. If
 * the object does not provide a AX.apply function, \f$\mathcal{E} = I \f$ is assumed.
 *
 * @sa mess_equation_As_apply
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_Es_apply_vector ( mess_equation eqn, mess_operation_t op, mess_vector in, mess_vector out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->EINV.apply == NULL ) {
        /* Essume Identity */
#ifdef DEBUG
        MSG_WARN("EINV apply is empty using E = I\n");
#endif

        ret = mess_vector_copy(in, out);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_copy);
        return 0;
    } else {
        mess_matrix input, output;
        ret = mess_matrix_init(&input);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&output);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);

        ret = mess_vector_tomatrix(in, input);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tomatrix);
        ret = eqn->EINV.apply(eqn, op, input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
        ret = mess_vector_frommatrix(output, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_frommatrix);
        mess_matrix_clear(&input);
        mess_matrix_clear(&output);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_Es_apply  ----- */


/**
 * @brief Initialization function of the linear mapping \f$x=(\mathcal{A}+p\mathcal{E})y\f$.
 * @param[in] eqn input  Equation object \f$\mathcal{A}+p\mathcal{E}\f$ belongs to
 * @param[in] parameters input Vector containing the shift parameters
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApE_pre function initializes the \f$x=(\mathcal{A}+\mathcal{E})y\f$ operation of a given
 * equation object. It is called normally once before the operator is uses the first time.
 * If the equation object does not provide an initialization functions the @ref mess_equation_ApE_pre
 * function returns immediately.
 *
 * @sa mess_equation_ApE_apply
 * @sa mess_equation_ApE_post
 *
 */
int mess_equation_ApE_pre ( mess_equation eqn, mess_vector parameters )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(eqn);

    if ( eqn->ApEX.generate == NULL) {
        return 0;
    }

    if ( eqn->ApEX.to_clear > 0 ) {
        MSG_INFO("The equation object's x = (A+pE)y function is already initialized.\n");
    }

    ret = eqn->ApEX.generate(eqn, parameters);
    eqn->ApEX.to_clear = 1;
    return ret;
}       /* -----  end of function mess_equation_ApE_pre  ----- */


/**
 * @brief Finalizes/Clears the linear mapping \f$x=(\mathcal{A}+\mathcal{E})y\f$.
 * @param[in]   eqn input     Equation object \f$x=(\mathcal{A}+\mathcal{E})y\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApE_post function finalizes the \f$x=(\mathcal{A}+\mathcal{E})y\f$ operator. That
 * means it calls the clear function of the equation object if it exits otherwise it will
 * return immediately.
 *
 * @sa mess_equation_ApE_pre
 * @sa mess_equation_ApE_apply
 *
 */
int mess_equation_ApE_post ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(eqn);

    if ( eqn->ApEX.to_clear == 0 ) {
        return 0;
    }

    ret = eqn->ApEX.clear(eqn);
    eqn->ApEX.to_clear = 0;
    return ret ;
}       /* -----  end of function mess_equation_ApE_post  ----- */



/**
 * @brief Compute \f$x =( \mathcal{A}+p\mathcal{E})y\f$.
 * @param[in] eqn input equation object \f$x =( \mathcal{A}+p\mathcal{E})y\f$
 * @param[in] op input  opertion to apply
 * @param[in] p input   current shift parameter
 * @param[in] p_idx input postion of the current parameter in the previously passed parameter set
 * @param[in] in input  matrix \f$y\f$
 * @param[in] out input matrix \f$x\f$
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApE_apply function computes \f$ x = op(\mathcal{A}+p\mathcal{E}) x\f$, where
 * x and y are matrices. For vectors use \ref mess_equation_ApE_apply_vector instead. If
 * the object does not provide a ApE.apply function, an error is returned. The parameter is given by
 * its value p and the postition in the parameter vector passed to the \ref mess_equation_ApE_pre function before.
 * If the index p_idx is lower than 0, the operator should be implemented such that the value of p is used.
 *
 * @sa mess_equation_A_apply
 * @sa mess_equation_E_apply
 *
 */
int mess_equation_ApE_apply ( mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_matrix in, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->ApEX.apply == NULL ) {
        MSG_ERROR("y = (Ap+E) not implemented.\n");
        return MESS_ERROR_NOSUPPORT;
    } else {
        ret = eqn->ApEX.apply(eqn, op, p, p_idx, in, out);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->ApEX.apply);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_ApE_apply  ----- */

/**
 * @brief Compute \f$x =( \mathcal{A}+p\mathcal{E})y\f$ with vectors..
 * @param[in] eqn input Equation object \f$x =( \mathcal{A}+p\mathcal{E})y\f$
 * @param[in] op input  Opertion to apply
 * @param[in] p input   Current shift parameter
 * @param[in] p_idx input Postiion of the current parameter in the previously passed parameter set
 * @param[in] in input  Input vector  y
 * @param[in] out input Ouput vector  x
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApE_apply_vector function computes \f$ x = op(\mathcal{A}+p\mathcal{E}) x\f$, where
 * x and y are vectors. For matrices  use \ref mess_equation_ApE_apply instead. If
 * the object does not provide a ApE.apply function, an error is returned. The parameter is given by
 * its value p and the postition in the parameter vector passed to the \ref mess_equation_ApE_pre function before.
 * If the index p_idx is lower than 0, the operator should be implemented such that the value of p is used.
 *
 * @sa mess_equation_A_apply_vector
 * @sa mess_equation_E_apply_vector
 *
 */
int mess_equation_ApE_apply_vector( mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_vector in, mess_vector out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->ApEX.apply == NULL ) {
        MSG_ERROR("y = (Ap+E) not implemented.\n");
        return MESS_ERROR_NOSUPPORT;
    } else {
        mess_matrix input, output;
        ret = mess_matrix_init(&input);             FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&output);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_vector_tomatrix(in, input);            FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tomatrix);
        ret = eqn->ApEX.apply(eqn, op, p, p_idx, input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
        ret = mess_vector_frommatrix(output, out);        FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_frommatrix);
        mess_matrix_clear(&input);
        mess_matrix_clear(&output);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_ApE_apply_vector  ----- */


/**
 * @brief Initialization function of the linear mapping \f$(\mathcal{A}+p\mathcal{E})x=y\f$.
 * @param[in] eqn input  Equation object \f$\mathcal{A}+p\mathcal{E}\f$ belongs to
 * @param[in] parameters input Vector containing the shift parameters
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApEs_pre function initializes the \f$(\mathcal{A}+\mathcal{E})x=y\f$ operation of a given
 * equation object. It is called normally once before the operator is uses the first time.
 * If the equation object does not provide an initialization functions the @ref mess_equation_ApEs_pre
 * function returns immediately.
 *
 * @sa mess_equation_ApEs_apply
 * @sa mess_equation_ApEs_post
 *
 */
int mess_equation_ApEs_pre ( mess_equation eqn, mess_vector parameters )
{
    MSG_FNAME(__func__);
    int ret = 0;

    mess_check_nullpointer(eqn);

    if ( eqn->ApEINV.generate == NULL) {
        return 0;
    }

    if ( eqn->ApEINV.to_clear > 0 ) {
        MSG_INFO("The equation object's x = (A+pE)y function is already initialized.\n");
    }

    ret = eqn->ApEINV.generate(eqn, parameters);
    eqn->ApEINV.to_clear = 1;
    return ret;
}       /* -----  end of function mess_equation_ApEs_pre  ----- */


/**
 * @brief Finalizes/Clears the linear mapping \f$(\mathcal{A}+\mathcal{E})x=y\f$.
 * @param[in]   eqn input     Equation object \f$(\mathcal{A}+\mathcal{E})x=y\f$ belongs to
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApEs_post function finalizes the \f$(\mathcal{A}+\mathcal{E})x=y\f$ operator. That
 * means it calls the clear function of the equation object if it exits otherwise it will
 * return immediately.
 *
 * @sa mess_equation_ApEs_pre
 * @sa mess_equation_ApEs_apply
 *
 */
int mess_equation_ApEs_post ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    int ret;

    mess_check_nullpointer(eqn);

    if ( eqn->ApEINV.to_clear == 0 ) {
        return 0;
    }

    ret = eqn->ApEINV.clear(eqn);
    eqn->ApEINV.to_clear = 0;
    return ret ;
}       /* -----  end of function mess_equation_ApEs_post  ----- */



/**
 * @brief Compute \f$( \mathcal{A}+p\mathcal{E})x=y\f$.
 * @param[in] eqn input equation object \f$( \mathcal{A}+p\mathcal{E})x=y\f$
 * @param[in] op input  opertion to apply
 * @param[in] p input   current shift parameter
 * @param[in] p_idx input postion of the current parameter in the previously passed parameter set
 * @param[in] in input  matrix \f$y\f$
 * @param[in] out input matrix \f$x\f$
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApEs_apply function computes \f$ op(\mathcal{A}+p\mathcal{E})x =y\f$, where
 * x and y are matrices. For vectors use \ref mess_equation_ApEs_apply_vector instead. If
 * the object does not provide a ApE.apply function, an error is returned. The parameter is given by
 * its value p and the postition in the parameter vector passed to the \ref mess_equation_ApEs_pre function before.
 * If the index p_idx is lower than 0, the operator should be implemented such that the value of p is used.
 *
 * @sa mess_equation_As_apply
 * @sa mess_equation_Es_apply
 *
 */
int mess_equation_ApEs_apply ( mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_matrix in, mess_matrix out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->ApEINV.apply == NULL ) {
        MSG_ERROR("y = (Ap+E) not implemented.\n");
        return MESS_ERROR_NOSUPPORT;
    } else {
        ret = eqn->ApEINV.apply(eqn, op, p, p_idx, in, out);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->ApEINV.apply);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_ApEs_apply  ----- */

/**
 * @brief Compute \f$( \mathcal{A}+p\mathcal{E})x=y\f$ with vectors..
 * @param[in] eqn input Equation object \f$( \mathcal{A}+p\mathcal{E})x=y\f$
 * @param[in] op input  Opertion to apply
 * @param[in] p input   Current shift parameter
 * @param[in] p_idx input Postion of the current parameter in the previously passed parameter set
 * @param[in] in input  vector  \f$y\f$
 * @param[out] out output vector  \f$x\f$
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_ApEs_apply_vector function computes \f$ op(\mathcal{A}+p\mathcal{E})x =y\f$, where
 * x and y are vectors. For matrices  use \ref mess_equation_ApEs_apply instead. If
 * the object does not provide a ApE.apply function, an error is returned. The parameter is given by
 * its value p and the postition in the parameter vector passed to the \ref mess_equation_ApEs_pre function before.
 * If the index p_idx is lower than 0, the operator should be implemented such that the value of p is used.
 *
 * @sa mess_equation_A_apply_vector
 * @sa mess_equation_E_apply_vector
 *
 */
int mess_equation_ApEs_apply_vector( mess_equation eqn, mess_operation_t op, mess_double_cpx_t p , mess_int_t p_idx, mess_vector in, mess_vector out )
{
    MSG_FNAME(__func__);
    int ret = 0;

    /* Check input   */
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(in);
    mess_check_nullpointer(out);
    mess_check_real_or_complex(in);

    if ( eqn->ApEINV.apply == NULL ) {
        MSG_ERROR("y = (Ap+E) not implemented.\n");
        return MESS_ERROR_NOSUPPORT;
    } else {
        mess_matrix input, output;
        ret = mess_matrix_init(&input);                                 FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_matrix_init(&output);                                FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
        ret = mess_vector_tomatrix(in, input);                          FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_tomatrix);
        ret = eqn->ApEINV.apply(eqn, op, p, p_idx, input, output);      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), eqn->AX.apply);
        ret = mess_vector_frommatrix(output, out);                      FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_vector_frommatrix);
        mess_matrix_clear(&input);
        mess_matrix_clear(&output);
        return 0;
    }
    return 0;
}       /* -----  end of function mess_equation_ApEs_apply_vector  ----- */


