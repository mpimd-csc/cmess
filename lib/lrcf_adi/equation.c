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
 * @file lib/lrcf_adi/equation.c
 * @brief Manage the @ref mess_equation structure.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"




/**
 * @brief Initialize an empty mess_equation object.
 * @param[in,out] eqn pointer to equation structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_init function initializes an empty mess_equation object. \n
 * In order to setup a proper equation use one of the generator functions from
 * \ref lrcfadi_eqn.
 *
 */
int mess_equation_init(mess_equation * eqn)
{
    MSG_FNAME(__func__);
    mess_equation e;

    mess_try_alloc( *eqn , struct mess_equation_st*, sizeof (struct mess_equation_st ));
    memset(*eqn, 0, sizeof(struct mess_equation_st));
    e = *eqn;

    e->eqn_type         = MESS_EQN_NONE;
    return 0;
}

/**
 * @brief Clean a mess_equation structure.
 * @param[in,out] eqn pointer to the mess_equation structure
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_clear function clears a mess_equation structure. \n
 * It calls all internal clean up handlers to remove all created data properly from
 * the memory. \n
 * If it is necessary to call the clean up handlers without removing
 * the whole equation from memory use \ref mess_equation_clean instead.
 *
 * \sa mess_equation_clean
 */
int mess_equation_clear(mess_equation * eqn){
    mess_equation e;
    if (eqn == NULL) return 0 ;
    if (*eqn == NULL) return 0;
    e = *eqn;

    if ( e->AX.to_clear){
        e->AX.clear(*eqn);
    }

    if ( e->EX.to_clear){
        e->EX.clear(*eqn);
    }

    if ( e->AINV.to_clear){
        e->AINV.clear(*eqn);
    }

    if ( e->EINV.to_clear){
        e->EINV.clear(*eqn);
    }

    if ( e->ApEX.to_clear){
        e->ApEX.clear(*eqn);
    }

    if ( e->ApEINV.to_clear){
        e->ApEINV.clear(*eqn);
    }

    if ( e->clearRHS ) {
        mess_matrix_clear(&(e->RHS));
    }

    if ( e->clearB ) {
        mess_matrix_clear(&(e->B));
    }

    if ( e->clearC ) {
        mess_matrix_clear(&(e->C));
    }

    if ((*eqn)->clear != NULL ) {
        (*eqn)->clear((*eqn));
    }
    mess_free(*eqn);
    *eqn = NULL;
    return 0;
}



/**
 * @brief Call cleanup functions for all internally generated data.
 * @param[in,out] eqn mess_equation object to clean
 * @return zero on success or a non zero error code
 *
 * The @ref mess_equation_clean function cleans the internally generated data of a
 * @ref mess_equation object, i.e. it calls the clear function for all internally
 * defined operations without removing the whole equation from memory. \n
 * This function is useful if the equation is used later on again but the memory footprint of
 * the internal data is too big. \n
 * If the equation should be removed completely
 * from memory use \ref mess_equation_clear instead.
 *
 * \sa mess_equation_clear
 */
int  mess_equation_clean ( mess_equation eqn )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);

    if ( eqn->AX.to_clear){
        eqn->AX.clear(eqn);
    }

    if ( eqn->EX.to_clear){
        eqn->EX.clear(eqn);
    }

    if ( eqn->AINV.to_clear){
        eqn->AINV.clear(eqn);
    }

    if ( eqn->EINV.to_clear){
        eqn->EINV.clear(eqn);
    }

    if ( eqn->ApEX.to_clear){
        eqn->ApEX.clear(eqn);
    }

    if ( eqn->ApEINV.to_clear){
        eqn->ApEINV.clear(eqn);
    }

    return 0;
}       /* -----  end of function mess_equation_clean  ----- */


/**
 * @brief Set a new right hand side in a mess_equation object.
 * @param[in,out] eqn   mess_equation object to modify
 * @param[in]   opt input   options structure of the current solver
 * @param[in]   rhs input     new right hand side
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_set_rhs function sets a new right hand side inside a
 * @ref mess_equation object. \n
 * The right hand side need to fit to the requirements of the equation type in \ref mess_options. \n
 * If @ref mess_options.type is @ref MESS_OP_NONE the right hand side need to be a tall and skinny matrix with the correct number of rows.
 * Otherwise it have to be the transposed.
 *
 * \sa mess_equation_set_rhs2
 *
 */
int mess_equation_set_rhs(mess_equation eqn, mess_options opt,  mess_matrix rhs){
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(rhs);
    mess_check_real(rhs);
    int ret = 0;


    switch ( eqn->eqn_type ) {
        /*-----------------------------------------------------------------------------
         *  Set new right hand side for the Lyapunov Equation
         *-----------------------------------------------------------------------------*/
        case MESS_EQN_LYAP:
        case MESS_EQN_GLYAP:
            if ( opt->type == MESS_OP_NONE &&  eqn->dim != rhs->rows) {
                MSG_ERROR("Equation and rhs have the wrong number of rows. eqn->dim = " MESS_PRINTF_INT ", rhs->rows = " MESS_PRINTF_INT"\n", eqn->dim, rhs->rows);
                return MESS_ERROR_DIMENSION;
            }
            if ( opt->type == MESS_OP_TRANSPOSE && eqn->dim != rhs->cols) {
                MSG_ERROR("Equation and rhs do not have fitting dimensions.  eqn->dim = " MESS_PRINTF_INT ", rhs->cols = " MESS_PRINTF_INT"\n", eqn->dim, rhs->cols);
                return MESS_ERROR_DIMENSION;
            }
            if ( opt->type == MESS_OP_NONE ) {
                if ( rhs->cols >= rhs->rows) {
                    MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(opt->type), (int) rhs->rows, (int) rhs->cols);
                }
            }
            if ( opt->type == MESS_OP_TRANSPOSE ) {
                if ( rhs->rows >= rhs->cols) {
                    MSG_WARN("Oversized right hand side factor for LRCF-ADI (op = %s, rows=%d, cols=%d)\n ", mess_operation_t_str(opt->type), (int) rhs->rows, (int) rhs->cols);
                }
            }
            if (eqn->clearRHS) {
                mess_matrix_clear(&eqn->RHS);
                eqn->RHS=NULL;
            }
            if ( opt->type == MESS_OP_NONE) {
                eqn->clearRHS = 0;
                eqn->RHS = rhs;
            } else {
                ret = mess_matrix_init(&eqn->RHS);  FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_init);
                eqn->clearRHS = 1;
                ret = mess_matrix_ctranspose(rhs, eqn->RHS); FUNCTION_FAILURE_HANDLE(ret, (ret!=0), mess_matrix_ctranspose);
            }
            break;
        default:
            if ( eqn->dim != rhs->rows){
                MSG_ERROR("The right hand side must have the same dimension as the equation.\n");
                return MESS_ERROR_DIMENSION;
            }
            if (eqn->clearRHS) {
                mess_matrix_clear(&eqn->RHS);
                eqn->RHS=NULL;
            }
            eqn->RHS = rhs;
    }

    return 0;
}

/**
 * @brief Set a new right hand side in a mess_equation object without additional checks.
 * @param[in,out] eqn   mess_equation object to modify
 * @param[in]   rhs input     new right hand side
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_equation_set_rhs2 function sets a new right hand side inside a
 * @ref mess_equation object without any sanity checks, except of checking for the
 * equation type.
 *
 * \sa mess_equation_set_rhs
 *
 */
int mess_equation_set_rhs2(mess_equation eqn, mess_matrix rhs){
    MSG_FNAME(__func__);
    mess_check_nullpointer(eqn);
    mess_check_nullpointer(rhs);
    mess_check_real(rhs);


    switch ( eqn->eqn_type ) {
        /*-----------------------------------------------------------------------------
         *  Set new right hand side for the Lyapunov Equation
         *-----------------------------------------------------------------------------*/
        case MESS_EQN_LYAP:
        case MESS_EQN_GLYAP:
            if (eqn->clearRHS) {
                mess_matrix_clear(&eqn->RHS);
                eqn->RHS=NULL;
            }
            eqn->clearRHS = 0;
            eqn->RHS = rhs;
        break;
        default:
            if (eqn->clearRHS) {
                mess_matrix_clear(&eqn->RHS);
                eqn->RHS=NULL;
            }
            eqn->RHS = rhs;
            eqn->clearRHS = 0;
    }

    return 0;
}


/**
 * @brief Return the dimension of the equation of object.
 * @param[in] eqn input   Equation object
 * @return The dimension or a negative error code
 *
 * The @ref mess_equation_dim function returns the dimension of the equation object.
 *
 */
mess_int_t  mess_equation_dim ( mess_equation eqn )
{
    if ( eqn == NULL) {
        return -MESS_ERROR_NULLPOINTER;
    }
    return eqn->dim;
}       /* -----  end of function mess_equation_dim  ----- */


/**
 * @brief Checks if \f$ x = \mathcal{A}y \f$ is implemented.
 * @param[in] eqn input
 * @return true if the operation is available, false otherwise inlcuding errors.
 *
 * The @ref mess_equation_has_A function returns true if the operator provides a
 * \f$ x = \mathcal{A}y \f$ function.
 * @sa mess_equation_has_As
 * @sa mess_equation_has_E
 * @sa mess_equation_has_Es
 * @sa mess_equation_has_ApE
 * @sa mess_equation_has_ApEs
 *
 */
mess_int_t  mess_equation_has_A ( mess_equation eqn )
{
    if ( eqn == NULL )
        return 0;
    if ( eqn->AX.apply == NULL)
        return 0;
    return 1;
}

/**
 * @brief Checks if \f$ \mathcal{A}x =y \f$ is implemented.
 * @param[in] eqn input
 * @return true if the operation is available, false otherwise inlcuding errors.
 *
 * The @ref mess_equation_has_As function returns true if the operator provides a
 * \f$ \mathcal{A}x = y \f$ function.
 * @sa mess_equation_has_A
 * @sa mess_equation_has_E
 * @sa mess_equation_has_Es
 * @sa mess_equation_has_ApE
 * @sa mess_equation_has_ApEs
 *
 */
mess_int_t  mess_equation_has_As ( mess_equation eqn )
{
    if ( eqn == NULL )
        return 0;
    if ( eqn->AINV.apply == NULL)
        return 0;
    return 1;
}


/**
 * @brief Checks if \f$ x = \mathcal{E}y \f$ is implemented.
 * @param[in] eqn input
 * @return true if the operation is available, false otherwise inlcuding errors.
 *
 * The @ref mess_equation_has_E function returns true if the operator provides a
 * \f$ x = \mathcal{E}y \f$ function.
 * @sa mess_equation_has_Es
 * @sa mess_equation_has_A
 * @sa mess_equation_has_As
 * @sa mess_equation_has_ApE
 * @sa mess_equation_has_ApEs
 *
 */
mess_int_t  mess_equation_has_E ( mess_equation eqn )
{
    if ( eqn == NULL )
        return 0;
    if ( eqn->EX.apply == NULL)
        return 0;
    return 1;
}

/**
 * @brief Checks if \f$ \mathcal{E}x =y \f$ is implemented.
 * @param[in] eqn input
 * @return true if the operation is available, false otherwise inlcuding errors.
 *
 * The @ref mess_equation_has_Es function returns true if the operator provides a
 * \f$ \mathcal{E}x = y \f$ function.
 * @sa mess_equation_has_E
 * @sa mess_equation_has_A
 * @sa mess_equation_has_As
 * @sa mess_equation_has_ApE
 * @sa mess_equation_has_ApEs
 *
 */
mess_int_t  mess_equation_has_Es ( mess_equation eqn )
{
    if ( eqn == NULL )
        return 0;
    if ( eqn->EINV.apply == NULL)
        return 0;
    return 1;
}

/**
 * @brief Checks if \f$ x = (\mathcal{A}+\mathcal{E})y \f$ is implemented.
 * @param[in] eqn input
 * @return true if the operation is available, false otherwise inlcuding errors.
 *
 * The @ref mess_equation_has_ApE function returns true if the operator provides a
 * \f$ x = (\mathcal{A}+\mathcal{E})y \f$ function.
 * @sa mess_equation_has_A
 * @sa mess_equation_has_As
 * @sa mess_equation_has_E
 * @sa mess_equation_has_Es
 * @sa mess_equation_has_ApEs
 *
 */
mess_int_t  mess_equation_has_ApE ( mess_equation eqn )
{
    if ( eqn == NULL )
        return 0;
    if ( eqn->ApEX.apply == NULL)
        return 0;
    return 1;
}

/**
 * @brief Checks if \f$ (\mathcal{A}+p\mathcal{E})x =y \f$ is implemented.
 * @param[in] eqn input
 * @return true if the operation is available, false otherwise inlcuding errors.
 *
 * The @ref mess_equation_has_ApEs function returns true if the operator provides a
 * \f$ (\mathcal{A}+p\mathcal{E})x = y \f$ function.
 * @sa mess_equation_has_A
 * @sa mess_equation_has_As
 * @sa mess_equation_has_E
 * @sa mess_equation_has_Es
 * @sa mess_equation_has_ApE
 *
 */
mess_int_t  mess_equation_has_ApEs ( mess_equation eqn )
{
    if ( eqn == NULL )
        return 0;
    if ( eqn->ApEINV.apply == NULL)
        return 0;
    return 1;
}





