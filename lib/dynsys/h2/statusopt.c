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
 * @file lib/dynsys/h2/statusopt.c
 * @brief Status and options for \f$ \mathcal{H}_2 \f$ things.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"


/**
 * @brief Initialize a mess_h2_options object.
 * @param[in,out] opt pointer to the object
 * @return always zero
 *
 * The @ref mess_h2_options_init function initializes a mess_h2_options object.
 *
 */
int mess_h2_options_init ( mess_h2_options *opt )
{
    MSG_FNAME(__func__);
    mess_check_nullpointer(opt);


    mess_try_alloc(*opt, mess_h2_options , sizeof(struct mess_h2_options_st));

    (*opt)->maxit = 10;
    (*opt)->tol = 1e-7;
    (*opt)->output = 1;
    (*opt)->calc_h2err = 0;
    (*opt)->calc_finalh2 = 1;
    (*opt)->stepdebug = NULL;
    (*opt)->stepdebug_data  = NULL;
    (*opt)->rdim = 0;
    return 0;
}       /* -----  end of function mess_h2_options_init  ----- */


/**
 * @brief Clean up a mess_h2_options object.
 * @param[in] opt  input pointer to the object
 * @return always zero
 *
 * The @ref mess_h2_options_clear function cleans up a mess_h2_options object
 *
 */
int  mess_h2_options_clear ( mess_h2_options *opt )
{
    if ( opt == NULL) return 0;
    if ( *opt == NULL) return 0;
    mess_free(*opt);
    *opt=NULL;
    return 0;
}       /* -----  end of function mess_h2_options_clear  ----- */


/**
 * @brief Print a mess_h2_options object to screen.
 * @param[in] opt  input mess_h2_options object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_h2_options_print function prints a mess_h2_options object to the screen.
 *
 */
int mess_h2_options_print ( mess_h2_options opt )
{
    MSG_FNAME(__func__);

    mess_check_nullpointer(opt);

    MSG_PRINT("maxit:         " MESS_PRINTF_INT "\n", opt->maxit);
    MSG_PRINT("tol:           %lg\n", opt->tol);
    MSG_PRINT("output:        %d\n", opt->output);
    MSG_PRINT("calc_h2err:    ");
    switch(opt->calc_h2err){
        case MESS_H2NORM_FULL:
            MSG_PRINT("FULL\n");
            break;
        case MESS_H2NORM_UPDATE:
            MSG_PRINT("UDPATE\n");
            break;
        default:
            MSG_PRINT("0\n");
            break;
    }
    MSG_PRINT("calc_finalh2:  %d\n", opt->calc_finalh2);
    MSG_PRINT("rdim:          " MESS_PRINTF_INT "\n", opt->rdim);
    return 0;
}       /* -----  end of function mess_h2_options_print  ----- */


/**
 * @brief Initialize a mess_h2_status object.
 * @param[in,out] status pointer to the status object
 * @return always zero
 *
 * The @ref mess_h2_status_init function initializes a mess_h2_status object.
 *
 */
int mess_h2_status_init ( mess_h2_status *status ){
    MSG_FNAME(__func__);
    int ret = 0;


    mess_check_nullpointer(status);
    mess_try_alloc(*status, mess_h2_status, sizeof(struct mess_h2_status_st));

    MESS_INIT_VECTORS(&((*status)->h2err), &((*status)->sigmadiff))
    ret = mess_vector_alloc(((*status)->h2err),1,MESS_REAL);                FUNCTION_FAILURE_HANDLE(ret, (ret!=0),mess_vector_alloc);
    ret = mess_vector_alloc(((*status)->sigmadiff),1,MESS_REAL);            FUNCTION_FAILURE_HANDLE(ret,(ret!=0), mess_vector_alloc);

    (*status)->it = 0;
    (*status)->finalh2 = 0;
    (*status)->finalsigma = 0;
    (*status)->time = 0;
    (*status)->cancel_sigma =0;
    (*status)->cancel_it = 0;
    return 0;
}       /* -----  end of function mess_h2_status_init  ----- */

/**
 * @brief Clean up a mess_h2_status object.
 * @param[in] status  input pointer to the object
 * @return always zero
 *
 * The @ref mess_h2_status_clear function cleans up a mess_h2_status object.
 *
 */
int mess_h2_status_clear ( mess_h2_status *status )
{
    if ( status == NULL)
        return 0;
    if ( *status == NULL)
        return 0;

    mess_vector_clear(&((*status)->h2err));
    mess_vector_clear(&((*status)->sigmadiff));

    mess_free(*status);
    *status = NULL;
    return 0;
}       /* -----  end of function mess_h2_status_clear  ----- */

/**
 * @brief Print a mess_h2_status object.
 * @param[in] status  input mess_h2_status object
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_h2_status_print function prints a mess_h2_status object on the screen.
 *
 */
int  mess_h2_status_print ( mess_h2_status status )
{
    MSG_FNAME(__func__);

    mess_check_nullpointer(status);

    MSG_PRINT("time:                 %lg\n", status->time);
    MSG_PRINT("number of iterations: " MESS_PRINTF_INT "\n", status->it);
    MSG_PRINT("final sigmadiff:      %lg\n", status->finalsigma);
    MSG_PRINT("final h2err:          %lg\n", status->finalh2);
    MSG_PRINT("cancel_sigma:         %d\n", status->cancel_sigma);
    MSG_PRINT("cancel_it:            %d\n", status->cancel_it);

    return 0;
}       /* -----  end of function mess_h2_status_print  ----- */

/**
 * @brief Print a mess_h2_status object (full output).
 * @param[in] status  input mess_h2_status object
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_h2_status_printfull function prints a mess_h2_status object on the screen.
 *
 */
int  mess_h2_status_printfull ( mess_h2_status status )
{
    MSG_FNAME(__func__);

    mess_check_nullpointer(status);

    MSG_PRINT("time:                 %lg\n", status->time);
    MSG_PRINT("number of iterations: " MESS_PRINTF_INT "\n", status->it);
    MSG_PRINT("final sigmadiff:      %lg\n", status->finalsigma);
    MSG_PRINT("final h2err:          %lg\n", status->finalh2);
    MSG_PRINT("cancel_sigma:         %d\n", status->cancel_sigma);
    MSG_PRINT("cancel_it:            %d\n", status->cancel_it);

    MSG_PRINT("h2-errors: \n");
    mess_vector_print(status->h2err);
    MSG_PRINT("sigma-diff:\n");
    mess_vector_print(status->sigmadiff);

    return 0;
}       /* -----  end of function mess_h2_status_print  ----- */


/**
 * @brief Write a mfile to plot in @matlab .
 * @param[in] status  input mess_h2_status object
 * @param[in] filename  input name of the mfile
 * @return zero on success or a non-zero error value
 *
 * The @ref mess_h2_status_write_mfile function writes a mfile to get the
 * status data easily to @matlab.
 *
 */
int  mess_h2_status_write_mfile ( mess_h2_status status, const char * filename )
{
    MSG_FNAME(__func__);
    mess_int_t i;
    int err;
    FILE *f;

    /*-----------------------------------------------------------------------------
     *  check input
     *-----------------------------------------------------------------------------*/
    mess_check_nullpointer(status);
    mess_check_nullpointer(filename);

    f = fopen(filename,"w");
    if ( f == NULL){
        err = errno;
        MSG_ERROR ("opening file: %s failed\n", filename);
        MSG_ERROR ("errno: %03d - %s\n", err, strerror(err));
        return MESS_ERROR_FILEIO;
    }
    fprintf(f,"finalsigma=%.15e;\n",status->finalsigma);
    fprintf(f,"finalh2=%.15e;\n",status->finalh2);
    fprintf(f,"time=%.15e;\n",status->time);
    fprintf(f,"h2err=[];\n");
    for (i=0;i < status->h2err->dim;i++){
        fprintf(f,"h2err(" MESS_PRINTF_INT ")=%.15e;\n",i+1, status->h2err->values[i]);
    }
    fprintf(f,"sigmadiff=[];\n");
    for (i=0;i < status->sigmadiff->dim;i++){
        fprintf(f,"sigmadiff(" MESS_PRINTF_INT ")=%.15e;\n",i+1, status->sigmadiff->values[i]);
    }

    fprintf(f,"figure;\n");
    fprintf(f,"semilogy(0:" MESS_PRINTF_INT ",sigmadiff,'r-',0:" MESS_PRINTF_INT ",h2err,'b-*');\n", status->sigmadiff->dim-1, status->h2err->dim-1);
    fprintf(f,"xlabel('iteration number');\n");
    fprintf(f,"title('convergence IRKA');\n");
    fprintf(f,"legend('||sigma-sigmaold||/||sigma||','H_2 error','Location','NorthEast');\n");



    fclose(f);

    return 0;
}       /* -----  end of function mess_h2_status_write_mfile  ----- */

