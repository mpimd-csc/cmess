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
 * @file include/mess/h2.h
 * @brief All functions for \f$ \mathcal{H}_2 \f$ model order reduction and balanced truncation.
 * @author @koehlerm
 */

#ifndef H2_H_
#define H2_H_

#ifdef  __cplusplus
extern "C" {
#endif
    /** @addtogroup dynsys_mor_h2
     * @{ */
    /**
     * @brief Structure to store \f$ \mathcal{H}_2 \f$ options.
     *
     * This structure stores the information for the \f$ \mathcal{H}_2 \f$ model order reduction
     * algorithms. It never be used directly but with the typedef mess_h2_options.
     */
    struct mess_h2_options_st {
        mess_int_t maxit;                                                                   /**< Maximum number of iteration */
        mess_int_t rdim;                                                                    /**< Dimension of the reduced order model.  */
        void (*stepdebug)(void *data, mess_int_t it, double sigmadiff, double h2err);       /**< Stepdebug function which is called every iteration, if wanted. */
        void *stepdebug_data;                                                               /**< Auxillary data for the debug function  */
        double tol;                                                                         /**< Cancellation tolerance */
        int output;                                                                         /**< Should output be created. */
        int calc_h2err;                                                                     /**< Should the \f$ \mathcal{H}_2 \f$ error be computed during the iteration.  */
        int calc_finalh2;                                                                   /**< Should a final  \f$ \mathcal{H}_2 \f$  error be computed. */
    };

    /**
     *  @brief Type definition for the \f$ \mathcal{H}_2 \f$ options.
     *
     * Type definition for better use of @ref mess_h2_options_st .
     */
    typedef struct mess_h2_options_st*  mess_h2_options;

    /**
     * @brief Structure to store \f$ \mathcal{H}_2 \f$ status.
     *
     * This structure stores the information for the status of the \f$ \mathcal{H}_2 \f$ error, the difference
     * between a system and its reduced one, the number of (needed) iterations and the time.
     */
    struct mess_h2_status_st {
        double finalsigma;          /**< Difference between a system and its reduced one in the end */
        double finalh2;             /**<  \f$ \mathcal{H}_2 \f$-error in the end */
        double time;                /**< Required time */
        mess_vector h2err;          /**< \f$ \mathcal{H}_2 \f$-error */
        mess_vector sigmadiff;      /**<  Difference between a system and its reduced one */
        mess_int_t it;              /**< Number of iterations */
        int cancel_sigma;           /**< Boolean value, 1 if tolerance is reached.  */
        int cancel_it;              /**< Number of iteration where the algorithm stops, 1 if tolerance is not reached. */
    };

    /**
     * @brief Type definition for the \f$ \mathcal{H}_2 \f$ status.
     *
     * The @ref mess_h2_status  is a type definition for the \f$ \mathcal{H}_2 \f$ status.
     *
     **/
    typedef struct mess_h2_status_st*  mess_h2_status;

    /** Do not calculate \f$ \mathcal{H}_2 \f$ norm. */
#define MESS_H2NORM_NO      0
    /** Calculate full \f$ \mathcal{H}_2 \f$ norm. */
#define MESS_H2NORM_FULL        1
    /** Calculate \f$ \mathcal{H}_2 \f$ norm with an update. */
#define MESS_H2NORM_UPDATE  2


    /*-----------------------------------------------------------------------------
     *  H2 norm and error
     *-----------------------------------------------------------------------------*/

    int mess_h2_norm_internal ( mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E, double *norm);
    int mess_h2_norm ( mess_dynsys LTI, double *norm );

    int mess_h2_error_internal ( mess_matrix A, mess_matrix B, mess_matrix C, mess_matrix E, mess_matrix Ar, mess_matrix Br, mess_matrix Cr , mess_matrix Er, double *norm);
    int mess_h2_error ( mess_dynsys ltiA, mess_dynsys ltiB,double *err );


    /*-----------------------------------------------------------------------------
     *  H2 Algorithmns
     *-----------------------------------------------------------------------------*/
    int mess_h2_irka_init ( mess_matrix A, mess_direct Asolver, mess_matrix E, mess_int_t* r0, mess_int_t kp, mess_int_t km, mess_vector sigma );
    int  mess_h2_irka_biorth (mess_dynsys orig, mess_vector sigma, mess_h2_options opt ,
            mess_dynsys reduced, mess_matrix V, mess_matrix W, mess_h2_status status);

    int mess_h2_tsia (  mess_dynsys orig, mess_vector sigmaEX, mess_h2_options opt,
            mess_dynsys redu, mess_matrix V, mess_matrix W, mess_h2_status status);

    int mess_h2_tsiag (     mess_dynsys orig, mess_vector sigmaEX, mess_h2_options opt,
            mess_dynsys redu, mess_matrix V, mess_matrix W, mess_h2_status status);

    /*-----------------------------------------------------------------------------
     *  H2 options
     *-----------------------------------------------------------------------------*/
    int mess_h2_options_init ( mess_h2_options *opt );
    int mess_h2_options_clear ( mess_h2_options *opt );
    int mess_h2_options_print ( mess_h2_options opt );


    /*-----------------------------------------------------------------------------
     *  H2 status
     *-----------------------------------------------------------------------------*/
    int mess_h2_status_init ( mess_h2_status *status );
    int mess_h2_status_clear ( mess_h2_status *status );
    int  mess_h2_status_print ( mess_h2_status status );
    int  mess_h2_status_printfull ( mess_h2_status status );
    int  mess_h2_status_write_mfile ( mess_h2_status status, const char * filename );

    /** @} */


    /** @addtogroup dynsys_mor_bt
     * @{ */
    /*-----------------------------------------------------------------------------
     * Balanced Truncation Stuff.
     *-----------------------------------------------------------------------------*/
    /** Type of second order system: position.  */
#define MESS_BT_POSITION 0x01
    /** Type of second order system: velocity.  */
#define MESS_BT_VELOCITY 0x02
    /** Type of second order system: position-velocity.  */
#define MESS_BT_POSITION_VELOCITY 0x03
    /** Type of second order system: velocity-position.  */
#define MESS_BT_VELOCITY_POSITION 0x04

    /**
     * @brief Type definition for the choose order function.
     *
     * To give the user the possibility to choose the order of the balanced truncation
     * reduced model as they want, they can use their own chooser function that
     * looks like this type definition. \n
     * In the case of success it needs to return \f$ 0 \f$. All
     * other return values are interpreted as a failure.
     */
    typedef int (*mess_bt_chooseorder_func) (mess_vector, double*, mess_int_t *, mess_int_t, void*);

    /**
     * @brief Options structure for balanced truncation .
     *
     * This structure stores the information for the balanced truncation  model order reduction
     * algorithms.
     */
    struct mess_bt_options_st {
        double tol;                             /**< Cut off tolerance for Hankel singular values */
        mess_bt_chooseorder_func chooseorder;   /**< Function to choose the reduced order in your own way  */
        void *chooseorder_aux;                  /**< Additional data for the chooseorder function. */
        mess_int_t rdim;                        /**< Maximum dimension of the reduced order model */
        unsigned short so_type;                 /**< Type of second order system */
    };

    /**
     * @brief Type definition for the balanced truncation options.
     *
     * Type definition for better use of @ref mess_bt_options . */
    typedef struct mess_bt_options_st*  mess_bt_options;

    /**
     * @brief Status structure for balanced truncation.
     *
     * This structure stores the information for the status of the balanced truncation algorithm.
     */
    struct mess_bt_status_st {
        double time;                /**< Overall time */
        double time_lyap;           /**< Time spent on solving Lyapunov Equations */
        double time_VW;             /**< Time spent on computing \f$ V \f$ and \f$ W \f$ */
        double esterror;            /**< Estimated error */
        mess_status statB;          /**<  LRCF-ADI status for \f$ AX+XA^T=-BB^T \f$ */
        mess_status statC;          /**<  LRCF-ADI status for \f$ A^TX+XA=-C^TC \f$ */
        mess_int_t rdim;            /**< Maximum dimension of the reduced order model */
    };

    /**
     * @brief Type definition for the balanced truncation status.
     *
     * The @ref mess_h2_status  is a type definition for the balanced truncation status.
     *
     **/
    typedef struct mess_bt_status_st*  mess_bt_status;


    /*-----------------------------------------------------------------------------
     *  BT options
     *-----------------------------------------------------------------------------*/
    int mess_bt_options_init(mess_bt_options *opt);
    int mess_bt_options_clear(mess_bt_options *opt);
    int mess_bt_options_print(mess_bt_options opt);
    const char *  mess_bt_getsotypestr( unsigned short so_type  );


    /*-----------------------------------------------------------------------------
     *  BT status
     *-----------------------------------------------------------------------------*/
    int mess_bt_status_init( mess_bt_status *status);
    int mess_bt_status_clear( mess_bt_status *status);
    int mess_bt_status_print( mess_bt_status status);

    /*-----------------------------------------------------------------------------
     *  BT algorithmns
     *-----------------------------------------------------------------------------*/
    int mess_bt_chooseorder_default(mess_vector SIGMA, double *tol, mess_int_t* maxr, mess_int_t FOM, void*aux);
    int mess_bt_chooseorder_minreal ( mess_vector SIGMA, double *tol, mess_int_t* maxr, mess_int_t FOM, void*aux);

    int mess_bt_lrsrm(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status);
    int mess_bt_gslrsrm(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status);
    int mess_bt_gglrsrm(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status);
    int mess_bt_lrsrm_somor(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status) ;
    int mess_bt_lrsrm_somor_so(mess_dynsys sys, mess_bt_options btopt, mess_options adiopt, mess_matrix V, mess_matrix W, mess_bt_status status);
    int mess_bt_lrsrm_somor_so_gram(mess_dynsys sys, mess_matrix ZB, mess_matrix ZC, mess_bt_options btopt, mess_matrix V, mess_matrix W, mess_bt_status status);

    /** @} */

#ifdef __cplusplus
}
#endif

#endif /* LRCF_ADI_H_ */
/** \}@ */
