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
 * @file lib/misc/time.c
 * @brief Functions for timers and time measuring.
 * @author @koehlerm
 *
 * This file contains some nice functions for easy time measuring. If your system
 * has POSIX-RT capabilities, they will be used.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>

#ifndef __USE_POSIX199309
#define __USE_POSIX199309
#include <time.h>
#undef __USE_POSIX199309
#else
#include <time.h>
#endif


#define MESS_ERROR
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"



/**
 * @brief Return a timer value.
 * @return time
 *
 * The @ref mess_wtime function returns a time as a @c double and can be used for easy runtime measuring. \n
 * The time measuring will works like the following:
 * \code{.c}
 *  double tStart, tEnd, tRun;
 *  tStart = mess_wtime();
 *  doSomething();
 *  tEnd = mess_wtime();
 *  tRun = tEnd-tStart;
 * \endcode
 *
 */
double mess_wtime ()
{
    struct timeval tv;
    gettimeofday (&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1e6;
}


/**
 * @brief Get the current CPU time.
 * @return a double value containing the enlasped CPU time
 *
 * The @ref mess_ctime function can be used for easy CPU time measuring. It returns a time as a @c double.
 * The time measuring will work like the following:
 * \code{.c}
 *  double tStart, tEnd, tRun;
 *  tStart = mess_ctime();
 *  doSomething();
 *  tEnd = mess_ctime();
 *  tRun = tEnd-tStart;
 * \endcode
 */
double mess_ctime() {
    clock_t  t = clock();
    return t / (double)(CLOCKS_PER_SEC);
}


/**
 * @brief Initialize a timer object.
 * @param[in,out] timer pointer to timer object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_timer_init function initializes a @ref mess_timer object and starts the timer.\n
 * It uses the POSIX Realtime Extensisons POSIX.1-2001.
 *
 * \sa mess_wtime
 */
int mess_timer_init(mess_timer  *timer){
    MSG_FNAME(__func__);
    mess_try_alloc((*timer), mess_timer, sizeof(struct mess_timer_st));
    mess_try_alloc((*timer)->timer, struct timespec *, sizeof(struct timespec));
#ifdef MESS_HAVE_CLOCK_GETTIME
    if ( clock_gettime(CLOCK_MONOTONIC,(*timer)->timer) != 0 ) {
        MSG_ERROR("An error occured during initializing the clock.\n");
        return (MESS_ERROR_MISC);
    }
#else
    (*timer)->dtimer = mess_wtime();
#endif
    return (0);
}


/**
 * @brief Clean up a mess_timer object.
 * @param[in,out] timer pointer to a mess_timer object
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_timer_clear function clears a @ref mess_timer object from memory.
 *
 */
int mess_timer_clear(mess_timer *timer){
    //  MESS_FN_HEADER;
    if ( timer == NULL ) return 0;
    if ( *timer == NULL ) return 0;
    mess_free((*timer)->timer);
    mess_free(*timer);
    *timer = NULL;
    return 0;
}

/**
 * @brief Get current elasped time since starting the timer.
 * @param[in] timer input  mess_timer object
 * @return elapsed time since starting the timer
 *
 * The @ref mess_timer_get function returns the elapsed time since the given timer
 * was started or reseted. \n
 * Time will be returned in seconds.
 *
 */
double mess_timer_get(mess_timer timer){
#ifdef MESS_HAVE_CLOCK_GETTIME
    MSG_FNAME(__func__);
    struct timespec endTimespec;
    if ( clock_gettime(CLOCK_MONOTONIC,&endTimespec) != 0 ) {
        MSG_ERROR("An error occured during initializing the clock.\n");
        return (-MESS_ERROR_MISC);

    }
    return ((endTimespec.tv_sec-timer->timer->tv_sec)+(endTimespec.tv_nsec-timer->timer->tv_nsec)*1e-9);
#else
    double t1 = mess_wtime();
    return  ( t1 - timer->dtimer);
#endif
}


/**
 * @brief Get elapsed time since starting the timer and reset.
 * @param[in,out] timer     a mess_timer object
 * @return elapsed time since start ing the  timer
 *
 * The @ref mess_timer_getandreset function returns the elapsed time since the timer was
 * started and restarts the timer.
 *
 */
double mess_timer_getandreset(mess_timer timer){
#ifdef MESS_HAVE_CLOCK_GETTIME
    MSG_FNAME(__func__);
    struct timespec endTimespec;
    if (clock_gettime(CLOCK_MONOTONIC,&endTimespec)!=0){
        MSG_ERROR("An error occured during initializing the clock.\n");
        return (-MESS_ERROR_MISC);
    }
    double result=(endTimespec.tv_sec-timer->timer->tv_sec)+ (endTimespec.tv_nsec-timer->timer->tv_nsec)*1e-9;
    *(timer->timer)=endTimespec;
    return (result);
#else
    double t1 = mess_wtime();
    double res = t1 - timer->dtimer;
    timer->dtimer = t1;
    return (res);
#endif
}


/**
 * @brief Get difference between two timers.
 * @param[in] timerStart input  first mess_timer object
 * @param[in] timerEnd input    second mess_timer object
 * @return difference between timers
 *
 * The @ref mess_timer_diff function returns the difference between two timers.
 *
 */
double mess_timer_diff(mess_timer timerStart, mess_timer timerEnd){
    MSG_FNAME(__func__);
    if ( timerStart == NULL ) {
        MSG_ERROR("timerStart points to NULL\n");
        return (-MESS_ERROR_NULLPOINTER);
    }
    if ( timerEnd == NULL ) {
        MSG_ERROR("timerEnd points to NULL\n");
        return (-MESS_ERROR_NULLPOINTER);
    }
#ifdef MESS_HAVE_CLOCK_GETTIME
    return ((timerEnd->timer->tv_sec-timerStart->timer->tv_sec)+(timerEnd->timer->tv_nsec-timerStart->timer->tv_nsec)*1e-9);
#else
    return (timerEnd->dtimer - timerStart->dtimer);
#endif
}

