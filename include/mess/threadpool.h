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
 * @file include/mess/threadpool.h
 * @brief Interface to handle a pool of threads.
 * @author @koehlerm
 */

#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif
    /** @addtogroup misc_thpool
     * @{ */
    /**
     * @brief Representation of a thread pool job.
     *
     * The @ref mess_threadpooljob_st structure contains all information to represent a job for the threadpool.
     */
    struct mess_threadpooljob_st {
        int exit;                               /**< Value to exit a job if it is done \f$ (-1) \f$*/
        mess_int_t id;                          /**< Identification number of the job */
        void *data;                             /**< Pointer to data of the job */
        void (*worker)(void *data);             /**< Pointer to the worker function. */
        struct mess_threadpooljob_st *next;     /**< Pointer to next job */
    };
    /**
     * @brief Type definition of a thread pool job.
     *
     * The @ref mess_threadpooljob type definition contains all information to represent a job for the threadpool.
     */
    typedef struct mess_threadpooljob_st *mess_threadpooljob;

    /**
     * @brief Representation of a thread pool.
     *
     * The @ref mess_threadpool_st structure contains all information to represent a threadpool.
     */
    struct mess_threadpool_st {
        int num_threads;                                /**< Number of created threads */
        int max_size;                                   /**< Maximum size of the queue holding indices to load */
        int cur_size;                                   /**< Number of current elements in the queue holding the indices to load*/
        mess_int_t nextid;                              /**< Identification number of the next job */
        pthread_t *threads;                             /**< Pointer to thread */
        mess_threadpooljob head, tail;                  /**< First and last job */
        pthread_mutex_t lock;                           /**< Mutex ensuring exclusive access to buffer */
        pthread_cond_t not_empty, not_full,job_done;    /**< Condition variable if a thread is not empty, not full or the job is done*/
        void *hashtable_done;                           /**< Pointer to a hashtable including done jobs.*/
    };

    /**
     * @brief Type definition of a thread pool.
     *
     * The @ref mess_threadpool type definition contains all information to represent a threadpool.
     */
    typedef struct mess_threadpool_st*  mess_threadpool;


    int mess_threadpool_init(mess_threadpool *pool, int threads, int max_size);
    void * mess_threadpool_worker(void *arg);
    int mess_threadpool_clear(mess_threadpool *pool);
    int mess_threadpool_insert(mess_threadpool pool, mess_threadpooljob job, mess_int_t *jobid);
    int mess_threadpool_isdone(mess_threadpool pool, mess_int_t id);
    int mess_threadpool_waitdone(mess_threadpool pool, mess_int_t id);
    int mess_threadpooljob_init(mess_threadpooljob *job);
    int mess_threadpooljob_clear(mess_threadpooljob *job);

    /** @}  */

#ifdef __cplusplus
}
#endif

#endif /* THREADPOOL_H_ */

/* @} */
