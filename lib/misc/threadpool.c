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
 * @file lib/misc/threadpool.c
 * @brief Basic functions for a loader pool.
 * @author @koehlerm
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <complex.h>
#include "mess/config.h"
#include "mess/mess.h"
#include "mess/error_macro.h"
#include "cscutils/ds.h"
#include <string.h>

// #define MESS_DEBUG_THREADPOOL

struct th_done_entry {
    int done;
    char *id;
    void *retdata;
};

static char * th_ht_getname(csc_ds_object_t obj) {
    struct th_done_entry *entry = (struct th_done_entry*) obj;
    return entry->id;
}

static void th_ht_free(csc_ds_object_t obj) {
    struct th_done_entry *entry = (struct th_done_entry*) obj;
    mess_free(entry->id);
    mess_free(entry);
}

static size_t hash (const char *name, size_t size){
    size_t sum = 0;
    size_t i = 0;
    size_t len = strlen(name);
    for (i = 0; i < len; i++) {
        sum += (size_t)name[i];
    }
    return sum % size;
}


/**
 * @brief Create a loader pool.
 * @param[out] pool         pointer to a threadpool object
 * @param[in] threads input     number of threads to create
 * @param[in] max_size input    maximal size of the queue
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_threadpool_init function initializes load threads
 * of the cache to disk property of a spmv. \n
 * Parameter \f$ threads \f$ describes the amount of loader threads which should be created, usually it should be one. \n
 * The size parameter \f$ max_{size} \f$ sets the maximum length of the queue which holds the indices to load, a good
 * size is a number of shifts but it depends on the usage of the shift parameters in the algorithm.
 *
 */
int mess_threadpool_init( mess_threadpool *pool,  int threads, int max_size) {
    MSG_FNAME(__func__);
    int i;
    pthread_attr_t tattr;
    mess_try_alloc( *pool ,  struct mess_threadpool_st* , sizeof(struct mess_threadpool_st));
    (*pool)->cur_size=0;
    (*pool)->max_size=max_size;
    (*pool)->num_threads = threads;
    (*pool)->nextid = 0;
    (*pool)->head = (*pool)->tail = NULL;
    (*pool)->hashtable_done = (void*) csc_ds_hashtable(257, th_ht_getname, hash, th_ht_free);
    if ( (*pool)->hashtable_done == NULL ) {
        MSG_ERROR("Failed to create hashtable for the results.\n");
        free(*pool);
        *pool = NULL;
        return MESS_ERROR_NULLPOINTER;
    }

    pthread_mutex_init(&((*pool)->lock),NULL);
    pthread_cond_init(&((*pool)->not_empty),NULL);
    pthread_cond_init(&((*pool)->not_full),NULL);
    pthread_cond_init(&((*pool)->job_done),NULL);
    pthread_attr_init(&tattr);
    pthread_attr_setschedpolicy(&tattr, SCHED_RR);
#ifdef MESS_DEBUG_THREADPOOL
    MSG_INFO("create threadpool with %d worker.\n", threads);
#endif
    mess_try_alloc((*pool)->threads, pthread_t *, sizeof(pthread_t) * threads);
    for ( i = 0; i < threads; i++ ){
        pthread_create(&((*pool)->threads[i]),&tattr, mess_threadpool_worker, (void*) (*pool));
    }

    return 0;
}

/**
 * @brief Implementation of the loader thread.
 * @param[in] arg input argument from the task pool
 *
 * The @ref mess_threadpool_worker function realizes the load algorithm for a
 * index specified in \f$ arg \f$. \n
 * If the values are loaded, a signal is send to a waiting thread.
 *
 */
void * mess_threadpool_worker(void *arg){
    MSG_FNAME(__func__);
    mess_threadpool pool = (mess_threadpool) arg;
    mess_threadpooljob job;
    char idstr[20];
    struct th_done_entry *entry;

    if ( arg == NULL ) {
        MSG_ERROR("no threadpool given.\n");
        pthread_exit(NULL);
        return NULL;
    }

    for (;;){
        pthread_mutex_lock(&(pool->lock));
        while (pool->cur_size == 0){
            pthread_cond_wait(&(pool->not_empty), &(pool->lock));
        }
        job = pool->head;
        pool->cur_size--;
        if ( pool->cur_size == 0) {
            pool->head=pool->tail = NULL;
        } else {
            pool->head = job->next;
        }
        if (pool->cur_size == pool->max_size-1){
            pthread_cond_broadcast(&(pool->not_full));
        }
        pthread_mutex_unlock(&(pool->lock));
        if ( job->exit){
#ifdef MESS_DEBUG_THREADPOOL
            MSG_INFO("exit worker thread.\n");
#endif
            mess_threadpooljob_clear(&job);
            return NULL;
        } else {
            job->worker(job->data);
            pthread_mutex_lock(&(pool->lock));
            snprintf(idstr, 19, "" MESS_PRINTF_INT "", job->id);
#ifdef MESS_DEBUG_THREADPOOL
            MSG_INFO("job, " MESS_PRINTF_INT ", done.\n", job->id);
#endif
            entry = (struct th_done_entry *) csc_ds_find((csc_ds_t *) pool->hashtable_done, idstr);
            if ( entry != NULL){
                entry -> done =1;
            }
            pthread_cond_signal(&(pool->job_done));
            pthread_mutex_unlock (&(pool->lock));
            mess_threadpooljob_clear(&job);
        }
    }
    return NULL;
}

/**
 *  @brief Stop and clean up a loader pool.
 *  @param[in,out] pool     pool to stop
 *  @return zero on success or a non-zero error value otherwise
 *
 *  The @ref mess_threadpool_clear function sends a stop job to all threads and
 *  cleans up all data structures.
 *
 */
int mess_threadpool_clear(mess_threadpool *pool){
    MSG_FNAME(__func__);
    mess_threadpooljob tmp;
    if ( pool == NULL || (*pool)== NULL){
        MSG_ERROR("pool points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    int i;
    for (i=0;i<(*pool)->num_threads;i++){
        mess_threadpooljob job;
        mess_threadpooljob_init(&job);
        job->exit = -1;
        mess_threadpool_insert(*pool, job,NULL);
    }
    for (i=0;i<(*pool)->num_threads;i++){
        pthread_join((*pool)->threads[i], NULL);
    }
    mess_free((*pool)->threads);
    while ((*pool)->head != NULL){
        tmp = (*pool)->head;
        (*pool)->head = (*pool)->head->next;
        mess_threadpooljob_clear(&tmp);
    }
    csc_ds_free((csc_ds_t*)(*pool)->hashtable_done);
    mess_free(*pool);
    *pool=NULL;
    return 0;
}

/**
 * @brief Insert a task in a loader pool.
 * @param[in,out] pool  threadpool object where the job should be insert
 * @param[in] job input     job object to put in the pool
 * @param[out] jobid    output id of the job
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_threadpool_insert function inserts a \f$ job \f$ into a loader pool \f$ pool \f$ getting a identification
 * \f$ jobid \f$.
 */
int mess_threadpool_insert(mess_threadpool pool, mess_threadpooljob job, mess_int_t *jobid){
    MSG_FNAME(__func__);

    csc_ds_t *ht = (csc_ds_t*) pool->hashtable_done;
    struct th_done_entry *entry;

    if ( pool == NULL) {
        MSG_ERROR("pool points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    if ( job == NULL) {
        MSG_ERROR("job points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    pthread_mutex_lock(&(pool->lock));
    while(pool->cur_size == pool->max_size){
        pthread_cond_wait(&(pool->not_full),&(pool->lock));
    }
    job->id = pool->nextid;
#ifdef MESS_DEBUG_THREADPOOL
    MSG_INFO("create job, id = " MESS_PRINTF_INT "\n", job->id);
#endif
    pool->nextid++;
    job->next = NULL;
    // insert in done table
    mess_try_alloc(entry, struct th_done_entry *, sizeof(struct th_done_entry));
    mess_try_alloc(entry->id, char *, sizeof(char)*20);
    memset(entry->id,0, 20);
    snprintf(entry->id,19, "" MESS_PRINTF_INT "",job->id);
    entry->retdata = NULL;
    entry->done = 0;
    csc_ds_insert(ht, (void*) entry);

    if (pool->cur_size== 0) {
        pool->head = pool->tail = job;
        pthread_cond_signal(&(pool->not_empty));
    } else {
        pool->tail->next = job;
        pool->tail = job;
        pthread_cond_signal(&(pool->not_empty));
    }
    pool->cur_size ++;
    pthread_mutex_unlock(&(pool->lock));
    if ( jobid != NULL){
        *jobid=job->id;
    }
    return 0;
}

/**
 * @brief Create a job for a loader pool.
 * @param[out] job  job object to create
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_threadpooljob_init function creates a job for a loader pool.
 *
 */
int mess_threadpooljob_init ( mess_threadpooljob *job) {
    MSG_FNAME(__func__);
    mess_try_alloc( *job, struct mess_threadpooljob_st *, sizeof(struct mess_threadpooljob_st));
    (*job)->exit = 0;
    (*job)->data = NULL;
    (*job)->worker = NULL;
    return 0;
}

/**
 * @brief Stop and clean up a loader pool job.
 * @param[in,out] job   job object to stop
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_threadpooljob_clear function stops and cleans up a job for a loader pool.
 */
int mess_threadpooljob_clear ( mess_threadpooljob *job){
    MSG_FNAME(__func__);
    if ( job == NULL || (*job) == NULL){
        MSG_ERROR("job points to NULL\n");
        return MESS_ERROR_NULLPOINTER;
    }
    mess_free(*job);
    return 0;
}

/**
 * @brief Check if a job is done.
 * @param[in] pool input loader pool
 * @param[in] id input  id to check
 * @return 1 if a job is not done, otherwise the job is done
 *
 * The @ref mess_threadpool_isdone function checks if a job
 * with identification \f$ id \f$ is done.
 *
 */
int  mess_threadpool_isdone ( mess_threadpool pool, mess_int_t id )
{
    MSG_FNAME(__func__);
    int ret  = 0;
    char idstr[20];
    struct th_done_entry *entry;
    mess_check_nullpointer(pool);
    snprintf(idstr, 19, "" MESS_PRINTF_INT "", id);
#ifdef MESS_DEBUG_THREADPOOL
    MSG_INFO("check if job, " MESS_PRINTF_INT ", is done.\n", id);
#endif
    pthread_mutex_lock (&(pool->lock));
    entry = (struct th_done_entry *) csc_ds_find((csc_ds_t*) pool->hashtable_done, idstr);
    if ( entry != NULL){
        ret = entry -> done;
    } else {
        ret =1;
    }
    pthread_mutex_unlock(&(pool->lock));
    return ret;
}       /* -----  end of function mess_threadpool_isdone  ----- */


/**
 * @brief Wait until a job is done.
 * @param[in] pool input loader pool
 * @param[in] id input id to wait
 * @return zero on success or a non-zero error value otherwise
 *
 * The @ref mess_threadpool_waitdone function waits until a job with identification \f$ id \f$ is done.
 *
 */
int mess_threadpool_waitdone ( mess_threadpool pool, mess_int_t id)
{
    MSG_FNAME(__func__);
    char idstr[20];
    struct th_done_entry *entry;

    mess_check_nullpointer(pool);
#ifdef MESS_DEBUG_THREADPOOL
    MSG_INFO("wait until " MESS_PRINTF_INT " is done.\n", id);
#endif
    snprintf(idstr, 19, "" MESS_PRINTF_INT "", id);

    pthread_mutex_lock(&(pool->lock));
    entry = (struct th_done_entry *) csc_ds_find((csc_ds_t* ) pool->hashtable_done, idstr);
    if ( entry == NULL){
        // MSG_ERROR("job " MESS_PRINTF_INT " not found.\n", id);
        pthread_mutex_unlock(&(pool->lock));
        // return MESS_ERROR_NULLPOINTER;
        return 0;
    }
    while (entry->done !=1 ){
        pthread_cond_wait(&(pool->job_done), &(pool->lock));
        entry = (struct th_done_entry *) csc_ds_find((csc_ds_t*) pool->hashtable_done, idstr);
    }
    pthread_mutex_unlock(&(pool->lock));
    return 0;
}       /* -----  end of function mess_threadpool_waitdone  ----- */

