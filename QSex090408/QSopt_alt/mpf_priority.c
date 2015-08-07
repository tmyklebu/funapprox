/****************************************************************************/
/*                                                                          */
/* This file is part of QSopt_ex.                                           */
/*                                                                          */
/* (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,       */
/* and Daniel Espinoza.  Sanjeeb Dash's ownership of copyright in           */
/* QSopt_ex is derived from his copyright in QSopt.                         */
/*                                                                          */
/* This code may be used under the terms of the GNU General Public License  */
/* (Version 2.1 or later) as published by the Free Software Foundation.     */
/*                                                                          */
/* Alternatively, use is granted for research purposes only.                */
/*                                                                          */
/* It is your choice of which of these two licenses you are operating       */
/* under.                                                                   */
/*                                                                          */
/* We make no guarantees about the correctness or usefulness of this code.  */
/*                                                                          */
/****************************************************************************/

/* RCSINFO $Id: mpf_priority.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/* */
/* This file is part of CONCORDE                                           */
/* */
/* (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/* Vasek Chvatal, and William Cook                                         */
/* */
/* Permission is granted for academic research use.  For other uses,       */
/* contact the authors for licensing options.                              */
/* */
/* Use at your own risk.  We make no guarantees about the                  */
/* correctness or usefulness of this code.                                 */
/* */
/****************************************************************************/

/****************************************************************************/
/* */
/* PRIORITY QUEUE ROUTINES                             */
/* */
/* */
/* TSP CODE                                    */
/* Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/* Date: March 3, 1997                                                     */
/* March 13, 2002 - Cook Modified for QS)                            */
/* Reference: R.E. Tarjan, Data Structures and Network Algorithms          */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/* int mpf_ILLutil_priority_init (mpf_ILLpriority *pri, int k)                     */
/* -h should point to a mpf_ILLpriority struct.                              */
/* -k an initial allocation for the priority queue.                      */
/* */
/* void mpf_ILLutil_priority_free (mpf_ILLpriority *pri)                           */
/* -frees the spaces allocated for the mpf_ILLpriority queue.                */
/* */
/* void mpf_ILLutil_priority_findmin (mpf_ILLpriority *pri, double *keyval         */
/* void **en)                                                          */
/* -en the entry with least key value (NULL if no entries in mpf_heap).      */
/* -if (keyval != NULL), *keyval will be the minimum key value.          */
/* */
/* int mpf_ILLutil_priority_insert (mpf_ILLpriority *pri, void *data,              */
/* double keyval, int *handle)                                         */
/* -adds (data, keyval) to h.                                            */
/* -handle returns a handle (>= 0) to use when deleting or changing the  */
/* entry                                                                */
/* */
/* void mpf_ILLutil_priority_delete (mpf_ILLpriority *pri, int handle)             */
/* -deletes an entry from the queue.  handle is the value returned by    */
/* mpf_ILLutil_priority_insert.                                             */
/* */
/* void mpf_ILLutil_priority_deletemin (mpf_ILLpriority *pri, double *keyval,      */
/* void **en)                                                         */
/* -like mpf_ILLutil_priority_findmin, but also deletes the entry.           */
/* */
/* void mpf_ILLutil_priority_changekey (mpf_ILLpriority *pri, int handle,          */
/* double newkey)                                                      */
/* -changes the key of an entry in the queue.  handle is the value       */
/* returned by mpf_ILLutil_priority_insert.                                 */
/* */
/****************************************************************************/

/****************************************************************************/
/* */
/* NOTES:                                                                  */
/* These priority queue routines use the mpf_ILLdheap routines to maintain */
/* the priority queue.                                                     */
/* */
/****************************************************************************/

#include "econfig.h"
#include "machdefs.h"
#include "mpf_priority.h"
#include "allocrus.h"
#include "except.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


int mpf_ILLutil_priority_init (mpf_ILLpriority * pri,
      int k)
{
    int i;
    int list;
    int rval = 0;

    pri->space = k;
    ILL_SAFE_MALLOC (pri->pri_info, k, union mpf_ILLpri_data);

    rval = mpf_ILLutil_dheap_init (&pri->mpf_heap, k);
    ILL_CLEANUP_IF (rval);

    list = -1;
    for (i = k - 1; i >= 0; i--) {
	pri->pri_info[i].next = list;
	list = i;
    }
    pri->freelist = list;

CLEANUP:

    if (rval) {
	ILL_IFFREE (pri->pri_info, union mpf_ILLpri_data);
    }
    return rval;
}

void mpf_ILLutil_priority_free (mpf_ILLpriority * pri)
{
    mpf_ILLutil_dheap_free (&pri->mpf_heap);
    ILL_IFFREE (pri->pri_info, union mpf_ILLpri_data);
    pri->space = 0;
}

void mpf_ILLutil_priority_findmin (mpf_ILLpriority * pri,
      mpf_t * keyval,
      void **en)
{
    int handle;

    mpf_ILLutil_dheap_findmin (&pri->mpf_heap, &handle);

    if (handle < 0) {
	*en = (void *) NULL;
    } else {
	if (keyval)
	    mpf_EGlpNumCopy (*keyval, pri->mpf_heap.key[handle]);
	*en = pri->pri_info[handle].data;
    }
}

int mpf_ILLutil_priority_insert (mpf_ILLpriority * pri,
      void *data,
      mpf_t * keyval,
      int *handle)
{
    int newsize;
    int i;
    int list;
    int rval = 0;

    if (pri->freelist == -1) {
	/* Change from 1.3 * pri->space to avoid a warning */
	newsize = pri->space + (pri->space / 3);
	if (newsize < pri->space + 1000)
	    newsize = pri->space + 1000;
	rval = mpf_ILLutil_dheap_resize (&pri->mpf_heap, newsize);
	ILL_CLEANUP_IF (rval);

	pri->pri_info = EGrealloc (pri->pri_info, sizeof (union mpf_ILLpri_data) * newsize);
	/*
			rval = ILLutil_reallocrus_count ((void **) &pri->pri_info, newsize, sizeof (union mpf_ILLpri_data));
			ILL_CLEANUP_IF (rval);
	*/

	list = -1;
	for (i = newsize - 1; i >= pri->space; i--) {
	    pri->pri_info[i].next = list;
	    list = i;
	}
	pri->space = newsize;
	pri->freelist = list;
    }
    i = pri->freelist;
    pri->freelist = pri->pri_info[i].next;
    pri->pri_info[i].data = data;
    mpf_EGlpNumCopy (pri->mpf_heap.key[i], *keyval);
    rval = mpf_ILLutil_dheap_insert (&pri->mpf_heap, i);
    ILL_CLEANUP_IF (rval);

    if (handle)
	*handle = i;

CLEANUP:

    return rval;
}

void mpf_ILLutil_priority_delete (mpf_ILLpriority * pri,
      int handle)
{
    mpf_ILLutil_dheap_delete (&pri->mpf_heap, handle);
    pri->pri_info[handle].next = pri->freelist;
    pri->freelist = handle;
}

void mpf_ILLutil_priority_deletemin (mpf_ILLpriority * pri,
      mpf_t * keyval,
      void **en)
{
    int handle;
    void *data;

    mpf_ILLutil_dheap_deletemin (&pri->mpf_heap, &handle);

    if (handle < 0) {
	*en = (void *) NULL;
    } else {
	if (keyval)
	    mpf_EGlpNumCopy (*keyval, pri->mpf_heap.key[handle]);
	data = pri->pri_info[handle].data;
	pri->pri_info[handle].next = pri->freelist;
	pri->freelist = handle;
	*en = data;
    }
}

void mpf_ILLutil_priority_changekey (mpf_ILLpriority * pri,
      int handle,
      mpf_t * newkey)
{
    mpf_ILLutil_dheap_changekey (&pri->mpf_heap, handle, newkey);
}
