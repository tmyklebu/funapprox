/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */ 
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

#ifndef mpf___PRIORITY_H__
#define mpf___PRIORITY_H__
#include "mpf_dheaps_i.h"
/****************************************************************************/
/*                                                                          */
/*                             mpf_priority.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct mpf_ILLpriority
{
	mpf_ILLdheap mpf_heap;
	union mpf_ILLpri_data
	{
		void *data;
		int next;
	}
	 *pri_info;
	int space;
	int freelist;
}
mpf_ILLpriority;

void mpf_ILLutil_priority_free (mpf_ILLpriority * pri),
  mpf_ILLutil_priority_delete (mpf_ILLpriority * pri,
													 int handle),
  mpf_ILLutil_priority_changekey (mpf_ILLpriority * pri,
															int handle,
															mpf_t * newkey),
  mpf_ILLutil_priority_findmin (mpf_ILLpriority * pri,
														mpf_t * keyval,
														void **en),
  mpf_ILLutil_priority_deletemin (mpf_ILLpriority * pri,
															mpf_t * keyval,
															void **en);

int mpf_ILLutil_priority_init (mpf_ILLpriority * pri,
													 int k),
  mpf_ILLutil_priority_insert (mpf_ILLpriority * pri,
													 void *data,
													 mpf_t * keyval,
													 int *handle);



#endif
