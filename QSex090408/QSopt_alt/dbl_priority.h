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

#ifndef dbl___PRIORITY_H__
#define dbl___PRIORITY_H__
#include "dbl_dheaps_i.h"
/****************************************************************************/
/*                                                                          */
/*                             dbl_priority.c                               */
/*                                                                          */
/****************************************************************************/

typedef struct dbl_ILLpriority
{
	dbl_ILLdheap dbl_heap;
	union dbl_ILLpri_data
	{
		void *data;
		int next;
	}
	 *pri_info;
	int space;
	int freelist;
}
dbl_ILLpriority;

void dbl_ILLutil_priority_free (dbl_ILLpriority * pri),
  dbl_ILLutil_priority_delete (dbl_ILLpriority * pri,
													 int handle),
  dbl_ILLutil_priority_changekey (dbl_ILLpriority * pri,
															int handle,
															double * newkey),
  dbl_ILLutil_priority_findmin (dbl_ILLpriority * pri,
														double * keyval,
														void **en),
  dbl_ILLutil_priority_deletemin (dbl_ILLpriority * pri,
															double * keyval,
															void **en);

int dbl_ILLutil_priority_init (dbl_ILLpriority * pri,
													 int k),
  dbl_ILLutil_priority_insert (dbl_ILLpriority * pri,
													 void *data,
													 double * keyval,
													 int *handle);



#endif
