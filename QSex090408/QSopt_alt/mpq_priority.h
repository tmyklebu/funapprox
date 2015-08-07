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

#ifndef mpq___PRIORITY_H__
#define mpq___PRIORITY_H__
#include "mpq_dheaps_i.h"
/****************************************************************************/
/*                                                                          */
/*                             mpq_priority.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct mpq_ILLpriority
{
	mpq_ILLdheap mpq_heap;
	union mpq_ILLpri_data
	{
		void *data;
		int next;
	}
	 *pri_info;
	int space;
	int freelist;
}
mpq_ILLpriority;

void mpq_ILLutil_priority_free (mpq_ILLpriority * pri),
  mpq_ILLutil_priority_delete (mpq_ILLpriority * pri,
													 int handle),
  mpq_ILLutil_priority_changekey (mpq_ILLpriority * pri,
															int handle,
															mpq_t * newkey),
  mpq_ILLutil_priority_findmin (mpq_ILLpriority * pri,
														mpq_t * keyval,
														void **en),
  mpq_ILLutil_priority_deletemin (mpq_ILLpriority * pri,
															mpq_t * keyval,
															void **en);

int mpq_ILLutil_priority_init (mpq_ILLpriority * pri,
													 int k),
  mpq_ILLutil_priority_insert (mpq_ILLpriority * pri,
													 void *data,
													 mpq_t * keyval,
													 int *handle);



#endif
