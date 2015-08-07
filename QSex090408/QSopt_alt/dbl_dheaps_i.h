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

#ifndef dbl___DHEAPS_I_H__
#define dbl___DHEAPS_I_H__
/****************************************************************************/
/*                                                                          */
/*                             dbl_dheaps_i.c                               */
/*                                                                          */
/****************************************************************************/

typedef struct dbl_ILLdheap
{
	double *key;
	int *entry;
	int *loc;
	int total_space;
	int size;
}
dbl_ILLdheap;

void dbl_ILLutil_dheap_free (dbl_ILLdheap * h),
  dbl_ILLutil_dheap_delete (dbl_ILLdheap * h,
												int i),
  dbl_ILLutil_dheap_changekey (dbl_ILLdheap * h,
													 int i,
													 double * newkey),
  dbl_ILLutil_dheap_findmin (dbl_ILLdheap * h,
												 int *i),
  dbl_ILLutil_dheap_deletemin (dbl_ILLdheap * h,
													 int *i);

int dbl_ILLutil_dheap_init (dbl_ILLdheap * h,
												int k),
  dbl_ILLutil_dheap_resize (dbl_ILLdheap * h,
												int newsize),
  dbl_ILLutil_dheap_insert (dbl_ILLdheap * h,
												int i);



#endif
