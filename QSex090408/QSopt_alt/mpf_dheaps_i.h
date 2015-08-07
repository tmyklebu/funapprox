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

#ifndef mpf___DHEAPS_I_H__
#define mpf___DHEAPS_I_H__
/****************************************************************************/
/*                                                                          */
/*                             mpf_dheaps_i.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct mpf_ILLdheap
{
	mpf_t *key;
	int *entry;
	int *loc;
	int total_space;
	int size;
}
mpf_ILLdheap;

void mpf_ILLutil_dheap_free (mpf_ILLdheap * h),
  mpf_ILLutil_dheap_delete (mpf_ILLdheap * h,
												int i),
  mpf_ILLutil_dheap_changekey (mpf_ILLdheap * h,
													 int i,
													 mpf_t * newkey),
  mpf_ILLutil_dheap_findmin (mpf_ILLdheap * h,
												 int *i),
  mpf_ILLutil_dheap_deletemin (mpf_ILLdheap * h,
													 int *i);

int mpf_ILLutil_dheap_init (mpf_ILLdheap * h,
												int k),
  mpf_ILLutil_dheap_resize (mpf_ILLdheap * h,
												int newsize),
  mpf_ILLutil_dheap_insert (mpf_ILLdheap * h,
												int i);



#endif
