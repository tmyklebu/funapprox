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

#ifndef mpq___DHEAPS_I_H__
#define mpq___DHEAPS_I_H__
/****************************************************************************/
/*                                                                          */
/*                             mpq_dheaps_i.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct mpq_ILLdheap
{
	mpq_t *key;
	int *entry;
	int *loc;
	int total_space;
	int size;
}
mpq_ILLdheap;

void mpq_ILLutil_dheap_free (mpq_ILLdheap * h),
  mpq_ILLutil_dheap_delete (mpq_ILLdheap * h,
												int i),
  mpq_ILLutil_dheap_changekey (mpq_ILLdheap * h,
													 int i,
													 mpq_t * newkey),
  mpq_ILLutil_dheap_findmin (mpq_ILLdheap * h,
												 int *i),
  mpq_ILLutil_dheap_deletemin (mpq_ILLdheap * h,
													 int *i);

int mpq_ILLutil_dheap_init (mpq_ILLdheap * h,
												int k),
  mpq_ILLutil_dheap_resize (mpq_ILLdheap * h,
												int newsize),
  mpq_ILLutil_dheap_insert (mpq_ILLdheap * h,
												int i);



#endif
