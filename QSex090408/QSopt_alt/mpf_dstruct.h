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

/* RCSINFO $Id: mpf_dstruct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*                           mpf_svector.h                                      */
/*                                                                          */
/****************************************************************************/

#ifndef mpf___SVECTOR_H
#define mpf___SVECTOR_H
#include "econfig.h"

typedef struct mpf_svector
{
	int nzcnt;
	int *indx;
	int size;
	mpf_t *coef;
}
mpf_svector;

void mpf_ILLsvector_init (mpf_svector * s),
  mpf_ILLsvector_free (mpf_svector * s);

int mpf_ILLsvector_alloc (mpf_svector * s,
											int nzcnt),
  mpf_ILLsvector_copy (const mpf_svector * s_in,
									 mpf_svector * s_out);

#endif /* mpf___SVECTOR_H */

/****************************************************************************/
/*                                                                          */
/*                           mpf_heap.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef mpf___HEAP_H
#define mpf___HEAP_H

typedef struct
{
	int *entry;
	int *loc;
	mpf_t *key;
	int hexist;
	int maxsize;
	int size;
}
mpf_heap;

void mpf_ILLheap_insert (mpf_heap * const h,
										 int const ix),
  mpf_ILLheap_modify (mpf_heap * const h,
									int const ix),
  mpf_ILLheap_delete (mpf_heap * const h,
									int const ix),
  mpf_ILLheap_init (mpf_heap * const h),
  mpf_ILLheap_free (mpf_heap * const h);

int mpf_ILLheap_findmin (mpf_heap * const h),
  mpf_ILLheap_build (mpf_heap * const h,
								 int const nelems,
								 mpf_t * key);

#endif /* mpf___HEAP_H */

/****************************************************************************/
/*                                                                          */
/*                         matrix.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef mpf___MATRIX_H
#define mpf___MATRIX_H

typedef struct mpf_ILLmatrix
{
	mpf_t *matval;						/* The coefficients.                       */
	int *matcnt;									/* Number of coefs in each col.            */
	int *matind;									/* The row indices of the coefs.           */
	int *matbeg;									/* The start of each col.                  */
	int matcols;									/* Number of columns.                      */
	int matrows;
	int matcolsize;								/* Length of matbeg and matcnt.            */
	int matsize;									/* Length of matind and matval.            */
	int matfree;									/* Free space at end of matind.            */
	/* Note: free elements marked by -1 in     */
	/* matind; we keep at least 1 free at end. */
}
mpf_ILLmatrix;

void mpf_ILLmatrix_init (mpf_ILLmatrix * A);
void mpf_ILLmatrix_free (mpf_ILLmatrix * A);
void mpf_ILLmatrix_prt (FILE * fd,
										mpf_ILLmatrix * A);

#endif /* mpf___MATRIX_H */
