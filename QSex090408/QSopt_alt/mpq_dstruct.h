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

/* RCSINFO $Id: mpq_dstruct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*                           mpq_svector.h                                      */
/*                                                                          */
/****************************************************************************/

#ifndef mpq___SVECTOR_H
#define mpq___SVECTOR_H
#include "econfig.h"

typedef struct mpq_svector
{
	int nzcnt;
	int *indx;
	int size;
	mpq_t *coef;
}
mpq_svector;

void mpq_ILLsvector_init (mpq_svector * s),
  mpq_ILLsvector_free (mpq_svector * s);

int mpq_ILLsvector_alloc (mpq_svector * s,
											int nzcnt),
  mpq_ILLsvector_copy (const mpq_svector * s_in,
									 mpq_svector * s_out);

#endif /* mpq___SVECTOR_H */

/****************************************************************************/
/*                                                                          */
/*                           mpq_heap.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef mpq___HEAP_H
#define mpq___HEAP_H

typedef struct
{
	int *entry;
	int *loc;
	mpq_t *key;
	int hexist;
	int maxsize;
	int size;
}
mpq_heap;

void mpq_ILLheap_insert (mpq_heap * const h,
										 int const ix),
  mpq_ILLheap_modify (mpq_heap * const h,
									int const ix),
  mpq_ILLheap_delete (mpq_heap * const h,
									int const ix),
  mpq_ILLheap_init (mpq_heap * const h),
  mpq_ILLheap_free (mpq_heap * const h);

int mpq_ILLheap_findmin (mpq_heap * const h),
  mpq_ILLheap_build (mpq_heap * const h,
								 int const nelems,
								 mpq_t * key);

#endif /* mpq___HEAP_H */

/****************************************************************************/
/*                                                                          */
/*                         matrix.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef mpq___MATRIX_H
#define mpq___MATRIX_H

typedef struct mpq_ILLmatrix
{
	mpq_t *matval;						/* The coefficients.                       */
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
mpq_ILLmatrix;

void mpq_ILLmatrix_init (mpq_ILLmatrix * A);
void mpq_ILLmatrix_free (mpq_ILLmatrix * A);
void mpq_ILLmatrix_prt (FILE * fd,
										mpq_ILLmatrix * A);

#endif /* mpq___MATRIX_H */
