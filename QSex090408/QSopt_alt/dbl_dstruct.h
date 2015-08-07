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

/* RCSINFO $Id: dbl_dstruct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*                       dbl_svector.h                                      */
/*                                                                          */
/****************************************************************************/

#ifndef dbl___SVECTOR_H
#define dbl___SVECTOR_H
#include "econfig.h"

typedef struct dbl_svector
{
	int nzcnt;
	int *indx;
	int size;
	double *coef;
}
dbl_svector;

void dbl_ILLsvector_init (dbl_svector * s),
  dbl_ILLsvector_free (dbl_svector * s);

int dbl_ILLsvector_alloc (dbl_svector * s,
											int nzcnt),
  dbl_ILLsvector_copy (const dbl_svector * s_in,
									 dbl_svector * s_out);

#endif /* dbl___SVECTOR_H */

/****************************************************************************/
/*                                                                          */
/*                           dbl_heap.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef dbl___HEAP_H
#define dbl___HEAP_H

typedef struct
{
	int *entry;
	int *loc;
	double *key;
	int hexist;
	int maxsize;
	int size;
}
dbl_heap;

void dbl_ILLheap_insert (dbl_heap * const h,
										 int const ix),
  dbl_ILLheap_modify (dbl_heap * const h,
									int const ix),
  dbl_ILLheap_delete (dbl_heap * const h,
									int const ix),
  dbl_ILLheap_init (dbl_heap * const h),
  dbl_ILLheap_free (dbl_heap * const h);

int dbl_ILLheap_findmin (dbl_heap * const h),
  dbl_ILLheap_build (dbl_heap * const h,
								 int const nelems,
								 double * key);

#endif /* dbl___HEAP_H */

/****************************************************************************/
/*                                                                          */
/*                         matrix.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef dbl___MATRIX_H
#define dbl___MATRIX_H

typedef struct dbl_ILLmatrix
{
	double *matval;						/* The coefficients.                       */
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
dbl_ILLmatrix;

void dbl_ILLmatrix_init (dbl_ILLmatrix * A);
void dbl_ILLmatrix_free (dbl_ILLmatrix * A);
void dbl_ILLmatrix_prt (FILE * fd,
										dbl_ILLmatrix * A);

#endif /* dbl___MATRIX_H */
