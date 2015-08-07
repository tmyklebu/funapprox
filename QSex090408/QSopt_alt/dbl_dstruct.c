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

/* RCS_INFO = "$RCSfile: dbl_dstruct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#include "config.h"
#include "dbl_iqsutil.h"
#include "dbl_dstruct.h"
#include "dbl_qsopt.h"
#include "dbl_lpdefs.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/****************************************************************************/
/* */
/* dbl_svector                                       */
/* */
/* Written by:  Applegate, Cook, Dash                                      */
/* Date:                                                                   */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/****************************************************************************/

void dbl_ILLsvector_init (dbl_svector * s)
{
    s->nzcnt = 0;
    s->indx = 0;
    s->coef = 0;
}

void dbl_ILLsvector_free (dbl_svector * s)
{
    ILL_IFFREE (s->indx, int);
    dbl_EGlpNumFreeArray (s->coef);
    s->nzcnt = 0;
}

int dbl_ILLsvector_alloc (dbl_svector * s,
      int nzcnt)
{
    int rval = 0;

    s->nzcnt = nzcnt;
    if (nzcnt == 0) {
	s->indx = 0;
	s->coef = 0;
    } else {
	ILL_SAFE_MALLOC (s->indx, nzcnt, int);
	s->coef = dbl_EGlpNumAllocArray (nzcnt);
    }
    return 0;
CLEANUP:
    ILL_IFFREE (s->indx, int);
    dbl_EGlpNumFreeArray (s->coef);
    ILL_RETURN (rval, "dbl_ILLsvector_alloc");
}

int dbl_ILLsvector_copy (const dbl_svector * s_in,
      dbl_svector * s_out)
{
    int i;
    int nzcnt = s_in->nzcnt;
    int rval = 0;

    rval = dbl_ILLsvector_alloc (s_out, nzcnt);
    ILL_CLEANUP_IF (rval);
    for (i = 0; i < nzcnt; i++) {
	s_out->indx[i] = s_in->indx[i];
	dbl_EGlpNumCopy (s_out->coef[i], s_in->coef[i]);
    }

CLEANUP:
    ILL_RETURN (rval, "dbl_ILLsvector_copy");
}

/****************************************************************************/
/* */
/* dbl_heap                                          */
/* */
/* Written by:  Applegate, Cook, Dash                                      */
/* Date:                                                                   */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/****************************************************************************/

#define dbl_DEBUG_HEAP 0

#define dbl_HEAP_D 3
#define dbl_HEAP_UP(x) (((x)-1)/dbl_HEAP_D)
#define dbl_HEAP_DOWN(x) (((x)*dbl_HEAP_D)+1)

static int dbl_siftup (dbl_heap * h,
      int hloc,
      int ix),
  dbl_siftdown (dbl_heap * h,
      int hloc,
      int ix),
  dbl_maxchild (dbl_heap * h,
      int hloc);

static int dbl_siftup (dbl_heap * h,
      int hloc,
      int ix)
{
    int i = hloc;
    int p = dbl_HEAP_UP (i);
    double val;
    dbl_EGlpNumInitVar (val);
    dbl_EGlpNumCopy (val, h->key[ix]);

    while (i > 0 && dbl_EGlpNumIsLess (h->key[h->entry[p]], val)) {
	h->entry[i] = h->entry[p];
	h->loc[h->entry[i]] = i;
	i = p;
	p = dbl_HEAP_UP (p);
    }
    h->entry[i] = ix;
    h->loc[ix] = i;
    ILL_IFTRACE2 ("%s:%la:%d:%d:%d\n", __func__, dbl_EGlpNumToLf (val), hloc, ix, i);
    dbl_EGlpNumClearVar (val);
    return i;
}

static int dbl_siftdown (dbl_heap * h,
      int hloc,
      int ix)
{
    int i = hloc;
    int c = dbl_maxchild (h, i);
    double val;
    dbl_EGlpNumInitVar (val);
    dbl_EGlpNumCopy (val, h->key[ix]);
    ILL_IFTRACE2 ("%s:%d:%d:%d:%la", __func__, hloc, ix, c, dbl_EGlpNumToLf (val));

    while (c != -1 && dbl_EGlpNumIsLess (val, h->key[h->entry[c]])) {
	h->entry[i] = h->entry[c];
	h->loc[h->entry[i]] = i;
	i = c;
	c = dbl_maxchild (h, c);
    }
    h->entry[i] = ix;
    h->loc[ix] = i;
    dbl_EGlpNumClearVar (val);
    ILL_IFTRACE2 ("%s:%d:%d\n", __func__, ix, i);
    return i;
}

/* extern double dbl_ILL_MINDOUBLE; */
static int dbl_maxchild (dbl_heap * h,
      int hloc)
{
    int i;
    int mc = -1;
    int hmin = dbl_HEAP_D * hloc + 1;
    int hmax = dbl_HEAP_D * hloc + dbl_HEAP_D;
    double val;
    dbl_EGlpNumInitVar (val);
    dbl_EGlpNumCopy (val, dbl_ILL_MINDOUBLE);
    ILL_IFTRACE2 (" %s:%d", __func__, hloc);

    for (i = hmin; i <= hmax && i < h->size; i++)
	if (dbl_EGlpNumIsLess (val, h->key[h->entry[i]])) {
	    dbl_EGlpNumCopy (val, h->key[h->entry[i]]);
	    mc = i;
	    ILL_IFTRACE2 (":%d:%la", mc, dbl_EGlpNumToLf (val));
	}
    dbl_EGlpNumClearVar (val);
    ILL_IFTRACE2 ("\n");
    return mc;
}

#if dbl_DEBUG_HEAP > 0

static void dbl_printheap (dbl_heap * h)
{
    int i;

    printf ("entry (%d): ", h->size);
    for (i = 0; i < h->size; i++)
	printf ("%d ", h->entry[i]);
    printf ("\n loc: ");
    for (i = 0; i < h->maxsize; i++)
	printf ("%d ", h->loc[i]);
    printf ("\n key: ");
    for (i = 0; i < h->maxsize; i++)
	printf ("%la ", dbl_EGlpNumToLf (h->key[i]));
    printf ("\n key(sorted): ");
    for (i = 0; i < h->size; i++)
	printf ("%la ", dbl_EGlpNumToLf (h->key[h->entry[i]]));
    printf ("\n");
}

static void dbl_heapcheck (dbl_heap * h)
{
    int i, tcnt = 0;

    for (i = 0; i < h->maxsize; i++) {
	if (h->loc[i] < -1)
	    printf ("error in dbl_heap\n");
	else if (h->loc[i] > -1)
	    tcnt++;
    }
    if (tcnt != h->size)
	printf ("error 3 in dbl_heap\n");

    for (i = 0; i < h->size; i++) {
	if (h->loc[h->entry[i]] != i)
	    printf ("error 1 in dbl_heap\n");
	if (dbl_EGlpNumIsEqqual (h->key[h->entry[i]], dbl_zeroLpNum))
	    printf ("error 2 in dbl_heap\n");
	if (dbl_EGlpNumIsLess (h->key[h->entry[dbl_HEAP_UP (i)]], h->key[h->entry[i]]))
	    printf ("error 4 in dbl_heap\n");
    }
}

#endif

void dbl_ILLheap_insert (dbl_heap * const h,
      int const ix)
{
    int i = h->size;
    ILL_IFTRACE ("%s:%d:%la\n", __func__, ix, dbl_EGlpNumToLf (h->key[ix]));

    i = dbl_siftup (h, i, ix);
    h->size++;

#if dbl_DEBUG_HEAP > 0
    dbl_heapcheck (h);
#endif
#if dbl_DEBUG_HEAP > 1
    dbl_printheap (h);
#endif
}

void dbl_ILLheap_modify (dbl_heap * const h,
      int const ix)
{
    int i = h->loc[ix];
    int pi = i;
    ILL_IFTRACE ("%s:%d\n", __func__, ix);

    if (h->loc[ix] == -1)
	return;
    i = dbl_siftup (h, i, ix);
    if (pi == i)
	i = dbl_siftdown (h, i, ix);

#if dbl_DEBUG_HEAP > 0
    dbl_heapcheck (h);
#endif
#if dbl_DEBUG_HEAP > 1
    dbl_printheap (h);
#endif
}

void dbl_ILLheap_delete (dbl_heap * const h,
      int const ix)
{
    int i = h->loc[ix];
    int pi = i;
    int nix = h->entry[h->size - 1];
    ILL_IFTRACE ("%s:%d:%d:%d\n", __func__, ix, nix, pi);

    h->loc[ix] = -1;
    h->size--;
    if (nix == ix) {
#if dbl_DEBUG_HEAP > 0
	dbl_heapcheck (h);
#endif
#if dbl_DEBUG_HEAP > 1
	dbl_printheap (h);
#endif
	return;
    }
    h->entry[i] = nix;
    h->loc[nix] = i;

    i = dbl_siftup (h, i, nix);
    ILL_IFTRACE ("%s:%d:%d:%d:%d\n", __func__, ix, nix, pi, i);
    if (pi == i)
	dbl_siftdown (h, i, nix);

#if dbl_DEBUG_HEAP > 0
    dbl_heapcheck (h);
#endif
#if dbl_DEBUG_HEAP > 1
    dbl_printheap (h);
#endif
}

int dbl_ILLheap_findmin (dbl_heap * const h)
{
    if (h->hexist == 0 || h->size <= 0)
	return -1;
    return h->entry[0];
}

void dbl_ILLheap_init (dbl_heap * const h)
{
    h->entry = NULL;
    h->loc = NULL;
    h->key = NULL;
    h->hexist = 0;
}

int dbl_ILLheap_build (dbl_heap * const h,
      int const nelems,
      double *key)
{
    int rval = 0;
    int i, n = 0;
    ILL_IFTRACE ("%s:%d\n", __func__, nelems);

    h->hexist = 1;
    h->size = 0;
    h->maxsize = nelems;
    h->key = key;
    ILL_SAFE_MALLOC (h->entry, nelems, int);
    ILL_SAFE_MALLOC (h->loc, nelems, int);

    for (i = 0; i < nelems; i++) {
	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, key[i])) {
	    h->entry[n] = i;
	    h->loc[i] = n;
	    n++;
	} else
	    h->loc[i] = -1;
    }
    h->size = n;
    for (i = n - 1; i >= 0; i--) {
	ILL_IFTRACE2 ("insert %la\n", dbl_EGlpNumToLf (h->key[h->entry[i]]));
	dbl_siftdown (h, i, h->entry[i]);
    }

#if dbl_DEBUG_HEAP > 0
    dbl_heapcheck (h);
#endif
#if dbl_DEBUG_HEAP > 1
    dbl_printheap (h);
#endif

CLEANUP:
    if (rval)
	dbl_ILLheap_free (h);
    ILL_RETURN (rval, "dbl_ILLheap_init");
}

void dbl_ILLheap_free (dbl_heap * const h)
{
    if (h->hexist) {
	ILL_IFFREE (h->entry, int);
	ILL_IFFREE (h->loc, int);
	h->hexist = 0;
	h->maxsize = 0;
	h->size = 0;
    }
}


/****************************************************************************/
/* */
/* matrix                                          */
/* */
/* Written by:  Applegate, Cook, Dash                                      */
/* Date:                                                                   */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/****************************************************************************/

void dbl_ILLmatrix_init (dbl_ILLmatrix * A)
{
    if (A) {
	A->matval = 0;
	A->matcnt = 0;
	A->matbeg = 0;
	A->matind = 0;
	A->matcols = 0;
	A->matcolsize = 0;
	A->matrows = 0;
	A->matsize = 0;
	A->matfree = 0;
    }
}

void dbl_ILLmatrix_free (dbl_ILLmatrix * A)
{
    if (A) {
	dbl_EGlpNumFreeArray (A->matval);
	ILL_IFFREE (A->matcnt, int);
	ILL_IFFREE (A->matbeg, int);
	ILL_IFFREE (A->matind, int);
	dbl_ILLmatrix_init (A);
    }
}

void dbl_ILLmatrix_prt (FILE * fd,
      dbl_ILLmatrix * A)
{
    int j, k;
    if (A == NULL) {
	fprintf (fd, "Matrix %p: empty\n", (void *) A);
    } else {
	fprintf (fd, "Matrix %p: nrows = %d ncols = %d\n",
	    (void *) A, A->matrows, A->matcols);
	for (j = 0; j < A->matcols; j++) {
	    fprintf (fd, "col %d: ", j);
	    for (k = A->matbeg[j]; k < A->matbeg[j] + A->matcnt[j]; k++) {
		fprintf (fd, "row %d=%.3f ", A->matind[k], dbl_EGlpNumToLf (A->matval[k]));
	    }
	    fprintf (fd, "\n");
	}
    }
}
