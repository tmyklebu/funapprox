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

/* RCS_INFO = "$RCSfile: mpq_dstruct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#include "config.h"
#include "mpq_iqsutil.h"
#include "mpq_dstruct.h"
#include "mpq_qsopt.h"
#include "mpq_lpdefs.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/****************************************************************************/
/* */
/* mpq_svector                                       */
/* */
/* Written by:  Applegate, Cook, Dash                                      */
/* Date:                                                                   */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/****************************************************************************/

void mpq_ILLsvector_init (mpq_svector * s)
{
    s->nzcnt = 0;
    s->indx = 0;
    s->coef = 0;
}

void mpq_ILLsvector_free (mpq_svector * s)
{
    ILL_IFFREE (s->indx, int);
    mpq_EGlpNumFreeArray (s->coef);
    s->nzcnt = 0;
}

int mpq_ILLsvector_alloc (mpq_svector * s,
      int nzcnt)
{
    int rval = 0;

    s->nzcnt = nzcnt;
    if (nzcnt == 0) {
	s->indx = 0;
	s->coef = 0;
    } else {
	ILL_SAFE_MALLOC (s->indx, nzcnt, int);
	s->coef = mpq_EGlpNumAllocArray (nzcnt);
    }
    return 0;
CLEANUP:
    ILL_IFFREE (s->indx, int);
    mpq_EGlpNumFreeArray (s->coef);
    ILL_RETURN (rval, "mpq_ILLsvector_alloc");
}

int mpq_ILLsvector_copy (const mpq_svector * s_in,
      mpq_svector * s_out)
{
    int i;
    int nzcnt = s_in->nzcnt;
    int rval = 0;

    rval = mpq_ILLsvector_alloc (s_out, nzcnt);
    ILL_CLEANUP_IF (rval);
    for (i = 0; i < nzcnt; i++) {
	s_out->indx[i] = s_in->indx[i];
	mpq_EGlpNumCopy (s_out->coef[i], s_in->coef[i]);
    }

CLEANUP:
    ILL_RETURN (rval, "mpq_ILLsvector_copy");
}

/****************************************************************************/
/* */
/* mpq_heap                                          */
/* */
/* Written by:  Applegate, Cook, Dash                                      */
/* Date:                                                                   */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/****************************************************************************/

#define mpq_DEBUG_HEAP 0

#define mpq_HEAP_D 3
#define mpq_HEAP_UP(x) (((x)-1)/mpq_HEAP_D)
#define mpq_HEAP_DOWN(x) (((x)*mpq_HEAP_D)+1)

static int mpq_siftup (mpq_heap * h,
      int hloc,
      int ix),
  mpq_siftdown (mpq_heap * h,
      int hloc,
      int ix),
  mpq_maxchild (mpq_heap * h,
      int hloc);

static int mpq_siftup (mpq_heap * h,
      int hloc,
      int ix)
{
    int i = hloc;
    int p = mpq_HEAP_UP (i);
    mpq_t val;
    mpq_EGlpNumInitVar (val);
    mpq_EGlpNumCopy (val, h->key[ix]);

    while (i > 0 && mpq_EGlpNumIsLess (h->key[h->entry[p]], val)) {
	h->entry[i] = h->entry[p];
	h->loc[h->entry[i]] = i;
	i = p;
	p = mpq_HEAP_UP (p);
    }
    h->entry[i] = ix;
    h->loc[ix] = i;
    ILL_IFTRACE2 ("%s:%la:%d:%d:%d\n", __func__, mpq_EGlpNumToLf (val), hloc, ix, i);
    mpq_EGlpNumClearVar (val);
    return i;
}

static int mpq_siftdown (mpq_heap * h,
      int hloc,
      int ix)
{
    int i = hloc;
    int c = mpq_maxchild (h, i);
    mpq_t val;
    mpq_EGlpNumInitVar (val);
    mpq_EGlpNumCopy (val, h->key[ix]);
    ILL_IFTRACE2 ("%s:%d:%d:%d:%la", __func__, hloc, ix, c, mpq_EGlpNumToLf (val));

    while (c != -1 && mpq_EGlpNumIsLess (val, h->key[h->entry[c]])) {
	h->entry[i] = h->entry[c];
	h->loc[h->entry[i]] = i;
	i = c;
	c = mpq_maxchild (h, c);
    }
    h->entry[i] = ix;
    h->loc[ix] = i;
    mpq_EGlpNumClearVar (val);
    ILL_IFTRACE2 ("%s:%d:%d\n", __func__, ix, i);
    return i;
}

/* extern mpq_t mpq_ILL_MINDOUBLE; */
static int mpq_maxchild (mpq_heap * h,
      int hloc)
{
    int i;
    int mc = -1;
    int hmin = mpq_HEAP_D * hloc + 1;
    int hmax = mpq_HEAP_D * hloc + mpq_HEAP_D;
    mpq_t val;
    mpq_EGlpNumInitVar (val);
    mpq_EGlpNumCopy (val, mpq_ILL_MINDOUBLE);
    ILL_IFTRACE2 (" %s:%d", __func__, hloc);

    for (i = hmin; i <= hmax && i < h->size; i++)
	if (mpq_EGlpNumIsLess (val, h->key[h->entry[i]])) {
	    mpq_EGlpNumCopy (val, h->key[h->entry[i]]);
	    mc = i;
	    ILL_IFTRACE2 (":%d:%la", mc, mpq_EGlpNumToLf (val));
	}
    mpq_EGlpNumClearVar (val);
    ILL_IFTRACE2 ("\n");
    return mc;
}

#if mpq_DEBUG_HEAP > 0

static void mpq_printheap (mpq_heap * h)
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
	printf ("%la ", mpq_EGlpNumToLf (h->key[i]));
    printf ("\n key(sorted): ");
    for (i = 0; i < h->size; i++)
	printf ("%la ", mpq_EGlpNumToLf (h->key[h->entry[i]]));
    printf ("\n");
}

static void mpq_heapcheck (mpq_heap * h)
{
    int i, tcnt = 0;

    for (i = 0; i < h->maxsize; i++) {
	if (h->loc[i] < -1)
	    printf ("error in mpq_heap\n");
	else if (h->loc[i] > -1)
	    tcnt++;
    }
    if (tcnt != h->size)
	printf ("error 3 in mpq_heap\n");

    for (i = 0; i < h->size; i++) {
	if (h->loc[h->entry[i]] != i)
	    printf ("error 1 in mpq_heap\n");
	if (mpq_EGlpNumIsEqqual (h->key[h->entry[i]], mpq_zeroLpNum))
	    printf ("error 2 in mpq_heap\n");
	if (mpq_EGlpNumIsLess (h->key[h->entry[mpq_HEAP_UP (i)]], h->key[h->entry[i]]))
	    printf ("error 4 in mpq_heap\n");
    }
}

#endif

void mpq_ILLheap_insert (mpq_heap * const h,
      int const ix)
{
    int i = h->size;
    ILL_IFTRACE ("%s:%d:%la\n", __func__, ix, mpq_EGlpNumToLf (h->key[ix]));

    i = mpq_siftup (h, i, ix);
    h->size++;

#if mpq_DEBUG_HEAP > 0
    mpq_heapcheck (h);
#endif
#if mpq_DEBUG_HEAP > 1
    mpq_printheap (h);
#endif
}

void mpq_ILLheap_modify (mpq_heap * const h,
      int const ix)
{
    int i = h->loc[ix];
    int pi = i;
    ILL_IFTRACE ("%s:%d\n", __func__, ix);

    if (h->loc[ix] == -1)
	return;
    i = mpq_siftup (h, i, ix);
    if (pi == i)
	i = mpq_siftdown (h, i, ix);

#if mpq_DEBUG_HEAP > 0
    mpq_heapcheck (h);
#endif
#if mpq_DEBUG_HEAP > 1
    mpq_printheap (h);
#endif
}

void mpq_ILLheap_delete (mpq_heap * const h,
      int const ix)
{
    int i = h->loc[ix];
    int pi = i;
    int nix = h->entry[h->size - 1];
    ILL_IFTRACE ("%s:%d:%d:%d\n", __func__, ix, nix, pi);

    h->loc[ix] = -1;
    h->size--;
    if (nix == ix) {
#if mpq_DEBUG_HEAP > 0
	mpq_heapcheck (h);
#endif
#if mpq_DEBUG_HEAP > 1
	mpq_printheap (h);
#endif
	return;
    }
    h->entry[i] = nix;
    h->loc[nix] = i;

    i = mpq_siftup (h, i, nix);
    ILL_IFTRACE ("%s:%d:%d:%d:%d\n", __func__, ix, nix, pi, i);
    if (pi == i)
	mpq_siftdown (h, i, nix);

#if mpq_DEBUG_HEAP > 0
    mpq_heapcheck (h);
#endif
#if mpq_DEBUG_HEAP > 1
    mpq_printheap (h);
#endif
}

int mpq_ILLheap_findmin (mpq_heap * const h)
{
    if (h->hexist == 0 || h->size <= 0)
	return -1;
    return h->entry[0];
}

void mpq_ILLheap_init (mpq_heap * const h)
{
    h->entry = NULL;
    h->loc = NULL;
    h->key = NULL;
    h->hexist = 0;
}

int mpq_ILLheap_build (mpq_heap * const h,
      int const nelems,
      mpq_t * key)
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
	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, key[i])) {
	    h->entry[n] = i;
	    h->loc[i] = n;
	    n++;
	} else
	    h->loc[i] = -1;
    }
    h->size = n;
    for (i = n - 1; i >= 0; i--) {
	ILL_IFTRACE2 ("insert %la\n", mpq_EGlpNumToLf (h->key[h->entry[i]]));
	mpq_siftdown (h, i, h->entry[i]);
    }

#if mpq_DEBUG_HEAP > 0
    mpq_heapcheck (h);
#endif
#if mpq_DEBUG_HEAP > 1
    mpq_printheap (h);
#endif

CLEANUP:
    if (rval)
	mpq_ILLheap_free (h);
    ILL_RETURN (rval, "mpq_ILLheap_init");
}

void mpq_ILLheap_free (mpq_heap * const h)
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

void mpq_ILLmatrix_init (mpq_ILLmatrix * A)
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

void mpq_ILLmatrix_free (mpq_ILLmatrix * A)
{
    if (A) {
	mpq_EGlpNumFreeArray (A->matval);
	ILL_IFFREE (A->matcnt, int);
	ILL_IFFREE (A->matbeg, int);
	ILL_IFFREE (A->matind, int);
	mpq_ILLmatrix_init (A);
    }
}

void mpq_ILLmatrix_prt (FILE * fd,
      mpq_ILLmatrix * A)
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
		fprintf (fd, "row %d=%.3f ", A->matind[k], mpq_EGlpNumToLf (A->matval[k]));
	    }
	    fprintf (fd, "\n");
	}
    }
}
