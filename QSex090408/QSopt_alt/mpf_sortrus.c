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

/* RCSINFO $Id: mpf_sortrus.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/* */
/* This file is part of CONCORDE                                           */
/* */
/* (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/* Vasek Chvatal, and William Cook                                         */
/* */
/* Permission is granted for academic research use.  For other uses,       */
/* contact the authors for licensing options.                              */
/* */
/* Use at your own risk.  We make no guarantees about the                  */
/* correctness or usefulness of this code.                                 */
/* */
/****************************************************************************/

/****************************************************************************/
/* */
/* SORTING ROUTINES                                 */
/* */
/* TSP CODE                                     */
/* */
/* Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/* DATE:  February 24, 1994                                               */
/* */
/* EXPORTED FUNCTIONS:                                                   */
/* */
/* char *mpf_ILLutil_linked_radixsort (char *data, char *datanext,              */
/* char *dataval, int valsize)                                         */
/* USAGE:                                                                */
/* head = (bar *) mpf_ILLutil_linked_radixsort ((char *) head,              */
/* (char *) &(head->next), (char *) &(head->val), sizeof (int));    */
/* Then head is the start of the linked list in increasing order of      */
/* val, with next as the field that links the bars.                      */
/* WARNING: DOES NOT HANDLE NEGATIVE NUMBERS PROPERLY.                   */
/* */
/* void mpf_ILLutil_int_array_quicksort (int *len, int n)                       */
/* len - the array to be sorted                                          */
/* n - the number of elements in len                                     */
/* Uses quicksort to put len in increasing order.                        */
/* */
/* void mpf_ILLutil_int_perm_quicksort (int *perm, int *len, int n)             */
/* perm - must be allocated and initialized by the calling routine,      */
/* it will be arranged in increasing order of len.                */
/* n - the number of elements in perm and len.                           */
/* */
/* void mpf_ILLutil_double_perm_quicksort (int *perm, double *len, int n)       */
/* perm - must be allocated and initialized by the calling routine,      */
/* it will be arranged in increasing order of len.                */
/* n - the number of elements in perm and len.                           */
/* */
/* void mpf_ILLutil_rselect (int *arr, int l, int r, int m,                     */
/* double *coord, ILLrandstate *rstate)                                 */
/* arr - permutation that will be rearranged                             */
/* l,r - specify the range of arr that we are interested in              */
/* m - is the index into l,r that is the break point for the perm        */
/* coord - gives the keys that determine the ordering                    */
/* */
/****************************************************************************/

#include "econfig.h"
#include "machdefs.h"
#include "mpf_util.h"
#include "except.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


#define mpf_BITS_PER_PASS (8)

#define mpf_NBINS (1<<mpf_BITS_PER_PASS)


static void mpf_select_split (int *arr,
      int n,
      double v,
      int *start,
      int *end,
      double *coord),
  mpf_select_sort (int *arr,
      int n,
      double *coord),
  mpf_select_sort_dsample (double *samp,
      int n),
  mpf_select_EGlpNum_split (int *arr,
      int n,
      mpf_t * v,
      int *start,
      int *end,
      mpf_t * coord),
  mpf_select_EGlpNum_sort (int *arr,
      int n,
      mpf_t * coord),
  mpf_select_EGlpNum_sort_dsample (mpf_t * samp,
      int n);



char *mpf_ILLutil_linked_radixsort (char *data,
      char *datanext,
      char *dataval,
      int valsize)
{
    size_t nextoff = datanext - data;
    size_t valoff = dataval - data;
    int i;
    char *head[mpf_NBINS];
    char **tail[mpf_NBINS];
    char *p;
    char **last;
    int j;
    int v;


    for (j = valsize - 1; j >= 0; j--) {
	for (i = 0; i < mpf_NBINS; i++) {
	    head[i] = (char *) NULL;
	    tail[i] = &head[i];
	}
	for (p = data; p; p = *(char **) (p + nextoff)) {
	    v = (unsigned char) p[valoff + j];
	    *tail[v] = p;
	    tail[v] = (char **) (p + nextoff);
	}
	last = &data;
	for (i = 0; i < mpf_NBINS; i++) {
	    if (head[i]) {
		*last = head[i];
		last = tail[i];
	    }
	}
	*last = (char *) NULL;
    }
    return data;
}

void mpf_ILLutil_int_array_quicksort (int *len,
      int n)
{
    int i, j, temp, t;

    if (n <= 1)
	return;

    mpf_ILL_SWAP (len[0], len[(n - 1) / 2], temp);

    i = 0;
    j = n;
    t = len[0];

    for (;;) {
	do
	    i++;
	while (i < n && len[i] < t);
	do
	    j--;
	while (len[j] > t);
	if (j < i)
	    break;
	mpf_ILL_SWAP (len[i], len[j], temp);
    }
    mpf_ILL_SWAP (len[0], len[j], temp);

    mpf_ILLutil_int_array_quicksort (len, j);
    mpf_ILLutil_int_array_quicksort (len + i, n - i);
}

void mpf_ILLutil_int_perm_quicksort (int *perm,
      int *len,
      int n)
{
    int i, j, temp, t;

    if (n <= 1)
	return;

    mpf_ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    for (;;) {
	do
	    i++;
	while (i < n && len[perm[i]] < t);
	do
	    j--;
	while (len[perm[j]] > t);
	if (j < i)
	    break;
	mpf_ILL_SWAP (perm[i], perm[j], temp);
    }
    mpf_ILL_SWAP (perm[0], perm[j], temp);

    mpf_ILLutil_int_perm_quicksort (perm, len, j);
    mpf_ILLutil_int_perm_quicksort (perm + i, len, n - i);
}


void mpf_ILLutil_EGlpNum_perm_quicksort (int *perm,
      mpf_t * len,
      int n)
{
    int i, j, temp;
    mpf_t t;
    if (n <= 1)
	return;

    mpf_EGlpNumInitVar (t);
    mpf_ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

    i = 0;
    j = n;
    mpf_EGlpNumCopy (t, len[perm[0]]);

    for (;;) {
	do
	    i++;
	while (i < n && mpf_EGlpNumIsLess (len[perm[i]], t));
	do
	    j--;
	while (mpf_EGlpNumIsLess (t, len[perm[j]]));
	if (j < i)
	    break;
	mpf_ILL_SWAP (perm[i], perm[j], temp);
    }
    mpf_ILL_SWAP (perm[0], perm[j], temp);

    mpf_EGlpNumClearVar (t);
    mpf_ILLutil_EGlpNum_perm_quicksort (perm, len, j);
    mpf_ILLutil_EGlpNum_perm_quicksort (perm + i, len, n - i);
}

void mpf_ILLutil_double_perm_quicksort (int *perm,
      double *len,
      int n)
{
    int i, j, temp;
    double t;

    if (n <= 1)
	return;

    mpf_ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    for (;;) {
	do
	    i++;
	while (i < n && len[perm[i]] < t);
	do
	    j--;
	while (t < len[perm[j]]);
	if (j < i)
	    break;
	mpf_ILL_SWAP (perm[i], perm[j], temp);
    }
    mpf_ILL_SWAP (perm[0], perm[j], temp);

    mpf_ILLutil_double_perm_quicksort (perm, len, j);
    mpf_ILLutil_double_perm_quicksort (perm + i, len, n - i);
}

void mpf_ILLutil_str_perm_quicksort (int *perm,
      char **len,
      int n)
{
    int i, j, temp;
    char *t;

    if (n <= 1)
	return;

    mpf_ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    for (;;) {
	do
	    i++;
	while (i < n && (strcmp (len[perm[i]], t) < 0));
	do
	    j--;
	while (strcmp (len[perm[j]], t) > 0);
	if (j < i)
	    break;
	mpf_ILL_SWAP (perm[i], perm[j], temp);
    }
    mpf_ILL_SWAP (perm[0], perm[j], temp);

    mpf_ILLutil_str_perm_quicksort (perm, len, j);
    mpf_ILLutil_str_perm_quicksort (perm + i, len, n - i);
}


/**********  Median - Select Routines **********/

/* mpf_NSAMPLES should be odd */
#define mpf_NSAMPLES 3
#define mpf_SORTSIZE 20


void mpf_ILLutil_EGlpNum_rselect (int *arr,
      int l,
      int r,
      int m,
      mpf_t * coord,
      ILLrandstate * rstate)
{
    mpf_t *samplevals = mpf_EGlpNumAllocArray (mpf_NSAMPLES);
    int i;
    int st, en;
    int n;

    arr += l;
    n = r - l + 1;
    m -= l;

    while (n > mpf_SORTSIZE) {
	for (i = 0; i < mpf_NSAMPLES; i++) {
	    mpf_EGlpNumCopy (samplevals[i], coord[arr[ILLutil_lprand (rstate) % n]]);
	}
	mpf_select_EGlpNum_sort_dsample (samplevals, mpf_NSAMPLES);
	mpf_select_EGlpNum_split (arr, n, &(samplevals[(mpf_NSAMPLES - 1) / 2]),
	    &st, &en, coord);
	if (st > m) {
	    n = st;
	} else if (en <= m) {
	    arr += en;
	    n -= en;
	    m -= en;
	} else {
	    return;
	}
    }

    mpf_select_EGlpNum_sort (arr, n, coord);
    mpf_EGlpNumFreeArray (samplevals);
    return;
}


void mpf_ILLutil_rselect (int *arr,
      int l,
      int r,
      int m,
      double *coord,
      ILLrandstate * rstate)
{
    double samplevals[mpf_NSAMPLES];
    int i;
    int st, en;
    int n;

    arr += l;
    n = r - l + 1;
    m -= l;

    while (n > mpf_SORTSIZE) {
	for (i = 0; i < mpf_NSAMPLES; i++) {
	    samplevals[i] = coord[arr[ILLutil_lprand (rstate) % n]];
	}
	mpf_select_sort_dsample (samplevals, mpf_NSAMPLES);
	mpf_select_split (arr, n, samplevals[(mpf_NSAMPLES - 1) / 2], &st, &en, coord);
	if (st > m) {
	    n = st;
	} else if (en <= m) {
	    arr += en;
	    n -= en;
	    m -= en;
	} else {
	    return;
	}
    }

    mpf_select_sort (arr, n, coord);
    return;
}

static void mpf_select_split (int *arr,
      int n,
      double v,
      int *start,
      int *end,
      double *coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j) {
	if (coord[arr[i]] < v) {
	    i++;
	} else if (coord[arr[i]] == v) {
	    j--;
	    mpf_ILL_SWAP (arr[i], arr[j], t);
	} else {
	    j--;
	    k--;
	    t = arr[i];
	    arr[i] = arr[j];
	    arr[j] = arr[k];
	    arr[k] = t;
	}
    }
    *start = j;
    *end = k;
    return;
}

static void mpf_select_sort (int *arr,
      int n,
      double *coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
	t = arr[i];
	for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--) {
	    arr[j] = arr[j - 1];
	}
	arr[j] = t;
    }
}

static void mpf_select_sort_dsample (double *samp,
      int n)
{
    int i, j;
    double t;

    for (i = 1; i < n; i++) {
	t = samp[i];
	for (j = i; j > 0 && samp[j - 1] > t; j--) {
	    samp[j] = samp[j - 1];
	}
	samp[j] = t;
    }
}

static void mpf_select_EGlpNum_split (int *arr,
      int n,
      mpf_t * v,
      int *start,
      int *end,
      mpf_t * coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j) {
	if (mpf_EGlpNumIsLess (coord[arr[i]], *v)) {
	    i++;
	} else if (mpf_EGlpNumIsEqqual (coord[arr[i]], *v)) {
	    j--;
	    mpf_ILL_SWAP (arr[i], arr[j], t);
	} else {
	    j--;
	    k--;
	    t = arr[i];
	    arr[i] = arr[j];
	    arr[j] = arr[k];
	    arr[k] = t;
	}
    }
    *start = j;
    *end = k;
    return;
}

static void mpf_select_EGlpNum_sort (int *arr,
      int n,
      mpf_t * coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
	t = arr[i];
	for (j = i; j > 0 && mpf_EGlpNumIsLess (coord[t], coord[arr[j - 1]]); j--) {
	    arr[j] = arr[j - 1];
	}
	arr[j] = t;
    }
}

static void mpf_select_EGlpNum_sort_dsample (mpf_t * samp,
      int n)
{
    int i, j;
    mpf_t t;
    mpf_EGlpNumInitVar (t);

    for (i = 1; i < n; i++) {
	mpf_EGlpNumCopy (t, samp[i]);
	for (j = i; j > 0 && mpf_EGlpNumIsLess (t, samp[j - 1]); j--) {
	    mpf_EGlpNumCopy (samp[j], samp[j - 1]);
	}
	mpf_EGlpNumCopy (samp[j], t);
    }
    mpf_EGlpNumClearVar (t);
}



#ifdef  TRY_CODE
static void mpf_perm_init (int *perm,
      int offset,
      int n)
{
    int i;
    for (i = 0; i < n; i++) {
	perm[i + offset] = i;
    }
}
int main (int argc,
      char **argv)
{
    ILLrandstate *rstate;
    double *dlen;
    int *perm;
    int i, n, m, fr;
    int rval = 0;

    argc--;
    ILL_NEW (rstate, ILLrandstate);
    ILLutil_sprand (10, rstate);
    ILL_SAFE_MALLOC (dlen, argc, double);
    ILL_SAFE_MALLOC (perm, argc, int);

    for (i = 0; i < argc; i++) {
	sscanf (argv[i + 1], "%lf", dlen + i);
    }
    for (fr = 0; fr <= argc / 2; fr++) {
	n = argc - 2 * fr;
	m = n / 2;

	mpf_perm_init (perm, 0, argc);
	fprintf (stdout, "rselect: fr=%d n=%d m=%d: ", fr, n, m);
	for (i = 0; i < n; i++) {
	    fprintf (stdout, "%.1f ", dlen[perm[fr + i]]);
	}

	mpf_ILLutil_rselect (perm, fr, fr + n - 1, m, dlen, rstate);

	fprintf (stdout, " -> ");
	for (i = 0; i < m; i++) {
	    fprintf (stdout, "%.1f ", dlen[perm[fr + i]]);
	}
	fprintf (stdout, "\n");
    }
CLEANUP:
    exit (rval);
}


#endif
