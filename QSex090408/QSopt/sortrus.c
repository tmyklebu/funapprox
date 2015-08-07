/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza.  Sanjeeb Dash's ownership of copyright in          */
/*  QSopt_ex is derived from his copyright in QSopt.                        */
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

/* RCSINFO $Id: sortrus.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  char *ILLutil_linked_radixsort (char *data, char *datanext,              */
/*      char *dataval, int valsize)                                         */
/*    USAGE:                                                                */
/*      head = (bar *) ILLutil_linked_radixsort ((char *) head,              */
/*         (char *) &(head->next), (char *) &(head->val), sizeof (int));    */
/*    Then head is the start of the linked list in increasing order of      */
/*    val, with next as the field that links the bars.                      */
/*    WARNING: DOES NOT HANDLE NEGATIVE NUMBERS PROPERLY.                   */
/*                                                                          */
/*  void ILLutil_int_array_quicksort (int *len, int n)                       */
/*    len - the array to be sorted                                          */
/*    n - the number of elements in len                                     */
/*    Uses quicksort to put len in increasing order.                        */
/*                                                                          */
/*  void ILLutil_int_perm_quicksort (int *perm, int *len, int n)             */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void ILLutil_double_perm_quicksort (int *perm, double *len, int n)       */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void ILLutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord, ILLrandstate *rstate)                                 */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

#include "econfig.h"
#include "machdefs.h"
#include "util.h"
#include "except.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


#define BITS_PER_PASS (8)

#define NBINS (1<<BITS_PER_PASS)


static void select_split (int *arr,
													int n,
													double v,
													int *start,
													int *end,
													double *coord),
  select_sort (int *arr,
							 int n,
							 double *coord),
  select_sort_dsample (double *samp,
											 int n),
  select_EGlpNum_split (int *arr,
												int n,
												EGlpNum_t * v,
												int *start,
												int *end,
												EGlpNum_t * coord),
  select_EGlpNum_sort (int *arr,
											 int n,
											 EGlpNum_t * coord),
  select_EGlpNum_sort_dsample (EGlpNum_t * samp,
															 int n);



char *ILLutil_linked_radixsort (char *data,
																char *datanext,
																char *dataval,
																int valsize)
{
	size_t nextoff = datanext - data;
	size_t valoff = dataval - data;
	int i;
	char *head[NBINS];
	char **tail[NBINS];
	char *p;
	char **last;
	int j;
	int v;


	for (j = valsize - 1; j >= 0; j--)
	{
		for (i = 0; i < NBINS; i++)
		{
			head[i] = (char *) NULL;
			tail[i] = &head[i];
		}
		for (p = data; p; p = *(char **) (p + nextoff))
		{
			v = (unsigned char) p[valoff + j];
			*tail[v] = p;
			tail[v] = (char **) (p + nextoff);
		}
		last = &data;
		for (i = 0; i < NBINS; i++)
		{
			if (head[i])
			{
				*last = head[i];
				last = tail[i];
			}
		}
		*last = (char *) NULL;
	}
	return data;
}

void ILLutil_int_array_quicksort (int *len,
																	int n)
{
	int i,
	  j,
	  temp,
	  t;

	if (n <= 1)
		return;

	ILL_SWAP (len[0], len[(n - 1) / 2], temp);

	i = 0;
	j = n;
	t = len[0];

	for (;;)
	{
		do
			i++;
		while (i < n && len[i] < t);
		do
			j--;
		while (len[j] > t);
		if (j < i)
			break;
		ILL_SWAP (len[i], len[j], temp);
	}
	ILL_SWAP (len[0], len[j], temp);

	ILLutil_int_array_quicksort (len, j);
	ILLutil_int_array_quicksort (len + i, n - i);
}

void ILLutil_int_perm_quicksort (int *perm,
																 int *len,
																 int n)
{
	int i,
	  j,
	  temp,
	  t;

	if (n <= 1)
		return;

	ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

	i = 0;
	j = n;
	t = len[perm[0]];

	for (;;)
	{
		do
			i++;
		while (i < n && len[perm[i]] < t);
		do
			j--;
		while (len[perm[j]] > t);
		if (j < i)
			break;
		ILL_SWAP (perm[i], perm[j], temp);
	}
	ILL_SWAP (perm[0], perm[j], temp);

	ILLutil_int_perm_quicksort (perm, len, j);
	ILLutil_int_perm_quicksort (perm + i, len, n - i);
}


void ILLutil_EGlpNum_perm_quicksort (int *perm,
																		 EGlpNum_t * len,
																		 int n)
{
	int i,
	  j,
	  temp;
	EGlpNum_t t;
	if (n <= 1)
		return;

	EGlpNumInitVar (t);
	ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

	i = 0;
	j = n;
	EGlpNumCopy (t, len[perm[0]]);

	for (;;)
	{
		do
			i++;
		while (i < n && EGlpNumIsLess (len[perm[i]], t));
		do
			j--;
		while (EGlpNumIsLess (t, len[perm[j]]));
		if (j < i)
			break;
		ILL_SWAP (perm[i], perm[j], temp);
	}
	ILL_SWAP (perm[0], perm[j], temp);

	EGlpNumClearVar (t);
	ILLutil_EGlpNum_perm_quicksort (perm, len, j);
	ILLutil_EGlpNum_perm_quicksort (perm + i, len, n - i);
}

void ILLutil_double_perm_quicksort (int *perm,
																		double *len,
																		int n)
{
	int i,
	  j,
	  temp;
	double t;

	if (n <= 1)
		return;

	ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

	i = 0;
	j = n;
	t = len[perm[0]];

	for (;;)
	{
		do
			i++;
		while (i < n && len[perm[i]] < t);
		do
			j--;
		while (t < len[perm[j]]);
		if (j < i)
			break;
		ILL_SWAP (perm[i], perm[j], temp);
	}
	ILL_SWAP (perm[0], perm[j], temp);

	ILLutil_double_perm_quicksort (perm, len, j);
	ILLutil_double_perm_quicksort (perm + i, len, n - i);
}

void ILLutil_str_perm_quicksort (int *perm,
																 char **len,
																 int n)
{
	int i,
	  j,
	  temp;
	char *t;

	if (n <= 1)
		return;

	ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

	i = 0;
	j = n;
	t = len[perm[0]];

	for (;;)
	{
		do
			i++;
		while (i < n && (strcmp (len[perm[i]], t) < 0));
		do
			j--;
		while (strcmp (len[perm[j]], t) > 0);
		if (j < i)
			break;
		ILL_SWAP (perm[i], perm[j], temp);
	}
	ILL_SWAP (perm[0], perm[j], temp);

	ILLutil_str_perm_quicksort (perm, len, j);
	ILLutil_str_perm_quicksort (perm + i, len, n - i);
}


/**********  Median - Select Routines **********/

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20


void ILLutil_EGlpNum_rselect (int *arr,
															int l,
															int r,
															int m,
															EGlpNum_t * coord,
															ILLrandstate * rstate)
{
	EGlpNum_t *samplevals = EGlpNumAllocArray (NSAMPLES);
	int i;
	int st,
	  en;
	int n;

	arr += l;
	n = r - l + 1;
	m -= l;

	while (n > SORTSIZE)
	{
		for (i = 0; i < NSAMPLES; i++)
		{
			EGlpNumCopy (samplevals[i], coord[arr[ILLutil_lprand (rstate) % n]]);
		}
		select_EGlpNum_sort_dsample (samplevals, NSAMPLES);
		select_EGlpNum_split (arr, n, &(samplevals[(NSAMPLES - 1) / 2]),
													&st, &en, coord);
		if (st > m)
		{
			n = st;
		}
		else if (en <= m)
		{
			arr += en;
			n -= en;
			m -= en;
		}
		else
		{
			return;
		}
	}

	select_EGlpNum_sort (arr, n, coord);
	EGlpNumFreeArray (samplevals);
	return;
}


void ILLutil_rselect (int *arr,
											int l,
											int r,
											int m,
											double *coord,
											ILLrandstate * rstate)
{
	double samplevals[NSAMPLES];
	int i;
	int st,
	  en;
	int n;

	arr += l;
	n = r - l + 1;
	m -= l;

	while (n > SORTSIZE)
	{
		for (i = 0; i < NSAMPLES; i++)
		{
			samplevals[i] = coord[arr[ILLutil_lprand (rstate) % n]];
		}
		select_sort_dsample (samplevals, NSAMPLES);
		select_split (arr, n, samplevals[(NSAMPLES - 1) / 2], &st, &en, coord);
		if (st > m)
		{
			n = st;
		}
		else if (en <= m)
		{
			arr += en;
			n -= en;
			m -= en;
		}
		else
		{
			return;
		}
	}

	select_sort (arr, n, coord);
	return;
}

static void select_split (int *arr,
													int n,
													double v,
													int *start,
													int *end,
													double *coord)
{
	int i,
	  j,
	  k;
	int t;

	i = 0;
	j = k = n;

	while (i < j)
	{
		if (coord[arr[i]] < v)
		{
			i++;
		}
		else if (coord[arr[i]] == v)
		{
			j--;
			ILL_SWAP (arr[i], arr[j], t);
		}
		else
		{
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

static void select_sort (int *arr,
												 int n,
												 double *coord)
{
	int i,
	  j;
	int t;

	for (i = 1; i < n; i++)
	{
		t = arr[i];
		for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--)
		{
			arr[j] = arr[j - 1];
		}
		arr[j] = t;
	}
}

static void select_sort_dsample (double *samp,
																 int n)
{
	int i,
	  j;
	double t;

	for (i = 1; i < n; i++)
	{
		t = samp[i];
		for (j = i; j > 0 && samp[j - 1] > t; j--)
		{
			samp[j] = samp[j - 1];
		}
		samp[j] = t;
	}
}

static void select_EGlpNum_split (int *arr,
																	int n,
																	EGlpNum_t * v,
																	int *start,
																	int *end,
																	EGlpNum_t * coord)
{
	int i,
	  j,
	  k;
	int t;

	i = 0;
	j = k = n;

	while (i < j)
	{
		if (EGlpNumIsLess (coord[arr[i]], *v))
		{
			i++;
		}
		else if (EGlpNumIsEqqual (coord[arr[i]], *v))
		{
			j--;
			ILL_SWAP (arr[i], arr[j], t);
		}
		else
		{
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

static void select_EGlpNum_sort (int *arr,
																 int n,
																 EGlpNum_t * coord)
{
	int i,
	  j;
	int t;

	for (i = 1; i < n; i++)
	{
		t = arr[i];
		for (j = i; j > 0 && EGlpNumIsLess (coord[t], coord[arr[j - 1]]); j--)
		{
			arr[j] = arr[j - 1];
		}
		arr[j] = t;
	}
}

static void select_EGlpNum_sort_dsample (EGlpNum_t * samp,
																				 int n)
{
	int i,
	  j;
	EGlpNum_t t;
	EGlpNumInitVar (t);

	for (i = 1; i < n; i++)
	{
		EGlpNumCopy (t, samp[i]);
		for (j = i; j > 0 && EGlpNumIsLess (t, samp[j - 1]); j--)
		{
			EGlpNumCopy (samp[j], samp[j - 1]);
		}
		EGlpNumCopy (samp[j], t);
	}
	EGlpNumClearVar (t);
}



#ifdef  TRY_CODE
static void perm_init (int *perm,
											 int offset,
											 int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		perm[i + offset] = i;
	}
}
int main (int argc,
					char **argv)
{
	ILLrandstate *rstate;
	double *dlen;
	int *perm;
	int i,
	  n,
	  m,
	  fr;
	int rval = 0;

	argc--;
	ILL_NEW (rstate, ILLrandstate);
	ILLutil_sprand (10, rstate);
	ILL_SAFE_MALLOC (dlen, argc, double);
	ILL_SAFE_MALLOC (perm, argc, int);

	for (i = 0; i < argc; i++)
	{
		sscanf (argv[i + 1], "%lf", dlen + i);
	}
	for (fr = 0; fr <= argc / 2; fr++)
	{
		n = argc - 2 * fr;
		m = n / 2;

		perm_init (perm, 0, argc);
		fprintf (stdout, "rselect: fr=%d n=%d m=%d: ", fr, n, m);
		for (i = 0; i < n; i++)
		{
			fprintf (stdout, "%.1f ", dlen[perm[fr + i]]);
		}

		ILLutil_rselect (perm, fr, fr + n - 1, m, dlen, rstate);

		fprintf (stdout, " -> ");
		for (i = 0; i < m; i++)
		{
			fprintf (stdout, "%.1f ", dlen[perm[fr + i]]);
		}
		fprintf (stdout, "\n");
	}
CLEANUP:
	exit (rval);
}


#endif
