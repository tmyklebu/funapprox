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

#ifndef mpf___SORTRUS_H__
#define mpf___SORTRUS_H__
/****************************************************************************/
/*                                                                          */
/*                             mpf_sortrus.c                                    */
/*                                                                          */
/****************************************************************************/
#include "econfig.h"
#include "urandom.h"
void mpf_ILLutil_int_array_quicksort (int *len, int n),
  mpf_ILLutil_int_perm_quicksort (int *perm, int *len, int n),
  mpf_ILLutil_double_perm_quicksort (int *perm, double *len, int n),
  mpf_ILLutil_EGlpNum_perm_quicksort (int *perm, mpf_t * len, int n),
  mpf_ILLutil_str_perm_quicksort (int *perm, char **len, int n),
  mpf_ILLutil_EGlpNum_rselect (int *arr, int l, int r, int m, mpf_t * coord,
      ILLrandstate * rstate),
  mpf_ILLutil_rselect (int *arr, int l, int r, int m, double *coord,
      ILLrandstate * rstate);

char *mpf_ILLutil_linked_radixsort (char *data, char *datanext, char *dataval,
      int valsize);

#endif
