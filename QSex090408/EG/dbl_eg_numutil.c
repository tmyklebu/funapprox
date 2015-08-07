/* EGlib "Efficient General Library" provides some basic structures and
   algorithms commons in many optimization algorithms.

Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.

This library is free software; you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

You should have received a copy of the GNU Lesser General Public License
   along with this library; if not, write to the Free Software Foundation,
   Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */
#include "dbl_eg_numutil.h"

/* ========================================================================= */
void dbl_EGutilPermSort (const size_t sz, int *const perm,
      const double *const elem)
{
    size_t i, j, temp;
    double t;
    if (sz <= 1)
	return;

    dbl_EGlpNumInitVar (t);
    EGswap (perm[0], perm[(sz - 1) / 2], temp);
    i = 0;
    j = sz;
    dbl_EGlpNumCopy (t, elem[perm[0]]);
    for (;;) {
	do
	    i++;
	while (i < sz && dbl_EGlpNumIsLess (elem[perm[i]], t));
	do
	    j--;
	while (dbl_EGlpNumIsLess (t, elem[perm[j]]));
	if (j < i)
	    break;
	EGswap (perm[i], perm[j], temp);
    }
    EGswap (perm[0], perm[j], temp);
    dbl_EGlpNumClearVar (t);
    dbl_EGutilPermSort (j, perm, elem);
    dbl_EGutilPermSort (sz - i, perm + i, elem);
}

/* ========================================================================= */
void dbl___EGlpNumInnProd (double *const rop, double *const __EGa1,
      double *const __EGa2, const size_t length)
{
    size_t __EGdim = length;
    dbl_EGlpNumZero ((*rop));
    while (__EGdim--)
	dbl_EGlpNumAddInnProdTo ((*rop), __EGa1[__EGdim], __EGa2[__EGdim]);
}

