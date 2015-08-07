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

/* RCSINFO $Id: dbl_util.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef dbl_ILL_UTIL_H
#define dbl_ILL_UTIL_H

#include "machdefs.h"
#include "config.h"

#ifdef _USRDLL

#ifdef QSLIB_EXPORTS
#define dbl_QSLIB_INTERFACE __declspec(dllexport)
#else
#define dbl_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else

#define dbl_QSLIB_INTERFACE extern

#endif

#ifdef WIN32
#define strcasecmp(s1, s2) 	stricmp(s1, s2)
#define strncasecmp(s1, s2, n) 	strnicmp(s1, s2, n)
#endif

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
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  EGLPNUMPT_SWAP                                                          */
/*  dbl_ILL_SWAP(a,b,t)                                                         */
/*    swaps a and b, using t as temporary space.  a, b, and t should all    */
/*    be the same type.                                                     */
/*                                                                          */
/*  dbl_ILL_OURABS(a)                                                           */
/*    returns the absolute value of a.                                      */
/*                                                                          */
/****************************************************************************/
typedef char dbl_ILLbool;
#define dbl_FALSE 0
#define dbl_TRUE  1

#define dbl_ILL_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))
#define dbl_EGLPNUM_SWAP(a,b,t) ((dbl_EGlpNumCopy(t,a)),(dbl_EGlpNumCopy(a,b)),(dbl_EGlpNumCopy(b,t)))

#define dbl_ILL_OURABS(a) (((a) >= 0) ? (a) : -(a))

#include "dbl_sortrus.h"
#include "allocrus.h"
#include "urandom.h"
#include "zeit.h"
/****************************************************************************/
/*                                                                          */
/*                             dbl_util.c                                       */
/*                                                                          */
/****************************************************************************/
#define dbl_ILL_UTIL_STR(new, str) \
    { new = dbl_ILLutil_str(str); \
      if (str != NULL) { ILL_CHECKnull(new, "out of memeory"); } }

extern char *dbl_ILLutil_str (const char *str);
	 /* allocates and returns a copy of s */

extern int dbl_ILLutil_array_index (char *list[],
																int n,
																const char *name);
	 /* returns index of name in list or -1  */

extern int dbl_ILLutil_index (const char *list[],
													const char *name);
	 /* returns index of name in list or -1  */

extern unsigned int dbl_ILLutil_nextprime (unsigned int x);

extern const char *dbl_ILLutil_strchr (const char *s,
																	 int c);

extern int dbl_ILLutil_strcasecmp (const char *s1,
															 const char *s2);
extern int dbl_ILLutil_strncasecmp (const char *s1,
																const char *s2,
																size_t n);


extern int dbl_ILLutil_our_gcd (int a,
														int b),
  dbl_ILLutil_our_lcm (int a,
									 int b),
  dbl_ILLutil_our_log2 (int a);

double dbl_ILLutil_our_floor (double x),
  dbl_ILLutil_our_ceil (double x),
  dbl_ILLutil_our_frac (double x),
  dbl_ILLutil_norm_sqr (double *v,
										int len);
#include "bgetopt.h"
#include "dbl_dheaps_i.h"
#include "dbl_priority.h"
#endif /* dbl_ILL_UTIL_H */
