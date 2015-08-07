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

/* RCSINFO $Id: dbl_basis.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef dbl___BASIS_H
#define dbl___BASIS_H

#include "config.h"
#include "dbl_dstruct.h"
#include "dbl_lpdefs.h"
#include "dbl_lpdata.h"

#if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE
extern double dbl_CB_PRI_RLIMIT;	/* = 0.25 */
extern double dbl_CB_INF_RATIO;	/* = 10.0 */
extern double dbl_CB_EPS;				/* = 0.001 */
#endif

typedef struct dbl_var_data
{
	int nartif;
	int nslacks;
	int nfree;
	int nbndone;
	int nbounded;
	int nfixed;
	double cmax;
}
dbl_var_data;

void dbl_ILLbasis_init_vardata (dbl_var_data * vd);
void dbl_ILLbasis_clear_vardata (dbl_var_data * vd);

int dbl_ILLbasis_build_basisinfo (dbl_lpinfo * lp),
  dbl_ILLbasis_get_initial (dbl_lpinfo * lp,
												int algorithm),
  dbl_ILLbasis_get_cinitial (dbl_lpinfo * lp,
												 int algorithm),
  dbl_ILLbasis_load (dbl_lpinfo * lp,
								 dbl_ILLlp_basis * B),
  dbl_ILLbasis_tableau_row (dbl_lpinfo * lp,
												int row,
												double * brow,
												double * trow,
												double * rhs,
												int strict),
  dbl_ILLbasis_factor (dbl_lpinfo * lp,
									 int *singular),
  dbl_ILLbasis_refactor (dbl_lpinfo * lp),
  dbl_ILLbasis_update (dbl_lpinfo * lp,
									 dbl_svector * y,
									 int lindex,
									 int *refactor,
									 int *singular);

void dbl_ILLbasis_column_solve (dbl_lpinfo * lp,
														dbl_svector * rhs,
														dbl_svector * soln),
  dbl_ILLbasis_column_solve_update (dbl_lpinfo * lp,
																dbl_svector * rhs,
																dbl_svector * upd,
																dbl_svector * soln),
  dbl_ILLbasis_row_solve (dbl_lpinfo * lp,
											dbl_svector * rhs,
											dbl_svector * soln),
  dbl_ILLbasis_free_basisinfo (dbl_lpinfo * lp),
  dbl_ILLbasis_free_fbasisinfo (dbl_lpinfo * lp),
  dbl_ILLbasis_init_basisinfo (dbl_lpinfo * lp);

#endif /* dbl___BASIS_H */
