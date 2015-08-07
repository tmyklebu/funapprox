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

/* RCSINFO $Id: mpf_basis.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpf___BASIS_H
#define mpf___BASIS_H

#include "config.h"
#include "mpf_dstruct.h"
#include "mpf_lpdefs.h"
#include "mpf_lpdata.h"

#if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE
extern mpf_t mpf_CB_PRI_RLIMIT;	/* = 0.25 */
extern mpf_t mpf_CB_INF_RATIO;	/* = 10.0 */
extern mpf_t mpf_CB_EPS;				/* = 0.001 */
#endif

typedef struct mpf_var_data
{
	int nartif;
	int nslacks;
	int nfree;
	int nbndone;
	int nbounded;
	int nfixed;
	mpf_t cmax;
}
mpf_var_data;

void mpf_ILLbasis_init_vardata (mpf_var_data * vd);
void mpf_ILLbasis_clear_vardata (mpf_var_data * vd);

int mpf_ILLbasis_build_basisinfo (mpf_lpinfo * lp),
  mpf_ILLbasis_get_initial (mpf_lpinfo * lp,
												int algorithm),
  mpf_ILLbasis_get_cinitial (mpf_lpinfo * lp,
												 int algorithm),
  mpf_ILLbasis_load (mpf_lpinfo * lp,
								 mpf_ILLlp_basis * B),
  mpf_ILLbasis_tableau_row (mpf_lpinfo * lp,
												int row,
												mpf_t * brow,
												mpf_t * trow,
												mpf_t * rhs,
												int strict),
  mpf_ILLbasis_factor (mpf_lpinfo * lp,
									 int *singular),
  mpf_ILLbasis_refactor (mpf_lpinfo * lp),
  mpf_ILLbasis_update (mpf_lpinfo * lp,
									 mpf_svector * y,
									 int lindex,
									 int *refactor,
									 int *singular);

void mpf_ILLbasis_column_solve (mpf_lpinfo * lp,
														mpf_svector * rhs,
														mpf_svector * soln),
  mpf_ILLbasis_column_solve_update (mpf_lpinfo * lp,
																mpf_svector * rhs,
																mpf_svector * upd,
																mpf_svector * soln),
  mpf_ILLbasis_row_solve (mpf_lpinfo * lp,
											mpf_svector * rhs,
											mpf_svector * soln),
  mpf_ILLbasis_free_basisinfo (mpf_lpinfo * lp),
  mpf_ILLbasis_free_fbasisinfo (mpf_lpinfo * lp),
  mpf_ILLbasis_init_basisinfo (mpf_lpinfo * lp);

#endif /* mpf___BASIS_H */
