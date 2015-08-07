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

/* RCSINFO $Id: mpq_basis.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpq___BASIS_H
#define mpq___BASIS_H

#include "config.h"
#include "mpq_dstruct.h"
#include "mpq_lpdefs.h"
#include "mpq_lpdata.h"

#if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE
extern mpq_t mpq_CB_PRI_RLIMIT;	/* = 0.25 */
extern mpq_t mpq_CB_INF_RATIO;	/* = 10.0 */
extern mpq_t mpq_CB_EPS;				/* = 0.001 */
#endif

typedef struct mpq_var_data
{
	int nartif;
	int nslacks;
	int nfree;
	int nbndone;
	int nbounded;
	int nfixed;
	mpq_t cmax;
}
mpq_var_data;

void mpq_ILLbasis_init_vardata (mpq_var_data * vd);
void mpq_ILLbasis_clear_vardata (mpq_var_data * vd);

int mpq_ILLbasis_build_basisinfo (mpq_lpinfo * lp),
  mpq_ILLbasis_get_initial (mpq_lpinfo * lp,
												int algorithm),
  mpq_ILLbasis_get_cinitial (mpq_lpinfo * lp,
												 int algorithm),
  mpq_ILLbasis_load (mpq_lpinfo * lp,
								 mpq_ILLlp_basis * B),
  mpq_ILLbasis_tableau_row (mpq_lpinfo * lp,
												int row,
												mpq_t * brow,
												mpq_t * trow,
												mpq_t * rhs,
												int strict),
  mpq_ILLbasis_factor (mpq_lpinfo * lp,
									 int *singular),
  mpq_ILLbasis_refactor (mpq_lpinfo * lp),
  mpq_ILLbasis_update (mpq_lpinfo * lp,
									 mpq_svector * y,
									 int lindex,
									 int *refactor,
									 int *singular);

void mpq_ILLbasis_column_solve (mpq_lpinfo * lp,
														mpq_svector * rhs,
														mpq_svector * soln),
  mpq_ILLbasis_column_solve_update (mpq_lpinfo * lp,
																mpq_svector * rhs,
																mpq_svector * upd,
																mpq_svector * soln),
  mpq_ILLbasis_row_solve (mpq_lpinfo * lp,
											mpq_svector * rhs,
											mpq_svector * soln),
  mpq_ILLbasis_free_basisinfo (mpq_lpinfo * lp),
  mpq_ILLbasis_free_fbasisinfo (mpq_lpinfo * lp),
  mpq_ILLbasis_init_basisinfo (mpq_lpinfo * lp);

#endif /* mpq___BASIS_H */
