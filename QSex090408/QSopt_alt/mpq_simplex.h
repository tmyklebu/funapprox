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

/* $RCSfile: mpq_simplex.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */

#ifndef mpq___SIMPLEX_H
#define mpq___SIMPLEX_H

#include "config.h"
#include "mpq_lpdata.h"
#include "basicdefs.h"
typedef struct mpq_param_info
{
	int origalgo;
	int pphaseI;
	int pphaseII;
	int dphaseI;
	int dphaseII;
	int p_strategy;
	int d_strategy;
}
mpq_param_info;

typedef struct mpq_iter_info
{
	int newphase;
	int nextphase;
	int nextstep;
	int sdisplay;
	int itercnt;
	int solstatus;
	int curtime;
	int rounds;
	int chkobj;
	int nosolve;
	int noprog;
	int inner;
	int algorithm;
	int resumeid;
	int pricetype;
	int n_restart;
	int n_pivot_fail;
	mpq_t prevobj;
	mpq_t objtol;
	mpq_param_info oldinfo;
}
mpq_iter_info;

void mpq_ILLsimplex_init_lpinfo (mpq_lpinfo * lp),
  mpq_ILLsimplex_free_lpinfo (mpq_lpinfo * lp),
  mpq_ILLsimplex_load_lpinfo (mpq_ILLlpdata * qslp,
													mpq_lpinfo * lp),
  mpq_ILLsimplex_set_bound (mpq_lpinfo * lp,
												const mpq_t * objbound,
												int sense);
void mpq_free_internal_lpinfo (mpq_lpinfo * lp);
void mpq_init_internal_lpinfo (mpq_lpinfo * lp);
int mpq_build_internal_lpinfo (mpq_lpinfo * lp);
int mpq_ILLsimplex_retest_psolution (mpq_lpinfo * lp,
																 mpq_price_info * p,
																 int phase,
																 mpq_feas_info * fs),
  mpq_ILLsimplex_retest_dsolution (mpq_lpinfo * lp,
															 mpq_price_info * p,
															 int phase,
															 mpq_feas_info * fs),
  mpq_ILLsimplex_solution (mpq_lpinfo * lp,
											 mpq_t * xz,
											 mpq_t * piz,
											 mpq_t * dz,
											 mpq_t * objval),
  mpq_ILLsimplex_infcertificate (mpq_lpinfo * lp,
														 mpq_t * pi),
  mpq_ILLsimplex (mpq_lpinfo * lp,
							int algorithm,
							mpq_ILLlp_basis * B,
							mpq_price_info * pinf,
							int *sol_status,
							int sdisplay),
  mpq_ILLsimplex_pivotin (mpq_lpinfo * lp,
											mpq_price_info * pinf,
											int rcnt,
											int *rlist,
											int pivot_opt,
											int *basis_mod);

#endif /* mpq___SIMPLEX_H */
