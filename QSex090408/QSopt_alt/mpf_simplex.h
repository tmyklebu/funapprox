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

/* $RCSfile: mpf_simplex.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */

#ifndef mpf___SIMPLEX_H
#define mpf___SIMPLEX_H

#include "config.h"
#include "mpf_lpdata.h"
#include "basicdefs.h"
typedef struct mpf_param_info
{
	int origalgo;
	int pphaseI;
	int pphaseII;
	int dphaseI;
	int dphaseII;
	int p_strategy;
	int d_strategy;
}
mpf_param_info;

typedef struct mpf_iter_info
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
	mpf_t prevobj;
	mpf_t objtol;
	mpf_param_info oldinfo;
}
mpf_iter_info;

void mpf_ILLsimplex_init_lpinfo (mpf_lpinfo * lp),
  mpf_ILLsimplex_free_lpinfo (mpf_lpinfo * lp),
  mpf_ILLsimplex_load_lpinfo (mpf_ILLlpdata * qslp,
													mpf_lpinfo * lp),
  mpf_ILLsimplex_set_bound (mpf_lpinfo * lp,
												const mpf_t * objbound,
												int sense);
void mpf_free_internal_lpinfo (mpf_lpinfo * lp);
void mpf_init_internal_lpinfo (mpf_lpinfo * lp);
int mpf_build_internal_lpinfo (mpf_lpinfo * lp);
int mpf_ILLsimplex_retest_psolution (mpf_lpinfo * lp,
																 mpf_price_info * p,
																 int phase,
																 mpf_feas_info * fs),
  mpf_ILLsimplex_retest_dsolution (mpf_lpinfo * lp,
															 mpf_price_info * p,
															 int phase,
															 mpf_feas_info * fs),
  mpf_ILLsimplex_solution (mpf_lpinfo * lp,
											 mpf_t * xz,
											 mpf_t * piz,
											 mpf_t * dz,
											 mpf_t * objval),
  mpf_ILLsimplex_infcertificate (mpf_lpinfo * lp,
														 mpf_t * pi),
  mpf_ILLsimplex (mpf_lpinfo * lp,
							int algorithm,
							mpf_ILLlp_basis * B,
							mpf_price_info * pinf,
							int *sol_status,
							int sdisplay),
  mpf_ILLsimplex_pivotin (mpf_lpinfo * lp,
											mpf_price_info * pinf,
											int rcnt,
											int *rlist,
											int pivot_opt,
											int *basis_mod);

#endif /* mpf___SIMPLEX_H */
