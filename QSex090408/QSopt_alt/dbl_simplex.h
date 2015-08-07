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

/* $RCSfile: dbl_simplex.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */

#ifndef dbl___SIMPLEX_H
#define dbl___SIMPLEX_H

#include "config.h"
#include "dbl_lpdata.h"
#include "basicdefs.h"
typedef struct dbl_param_info
{
	int origalgo;
	int pphaseI;
	int pphaseII;
	int dphaseI;
	int dphaseII;
	int p_strategy;
	int d_strategy;
}
dbl_param_info;

typedef struct dbl_iter_info
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
	double prevobj;
	double objtol;
	dbl_param_info oldinfo;
}
dbl_iter_info;

void dbl_ILLsimplex_init_lpinfo (dbl_lpinfo * lp),
  dbl_ILLsimplex_free_lpinfo (dbl_lpinfo * lp),
  dbl_ILLsimplex_load_lpinfo (dbl_ILLlpdata * qslp,
													dbl_lpinfo * lp),
  dbl_ILLsimplex_set_bound (dbl_lpinfo * lp,
												const double * objbound,
												int sense);
void dbl_free_internal_lpinfo (dbl_lpinfo * lp);
void dbl_init_internal_lpinfo (dbl_lpinfo * lp);
int dbl_build_internal_lpinfo (dbl_lpinfo * lp);
int dbl_ILLsimplex_retest_psolution (dbl_lpinfo * lp,
																 dbl_price_info * p,
																 int phase,
																 dbl_feas_info * fs),
  dbl_ILLsimplex_retest_dsolution (dbl_lpinfo * lp,
															 dbl_price_info * p,
															 int phase,
															 dbl_feas_info * fs),
  dbl_ILLsimplex_solution (dbl_lpinfo * lp,
											 double * xz,
											 double * piz,
											 double * dz,
											 double * objval),
  dbl_ILLsimplex_infcertificate (dbl_lpinfo * lp,
														 double * pi),
  dbl_ILLsimplex (dbl_lpinfo * lp,
							int algorithm,
							dbl_ILLlp_basis * B,
							dbl_price_info * pinf,
							int *sol_status,
							int sdisplay),
  dbl_ILLsimplex_pivotin (dbl_lpinfo * lp,
											dbl_price_info * pinf,
											int rcnt,
											int *rlist,
											int pivot_opt,
											int *basis_mod);

#endif /* dbl___SIMPLEX_H */
