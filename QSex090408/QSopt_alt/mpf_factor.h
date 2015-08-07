/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
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

/* RCSINFO $Id: mpf_factor.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpf___QS_FACTOR_H_
#define mpf___QS_FACTOR_H_
#include "basicdefs.h"
#include "econfig.h"
#include "mpf_dstruct.h"

typedef char mpf_QSbool;

typedef struct mpf_uc_info
{
	int cbeg;
	int nzcnt;
	int next;
	int prev;
	int delay;
}
mpf_uc_info;

typedef struct mpf_ur_info
{
	mpf_t max;
	int rbeg;
	int nzcnt;
	int pivcnt;
	int next;
	int prev;
	int delay;
}
mpf_ur_info;

typedef struct mpf_lc_info
{
	int cbeg;
	int nzcnt;
	int c;
	int crank;
	int delay;
}
mpf_lc_info;

typedef struct mpf_lr_info
{
	int rbeg;
	int nzcnt;
	int r;
	int rrank;
	int delay;
}
mpf_lr_info;

typedef struct mpf_er_info
{
	int rbeg;
	int nzcnt;
	int r;
}
mpf_er_info;

typedef struct mpf_factor_work
{
	int max_k;
	mpf_t fzero_tol;
	mpf_t szero_tol;
	mpf_t partial_tol;
	double ur_space_mul;
	double uc_space_mul;
	double lc_space_mul;
	double lr_space_mul;
	double er_space_mul;
	double grow_mul;
	int p;
	int etamax;
	double minmult;
	double maxmult;
	double updmaxmult;
	double dense_fract;
	int dense_min;

	mpf_t maxelem_orig;
	int nzcnt_orig;
	mpf_t maxelem_factor;
	int nzcnt_factor;
	mpf_t maxelem_cur;
	int nzcnt_cur;

	mpf_t partial_cur;

	int dim;
	int stage;
	int nstages;
	int etacnt;
	mpf_t *work_coef;
	int *work_indx;
	mpf_uc_info *uc_inf;
	mpf_ur_info *ur_inf;
	mpf_lc_info *lc_inf;
	mpf_lr_info *lr_inf;
	mpf_er_info *er_inf;
	int *ucindx;									/* row index for column data */
	int *ucrind;									/* index of column in row data */
	mpf_t *uccoef;						/* coefficient for column data */
	int *urindx;									/* col index for row data */
	int *urcind;									/* index of row in column data */
	mpf_t *urcoef;						/* coefficient for row data */
	int *lcindx;									/* row index for L data */
	mpf_t *lccoef;						/* coefficient for L row data */
	int *lrindx;									/* col index for L data */
	mpf_t *lrcoef;						/* coefficient for L col data */
	int *erindx;									/* col index for eta data */
	mpf_t *ercoef;						/* coefficient for eta data */
	int *rperm;
	int *rrank;
	int *cperm;
	int *crank;
	mpf_svector xtmp;
	int ur_freebeg;
	int ur_space;
	int uc_freebeg;
	int uc_space;
	int lc_freebeg;
	int lc_space;
	int lr_freebeg;
	int lr_space;
	int er_freebeg;
	int er_space;

	int *p_nsing;
	int **p_singr;
	int **p_singc;

	mpf_t *dmat;
	int drows;
	int dcols;
	int dense_base;
}
mpf_factor_work;

void mpf_ILLfactor_init_factor_work (mpf_factor_work * f),
  mpf_ILLfactor_free_factor_work (mpf_factor_work * f),
  mpf_ILLfactor_ftran (mpf_factor_work * f,
									 mpf_svector * a,
									 mpf_svector * x),
  mpf_ILLfactor_ftran_update (mpf_factor_work * f,
													mpf_svector * a,
													mpf_svector * upd,
													mpf_svector * x),
  mpf_ILLfactor_btran (mpf_factor_work * f,
									 mpf_svector * a,
									 mpf_svector * x);

int mpf_ILLfactor_create_factor_work (mpf_factor_work * f,
																	int dim),
  mpf_ILLfactor_set_factor_iparam (mpf_factor_work * f,
															 int param,
															 int val),
  mpf_ILLfactor_set_factor_dparam (mpf_factor_work * f,
															 int param,
															 mpf_t val),
  mpf_ILLfactor (mpf_factor_work * f,
						 int *basis,
						 int *cbeg,
						 int *clen,
						 int *cindx,
						 mpf_t * ccoef,
						 int *p_nsing,
						 int **p_singr,
						 int **p_singc),
  mpf_ILLfactor_update (mpf_factor_work * f,
										mpf_svector * a,
										int col,
										int *p_refact);

#endif /* mpf___QS_FACTOR_H_ */
