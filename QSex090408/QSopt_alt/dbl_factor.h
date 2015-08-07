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

/* RCSINFO $Id: dbl_factor.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef dbl___QS_FACTOR_H_
#define dbl___QS_FACTOR_H_
#include "basicdefs.h"
#include "econfig.h"
#include "dbl_dstruct.h"

typedef char dbl_QSbool;

typedef struct dbl_uc_info
{
	int cbeg;
	int nzcnt;
	int next;
	int prev;
	int delay;
}
dbl_uc_info;

typedef struct dbl_ur_info
{
	double max;
	int rbeg;
	int nzcnt;
	int pivcnt;
	int next;
	int prev;
	int delay;
}
dbl_ur_info;

typedef struct dbl_lc_info
{
	int cbeg;
	int nzcnt;
	int c;
	int crank;
	int delay;
}
dbl_lc_info;

typedef struct dbl_lr_info
{
	int rbeg;
	int nzcnt;
	int r;
	int rrank;
	int delay;
}
dbl_lr_info;

typedef struct dbl_er_info
{
	int rbeg;
	int nzcnt;
	int r;
}
dbl_er_info;

typedef struct dbl_factor_work
{
	int max_k;
	double fzero_tol;
	double szero_tol;
	double partial_tol;
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

	double maxelem_orig;
	int nzcnt_orig;
	double maxelem_factor;
	int nzcnt_factor;
	double maxelem_cur;
	int nzcnt_cur;

	double partial_cur;

	int dim;
	int stage;
	int nstages;
	int etacnt;
	double *work_coef;
	int *work_indx;
	dbl_uc_info *uc_inf;
	dbl_ur_info *ur_inf;
	dbl_lc_info *lc_inf;
	dbl_lr_info *lr_inf;
	dbl_er_info *er_inf;
	int *ucindx;									/* row index for column data */
	int *ucrind;									/* index of column in row data */
	double *uccoef;						/* coefficient for column data */
	int *urindx;									/* col index for row data */
	int *urcind;									/* index of row in column data */
	double *urcoef;						/* coefficient for row data */
	int *lcindx;									/* row index for L data */
	double *lccoef;						/* coefficient for L row data */
	int *lrindx;									/* col index for L data */
	double *lrcoef;						/* coefficient for L col data */
	int *erindx;									/* col index for eta data */
	double *ercoef;						/* coefficient for eta data */
	int *rperm;
	int *rrank;
	int *cperm;
	int *crank;
	dbl_svector xtmp;
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

	double *dmat;
	int drows;
	int dcols;
	int dense_base;
}
dbl_factor_work;

void dbl_ILLfactor_init_factor_work (dbl_factor_work * f),
  dbl_ILLfactor_free_factor_work (dbl_factor_work * f),
  dbl_ILLfactor_ftran (dbl_factor_work * f,
									 dbl_svector * a,
									 dbl_svector * x),
  dbl_ILLfactor_ftran_update (dbl_factor_work * f,
													dbl_svector * a,
													dbl_svector * upd,
													dbl_svector * x),
  dbl_ILLfactor_btran (dbl_factor_work * f,
									 dbl_svector * a,
									 dbl_svector * x);

int dbl_ILLfactor_create_factor_work (dbl_factor_work * f,
																	int dim),
  dbl_ILLfactor_set_factor_iparam (dbl_factor_work * f,
															 int param,
															 int val),
  dbl_ILLfactor_set_factor_dparam (dbl_factor_work * f,
															 int param,
															 double val),
  dbl_ILLfactor (dbl_factor_work * f,
						 int *basis,
						 int *cbeg,
						 int *clen,
						 int *cindx,
						 double * ccoef,
						 int *p_nsing,
						 int **p_singr,
						 int **p_singc),
  dbl_ILLfactor_update (dbl_factor_work * f,
										dbl_svector * a,
										int col,
										int *p_refact);

#endif /* dbl___QS_FACTOR_H_ */
