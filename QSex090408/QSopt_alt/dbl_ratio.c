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

/* "$RCSfile: dbl_ratio.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#include "config.h"
#include "dbl_sortrus.h"
#include "stddefs.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdefs.h"
#include "dbl_ratio.h"
#include "dbl_fct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

void dbl_ILLratio_pI_test (dbl_lpinfo * lp,
      int eindex,
      int dir,
      dbl_ratio_res * rs)
{
    int i = 0, k = 0;
    int col, ecol;
    int cbnd, indx = 0;
    int tctr = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    double *pivtol = &(lp->tol->pivot_tol);
    double *dftol = &(lp->tol->id_tol);
     /* HHH */ double *t = lp->upd.t;
    double t_i, delta, y_ij, rcost, nrcost, ntmp;
    double *x, *l, *u;
     /* HHH */ dbl_EGlpNumInitVar (t_i);
    dbl_EGlpNumInitVar (delta);
    dbl_EGlpNumInitVar (y_ij);
    dbl_EGlpNumInitVar (rcost);
    dbl_EGlpNumInitVar (nrcost);
    dbl_EGlpNumInitVar (ntmp);
    dbl_EGlpNumZero (t_i);
    dbl_EGlpNumZero (y_ij);
    dbl_EGlpNumZero (delta);
    rs->lindex = -1;
    dbl_EGlpNumZero (rs->tz);
    dbl_EGlpNumZero (rs->pivotval);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    ecol = lp->nbaz[eindex];
    ILL_IFTRACE2 ("%s:%d:%d:%d:%d", __func__, eindex, dir, ecol,
	(VBOUNDED == lp->vtype[ecol]));
    if (lp->vtype[ecol] == VBOUNDED) {
	dbl_EGlpNumCopyDiff (t[0], lp->uz[ecol], lp->lz[ecol]);
	ix[0] = BBOUND;
	ILL_IFTRACE2 (":%d[%d](%la,%la,%la)\n", ix[tctr], tctr,
	    dbl_EGlpNumToLf (t[tctr]), dbl_EGlpNumToLf (lp->uz[ecol]),
	    dbl_EGlpNumToLf (lp->lz[ecol]));
	tctr++;
    }
    ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	dbl_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	if (dbl_EGlpNumIsEqual (y_ij, dbl_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && dbl_EGlpNumIsLess (dbl_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum))) {
	    if (dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum))
		dbl_EGlpNumSign (y_ij);
	    ILL_IFTRACE2 (":%d", lp->bfeas[i]);
	    if (lp->bfeas[i] > 0) {
		dbl_EGlpNumCopyDiffRatio (t[tctr], *x, *u, y_ij);
		ix[tctr] = 10 * k + BATOUPPER;
		ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, dbl_EGlpNumToLf (t[tctr]));
		tctr++;
		if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY)) {
		    dbl_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
		    ix[tctr] = 10 * k + BATOLOWER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			dbl_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    } else if (lp->bfeas[i] == 0) {
		if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY)) {
		    dbl_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
		    ix[tctr] = 10 * k + BATOLOWER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			dbl_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    }
	} else if ((dir == VINCREASE && dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLess (dbl_zeroLpNum, y_ij))) {
	    if (dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum))
		dbl_EGlpNumSign (y_ij);
	    ILL_IFTRACE2 (":%d", lp->bfeas[i]);
	    if (lp->bfeas[i] < 0) {
		dbl_EGlpNumCopyDiffRatio (t[tctr], *l, *x, y_ij);
		ix[tctr] = 10 * k + BBTOLOWER;
		ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, dbl_EGlpNumToLf (t[tctr]));
		tctr++;
		if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY)) {
		    dbl_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			dbl_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    } else if (lp->bfeas[i] == 0) {
		if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY)) {
		    dbl_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			dbl_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    }
	}
    }
    if (tctr == 0) {
	rs->ratio_stat = RATIO_FAILED;
	ILL_CLEANUP;
    }
    for (i = 0; i < tctr; i++)
	perm[i] = i;
    dbl_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

    dbl_EGlpNumZero (lp->upd.c_obj);
    dbl_EGlpNumCopy (rcost, lp->pIdz[eindex]);
    ILL_IFTRACE2 ("\n%s:%d:%lf", __func__, tctr, dbl_EGlpNumToLf (rcost));
    for (i = 0; i < tctr; i++) {
	dbl_EGlpNumCopy (t_i, t[perm[i]]);
	dbl_EGlpNumCopy (ntmp, t_i);
	dbl_EGlpNumSubTo (ntmp, delta);
	dbl_EGlpNumAddInnProdTo (lp->upd.c_obj, ntmp, rcost);
	dbl_EGlpNumCopy (delta, t_i);
	ILL_IFTRACE2 (":%d:%lf", perm[i], dbl_EGlpNumToLf (delta));
	 /* HHH */ cbnd = ix[perm[i]] % 10;
	if (cbnd != BBOUND) {
	    k = ix[perm[i]] / 10;
	    dbl_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	    indx = lp->yjz.indx[k];
	    ILL_IFTRACE2 (":%d", indx);
	}
	switch (cbnd) {
	case BBOUND:
	    rs->ratio_stat = RATIO_NOBCHANGE;
	    dbl_EGlpNumCopy (rs->tz, t_i);
	    if (dir != VINCREASE)
		dbl_EGlpNumSign (rs->tz);
	    ILL_CLEANUP;

	case BATOLOWER:
	case BATOUPPER:
	    dbl_EGlpNumAddTo (rcost, y_ij);
	    break;
	case BBTOLOWER:
	case BBTOUPPER:
	    dbl_EGlpNumSubTo (rcost, y_ij);
	    break;
	}
	dbl_EGlpNumCopyNeg (nrcost, rcost);
	if ((dir == VINCREASE && dbl_EGlpNumIsLeq (nrcost, *dftol)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLeq (rcost, *dftol))) {
	    /* change 5 to -1 if t_i > 0 is required below */
	    if (dbl_EGlpNumIsLess (t_i, dbl_zeroLpNum) && i > 5) {
		/* printf ("pIhell %.5f %d\n", t_i, i); */
		dbl_EGlpNumDivUiTo (t_i, 2);
		rs->ratio_stat = RATIO_NEGATIVE;
		dbl_EGlpNumZero (rs->tz);
		ILL_CLEANUP;
	    }
	    rs->lindex = indx;
	    rs->ratio_stat = RATIO_BCHANGE;
	    if (cbnd == BATOLOWER || cbnd == BBTOLOWER)
		rs->lvstat = STAT_LOWER;
	    else
		rs->lvstat = STAT_UPPER;

	    dbl_EGlpNumCopy (rs->pivotval, y_ij);
	    dbl_EGlpNumCopy (rs->tz, t_i);
	    if (dir != VINCREASE)
		dbl_EGlpNumSign (rs->tz);
	    ILL_CLEANUP;
	}
    }

CLEANUP:
    dbl_ILLfct_update_counts (lp, CNT_PIPIV, 0, rs->pivotval);
    ILL_IFTRACE2 (":tctr %d:%d\n", tctr, rs->ratio_stat);
    lp->upd.tctr = tctr;
    lp->upd.i = i;
    dbl_EGlpNumCopy (lp->upd.tz, t_i);
    dbl_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    if (dir == VDECREASE)
	dbl_EGlpNumSign (lp->upd.c_obj);
    if (rs->lindex != -1)
	lp->upd.fs = lp->bfeas[rs->lindex];
    dbl_EGlpNumClearVar (t_i);
    dbl_EGlpNumClearVar (delta);
    dbl_EGlpNumClearVar (y_ij);
    dbl_EGlpNumClearVar (rcost);
    dbl_EGlpNumClearVar (nrcost);
    dbl_EGlpNumClearVar (ntmp);
}

void dbl_ILLratio_pII_test (dbl_lpinfo * lp,
      int eindex,
      int dir,
      dbl_ratio_res * rs)
{
    int i, k, indx, col, ecol;
    double *x, *l, *u, t_max, ayi_max, yi_max, ay_ij, y_ij, t_i, t_z;
    double *pivtol = &(lp->tol->pivot_tol);
    double *pftol = &(lp->tol->pfeas_tol);
    dbl_EGlpNumInitVar (y_ij);
    dbl_EGlpNumInitVar (ay_ij);
    dbl_EGlpNumInitVar (t_i);
    dbl_EGlpNumInitVar (t_z);
    dbl_EGlpNumInitVar (t_max);
    dbl_EGlpNumInitVar (yi_max);
    dbl_EGlpNumInitVar (ayi_max);
     /* HHH */ rs->boundch = 0;
    rs->lindex = -1;
    dbl_EGlpNumZero (rs->tz);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    dbl_EGlpNumZero (rs->pivotval);
    dbl_EGlpNumZero (rs->lbound);
    ecol = lp->nbaz[eindex];

    for (k = 0, dbl_EGlpNumCopy (t_max, dbl_INFTY); k < lp->yjz.nzcnt; k++) {
	dbl_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	dbl_EGlpNumCopyAbs (ay_ij, y_ij);
	if (dbl_EGlpNumIsEqual (y_ij, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_i, dbl_INFTY);
	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && dbl_EGlpNumIsLess (dbl_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum))) {
	    if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY)) {
		dbl_EGlpNumCopyDiff (t_i, *x, *l);
		dbl_EGlpNumAddTo (t_i, *pftol);
		dbl_EGlpNumDivTo (t_i, ay_ij);
	    }
	} else if ((dir == VINCREASE && dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLess (dbl_zeroLpNum, y_ij))) {
	    if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY)) {
		dbl_EGlpNumCopySum (t_i, *u, *pftol);
		dbl_EGlpNumSubTo (t_i, *x);
		dbl_EGlpNumDivTo (t_i, ay_ij);
	    }
	}
	if (dbl_EGlpNumIsEqqual (t_i, dbl_INFTY))
	    continue;

	if (dbl_EGlpNumIsLess (t_i, t_max)) {
	    /* HHH tind = i; yval = fabs (y_ij); tval = t_i -
	       pftol/fabs(y_ij); */
	    dbl_EGlpNumCopy (t_max, t_i);
	}
    }
    /* we use yi_max as temporal variable here */
    dbl_EGlpNumCopyDiff (yi_max, lp->uz[ecol], lp->lz[ecol]);
    if (lp->vtype[ecol] == VBOUNDED && dbl_EGlpNumIsLeq (yi_max, t_max)) {

	dbl_EGlpNumCopy (t_max, yi_max);
	rs->ratio_stat = RATIO_NOBCHANGE;
	dbl_EGlpNumCopy (rs->tz, t_max);
	if (dir != VINCREASE)
	    dbl_EGlpNumSign (rs->tz);
	ILL_CLEANUP;
    }
    if (dbl_EGlpNumIsLeq (dbl_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    /* if (dbl_EGlpNumIsLess (t_max, dbl_zeroLpNum)) printf ("pIIhell\n"); */
    indx = -1;
    dbl_EGlpNumZero (t_z);
    dbl_EGlpNumZero (yi_max);
    dbl_EGlpNumZero (ayi_max);
    ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	dbl_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	dbl_EGlpNumCopyAbs (ay_ij, y_ij);
	if (dbl_EGlpNumIsEqual (y_ij, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_i, dbl_INFTY);
	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && dbl_EGlpNumIsLess (dbl_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum))) {
	    if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY))
		dbl_EGlpNumCopyDiffRatio (t_i, *x, *l, ay_ij);
	} else if ((dir == VINCREASE && dbl_EGlpNumIsLess (y_ij, dbl_zeroLpNum)) ||
	    (dir == VDECREASE && dbl_EGlpNumIsLess (dbl_zeroLpNum, y_ij))) {
	    if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY))
		dbl_EGlpNumCopyDiffRatio (t_i, *u, *x, ay_ij);
	}
	if (dbl_EGlpNumIsLeq (t_i, t_max)) {
	    if (dbl_EGlpNumIsLess (ayi_max, ay_ij)) {
		dbl_EGlpNumCopy (yi_max, y_ij);
		dbl_EGlpNumCopy (ayi_max, ay_ij);
		indx = i;
		dbl_EGlpNumCopy (t_z, t_i);
		ILL_IFTRACE2 (":%d:%lf:%lf:%lf:%lf", indx, dbl_EGlpNumToLf (t_i),
		    dbl_EGlpNumToLf (t_max), dbl_EGlpNumToLf (ayi_max),
		    dbl_EGlpNumToLf (ay_ij));
	    }
	}
    }

    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	/*
         * if (tind != rs->lindex){
         * HHHprintf ("tmax %e tval = %e yval = %e tind = %d\n", t_max, tval, yval, tind);
         * HHHprintf ("h tval = %e yval = %e tind = %d\n",rs->tz, yi_max, rs->lindex);
         * }
         */
	ILL_IFTRACE2 (":%d", indx);
	rs->lindex = indx;
	dbl_EGlpNumCopy (rs->tz, t_z);
	dbl_EGlpNumCopy (rs->pivotval, yi_max);
	rs->ratio_stat = RATIO_BCHANGE;

	if (dir == VINCREASE)
	    rs->lvstat =
		(dbl_EGlpNumIsLess (dbl_zeroLpNum, yi_max)) ? STAT_LOWER : STAT_UPPER;
	else
	    rs->lvstat =
		(dbl_EGlpNumIsLess (dbl_zeroLpNum, yi_max)) ? STAT_UPPER : STAT_LOWER;

	if (dbl_EGlpNumIsLess (rs->tz, dbl_zeroLpNum)) {
	    ILL_IFTRACE2 ("need to change bound, tz=%la\n", dbl_EGlpNumToLf (rs->tz));
	    dbl_EGlpNumCopyAbs (rs->tz, t_max);
	    dbl_EGlpNumDivUiTo (rs->tz, 10);
	    rs->boundch = 1;
	    dbl_EGlpNumCopy (rs->lbound, lp->xbz[rs->lindex]);
	    if (rs->lvstat == STAT_LOWER)
		dbl_EGlpNumSubInnProdTo (rs->lbound, rs->tz, ayi_max);
	    else
		dbl_EGlpNumAddInnProdTo (rs->lbound, rs->tz, ayi_max);
	}
	if (dir == VDECREASE)
	    dbl_EGlpNumSign (rs->tz);
    }
CLEANUP:
    dbl_ILLfct_update_counts (lp, CNT_PIIPIV, 0, rs->pivotval);
    dbl_EGlpNumClearVar (y_ij);
    dbl_EGlpNumClearVar (ay_ij);
    dbl_EGlpNumClearVar (t_i);
    dbl_EGlpNumClearVar (t_z);
    dbl_EGlpNumClearVar (t_max);
    dbl_EGlpNumClearVar (yi_max);
    dbl_EGlpNumClearVar (ayi_max);
}

#define dbl_GET_XY_DRATIOTEST \
      if (lp->vstat[col] == STAT_UPPER){ \
                dbl_EGlpNumCopyNeg(x,lp->dz[j]);\
        dbl_EGlpNumCopy(y, *zAj);\
      } \
      else{ \
         dbl_EGlpNumCopy(x, lp->dz[j]); \
         dbl_EGlpNumCopyNeg(y, *zAj);\
      } \
      if (lvstat == STAT_UPPER) \
         dbl_EGlpNumSign(y);


void dbl_ILLratio_dI_test (dbl_lpinfo * lp,
      int lindex,
      int lvstat,
      dbl_ratio_res * rs)
{
    int j = 0, k;
    int col;
    int cbnd, indx;
    int tctr = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    double *t = lp->upd.t;
    double *zAj, x, y, t_j, theta, rcost, delta;
    double *pftol = &(lp->tol->ip_tol);
    double *pivtol = &(lp->tol->pivot_tol);
    dbl_EGlpNumInitVar (x);
    dbl_EGlpNumInitVar (y);
    dbl_EGlpNumInitVar (t_j);
    dbl_EGlpNumInitVar (theta);
    dbl_EGlpNumInitVar (rcost);
    dbl_EGlpNumInitVar (delta);
    dbl_EGlpNumZero (delta);
    dbl_EGlpNumZero (t_j);
    dbl_EGlpNumZero (rs->tz);
     /* HHH */ rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;
    dbl_EGlpNumZero (rs->pivotval);

    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (dbl_EGlpNumIsEqual (*zAj, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_j, dbl_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	dbl_GET_XY_DRATIOTEST;

	if (dbl_EGlpNumIsLess (y, dbl_zeroLpNum)) {
	    if (lp->dfeas[j] != 0 && lp->vstat[col] != STAT_ZERO) {
		dbl_EGlpNumCopyFrac (t[tctr], x, y);
		ix[tctr] = 10 * k + BBTOLOWER;
		tctr++;
	    } else if (lp->vstat[col] == STAT_ZERO) {
		if (lp->dfeas[j] < 0) {
		    dbl_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BBTOLOWER;
		    tctr++;
		}
		if (lp->dfeas[j] <= 0) {
		    dbl_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    tctr++;
		}
	    }
	} else {
	    if (lp->dfeas[j] > 0) {
		if (lp->vstat[col] == STAT_ZERO) {
		    dbl_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BATOUPPER;
		    tctr++;
		    dbl_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BATOLOWER;
		    tctr++;
		}
	    } else if (lp->dfeas[j] == 0) {
		dbl_EGlpNumCopyFrac (t[tctr], x, y);
		if (lp->vtype[col] == VBOUNDED)
		    ix[tctr] = 10 * k + BSKIP;
		else
		    ix[tctr] = 10 * k + BATOLOWER;
		tctr++;
	    }
	}
    }

    if (tctr == 0) {
	rs->ratio_stat = RATIO_FAILED;
	ILL_CLEANUP;
    }
    for (j = 0; j < tctr; j++)
	perm[j] = j;
    dbl_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

    dbl_EGlpNumZero (lp->upd.c_obj);
    dbl_EGlpNumCopy (rcost, lp->xbz[lindex]);
    if (lvstat == STAT_LOWER)
	dbl_EGlpNumSign (rcost);
    for (j = 0; j < tctr; j++) {
	cbnd = ix[perm[j]] % 10;
	if (cbnd == BSKIP)
	    continue;

	dbl_EGlpNumCopy (t_j, t[perm[j]]);
	dbl_EGlpNumCopy (x, t_j);
	dbl_EGlpNumSubTo (x, delta);
	dbl_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	dbl_EGlpNumCopy (delta, t_j);
	k = ix[perm[j]] / 10;
	zAj = &(lp->zA.coef[k]);
	indx = lp->zA.indx[k];

	if (lp->vstat[lp->nbaz[indx]] == STAT_LOWER
	    || lp->vstat[lp->nbaz[indx]] == STAT_ZERO)
	    dbl_EGlpNumCopyNeg (theta, *zAj);
	else
	    dbl_EGlpNumCopy (theta, *zAj);

	if (lvstat == STAT_UPPER)
	    dbl_EGlpNumSign (theta);

	switch (cbnd) {
	case BATOLOWER:
	case BATOUPPER:
	    dbl_EGlpNumSubTo (rcost, theta);
	    break;
	case BBTOLOWER:
	case BBTOUPPER:
	    dbl_EGlpNumAddTo (rcost, theta);
	    break;
	}
	if (dbl_EGlpNumIsLeq (rcost, *pftol)) {
	    /* if (t_j < 0.0) printf ("dIhell\n"); */
	    rs->eindex = indx;
	    dbl_EGlpNumCopy (rs->tz, t_j);
	    dbl_EGlpNumCopy (rs->pivotval, *zAj);
	    rs->ratio_stat = RATIO_BCHANGE;
	    ILL_CLEANUP;
	}
    }

CLEANUP:
    dbl_ILLfct_update_counts (lp, CNT_DIPIV, 0, rs->pivotval);
    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
    lp->upd.tctr = tctr;
    lp->upd.i = j;
    dbl_EGlpNumCopyAbs (lp->upd.tz, t_j);
    dbl_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    if (rs->eindex != -1)
	lp->upd.fs = lp->dfeas[rs->eindex];
    dbl_EGlpNumClearVar (x);
    dbl_EGlpNumClearVar (y);
    dbl_EGlpNumClearVar (t_j);
    dbl_EGlpNumClearVar (theta);
    dbl_EGlpNumClearVar (rcost);
    dbl_EGlpNumClearVar (delta);
}

void dbl_ILLratio_dII_test (dbl_lpinfo * lp, int lvstat, dbl_ratio_res * rs)
{
    int j, k, indx;
    int col, ecol;
    double *zAj, azAj, az_max, x, y, t_j, z_max, t_max, t_z;
    double *dftol = &(lp->tol->dfeas_tol);
    double *pivtol = &(lp->tol->pivot_tol);
    dbl_EGlpNumInitVar (x);
    dbl_EGlpNumInitVar (y);
    dbl_EGlpNumInitVar (t_j);
    dbl_EGlpNumInitVar (z_max);
    dbl_EGlpNumInitVar (t_max);
    dbl_EGlpNumInitVar (az_max);
    dbl_EGlpNumInitVar (azAj);
    dbl_EGlpNumInitVar (t_z);
    dbl_EGlpNumZero (t_j);
    rs->coeffch = 0;
    dbl_EGlpNumZero (rs->ecoeff);
    rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;
    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
    lp->upd.tctr = 0;
    dbl_EGlpNumZero (lp->upd.dty);
    for (k = 0, dbl_EGlpNumCopy (t_max, dbl_INFTY); k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (dbl_EGlpNumIsEqual (*zAj, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_j, dbl_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	dbl_GET_XY_DRATIOTEST;

	/* #warning adding/substracting tolerances to used value, is it
	   rigght? */
	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, y)) {
	    /* t_j = (x + dftol) / y; */
	    dbl_EGlpNumCopySum (t_j, x, *dftol);
	    dbl_EGlpNumDivTo (t_j, y);
	} else {
	    /* #warning adding/substracting tolerances to used value, is it
	       rigght? */
	    if (lp->vstat[col] == STAT_ZERO)
		dbl_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
	}
	/* if (t_j == dbl_INFTY) */
	if (dbl_EGlpNumIsEqqual (t_j, dbl_INFTY))
	    continue;

	if (dbl_EGlpNumIsLess (t_j, t_max))
	    dbl_EGlpNumCopy (t_max, t_j);
    }

    if (dbl_EGlpNumIsLeq (dbl_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    /* if (t_max < 0.0) printf ("dIIhell\n"); */

    indx = -1;
    dbl_EGlpNumZero (t_z);
    dbl_EGlpNumZero (z_max);
    dbl_EGlpNumZero (az_max);

    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	dbl_EGlpNumCopyAbs (azAj, *zAj);
	if (dbl_EGlpNumIsEqual (*zAj, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_j, dbl_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	dbl_GET_XY_DRATIOTEST;

	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, y) || lp->vstat[col] == STAT_ZERO)
	    dbl_EGlpNumCopyFrac (t_j, x, y);

	if (dbl_EGlpNumIsLeq (t_j, t_max) && (dbl_EGlpNumIsLess (az_max, azAj))) {
	    dbl_EGlpNumCopy (z_max, *zAj);
	    dbl_EGlpNumCopy (az_max, azAj);
	    indx = j;
	    dbl_EGlpNumCopy (t_z, t_j);
	}
    }


    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	rs->eindex = indx;
	dbl_EGlpNumCopy (rs->tz, t_z);
	dbl_EGlpNumCopy (rs->pivotval, z_max);
	rs->ratio_stat = RATIO_BCHANGE;

	if (dbl_EGlpNumIsLess (rs->tz, dbl_zeroLpNum)) {
	    dbl_EGlpNumCopyAbs (rs->tz, t_max);
	    dbl_EGlpNumDivUiTo (rs->tz, 20);
	    rs->coeffch = 1;
	    ecol = lp->nbaz[indx];
	    dbl_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
	    switch (lp->vstat[ecol]) {
	    case STAT_LOWER:
		dbl_EGlpNumAddInnProdTo (rs->ecoeff, rs->tz, az_max);
		break;
	    case STAT_UPPER:
		dbl_EGlpNumSubInnProdTo (rs->ecoeff, rs->tz, az_max);
		break;
	    default:
		dbl_EGlpNumZero (rs->tz);
		break;
	    }
	}
    }

CLEANUP:
    dbl_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
    dbl_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    dbl_EGlpNumClearVar (x);
    dbl_EGlpNumClearVar (y);
    dbl_EGlpNumClearVar (t_j);
    dbl_EGlpNumClearVar (z_max);
    dbl_EGlpNumClearVar (t_max);
    dbl_EGlpNumClearVar (t_z);
    dbl_EGlpNumClearVar (az_max);
    dbl_EGlpNumClearVar (azAj);
}

void dbl_ILLratio_longdII_test (dbl_lpinfo * lp,
      int lindex,
      int lvstat,
      dbl_ratio_res * rs)
{
    int j, k, indx = 0, tctr = 0;
    int col, ecol;
    int vs, bnd_exist = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int b_indx = -1;
    double *t = lp->upd.t;
    double *l, *u, *xb, *zAj = 0, x, y, t_j, z_max, t_max, t_z, theta,
      rcost, delta, zb_val, tb_val, az_max, azb_val, azAj;
    double *pftol = &(lp->tol->pfeas_tol);
    double *dftol = &(lp->tol->dfeas_tol);
    double *pivtol = &(lp->tol->pivot_tol);
    dbl_EGlpNumInitVar (x);
    dbl_EGlpNumInitVar (azAj);
    dbl_EGlpNumInitVar (y);
    dbl_EGlpNumInitVar (t_j);
    dbl_EGlpNumInitVar (z_max);
    dbl_EGlpNumInitVar (az_max);
    dbl_EGlpNumInitVar (t_max);
    dbl_EGlpNumInitVar (t_z);
    dbl_EGlpNumInitVar (theta);
    dbl_EGlpNumInitVar (rcost);
    dbl_EGlpNumInitVar (delta);
    dbl_EGlpNumInitVar (zb_val);
    dbl_EGlpNumInitVar (azb_val);
    dbl_EGlpNumInitVar (tb_val);
    dbl_EGlpNumZero (t_j);
    dbl_EGlpNumZero (delta);
    dbl_EGlpNumZero (zb_val);
    dbl_EGlpNumZero (azb_val);
    dbl_EGlpNumCopy (tb_val, dbl_NINFTY);
    /* #warning not sure about THIS line */
    dbl_EGlpNumZero (rs->pivotval);

    rs->coeffch = 0;
    rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;

    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
    lp->upd.tctr = 0;
    lp->upd.i = 0;
    dbl_EGlpNumZero (lp->upd.tz);
    dbl_EGlpNumZero (lp->upd.piv);
    dbl_EGlpNumZero (lp->upd.c_obj);
    dbl_EGlpNumZero (lp->upd.dty);

    xb = &(lp->xbz[lindex]);
    col = lp->baz[lindex];
    l = &(lp->lz[col]);
    u = &(lp->uz[col]);
    /* rcost = (lvstat == STAT_LOWER) ? l - xb : xb - u; */
    if (lvstat == STAT_LOWER)
	dbl_EGlpNumCopyDiff (rcost, *l, *xb);
    else
	dbl_EGlpNumCopyDiff (rcost, *xb, *u);

    for (k = 0, dbl_EGlpNumCopy (t_max, dbl_INFTY); k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (dbl_EGlpNumIsEqual (*zAj, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_j, dbl_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;
	if (lp->vtype[col] == VBOUNDED) {
	    bnd_exist++;
	    continue;
	}
	dbl_GET_XY_DRATIOTEST;

	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, y)) {
	    /* t_j = (x + dftol) / y; */
	    /* #warning Using tolerances to add to result, is it right? */
	    dbl_EGlpNumCopySum (t_j, x, *dftol);
	    dbl_EGlpNumDivTo (t_j, y);
	} else {
	    if (lp->vstat[col] == STAT_ZERO)
		dbl_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
	}
	if (dbl_EGlpNumIsEqqual (t_j, dbl_INFTY))
	    continue;

	if (dbl_EGlpNumIsLess (t_j, t_max))
	    dbl_EGlpNumCopy (t_max, t_j);
    }
    if (dbl_EGlpNumIsLess (t_max, dbl_zeroLpNum)) {
	/* printf ("dIIhell, %.4f\n", t_max); */
	rs->ratio_stat = RATIO_NEGATIVE;
	ILL_CLEANUP;
    }
    if (bnd_exist == 0 && dbl_EGlpNumIsLeq (dbl_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	/*
         * printf ("x = %.8f, b = %.2f \n", lp->xbz[lindex], (lvstat == STAT_LOWER ) ? lp->lz[lp->baz[lindex]] : lp->uz[lp->baz[lindex]]);
         */
	ILL_CLEANUP;
    }
    if (bnd_exist != 0) {
	for (k = 0; k < lp->zA.nzcnt; k++) {
	    zAj = &(lp->zA.coef[k]);
	    if (dbl_EGlpNumIsEqual (*zAj, dbl_zeroLpNum, *pivtol))
		continue;

	    dbl_EGlpNumCopy (t_j, dbl_INFTY);
	    j = lp->zA.indx[k];
	    col = lp->nbaz[j];

	    if (lp->vtype[col] != VBOUNDED)
		continue;

	    dbl_GET_XY_DRATIOTEST;

	    if (dbl_EGlpNumIsLess (dbl_zeroLpNum, y)) {
		dbl_EGlpNumCopyFrac (t_j, x, y);
		if (dbl_EGlpNumIsLeq (t_j, t_max)) {
		    dbl_EGlpNumCopy (t[tctr], t_j);
		    ix[tctr] = k;
		    tctr++;
		}
	    }
	}
    }
    if (tctr != 0) {
	for (j = 0; j < tctr; j++)
	    perm[j] = j;
	dbl_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	for (j = 0; j < tctr; j++) {

	    dbl_EGlpNumCopy (t_j, t[perm[j]]);
	    /* we use x as temporal storage */
	    /* lp->upd.c_obj += (t_j - delta) * rcost; */
	    dbl_EGlpNumCopy (x, t_j);
	    dbl_EGlpNumSubTo (x, delta);
	    dbl_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	    dbl_EGlpNumCopy (delta, t_j);
	     /* HHH */ k = ix[perm[j]];
	    zAj = &(lp->zA.coef[k]);
	    indx = lp->zA.indx[k];
	    col = lp->nbaz[indx];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);
	    vs = lp->vstat[col];
	    /* theta = (vs == STAT_UPPER) ? (l - u) * zAj : (u - l) * zAj; */
	    dbl_EGlpNumCopyDiff (theta, *l, *u);
	    dbl_EGlpNumMultTo (theta, *zAj);
	    if (vs != STAT_UPPER)
		dbl_EGlpNumSign (theta);
	    if (lvstat == STAT_LOWER)
		dbl_EGlpNumAddTo (rcost, theta);
	    else
		dbl_EGlpNumSubTo (rcost, theta);

	    if (dbl_EGlpNumIsLeq (rcost, *pftol)) {
		rs->eindex = indx;
		dbl_EGlpNumCopy (rs->tz, t_j);
		dbl_EGlpNumCopy (rs->pivotval, *zAj);
		rs->ratio_stat = RATIO_BCHANGE;

		if (dbl_EGlpNumIsLess (rs->tz, dbl_zeroLpNum)) {
		    dbl_EGlpNumZero (rs->tz);
		    rs->coeffch = 1;
		    /* rs->ecoeff = lp->cz[col] - lp->dz[indx]; */
		    dbl_EGlpNumCopyDiff (rs->ecoeff, lp->cz[col], lp->dz[indx]);
		    /* lp->upd.c_obj += (rs->tz - delta) * rcost; note ts->tz
		       == 0; */
		    dbl_EGlpNumSubInnProdTo (lp->upd.c_obj, delta, rcost);
		}
		ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
		lp->upd.tctr = tctr;
		lp->upd.i = j;
		dbl_EGlpNumCopy (lp->upd.tz, rs->tz);
		ILL_CLEANUP;
	    }
	}
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
	lp->upd.tctr = tctr;
	lp->upd.i = tctr;
	dbl_EGlpNumCopy (lp->upd.tz, t_j);
	dbl_EGlpNumCopy (zb_val, *zAj);
	dbl_EGlpNumCopyAbs (azb_val, zb_val);
	dbl_EGlpNumCopy (tb_val, t_j);
	b_indx = indx;
    }
    if (bnd_exist != 0 && dbl_EGlpNumIsLeq (dbl_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	/* printf ("rcost: %.8f\n", rcost); */
	ILL_CLEANUP;
    }
    dbl_EGlpNumZero (z_max);
    dbl_EGlpNumZero (az_max);
    indx = -1;
    dbl_EGlpNumZero (t_z);
    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	dbl_EGlpNumCopyAbs (azAj, *zAj);
	if (dbl_EGlpNumIsEqual (*zAj, dbl_zeroLpNum, *pivtol))
	    continue;

	dbl_EGlpNumCopy (t_j, dbl_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED ||
	    lp->vtype[col] == VBOUNDED)
	    continue;

	dbl_GET_XY_DRATIOTEST;

	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, y) || lp->vstat[col] == STAT_ZERO)
	    dbl_EGlpNumCopyFrac (t_j, x, y);

	if (dbl_EGlpNumIsLeq (t_j, t_max)) {
	    if (dbl_EGlpNumIsLess (az_max, azAj)) {
		dbl_EGlpNumCopy (z_max, *zAj);
		dbl_EGlpNumCopy (az_max, azAj);
		indx = j;
		dbl_EGlpNumCopy (t_z, t_j);
	    }
	}
    }

    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
	ILL_CLEANUP;
    }
    if ((tctr == 0) || (dbl_EGlpNumIsLess (tb_val, dbl_zeroLpNum)) ||
	(tctr != 0 && dbl_EGlpNumIsLeq (tb_val, t_z) &&
	    dbl_EGlpNumIsLeq (azb_val, az_max))) {
	/* we use x as temporal vvariable */
	/* lp->upd.c_obj += (t_z - delta) * rcost; */
	dbl_EGlpNumCopyDiff (x, t_z, delta);
	dbl_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	dbl_EGlpNumCopy (delta, t_z);
	rs->eindex = indx;
	dbl_EGlpNumCopy (rs->tz, t_z);
	dbl_EGlpNumCopy (rs->pivotval, z_max);
	rs->ratio_stat = RATIO_BCHANGE;
    }
    /* For now */
    else if (tctr != 0) {
	rs->eindex = b_indx;
	dbl_EGlpNumCopy (rs->tz, tb_val);
	dbl_EGlpNumCopy (rs->pivotval, zb_val);
	rs->ratio_stat = RATIO_BCHANGE;
	lp->upd.i -= 1;
    }
    if (dbl_EGlpNumIsLess (rs->tz, dbl_zeroLpNum)) {
	/* if (tctr != 0) printf ("despite long step\n"); */
	/* rs->tz = fabs (t_max / 20.0); */
	dbl_EGlpNumCopyAbs (rs->tz, t_max);
	dbl_EGlpNumDivUiTo (rs->tz, 20);
	rs->coeffch = 1;

	ecol = lp->nbaz[indx];
	if (lp->vstat[ecol] == STAT_LOWER) {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs
	       (z_max); */
	    dbl_EGlpNumCopy (rs->ecoeff, az_max);
	    dbl_EGlpNumMultTo (rs->ecoeff, rs->tz);
	    dbl_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
	    dbl_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
	} else if (lp->vstat[ecol] == STAT_UPPER) {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs
	       (z_max); */
	    dbl_EGlpNumCopy (rs->ecoeff, az_max);
	    dbl_EGlpNumMultTo (rs->ecoeff, rs->tz);
	    dbl_EGlpNumSign (rs->ecoeff);
	    dbl_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
	    dbl_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
	} else {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx]; */
	    dbl_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
	    dbl_EGlpNumZero (rs->tz);
	}
	/* we use x as temporal storage */
	/* lp->upd.c_obj += (rs->tz - delta) * rcost; */
	dbl_EGlpNumCopy (x, rs->tz);
	dbl_EGlpNumSubTo (x, delta);
	dbl_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
    }
CLEANUP:
    dbl_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
    dbl_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    dbl_EGlpNumClearVar (x);
    dbl_EGlpNumClearVar (y);
    dbl_EGlpNumClearVar (t_j);
    dbl_EGlpNumClearVar (z_max);
    dbl_EGlpNumClearVar (az_max);
    dbl_EGlpNumClearVar (t_max);
    dbl_EGlpNumClearVar (t_z);
    dbl_EGlpNumClearVar (theta);
    dbl_EGlpNumClearVar (rcost);
    dbl_EGlpNumClearVar (delta);
    dbl_EGlpNumClearVar (zb_val);
    dbl_EGlpNumClearVar (azb_val);
    dbl_EGlpNumClearVar (tb_val);
    dbl_EGlpNumClearVar (azAj);
}

void dbl_ILLratio_pivotin_test (dbl_lpinfo * lp,
      int *rlist,
      int rcnt,
      dbl_ratio_res * rs)
{
    int i, k, col;
    double *x, *l, *u;
    double ay_ij, at_i, at_l, at_u, ayi_max, y_ij, t_i, t_l, t_u, t_max,
      yi_max;
    double *pivtol = &(lp->tol->pivot_tol);
    if (rcnt <= 0 || rs == NULL)
	return;
    dbl_EGlpNumInitVar (ay_ij);
    dbl_EGlpNumInitVar (at_i);
    dbl_EGlpNumInitVar (at_l);
    dbl_EGlpNumInitVar (at_u);
    dbl_EGlpNumInitVar (ayi_max);
    dbl_EGlpNumInitVar (t_max);
    dbl_EGlpNumInitVar (y_ij);
    dbl_EGlpNumInitVar (t_i);
    dbl_EGlpNumInitVar (t_l);
    dbl_EGlpNumInitVar (t_u);
    dbl_EGlpNumInitVar (yi_max);
    rs->boundch = 0;
    rs->lindex = -1;
    dbl_EGlpNumZero (rs->tz);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    dbl_EGlpNumZero (rs->pivotval);
    dbl_EGlpNumZero (rs->lbound);

    for (i = 0; i < rcnt; i++)
	lp->iwork[rlist[i]] = 1;

    for (k = 0, dbl_EGlpNumCopy (t_max, dbl_INFTY); k < lp->yjz.nzcnt; k++) {
	dbl_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	if (dbl_EGlpNumIsEqual (y_ij, dbl_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	if (lp->iwork[lp->baz[i]] == 1)
	    continue;
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	dbl_EGlpNumCopy (t_u, dbl_INFTY);
	dbl_EGlpNumCopy (at_u, dbl_INFTY);
	dbl_EGlpNumCopy (t_l, dbl_NINFTY);
	dbl_EGlpNumCopy (at_l, dbl_INFTY);

	if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY)) {
	    dbl_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
	    dbl_EGlpNumCopyAbs (at_l, t_l);
	    if (dbl_EGlpNumIsLess (at_l, t_max))
		dbl_EGlpNumCopy (t_max, at_l);
	}
	if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY)) {
	    dbl_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
	    dbl_EGlpNumCopyAbs (at_u, t_u);
	    if (dbl_EGlpNumIsLess (at_u, t_max))
		dbl_EGlpNumCopy (t_max, at_u);
	}
    }

    if (dbl_EGlpNumIsLeq (dbl_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    dbl_EGlpNumZero (yi_max);
    dbl_EGlpNumZero (ayi_max);
    dbl_EGlpNumMultUiTo (t_max, 101);
    dbl_EGlpNumDivUiTo (t_max, 100);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	dbl_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	dbl_EGlpNumCopyAbs (ay_ij, y_ij);
	if (dbl_EGlpNumIsEqual (y_ij, dbl_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	if (lp->iwork[lp->baz[i]] == 1)
	    continue;
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	dbl_EGlpNumCopy (t_u, dbl_INFTY);
	dbl_EGlpNumCopy (at_u, t_u);
	dbl_EGlpNumCopy (t_l, dbl_NINFTY);
	dbl_EGlpNumCopy (at_l, t_u);
	if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY)) {
	    dbl_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
	    dbl_EGlpNumCopyAbs (at_l, t_l);
	}
	if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY)) {
	    dbl_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
	    dbl_EGlpNumCopyAbs (at_u, t_u);
	}
	/* t_i = (fabs (t_l) < fabs (t_u)) ? t_l : t_u; */
	if (dbl_EGlpNumIsLess (at_l, at_u)) {
	    dbl_EGlpNumCopy (t_i, t_l);
	    dbl_EGlpNumCopy (at_i, at_l);
	} else {
	    dbl_EGlpNumCopy (t_i, t_u);
	    dbl_EGlpNumCopy (at_i, at_u);
	}
	/* if (fabs (t_i) <= t_max + t_max * (1.0e-2)) */
	if (dbl_EGlpNumIsLeq (at_i, t_max)) {
	    if (dbl_EGlpNumIsLess (ayi_max, ay_ij)) {
		dbl_EGlpNumCopy (yi_max, y_ij);
		dbl_EGlpNumCopy (ayi_max, ay_ij);
		rs->lindex = i;
		dbl_EGlpNumCopy (rs->tz, t_i);
		rs->lvstat = (dbl_EGlpNumIsLess (at_l, at_u)) ? STAT_LOWER : STAT_UPPER;
	    }
	}
    }

    if (rs->lindex < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	rs->ratio_stat = RATIO_BCHANGE;
	dbl_EGlpNumCopy (rs->pivotval, yi_max);
    }
CLEANUP:
    for (i = 0; i < rcnt; i++)
	lp->iwork[rlist[i]] = 0;
    dbl_EGlpNumClearVar (t_max);
    dbl_EGlpNumClearVar (ay_ij);
    dbl_EGlpNumClearVar (at_i);
    dbl_EGlpNumClearVar (at_l);
    dbl_EGlpNumClearVar (at_u);
    dbl_EGlpNumClearVar (ayi_max);
    dbl_EGlpNumClearVar (y_ij);
    dbl_EGlpNumClearVar (t_i);
    dbl_EGlpNumClearVar (t_l);
    dbl_EGlpNumClearVar (t_u);
    dbl_EGlpNumClearVar (yi_max);
    return;
}
