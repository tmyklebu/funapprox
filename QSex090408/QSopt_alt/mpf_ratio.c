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

/* "$RCSfile: mpf_ratio.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#include "config.h"
#include "mpf_sortrus.h"
#include "stddefs.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdefs.h"
#include "mpf_ratio.h"
#include "mpf_fct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

void mpf_ILLratio_pI_test (mpf_lpinfo * lp,
      int eindex,
      int dir,
      mpf_ratio_res * rs)
{
    int i = 0, k = 0;
    int col, ecol;
    int cbnd, indx = 0;
    int tctr = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    mpf_t *pivtol = &(lp->tol->pivot_tol);
    mpf_t *dftol = &(lp->tol->id_tol);
     /* HHH */ mpf_t *t = lp->upd.t;
    mpf_t t_i, delta, y_ij, rcost, nrcost, ntmp;
    mpf_t *x, *l, *u;
     /* HHH */ mpf_EGlpNumInitVar (t_i);
    mpf_EGlpNumInitVar (delta);
    mpf_EGlpNumInitVar (y_ij);
    mpf_EGlpNumInitVar (rcost);
    mpf_EGlpNumInitVar (nrcost);
    mpf_EGlpNumInitVar (ntmp);
    mpf_EGlpNumZero (t_i);
    mpf_EGlpNumZero (y_ij);
    mpf_EGlpNumZero (delta);
    rs->lindex = -1;
    mpf_EGlpNumZero (rs->tz);
    mpf_EGlpNumZero (rs->pivotval);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    ecol = lp->nbaz[eindex];
    ILL_IFTRACE2 ("%s:%d:%d:%d:%d", __func__, eindex, dir, ecol,
	(VBOUNDED == lp->vtype[ecol]));
    if (lp->vtype[ecol] == VBOUNDED) {
	mpf_EGlpNumCopyDiff (t[0], lp->uz[ecol], lp->lz[ecol]);
	ix[0] = BBOUND;
	ILL_IFTRACE2 (":%d[%d](%la,%la,%la)\n", ix[tctr], tctr,
	    mpf_EGlpNumToLf (t[tctr]), mpf_EGlpNumToLf (lp->uz[ecol]),
	    mpf_EGlpNumToLf (lp->lz[ecol]));
	tctr++;
    }
    ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	mpf_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	if (mpf_EGlpNumIsEqual (y_ij, mpf_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && mpf_EGlpNumIsLess (mpf_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum))) {
	    if (mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum))
		mpf_EGlpNumSign (y_ij);
	    ILL_IFTRACE2 (":%d", lp->bfeas[i]);
	    if (lp->bfeas[i] > 0) {
		mpf_EGlpNumCopyDiffRatio (t[tctr], *x, *u, y_ij);
		ix[tctr] = 10 * k + BATOUPPER;
		ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, mpf_EGlpNumToLf (t[tctr]));
		tctr++;
		if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY)) {
		    mpf_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
		    ix[tctr] = 10 * k + BATOLOWER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpf_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    } else if (lp->bfeas[i] == 0) {
		if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY)) {
		    mpf_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
		    ix[tctr] = 10 * k + BATOLOWER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpf_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    }
	} else if ((dir == VINCREASE && mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLess (mpf_zeroLpNum, y_ij))) {
	    if (mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum))
		mpf_EGlpNumSign (y_ij);
	    ILL_IFTRACE2 (":%d", lp->bfeas[i]);
	    if (lp->bfeas[i] < 0) {
		mpf_EGlpNumCopyDiffRatio (t[tctr], *l, *x, y_ij);
		ix[tctr] = 10 * k + BBTOLOWER;
		ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, mpf_EGlpNumToLf (t[tctr]));
		tctr++;
		if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY)) {
		    mpf_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpf_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    } else if (lp->bfeas[i] == 0) {
		if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY)) {
		    mpf_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpf_EGlpNumToLf (t[tctr]));
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
    mpf_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

    mpf_EGlpNumZero (lp->upd.c_obj);
    mpf_EGlpNumCopy (rcost, lp->pIdz[eindex]);
    ILL_IFTRACE2 ("\n%s:%d:%lf", __func__, tctr, mpf_EGlpNumToLf (rcost));
    for (i = 0; i < tctr; i++) {
	mpf_EGlpNumCopy (t_i, t[perm[i]]);
	mpf_EGlpNumCopy (ntmp, t_i);
	mpf_EGlpNumSubTo (ntmp, delta);
	mpf_EGlpNumAddInnProdTo (lp->upd.c_obj, ntmp, rcost);
	mpf_EGlpNumCopy (delta, t_i);
	ILL_IFTRACE2 (":%d:%lf", perm[i], mpf_EGlpNumToLf (delta));
	 /* HHH */ cbnd = ix[perm[i]] % 10;
	if (cbnd != BBOUND) {
	    k = ix[perm[i]] / 10;
	    mpf_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	    indx = lp->yjz.indx[k];
	    ILL_IFTRACE2 (":%d", indx);
	}
	switch (cbnd) {
	case BBOUND:
	    rs->ratio_stat = RATIO_NOBCHANGE;
	    mpf_EGlpNumCopy (rs->tz, t_i);
	    if (dir != VINCREASE)
		mpf_EGlpNumSign (rs->tz);
	    ILL_CLEANUP;

	case BATOLOWER:
	case BATOUPPER:
	    mpf_EGlpNumAddTo (rcost, y_ij);
	    break;
	case BBTOLOWER:
	case BBTOUPPER:
	    mpf_EGlpNumSubTo (rcost, y_ij);
	    break;
	}
	mpf_EGlpNumCopyNeg (nrcost, rcost);
	if ((dir == VINCREASE && mpf_EGlpNumIsLeq (nrcost, *dftol)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLeq (rcost, *dftol))) {
	    /* change 5 to -1 if t_i > 0 is required below */
	    if (mpf_EGlpNumIsLess (t_i, mpf_zeroLpNum) && i > 5) {
		/* printf ("pIhell %.5f %d\n", t_i, i); */
		mpf_EGlpNumDivUiTo (t_i, 2);
		rs->ratio_stat = RATIO_NEGATIVE;
		mpf_EGlpNumZero (rs->tz);
		ILL_CLEANUP;
	    }
	    rs->lindex = indx;
	    rs->ratio_stat = RATIO_BCHANGE;
	    if (cbnd == BATOLOWER || cbnd == BBTOLOWER)
		rs->lvstat = STAT_LOWER;
	    else
		rs->lvstat = STAT_UPPER;

	    mpf_EGlpNumCopy (rs->pivotval, y_ij);
	    mpf_EGlpNumCopy (rs->tz, t_i);
	    if (dir != VINCREASE)
		mpf_EGlpNumSign (rs->tz);
	    ILL_CLEANUP;
	}
    }

CLEANUP:
    mpf_ILLfct_update_counts (lp, CNT_PIPIV, 0, rs->pivotval);
    ILL_IFTRACE2 (":tctr %d:%d\n", tctr, rs->ratio_stat);
    lp->upd.tctr = tctr;
    lp->upd.i = i;
    mpf_EGlpNumCopy (lp->upd.tz, t_i);
    mpf_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    if (dir == VDECREASE)
	mpf_EGlpNumSign (lp->upd.c_obj);
    if (rs->lindex != -1)
	lp->upd.fs = lp->bfeas[rs->lindex];
    mpf_EGlpNumClearVar (t_i);
    mpf_EGlpNumClearVar (delta);
    mpf_EGlpNumClearVar (y_ij);
    mpf_EGlpNumClearVar (rcost);
    mpf_EGlpNumClearVar (nrcost);
    mpf_EGlpNumClearVar (ntmp);
}

void mpf_ILLratio_pII_test (mpf_lpinfo * lp,
      int eindex,
      int dir,
      mpf_ratio_res * rs)
{
    int i, k, indx, col, ecol;
    mpf_t *x, *l, *u, t_max, ayi_max, yi_max, ay_ij, y_ij, t_i, t_z;
    mpf_t *pivtol = &(lp->tol->pivot_tol);
    mpf_t *pftol = &(lp->tol->pfeas_tol);
    mpf_EGlpNumInitVar (y_ij);
    mpf_EGlpNumInitVar (ay_ij);
    mpf_EGlpNumInitVar (t_i);
    mpf_EGlpNumInitVar (t_z);
    mpf_EGlpNumInitVar (t_max);
    mpf_EGlpNumInitVar (yi_max);
    mpf_EGlpNumInitVar (ayi_max);
     /* HHH */ rs->boundch = 0;
    rs->lindex = -1;
    mpf_EGlpNumZero (rs->tz);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    mpf_EGlpNumZero (rs->pivotval);
    mpf_EGlpNumZero (rs->lbound);
    ecol = lp->nbaz[eindex];

    for (k = 0, mpf_EGlpNumCopy (t_max, mpf_INFTY); k < lp->yjz.nzcnt; k++) {
	mpf_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	mpf_EGlpNumCopyAbs (ay_ij, y_ij);
	if (mpf_EGlpNumIsEqual (y_ij, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_i, mpf_INFTY);
	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && mpf_EGlpNumIsLess (mpf_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum))) {
	    if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY)) {
		mpf_EGlpNumCopyDiff (t_i, *x, *l);
		mpf_EGlpNumAddTo (t_i, *pftol);
		mpf_EGlpNumDivTo (t_i, ay_ij);
	    }
	} else if ((dir == VINCREASE && mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLess (mpf_zeroLpNum, y_ij))) {
	    if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY)) {
		mpf_EGlpNumCopySum (t_i, *u, *pftol);
		mpf_EGlpNumSubTo (t_i, *x);
		mpf_EGlpNumDivTo (t_i, ay_ij);
	    }
	}
	if (mpf_EGlpNumIsEqqual (t_i, mpf_INFTY))
	    continue;

	if (mpf_EGlpNumIsLess (t_i, t_max)) {
	    /* HHH tind = i; yval = fabs (y_ij); tval = t_i -
	       pftol/fabs(y_ij); */
	    mpf_EGlpNumCopy (t_max, t_i);
	}
    }
    /* we use yi_max as temporal variable here */
    mpf_EGlpNumCopyDiff (yi_max, lp->uz[ecol], lp->lz[ecol]);
    if (lp->vtype[ecol] == VBOUNDED && mpf_EGlpNumIsLeq (yi_max, t_max)) {

	mpf_EGlpNumCopy (t_max, yi_max);
	rs->ratio_stat = RATIO_NOBCHANGE;
	mpf_EGlpNumCopy (rs->tz, t_max);
	if (dir != VINCREASE)
	    mpf_EGlpNumSign (rs->tz);
	ILL_CLEANUP;
    }
    if (mpf_EGlpNumIsLeq (mpf_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    /* if (mpf_EGlpNumIsLess (t_max, mpf_zeroLpNum)) printf ("pIIhell\n"); */
    indx = -1;
    mpf_EGlpNumZero (t_z);
    mpf_EGlpNumZero (yi_max);
    mpf_EGlpNumZero (ayi_max);
    ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	mpf_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	mpf_EGlpNumCopyAbs (ay_ij, y_ij);
	if (mpf_EGlpNumIsEqual (y_ij, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_i, mpf_INFTY);
	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && mpf_EGlpNumIsLess (mpf_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum))) {
	    if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY))
		mpf_EGlpNumCopyDiffRatio (t_i, *x, *l, ay_ij);
	} else if ((dir == VINCREASE && mpf_EGlpNumIsLess (y_ij, mpf_zeroLpNum)) ||
	    (dir == VDECREASE && mpf_EGlpNumIsLess (mpf_zeroLpNum, y_ij))) {
	    if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY))
		mpf_EGlpNumCopyDiffRatio (t_i, *u, *x, ay_ij);
	}
	if (mpf_EGlpNumIsLeq (t_i, t_max)) {
	    if (mpf_EGlpNumIsLess (ayi_max, ay_ij)) {
		mpf_EGlpNumCopy (yi_max, y_ij);
		mpf_EGlpNumCopy (ayi_max, ay_ij);
		indx = i;
		mpf_EGlpNumCopy (t_z, t_i);
		ILL_IFTRACE2 (":%d:%lf:%lf:%lf:%lf", indx, mpf_EGlpNumToLf (t_i),
		    mpf_EGlpNumToLf (t_max), mpf_EGlpNumToLf (ayi_max),
		    mpf_EGlpNumToLf (ay_ij));
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
	mpf_EGlpNumCopy (rs->tz, t_z);
	mpf_EGlpNumCopy (rs->pivotval, yi_max);
	rs->ratio_stat = RATIO_BCHANGE;

	if (dir == VINCREASE)
	    rs->lvstat =
		(mpf_EGlpNumIsLess (mpf_zeroLpNum, yi_max)) ? STAT_LOWER : STAT_UPPER;
	else
	    rs->lvstat =
		(mpf_EGlpNumIsLess (mpf_zeroLpNum, yi_max)) ? STAT_UPPER : STAT_LOWER;

	if (mpf_EGlpNumIsLess (rs->tz, mpf_zeroLpNum)) {
	    ILL_IFTRACE2 ("need to change bound, tz=%la\n", mpf_EGlpNumToLf (rs->tz));
	    mpf_EGlpNumCopyAbs (rs->tz, t_max);
	    mpf_EGlpNumDivUiTo (rs->tz, 10);
	    rs->boundch = 1;
	    mpf_EGlpNumCopy (rs->lbound, lp->xbz[rs->lindex]);
	    if (rs->lvstat == STAT_LOWER)
		mpf_EGlpNumSubInnProdTo (rs->lbound, rs->tz, ayi_max);
	    else
		mpf_EGlpNumAddInnProdTo (rs->lbound, rs->tz, ayi_max);
	}
	if (dir == VDECREASE)
	    mpf_EGlpNumSign (rs->tz);
    }
CLEANUP:
    mpf_ILLfct_update_counts (lp, CNT_PIIPIV, 0, rs->pivotval);
    mpf_EGlpNumClearVar (y_ij);
    mpf_EGlpNumClearVar (ay_ij);
    mpf_EGlpNumClearVar (t_i);
    mpf_EGlpNumClearVar (t_z);
    mpf_EGlpNumClearVar (t_max);
    mpf_EGlpNumClearVar (yi_max);
    mpf_EGlpNumClearVar (ayi_max);
}

#define mpf_GET_XY_DRATIOTEST \
      if (lp->vstat[col] == STAT_UPPER){ \
                mpf_EGlpNumCopyNeg(x,lp->dz[j]);\
        mpf_EGlpNumCopy(y, *zAj);\
      } \
      else{ \
         mpf_EGlpNumCopy(x, lp->dz[j]); \
         mpf_EGlpNumCopyNeg(y, *zAj);\
      } \
      if (lvstat == STAT_UPPER) \
         mpf_EGlpNumSign(y);


void mpf_ILLratio_dI_test (mpf_lpinfo * lp,
      int lindex,
      int lvstat,
      mpf_ratio_res * rs)
{
    int j = 0, k;
    int col;
    int cbnd, indx;
    int tctr = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    mpf_t *t = lp->upd.t;
    mpf_t *zAj, x, y, t_j, theta, rcost, delta;
    mpf_t *pftol = &(lp->tol->ip_tol);
    mpf_t *pivtol = &(lp->tol->pivot_tol);
    mpf_EGlpNumInitVar (x);
    mpf_EGlpNumInitVar (y);
    mpf_EGlpNumInitVar (t_j);
    mpf_EGlpNumInitVar (theta);
    mpf_EGlpNumInitVar (rcost);
    mpf_EGlpNumInitVar (delta);
    mpf_EGlpNumZero (delta);
    mpf_EGlpNumZero (t_j);
    mpf_EGlpNumZero (rs->tz);
     /* HHH */ rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;
    mpf_EGlpNumZero (rs->pivotval);

    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (mpf_EGlpNumIsEqual (*zAj, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_j, mpf_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	mpf_GET_XY_DRATIOTEST;

	if (mpf_EGlpNumIsLess (y, mpf_zeroLpNum)) {
	    if (lp->dfeas[j] != 0 && lp->vstat[col] != STAT_ZERO) {
		mpf_EGlpNumCopyFrac (t[tctr], x, y);
		ix[tctr] = 10 * k + BBTOLOWER;
		tctr++;
	    } else if (lp->vstat[col] == STAT_ZERO) {
		if (lp->dfeas[j] < 0) {
		    mpf_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BBTOLOWER;
		    tctr++;
		}
		if (lp->dfeas[j] <= 0) {
		    mpf_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    tctr++;
		}
	    }
	} else {
	    if (lp->dfeas[j] > 0) {
		if (lp->vstat[col] == STAT_ZERO) {
		    mpf_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BATOUPPER;
		    tctr++;
		    mpf_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BATOLOWER;
		    tctr++;
		}
	    } else if (lp->dfeas[j] == 0) {
		mpf_EGlpNumCopyFrac (t[tctr], x, y);
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
    mpf_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

    mpf_EGlpNumZero (lp->upd.c_obj);
    mpf_EGlpNumCopy (rcost, lp->xbz[lindex]);
    if (lvstat == STAT_LOWER)
	mpf_EGlpNumSign (rcost);
    for (j = 0; j < tctr; j++) {
	cbnd = ix[perm[j]] % 10;
	if (cbnd == BSKIP)
	    continue;

	mpf_EGlpNumCopy (t_j, t[perm[j]]);
	mpf_EGlpNumCopy (x, t_j);
	mpf_EGlpNumSubTo (x, delta);
	mpf_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	mpf_EGlpNumCopy (delta, t_j);
	k = ix[perm[j]] / 10;
	zAj = &(lp->zA.coef[k]);
	indx = lp->zA.indx[k];

	if (lp->vstat[lp->nbaz[indx]] == STAT_LOWER
	    || lp->vstat[lp->nbaz[indx]] == STAT_ZERO)
	    mpf_EGlpNumCopyNeg (theta, *zAj);
	else
	    mpf_EGlpNumCopy (theta, *zAj);

	if (lvstat == STAT_UPPER)
	    mpf_EGlpNumSign (theta);

	switch (cbnd) {
	case BATOLOWER:
	case BATOUPPER:
	    mpf_EGlpNumSubTo (rcost, theta);
	    break;
	case BBTOLOWER:
	case BBTOUPPER:
	    mpf_EGlpNumAddTo (rcost, theta);
	    break;
	}
	if (mpf_EGlpNumIsLeq (rcost, *pftol)) {
	    /* if (t_j < 0.0) printf ("dIhell\n"); */
	    rs->eindex = indx;
	    mpf_EGlpNumCopy (rs->tz, t_j);
	    mpf_EGlpNumCopy (rs->pivotval, *zAj);
	    rs->ratio_stat = RATIO_BCHANGE;
	    ILL_CLEANUP;
	}
    }

CLEANUP:
    mpf_ILLfct_update_counts (lp, CNT_DIPIV, 0, rs->pivotval);
    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
    lp->upd.tctr = tctr;
    lp->upd.i = j;
    mpf_EGlpNumCopyAbs (lp->upd.tz, t_j);
    mpf_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    if (rs->eindex != -1)
	lp->upd.fs = lp->dfeas[rs->eindex];
    mpf_EGlpNumClearVar (x);
    mpf_EGlpNumClearVar (y);
    mpf_EGlpNumClearVar (t_j);
    mpf_EGlpNumClearVar (theta);
    mpf_EGlpNumClearVar (rcost);
    mpf_EGlpNumClearVar (delta);
}

void mpf_ILLratio_dII_test (mpf_lpinfo * lp, int lvstat, mpf_ratio_res * rs)
{
    int j, k, indx;
    int col, ecol;
    mpf_t *zAj, azAj, az_max, x, y, t_j, z_max, t_max, t_z;
    mpf_t *dftol = &(lp->tol->dfeas_tol);
    mpf_t *pivtol = &(lp->tol->pivot_tol);
    mpf_EGlpNumInitVar (x);
    mpf_EGlpNumInitVar (y);
    mpf_EGlpNumInitVar (t_j);
    mpf_EGlpNumInitVar (z_max);
    mpf_EGlpNumInitVar (t_max);
    mpf_EGlpNumInitVar (az_max);
    mpf_EGlpNumInitVar (azAj);
    mpf_EGlpNumInitVar (t_z);
    mpf_EGlpNumZero (t_j);
    rs->coeffch = 0;
    mpf_EGlpNumZero (rs->ecoeff);
    rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;
    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
    lp->upd.tctr = 0;
    mpf_EGlpNumZero (lp->upd.dty);
    for (k = 0, mpf_EGlpNumCopy (t_max, mpf_INFTY); k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (mpf_EGlpNumIsEqual (*zAj, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_j, mpf_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	mpf_GET_XY_DRATIOTEST;

	/* warning adding/substracting tolerances to used value, is it right? */
	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, y)) {
	    /* t_j = (x + dftol) / y; */
	    mpf_EGlpNumCopySum (t_j, x, *dftol);
	    mpf_EGlpNumDivTo (t_j, y);
	} else {
	    /* warning adding/substracting tolerances to used value, is it
	       right? */
	    if (lp->vstat[col] == STAT_ZERO)
		mpf_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
	}
	/* if (t_j == mpf_INFTY) */
	if (mpf_EGlpNumIsEqqual (t_j, mpf_INFTY))
	    continue;

	if (mpf_EGlpNumIsLess (t_j, t_max))
	    mpf_EGlpNumCopy (t_max, t_j);
    }

    if (mpf_EGlpNumIsLeq (mpf_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    /* if (t_max < 0.0) printf ("dIIhell\n"); */

    indx = -1;
    mpf_EGlpNumZero (t_z);
    mpf_EGlpNumZero (z_max);
    mpf_EGlpNumZero (az_max);

    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	mpf_EGlpNumCopyAbs (azAj, *zAj);
	if (mpf_EGlpNumIsEqual (*zAj, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_j, mpf_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	mpf_GET_XY_DRATIOTEST;

	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, y) || lp->vstat[col] == STAT_ZERO)
	    mpf_EGlpNumCopyFrac (t_j, x, y);

	if (mpf_EGlpNumIsLeq (t_j, t_max) && (mpf_EGlpNumIsLess (az_max, azAj))) {
	    mpf_EGlpNumCopy (z_max, *zAj);
	    mpf_EGlpNumCopy (az_max, azAj);
	    indx = j;
	    mpf_EGlpNumCopy (t_z, t_j);
	}
    }


    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	rs->eindex = indx;
	mpf_EGlpNumCopy (rs->tz, t_z);
	mpf_EGlpNumCopy (rs->pivotval, z_max);
	rs->ratio_stat = RATIO_BCHANGE;

	if (mpf_EGlpNumIsLess (rs->tz, mpf_zeroLpNum)) {
	    mpf_EGlpNumCopyAbs (rs->tz, t_max);
	    mpf_EGlpNumDivUiTo (rs->tz, 20);
	    rs->coeffch = 1;
	    ecol = lp->nbaz[indx];
	    mpf_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
	    switch (lp->vstat[ecol]) {
	    case STAT_LOWER:
		mpf_EGlpNumAddInnProdTo (rs->ecoeff, rs->tz, az_max);
		break;
	    case STAT_UPPER:
		mpf_EGlpNumSubInnProdTo (rs->ecoeff, rs->tz, az_max);
		break;
	    default:
		mpf_EGlpNumZero (rs->tz);
		break;
	    }
	}
    }

CLEANUP:
    mpf_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
    mpf_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    mpf_EGlpNumClearVar (x);
    mpf_EGlpNumClearVar (y);
    mpf_EGlpNumClearVar (t_j);
    mpf_EGlpNumClearVar (z_max);
    mpf_EGlpNumClearVar (t_max);
    mpf_EGlpNumClearVar (t_z);
    mpf_EGlpNumClearVar (az_max);
    mpf_EGlpNumClearVar (azAj);
}

void mpf_ILLratio_longdII_test (mpf_lpinfo * lp,
      int lindex,
      int lvstat,
      mpf_ratio_res * rs)
{
    int j, k, indx = 0, tctr = 0;
    int col, ecol;
    int vs, bnd_exist = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int b_indx = -1;
    mpf_t *t = lp->upd.t;
    mpf_t *l, *u, *xb, *zAj = 0, x, y, t_j, z_max, t_max, t_z, theta, rcost,
      delta, zb_val, tb_val, az_max, azb_val, azAj;
    mpf_t *pftol = &(lp->tol->pfeas_tol);
    mpf_t *dftol = &(lp->tol->dfeas_tol);
    mpf_t *pivtol = &(lp->tol->pivot_tol);
    mpf_EGlpNumInitVar (x);
    mpf_EGlpNumInitVar (azAj);
    mpf_EGlpNumInitVar (y);
    mpf_EGlpNumInitVar (t_j);
    mpf_EGlpNumInitVar (z_max);
    mpf_EGlpNumInitVar (az_max);
    mpf_EGlpNumInitVar (t_max);
    mpf_EGlpNumInitVar (t_z);
    mpf_EGlpNumInitVar (theta);
    mpf_EGlpNumInitVar (rcost);
    mpf_EGlpNumInitVar (delta);
    mpf_EGlpNumInitVar (zb_val);
    mpf_EGlpNumInitVar (azb_val);
    mpf_EGlpNumInitVar (tb_val);
    mpf_EGlpNumZero (t_j);
    mpf_EGlpNumZero (delta);
    mpf_EGlpNumZero (zb_val);
    mpf_EGlpNumZero (azb_val);
    mpf_EGlpNumCopy (tb_val, mpf_NINFTY);
    /* #warning not sure about THIS line */
    mpf_EGlpNumZero (rs->pivotval);

    rs->coeffch = 0;
    rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;

    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
    lp->upd.tctr = 0;
    lp->upd.i = 0;
    mpf_EGlpNumZero (lp->upd.tz);
    mpf_EGlpNumZero (lp->upd.piv);
    mpf_EGlpNumZero (lp->upd.c_obj);
    mpf_EGlpNumZero (lp->upd.dty);

    xb = &(lp->xbz[lindex]);
    col = lp->baz[lindex];
    l = &(lp->lz[col]);
    u = &(lp->uz[col]);
    /* rcost = (lvstat == STAT_LOWER) ? l - xb : xb - u; */
    if (lvstat == STAT_LOWER)
	mpf_EGlpNumCopyDiff (rcost, *l, *xb);
    else
	mpf_EGlpNumCopyDiff (rcost, *xb, *u);

    for (k = 0, mpf_EGlpNumCopy (t_max, mpf_INFTY); k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (mpf_EGlpNumIsEqual (*zAj, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_j, mpf_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;
	if (lp->vtype[col] == VBOUNDED) {
	    bnd_exist++;
	    continue;
	}
	mpf_GET_XY_DRATIOTEST;

	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, y)) {
	    /* t_j = (x + dftol) / y; */
	    /* #warning Using tolerances to add to result, is it right? */
	    mpf_EGlpNumCopySum (t_j, x, *dftol);
	    mpf_EGlpNumDivTo (t_j, y);
	} else {
	    if (lp->vstat[col] == STAT_ZERO)
		mpf_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
	}
	if (mpf_EGlpNumIsEqqual (t_j, mpf_INFTY))
	    continue;

	if (mpf_EGlpNumIsLess (t_j, t_max))
	    mpf_EGlpNumCopy (t_max, t_j);
    }
    if (mpf_EGlpNumIsLess (t_max, mpf_zeroLpNum)) {
	/* printf ("dIIhell, %.4f\n", t_max); */
	rs->ratio_stat = RATIO_NEGATIVE;
	ILL_CLEANUP;
    }
    if (bnd_exist == 0 && mpf_EGlpNumIsLeq (mpf_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	/*
         * printf ("x = %.8f, b = %.2f \n", lp->xbz[lindex], (lvstat == STAT_LOWER ) ? lp->lz[lp->baz[lindex]] : lp->uz[lp->baz[lindex]]);
         */
	ILL_CLEANUP;
    }
    if (bnd_exist != 0) {
	for (k = 0; k < lp->zA.nzcnt; k++) {
	    zAj = &(lp->zA.coef[k]);
	    if (mpf_EGlpNumIsEqual (*zAj, mpf_zeroLpNum, *pivtol))
		continue;

	    mpf_EGlpNumCopy (t_j, mpf_INFTY);
	    j = lp->zA.indx[k];
	    col = lp->nbaz[j];

	    if (lp->vtype[col] != VBOUNDED)
		continue;

	    mpf_GET_XY_DRATIOTEST;

	    if (mpf_EGlpNumIsLess (mpf_zeroLpNum, y)) {
		mpf_EGlpNumCopyFrac (t_j, x, y);
		if (mpf_EGlpNumIsLeq (t_j, t_max)) {
		    mpf_EGlpNumCopy (t[tctr], t_j);
		    ix[tctr] = k;
		    tctr++;
		}
	    }
	}
    }
    if (tctr != 0) {
	for (j = 0; j < tctr; j++)
	    perm[j] = j;
	mpf_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	for (j = 0; j < tctr; j++) {

	    mpf_EGlpNumCopy (t_j, t[perm[j]]);
	    /* we use x as temporal storage */
	    /* lp->upd.c_obj += (t_j - delta) * rcost; */
	    mpf_EGlpNumCopy (x, t_j);
	    mpf_EGlpNumSubTo (x, delta);
	    mpf_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	    mpf_EGlpNumCopy (delta, t_j);
	     /* HHH */ k = ix[perm[j]];
	    zAj = &(lp->zA.coef[k]);
	    indx = lp->zA.indx[k];
	    col = lp->nbaz[indx];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);
	    vs = lp->vstat[col];
	    /* theta = (vs == STAT_UPPER) ? (l - u) * zAj : (u - l) * zAj; */
	    mpf_EGlpNumCopyDiff (theta, *l, *u);
	    mpf_EGlpNumMultTo (theta, *zAj);
	    if (vs != STAT_UPPER)
		mpf_EGlpNumSign (theta);
	    if (lvstat == STAT_LOWER)
		mpf_EGlpNumAddTo (rcost, theta);
	    else
		mpf_EGlpNumSubTo (rcost, theta);

	    if (mpf_EGlpNumIsLeq (rcost, *pftol)) {
		rs->eindex = indx;
		mpf_EGlpNumCopy (rs->tz, t_j);
		mpf_EGlpNumCopy (rs->pivotval, *zAj);
		rs->ratio_stat = RATIO_BCHANGE;

		if (mpf_EGlpNumIsLess (rs->tz, mpf_zeroLpNum)) {
		    mpf_EGlpNumZero (rs->tz);
		    rs->coeffch = 1;
		    /* rs->ecoeff = lp->cz[col] - lp->dz[indx]; */
		    mpf_EGlpNumCopyDiff (rs->ecoeff, lp->cz[col], lp->dz[indx]);
		    /* lp->upd.c_obj += (rs->tz - delta) * rcost; note ts->tz
		       == 0; */
		    mpf_EGlpNumSubInnProdTo (lp->upd.c_obj, delta, rcost);
		}
		ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
		lp->upd.tctr = tctr;
		lp->upd.i = j;
		mpf_EGlpNumCopy (lp->upd.tz, rs->tz);
		ILL_CLEANUP;
	    }
	}
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
	lp->upd.tctr = tctr;
	lp->upd.i = tctr;
	mpf_EGlpNumCopy (lp->upd.tz, t_j);
	mpf_EGlpNumCopy (zb_val, *zAj);
	mpf_EGlpNumCopyAbs (azb_val, zb_val);
	mpf_EGlpNumCopy (tb_val, t_j);
	b_indx = indx;
    }
    if (bnd_exist != 0 && mpf_EGlpNumIsLeq (mpf_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	/* printf ("rcost: %.8f\n", rcost); */
	ILL_CLEANUP;
    }
    mpf_EGlpNumZero (z_max);
    mpf_EGlpNumZero (az_max);
    indx = -1;
    mpf_EGlpNumZero (t_z);
    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	mpf_EGlpNumCopyAbs (azAj, *zAj);
	if (mpf_EGlpNumIsEqual (*zAj, mpf_zeroLpNum, *pivtol))
	    continue;

	mpf_EGlpNumCopy (t_j, mpf_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED ||
	    lp->vtype[col] == VBOUNDED)
	    continue;

	mpf_GET_XY_DRATIOTEST;

	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, y) || lp->vstat[col] == STAT_ZERO)
	    mpf_EGlpNumCopyFrac (t_j, x, y);

	if (mpf_EGlpNumIsLeq (t_j, t_max)) {
	    if (mpf_EGlpNumIsLess (az_max, azAj)) {
		mpf_EGlpNumCopy (z_max, *zAj);
		mpf_EGlpNumCopy (az_max, azAj);
		indx = j;
		mpf_EGlpNumCopy (t_z, t_j);
	    }
	}
    }

    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
	ILL_CLEANUP;
    }
    if ((tctr == 0) || (mpf_EGlpNumIsLess (tb_val, mpf_zeroLpNum)) ||
	(tctr != 0 && mpf_EGlpNumIsLeq (tb_val, t_z) &&
	    mpf_EGlpNumIsLeq (azb_val, az_max))) {
	/* we use x as temporal vvariable */
	/* lp->upd.c_obj += (t_z - delta) * rcost; */
	mpf_EGlpNumCopyDiff (x, t_z, delta);
	mpf_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	mpf_EGlpNumCopy (delta, t_z);
	rs->eindex = indx;
	mpf_EGlpNumCopy (rs->tz, t_z);
	mpf_EGlpNumCopy (rs->pivotval, z_max);
	rs->ratio_stat = RATIO_BCHANGE;
    }
    /* For now */
    else if (tctr != 0) {
	rs->eindex = b_indx;
	mpf_EGlpNumCopy (rs->tz, tb_val);
	mpf_EGlpNumCopy (rs->pivotval, zb_val);
	rs->ratio_stat = RATIO_BCHANGE;
	lp->upd.i -= 1;
    }
    if (mpf_EGlpNumIsLess (rs->tz, mpf_zeroLpNum)) {
	/* if (tctr != 0) printf ("despite long step\n"); */
	/* rs->tz = fabs (t_max / 20.0); */
	mpf_EGlpNumCopyAbs (rs->tz, t_max);
	mpf_EGlpNumDivUiTo (rs->tz, 20);
	rs->coeffch = 1;

	ecol = lp->nbaz[indx];
	if (lp->vstat[ecol] == STAT_LOWER) {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs
	       (z_max); */
	    mpf_EGlpNumCopy (rs->ecoeff, az_max);
	    mpf_EGlpNumMultTo (rs->ecoeff, rs->tz);
	    mpf_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
	    mpf_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
	} else if (lp->vstat[ecol] == STAT_UPPER) {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs
	       (z_max); */
	    mpf_EGlpNumCopy (rs->ecoeff, az_max);
	    mpf_EGlpNumMultTo (rs->ecoeff, rs->tz);
	    mpf_EGlpNumSign (rs->ecoeff);
	    mpf_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
	    mpf_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
	} else {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx]; */
	    mpf_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
	    mpf_EGlpNumZero (rs->tz);
	}
	/* we use x as temporal storage */
	/* lp->upd.c_obj += (rs->tz - delta) * rcost; */
	mpf_EGlpNumCopy (x, rs->tz);
	mpf_EGlpNumSubTo (x, delta);
	mpf_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
    }
CLEANUP:
    mpf_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
    mpf_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    mpf_EGlpNumClearVar (x);
    mpf_EGlpNumClearVar (y);
    mpf_EGlpNumClearVar (t_j);
    mpf_EGlpNumClearVar (z_max);
    mpf_EGlpNumClearVar (az_max);
    mpf_EGlpNumClearVar (t_max);
    mpf_EGlpNumClearVar (t_z);
    mpf_EGlpNumClearVar (theta);
    mpf_EGlpNumClearVar (rcost);
    mpf_EGlpNumClearVar (delta);
    mpf_EGlpNumClearVar (zb_val);
    mpf_EGlpNumClearVar (azb_val);
    mpf_EGlpNumClearVar (tb_val);
    mpf_EGlpNumClearVar (azAj);
}

void mpf_ILLratio_pivotin_test (mpf_lpinfo * lp,
      int *rlist,
      int rcnt,
      mpf_ratio_res * rs)
{
    int i, k, col;
    mpf_t *x, *l, *u;
    mpf_t ay_ij, at_i, at_l, at_u, ayi_max, y_ij, t_i, t_l, t_u, t_max,
      yi_max;
    mpf_t *pivtol = &(lp->tol->pivot_tol);
    if (rcnt <= 0 || rs == NULL)
	return;
    mpf_EGlpNumInitVar (ay_ij);
    mpf_EGlpNumInitVar (at_i);
    mpf_EGlpNumInitVar (at_l);
    mpf_EGlpNumInitVar (at_u);
    mpf_EGlpNumInitVar (ayi_max);
    mpf_EGlpNumInitVar (t_max);
    mpf_EGlpNumInitVar (y_ij);
    mpf_EGlpNumInitVar (t_i);
    mpf_EGlpNumInitVar (t_l);
    mpf_EGlpNumInitVar (t_u);
    mpf_EGlpNumInitVar (yi_max);
    rs->boundch = 0;
    rs->lindex = -1;
    mpf_EGlpNumZero (rs->tz);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    mpf_EGlpNumZero (rs->pivotval);
    mpf_EGlpNumZero (rs->lbound);

    for (i = 0; i < rcnt; i++)
	lp->iwork[rlist[i]] = 1;

    for (k = 0, mpf_EGlpNumCopy (t_max, mpf_INFTY); k < lp->yjz.nzcnt; k++) {
	mpf_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	if (mpf_EGlpNumIsEqual (y_ij, mpf_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	if (lp->iwork[lp->baz[i]] == 1)
	    continue;
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	mpf_EGlpNumCopy (t_u, mpf_INFTY);
	mpf_EGlpNumCopy (at_u, mpf_INFTY);
	mpf_EGlpNumCopy (t_l, mpf_NINFTY);
	mpf_EGlpNumCopy (at_l, mpf_INFTY);

	if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY)) {
	    mpf_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
	    mpf_EGlpNumCopyAbs (at_l, t_l);
	    if (mpf_EGlpNumIsLess (at_l, t_max))
		mpf_EGlpNumCopy (t_max, at_l);
	}
	if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY)) {
	    mpf_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
	    mpf_EGlpNumCopyAbs (at_u, t_u);
	    if (mpf_EGlpNumIsLess (at_u, t_max))
		mpf_EGlpNumCopy (t_max, at_u);
	}
    }

    if (mpf_EGlpNumIsLeq (mpf_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    mpf_EGlpNumZero (yi_max);
    mpf_EGlpNumZero (ayi_max);
    mpf_EGlpNumMultUiTo (t_max, 101);
    mpf_EGlpNumDivUiTo (t_max, 100);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	mpf_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	mpf_EGlpNumCopyAbs (ay_ij, y_ij);
	if (mpf_EGlpNumIsEqual (y_ij, mpf_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	if (lp->iwork[lp->baz[i]] == 1)
	    continue;
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	mpf_EGlpNumCopy (t_u, mpf_INFTY);
	mpf_EGlpNumCopy (at_u, t_u);
	mpf_EGlpNumCopy (t_l, mpf_NINFTY);
	mpf_EGlpNumCopy (at_l, t_u);
	if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY)) {
	    mpf_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
	    mpf_EGlpNumCopyAbs (at_l, t_l);
	}
	if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY)) {
	    mpf_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
	    mpf_EGlpNumCopyAbs (at_u, t_u);
	}
	/* t_i = (fabs (t_l) < fabs (t_u)) ? t_l : t_u; */
	if (mpf_EGlpNumIsLess (at_l, at_u)) {
	    mpf_EGlpNumCopy (t_i, t_l);
	    mpf_EGlpNumCopy (at_i, at_l);
	} else {
	    mpf_EGlpNumCopy (t_i, t_u);
	    mpf_EGlpNumCopy (at_i, at_u);
	}
	/* if (fabs (t_i) <= t_max + t_max * (1.0e-2)) */
	if (mpf_EGlpNumIsLeq (at_i, t_max)) {
	    if (mpf_EGlpNumIsLess (ayi_max, ay_ij)) {
		mpf_EGlpNumCopy (yi_max, y_ij);
		mpf_EGlpNumCopy (ayi_max, ay_ij);
		rs->lindex = i;
		mpf_EGlpNumCopy (rs->tz, t_i);
		rs->lvstat = (mpf_EGlpNumIsLess (at_l, at_u)) ? STAT_LOWER : STAT_UPPER;
	    }
	}
    }

    if (rs->lindex < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	rs->ratio_stat = RATIO_BCHANGE;
	mpf_EGlpNumCopy (rs->pivotval, yi_max);
    }
CLEANUP:
    for (i = 0; i < rcnt; i++)
	lp->iwork[rlist[i]] = 0;
    mpf_EGlpNumClearVar (t_max);
    mpf_EGlpNumClearVar (ay_ij);
    mpf_EGlpNumClearVar (at_i);
    mpf_EGlpNumClearVar (at_l);
    mpf_EGlpNumClearVar (at_u);
    mpf_EGlpNumClearVar (ayi_max);
    mpf_EGlpNumClearVar (y_ij);
    mpf_EGlpNumClearVar (t_i);
    mpf_EGlpNumClearVar (t_l);
    mpf_EGlpNumClearVar (t_u);
    mpf_EGlpNumClearVar (yi_max);
    return;
}
