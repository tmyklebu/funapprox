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

/* "$RCSfile: mpq_ratio.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#include "config.h"
#include "mpq_sortrus.h"
#include "stddefs.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdefs.h"
#include "mpq_ratio.h"
#include "mpq_fct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

void mpq_ILLratio_pI_test (mpq_lpinfo * lp, int eindex, int dir,
      mpq_ratio_res * rs)
{
    int i = 0, k = 0;
    int col, ecol;
    int cbnd, indx = 0;
    int tctr = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    mpq_t *dftol = &(lp->tol->id_tol);
    mpq_t *t = lp->upd.t;
    mpq_t t_i, delta, y_ij, rcost, nrcost, ntmp;
    mpq_t *x, *l, *u;

    mpq_EGlpNumInitVar (t_i);
    mpq_EGlpNumInitVar (delta);
    mpq_EGlpNumInitVar (y_ij);
    mpq_EGlpNumInitVar (rcost);
    mpq_EGlpNumInitVar (nrcost);
    mpq_EGlpNumInitVar (ntmp);
    mpq_EGlpNumZero (t_i);
    mpq_EGlpNumZero (y_ij);
    mpq_EGlpNumZero (delta);
    rs->lindex = -1;
    mpq_EGlpNumZero (rs->tz);
    mpq_EGlpNumZero (rs->pivotval);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    ecol = lp->nbaz[eindex];
    ILL_IFTRACE2 ("%s:%d:%d:%d:%d", __func__, eindex, dir, ecol,
	(VBOUNDED == lp->vtype[ecol]));
    if (lp->vtype[ecol] == VBOUNDED) {
	mpq_EGlpNumCopyDiff (t[0], lp->uz[ecol], lp->lz[ecol]);
	ix[0] = BBOUND;
	ILL_IFTRACE2 (":%d[%d](%la,%la,%la)\n", ix[tctr], tctr,
	    mpq_EGlpNumToLf (t[tctr]), mpq_EGlpNumToLf (lp->uz[ecol]),
	    mpq_EGlpNumToLf (lp->lz[ecol]));
	tctr++;
    }
    ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	mpq_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	if (mpq_EGlpNumIsEqual (y_ij, mpq_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && mpq_EGlpNumIsLess (mpq_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum))) {
	    if (mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum))
		mpq_EGlpNumSign (y_ij);
	    ILL_IFTRACE2 (":%d", lp->bfeas[i]);
	    if (lp->bfeas[i] > 0) {
		mpq_EGlpNumCopyDiffRatio (t[tctr], *x, *u, y_ij);
		ix[tctr] = 10 * k + BATOUPPER;
		ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, mpq_EGlpNumToLf (t[tctr]));
		tctr++;
		if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY)) {
		    mpq_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
		    ix[tctr] = 10 * k + BATOLOWER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpq_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    } else if (lp->bfeas[i] == 0) {
		if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY)) {
		    mpq_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
		    ix[tctr] = 10 * k + BATOLOWER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpq_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    }
	} else if ((dir == VINCREASE && mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLess (mpq_zeroLpNum, y_ij))) {
	    if (mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum))
		mpq_EGlpNumSign (y_ij);
	    ILL_IFTRACE2 (":%d", lp->bfeas[i]);
	    if (lp->bfeas[i] < 0) {
		mpq_EGlpNumCopyDiffRatio (t[tctr], *l, *x, y_ij);
		ix[tctr] = 10 * k + BBTOLOWER;
		ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, mpq_EGlpNumToLf (t[tctr]));
		tctr++;
		if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY)) {
		    mpq_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpq_EGlpNumToLf (t[tctr]));
		    tctr++;
		}
	    } else if (lp->bfeas[i] == 0) {
		if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY)) {
		    mpq_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
			mpq_EGlpNumToLf (t[tctr]));
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
    mpq_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

    mpq_EGlpNumZero (lp->upd.c_obj);
    mpq_EGlpNumCopy (rcost, lp->pIdz[eindex]);
    ILL_IFTRACE2 ("\n%s:%d:%lf", __func__, tctr, mpq_EGlpNumToLf (rcost));
    for (i = 0; i < tctr; i++) {
	mpq_EGlpNumCopy (t_i, t[perm[i]]);
	mpq_EGlpNumCopy (ntmp, t_i);
	mpq_EGlpNumSubTo (ntmp, delta);
	mpq_EGlpNumAddInnProdTo (lp->upd.c_obj, ntmp, rcost);
	mpq_EGlpNumCopy (delta, t_i);
	ILL_IFTRACE2 (":%d:%lf", perm[i], mpq_EGlpNumToLf (delta));
	 /* HHH */ cbnd = ix[perm[i]] % 10;
	if (cbnd != BBOUND) {
	    k = ix[perm[i]] / 10;
	    mpq_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	    indx = lp->yjz.indx[k];
	    ILL_IFTRACE2 (":%d", indx);
	}
	switch (cbnd) {
	case BBOUND:
	    rs->ratio_stat = RATIO_NOBCHANGE;
	    mpq_EGlpNumCopy (rs->tz, t_i);
	    if (dir != VINCREASE)
		mpq_EGlpNumSign (rs->tz);
	    ILL_CLEANUP;

	case BATOLOWER:
	case BATOUPPER:
	    mpq_EGlpNumAddTo (rcost, y_ij);
	    break;
	case BBTOLOWER:
	case BBTOUPPER:
	    mpq_EGlpNumSubTo (rcost, y_ij);
	    break;
	}
	mpq_EGlpNumCopyNeg (nrcost, rcost);
	if ((dir == VINCREASE && mpq_EGlpNumIsLeq (nrcost, *dftol)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLeq (rcost, *dftol))) {
	    /* change 5 to -1 if t_i > 0 is required below */
	    if (mpq_EGlpNumIsLess (t_i, mpq_zeroLpNum) && i > 5) {
		/* printf ("pIhell %.5f %d\n", t_i, i); */
		mpq_EGlpNumDivUiTo (t_i, 2);
		rs->ratio_stat = RATIO_NEGATIVE;
		mpq_EGlpNumZero (rs->tz);
		ILL_CLEANUP;
	    }
	    rs->lindex = indx;
	    rs->ratio_stat = RATIO_BCHANGE;
	    if (cbnd == BATOLOWER || cbnd == BBTOLOWER)
		rs->lvstat = STAT_LOWER;
	    else
		rs->lvstat = STAT_UPPER;

	    mpq_EGlpNumCopy (rs->pivotval, y_ij);
	    mpq_EGlpNumCopy (rs->tz, t_i);
	    if (dir != VINCREASE)
		mpq_EGlpNumSign (rs->tz);
	    ILL_CLEANUP;
	}
    }

CLEANUP:
    mpq_ILLfct_update_counts (lp, CNT_PIPIV, 0, rs->pivotval);
    ILL_IFTRACE2 (":tctr %d:%d\n", tctr, rs->ratio_stat);
    lp->upd.tctr = tctr;
    lp->upd.i = i;
    mpq_EGlpNumCopy (lp->upd.tz, t_i);
    mpq_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    if (dir == VDECREASE)
	mpq_EGlpNumSign (lp->upd.c_obj);
    if (rs->lindex != -1)
	lp->upd.fs = lp->bfeas[rs->lindex];
    mpq_EGlpNumClearVar (t_i);
    mpq_EGlpNumClearVar (delta);
    mpq_EGlpNumClearVar (y_ij);
    mpq_EGlpNumClearVar (rcost);
    mpq_EGlpNumClearVar (nrcost);
    mpq_EGlpNumClearVar (ntmp);
}

void mpq_ILLratio_pII_test (mpq_lpinfo * lp, int eindex, int dir,
      mpq_ratio_res * rs)
{
    int i, k, indx, col, ecol;
    mpq_t *x, *l, *u, t_max, ayi_max, yi_max, ay_ij, y_ij, t_i, t_z;
    mpq_t *pftol = &(lp->tol->pfeas_tol);
    mpq_EGlpNumInitVar (y_ij);
    mpq_EGlpNumInitVar (ay_ij);
    mpq_EGlpNumInitVar (t_i);
    mpq_EGlpNumInitVar (t_z);
    mpq_EGlpNumInitVar (t_max);
    mpq_EGlpNumInitVar (yi_max);
    mpq_EGlpNumInitVar (ayi_max);
     /* HHH */ rs->boundch = 0;
    rs->lindex = -1;
    mpq_EGlpNumZero (rs->tz);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    mpq_EGlpNumZero (rs->pivotval);
    mpq_EGlpNumZero (rs->lbound);
    ecol = lp->nbaz[eindex];

    for (k = 0, mpq_EGlpNumCopy (t_max, mpq_INFTY); k < lp->yjz.nzcnt; k++) {
	mpq_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	mpq_EGlpNumCopyAbs (ay_ij, y_ij);
	if (mpq_EGlpNumIsEqual (y_ij, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_i, mpq_INFTY);
	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && mpq_EGlpNumIsLess (mpq_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum))) {
	    if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY)) {
		mpq_EGlpNumCopyDiff (t_i, *x, *l);
		mpq_EGlpNumAddTo (t_i, *pftol);
		mpq_EGlpNumDivTo (t_i, ay_ij);
	    }
	} else if ((dir == VINCREASE && mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLess (mpq_zeroLpNum, y_ij))) {
	    if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY)) {
		mpq_EGlpNumCopySum (t_i, *u, *pftol);
		mpq_EGlpNumSubTo (t_i, *x);
		mpq_EGlpNumDivTo (t_i, ay_ij);
	    }
	}
	if (mpq_EGlpNumIsEqqual (t_i, mpq_INFTY))
	    continue;

	if (mpq_EGlpNumIsLess (t_i, t_max)) {
	    /* HHH tind = i; yval = fabs (y_ij); tval = t_i -
	       pftol/fabs(y_ij); */
	    mpq_EGlpNumCopy (t_max, t_i);
	}
    }
    /* we use yi_max as temporal variable here */
    mpq_EGlpNumCopyDiff (yi_max, lp->uz[ecol], lp->lz[ecol]);
    if (lp->vtype[ecol] == VBOUNDED && mpq_EGlpNumIsLeq (yi_max, t_max)) {

	mpq_EGlpNumCopy (t_max, yi_max);
	rs->ratio_stat = RATIO_NOBCHANGE;
	mpq_EGlpNumCopy (rs->tz, t_max);
	if (dir != VINCREASE)
	    mpq_EGlpNumSign (rs->tz);
	ILL_CLEANUP;
    }
    if (mpq_EGlpNumIsLeq (mpq_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    /* if (mpq_EGlpNumIsLess (t_max, mpq_zeroLpNum)) printf ("pIIhell\n"); */
    indx = -1;
    mpq_EGlpNumZero (t_z);
    mpq_EGlpNumZero (yi_max);
    mpq_EGlpNumZero (ayi_max);
    ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	mpq_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	mpq_EGlpNumCopyAbs (ay_ij, y_ij);
	if (mpq_EGlpNumIsEqual (y_ij, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_i, mpq_INFTY);
	i = lp->yjz.indx[k];
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if ((dir == VINCREASE && mpq_EGlpNumIsLess (mpq_zeroLpNum, y_ij)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum))) {
	    if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY))
		mpq_EGlpNumCopyDiffRatio (t_i, *x, *l, ay_ij);
	} else if ((dir == VINCREASE && mpq_EGlpNumIsLess (y_ij, mpq_zeroLpNum)) ||
	    (dir == VDECREASE && mpq_EGlpNumIsLess (mpq_zeroLpNum, y_ij))) {
	    if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY))
		mpq_EGlpNumCopyDiffRatio (t_i, *u, *x, ay_ij);
	}
	if (mpq_EGlpNumIsLeq (t_i, t_max)) {
	    if (mpq_EGlpNumIsLess (ayi_max, ay_ij)) {
		mpq_EGlpNumCopy (yi_max, y_ij);
		mpq_EGlpNumCopy (ayi_max, ay_ij);
		indx = i;
		mpq_EGlpNumCopy (t_z, t_i);
		ILL_IFTRACE2 (":%d:%lf:%lf:%lf:%lf", indx, mpq_EGlpNumToLf (t_i),
		    mpq_EGlpNumToLf (t_max), mpq_EGlpNumToLf (ayi_max),
		    mpq_EGlpNumToLf (ay_ij));
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
	mpq_EGlpNumCopy (rs->tz, t_z);
	mpq_EGlpNumCopy (rs->pivotval, yi_max);
	rs->ratio_stat = RATIO_BCHANGE;

	if (dir == VINCREASE)
	    rs->lvstat =
		(mpq_EGlpNumIsLess (mpq_zeroLpNum, yi_max)) ? STAT_LOWER : STAT_UPPER;
	else
	    rs->lvstat =
		(mpq_EGlpNumIsLess (mpq_zeroLpNum, yi_max)) ? STAT_UPPER : STAT_LOWER;

	if (mpq_EGlpNumIsLess (rs->tz, mpq_zeroLpNum)) {
	    ILL_IFTRACE2 ("need to change bound, tz=%la\n", mpq_EGlpNumToLf (rs->tz));
	    mpq_EGlpNumCopyAbs (rs->tz, t_max);
	    mpq_EGlpNumDivUiTo (rs->tz, 10);
	    rs->boundch = 1;
	    mpq_EGlpNumCopy (rs->lbound, lp->xbz[rs->lindex]);
	    if (rs->lvstat == STAT_LOWER)
		mpq_EGlpNumSubInnProdTo (rs->lbound, rs->tz, ayi_max);
	    else
		mpq_EGlpNumAddInnProdTo (rs->lbound, rs->tz, ayi_max);
	}
	if (dir == VDECREASE)
	    mpq_EGlpNumSign (rs->tz);
    }
CLEANUP:
    mpq_ILLfct_update_counts (lp, CNT_PIIPIV, 0, rs->pivotval);
    mpq_EGlpNumClearVar (y_ij);
    mpq_EGlpNumClearVar (ay_ij);
    mpq_EGlpNumClearVar (t_i);
    mpq_EGlpNumClearVar (t_z);
    mpq_EGlpNumClearVar (t_max);
    mpq_EGlpNumClearVar (yi_max);
    mpq_EGlpNumClearVar (ayi_max);
}

#define mpq_GET_XY_DRATIOTEST \
      if (lp->vstat[col] == STAT_UPPER){ \
                mpq_EGlpNumCopyNeg(x,lp->dz[j]);\
        mpq_EGlpNumCopy(y, *zAj);\
      } \
      else{ \
         mpq_EGlpNumCopy(x, lp->dz[j]); \
         mpq_EGlpNumCopyNeg(y, *zAj);\
      } \
      if (lvstat == STAT_UPPER) \
         mpq_EGlpNumSign(y);


void mpq_ILLratio_dI_test (mpq_lpinfo * lp, int lindex, int lvstat,
      mpq_ratio_res * rs)
{
    int j = 0, k;
    int col;
    int cbnd, indx;
    int tctr = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    mpq_t *t = lp->upd.t;
    mpq_t *zAj, x, y, t_j, theta, rcost, delta;
    mpq_t *pftol = &(lp->tol->ip_tol);
    mpq_EGlpNumInitVar (x);
    mpq_EGlpNumInitVar (y);
    mpq_EGlpNumInitVar (t_j);
    mpq_EGlpNumInitVar (theta);
    mpq_EGlpNumInitVar (rcost);
    mpq_EGlpNumInitVar (delta);
    mpq_EGlpNumZero (delta);
    mpq_EGlpNumZero (t_j);
    mpq_EGlpNumZero (rs->tz);
     /* HHH */ rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;
    mpq_EGlpNumZero (rs->pivotval);

    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (mpq_EGlpNumIsEqual (*zAj, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_j, mpq_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	mpq_GET_XY_DRATIOTEST;

	if (mpq_EGlpNumIsLess (y, mpq_zeroLpNum)) {
	    if (lp->dfeas[j] != 0 && lp->vstat[col] != STAT_ZERO) {
		mpq_EGlpNumCopyFrac (t[tctr], x, y);
		ix[tctr] = 10 * k + BBTOLOWER;
		tctr++;
	    } else if (lp->vstat[col] == STAT_ZERO) {
		if (lp->dfeas[j] < 0) {
		    mpq_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BBTOLOWER;
		    tctr++;
		}
		if (lp->dfeas[j] <= 0) {
		    mpq_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BBTOUPPER;
		    tctr++;
		}
	    }
	} else {
	    if (lp->dfeas[j] > 0) {
		if (lp->vstat[col] == STAT_ZERO) {
		    mpq_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BATOUPPER;
		    tctr++;
		    mpq_EGlpNumCopyFrac (t[tctr], x, y);
		    ix[tctr] = 10 * k + BATOLOWER;
		    tctr++;
		}
	    } else if (lp->dfeas[j] == 0) {
		mpq_EGlpNumCopyFrac (t[tctr], x, y);
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
    mpq_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

    mpq_EGlpNumZero (lp->upd.c_obj);
    mpq_EGlpNumCopy (rcost, lp->xbz[lindex]);
    if (lvstat == STAT_LOWER)
	mpq_EGlpNumSign (rcost);
    for (j = 0; j < tctr; j++) {
	cbnd = ix[perm[j]] % 10;
	if (cbnd == BSKIP)
	    continue;

	mpq_EGlpNumCopy (t_j, t[perm[j]]);
	mpq_EGlpNumCopy (x, t_j);
	mpq_EGlpNumSubTo (x, delta);
	mpq_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	mpq_EGlpNumCopy (delta, t_j);
	k = ix[perm[j]] / 10;
	zAj = &(lp->zA.coef[k]);
	indx = lp->zA.indx[k];

	if (lp->vstat[lp->nbaz[indx]] == STAT_LOWER
	    || lp->vstat[lp->nbaz[indx]] == STAT_ZERO)
	    mpq_EGlpNumCopyNeg (theta, *zAj);
	else
	    mpq_EGlpNumCopy (theta, *zAj);

	if (lvstat == STAT_UPPER)
	    mpq_EGlpNumSign (theta);

	switch (cbnd) {
	case BATOLOWER:
	case BATOUPPER:
	    mpq_EGlpNumSubTo (rcost, theta);
	    break;
	case BBTOLOWER:
	case BBTOUPPER:
	    mpq_EGlpNumAddTo (rcost, theta);
	    break;
	}
	if (mpq_EGlpNumIsLeq (rcost, *pftol)) {
	    /* if (t_j < 0.0) printf ("dIhell\n"); */
	    rs->eindex = indx;
	    mpq_EGlpNumCopy (rs->tz, t_j);
	    mpq_EGlpNumCopy (rs->pivotval, *zAj);
	    rs->ratio_stat = RATIO_BCHANGE;
	    ILL_CLEANUP;
	}
    }

CLEANUP:
    mpq_ILLfct_update_counts (lp, CNT_DIPIV, 0, rs->pivotval);
    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
    lp->upd.tctr = tctr;
    lp->upd.i = j;
    mpq_EGlpNumCopyAbs (lp->upd.tz, t_j);
    mpq_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    if (rs->eindex != -1)
	lp->upd.fs = lp->dfeas[rs->eindex];
    mpq_EGlpNumClearVar (x);
    mpq_EGlpNumClearVar (y);
    mpq_EGlpNumClearVar (t_j);
    mpq_EGlpNumClearVar (theta);
    mpq_EGlpNumClearVar (rcost);
    mpq_EGlpNumClearVar (delta);
}

void mpq_ILLratio_dII_test (mpq_lpinfo * lp, int lvstat, mpq_ratio_res * rs)
{
    int j, k, indx;
    int col, ecol;
    mpq_t *zAj, azAj, az_max, x, y, t_j, z_max, t_max, t_z;
    mpq_t *dftol = &(lp->tol->dfeas_tol);
    mpq_EGlpNumInitVar (x);
    mpq_EGlpNumInitVar (y);
    mpq_EGlpNumInitVar (t_j);
    mpq_EGlpNumInitVar (z_max);
    mpq_EGlpNumInitVar (t_max);
    mpq_EGlpNumInitVar (az_max);
    mpq_EGlpNumInitVar (azAj);
    mpq_EGlpNumInitVar (t_z);
    mpq_EGlpNumZero (t_j);
    rs->coeffch = 0;
    mpq_EGlpNumZero (rs->ecoeff);
    rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;
    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
    lp->upd.tctr = 0;
    mpq_EGlpNumZero (lp->upd.dty);
    for (k = 0, mpq_EGlpNumCopy (t_max, mpq_INFTY); k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (mpq_EGlpNumIsEqual (*zAj, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_j, mpq_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	mpq_GET_XY_DRATIOTEST;

	/* warning adding/substracting tolerances to used value, is it right? */
	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, y)) {
	    /* t_j = (x + dftol) / y; */
	    mpq_EGlpNumCopySum (t_j, x, *dftol);
	    mpq_EGlpNumDivTo (t_j, y);
	} else {
	    /* warning adding/substracting tolerances to used value, is it
	       right? */
	    if (lp->vstat[col] == STAT_ZERO)
		mpq_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
	}
	/* if (t_j == mpq_INFTY) */
	if (mpq_EGlpNumIsEqqual (t_j, mpq_INFTY))
	    continue;

	if (mpq_EGlpNumIsLess (t_j, t_max))
	    mpq_EGlpNumCopy (t_max, t_j);
    }

    if (mpq_EGlpNumIsLeq (mpq_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    /* if (t_max < 0.0) printf ("dIIhell\n"); */

    indx = -1;
    mpq_EGlpNumZero (t_z);
    mpq_EGlpNumZero (z_max);
    mpq_EGlpNumZero (az_max);

    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	mpq_EGlpNumCopyAbs (azAj, *zAj);
	if (mpq_EGlpNumIsEqual (*zAj, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_j, mpq_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	mpq_GET_XY_DRATIOTEST;

	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, y) || lp->vstat[col] == STAT_ZERO)
	    mpq_EGlpNumCopyFrac (t_j, x, y);

	if (mpq_EGlpNumIsLeq (t_j, t_max) && (mpq_EGlpNumIsLess (az_max, azAj))) {
	    mpq_EGlpNumCopy (z_max, *zAj);
	    mpq_EGlpNumCopy (az_max, azAj);
	    indx = j;
	    mpq_EGlpNumCopy (t_z, t_j);
	}
    }


    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	rs->eindex = indx;
	mpq_EGlpNumCopy (rs->tz, t_z);
	mpq_EGlpNumCopy (rs->pivotval, z_max);
	rs->ratio_stat = RATIO_BCHANGE;

	if (mpq_EGlpNumIsLess (rs->tz, mpq_zeroLpNum)) {
	    mpq_EGlpNumCopyAbs (rs->tz, t_max);
	    mpq_EGlpNumDivUiTo (rs->tz, 20);
	    rs->coeffch = 1;
	    ecol = lp->nbaz[indx];
	    mpq_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
	    switch (lp->vstat[ecol]) {
	    case STAT_LOWER:
		mpq_EGlpNumAddInnProdTo (rs->ecoeff, rs->tz, az_max);
		break;
	    case STAT_UPPER:
		mpq_EGlpNumSubInnProdTo (rs->ecoeff, rs->tz, az_max);
		break;
	    default:
		mpq_EGlpNumZero (rs->tz);
		break;
	    }
	}
    }

CLEANUP:
    mpq_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
    mpq_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    mpq_EGlpNumClearVar (x);
    mpq_EGlpNumClearVar (y);
    mpq_EGlpNumClearVar (t_j);
    mpq_EGlpNumClearVar (z_max);
    mpq_EGlpNumClearVar (t_max);
    mpq_EGlpNumClearVar (t_z);
    mpq_EGlpNumClearVar (az_max);
    mpq_EGlpNumClearVar (azAj);
}

void mpq_ILLratio_longdII_test (mpq_lpinfo * lp, int lindex, int lvstat,
      mpq_ratio_res * rs)
{
    int j, k, indx = 0, tctr = 0;
    int col, ecol;
    int vs, bnd_exist = 0;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int b_indx = -1;
    mpq_t *t = lp->upd.t;
    mpq_t *l, *u, *xb, *zAj = 0, x, y, t_j, z_max, t_max, t_z, theta, rcost, delta,
      zb_val, tb_val, az_max, azb_val, azAj;
    mpq_t *pftol = &(lp->tol->pfeas_tol);
    mpq_t *dftol = &(lp->tol->dfeas_tol);
    mpq_EGlpNumInitVar (x);
    mpq_EGlpNumInitVar (azAj);
    mpq_EGlpNumInitVar (y);
    mpq_EGlpNumInitVar (t_j);
    mpq_EGlpNumInitVar (z_max);
    mpq_EGlpNumInitVar (az_max);
    mpq_EGlpNumInitVar (t_max);
    mpq_EGlpNumInitVar (t_z);
    mpq_EGlpNumInitVar (theta);
    mpq_EGlpNumInitVar (rcost);
    mpq_EGlpNumInitVar (delta);
    mpq_EGlpNumInitVar (zb_val);
    mpq_EGlpNumInitVar (azb_val);
    mpq_EGlpNumInitVar (tb_val);
    mpq_EGlpNumZero (t_j);
    mpq_EGlpNumZero (delta);
    mpq_EGlpNumZero (zb_val);
    mpq_EGlpNumZero (azb_val);
    mpq_EGlpNumCopy (tb_val, mpq_NINFTY);
    /* #warning not sure about THIS line */
    mpq_EGlpNumZero (rs->pivotval);

    rs->coeffch = 0;
    rs->eindex = -1;
    rs->ratio_stat = RATIO_FAILED;

    ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
    lp->upd.tctr = 0;
    lp->upd.i = 0;
    mpq_EGlpNumZero (lp->upd.tz);
    mpq_EGlpNumZero (lp->upd.piv);
    mpq_EGlpNumZero (lp->upd.c_obj);
    mpq_EGlpNumZero (lp->upd.dty);

    xb = &(lp->xbz[lindex]);
    col = lp->baz[lindex];
    l = &(lp->lz[col]);
    u = &(lp->uz[col]);
    /* rcost = (lvstat == STAT_LOWER) ? l - xb : xb - u; */
    if (lvstat == STAT_LOWER)
	mpq_EGlpNumCopyDiff (rcost, *l, *xb);
    else
	mpq_EGlpNumCopyDiff (rcost, *xb, *u);

    for (k = 0, mpq_EGlpNumCopy (t_max, mpq_INFTY); k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	if (mpq_EGlpNumIsEqual (*zAj, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_j, mpq_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;
	if (lp->vtype[col] == VBOUNDED) {
	    bnd_exist++;
	    continue;
	}
	mpq_GET_XY_DRATIOTEST;

	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, y)) {
	    /* t_j = (x + dftol) / y; */
	    /* #warning Using tolerances to add to result, is it right? */
	    mpq_EGlpNumCopySum (t_j, x, *dftol);
	    mpq_EGlpNumDivTo (t_j, y);
	} else {
	    if (lp->vstat[col] == STAT_ZERO)
		mpq_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
	}
	if (mpq_EGlpNumIsEqqual (t_j, mpq_INFTY))
	    continue;

	if (mpq_EGlpNumIsLess (t_j, t_max))
	    mpq_EGlpNumCopy (t_max, t_j);
    }
    if (mpq_EGlpNumIsLess (t_max, mpq_zeroLpNum)) {
	/* printf ("dIIhell, %.4f\n", t_max); */
	rs->ratio_stat = RATIO_NEGATIVE;
	ILL_CLEANUP;
    }
    if (bnd_exist == 0 && mpq_EGlpNumIsLeq (mpq_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	/*
         * printf ("x = %.8f, b = %.2f \n", lp->xbz[lindex], (lvstat == STAT_LOWER ) ? lp->lz[lp->baz[lindex]] : lp->uz[lp->baz[lindex]]);
         */
	ILL_CLEANUP;
    }
    if (bnd_exist != 0) {
	for (k = 0; k < lp->zA.nzcnt; k++) {
	    zAj = &(lp->zA.coef[k]);
	    if (mpq_EGlpNumIsEqual (*zAj, mpq_zeroLpNum, *pivtol))
		continue;

	    mpq_EGlpNumCopy (t_j, mpq_INFTY);
	    j = lp->zA.indx[k];
	    col = lp->nbaz[j];

	    if (lp->vtype[col] != VBOUNDED)
		continue;

	    mpq_GET_XY_DRATIOTEST;

	    if (mpq_EGlpNumIsLess (mpq_zeroLpNum, y)) {
		mpq_EGlpNumCopyFrac (t_j, x, y);
		if (mpq_EGlpNumIsLeq (t_j, t_max)) {
		    mpq_EGlpNumCopy (t[tctr], t_j);
		    ix[tctr] = k;
		    tctr++;
		}
	    }
	}
    }
    if (tctr != 0) {
	for (j = 0; j < tctr; j++)
	    perm[j] = j;
	mpq_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	for (j = 0; j < tctr; j++) {

	    mpq_EGlpNumCopy (t_j, t[perm[j]]);
	    /* we use x as temporal storage */
	    /* lp->upd.c_obj += (t_j - delta) * rcost; */
	    mpq_EGlpNumCopy (x, t_j);
	    mpq_EGlpNumSubTo (x, delta);
	    mpq_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	    mpq_EGlpNumCopy (delta, t_j);
	     /* HHH */ k = ix[perm[j]];
	    zAj = &(lp->zA.coef[k]);
	    indx = lp->zA.indx[k];
	    col = lp->nbaz[indx];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);
	    vs = lp->vstat[col];
	    /* theta = (vs == STAT_UPPER) ? (l - u) * zAj : (u - l) * zAj; */
	    mpq_EGlpNumCopyDiff (theta, *l, *u);
	    mpq_EGlpNumMultTo (theta, *zAj);
	    if (vs != STAT_UPPER)
		mpq_EGlpNumSign (theta);
	    if (lvstat == STAT_LOWER)
		mpq_EGlpNumAddTo (rcost, theta);
	    else
		mpq_EGlpNumSubTo (rcost, theta);

	    if (mpq_EGlpNumIsLeq (rcost, *pftol)) {
		rs->eindex = indx;
		mpq_EGlpNumCopy (rs->tz, t_j);
		mpq_EGlpNumCopy (rs->pivotval, *zAj);
		rs->ratio_stat = RATIO_BCHANGE;

		if (mpq_EGlpNumIsLess (rs->tz, mpq_zeroLpNum)) {
		    mpq_EGlpNumZero (rs->tz);
		    rs->coeffch = 1;
		    /* rs->ecoeff = lp->cz[col] - lp->dz[indx]; */
		    mpq_EGlpNumCopyDiff (rs->ecoeff, lp->cz[col], lp->dz[indx]);
		    /* lp->upd.c_obj += (rs->tz - delta) * rcost; note ts->tz
		       == 0; */
		    mpq_EGlpNumSubInnProdTo (lp->upd.c_obj, delta, rcost);
		}
		ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
		lp->upd.tctr = tctr;
		lp->upd.i = j;
		mpq_EGlpNumCopy (lp->upd.tz, rs->tz);
		ILL_CLEANUP;
	    }
	}
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
	lp->upd.tctr = tctr;
	lp->upd.i = tctr;
	mpq_EGlpNumCopy (lp->upd.tz, t_j);
	mpq_EGlpNumCopy (zb_val, *zAj);
	mpq_EGlpNumCopyAbs (azb_val, zb_val);
	mpq_EGlpNumCopy (tb_val, t_j);
	b_indx = indx;
    }
    if (bnd_exist != 0 && mpq_EGlpNumIsLeq (mpq_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	/* printf ("rcost: %.8f\n", rcost); */
	ILL_CLEANUP;
    }
    mpq_EGlpNumZero (z_max);
    mpq_EGlpNumZero (az_max);
    indx = -1;
    mpq_EGlpNumZero (t_z);
    for (k = 0; k < lp->zA.nzcnt; k++) {
	zAj = &(lp->zA.coef[k]);
	mpq_EGlpNumCopyAbs (azAj, *zAj);
	if (mpq_EGlpNumIsEqual (*zAj, mpq_zeroLpNum, *pivtol))
	    continue;

	mpq_EGlpNumCopy (t_j, mpq_INFTY);
	j = lp->zA.indx[k];
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED ||
	    lp->vtype[col] == VBOUNDED)
	    continue;

	mpq_GET_XY_DRATIOTEST;

	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, y) || lp->vstat[col] == STAT_ZERO)
	    mpq_EGlpNumCopyFrac (t_j, x, y);

	if (mpq_EGlpNumIsLeq (t_j, t_max)) {
	    if (mpq_EGlpNumIsLess (az_max, azAj)) {
		mpq_EGlpNumCopy (z_max, *zAj);
		mpq_EGlpNumCopy (az_max, azAj);
		indx = j;
		mpq_EGlpNumCopy (t_z, t_j);
	    }
	}
    }

    if (indx < 0) {
	rs->ratio_stat = RATIO_FAILED;
	ILL_CLEANUP;
    }
    if ((tctr == 0) || (mpq_EGlpNumIsLess (tb_val, mpq_zeroLpNum)) ||
	(tctr != 0 && mpq_EGlpNumIsLeq (tb_val, t_z) &&
	    mpq_EGlpNumIsLeq (azb_val, az_max))) {
	/* we use x as temporal vvariable */
	/* lp->upd.c_obj += (t_z - delta) * rcost; */
	mpq_EGlpNumCopyDiff (x, t_z, delta);
	mpq_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	mpq_EGlpNumCopy (delta, t_z);
	rs->eindex = indx;
	mpq_EGlpNumCopy (rs->tz, t_z);
	mpq_EGlpNumCopy (rs->pivotval, z_max);
	rs->ratio_stat = RATIO_BCHANGE;
    }
    /* For now */
    else if (tctr != 0) {
	rs->eindex = b_indx;
	mpq_EGlpNumCopy (rs->tz, tb_val);
	mpq_EGlpNumCopy (rs->pivotval, zb_val);
	rs->ratio_stat = RATIO_BCHANGE;
	lp->upd.i -= 1;
    }
    if (mpq_EGlpNumIsLess (rs->tz, mpq_zeroLpNum)) {
	/* if (tctr != 0) printf ("despite long step\n"); */
	/* rs->tz = fabs (t_max / 20.0); */
	mpq_EGlpNumCopyAbs (rs->tz, t_max);
	mpq_EGlpNumDivUiTo (rs->tz, 20);
	rs->coeffch = 1;

	ecol = lp->nbaz[indx];
	if (lp->vstat[ecol] == STAT_LOWER) {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs
	       (z_max); */
	    mpq_EGlpNumCopy (rs->ecoeff, az_max);
	    mpq_EGlpNumMultTo (rs->ecoeff, rs->tz);
	    mpq_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
	    mpq_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
	} else if (lp->vstat[ecol] == STAT_UPPER) {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs
	       (z_max); */
	    mpq_EGlpNumCopy (rs->ecoeff, az_max);
	    mpq_EGlpNumMultTo (rs->ecoeff, rs->tz);
	    mpq_EGlpNumSign (rs->ecoeff);
	    mpq_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
	    mpq_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
	} else {
	    /* rs->ecoeff = lp->cz[ecol] - lp->dz[indx]; */
	    mpq_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
	    mpq_EGlpNumZero (rs->tz);
	}
	/* we use x as temporal storage */
	/* lp->upd.c_obj += (rs->tz - delta) * rcost; */
	mpq_EGlpNumCopy (x, rs->tz);
	mpq_EGlpNumSubTo (x, delta);
	mpq_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
    }
CLEANUP:
    mpq_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
    mpq_EGlpNumCopy (lp->upd.piv, rs->pivotval);
    mpq_EGlpNumClearVar (x);
    mpq_EGlpNumClearVar (y);
    mpq_EGlpNumClearVar (t_j);
    mpq_EGlpNumClearVar (z_max);
    mpq_EGlpNumClearVar (az_max);
    mpq_EGlpNumClearVar (t_max);
    mpq_EGlpNumClearVar (t_z);
    mpq_EGlpNumClearVar (theta);
    mpq_EGlpNumClearVar (rcost);
    mpq_EGlpNumClearVar (delta);
    mpq_EGlpNumClearVar (zb_val);
    mpq_EGlpNumClearVar (azb_val);
    mpq_EGlpNumClearVar (tb_val);
    mpq_EGlpNumClearVar (azAj);
}

void mpq_ILLratio_pivotin_test (mpq_lpinfo * lp, int *rlist, int rcnt,
      mpq_ratio_res * rs)
{
    int i, k, col;
    mpq_t *x, *l, *u;
    mpq_t ay_ij, at_i, at_l, at_u, ayi_max, y_ij, t_i, t_l, t_u, t_max, yi_max;
    if (rcnt <= 0 || rs == NULL)
	return;
    mpq_EGlpNumInitVar (ay_ij);
    mpq_EGlpNumInitVar (at_i);
    mpq_EGlpNumInitVar (at_l);
    mpq_EGlpNumInitVar (at_u);
    mpq_EGlpNumInitVar (ayi_max);
    mpq_EGlpNumInitVar (t_max);
    mpq_EGlpNumInitVar (y_ij);
    mpq_EGlpNumInitVar (t_i);
    mpq_EGlpNumInitVar (t_l);
    mpq_EGlpNumInitVar (t_u);
    mpq_EGlpNumInitVar (yi_max);
    rs->boundch = 0;
    rs->lindex = -1;
    mpq_EGlpNumZero (rs->tz);
    rs->ratio_stat = RATIO_FAILED;
    rs->lvstat = -1;
    mpq_EGlpNumZero (rs->pivotval);
    mpq_EGlpNumZero (rs->lbound);

    for (i = 0; i < rcnt; i++)
	lp->iwork[rlist[i]] = 1;

    for (k = 0, mpq_EGlpNumCopy (t_max, mpq_INFTY); k < lp->yjz.nzcnt; k++) {
	mpq_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	if (mpq_EGlpNumIsEqual (y_ij, mpq_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	if (lp->iwork[lp->baz[i]] == 1)
	    continue;
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	mpq_EGlpNumCopy (t_u, mpq_INFTY);
	mpq_EGlpNumCopy (at_u, mpq_INFTY);
	mpq_EGlpNumCopy (t_l, mpq_NINFTY);
	mpq_EGlpNumCopy (at_l, mpq_INFTY);

	if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY)) {
	    mpq_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
	    mpq_EGlpNumCopyAbs (at_l, t_l);
	    if (mpq_EGlpNumIsLess (at_l, t_max))
		mpq_EGlpNumCopy (t_max, at_l);
	}
	if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY)) {
	    mpq_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
	    mpq_EGlpNumCopyAbs (at_u, t_u);
	    if (mpq_EGlpNumIsLess (at_u, t_max))
		mpq_EGlpNumCopy (t_max, at_u);
	}
    }

    if (mpq_EGlpNumIsLeq (mpq_INFTY, t_max)) {
	rs->ratio_stat = RATIO_UNBOUNDED;
	ILL_CLEANUP;
    }
    mpq_EGlpNumZero (yi_max);
    mpq_EGlpNumZero (ayi_max);
    mpq_EGlpNumMultUiTo (t_max, 101);
    mpq_EGlpNumDivUiTo (t_max, 100);
    for (k = 0; k < lp->yjz.nzcnt; k++) {
	mpq_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
	mpq_EGlpNumCopyAbs (ay_ij, y_ij);
	if (mpq_EGlpNumIsEqual (y_ij, mpq_zeroLpNum, *pivtol))
	    continue;

	i = lp->yjz.indx[k];
	if (lp->iwork[lp->baz[i]] == 1)
	    continue;
	x = &(lp->xbz[i]);
	col = lp->baz[i];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	mpq_EGlpNumCopy (t_u, mpq_INFTY);
	mpq_EGlpNumCopy (at_u, t_u);
	mpq_EGlpNumCopy (t_l, mpq_NINFTY);
	mpq_EGlpNumCopy (at_l, t_u);
	if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY)) {
	    mpq_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
	    mpq_EGlpNumCopyAbs (at_l, t_l);
	}
	if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY)) {
	    mpq_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
	    mpq_EGlpNumCopyAbs (at_u, t_u);
	}
	/* t_i = (fabs (t_l) < fabs (t_u)) ? t_l : t_u; */
	if (mpq_EGlpNumIsLess (at_l, at_u)) {
	    mpq_EGlpNumCopy (t_i, t_l);
	    mpq_EGlpNumCopy (at_i, at_l);
	} else {
	    mpq_EGlpNumCopy (t_i, t_u);
	    mpq_EGlpNumCopy (at_i, at_u);
	}
	/* if (fabs (t_i) <= t_max + t_max * (1.0e-2)) */
	if (mpq_EGlpNumIsLeq (at_i, t_max)) {
	    if (mpq_EGlpNumIsLess (ayi_max, ay_ij)) {
		mpq_EGlpNumCopy (yi_max, y_ij);
		mpq_EGlpNumCopy (ayi_max, ay_ij);
		rs->lindex = i;
		mpq_EGlpNumCopy (rs->tz, t_i);
		rs->lvstat = (mpq_EGlpNumIsLess (at_l, at_u)) ? STAT_LOWER : STAT_UPPER;
	    }
	}
    }

    if (rs->lindex < 0) {
	rs->ratio_stat = RATIO_FAILED;
    } else {
	rs->ratio_stat = RATIO_BCHANGE;
	mpq_EGlpNumCopy (rs->pivotval, yi_max);
    }
CLEANUP:
    for (i = 0; i < rcnt; i++)
	lp->iwork[rlist[i]] = 0;
    mpq_EGlpNumClearVar (t_max);
    mpq_EGlpNumClearVar (ay_ij);
    mpq_EGlpNumClearVar (at_i);
    mpq_EGlpNumClearVar (at_l);
    mpq_EGlpNumClearVar (at_u);
    mpq_EGlpNumClearVar (ayi_max);
    mpq_EGlpNumClearVar (y_ij);
    mpq_EGlpNumClearVar (t_i);
    mpq_EGlpNumClearVar (t_l);
    mpq_EGlpNumClearVar (t_u);
    mpq_EGlpNumClearVar (yi_max);
    return;
}
