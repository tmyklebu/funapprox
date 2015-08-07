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

/* RCS_INFO = "$RCSfile: dbl_fct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#define dbl_FCT_DEBUG 0

#include "econfig.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdefs.h"
#include "stddefs.h"
#include "dbl_basis.h"
#include "dbl_fct.h"
#include "dbl_price.h"
#include "dbl_ratio.h"
#include "dbl_dstruct.h"

dbl_bndinfo *dbl_ILLfct_new_bndinfo (void)
{
    dbl_bndinfo *nbnd = (dbl_bndinfo *) malloc (sizeof (dbl_bndinfo));
    if (!nbnd) {
	fprintf (stderr, "not enough memory, in %s\n", __func__);
	exit (1);
    }
    dbl_EGlpNumInitVar ((nbnd->pbound));
    dbl_EGlpNumInitVar ((nbnd->cbound));
    return nbnd;
}

void dbl_ILLfct_free_bndinfo (dbl_bndinfo * binfo)
{
    dbl_EGlpNumClearVar ((binfo->pbound));
    dbl_EGlpNumClearVar ((binfo->cbound));
    ILL_IFFREE (binfo, dbl_bndinfo);
    return;
}

static int dbl_compute_zA1 (dbl_lpinfo * lp,
      dbl_svector * z,
      dbl_svector * zA,
      double ztoler),
/*
  compute_zA2 (dbl_lpinfo * lp,
                             dbl_svector * z,
                             dbl_svector * zA,
                             const double* ztoler), */
  dbl_compute_zA3 (dbl_lpinfo * lp,
      dbl_svector * z,
      dbl_svector * zA,
      double ztoler),
  dbl_expand_var_bounds (dbl_lpinfo * lp,
      double ftol,
      int *chgb),
  dbl_expand_var_coefs (dbl_lpinfo * lp,
      double ftol,
      int *chgc);

static void dbl_update_piv_values (dbl_count_struct * c,
      int phase,
      double piv),
/* copy_vectors (dbl_svector * a, dbl_svector * b), */
  dbl_add_vectors (dbl_lpinfo * lp,
      dbl_svector * a,
      dbl_svector * b,
      dbl_svector * c,
      double t);

static double dbl_my_rand (int bound,
      ILLrandstate * r);


void dbl_ILLfct_load_workvector (dbl_lpinfo * lp,
      dbl_svector * s)
{
    int i;

    for (i = 0; i < s->nzcnt; i++) {
	lp->work.indx[i] = s->indx[i];
	dbl_EGlpNumCopy (lp->work.coef[s->indx[i]], s->coef[i]);
    }
    lp->work.nzcnt = s->nzcnt;
}

void dbl_ILLfct_zero_workvector (dbl_lpinfo * lp)
{
    int i;

    for (i = 0; i < lp->work.nzcnt; i++)
	dbl_EGlpNumZero (lp->work.coef[lp->work.indx[i]]);
    lp->work.nzcnt = 0;
}

void dbl_ILLfct_set_variable_type (dbl_lpinfo * lp)
{
    int j;

    for (j = 0; j < lp->ncols; j++) {

	if (lp->matcnt[j] == 1 && lp->O->rowmap[lp->matind[lp->matbeg[j]]] == j)
	    lp->vclass[j] = CLASS_LOGICAL;
	else
	    lp->vclass[j] = CLASS_STRUCT;
	switch ((dbl_EGlpNumIsEqqual (lp->uz[j], dbl_INFTY) ? 1U : 0U) |
	    (dbl_EGlpNumIsEqqual (lp->lz[j], dbl_NINFTY) ? 2U : 0U)) {
	case 0:
	    if (dbl_EGlpNumIsLess (lp->lz[j], lp->uz[j]))
		lp->vtype[j] = VBOUNDED;
	    else if (dbl_EGlpNumIsEqqual (lp->lz[j], dbl_zeroLpNum) &&
		(lp->vclass[j] == CLASS_LOGICAL))
		lp->vtype[j] = VARTIFICIAL;
	    else
		lp->vtype[j] = VFIXED;
	    break;
	case 3:
	    lp->vtype[j] = VFREE;
	    break;
	case 1:
	    lp->vtype[j] = VLOWER;
	    break;
	case 2:
	    lp->vtype[j] = VUPPER;
	    break;
	}
    }
}

/* compute various vectors */

void dbl_ILLfct_compute_pobj (dbl_lpinfo * lp)
{
    int i, j;
    int col;
    double sum;
    dbl_EGlpNumInitVar (sum);
    dbl_EGlpNumZero (sum);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumAddInnProdTo (sum, lp->cz[lp->baz[i]], lp->xbz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	if (lp->vstat[col] == STAT_UPPER)
	    dbl_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER)
	    dbl_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->lz[col]);
    }
    dbl_EGlpNumCopy (lp->pobjval, sum);
    dbl_EGlpNumCopy (lp->objval, sum);
    dbl_EGlpNumClearVar (sum);
}

void dbl_ILLfct_compute_dobj (dbl_lpinfo * lp)
{
    int i, j;
    int col;
    double sum;
    dbl_EGlpNumInitVar (sum);
    dbl_EGlpNumZero (sum);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumAddInnProdTo (sum, lp->piz[i], lp->bz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	if (lp->vstat[col] == STAT_UPPER)
	    dbl_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER)
	    dbl_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->lz[col]);
    }
    dbl_EGlpNumCopy (lp->dobjval, sum);
    dbl_EGlpNumCopy (lp->objval, sum);
    dbl_EGlpNumClearVar (sum);
}

void dbl_ILLfct_compute_xbz (dbl_lpinfo * lp)
{
    int i, j, r;
    int col, mcnt, mbeg;
    dbl_svector *srhs = &(lp->srhs);
    dbl_svector *ssoln = &(lp->ssoln);
    double xval;
    dbl_EGlpNumInitVar (xval);

    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumZero (lp->xbz[i]);
	dbl_EGlpNumCopy (srhs->coef[i], lp->bz[i]);
    }
    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	dbl_EGlpNumZero (xval);
	if (lp->vstat[col] == STAT_UPPER && dbl_EGlpNumIsNeqqZero (lp->uz[col]))
	    dbl_EGlpNumCopy (xval, lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER && dbl_EGlpNumIsNeqqZero (lp->lz[col]))
	    dbl_EGlpNumCopy (xval, lp->lz[col]);

	if (dbl_EGlpNumIsNeqqZero (xval)) {
	    mcnt = lp->matcnt[col];
	    mbeg = lp->matbeg[col];
	    for (i = 0; i < mcnt; i++)
		dbl_EGlpNumSubInnProdTo (srhs->coef[lp->matind[mbeg + i]], xval,
		    lp->matval[mbeg + i]);
	}
    }
    for (i = 0, r = 0; i < lp->nrows; i++)
	if (dbl_EGlpNumIsNeqqZero (srhs->coef[i])) {
	    dbl_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
	    srhs->indx[r] = i;
	    r++;
	}
    srhs->nzcnt = r;

    dbl_ILLbasis_column_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	dbl_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
    dbl_EGlpNumClearVar (xval);
}

void dbl_ILLfct_compute_piz (dbl_lpinfo * lp)
{
    int i, r;
    dbl_svector *srhs = &(lp->srhs);
    dbl_svector *ssoln = &(lp->ssoln);

    for (i = 0, r = 0; i < lp->nrows; i++) {
	dbl_EGlpNumZero (lp->piz[i]);
	if (dbl_EGlpNumIsNeqqZero (lp->cz[lp->baz[i]])) {
	    srhs->indx[r] = i;
	    dbl_EGlpNumCopy (srhs->coef[r], lp->cz[lp->baz[i]]);
	    r++;
	}
    }
    srhs->nzcnt = r;

    dbl_ILLbasis_row_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	dbl_EGlpNumCopy (lp->piz[ssoln->indx[i]], ssoln->coef[i]);
}

void dbl_ILLfct_compute_dz (dbl_lpinfo * lp)
{
    int i, j;
    int col;
    int mcnt, mbeg;
    double sum;
    dbl_EGlpNumInitVar (sum);

    for (j = 0; j < lp->nnbasic; j++) {
	dbl_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    dbl_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	dbl_EGlpNumCopyDiff (lp->dz[j], lp->cz[col], sum);
    }
    dbl_EGlpNumClearVar (sum);
}

void dbl_ILLfct_compute_phaseI_xbz (dbl_lpinfo * lp)
{
    int i, j, r;
    int col, mcnt, mbeg;
    dbl_svector *srhs = &(lp->srhs);
    dbl_svector *ssoln = &(lp->ssoln);

    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumZero (lp->xbz[i]);
	dbl_EGlpNumZero (srhs->coef[i]);
    }
    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];

	if (lp->dfeas[j]) {
	    mcnt = lp->matcnt[col];
	    mbeg = lp->matbeg[col];
	    if (lp->dfeas[j] == -1)
		for (i = 0; i < mcnt; i++)
		    dbl_EGlpNumSubTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
	    else
		for (i = 0; i < mcnt; i++)
		    dbl_EGlpNumAddTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
	}
    }
    for (i = 0, r = 0; i < lp->nrows; i++)
	if (dbl_EGlpNumIsNeqqZero (srhs->coef[i])) {
	    dbl_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
	    srhs->indx[r] = i;
	    r++;
	}
    srhs->nzcnt = r;

    dbl_ILLbasis_column_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	dbl_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
}

void dbl_ILLfct_compute_phaseI_piz (dbl_lpinfo * lp)
{
    int i, r;
    dbl_svector *srhs = &(lp->srhs);
    dbl_svector *ssoln = &(lp->ssoln);

    for (i = 0, r = 0; i < lp->nrows; i++) {
	dbl_EGlpNumZero (lp->pIpiz[i]);
	if (lp->bfeas[i] != 0) {
	    srhs->indx[r] = i;
	    dbl_EGlpNumSet (srhs->coef[r], (double) lp->bfeas[i]);
	    r++;
	}
    }
    srhs->nzcnt = r;

    dbl_ILLbasis_row_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	dbl_EGlpNumCopy (lp->pIpiz[ssoln->indx[i]], ssoln->coef[i]);
    dbl_ILLfct_update_counts (lp, CNT_P1PINZ, ssoln->nzcnt, dbl_zeroLpNum);
}

void dbl_ILLfct_compute_phaseI_dz (dbl_lpinfo * lp)
{
    int i, j;
    int col;
    int mcnt, mbeg;
    double sum;
    dbl_EGlpNumInitVar (sum);
    ILL_IFTRACE ("%s\n", __func__);

    for (j = 0; j < lp->nnbasic; j++) {
	dbl_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    dbl_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	dbl_EGlpNumCopyNeg (lp->pIdz[j], sum);
	ILL_IFTRACE ("%d:%d:%lf:%la\n", j, col, dbl_EGlpNumToLf (sum),
	    dbl_EGlpNumToLf (sum));
    }
    dbl_EGlpNumClearVar (sum);
}

void dbl_ILLfct_compute_yz (dbl_lpinfo * lp,
      dbl_svector * yz,
      dbl_svector * updz,
      int col)
{
    dbl_svector a;

    a.nzcnt = lp->matcnt[col];
    a.indx = &(lp->matind[lp->matbeg[col]]);
    a.coef = &(lp->matval[lp->matbeg[col]]);

    dbl_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, dbl_PIVZ_TOLER);
    if (updz)
	dbl_ILLbasis_column_solve_update (lp, &a, updz, yz);
    else
	dbl_ILLbasis_column_solve (lp, &a, yz);
    dbl_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, dbl_SZERO_TOLER);
}

void dbl_ILLfct_compute_zz (dbl_lpinfo * lp,
      dbl_svector * zz,
      int row)
{
    dbl_ILLfct_compute_binvrow (lp, zz, row, dbl_PIVZ_TOLER);
}

void dbl_ILLfct_compute_binvrow (dbl_lpinfo * lp,
      dbl_svector * zz,
      int row,
      double ztoler)
{
    dbl_svector a;
    double e;
    dbl_EGlpNumInitVar (e);
    dbl_EGlpNumOne (e);

    a.nzcnt = 1;
    a.coef = &e;
    a.indx = &row;

    if (dbl_EGlpNumIsLess (dbl_zeroLpNum, ztoler))
	dbl_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, ztoler);
    dbl_ILLbasis_row_solve (lp, &a, zz);
    if (dbl_EGlpNumIsLess (dbl_zeroLpNum, ztoler))
	dbl_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, dbl_SZERO_TOLER);
    dbl_EGlpNumClearVar (e);
}

void dbl_ILLfct_compute_psteep_upv (dbl_lpinfo * lp,
      dbl_svector * swz)
{
    dbl_ILLbasis_row_solve (lp, &(lp->yjz), swz);
}

void dbl_ILLfct_compute_dsteep_upv (dbl_lpinfo * lp,
      dbl_svector * swz)
{
    dbl_ILLbasis_column_solve (lp, &(lp->zz), swz);
}

static int dbl_compute_zA1 (dbl_lpinfo * lp,
      dbl_svector * z,
      dbl_svector * zA,
      double ztoler)
{
    int rval = 0;
    int i, j, nz = 0;
    int col, mcnt, mbeg;
    double sum;
    double *v = 0;
    dbl_EGlpNumInitVar (sum);
    v = dbl_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumZero (v[i]);
    for (i = 0; i < z->nzcnt; i++)
	dbl_EGlpNumCopy (v[z->indx[i]], z->coef[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	dbl_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    dbl_EGlpNumAddInnProdTo (sum, v[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

	if (dbl_EGlpNumIsNeqZero (sum, ztoler)) {
	    dbl_EGlpNumCopy (zA->coef[nz], sum);
	    zA->indx[nz] = j;
	    nz++;
	}
    }
    zA->nzcnt = nz;

    dbl_EGlpNumClearVar (sum);
    dbl_EGlpNumFreeArray (v);
    ILL_RETURN (rval, "dbl_compute_zA1");
}


static int dbl_compute_zA3 (dbl_lpinfo * lp,
      dbl_svector * z,
      dbl_svector * zA,
      double ztoler)
{
    int rval = 0;
    int i, j, k, ix;
    int nz = 0;
    int row, col;
    int rcnt, rbeg;
    double val;
    dbl_EGlpNumInitVar (val);
    k = 0;
    for (i = 0; i < z->nzcnt; i++) {
	row = z->indx[i];
	dbl_EGlpNumCopy (val, z->coef[i]);
	rcnt = lp->rowcnt[row];
	rbeg = lp->rowbeg[row];
	for (j = 0; j < rcnt; j++) {
	    col = lp->rowind[rbeg + j];
	    if (lp->vstat[col] != STAT_BASIC) {
		ix = lp->vindex[col];
		if (lp->iwork[ix] == 0) {
		    lp->iwork[ix] = 1;
		    lp->work.indx[k++] = ix;
		}
		dbl_EGlpNumAddInnProdTo (lp->work.coef[ix], val, lp->rowval[rbeg + j]);
	    }
	}
    }
    for (j = 0; j < k; j++) {
	ix = lp->work.indx[j];
	dbl_EGlpNumCopy (val, lp->work.coef[ix]);
	dbl_EGlpNumZero (lp->work.coef[ix]);
	lp->iwork[ix] = 0;
	if (dbl_EGlpNumIsNeqZero (val, ztoler)) {
	    dbl_EGlpNumCopy (zA->coef[nz], val);
	    zA->indx[nz] = ix;
	    nz++;
	}
    }
    zA->nzcnt = nz;
    dbl_EGlpNumClearVar (val);
    ILL_RETURN (rval, "dbl_compute_zA3");
}

int dbl_ILLfct_compute_zA (dbl_lpinfo * lp,
      dbl_svector * z,
      dbl_svector * zA)
{
    if (z->nzcnt < lp->nrows / 2)
	return dbl_compute_zA3 (lp, z, zA, dbl_PIVZ_TOLER);
    else
	return dbl_compute_zA1 (lp, z, zA, dbl_PIVZ_TOLER);
}

/* compute v^T A */
void dbl_ILLfct_compute_vA (dbl_lpinfo * lp,
      dbl_svector * v,
      double *vA)
{
    int i, j;
    int row, col;
    int rcnt, rbeg;
    double val;
    dbl_EGlpNumInitVar (val);

    for (j = 0; j < lp->ncols; j++)
	dbl_EGlpNumZero (vA[j]);

    for (i = 0; i < v->nzcnt; i++) {
	row = v->indx[i];
	dbl_EGlpNumCopy (val, v->coef[i]);
	rcnt = lp->rowcnt[row];
	rbeg = lp->rowbeg[row];
	for (j = 0; j < rcnt; j++) {
	    col = lp->rowind[rbeg + j];
	    dbl_EGlpNumAddInnProdTo (vA[col], val, lp->rowval[rbeg + j]);
	}
    }

    for (j = 0; j < lp->ncols; j++)
	if (dbl_EGlpNumIsEqual (vA[j], dbl_zeroLpNum, dbl_SZERO_TOLER))
	    dbl_EGlpNumZero (vA[j]);

    dbl_EGlpNumClearVar (val);
    return;
}

/* update information */

/*
1) lvstat - new status of leaving var.
*/
void dbl_ILLfct_update_basis_info (dbl_lpinfo * lp,
      int eindex,
      int lindex,
      int lvstat)
{
    int evar;
    int lvar;

    evar = lp->nbaz[eindex];

    if (lindex >= 0) {		/* variable leaves basis */
	lvar = lp->baz[lindex];
	lp->vstat[evar] = STAT_BASIC;
	lp->vstat[lvar] = lvstat;
	lp->vindex[evar] = lindex;
	lp->vindex[lvar] = eindex;
	lp->baz[lindex] = evar;
	lp->nbaz[eindex] = lvar;
	(lp->basisid)++;
    } else {
	lp->vstat[evar] = (lp->vstat[evar] == STAT_LOWER) ? STAT_UPPER : STAT_LOWER;
    }
}

void dbl_ILLfct_update_xz (dbl_lpinfo * lp,
      double tz,
      int eindex,
      int lindex)
{
    int i, evar, estat;
    ILL_IFTRACE ("%s:%la:%d:%d:%d\n", __func__, dbl_EGlpNumToLf (tz), eindex,
	lindex, lp->yjz.nzcnt);

    if (dbl_EGlpNumIsNeqqZero (tz))
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    dbl_EGlpNumSubInnProdTo (lp->xbz[lp->yjz.indx[i]], tz, lp->yjz.coef[i]);

    if (lindex >= 0) {		/* variable leaves basis */
	evar = lp->nbaz[eindex];
	estat = lp->vstat[evar];
	if (estat == STAT_LOWER)
	    dbl_EGlpNumCopySum (lp->xbz[lindex], lp->lz[evar], tz);
	else if (estat == STAT_UPPER)
	    dbl_EGlpNumCopySum (lp->xbz[lindex], lp->uz[evar], tz);
	else if (estat == STAT_ZERO)
	    dbl_EGlpNumCopy (lp->xbz[lindex], tz);
    }
}

void dbl_ILLfct_update_piz (dbl_lpinfo * lp,
      double alpha)
{
    int i;

    for (i = 0; i < lp->zz.nzcnt; i++)
	dbl_EGlpNumAddInnProdTo (lp->piz[lp->zz.indx[i]], alpha, lp->zz.coef[i]);
}

void dbl_ILLfct_update_pIpiz (dbl_lpinfo * lp,
      dbl_svector * z,
      double alpha)
{
    int i;
    if (dbl_EGlpNumIsEqqual (alpha, dbl_zeroLpNum))
	return;
    if (dbl_EGlpNumIsEqqual (alpha, dbl_oneLpNum)) {
	for (i = 0; i < z->nzcnt; i++)
	    dbl_EGlpNumAddTo (lp->pIpiz[z->indx[i]], z->coef[i]);
    } else {
	for (i = 0; i < z->nzcnt; i++)
	    dbl_EGlpNumAddInnProdTo (lp->pIpiz[z->indx[i]], alpha, z->coef[i]);
    }
}

void dbl_ILLfct_update_dz (dbl_lpinfo * lp,
      int eindex,
      double alpha)
{
    int i;

    for (i = 0; i < lp->zA.nzcnt; i++)
	dbl_EGlpNumSubInnProdTo (lp->dz[lp->zA.indx[i]], alpha, lp->zA.coef[i]);
    dbl_EGlpNumCopyNeg (lp->dz[eindex], alpha);
}

void dbl_ILLfct_update_pIdz (dbl_lpinfo * lp,
      dbl_svector * zA,
      int eindex,
      double alpha)
{
    int i;
    if (dbl_EGlpNumIsEqqual (alpha, dbl_zeroLpNum))
	return;

    if (dbl_EGlpNumIsEqqual (alpha, dbl_oneLpNum)) {
	for (i = 0; i < zA->nzcnt; i++)
	    dbl_EGlpNumSubTo (lp->pIdz[zA->indx[i]], zA->coef[i]);
    } else {
	for (i = 0; i < zA->nzcnt; i++)
	    dbl_EGlpNumSubInnProdTo (lp->pIdz[zA->indx[i]], alpha, zA->coef[i]);
    }
    if (eindex > -1)
	dbl_EGlpNumCopyNeg (lp->pIdz[eindex], alpha);
}

/* bound and coef shift routines */

/* scale bound in dbl_my_rand to get more random digits, unless bound is
   large */
static double dbl_my_rand (int bound,
      ILLrandstate * r)
{
    int k = bound, scale = 1;
    double v = 0.0;

    if (bound < 100000) {
	k = 20000 * bound;
	scale = 20000;
    }
    v = 1 + (ILLutil_lprand (r) % (k));
    return v / (double) scale;
}

static int dbl_expand_var_bounds (dbl_lpinfo * lp,
      double ftol,
      int *chgb)
{
    int rval = 0;
    int i, col, nchg = 0;
    double newb, cftol;
    double *x, *l, *u;
    ILLrandstate r;
    dbl_EGlpNumInitVar (newb);
    dbl_EGlpNumInitVar (cftol);
    dbl_EGlpNumCopyAbs (cftol, ftol);
    dbl_EGlpNumDivUiTo (cftol, 10);

    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
	    continue;
	x = &(lp->xbz[i]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	/* we use newb as temporal variable outside the if's scope */
	dbl_EGlpNumCopyDiff (newb, *x, ftol);
	if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY) && dbl_EGlpNumIsLess (newb, *l)) {
	    dbl_EGlpNumSet (newb, -1.0 * (dbl_my_rand (50, &(lp->rstate)) + 1.0));
	    dbl_EGlpNumMultTo (newb, cftol);
	    if (dbl_EGlpNumIsLess (*x, *l))
		dbl_EGlpNumAddTo (newb, *x);
	    else
		dbl_EGlpNumAddTo (newb, *l);
	    rval = dbl_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
	dbl_EGlpNumCopySum (newb, *x, ftol);
	if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY) && dbl_EGlpNumIsLess (*u, newb)) {
	    dbl_EGlpNumSet (newb, dbl_my_rand (50, &(lp->rstate)) + 1.0);
	    dbl_EGlpNumMultTo (newb, cftol);
	    if (dbl_EGlpNumIsLess (*x, *u))
		dbl_EGlpNumAddTo (newb, *u);
	    else
		dbl_EGlpNumAddTo (newb, *x);
	    rval = dbl_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
    }
    *chgb = nchg;

CLEANUP:
    dbl_EGlpNumClearVar (newb);
    dbl_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "dbl_expand_var_bounds");
}

static int dbl_expand_phaseI_bounds (dbl_lpinfo * lp,
      int *chgb)
{
    int rval = 0;
    int i, col, nchg = 0;
    double newb, cftol;
    double *u, *l, *x;
    ILLrandstate r;
    dbl_EGlpNumInitVar (newb);
    dbl_EGlpNumInitVar (cftol);
    dbl_EGlpNumCopyAbs (cftol, lp->tol->ip_tol);
    dbl_EGlpNumDivUiTo (cftol, 10);
    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
	    continue;
	x = &(lp->xbz[i]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY) && dbl_EGlpNumIsEqual (*x, *l, cftol)) {
	    dbl_EGlpNumSet (newb, dbl_my_rand (50, &(lp->rstate)) + 1.0);
	    dbl_EGlpNumMultTo (newb, cftol);
	    dbl_EGlpNumSign (newb);
	    dbl_EGlpNumAddTo (newb, *l);
	    rval = dbl_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
	if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY) && dbl_EGlpNumIsEqual (*x, *u, cftol)) {
	    dbl_EGlpNumSet (newb, dbl_my_rand (50, &(lp->rstate)) + 1.0);
	    dbl_EGlpNumMultTo (newb, cftol);
	    dbl_EGlpNumAddTo (newb, *u);
	    rval = dbl_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
    }
    *chgb = nchg;

CLEANUP:
    dbl_EGlpNumClearVar (newb);
    dbl_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "dbl_expand_phaseI_bounds");
}

int dbl_ILLfct_adjust_viol_bounds (dbl_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;
    double tol;
    dbl_EGlpNumInitVar (tol);
    dbl_EGlpNumCopyNeg (tol, lp->tol->pfeas_tol);
    rval = dbl_expand_var_bounds (lp, tol, &chgb);
#if dbl_FCT_DEBUG > 0
    if (rval == 0)
	printf ("adjusting %d bounds\n", chgb);
#endif
    dbl_EGlpNumClearVar (tol);
    ILL_RETURN (rval, "dbl_ILLfct_adjust_viol_bounds");
}

int dbl_ILLfct_perturb_bounds (dbl_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;

    rval = dbl_expand_var_bounds (lp, lp->tol->ip_tol, &chgb);
#if dbl_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d bounds\n", chgb);
#endif
    ILL_RETURN (rval, "dbl_ILLfct_perturb_bounds");
}

int dbl_ILLfct_perturb_phaseI_bounds (dbl_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;

    rval = dbl_expand_phaseI_bounds (lp, &chgb);
#if dbl_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d phase I bounds\n", chgb);
#endif
    ILL_RETURN (rval, "dbl_ILLfct_perturb_phaseI_bounds");
}

int dbl_ILLfct_bound_shift (dbl_lpinfo * lp,
      int col,
      int bndtype,
      double newbnd)
{
    int rval = 0;
    dbl_bndinfo *nbnd = 0;
    ILL_IFTRACE ("\n%s:%d:%d:%la", __func__, col, bndtype, dbl_EGlpNumToLf (newbnd));
    nbnd = dbl_ILLfct_new_bndinfo ();

    nbnd->varnum = col;
    nbnd->btype = bndtype;
    if (bndtype == BOUND_LOWER) {
	dbl_EGlpNumCopy (nbnd->pbound, lp->lz[col]);
	dbl_EGlpNumCopy (nbnd->cbound, newbnd);
	dbl_EGlpNumCopy (lp->lz[col], newbnd);
    } else {
	dbl_EGlpNumCopy (nbnd->pbound, lp->uz[col]);
	dbl_EGlpNumCopy (nbnd->cbound, newbnd);
	dbl_EGlpNumCopy (lp->uz[col], newbnd);
    }
    ILL_IFTRACE (":%la", dbl_EGlpNumToLf (nbnd->pbound));
    if (lp->vtype[col] == VFIXED || lp->vtype[col] == VARTIFICIAL) {
	/* printf ("changing f/a bound\n"); */
	if (dbl_EGlpNumIsLess (lp->lz[col], lp->uz[col]))
	    lp->vtype[col] = VBOUNDED;
    }
    nbnd->next = lp->bchanges;
    lp->bchanges = nbnd;
    lp->nbchange++;

    /* CLEANUP: */
    if (rval)
	dbl_ILLfct_free_bndinfo (nbnd);
    ILL_IFTRACE ("\n");
    ILL_RETURN (rval, "dbl_ILLfct_bound_shift");
}

void dbl_ILLfct_unroll_bound_change (dbl_lpinfo * lp)
{
    int col;
    int changex = 0;
    dbl_bndinfo *bptr = lp->bchanges;
    dbl_bndinfo *nptr = 0;
    ILL_IFTRACE ("%s:", __func__);

    while (lp->nbchange != 0) {
	col = bptr->varnum;
	ILL_IFTRACE (":%d", col);

	if (bptr->btype == BOUND_UPPER)
	    dbl_EGlpNumCopy (lp->uz[col], bptr->pbound);
	else
	    dbl_EGlpNumCopy (lp->lz[col], bptr->pbound);

	if (lp->vtype[col] == VBOUNDED) {
	    if (dbl_EGlpNumIsEqqual (lp->lz[col], lp->uz[col]))
		lp->vtype[col] = (dbl_EGlpNumIsEqqual (lp->lz[col], dbl_zeroLpNum)) ?
		    VARTIFICIAL : VFIXED;
	}
	if (lp->vstat[col] != STAT_BASIC) {
	    if ((bptr->btype == BOUND_UPPER && lp->vstat[col] == STAT_UPPER) ||
		(bptr->btype == BOUND_LOWER && lp->vstat[col] == STAT_LOWER))
		changex++;
	}
	nptr = bptr->next;
	dbl_EGlpNumClearVar ((bptr->cbound));
	dbl_EGlpNumClearVar ((bptr->pbound));
	ILL_IFFREE (bptr, dbl_bndinfo);
	bptr = nptr;
	lp->nbchange--;
    }
    lp->bchanges = bptr;
    ILL_IFTRACE ("\n");
    if (changex)
	dbl_ILLfct_compute_xbz (lp);
}

static int dbl_expand_var_coefs (dbl_lpinfo * lp,
      double ftol,
      int *chgc)
{
    int rval = 0;
    int i, col, vs, vt;
    int nchg = 0;
    double newc, cftol, mftol[1];
    double *c, *dj;
    ILLrandstate r;
    dbl_EGlpNumInitVar (newc);
    dbl_EGlpNumInitVar (cftol);
    dbl_EGlpNumInitVar (mftol[0]);
    dbl_EGlpNumCopyAbs (cftol, ftol);
    dbl_EGlpNumDivUiTo (cftol, 10);
    dbl_EGlpNumCopyNeg (mftol[0], ftol);
    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nnbasic; i++) {
	dj = &(lp->dz[i]);
	col = lp->nbaz[i];
	c = &(lp->cz[col]);
	vs = lp->vstat[col];
	vt = lp->vtype[col];

	if (vt == VARTIFICIAL || vt == VFIXED)
	    continue;
	switch (vs) {
	case STAT_ZERO:
	    dbl_EGlpNumCopyDiff (newc, *c, *dj);
	    rval = dbl_ILLfct_coef_shift (lp, col, newc);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	    break;
	case STAT_LOWER:
	    if (dbl_EGlpNumIsLess (*dj, ftol)) {
		dbl_EGlpNumSet (newc, dbl_my_rand (50, &(lp->rstate)) + 1.0);
		dbl_EGlpNumMultTo (newc, cftol);
		dbl_EGlpNumAddTo (newc, *c);
		if (dbl_EGlpNumIsLess (*dj, dbl_zeroLpNum))
		    dbl_EGlpNumSubTo (newc, *dj);
		rval = dbl_ILLfct_coef_shift (lp, col, newc);
		ILL_CLEANUP_IF (rval);
		nchg++;
	    }
	    break;
	case STAT_UPPER:
	    if (dbl_EGlpNumIsLess (mftol[0], *dj)) {
		dbl_EGlpNumSet (newc, dbl_my_rand (50, &(lp->rstate)) + 1.0);
		dbl_EGlpNumMultTo (newc, cftol);
		dbl_EGlpNumSign (newc);
		dbl_EGlpNumAddTo (newc, *c);
		if (dbl_EGlpNumIsLess (dbl_zeroLpNum, *dj))
		    dbl_EGlpNumSubTo (newc, *dj);
		rval = dbl_ILLfct_coef_shift (lp, col, newc);
		ILL_CLEANUP_IF (rval);
		nchg++;
	    }
	    break;
	default:
	    break;
	}
    }
    *chgc = nchg;

CLEANUP:
    dbl_EGlpNumClearVar (mftol[0]);
    dbl_EGlpNumClearVar (newc);
    dbl_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "dbl_expand_var_coefs");
}

int dbl_ILLfct_adjust_viol_coefs (dbl_lpinfo * lp)
{
    int rval = 0;
    int chgc = 0;
    double tol;
    dbl_EGlpNumInitVar (tol);
    dbl_EGlpNumCopyNeg (tol, lp->tol->dfeas_tol);

    rval = dbl_expand_var_coefs (lp, tol, &chgc);
#if dbl_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d coefs\n", chgc);
#endif
    dbl_EGlpNumClearVar (tol);
    ILL_RETURN (rval, "dbl_ILLfct_adjust_viol_coefs");
}

int dbl_ILLfct_perturb_coefs (dbl_lpinfo * lp)
{
    int rval = 0;
    int chgc = 0;

    rval = dbl_expand_var_coefs (lp, lp->tol->id_tol, &chgc);
#if dbl_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d coefs\n", chgc);
#endif
    ILL_RETURN (rval, "dbl_ILLfct_perturb_coefs");
}

int dbl_ILLfct_coef_shift (dbl_lpinfo * lp,
      int col,
      double newcoef)
{
    int rval = 0;
    dbl_coefinfo *ncoef = 0;

    ILL_SAFE_MALLOC (ncoef, 1, dbl_coefinfo);
    dbl_EGlpNumInitVar ((ncoef->pcoef));
    dbl_EGlpNumInitVar ((ncoef->ccoef));

    ncoef->varnum = col;
    dbl_EGlpNumCopy (ncoef->pcoef, lp->cz[col]);
    dbl_EGlpNumCopy (ncoef->ccoef, newcoef);
    dbl_EGlpNumCopy (lp->cz[col], newcoef);
    ncoef->next = lp->cchanges;
    lp->cchanges = ncoef;
    dbl_EGlpNumAddTo (lp->dz[lp->vindex[col]], ncoef->ccoef);
    dbl_EGlpNumSubTo (lp->dz[lp->vindex[col]], ncoef->pcoef);
    lp->ncchange++;

CLEANUP:
    if (rval) {
	dbl_EGlpNumClearVar ((ncoef->pcoef));
	dbl_EGlpNumClearVar ((ncoef->ccoef));
	ILL_IFFREE (ncoef, dbl_coefinfo);
    }
    ILL_RETURN (rval, "dbl_ILLfct_coef_shift");
}

void dbl_ILLfct_unroll_coef_change (dbl_lpinfo * lp)
{
    int bascoef = 0;
    dbl_coefinfo *cptr = (dbl_coefinfo *) lp->cchanges;
    dbl_coefinfo *nptr = 0;

    while (lp->ncchange != 0) {
	dbl_EGlpNumCopy (lp->cz[cptr->varnum], cptr->pcoef);
	if (lp->vstat[cptr->varnum] != STAT_BASIC) {
	    dbl_EGlpNumAddTo (lp->dz[lp->vindex[cptr->varnum]], cptr->pcoef);
	    dbl_EGlpNumSubTo (lp->dz[lp->vindex[cptr->varnum]], cptr->ccoef);
	} else
	    bascoef++;

	nptr = cptr->next;
	dbl_EGlpNumClearVar ((cptr->pcoef));
	dbl_EGlpNumClearVar ((cptr->ccoef));
	ILL_IFFREE (cptr, dbl_coefinfo);
	cptr = nptr;
	lp->ncchange--;
    }
    lp->cchanges = cptr;
    if (bascoef) {
	dbl_ILLfct_compute_piz (lp);
	dbl_ILLfct_compute_dz (lp);
    }
}

/* feasibility routines */
void dbl_ILLfct_check_pfeasible (dbl_lpinfo * lp,
      dbl_feas_info * fs,
      double ftol)
{
    int i, col;
    double infeas, err1, err2;
    dbl_EGlpNumInitVar (infeas);
    dbl_EGlpNumInitVar (err1);
    dbl_EGlpNumInitVar (err2);
    dbl_EGlpNumZero (infeas);
    fs->pstatus = PRIMAL_FEASIBLE;
    dbl_EGlpNumZero (fs->totinfeas);
    ILL_IFTRACE ("%s:tol %la\n", __func__, dbl_EGlpNumToLf (ftol));

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	dbl_EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
	dbl_EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
	if (dbl_EGlpNumIsLess (ftol, err1)
	    && dbl_EGlpNumIsNeq (lp->uz[col], dbl_INFTY, dbl_oneLpNum)) {
	    dbl_EGlpNumAddTo (infeas, err1);
	    WARNING (dbl_EGlpNumIsLess (dbl_INFTY, err1),
		"This is imposible lu = %15lg xbz = %15lg" " dbl_INFTY = %15lg",
		dbl_EGlpNumToLf (lp->uz[col]), dbl_EGlpNumToLf (lp->xbz[i]),
		dbl_EGlpNumToLf (dbl_INFTY));
	    lp->bfeas[i] = 1;
	} else if (dbl_EGlpNumIsLess (ftol, err2)
	    && dbl_EGlpNumIsNeq (lp->lz[col], dbl_NINFTY, dbl_oneLpNum)) {
	    dbl_EGlpNumAddTo (infeas, err2);
	    WARNING (dbl_EGlpNumIsLess (dbl_INFTY, err2),
		"This is imposible lz = %15lg xbz = %15lg" " dbl_NINFTY = %15lg",
		dbl_EGlpNumToLf (lp->lz[col]), dbl_EGlpNumToLf (lp->xbz[i]),
		dbl_EGlpNumToLf (dbl_NINFTY));
	    lp->bfeas[i] = -1;
	} else
	    lp->bfeas[i] = 0;
    }
    if (dbl_EGlpNumIsNeqqZero (infeas)) {
	fs->pstatus = PRIMAL_INFEASIBLE;
	dbl_EGlpNumCopy (fs->totinfeas, infeas);
	ILL_IFTRACE ("%s:inf %la\n", __func__, dbl_EGlpNumToLf (infeas));
	if (dbl_EGlpNumIsLess (fs->totinfeas, dbl_zeroLpNum)) {
	    printf ("Negative infeasibility, Imposible! %lf %la\n",
		dbl_EGlpNumToLf (infeas), dbl_EGlpNumToLf (infeas));
	}
    }
    dbl_EGlpNumCopy (lp->pinfeas, infeas);
    dbl_EGlpNumClearVar (infeas);
    dbl_EGlpNumClearVar (err1);
    dbl_EGlpNumClearVar (err2);
}

/* feasibility routines */
void dbl_ILLfct_check_pIpfeasible (dbl_lpinfo * lp,
      dbl_feas_info * fs,
      double ftol)
{
    int i, col;
    int ninf = 0;

    fs->pstatus = PRIMAL_FEASIBLE;
    dbl_EGlpNumZero (fs->totinfeas);

    for (i = 0; i < lp->nrows; i++) {
	if (dbl_EGlpNumIsEqual (lp->xbz[i], dbl_zeroLpNum, ftol))
	    continue;
	col = lp->baz[i];
	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, lp->xbz[i]) &&
	    dbl_EGlpNumIsNeqq (lp->uz[col], dbl_INFTY)) {
	    ninf++;
	} else if (dbl_EGlpNumIsLess (lp->xbz[i], dbl_zeroLpNum) &&
	    dbl_EGlpNumIsNeqq (lp->lz[col], dbl_NINFTY)) {
	    ninf++;
	}
    }
    if (ninf != 0)
	fs->pstatus = PRIMAL_INFEASIBLE;
}

void dbl_ILLfct_check_dfeasible (dbl_lpinfo * lp,
      dbl_feas_info * fs,
      double ftol)
{
    int j, col;
    double infeas;
    dbl_EGlpNumInitVar (infeas);
    dbl_EGlpNumZero (infeas);
    fs->dstatus = DUAL_FEASIBLE;
    dbl_EGlpNumZero (fs->totinfeas);

    for (j = 0; j < lp->nnbasic; j++) {
	lp->dfeas[j] = 0;
	if (dbl_EGlpNumIsEqual (lp->dz[j], dbl_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	if (dbl_EGlpNumIsLess (lp->dz[j], dbl_zeroLpNum) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)) {
	    dbl_EGlpNumSubTo (infeas, lp->dz[j]);
	    lp->dfeas[j] = -1;
	} else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, lp->dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)) {
	    dbl_EGlpNumAddTo (infeas, lp->dz[j]);
	    lp->dfeas[j] = 1;
	}
    }

    if (dbl_EGlpNumIsNeqqZero (infeas)) {
	dbl_EGlpNumCopy (fs->totinfeas, infeas);
	fs->dstatus = DUAL_INFEASIBLE;
	ILL_IFTRACE ("%s:inf %la\n", __func__, dbl_EGlpNumToLf (infeas));
	if (dbl_EGlpNumIsLess (fs->totinfeas, dbl_zeroLpNum)) {
	    printf ("Negative infeasibility, Imposible! %lf %la\n",
		dbl_EGlpNumToLf (infeas), dbl_EGlpNumToLf (infeas));
	}
    }
    dbl_EGlpNumCopy (lp->dinfeas, infeas);
    dbl_EGlpNumClearVar (infeas);
}

void dbl_ILLfct_check_pIdfeasible (dbl_lpinfo * lp,
      dbl_feas_info * fs,
      double ftol)
{
    int j, col;
    int ninf = 0;
    double *dz = lp->pIdz;

    fs->dstatus = DUAL_FEASIBLE;

    for (j = 0; j < lp->nnbasic; j++) {
	if (dbl_EGlpNumIsEqual (dz[j], dbl_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	if (dbl_EGlpNumIsLess (dz[j], dbl_zeroLpNum) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
	    ninf++;
	else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
	    ninf++;
    }

    if (ninf != 0)
	fs->dstatus = DUAL_INFEASIBLE;
}

void dbl_ILLfct_dual_adjust (dbl_lpinfo * lp,
      double ftol)
{
    int j, col;

    for (j = 0; j < lp->nnbasic; j++) {
	if (dbl_EGlpNumIsEqual (lp->dz[j], dbl_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (dbl_EGlpNumIsLess (lp->dz[j], dbl_zeroLpNum) &&
	    dbl_EGlpNumIsNeqq (lp->uz[col], dbl_INFTY))
	    lp->vstat[col] = STAT_UPPER;
	else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, lp->dz[j]) &&
	    dbl_EGlpNumIsNeqq (lp->lz[col], dbl_NINFTY))
	    lp->vstat[col] = STAT_LOWER;
    }
}

void dbl_ILLfct_dphaseI_simple_update (dbl_lpinfo * lp,
      double ftol)
{
    int j, col;

    for (j = 0; j < lp->nnbasic; j++) {
	if (dbl_EGlpNumIsEqual (lp->dz[j], dbl_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (dbl_EGlpNumIsLess (lp->dz[j], dbl_zeroLpNum) && lp->vtype[col] == VBOUNDED)
	    lp->vstat[col] = STAT_UPPER;
	else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, lp->dz[j]) && lp->vtype[col] == VBOUNDED)
	    lp->vstat[col] = STAT_LOWER;
    }
}

/* set status values */
void dbl_ILLfct_set_status_values (dbl_lpinfo * lp,
      int pstatus,
      int dstatus,
      int ptype,
      int dtype)
{
    if (dstatus == DUAL_FEASIBLE && dtype == PHASEII) {
	if (!lp->ncchange) {
	    lp->probstat.dual_feasible = 1;
	    lp->basisstat.dual_feasible = 1;
	    lp->basisstat.dual_infeasible = 0;
	}
    }
    if (dstatus == DUAL_INFEASIBLE && dtype == PHASEII) {
	if (!lp->ncchange) {
	    lp->basisstat.dual_feasible = 0;
	    lp->basisstat.dual_infeasible = 1;
	}
	if (pstatus == PRIMAL_FEASIBLE && ptype == PHASEI)
	    if (!lp->ncchange)
		lp->probstat.dual_infeasible = 1;
    }
    if (pstatus == PRIMAL_FEASIBLE && ptype == PHASEII) {
	if (!lp->nbchange) {
	    lp->probstat.primal_feasible = 1;
	    lp->basisstat.primal_feasible = 1;
	    lp->basisstat.primal_infeasible = 0;
	}
    }
    if (pstatus == PRIMAL_INFEASIBLE && ptype == PHASEII) {
	lp->basisstat.primal_feasible = 0;
	lp->basisstat.primal_infeasible = 1;

	if (dstatus == DUAL_FEASIBLE && dtype == PHASEI)
	    lp->probstat.primal_infeasible = 1;
    }
    if (pstatus == PRIMAL_UNBOUNDED) {
	if (!lp->nbchange) {
	    lp->probstat.primal_unbounded = 1;
	    lp->basisstat.primal_unbounded = 1;
	    lp->probstat.dual_infeasible = 1;
	    lp->basisstat.dual_infeasible = 1;
	    lp->basisstat.dual_feasible = 0;
	}
    }
    if (dstatus == DUAL_UNBOUNDED) {
	if (!lp->ncchange) {
	    lp->probstat.dual_unbounded = 1;
	    lp->basisstat.dual_unbounded = 1;
	    lp->probstat.primal_infeasible = 1;
	    lp->basisstat.primal_infeasible = 1;
	    lp->basisstat.primal_feasible = 0;
	}
    }
    if (lp->probstat.primal_feasible && lp->probstat.dual_feasible)
	lp->probstat.optimal = 1;

    if (lp->basisstat.primal_feasible && lp->basisstat.dual_feasible)
	lp->basisstat.optimal = 1;
    else
	lp->basisstat.optimal = 0;
}

void dbl_ILLfct_init_counts (dbl_lpinfo * lp)
{
    int i;
    dbl_count_struct *c = lp->cnts;
#define dbl_C_VALUE(a) (1.0+(double)(a)/(PARAM_HEAP_RATIO*dbl_ILLutil_our_log2(a)))
    dbl_EGlpNumSet (c->y_ravg, dbl_C_VALUE (lp->nrows));
    dbl_EGlpNumSet (c->za_ravg, dbl_C_VALUE (lp->nnbasic));
    ILL_IFTRACE ("%s:%la\n", __func__, dbl_EGlpNumToLf (c->za_ravg));
#undef dbl_C_VALUE
    c->ynz_cnt = 0;
    c->num_y = 0;
    c->znz_cnt = 0;
    c->num_z = 0;
    c->zanz_cnt = 0;
    c->num_za = 0;
    c->pnorm_cnt = 0;
    c->dnorm_cnt = 0;
    c->pinz_cnt = 0;
    c->num_pi = 0;
    c->pi1nz_cnt = 0;
    c->num_pi1 = 0;
    c->upnz_cnt = 0;
    c->num_up = 0;
    c->pupv_cnt = 0;
    c->dupv_cnt = 0;
    c->pI_iter = 0;
    c->pII_iter = 0;
    c->dI_iter = 0;
    c->dII_iter = 0;
    c->tot_iter = 0;
    for (i = 0; i < 10; i++) {
	c->pivpI[i] = 0;
	c->pivpII[i] = 0;
	c->pivdI[i] = 0;
	c->pivdII[i] = 0;
    }
}

static void dbl_update_piv_values (dbl_count_struct * c,
      int phase,
      double piv2)
{
    int i = 0;
    double v, piv;

    if (dbl_EGlpNumIsEqqual (piv2, dbl_zeroLpNum))
	return;
    dbl_EGlpNumInitVar (v);
    dbl_EGlpNumInitVar (piv);
    dbl_EGlpNumCopyAbs (piv, piv2);
    dbl_EGlpNumOne (v);
    while (dbl_EGlpNumIsLess (piv, v) && i < 9) {
	dbl_EGlpNumDivUiTo (v, 10);
	i++;
    }
    switch (phase) {
    case PRIMAL_PHASEI:
	c->pivpI[i]++;
	break;
    case PRIMAL_PHASEII:
	c->pivpII[i]++;
	break;
    case DUAL_PHASEI:
	c->pivdI[i]++;
	break;
    case DUAL_PHASEII:
	c->pivdII[i]++;
	break;
    default:
	break;
    }
    dbl_EGlpNumClearVar (v);
    dbl_EGlpNumClearVar (piv);
}

void dbl_ILLfct_update_counts (dbl_lpinfo * lp,
      int f,
      int upi,
      double upd)
{
    dbl_count_struct *c = lp->cnts;

    switch (f) {
    case CNT_PPHASE1ITER:
	c->pI_iter++;
	c->tot_iter++;
	break;
    case CNT_PPHASE2ITER:
	c->pII_iter++;
	c->tot_iter++;
	break;
    case CNT_DPHASE1ITER:
	c->dI_iter++;
	c->tot_iter++;
	break;
    case CNT_DPHASE2ITER:
	c->dII_iter++;
	c->tot_iter++;
	break;
    case CNT_YNZ:
	c->ynz_cnt += upi;
	c->num_y++;
	break;
    case CNT_ZANZ:
	c->zanz_cnt += upi;
	c->num_za++;
	break;
    case CNT_PINZ:
	c->pinz_cnt += upi;
	c->num_pi++;
	break;
    case CNT_P1PINZ:
	c->pi1nz_cnt += upi;
	c->num_pi1++;
	break;
    case CNT_UPNZ:
	c->upnz_cnt += upi;
	c->num_up++;
	break;
    case CNT_PIPIV:
	dbl_update_piv_values (c, PRIMAL_PHASEI, upd);
	break;
    case CNT_PIIPIV:
	dbl_update_piv_values (c, PRIMAL_PHASEII, upd);
	break;
    case CNT_DIPIV:
	dbl_update_piv_values (c, DUAL_PHASEI, upd);
	break;
    case CNT_DIIPIV:
	dbl_update_piv_values (c, DUAL_PHASEII, upd);
	break;
    case CNT_YRAVG:
	dbl_EGlpNumMultUiTo (c->y_ravg, c->tot_iter);
	dbl_EGlpNumAddUiTo (c->y_ravg, upi);
	dbl_EGlpNumDivUiTo (c->y_ravg, c->tot_iter + 1);
	break;
    case CNT_ZARAVG:
	ILL_IFTRACE ("%s:%d:%d:%d:%la:%la", __func__, f, c->tot_iter, upi,
	    dbl_EGlpNumToLf (upd), dbl_EGlpNumToLf (c->za_ravg));
	dbl_EGlpNumMultUiTo (c->za_ravg, c->tot_iter);
	dbl_EGlpNumAddUiTo (c->za_ravg, upi);
	dbl_EGlpNumDivUiTo (c->za_ravg, c->tot_iter + 1);
	ILL_IFTRACE (":%la\n", dbl_EGlpNumToLf (c->za_ravg));
	break;
    }
}

void dbl_ILLfct_print_counts (dbl_lpinfo * lp)
{
    int i, niter;
    dbl_count_struct *c = lp->cnts;

    c->tot_iter = c->pI_iter + c->pII_iter + c->dI_iter + c->dII_iter;
    niter = (c->tot_iter == 0) ? 1 : c->tot_iter;
    printf ("Counts for problem %s\n", lp->O->probname);
    if (c->num_y != 0)
	printf ("avg ynz = %.2f\n", (double) c->ynz_cnt / c->num_y);
    if (c->num_z != 0)
	printf ("avg znz = %.2f\n", (double) c->znz_cnt / c->num_z);
    if (c->num_za != 0)
	printf ("avg zanz = %.2f\n", (double) c->zanz_cnt / c->num_za);
    printf ("avg pnorm = %.2f\n", (double) c->pnorm_cnt / lp->nnbasic);
    printf ("avg dnorm = %.2f\n", (double) c->dnorm_cnt / lp->nrows);
    if (c->num_pi != 0)
	printf ("avg pinz = %.2f\n", (double) c->pinz_cnt / c->num_pi);
    if (c->num_pi1 != 0)
	printf ("avg piInz = %.2f\n", (double) c->pi1nz_cnt / c->num_pi1);
    if (c->num_up != 0)
	printf ("avg upnz = %.2f\n", (double) c->upnz_cnt / c->num_up);

    for (i = 0; i < 10; i++)
	printf ("piv 1.0e-%d : %d %d %d %d\n",
	    i, c->pivpI[i], c->pivpII[i], c->pivdI[i], c->pivdII[i]);
}


/* c <- a + t*b */
static void dbl_add_vectors (dbl_lpinfo * lp,
      dbl_svector * a,
      dbl_svector * b,
      dbl_svector * c,
      double t)
{
    int i, r, l;
    dbl_svector *w = &(lp->work);

    for (i = 0; i < b->nzcnt; i++) {
	r = b->indx[i];
	w->indx[i] = r;
	dbl_EGlpNumCopy (w->coef[r], t);
	dbl_EGlpNumMultTo (w->coef[r], b->coef[i]);
	lp->iwork[r] = 1;
    }
    l = b->nzcnt;

    for (i = 0; i < a->nzcnt; i++) {
	r = a->indx[i];
	if (lp->iwork[r] == 0)
	    w->indx[l++] = r;
	dbl_EGlpNumAddTo (w->coef[r], a->coef[i]);
    }
    for (i = 0; i < l; i++) {
	r = w->indx[i];
	c->indx[i] = r;
	dbl_EGlpNumCopy (c->coef[i], w->coef[r]);
	dbl_EGlpNumZero (w->coef[r]);
	lp->iwork[r] = 0;
    }
    w->nzcnt = 0;
    c->nzcnt = l;
}

void dbl_ILLfct_update_pfeas (dbl_lpinfo * lp,
      int lindex,
      dbl_svector * srhs)
{
    int i, k, r;
    int col, nz = 0;
    int cbnd, f;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    double *t = lp->upd.t;
    double tz, *dty, ntmp;
    double *l, *x, *u, *pftol = &(lp->tol->ip_tol);
    dbl_EGlpNumInitVar (tz);
    dbl_EGlpNumInitVar (ntmp);
    dty = &(lp->upd.dty);
    dbl_EGlpNumZero (*dty);
    dbl_EGlpNumCopyAbs (tz, lp->upd.tz);
    dbl_EGlpNumDivUiTo (tz, 100);
    dbl_EGlpNumAddTo (tz, lp->upd.tz);
    ILL_IFTRACE ("%s:%d", __func__, tctr);
    for (i = 0; i < tctr && dbl_EGlpNumIsLeq (t[perm[i]], tz); i++) {
	cbnd = ix[perm[i]] % 10;
	ILL_IFTRACE (":%d", cbnd);
	if (cbnd == BBOUND)
	    continue;
	k = ix[perm[i]] / 10;
	r = lp->yjz.indx[k];
	ILL_IFTRACE (":%d:%d:%d", k, r, lp->iwork[r]);

	if (lp->iwork[r] != 1) {
	    lp->iwork[r] = 1;
	    x = &(lp->xbz[r]);
	    col = lp->baz[r];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);

	    if (r != lindex) {
		f = 0;
		dbl_EGlpNumCopyDiff (ntmp, *l, *x);
		if (dbl_EGlpNumIsNeqq (*l, dbl_NINFTY) && dbl_EGlpNumIsLess (*pftol, ntmp))
		    f = -1;
		else {
		    dbl_EGlpNumCopyDiff (ntmp, *x, *u);
		    if (dbl_EGlpNumIsNeqq (*u, dbl_INFTY) && dbl_EGlpNumIsLess (*pftol, ntmp))
			f = 1;
		}

		ILL_IFTRACE (":%d:%d", f, lp->bfeas[r]);
		if (f != lp->bfeas[r]) {
		    srhs->indx[nz] = r;
		    dbl_EGlpNumSet (srhs->coef[nz], (double) (f - lp->bfeas[r]));
		    dbl_EGlpNumAddInnProdTo (*dty, srhs->coef[nz], lp->yjz.coef[k]);
		    nz++;
		    lp->bfeas[r] = f;
		}
	    } else {
		lp->bfeas[r] = 0;
	    }
	}
    }
    while (--i >= 0) {
	cbnd = ix[perm[i]] % 10;
	if (cbnd == BBOUND)
	    continue;
	k = ix[perm[i]] / 10;
	r = lp->yjz.indx[k];
	lp->iwork[r] = 0;
    }
    srhs->nzcnt = nz;
    ILL_IFTRACE (":%d\n", nz);
    dbl_EGlpNumClearVar (tz);
    dbl_EGlpNumClearVar (ntmp);
}

void dbl_ILLfct_compute_ppIzz (dbl_lpinfo * lp,
      dbl_svector * srhs,
      dbl_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	ILL_IFTRACE ("%s:\n", __func__);
	dbl_ILLbasis_row_solve (lp, srhs, ssoln);
    }
}

void dbl_ILLfct_update_ppI_prices (dbl_lpinfo * lp,
      dbl_price_info * pinf,
      dbl_svector * srhs,
      dbl_svector * ssoln,
      int eindex,
      int lindex,
      double alpha)
{
    double ntmp;
    dbl_EGlpNumInitVar (ntmp);
    dbl_EGlpNumCopy (ntmp, alpha);
    ILL_IFTRACE ("%s:\n", __func__);
    if (lindex == -1) {
	if (srhs->nzcnt != 0) {
	    dbl_ILLfct_update_pIpiz (lp, ssoln, dbl_oneLpNum);
	    if (pinf->p_strategy == COMPLETE_PRICING) {
		dbl_ILLfct_compute_zA (lp, ssoln, &(lp->zA));
		dbl_ILLfct_update_pIdz (lp, &(lp->zA), -1, dbl_oneLpNum);
	    }
	} else {
	    if (pinf->p_strategy == COMPLETE_PRICING)
		dbl_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
	    else
		dbl_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
	    dbl_EGlpNumClearVar (ntmp);
	    return;
	}
    } else {
	if (srhs->nzcnt == 0) {
	    dbl_ILLfct_update_pIpiz (lp, &(lp->zz), ntmp);
	    if (pinf->p_strategy == COMPLETE_PRICING)
		dbl_ILLfct_update_pIdz (lp, &(lp->zA), eindex, ntmp);
	} else {
	    dbl_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
	    dbl_EGlpNumSubTo (ntmp, alpha);
	    dbl_EGlpNumSign (ntmp);
	    dbl_add_vectors (lp, ssoln, &(lp->zz), &(lp->zz), ntmp);
	    dbl_ILLfct_update_pIpiz (lp, &(lp->zz), dbl_oneLpNum);
	    if (pinf->p_strategy == COMPLETE_PRICING) {
		dbl_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
		dbl_ILLfct_update_pIdz (lp, &(lp->zA), eindex, dbl_oneLpNum);
	    }
	}
	dbl_EGlpNumSet (lp->pIdz[eindex], (double) (lp->upd.fs));
	dbl_EGlpNumAddTo (lp->pIdz[eindex], ntmp);
	dbl_EGlpNumSign (lp->pIdz[eindex]);
    }
    if (pinf->p_strategy == COMPLETE_PRICING) {
	dbl_ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
	    PRIMAL_PHASEI);
	if (eindex > -1)
	    dbl_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
	dbl_ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, dbl_zeroLpNum);
    } else
	dbl_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
    dbl_EGlpNumClearVar (ntmp);
    return;
}

void dbl_ILLfct_update_dfeas (dbl_lpinfo * lp,
      int eindex,
      dbl_svector * srhs)
{
    int i, j, k, c;
    int cbnd, col, nz = 0;
    int vs, vt, f;
    int delta;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    int mcnt, mbeg;
    double *t = lp->upd.t;
    double *w = lp->work.coef;
    double tz;
    double *dty = &(lp->upd.dty);
    double *dftol = &(lp->tol->id_tol);
    double dj;
    dbl_EGlpNumInitVar (dj);
    dbl_EGlpNumInitVar (tz);
    dbl_EGlpNumZero (*dty);
    dbl_EGlpNumCopy (tz, lp->upd.tz);
    dbl_EGlpNumMultUiTo (tz, 101);
    dbl_EGlpNumDivUiTo (tz, 100);

    for (j = 0; j < tctr && dbl_EGlpNumIsLeq (t[perm[j]], tz); j++) {
	k = ix[perm[j]] / 10;
	c = lp->zA.indx[k];

	if (lp->iwork[c] != 1) {
	    lp->iwork[c] = 1;
	    cbnd = ix[perm[j]] % 10;
	    col = lp->nbaz[c];
	    dbl_EGlpNumCopy (dj, lp->dz[c]);
	    vs = lp->vstat[col];
	    vt = lp->vtype[col];

	    if (cbnd == BSKIP) {
		if (dbl_EGlpNumIsEqual (dj, dbl_zeroLpNum, *dftol));
		else if (dbl_EGlpNumIsLess (dj, dbl_zeroLpNum) && vs == STAT_LOWER)
		    lp->vstat[col] = STAT_UPPER;
		else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, dj) && vs == STAT_UPPER)
		    lp->vstat[col] = STAT_LOWER;
	    } else if (c != eindex) {
		if (dbl_EGlpNumIsEqual (dj, dbl_zeroLpNum, *dftol))
		    f = 0;
		else if (dbl_EGlpNumIsLess (dj, dbl_zeroLpNum) &&
		    (vs == STAT_LOWER || vs == STAT_ZERO))
		    f = -1;
		else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, dj) &&
		    (vs == STAT_UPPER || vs == STAT_ZERO))
		    f = 1;
		else
		    f = 0;

		if (f != lp->dfeas[c]) {
		    delta = f - lp->dfeas[c];
		    mcnt = lp->matcnt[col];
		    mbeg = lp->matbeg[col];
		    dbl_EGlpNumSet (dj, (double) (delta));
		    for (i = 0; i < mcnt; i++)
			dbl_EGlpNumAddInnProdTo (w[lp->matind[mbeg + i]], dj,
			    lp->matval[mbeg + i]);
		    dbl_EGlpNumAddInnProdTo (*dty, dj, lp->zA.coef[k]);
		    nz = 1;
		    lp->dfeas[c] = f;
		}
	    } else {
		lp->dfeas[c] = 0;
	    }
	}
    }
    while (--j >= 0) {
	k = ix[perm[j]] / 10;
	c = lp->zA.indx[k];
	lp->iwork[c] = 0;
    }

    if (nz) {
	for (i = 0, nz = 0; i < lp->nrows; i++)
	    if (dbl_EGlpNumIsNeqqZero (w[i])) {
		dbl_EGlpNumCopy (srhs->coef[nz], w[i]);
		srhs->indx[nz] = i;
		nz++;
		dbl_EGlpNumZero (w[i]);
	    }
    }
    srhs->nzcnt = nz;
    dbl_EGlpNumClearVar (dj);
    dbl_EGlpNumClearVar (tz);
}

void dbl_ILLfct_compute_dpIy (dbl_lpinfo * lp,
      dbl_svector * srhs,
      dbl_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	dbl_ILLbasis_column_solve (lp, srhs, ssoln);
    }
}

void dbl_ILLfct_update_dpI_prices (dbl_lpinfo * lp,
      dbl_price_info * pinf,
      dbl_svector * srhs,
      dbl_svector * ssoln,
      int lindex,
      double alpha)
{
    int i;
    double ntmp;
    dbl_EGlpNumInitVar (ntmp);
    dbl_EGlpNumZero (ntmp);

    if (srhs->nzcnt == 0) {
	dbl_ILLfct_update_xz (lp, alpha, -1, -1);
    } else {
	dbl_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
	dbl_EGlpNumAddTo (ntmp, alpha);
	dbl_EGlpNumSign (ntmp);
	dbl_add_vectors (lp, ssoln, &(lp->yjz), &(lp->yjz), ntmp);
	dbl_EGlpNumSign (ntmp);
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    dbl_EGlpNumAddTo (lp->xbz[lp->yjz.indx[i]], lp->yjz.coef[i]);
    }
    dbl_EGlpNumSet (lp->xbz[lindex], ((double) (-lp->upd.fs)));
    dbl_EGlpNumAddTo (lp->xbz[lindex], ntmp);

    if (pinf->d_strategy == COMPLETE_PRICING) {
	dbl_ILLprice_compute_primal_inf (lp, pinf, lp->yjz.indx, lp->yjz.nzcnt,
	    DUAL_PHASEI);
	dbl_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEI);
	dbl_ILLfct_update_counts (lp, CNT_YRAVG, lp->yjz.nzcnt, dbl_zeroLpNum);
    } else
	dbl_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEI, ROW_PRICING);
    dbl_EGlpNumClearVar (ntmp);
}

void dbl_ILLfct_update_dIIfeas (dbl_lpinfo * lp,
      int eindex,
      dbl_svector * srhs)
{
    int j, k;
    int col, indx, vs;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    double *zAj, *l, *u;
    double *dty = &(lp->upd.dty);
    double *t_max = &(lp->upd.tz);
    double *t = lp->upd.t;
    double delta;
    dbl_svector a;
    dbl_EGlpNumInitVar (delta);
    dbl_EGlpNumZero (delta);
    dbl_EGlpNumZero (*dty);

    srhs->nzcnt = 0;
    for (j = 0; j < tctr && dbl_EGlpNumIsLeq (t[perm[j]], *t_max); j++) {
	k = ix[perm[j]];
	indx = lp->zA.indx[k];

	if (indx != eindex) {
	    zAj = &(lp->zA.coef[k]);
	    col = lp->nbaz[indx];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);
	    vs = lp->vstat[col];
	    if (vs == STAT_UPPER)
		dbl_EGlpNumCopyDiff (delta, *l, *u);
	    else
		dbl_EGlpNumCopyDiff (delta, *u, *l);
	    dbl_EGlpNumAddInnProdTo (*dty, delta, *zAj);
	    lp->vstat[col] = (vs == STAT_UPPER) ? STAT_LOWER : STAT_UPPER;

	    a.nzcnt = lp->matcnt[col];
	    a.indx = &(lp->matind[lp->matbeg[col]]);
	    a.coef = &(lp->matval[lp->matbeg[col]]);
	    dbl_add_vectors (lp, srhs, &a, srhs, delta);
	}
    }
    dbl_EGlpNumClearVar (delta);
}

void dbl_ILLfct_compute_dpIIy (dbl_lpinfo * lp,
      dbl_svector * srhs,
      dbl_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	dbl_ILLbasis_column_solve (lp, srhs, ssoln);
    }
}

void dbl_ILLfct_update_dpII_prices (dbl_lpinfo * lp,
      dbl_price_info * pinf, dbl_svector * srhs, dbl_svector * ssoln,
      int lindex, double eval, double alpha)
{
    int i;
    dbl_svector *u;

    if (srhs->nzcnt == 0) {
	dbl_ILLfct_update_xz (lp, alpha, -1, -1);
	u = &(lp->yjz);
    } else {
	if (ssoln->nzcnt != 0)
	    for (i = 0; i < ssoln->nzcnt; i++)
		dbl_EGlpNumSubTo (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
	dbl_ILLfct_update_xz (lp, alpha, -1, -1);
	dbl_add_vectors (lp, ssoln, &(lp->yjz), ssoln, dbl_oneLpNum);
	u = ssoln;
    }
    dbl_EGlpNumCopySum (lp->xbz[lindex], eval, alpha);

    if (pinf->d_strategy == COMPLETE_PRICING) {
	dbl_ILLprice_compute_primal_inf (lp, pinf, u->indx, u->nzcnt, DUAL_PHASEII);
	dbl_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEII);
	dbl_ILLfct_update_counts (lp, CNT_YRAVG, u->nzcnt, dbl_zeroLpNum);
    } else
	dbl_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEII, ROW_PRICING);
}

int dbl_ILLfct_test_pivot (dbl_lpinfo * lp, int indx, int indxtype,
      double piv_val)
{
    int i;
    double pval, ntmp;
    dbl_EGlpNumInitVar (pval);
    dbl_EGlpNumInitVar (ntmp);
    dbl_EGlpNumZero (pval);

    if (indxtype == ROW_PIVOT) {
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    if (lp->yjz.indx[i] == indx) {
		dbl_EGlpNumCopy (pval, lp->yjz.coef[i]);
		break;
	    }
    } else {
	for (i = 0; i < lp->zA.nzcnt; i++)
	    if (lp->zA.indx[i] == indx) {
		dbl_EGlpNumCopy (pval, lp->zA.coef[i]);
		break;
	    }
    }
    dbl_EGlpNumCopyDiff (ntmp, pval, piv_val);
    dbl_EGlpNumDivTo (ntmp, piv_val);
    if (dbl_EGlpNumIsLess (ntmp, dbl_zeroLpNum))
	dbl_EGlpNumSign (ntmp);
    if (dbl_EGlpNumIsLess (dbl_ALTPIV_TOLER, ntmp)) {
#if dbl_FCT_DEBUG > 1
	if (indxtype == ROW_PIVOT)
	    printf ("y_i = %.8f, z_j = %.8f %la %la\n", dbl_EGlpNumToLf (pval),
		dbl_EGlpNumToLf (piv_val), dbl_EGlpNumToLf (dbl_ALTPIV_TOLER),
		dbl_EGlpNumToLf (ntmp));
	else
	    printf ("z_j = %.8f, y_i = %.8f\n", dbl_EGlpNumToLf (pval),
		dbl_EGlpNumToLf (piv_val));
#endif
	dbl_EGlpNumClearVar (ntmp);
	dbl_EGlpNumClearVar (pval);
	return 1;
    }
    dbl_EGlpNumClearVar (pval);
    dbl_EGlpNumClearVar (ntmp);
    return 0;
}

#if dbl_FCT_DEBUG > 0

void dbl_fct_test_workvector (dbl_lpinfo * lp)
{
    int i, err = 0;
    for (i = 0; i < lp->ncols; i++) {
	if (dbl_EGlpNumIsNeqqZero (lp->work.coef[i])) {
	    err++;
	    dbl_EGlpNumZero (lp->work.coef[i]);
	}
	if (lp->iwork[i] != 0) {
	    err++;
	    lp->iwork[i] = 0;
	}
    }
    if (err)
	printf ("bad work vector, err=%d\n", err);
}

void dbl_fct_test_pfeasible (dbl_lpinfo * lp)
{
    int i, col;
    int err = 0;
    double *ftol = &(lp->tol->pfeas_tol);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];

	if (dbl_EGlpNumIsNeqq (lp->uz[col], dbl_INFTY)
	    && dbl_EGlpNumIsSumLess (*ftol, lp->uz[col], lp->xbz[i])) {
	    if (lp->bfeas[i] != 1) {
		err++;
		lp->bfeas[i] = 1;
	    }
	} else if (dbl_EGlpNumIsNeqq (lp->lz[col], dbl_NINFTY)
	    && dbl_EGlpNumIsSumLess (lp->xbz[i], *ftol, lp->lz[col])) {
	    if (lp->bfeas[i] != -1) {
		err++;
		lp->bfeas[i] = -1;
	    }
	}
	/* else if (lp->bfeas[i] != 0) {err++; lp->bfeas[i] = 0;} */
    }
    if (err != 0)
	printf ("test_pfeas err =%d\n", err);
}

void dbl_fct_test_dfeasible (dbl_lpinfo * lp)
{
    int j, col;
    int err = 0;
    double *ftol = &(lp->tol->dfeas_tol);
    double mftol[1];
    dbl_EGlpNumInitVar (mftol[0]);
    dbl_EGlpNumCopyNeg (mftol[0], *ftol);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;
	if (dbl_EGlpNumIsLess (lp->dz[j], mftol[0]) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)) {
	    if (lp->dfeas[j] != -1) {
		err++;
		lp->dfeas[j] = -1;
	    }
	}
	if (dbl_EGlpNumIsLess (*ftol, lp->dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)) {
	    if (lp->dfeas[j] != 1) {
		err++;
		lp->dfeas[j] = 1;
	    }
	}
	/* else if (lp->dfeas[j] != 0) {err++; lp->dfeas[j] = 0;} */
    }
    if (err != 0)
	printf ("test_dfeas err =%d\n", err);
}

void dbl_fct_test_pI_x (dbl_lpinfo * lp,
      dbl_price_info * p)
{
    int i;
    int ern = 0;
    double *x;
    double err, diff;
    dbl_EGlpNumInitVar (err);
    dbl_EGlpNumInitVar (diff);
    dbl_EGlpNumZero (err);
    x = dbl_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumCopy (x[i], lp->xbz[i]);
    dbl_ILLfct_compute_phaseI_xbz (lp);
    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_PFEAS_TOLER, diff)) {
	    dbl_EGlpNumAddTo (err, diff);
	    ern++;
	    printf ("bad i = %d\n", i);
	}
    }
    if (dbl_EGlpNumIsNeqqZero (err))
	printf ("dI x err = %.7f, ern = %d\n", dbl_EGlpNumToLf (err), ern);
    dbl_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
    dbl_EGlpNumFreeArray (x);
    dbl_EGlpNumClearVar (diff);
    dbl_EGlpNumClearVar (err);
}

void dbl_fct_test_pII_x (dbl_lpinfo * lp,
      dbl_price_info * p)
{
    int i;
    int ern = 0;
    double *x;
    double err, diff;
    dbl_EGlpNumInitVar (err);
    dbl_EGlpNumInitVar (diff);
    dbl_EGlpNumZero (err);
    x = dbl_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumCopy (x[i], lp->xbz[i]);
    dbl_ILLfct_compute_xbz (lp);
    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_PFEAS_TOLER, diff)) {
	    dbl_EGlpNumAddTo (err, diff);
	    ern++;
	    printf ("bad i = %d\n", i);
	}
    }
    if (dbl_EGlpNumIsNeqqZero (err))
	printf ("dII x err = %.7f, ern = %d\n", dbl_EGlpNumToLf (err), ern);
    dbl_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
    dbl_EGlpNumFreeArray (x);
    dbl_EGlpNumClearVar (diff);
    dbl_EGlpNumClearVar (err);
}

void dbl_fct_test_pI_pi_dz (dbl_lpinfo * lp,
      dbl_price_info * p)
{
    int i;
    int ern = 0;
    double *pidz;
    double err, diff;
    dbl_EGlpNumInitVar (err);
    dbl_EGlpNumInitVar (diff);
    pidz = dbl_EGlpNumAllocArray (lp->ncols);
    dbl_EGlpNumZero (err);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumCopy (pidz[i], lp->pIpiz[i]);
    dbl_ILLfct_compute_phaseI_piz (lp);
    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumCopyDiff (diff, pidz[i], lp->pIpiz[i]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_DFEAS_TOLER, diff)) {
	    dbl_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (dbl_EGlpNumIsNeqqZero (err))
	printf ("pI pi err = %.7f, ern = %d\n", dbl_EGlpNumToLf (err), ern);

    dbl_EGlpNumZero (err);
    ern = 0;
    for (i = 0; i < lp->nnbasic; i++)
	dbl_EGlpNumCopy (pidz[i], lp->pIdz[i]);
    dbl_ILLfct_compute_phaseI_dz (lp);
    for (i = 0; i < lp->nnbasic; i++) {
	dbl_EGlpNumCopyDiff (diff, pidz[i], lp->pIdz[i]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_DFEAS_TOLER, diff)) {
	    dbl_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (dbl_EGlpNumIsNeqqZero (err))
	printf ("pI dz err = %.7f, ern = %d\n", dbl_EGlpNumToLf (err), ern);
    dbl_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
    dbl_EGlpNumClearVar (err);
    dbl_EGlpNumClearVar (diff);
    dbl_EGlpNumFreeArray (pidz);
}

void dbl_fct_test_pII_pi_dz (dbl_lpinfo * lp,
      dbl_price_info * p)
{
    int i;
    int ern = 0;
    double *pidz;
    double err, diff;
    dbl_EGlpNumInitVar (err);
    dbl_EGlpNumInitVar (diff);
    dbl_EGlpNumZero (err);
    pidz = dbl_EGlpNumAllocArray (lp->ncols);

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumCopy (pidz[i], lp->piz[i]);
    dbl_ILLfct_compute_piz (lp);
    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumCopyDiff (diff, pidz[i], lp->piz[i]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_DFEAS_TOLER, diff)) {
	    dbl_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (dbl_EGlpNumIsNeqqZero (err))
	printf ("pII pi err = %.7f, ern = %d\n", dbl_EGlpNumToLf (err), ern);

    dbl_EGlpNumZero (err);
    ern = 0;
    for (i = 0; i < lp->nnbasic; i++)
	dbl_EGlpNumCopy (pidz[i], lp->dz[i]);
    dbl_ILLfct_compute_dz (lp);
    for (i = 0; i < lp->nnbasic; i++) {
	dbl_EGlpNumCopyDiff (diff, pidz[i], lp->dz[i]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_DFEAS_TOLER, diff)) {
	    dbl_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (dbl_EGlpNumIsNeqqZero (err))
	printf ("pII dz err = %.7f, ern = %d\n", dbl_EGlpNumToLf (err), ern);
    /*
     * dbl_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
     */
    dbl_EGlpNumClearVar (err);
    dbl_EGlpNumClearVar (diff);
    dbl_EGlpNumFreeArray (pidz);
}

#endif
