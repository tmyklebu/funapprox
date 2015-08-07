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

/* RCS_INFO = "$RCSfile: mpq_fct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#define mpq_FCT_DEBUG 0

#include "econfig.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdefs.h"
#include "stddefs.h"
#include "mpq_basis.h"
#include "mpq_fct.h"
#include "mpq_price.h"
#include "mpq_ratio.h"
#include "mpq_dstruct.h"

mpq_bndinfo *mpq_ILLfct_new_bndinfo (void)
{
    mpq_bndinfo *nbnd = (mpq_bndinfo *) malloc (sizeof (mpq_bndinfo));
    if (!nbnd) {
	fprintf (stderr, "not enough memory, in %s\n", __func__);
	exit (1);
    }
    mpq_EGlpNumInitVar ((nbnd->pbound));
    mpq_EGlpNumInitVar ((nbnd->cbound));
    return nbnd;
}

void mpq_ILLfct_free_bndinfo (mpq_bndinfo * binfo)
{
    mpq_EGlpNumClearVar ((binfo->pbound));
    mpq_EGlpNumClearVar ((binfo->cbound));
    ILL_IFFREE (binfo, mpq_bndinfo);
    return;
}

static int
  mpq_compute_zA1 (mpq_lpinfo * lp, mpq_svector * z, mpq_svector * zA),
  mpq_compute_zA3 (mpq_lpinfo * lp, mpq_svector * z, mpq_svector * zA),
  mpq_expand_var_bounds (mpq_lpinfo * lp, mpq_t ftol, int *chgb),
  mpq_expand_var_coefs (mpq_lpinfo * lp, mpq_t ftol, int *chgc);

static void mpq_update_piv_values (mpq_count_struct * c,
      int phase,
      mpq_t piv),
/* copy_vectors (mpq_svector * a, mpq_svector * b), */
  mpq_add_vectors (mpq_lpinfo * lp,
      mpq_svector * a,
      mpq_svector * b,
      mpq_svector * c,
      mpq_t t);

static double mpq_my_rand (int bound,
      ILLrandstate * r);


void mpq_ILLfct_load_workvector (mpq_lpinfo * lp,
      mpq_svector * s)
{
    int i;

    for (i = 0; i < s->nzcnt; i++) {
	lp->work.indx[i] = s->indx[i];
	mpq_EGlpNumCopy (lp->work.coef[s->indx[i]], s->coef[i]);
    }
    lp->work.nzcnt = s->nzcnt;
}

void mpq_ILLfct_zero_workvector (mpq_lpinfo * lp)
{
    int i;

    for (i = 0; i < lp->work.nzcnt; i++)
	mpq_EGlpNumZero (lp->work.coef[lp->work.indx[i]]);
    lp->work.nzcnt = 0;
}

void mpq_ILLfct_set_variable_type (mpq_lpinfo * lp)
{
    int j;

    for (j = 0; j < lp->ncols; j++) {

	if (lp->matcnt[j] == 1 && lp->O->rowmap[lp->matind[lp->matbeg[j]]] == j)
	    lp->vclass[j] = CLASS_LOGICAL;
	else
	    lp->vclass[j] = CLASS_STRUCT;
	switch ((mpq_EGlpNumIsEqqual (lp->uz[j], mpq_INFTY) ? 1U : 0U) |
	    (mpq_EGlpNumIsEqqual (lp->lz[j], mpq_NINFTY) ? 2U : 0U)) {
	case 0:
	    if (mpq_EGlpNumIsLess (lp->lz[j], lp->uz[j]))
		lp->vtype[j] = VBOUNDED;
	    else if (mpq_EGlpNumIsEqqual (lp->lz[j], mpq_zeroLpNum) &&
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

void mpq_ILLfct_compute_pobj (mpq_lpinfo * lp)
{
    int i, j;
    int col;
    mpq_t sum;
    mpq_EGlpNumInitVar (sum);
    mpq_EGlpNumZero (sum);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumAddInnProdTo (sum, lp->cz[lp->baz[i]], lp->xbz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	if (lp->vstat[col] == STAT_UPPER)
	    mpq_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER)
	    mpq_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->lz[col]);
    }
    mpq_EGlpNumCopy (lp->pobjval, sum);
    mpq_EGlpNumCopy (lp->objval, sum);
    mpq_EGlpNumClearVar (sum);
}

void mpq_ILLfct_compute_dobj (mpq_lpinfo * lp)
{
    int i, j;
    int col;
    mpq_t sum;
    mpq_EGlpNumInitVar (sum);
    mpq_EGlpNumZero (sum);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumAddInnProdTo (sum, lp->piz[i], lp->bz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	if (lp->vstat[col] == STAT_UPPER)
	    mpq_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER)
	    mpq_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->lz[col]);
    }
    mpq_EGlpNumCopy (lp->dobjval, sum);
    mpq_EGlpNumCopy (lp->objval, sum);
    mpq_EGlpNumClearVar (sum);
}

void mpq_ILLfct_compute_xbz (mpq_lpinfo * lp)
{
    int i, j, r;
    int col, mcnt, mbeg;
    mpq_svector *srhs = &(lp->srhs);
    mpq_svector *ssoln = &(lp->ssoln);
    mpq_t xval;
    mpq_EGlpNumInitVar (xval);

    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumZero (lp->xbz[i]);
	mpq_EGlpNumCopy (srhs->coef[i], lp->bz[i]);
    }
    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	mpq_EGlpNumZero (xval);
	if (lp->vstat[col] == STAT_UPPER && mpq_EGlpNumIsNeqqZero (lp->uz[col]))
	    mpq_EGlpNumCopy (xval, lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER && mpq_EGlpNumIsNeqqZero (lp->lz[col]))
	    mpq_EGlpNumCopy (xval, lp->lz[col]);

	if (mpq_EGlpNumIsNeqqZero (xval)) {
	    mcnt = lp->matcnt[col];
	    mbeg = lp->matbeg[col];
	    for (i = 0; i < mcnt; i++)
		mpq_EGlpNumSubInnProdTo (srhs->coef[lp->matind[mbeg + i]], xval,
		    lp->matval[mbeg + i]);
	}
    }
    for (i = 0, r = 0; i < lp->nrows; i++)
	if (mpq_EGlpNumIsNeqqZero (srhs->coef[i])) {
	    mpq_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
	    srhs->indx[r] = i;
	    r++;
	}
    srhs->nzcnt = r;

    mpq_ILLbasis_column_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpq_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
    mpq_EGlpNumClearVar (xval);
}

void mpq_ILLfct_compute_piz (mpq_lpinfo * lp)
{
    int i, r;
    mpq_svector *srhs = &(lp->srhs);
    mpq_svector *ssoln = &(lp->ssoln);

    for (i = 0, r = 0; i < lp->nrows; i++) {
	mpq_EGlpNumZero (lp->piz[i]);
	if (mpq_EGlpNumIsNeqqZero (lp->cz[lp->baz[i]])) {
	    srhs->indx[r] = i;
	    mpq_EGlpNumCopy (srhs->coef[r], lp->cz[lp->baz[i]]);
	    r++;
	}
    }
    srhs->nzcnt = r;

    mpq_ILLbasis_row_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpq_EGlpNumCopy (lp->piz[ssoln->indx[i]], ssoln->coef[i]);
}

void mpq_ILLfct_compute_dz (mpq_lpinfo * lp)
{
    int i, j;
    int col;
    int mcnt, mbeg;
    mpq_t sum;
    mpq_EGlpNumInitVar (sum);

    for (j = 0; j < lp->nnbasic; j++) {
	mpq_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    mpq_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpq_EGlpNumCopyDiff (lp->dz[j], lp->cz[col], sum);
    }
    mpq_EGlpNumClearVar (sum);
}

void mpq_ILLfct_compute_phaseI_xbz (mpq_lpinfo * lp)
{
    int i, j, r;
    int col, mcnt, mbeg;
    mpq_svector *srhs = &(lp->srhs);
    mpq_svector *ssoln = &(lp->ssoln);

    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumZero (lp->xbz[i]);
	mpq_EGlpNumZero (srhs->coef[i]);
    }
    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];

	if (lp->dfeas[j]) {
	    mcnt = lp->matcnt[col];
	    mbeg = lp->matbeg[col];
	    if (lp->dfeas[j] == -1)
		for (i = 0; i < mcnt; i++)
		    mpq_EGlpNumSubTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
	    else
		for (i = 0; i < mcnt; i++)
		    mpq_EGlpNumAddTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
	}
    }
    for (i = 0, r = 0; i < lp->nrows; i++)
	if (mpq_EGlpNumIsNeqqZero (srhs->coef[i])) {
	    mpq_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
	    srhs->indx[r] = i;
	    r++;
	}
    srhs->nzcnt = r;

    mpq_ILLbasis_column_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpq_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
}

void mpq_ILLfct_compute_phaseI_piz (mpq_lpinfo * lp)
{
    int i, r;
    mpq_svector *srhs = &(lp->srhs);
    mpq_svector *ssoln = &(lp->ssoln);

    for (i = 0, r = 0; i < lp->nrows; i++) {
	mpq_EGlpNumZero (lp->pIpiz[i]);
	if (lp->bfeas[i] != 0) {
	    srhs->indx[r] = i;
	    mpq_EGlpNumSet (srhs->coef[r], (double) lp->bfeas[i]);
	    r++;
	}
    }
    srhs->nzcnt = r;

    mpq_ILLbasis_row_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpq_EGlpNumCopy (lp->pIpiz[ssoln->indx[i]], ssoln->coef[i]);
    mpq_ILLfct_update_counts (lp, CNT_P1PINZ, ssoln->nzcnt, mpq_zeroLpNum);
}

void mpq_ILLfct_compute_phaseI_dz (mpq_lpinfo * lp)
{
    int i, j;
    int col;
    int mcnt, mbeg;
    mpq_t sum;
    mpq_EGlpNumInitVar (sum);
    ILL_IFTRACE ("%s\n", __func__);

    for (j = 0; j < lp->nnbasic; j++) {
	mpq_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    mpq_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpq_EGlpNumCopyNeg (lp->pIdz[j], sum);
	ILL_IFTRACE ("%d:%d:%lf:%la\n", j, col, mpq_EGlpNumToLf (sum),
	    mpq_EGlpNumToLf (sum));
    }
    mpq_EGlpNumClearVar (sum);
}

void mpq_ILLfct_compute_yz (mpq_lpinfo * lp,
      mpq_svector * yz,
      mpq_svector * updz,
      int col)
{
    mpq_svector a;

    a.nzcnt = lp->matcnt[col];
    a.indx = &(lp->matind[lp->matbeg[col]]);
    a.coef = &(lp->matval[lp->matbeg[col]]);

    mpq_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, mpq_PIVZ_TOLER);
    if (updz)
	mpq_ILLbasis_column_solve_update (lp, &a, updz, yz);
    else
	mpq_ILLbasis_column_solve (lp, &a, yz);
    mpq_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, mpq_SZERO_TOLER);
}

void mpq_ILLfct_compute_zz (mpq_lpinfo * lp,
      mpq_svector * zz,
      int row)
{
    mpq_ILLfct_compute_binvrow (lp, zz, row, mpq_PIVZ_TOLER);
}

void mpq_ILLfct_compute_binvrow (mpq_lpinfo * lp,
      mpq_svector * zz,
      int row,
      mpq_t ztoler)
{
    mpq_svector a;
    mpq_t e;
    mpq_EGlpNumInitVar (e);
    mpq_EGlpNumOne (e);

    a.nzcnt = 1;
    a.coef = &e;
    a.indx = &row;

    if (mpq_EGlpNumIsLess (mpq_zeroLpNum, ztoler))
	mpq_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, ztoler);
    mpq_ILLbasis_row_solve (lp, &a, zz);
    if (mpq_EGlpNumIsLess (mpq_zeroLpNum, ztoler))
	mpq_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, mpq_SZERO_TOLER);
    mpq_EGlpNumClearVar (e);
}

void mpq_ILLfct_compute_psteep_upv (mpq_lpinfo * lp,
      mpq_svector * swz)
{
    mpq_ILLbasis_row_solve (lp, &(lp->yjz), swz);
}

void mpq_ILLfct_compute_dsteep_upv (mpq_lpinfo * lp,
      mpq_svector * swz)
{
    mpq_ILLbasis_column_solve (lp, &(lp->zz), swz);
}

static int mpq_compute_zA1 (mpq_lpinfo * lp, mpq_svector * z,
      mpq_svector * zA)
{
    int rval = 0;
    int i, j, nz = 0;
    int col, mcnt, mbeg;
    mpq_t sum;
    mpq_t *v = 0;
    mpq_EGlpNumInitVar (sum);
    v = mpq_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumZero (v[i]);
    for (i = 0; i < z->nzcnt; i++)
	mpq_EGlpNumCopy (v[z->indx[i]], z->coef[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	mpq_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    mpq_EGlpNumAddInnProdTo (sum, v[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

	if (mpq_EGlpNumIsNeqZero (sum, ztoler)) {
	    mpq_EGlpNumCopy (zA->coef[nz], sum);
	    zA->indx[nz] = j;
	    nz++;
	}
    }
    zA->nzcnt = nz;

    mpq_EGlpNumClearVar (sum);
    mpq_EGlpNumFreeArray (v);
    ILL_RETURN (rval, "mpq_compute_zA1");
}


static int mpq_compute_zA3 (mpq_lpinfo * lp, mpq_svector * z, mpq_svector * zA)
{
    int rval = 0;
    int i, j, k, ix;
    int nz = 0;
    int row, col;
    int rcnt, rbeg;
    mpq_t val;
    mpq_EGlpNumInitVar (val);
    k = 0;
    for (i = 0; i < z->nzcnt; i++) {
	row = z->indx[i];
	mpq_EGlpNumCopy (val, z->coef[i]);
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
		mpq_EGlpNumAddInnProdTo (lp->work.coef[ix], val, lp->rowval[rbeg + j]);
	    }
	}
    }
    for (j = 0; j < k; j++) {
	ix = lp->work.indx[j];
	mpq_EGlpNumCopy (val, lp->work.coef[ix]);
	mpq_EGlpNumZero (lp->work.coef[ix]);
	lp->iwork[ix] = 0;
	if (mpq_EGlpNumIsNeqZero (val, ztoler)) {
	    mpq_EGlpNumCopy (zA->coef[nz], val);
	    zA->indx[nz] = ix;
	    nz++;
	}
    }
    zA->nzcnt = nz;
    mpq_EGlpNumClearVar (val);
    ILL_RETURN (rval, "mpq_compute_zA3");
}

int mpq_ILLfct_compute_zA (mpq_lpinfo * lp,
      mpq_svector * z,
      mpq_svector * zA)
{
    if (z->nzcnt < lp->nrows / 2)
	return mpq_compute_zA3 (lp, z, zA);
    else
	return mpq_compute_zA1 (lp, z, zA);
}

/* compute v^T A */
void mpq_ILLfct_compute_vA (mpq_lpinfo * lp,
      mpq_svector * v,
      mpq_t * vA)
{
    int i, j;
    int row, col;
    int rcnt, rbeg;
    mpq_t val;
    mpq_EGlpNumInitVar (val);

    for (j = 0; j < lp->ncols; j++)
	mpq_EGlpNumZero (vA[j]);

    for (i = 0; i < v->nzcnt; i++) {
	row = v->indx[i];
	mpq_EGlpNumCopy (val, v->coef[i]);
	rcnt = lp->rowcnt[row];
	rbeg = lp->rowbeg[row];
	for (j = 0; j < rcnt; j++) {
	    col = lp->rowind[rbeg + j];
	    mpq_EGlpNumAddInnProdTo (vA[col], val, lp->rowval[rbeg + j]);
	}
    }

    for (j = 0; j < lp->ncols; j++)
	if (mpq_EGlpNumIsEqual (vA[j], mpq_zeroLpNum, mpq_SZERO_TOLER))
	    mpq_EGlpNumZero (vA[j]);

    mpq_EGlpNumClearVar (val);
    return;
}

/* update information */

/*
1) lvstat - new status of leaving var.
*/
void mpq_ILLfct_update_basis_info (mpq_lpinfo * lp,
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

void mpq_ILLfct_update_xz (mpq_lpinfo * lp,
      mpq_t tz,
      int eindex,
      int lindex)
{
    int i, evar, estat;
    ILL_IFTRACE ("%s:%la:%d:%d:%d\n", __func__, mpq_EGlpNumToLf (tz), eindex,
	lindex, lp->yjz.nzcnt);

    if (mpq_EGlpNumIsNeqqZero (tz))
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    mpq_EGlpNumSubInnProdTo (lp->xbz[lp->yjz.indx[i]], tz, lp->yjz.coef[i]);

    if (lindex >= 0) {		/* variable leaves basis */
	evar = lp->nbaz[eindex];
	estat = lp->vstat[evar];
	if (estat == STAT_LOWER)
	    mpq_EGlpNumCopySum (lp->xbz[lindex], lp->lz[evar], tz);
	else if (estat == STAT_UPPER)
	    mpq_EGlpNumCopySum (lp->xbz[lindex], lp->uz[evar], tz);
	else if (estat == STAT_ZERO)
	    mpq_EGlpNumCopy (lp->xbz[lindex], tz);
    }
}

void mpq_ILLfct_update_piz (mpq_lpinfo * lp,
      mpq_t alpha)
{
    int i;

    for (i = 0; i < lp->zz.nzcnt; i++)
	mpq_EGlpNumAddInnProdTo (lp->piz[lp->zz.indx[i]], alpha, lp->zz.coef[i]);
}

void mpq_ILLfct_update_pIpiz (mpq_lpinfo * lp,
      mpq_svector * z,
      mpq_t alpha)
{
    int i;
    if (mpq_EGlpNumIsEqqual (alpha, mpq_zeroLpNum))
	return;
    if (mpq_EGlpNumIsEqqual (alpha, mpq_oneLpNum)) {
	for (i = 0; i < z->nzcnt; i++)
	    mpq_EGlpNumAddTo (lp->pIpiz[z->indx[i]], z->coef[i]);
    } else {
	for (i = 0; i < z->nzcnt; i++)
	    mpq_EGlpNumAddInnProdTo (lp->pIpiz[z->indx[i]], alpha, z->coef[i]);
    }
}

void mpq_ILLfct_update_dz (mpq_lpinfo * lp,
      int eindex,
      mpq_t alpha)
{
    int i;

    for (i = 0; i < lp->zA.nzcnt; i++)
	mpq_EGlpNumSubInnProdTo (lp->dz[lp->zA.indx[i]], alpha, lp->zA.coef[i]);
    mpq_EGlpNumCopyNeg (lp->dz[eindex], alpha);
}

void mpq_ILLfct_update_pIdz (mpq_lpinfo * lp,
      mpq_svector * zA,
      int eindex,
      mpq_t alpha)
{
    int i;
    if (mpq_EGlpNumIsEqqual (alpha, mpq_zeroLpNum))
	return;

    if (mpq_EGlpNumIsEqqual (alpha, mpq_oneLpNum)) {
	for (i = 0; i < zA->nzcnt; i++)
	    mpq_EGlpNumSubTo (lp->pIdz[zA->indx[i]], zA->coef[i]);
    } else {
	for (i = 0; i < zA->nzcnt; i++)
	    mpq_EGlpNumSubInnProdTo (lp->pIdz[zA->indx[i]], alpha, zA->coef[i]);
    }
    if (eindex > -1)
	mpq_EGlpNumCopyNeg (lp->pIdz[eindex], alpha);
}

/* bound and coef shift routines */

/* scale bound in mpq_my_rand to get more random digits, unless bound is
   large */
static double mpq_my_rand (int bound,
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

static int mpq_expand_var_bounds (mpq_lpinfo * lp,
      mpq_t ftol,
      int *chgb)
{
    int rval = 0;
    int i, col, nchg = 0;
    mpq_t newb, cftol;
    mpq_t *x, *l, *u;
    ILLrandstate r;
    mpq_EGlpNumInitVar (newb);
    mpq_EGlpNumInitVar (cftol);
    mpq_EGlpNumCopyAbs (cftol, ftol);
    mpq_EGlpNumDivUiTo (cftol, 10);

    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
	    continue;
	x = &(lp->xbz[i]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	/* we use newb as temporal variable outside the if's scope */
	mpq_EGlpNumCopyDiff (newb, *x, ftol);
	if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY) && mpq_EGlpNumIsLess (newb, *l)) {
	    mpq_EGlpNumSet (newb, -1.0 * (mpq_my_rand (50, &(lp->rstate)) + 1.0));
	    mpq_EGlpNumMultTo (newb, cftol);
	    if (mpq_EGlpNumIsLess (*x, *l))
		mpq_EGlpNumAddTo (newb, *x);
	    else
		mpq_EGlpNumAddTo (newb, *l);
	    rval = mpq_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
	mpq_EGlpNumCopySum (newb, *x, ftol);
	if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY) && mpq_EGlpNumIsLess (*u, newb)) {
	    mpq_EGlpNumSet (newb, mpq_my_rand (50, &(lp->rstate)) + 1.0);
	    mpq_EGlpNumMultTo (newb, cftol);
	    if (mpq_EGlpNumIsLess (*x, *u))
		mpq_EGlpNumAddTo (newb, *u);
	    else
		mpq_EGlpNumAddTo (newb, *x);
	    rval = mpq_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
    }
    *chgb = nchg;

CLEANUP:
    mpq_EGlpNumClearVar (newb);
    mpq_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "mpq_expand_var_bounds");
}

static int mpq_expand_phaseI_bounds (mpq_lpinfo * lp,
      int *chgb)
{
    int rval = 0;
    int i, col, nchg = 0;
    mpq_t newb, cftol;
    mpq_t *u, *l, *x;
    ILLrandstate r;
    mpq_EGlpNumInitVar (newb);
    mpq_EGlpNumInitVar (cftol);
    mpq_EGlpNumCopyAbs (cftol, lp->tol->ip_tol);
    mpq_EGlpNumDivUiTo (cftol, 10);
    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
	    continue;
	x = &(lp->xbz[i]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY) && mpq_EGlpNumIsEqual (*x, *l, cftol)) {
	    mpq_EGlpNumSet (newb, mpq_my_rand (50, &(lp->rstate)) + 1.0);
	    mpq_EGlpNumMultTo (newb, cftol);
	    mpq_EGlpNumSign (newb);
	    mpq_EGlpNumAddTo (newb, *l);
	    rval = mpq_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
	if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY) && mpq_EGlpNumIsEqual (*x, *u, cftol)) {
	    mpq_EGlpNumSet (newb, mpq_my_rand (50, &(lp->rstate)) + 1.0);
	    mpq_EGlpNumMultTo (newb, cftol);
	    mpq_EGlpNumAddTo (newb, *u);
	    rval = mpq_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
    }
    *chgb = nchg;

CLEANUP:
    mpq_EGlpNumClearVar (newb);
    mpq_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "mpq_expand_phaseI_bounds");
}

int mpq_ILLfct_adjust_viol_bounds (mpq_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;
    mpq_t tol;
    mpq_EGlpNumInitVar (tol);
    mpq_EGlpNumCopyNeg (tol, lp->tol->pfeas_tol);
    rval = mpq_expand_var_bounds (lp, tol, &chgb);
#if mpq_FCT_DEBUG > 0
    if (rval == 0)
	printf ("adjusting %d bounds\n", chgb);
#endif
    mpq_EGlpNumClearVar (tol);
    ILL_RETURN (rval, "mpq_ILLfct_adjust_viol_bounds");
}

int mpq_ILLfct_perturb_bounds (mpq_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;

    rval = mpq_expand_var_bounds (lp, lp->tol->ip_tol, &chgb);
#if mpq_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d bounds\n", chgb);
#endif
    ILL_RETURN (rval, "mpq_ILLfct_perturb_bounds");
}

int mpq_ILLfct_perturb_phaseI_bounds (mpq_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;

    rval = mpq_expand_phaseI_bounds (lp, &chgb);
#if mpq_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d phase I bounds\n", chgb);
#endif
    ILL_RETURN (rval, "mpq_ILLfct_perturb_phaseI_bounds");
}

int mpq_ILLfct_bound_shift (mpq_lpinfo * lp,
      int col,
      int bndtype,
      mpq_t newbnd)
{
    int rval = 0;
    mpq_bndinfo *nbnd = 0;
    ILL_IFTRACE ("\n%s:%d:%d:%la", __func__, col, bndtype, mpq_EGlpNumToLf (newbnd));
    nbnd = mpq_ILLfct_new_bndinfo ();

    nbnd->varnum = col;
    nbnd->btype = bndtype;
    if (bndtype == BOUND_LOWER) {
	mpq_EGlpNumCopy (nbnd->pbound, lp->lz[col]);
	mpq_EGlpNumCopy (nbnd->cbound, newbnd);
	mpq_EGlpNumCopy (lp->lz[col], newbnd);
    } else {
	mpq_EGlpNumCopy (nbnd->pbound, lp->uz[col]);
	mpq_EGlpNumCopy (nbnd->cbound, newbnd);
	mpq_EGlpNumCopy (lp->uz[col], newbnd);
    }
    ILL_IFTRACE (":%la", mpq_EGlpNumToLf (nbnd->pbound));
    if (lp->vtype[col] == VFIXED || lp->vtype[col] == VARTIFICIAL) {
	/* printf ("changing f/a bound\n"); */
	if (mpq_EGlpNumIsLess (lp->lz[col], lp->uz[col]))
	    lp->vtype[col] = VBOUNDED;
    }
    nbnd->next = lp->bchanges;
    lp->bchanges = nbnd;
    lp->nbchange++;

    /* CLEANUP: */
    if (rval)
	mpq_ILLfct_free_bndinfo (nbnd);
    ILL_IFTRACE ("\n");
    ILL_RETURN (rval, "mpq_ILLfct_bound_shift");
}

void mpq_ILLfct_unroll_bound_change (mpq_lpinfo * lp)
{
    int col;
    int changex = 0;
    mpq_bndinfo *bptr = lp->bchanges;
    mpq_bndinfo *nptr = 0;
    ILL_IFTRACE ("%s:", __func__);

    while (lp->nbchange != 0) {
	col = bptr->varnum;
	ILL_IFTRACE (":%d", col);

	if (bptr->btype == BOUND_UPPER)
	    mpq_EGlpNumCopy (lp->uz[col], bptr->pbound);
	else
	    mpq_EGlpNumCopy (lp->lz[col], bptr->pbound);

	if (lp->vtype[col] == VBOUNDED) {
	    if (mpq_EGlpNumIsEqqual (lp->lz[col], lp->uz[col]))
		lp->vtype[col] = (mpq_EGlpNumIsEqqual (lp->lz[col], mpq_zeroLpNum)) ?
		    VARTIFICIAL : VFIXED;
	}
	if (lp->vstat[col] != STAT_BASIC) {
	    if ((bptr->btype == BOUND_UPPER && lp->vstat[col] == STAT_UPPER) ||
		(bptr->btype == BOUND_LOWER && lp->vstat[col] == STAT_LOWER))
		changex++;
	}
	nptr = bptr->next;
	mpq_EGlpNumClearVar ((bptr->cbound));
	mpq_EGlpNumClearVar ((bptr->pbound));
	ILL_IFFREE (bptr, mpq_bndinfo);
	bptr = nptr;
	lp->nbchange--;
    }
    lp->bchanges = bptr;
    ILL_IFTRACE ("\n");
    if (changex)
	mpq_ILLfct_compute_xbz (lp);
}

static int mpq_expand_var_coefs (mpq_lpinfo * lp,
      mpq_t ftol,
      int *chgc)
{
    int rval = 0;
    int i, col, vs, vt;
    int nchg = 0;
    mpq_t newc, cftol, mftol[1];
    mpq_t *c, *dj;
    ILLrandstate r;
    mpq_EGlpNumInitVar (newc);
    mpq_EGlpNumInitVar (cftol);
    mpq_EGlpNumInitVar (mftol[0]);
    mpq_EGlpNumCopyAbs (cftol, ftol);
    mpq_EGlpNumDivUiTo (cftol, 10);
    mpq_EGlpNumCopyNeg (mftol[0], ftol);
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
	    mpq_EGlpNumCopyDiff (newc, *c, *dj);
	    rval = mpq_ILLfct_coef_shift (lp, col, newc);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	    break;
	case STAT_LOWER:
	    if (mpq_EGlpNumIsLess (*dj, ftol)) {
		mpq_EGlpNumSet (newc, mpq_my_rand (50, &(lp->rstate)) + 1.0);
		mpq_EGlpNumMultTo (newc, cftol);
		mpq_EGlpNumAddTo (newc, *c);
		if (mpq_EGlpNumIsLess (*dj, mpq_zeroLpNum))
		    mpq_EGlpNumSubTo (newc, *dj);
		rval = mpq_ILLfct_coef_shift (lp, col, newc);
		ILL_CLEANUP_IF (rval);
		nchg++;
	    }
	    break;
	case STAT_UPPER:
	    if (mpq_EGlpNumIsLess (mftol[0], *dj)) {
		mpq_EGlpNumSet (newc, mpq_my_rand (50, &(lp->rstate)) + 1.0);
		mpq_EGlpNumMultTo (newc, cftol);
		mpq_EGlpNumSign (newc);
		mpq_EGlpNumAddTo (newc, *c);
		if (mpq_EGlpNumIsLess (mpq_zeroLpNum, *dj))
		    mpq_EGlpNumSubTo (newc, *dj);
		rval = mpq_ILLfct_coef_shift (lp, col, newc);
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
    mpq_EGlpNumClearVar (mftol[0]);
    mpq_EGlpNumClearVar (newc);
    mpq_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "mpq_expand_var_coefs");
}

int mpq_ILLfct_adjust_viol_coefs (mpq_lpinfo * lp)
{
    int rval = 0;
    int chgc = 0;
    mpq_t tol;
    mpq_EGlpNumInitVar (tol);
    mpq_EGlpNumCopyNeg (tol, lp->tol->dfeas_tol);

    rval = mpq_expand_var_coefs (lp, tol, &chgc);
#if mpq_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d coefs\n", chgc);
#endif
    mpq_EGlpNumClearVar (tol);
    ILL_RETURN (rval, "mpq_ILLfct_adjust_viol_coefs");
}

int mpq_ILLfct_perturb_coefs (mpq_lpinfo * lp)
{
    int rval = 0;
    int chgc = 0;

    rval = mpq_expand_var_coefs (lp, lp->tol->id_tol, &chgc);
#if mpq_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d coefs\n", chgc);
#endif
    ILL_RETURN (rval, "mpq_ILLfct_perturb_coefs");
}

int mpq_ILLfct_coef_shift (mpq_lpinfo * lp,
      int col,
      mpq_t newcoef)
{
    int rval = 0;
    mpq_coefinfo *ncoef = 0;

    ILL_SAFE_MALLOC (ncoef, 1, mpq_coefinfo);
    mpq_EGlpNumInitVar ((ncoef->pcoef));
    mpq_EGlpNumInitVar ((ncoef->ccoef));

    ncoef->varnum = col;
    mpq_EGlpNumCopy (ncoef->pcoef, lp->cz[col]);
    mpq_EGlpNumCopy (ncoef->ccoef, newcoef);
    mpq_EGlpNumCopy (lp->cz[col], newcoef);
    ncoef->next = lp->cchanges;
    lp->cchanges = ncoef;
    mpq_EGlpNumAddTo (lp->dz[lp->vindex[col]], ncoef->ccoef);
    mpq_EGlpNumSubTo (lp->dz[lp->vindex[col]], ncoef->pcoef);
    lp->ncchange++;

CLEANUP:
    if (rval) {
	mpq_EGlpNumClearVar ((ncoef->pcoef));
	mpq_EGlpNumClearVar ((ncoef->ccoef));
	ILL_IFFREE (ncoef, mpq_coefinfo);
    }
    ILL_RETURN (rval, "mpq_ILLfct_coef_shift");
}

void mpq_ILLfct_unroll_coef_change (mpq_lpinfo * lp)
{
    int bascoef = 0;
    mpq_coefinfo *cptr = (mpq_coefinfo *) lp->cchanges;
    mpq_coefinfo *nptr = 0;

    while (lp->ncchange != 0) {
	mpq_EGlpNumCopy (lp->cz[cptr->varnum], cptr->pcoef);
	if (lp->vstat[cptr->varnum] != STAT_BASIC) {
	    mpq_EGlpNumAddTo (lp->dz[lp->vindex[cptr->varnum]], cptr->pcoef);
	    mpq_EGlpNumSubTo (lp->dz[lp->vindex[cptr->varnum]], cptr->ccoef);
	} else
	    bascoef++;

	nptr = cptr->next;
	mpq_EGlpNumClearVar ((cptr->pcoef));
	mpq_EGlpNumClearVar ((cptr->ccoef));
	ILL_IFFREE (cptr, mpq_coefinfo);
	cptr = nptr;
	lp->ncchange--;
    }
    lp->cchanges = cptr;
    if (bascoef) {
	mpq_ILLfct_compute_piz (lp);
	mpq_ILLfct_compute_dz (lp);
    }
}

/* feasibility routines */
void mpq_ILLfct_check_pfeasible (mpq_lpinfo * lp,
      mpq_feas_info * fs,
      mpq_t ftol)
{
    int i, col;
    mpq_t infeas, err1, err2;
    mpq_EGlpNumInitVar (infeas);
    mpq_EGlpNumInitVar (err1);
    mpq_EGlpNumInitVar (err2);
    mpq_EGlpNumZero (infeas);
    fs->pstatus = PRIMAL_FEASIBLE;
    mpq_EGlpNumZero (fs->totinfeas);
    ILL_IFTRACE ("%s:tol %la\n", __func__, mpq_EGlpNumToLf (ftol));

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	mpq_EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
	mpq_EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
	if (mpq_EGlpNumIsLess (ftol, err1)
	    && mpq_EGlpNumIsNeq (lp->uz[col], mpq_INFTY, mpq_oneLpNum)) {
	    mpq_EGlpNumAddTo (infeas, err1);
	    WARNING (mpq_EGlpNumIsLess (mpq_INFTY, err1),
		"This is imposible lu = %15lg xbz = %15lg" " mpq_INFTY = %15lg",
		mpq_EGlpNumToLf (lp->uz[col]), mpq_EGlpNumToLf (lp->xbz[i]),
		mpq_EGlpNumToLf (mpq_INFTY));
	    lp->bfeas[i] = 1;
	} else if (mpq_EGlpNumIsLess (ftol, err2)
	    && mpq_EGlpNumIsNeq (lp->lz[col], mpq_NINFTY, mpq_oneLpNum)) {
	    mpq_EGlpNumAddTo (infeas, err2);
	    WARNING (mpq_EGlpNumIsLess (mpq_INFTY, err2),
		"This is imposible lz = %15lg xbz = %15lg" " mpq_NINFTY = %15lg",
		mpq_EGlpNumToLf (lp->lz[col]), mpq_EGlpNumToLf (lp->xbz[i]),
		mpq_EGlpNumToLf (mpq_NINFTY));
	    lp->bfeas[i] = -1;
	} else
	    lp->bfeas[i] = 0;
    }
    if (mpq_EGlpNumIsNeqqZero (infeas)) {
	fs->pstatus = PRIMAL_INFEASIBLE;
	mpq_EGlpNumCopy (fs->totinfeas, infeas);
	ILL_IFTRACE ("%s:inf %la\n", __func__, mpq_EGlpNumToLf (infeas));
	if (mpq_EGlpNumIsLess (fs->totinfeas, mpq_zeroLpNum)) {
	    printf ("Negative infeasibility, Imposible! %lf %la\n",
		mpq_EGlpNumToLf (infeas), mpq_EGlpNumToLf (infeas));
	}
    }
    mpq_EGlpNumCopy (lp->pinfeas, infeas);
    mpq_EGlpNumClearVar (infeas);
    mpq_EGlpNumClearVar (err1);
    mpq_EGlpNumClearVar (err2);
}

/* feasibility routines */
void mpq_ILLfct_check_pIpfeasible (mpq_lpinfo * lp, mpq_feas_info * fs)
{
    int i, col;
    int ninf = 0;

    fs->pstatus = PRIMAL_FEASIBLE;
    mpq_EGlpNumZero (fs->totinfeas);

    for (i = 0; i < lp->nrows; i++) {
	if (mpq_EGlpNumIsEqual (lp->xbz[i], mpq_zeroLpNum, ftol))
	    continue;
	col = lp->baz[i];
	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, lp->xbz[i]) &&
	    mpq_EGlpNumIsNeqq (lp->uz[col], mpq_INFTY)) {
	    ninf++;
	} else if (mpq_EGlpNumIsLess (lp->xbz[i], mpq_zeroLpNum) &&
	    mpq_EGlpNumIsNeqq (lp->lz[col], mpq_NINFTY)) {
	    ninf++;
	}
    }
    if (ninf != 0)
	fs->pstatus = PRIMAL_INFEASIBLE;
}

void mpq_ILLfct_check_dfeasible (mpq_lpinfo * lp, mpq_feas_info * fs)
{
    int j, col;
    mpq_t infeas;
    mpq_EGlpNumInitVar (infeas);
    mpq_EGlpNumZero (infeas);
    fs->dstatus = DUAL_FEASIBLE;
    mpq_EGlpNumZero (fs->totinfeas);

    for (j = 0; j < lp->nnbasic; j++) {
	lp->dfeas[j] = 0;
	if (mpq_EGlpNumIsEqual (lp->dz[j], mpq_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	if (mpq_EGlpNumIsLess (lp->dz[j], mpq_zeroLpNum) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)) {
	    mpq_EGlpNumSubTo (infeas, lp->dz[j]);
	    lp->dfeas[j] = -1;
	} else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, lp->dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)) {
	    mpq_EGlpNumAddTo (infeas, lp->dz[j]);
	    lp->dfeas[j] = 1;
	}
    }

    if (mpq_EGlpNumIsNeqqZero (infeas)) {
	mpq_EGlpNumCopy (fs->totinfeas, infeas);
	fs->dstatus = DUAL_INFEASIBLE;
	ILL_IFTRACE ("%s:inf %la\n", __func__, mpq_EGlpNumToLf (infeas));
	if (mpq_EGlpNumIsLess (fs->totinfeas, mpq_zeroLpNum)) {
	    printf ("Negative infeasibility, Imposible! %lf %la\n",
		mpq_EGlpNumToLf (infeas), mpq_EGlpNumToLf (infeas));
	}
    }
    mpq_EGlpNumCopy (lp->dinfeas, infeas);
    mpq_EGlpNumClearVar (infeas);
}

void mpq_ILLfct_check_pIdfeasible (mpq_lpinfo * lp, mpq_feas_info * fs)
{
    int j, col;
    int ninf = 0;
    mpq_t *dz = lp->pIdz;

    fs->dstatus = DUAL_FEASIBLE;

    for (j = 0; j < lp->nnbasic; j++) {
	if (mpq_EGlpNumIsEqual (dz[j], mpq_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	if (mpq_EGlpNumIsLess (dz[j], mpq_zeroLpNum) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
	    ninf++;
	else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
	    ninf++;
    }

    if (ninf != 0)
	fs->dstatus = DUAL_INFEASIBLE;
}

void mpq_ILLfct_dual_adjust (mpq_lpinfo * lp)
{
    int j, col;

    for (j = 0; j < lp->nnbasic; j++) {
	if (mpq_EGlpNumIsEqual (lp->dz[j], mpq_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (mpq_EGlpNumIsLess (lp->dz[j], mpq_zeroLpNum) &&
	    mpq_EGlpNumIsNeqq (lp->uz[col], mpq_INFTY))
	    lp->vstat[col] = STAT_UPPER;
	else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, lp->dz[j]) &&
	    mpq_EGlpNumIsNeqq (lp->lz[col], mpq_NINFTY))
	    lp->vstat[col] = STAT_LOWER;
    }
}

void mpq_ILLfct_dphaseI_simple_update (mpq_lpinfo * lp)
{
    int j, col;

    for (j = 0; j < lp->nnbasic; j++) {
	if (mpq_EGlpNumIsEqual (lp->dz[j], mpq_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (mpq_EGlpNumIsLess (lp->dz[j], mpq_zeroLpNum) && lp->vtype[col] == VBOUNDED)
	    lp->vstat[col] = STAT_UPPER;
	else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, lp->dz[j]) && lp->vtype[col] == VBOUNDED)
	    lp->vstat[col] = STAT_LOWER;
    }
}

/* set status values */
void mpq_ILLfct_set_status_values (mpq_lpinfo * lp,
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

void mpq_ILLfct_init_counts (mpq_lpinfo * lp)
{
    int i;
    mpq_count_struct *c = lp->cnts;
#define mpq_C_VALUE(a) (1.0+(double)(a)/(PARAM_HEAP_RATIO*mpq_ILLutil_our_log2(a)))
    mpq_EGlpNumSet (c->y_ravg, mpq_C_VALUE (lp->nrows));
    mpq_EGlpNumSet (c->za_ravg, mpq_C_VALUE (lp->nnbasic));
    ILL_IFTRACE ("%s:%la\n", __func__, mpq_EGlpNumToLf (c->za_ravg));
#undef mpq_C_VALUE
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

static void mpq_update_piv_values (mpq_count_struct * c,
      int phase,
      mpq_t piv2)
{
    int i = 0;
    mpq_t v, piv;

    if (mpq_EGlpNumIsEqqual (piv2, mpq_zeroLpNum))
	return;
    mpq_EGlpNumInitVar (v);
    mpq_EGlpNumInitVar (piv);
    mpq_EGlpNumCopyAbs (piv, piv2);
    mpq_EGlpNumOne (v);
    while (mpq_EGlpNumIsLess (piv, v) && i < 9) {
	mpq_EGlpNumDivUiTo (v, 10);
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
    mpq_EGlpNumClearVar (v);
    mpq_EGlpNumClearVar (piv);
}

void mpq_ILLfct_update_counts (mpq_lpinfo * lp,
      int f,
      int upi,
      mpq_t upd)
{
    mpq_count_struct *c = lp->cnts;

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
	mpq_update_piv_values (c, PRIMAL_PHASEI, upd);
	break;
    case CNT_PIIPIV:
	mpq_update_piv_values (c, PRIMAL_PHASEII, upd);
	break;
    case CNT_DIPIV:
	mpq_update_piv_values (c, DUAL_PHASEI, upd);
	break;
    case CNT_DIIPIV:
	mpq_update_piv_values (c, DUAL_PHASEII, upd);
	break;
    case CNT_YRAVG:
	mpq_EGlpNumMultUiTo (c->y_ravg, c->tot_iter);
	mpq_EGlpNumAddUiTo (c->y_ravg, upi);
	mpq_EGlpNumDivUiTo (c->y_ravg, c->tot_iter + 1);
	break;
    case CNT_ZARAVG:
	ILL_IFTRACE ("%s:%d:%d:%d:%la:%la", __func__, f, c->tot_iter, upi,
	    mpq_EGlpNumToLf (upd), mpq_EGlpNumToLf (c->za_ravg));
	mpq_EGlpNumMultUiTo (c->za_ravg, c->tot_iter);
	mpq_EGlpNumAddUiTo (c->za_ravg, upi);
	mpq_EGlpNumDivUiTo (c->za_ravg, c->tot_iter + 1);
	ILL_IFTRACE (":%la\n", mpq_EGlpNumToLf (c->za_ravg));
	break;
    }
}

void mpq_ILLfct_print_counts (mpq_lpinfo * lp)
{
    int i, niter;
    mpq_count_struct *c = lp->cnts;

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
static void mpq_add_vectors (mpq_lpinfo * lp,
      mpq_svector * a,
      mpq_svector * b,
      mpq_svector * c,
      mpq_t t)
{
    int i, r, l;
    mpq_svector *w = &(lp->work);

    for (i = 0; i < b->nzcnt; i++) {
	r = b->indx[i];
	w->indx[i] = r;
	mpq_EGlpNumCopy (w->coef[r], t);
	mpq_EGlpNumMultTo (w->coef[r], b->coef[i]);
	lp->iwork[r] = 1;
    }
    l = b->nzcnt;

    for (i = 0; i < a->nzcnt; i++) {
	r = a->indx[i];
	if (lp->iwork[r] == 0)
	    w->indx[l++] = r;
	mpq_EGlpNumAddTo (w->coef[r], a->coef[i]);
    }
    for (i = 0; i < l; i++) {
	r = w->indx[i];
	c->indx[i] = r;
	mpq_EGlpNumCopy (c->coef[i], w->coef[r]);
	mpq_EGlpNumZero (w->coef[r]);
	lp->iwork[r] = 0;
    }
    w->nzcnt = 0;
    c->nzcnt = l;
}

void mpq_ILLfct_update_pfeas (mpq_lpinfo * lp,
      int lindex,
      mpq_svector * srhs)
{
    int i, k, r;
    int col, nz = 0;
    int cbnd, f;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    mpq_t *t = lp->upd.t;
    mpq_t tz, *dty, ntmp;
    mpq_t *l, *x, *u, *pftol = &(lp->tol->ip_tol);
    mpq_EGlpNumInitVar (tz);
    mpq_EGlpNumInitVar (ntmp);
    dty = &(lp->upd.dty);
    mpq_EGlpNumZero (*dty);
    mpq_EGlpNumCopyAbs (tz, lp->upd.tz);
    mpq_EGlpNumDivUiTo (tz, 100);
    mpq_EGlpNumAddTo (tz, lp->upd.tz);
    ILL_IFTRACE ("%s:%d", __func__, tctr);
    for (i = 0; i < tctr && mpq_EGlpNumIsLeq (t[perm[i]], tz); i++) {
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
		mpq_EGlpNumCopyDiff (ntmp, *l, *x);
		if (mpq_EGlpNumIsNeqq (*l, mpq_NINFTY) && mpq_EGlpNumIsLess (*pftol, ntmp))
		    f = -1;
		else {
		    mpq_EGlpNumCopyDiff (ntmp, *x, *u);
		    if (mpq_EGlpNumIsNeqq (*u, mpq_INFTY) && mpq_EGlpNumIsLess (*pftol, ntmp))
			f = 1;
		}

		ILL_IFTRACE (":%d:%d", f, lp->bfeas[r]);
		if (f != lp->bfeas[r]) {
		    srhs->indx[nz] = r;
		    mpq_EGlpNumSet (srhs->coef[nz], (double) (f - lp->bfeas[r]));
		    mpq_EGlpNumAddInnProdTo (*dty, srhs->coef[nz], lp->yjz.coef[k]);
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
    mpq_EGlpNumClearVar (tz);
    mpq_EGlpNumClearVar (ntmp);
}

void mpq_ILLfct_compute_ppIzz (mpq_lpinfo * lp,
      mpq_svector * srhs,
      mpq_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	ILL_IFTRACE ("%s:\n", __func__);
	mpq_ILLbasis_row_solve (lp, srhs, ssoln);
    }
}

void mpq_ILLfct_update_ppI_prices (mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_svector * srhs,
      mpq_svector * ssoln,
      int eindex,
      int lindex,
      mpq_t alpha)
{
    mpq_t ntmp;
    mpq_EGlpNumInitVar (ntmp);
    mpq_EGlpNumCopy (ntmp, alpha);
    ILL_IFTRACE ("%s:\n", __func__);
    if (lindex == -1) {
	if (srhs->nzcnt != 0) {
	    mpq_ILLfct_update_pIpiz (lp, ssoln, mpq_oneLpNum);
	    if (pinf->p_strategy == COMPLETE_PRICING) {
		mpq_ILLfct_compute_zA (lp, ssoln, &(lp->zA));
		mpq_ILLfct_update_pIdz (lp, &(lp->zA), -1, mpq_oneLpNum);
	    }
	} else {
	    if (pinf->p_strategy == COMPLETE_PRICING)
		mpq_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
	    else
		mpq_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
	    mpq_EGlpNumClearVar (ntmp);
	    return;
	}
    } else {
	if (srhs->nzcnt == 0) {
	    mpq_ILLfct_update_pIpiz (lp, &(lp->zz), ntmp);
	    if (pinf->p_strategy == COMPLETE_PRICING)
		mpq_ILLfct_update_pIdz (lp, &(lp->zA), eindex, ntmp);
	} else {
	    mpq_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
	    mpq_EGlpNumSubTo (ntmp, alpha);
	    mpq_EGlpNumSign (ntmp);
	    mpq_add_vectors (lp, ssoln, &(lp->zz), &(lp->zz), ntmp);
	    mpq_ILLfct_update_pIpiz (lp, &(lp->zz), mpq_oneLpNum);
	    if (pinf->p_strategy == COMPLETE_PRICING) {
		mpq_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
		mpq_ILLfct_update_pIdz (lp, &(lp->zA), eindex, mpq_oneLpNum);
	    }
	}
	mpq_EGlpNumSet (lp->pIdz[eindex], (double) (lp->upd.fs));
	mpq_EGlpNumAddTo (lp->pIdz[eindex], ntmp);
	mpq_EGlpNumSign (lp->pIdz[eindex]);
    }
    if (pinf->p_strategy == COMPLETE_PRICING) {
	mpq_ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
	    PRIMAL_PHASEI);
	if (eindex > -1)
	    mpq_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
	mpq_ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, mpq_zeroLpNum);
    } else
	mpq_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
    mpq_EGlpNumClearVar (ntmp);
    return;
}

void mpq_ILLfct_update_dfeas (mpq_lpinfo * lp, int eindex, mpq_svector * srhs)
{
    int i, j, k, c;
    int cbnd, col, nz = 0;
    int vs, vt, f;
    int delta;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    int mcnt, mbeg;
    mpq_t *t = lp->upd.t;
    mpq_t *w = lp->work.coef;
    mpq_t tz;
    mpq_t *dty = &(lp->upd.dty);
    mpq_t dj;
    mpq_EGlpNumInitVar (dj);
    mpq_EGlpNumInitVar (tz);
    mpq_EGlpNumZero (*dty);
    mpq_EGlpNumCopy (tz, lp->upd.tz);
    mpq_EGlpNumMultUiTo (tz, 101);
    mpq_EGlpNumDivUiTo (tz, 100);

    for (j = 0; j < tctr && mpq_EGlpNumIsLeq (t[perm[j]], tz); j++) {
	k = ix[perm[j]] / 10;
	c = lp->zA.indx[k];

	if (lp->iwork[c] != 1) {
	    lp->iwork[c] = 1;
	    cbnd = ix[perm[j]] % 10;
	    col = lp->nbaz[c];
	    mpq_EGlpNumCopy (dj, lp->dz[c]);
	    vs = lp->vstat[col];
	    vt = lp->vtype[col];

	    if (cbnd == BSKIP) {
		if (mpq_EGlpNumIsEqual (dj, mpq_zeroLpNum, *dftol));
		else if (mpq_EGlpNumIsLess (dj, mpq_zeroLpNum) && vs == STAT_LOWER)
		    lp->vstat[col] = STAT_UPPER;
		else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, dj) && vs == STAT_UPPER)
		    lp->vstat[col] = STAT_LOWER;
	    } else if (c != eindex) {
		if (mpq_EGlpNumIsEqual (dj, mpq_zeroLpNum, *dftol))
		    f = 0;
		else if (mpq_EGlpNumIsLess (dj, mpq_zeroLpNum) &&
		    (vs == STAT_LOWER || vs == STAT_ZERO))
		    f = -1;
		else if (mpq_EGlpNumIsLess (mpq_zeroLpNum, dj) &&
		    (vs == STAT_UPPER || vs == STAT_ZERO))
		    f = 1;
		else
		    f = 0;

		if (f != lp->dfeas[c]) {
		    delta = f - lp->dfeas[c];
		    mcnt = lp->matcnt[col];
		    mbeg = lp->matbeg[col];
		    mpq_EGlpNumSet (dj, (double) (delta));
		    for (i = 0; i < mcnt; i++)
			mpq_EGlpNumAddInnProdTo (w[lp->matind[mbeg + i]], dj,
			    lp->matval[mbeg + i]);
		    mpq_EGlpNumAddInnProdTo (*dty, dj, lp->zA.coef[k]);
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
	    if (mpq_EGlpNumIsNeqqZero (w[i])) {
		mpq_EGlpNumCopy (srhs->coef[nz], w[i]);
		srhs->indx[nz] = i;
		nz++;
		mpq_EGlpNumZero (w[i]);
	    }
    }
    srhs->nzcnt = nz;
    mpq_EGlpNumClearVar (dj);
    mpq_EGlpNumClearVar (tz);
}

void mpq_ILLfct_compute_dpIy (mpq_lpinfo * lp,
      mpq_svector * srhs,
      mpq_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	mpq_ILLbasis_column_solve (lp, srhs, ssoln);
    }
}

void mpq_ILLfct_update_dpI_prices (mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_svector * srhs,
      mpq_svector * ssoln,
      int lindex,
      mpq_t alpha)
{
    int i;
    mpq_t ntmp;
    mpq_EGlpNumInitVar (ntmp);
    mpq_EGlpNumZero (ntmp);

    if (srhs->nzcnt == 0) {
	mpq_ILLfct_update_xz (lp, alpha, -1, -1);
    } else {
	mpq_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
	mpq_EGlpNumAddTo (ntmp, alpha);
	mpq_EGlpNumSign (ntmp);
	mpq_add_vectors (lp, ssoln, &(lp->yjz), &(lp->yjz), ntmp);
	mpq_EGlpNumSign (ntmp);
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    mpq_EGlpNumAddTo (lp->xbz[lp->yjz.indx[i]], lp->yjz.coef[i]);
    }
    mpq_EGlpNumSet (lp->xbz[lindex], ((double) (-lp->upd.fs)));
    mpq_EGlpNumAddTo (lp->xbz[lindex], ntmp);

    if (pinf->d_strategy == COMPLETE_PRICING) {
	mpq_ILLprice_compute_primal_inf (lp, pinf, lp->yjz.indx, lp->yjz.nzcnt,
	    DUAL_PHASEI);
	mpq_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEI);
	mpq_ILLfct_update_counts (lp, CNT_YRAVG, lp->yjz.nzcnt, mpq_zeroLpNum);
    } else
	mpq_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEI, ROW_PRICING);
    mpq_EGlpNumClearVar (ntmp);
}

void mpq_ILLfct_update_dIIfeas (mpq_lpinfo * lp,
      int eindex,
      mpq_svector * srhs)
{
    int j, k;
    int col, indx, vs;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    mpq_t *zAj, *l, *u;
    mpq_t *dty = &(lp->upd.dty);
    mpq_t *t_max = &(lp->upd.tz);
    mpq_t *t = lp->upd.t;
    mpq_t delta;
    mpq_svector a;
    mpq_EGlpNumInitVar (delta);
    mpq_EGlpNumZero (delta);
    mpq_EGlpNumZero (*dty);

    srhs->nzcnt = 0;
    for (j = 0; j < tctr && mpq_EGlpNumIsLeq (t[perm[j]], *t_max); j++) {
	k = ix[perm[j]];
	indx = lp->zA.indx[k];

	if (indx != eindex) {
	    zAj = &(lp->zA.coef[k]);
	    col = lp->nbaz[indx];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);
	    vs = lp->vstat[col];
	    if (vs == STAT_UPPER)
		mpq_EGlpNumCopyDiff (delta, *l, *u);
	    else
		mpq_EGlpNumCopyDiff (delta, *u, *l);
	    mpq_EGlpNumAddInnProdTo (*dty, delta, *zAj);
	    lp->vstat[col] = (vs == STAT_UPPER) ? STAT_LOWER : STAT_UPPER;

	    a.nzcnt = lp->matcnt[col];
	    a.indx = &(lp->matind[lp->matbeg[col]]);
	    a.coef = &(lp->matval[lp->matbeg[col]]);
	    mpq_add_vectors (lp, srhs, &a, srhs, delta);
	}
    }
    mpq_EGlpNumClearVar (delta);
}

void mpq_ILLfct_compute_dpIIy (mpq_lpinfo * lp,
      mpq_svector * srhs,
      mpq_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	mpq_ILLbasis_column_solve (lp, srhs, ssoln);
    }
}

void mpq_ILLfct_update_dpII_prices (mpq_lpinfo * lp, mpq_price_info * pinf,
      mpq_svector * srhs, mpq_svector * ssoln, int lindex, mpq_t eval,
      mpq_t alpha)
{
    int i;
    mpq_svector *u;

    if (srhs->nzcnt == 0) {
	mpq_ILLfct_update_xz (lp, alpha, -1, -1);
	u = &(lp->yjz);
    } else {
	if (ssoln->nzcnt != 0)
	    for (i = 0; i < ssoln->nzcnt; i++)
		mpq_EGlpNumSubTo (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
	mpq_ILLfct_update_xz (lp, alpha, -1, -1);
	mpq_add_vectors (lp, ssoln, &(lp->yjz), ssoln, mpq_oneLpNum);
	u = ssoln;
    }
    mpq_EGlpNumCopySum (lp->xbz[lindex], eval, alpha);

    if (pinf->d_strategy == COMPLETE_PRICING) {
	mpq_ILLprice_compute_primal_inf (lp, pinf, u->indx, u->nzcnt, DUAL_PHASEII);
	mpq_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEII);
	mpq_ILLfct_update_counts (lp, CNT_YRAVG, u->nzcnt, mpq_zeroLpNum);
    } else
	mpq_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEII, ROW_PRICING);
}

int mpq_ILLfct_test_pivot (mpq_lpinfo * lp,
      int indx,
      int indxtype,
      mpq_t piv_val)
{
    int i;
    mpq_t pval, ntmp;
    mpq_EGlpNumInitVar (pval);
    mpq_EGlpNumInitVar (ntmp);
    mpq_EGlpNumZero (pval);

    if (indxtype == ROW_PIVOT) {
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    if (lp->yjz.indx[i] == indx) {
		mpq_EGlpNumCopy (pval, lp->yjz.coef[i]);
		break;
	    }
    } else {
	for (i = 0; i < lp->zA.nzcnt; i++)
	    if (lp->zA.indx[i] == indx) {
		mpq_EGlpNumCopy (pval, lp->zA.coef[i]);
		break;
	    }
    }
    mpq_EGlpNumCopyDiff (ntmp, pval, piv_val);
    mpq_EGlpNumDivTo (ntmp, piv_val);
    if (mpq_EGlpNumIsLess (ntmp, mpq_zeroLpNum))
	mpq_EGlpNumSign (ntmp);
    if (mpq_EGlpNumIsLess (mpq_ALTPIV_TOLER, ntmp)) {
#if mpq_FCT_DEBUG > 1
	if (indxtype == ROW_PIVOT)
	    printf ("y_i = %.8f, z_j = %.8f %la %la\n", mpq_EGlpNumToLf (pval),
		mpq_EGlpNumToLf (piv_val), mpq_EGlpNumToLf (mpq_ALTPIV_TOLER),
		mpq_EGlpNumToLf (ntmp));
	else
	    printf ("z_j = %.8f, y_i = %.8f\n", mpq_EGlpNumToLf (pval),
		mpq_EGlpNumToLf (piv_val));
#endif
	mpq_EGlpNumClearVar (ntmp);
	mpq_EGlpNumClearVar (pval);
	return 1;
    }
    mpq_EGlpNumClearVar (pval);
    mpq_EGlpNumClearVar (ntmp);
    return 0;
}

#if mpq_FCT_DEBUG > 0

void mpq_fct_test_workvector (mpq_lpinfo * lp)
{
    int i, err = 0;
    for (i = 0; i < lp->ncols; i++) {
	if (mpq_EGlpNumIsNeqqZero (lp->work.coef[i])) {
	    err++;
	    mpq_EGlpNumZero (lp->work.coef[i]);
	}
	if (lp->iwork[i] != 0) {
	    err++;
	    lp->iwork[i] = 0;
	}
    }
    if (err)
	printf ("bad work vector, err=%d\n", err);
}

void mpq_fct_test_pfeasible (mpq_lpinfo * lp)
{
    int i, col;
    int err = 0;
    mpq_t *ftol = &(lp->tol->pfeas_tol);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];

	if (mpq_EGlpNumIsNeqq (lp->uz[col], mpq_INFTY)
	    && mpq_EGlpNumIsSumLess (*ftol, lp->uz[col], lp->xbz[i])) {
	    if (lp->bfeas[i] != 1) {
		err++;
		lp->bfeas[i] = 1;
	    }
	} else if (mpq_EGlpNumIsNeqq (lp->lz[col], mpq_NINFTY)
	    && mpq_EGlpNumIsSumLess (lp->xbz[i], *ftol, lp->lz[col])) {
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

void mpq_fct_test_dfeasible (mpq_lpinfo * lp)
{
    int j, col;
    int err = 0;
    mpq_t *ftol = &(lp->tol->dfeas_tol);
    mpq_t mftol[1];
    mpq_EGlpNumInitVar (mftol[0]);
    mpq_EGlpNumCopyNeg (mftol[0], *ftol);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;
	if (mpq_EGlpNumIsLess (lp->dz[j], mftol[0]) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)) {
	    if (lp->dfeas[j] != -1) {
		err++;
		lp->dfeas[j] = -1;
	    }
	}
	if (mpq_EGlpNumIsLess (*ftol, lp->dz[j]) &&
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

void mpq_fct_test_pI_x (mpq_lpinfo * lp,
      mpq_price_info * p)
{
    int i;
    int ern = 0;
    mpq_t *x;
    mpq_t err, diff;
    mpq_EGlpNumInitVar (err);
    mpq_EGlpNumInitVar (diff);
    mpq_EGlpNumZero (err);
    x = mpq_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumCopy (x[i], lp->xbz[i]);
    mpq_ILLfct_compute_phaseI_xbz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_PFEAS_TOLER, diff)) {
	    mpq_EGlpNumAddTo (err, diff);
	    ern++;
	    printf ("bad i = %d\n", i);
	}
    }
    if (mpq_EGlpNumIsNeqqZero (err))
	printf ("dI x err = %.7f, ern = %d\n", mpq_EGlpNumToLf (err), ern);
    mpq_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
    mpq_EGlpNumFreeArray (x);
    mpq_EGlpNumClearVar (diff);
    mpq_EGlpNumClearVar (err);
}

void mpq_fct_test_pII_x (mpq_lpinfo * lp,
      mpq_price_info * p)
{
    int i;
    int ern = 0;
    mpq_t *x;
    mpq_t err, diff;
    mpq_EGlpNumInitVar (err);
    mpq_EGlpNumInitVar (diff);
    mpq_EGlpNumZero (err);
    x = mpq_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumCopy (x[i], lp->xbz[i]);
    mpq_ILLfct_compute_xbz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_PFEAS_TOLER, diff)) {
	    mpq_EGlpNumAddTo (err, diff);
	    ern++;
	    printf ("bad i = %d\n", i);
	}
    }
    if (mpq_EGlpNumIsNeqqZero (err))
	printf ("dII x err = %.7f, ern = %d\n", mpq_EGlpNumToLf (err), ern);
    mpq_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
    mpq_EGlpNumFreeArray (x);
    mpq_EGlpNumClearVar (diff);
    mpq_EGlpNumClearVar (err);
}

void mpq_fct_test_pI_pi_dz (mpq_lpinfo * lp,
      mpq_price_info * p)
{
    int i;
    int ern = 0;
    mpq_t *pidz;
    mpq_t err, diff;
    mpq_EGlpNumInitVar (err);
    mpq_EGlpNumInitVar (diff);
    pidz = mpq_EGlpNumAllocArray (lp->ncols);
    mpq_EGlpNumZero (err);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumCopy (pidz[i], lp->pIpiz[i]);
    mpq_ILLfct_compute_phaseI_piz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumCopyDiff (diff, pidz[i], lp->pIpiz[i]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_DFEAS_TOLER, diff)) {
	    mpq_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpq_EGlpNumIsNeqqZero (err))
	printf ("pI pi err = %.7f, ern = %d\n", mpq_EGlpNumToLf (err), ern);

    mpq_EGlpNumZero (err);
    ern = 0;
    for (i = 0; i < lp->nnbasic; i++)
	mpq_EGlpNumCopy (pidz[i], lp->pIdz[i]);
    mpq_ILLfct_compute_phaseI_dz (lp);
    for (i = 0; i < lp->nnbasic; i++) {
	mpq_EGlpNumCopyDiff (diff, pidz[i], lp->pIdz[i]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_DFEAS_TOLER, diff)) {
	    mpq_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpq_EGlpNumIsNeqqZero (err))
	printf ("pI dz err = %.7f, ern = %d\n", mpq_EGlpNumToLf (err), ern);
    mpq_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
    mpq_EGlpNumClearVar (err);
    mpq_EGlpNumClearVar (diff);
    mpq_EGlpNumFreeArray (pidz);
}

void mpq_fct_test_pII_pi_dz (mpq_lpinfo * lp,
      mpq_price_info * p)
{
    int i;
    int ern = 0;
    mpq_t *pidz;
    mpq_t err, diff;
    mpq_EGlpNumInitVar (err);
    mpq_EGlpNumInitVar (diff);
    mpq_EGlpNumZero (err);
    pidz = mpq_EGlpNumAllocArray (lp->ncols);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumCopy (pidz[i], lp->piz[i]);
    mpq_ILLfct_compute_piz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumCopyDiff (diff, pidz[i], lp->piz[i]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_DFEAS_TOLER, diff)) {
	    mpq_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpq_EGlpNumIsNeqqZero (err))
	printf ("pII pi err = %.7f, ern = %d\n", mpq_EGlpNumToLf (err), ern);

    mpq_EGlpNumZero (err);
    ern = 0;
    for (i = 0; i < lp->nnbasic; i++)
	mpq_EGlpNumCopy (pidz[i], lp->dz[i]);
    mpq_ILLfct_compute_dz (lp);
    for (i = 0; i < lp->nnbasic; i++) {
	mpq_EGlpNumCopyDiff (diff, pidz[i], lp->dz[i]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_DFEAS_TOLER, diff)) {
	    mpq_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpq_EGlpNumIsNeqqZero (err))
	printf ("pII dz err = %.7f, ern = %d\n", mpq_EGlpNumToLf (err), ern);
    /*
     * mpq_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
     */
    mpq_EGlpNumClearVar (err);
    mpq_EGlpNumClearVar (diff);
    mpq_EGlpNumFreeArray (pidz);
}

#endif
