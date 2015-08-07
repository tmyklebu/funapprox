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

/* RCS_INFO = "$RCSfile: mpf_fct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#define mpf_FCT_DEBUG 0

#include "econfig.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdefs.h"
#include "stddefs.h"
#include "mpf_basis.h"
#include "mpf_fct.h"
#include "mpf_price.h"
#include "mpf_ratio.h"
#include "mpf_dstruct.h"

mpf_bndinfo *mpf_ILLfct_new_bndinfo (void)
{
    mpf_bndinfo *nbnd = (mpf_bndinfo *) malloc (sizeof (mpf_bndinfo));
    if (!nbnd) {
	fprintf (stderr, "not enough memory, in %s\n", __func__);
	exit (1);
    }
    mpf_EGlpNumInitVar ((nbnd->pbound));
    mpf_EGlpNumInitVar ((nbnd->cbound));
    return nbnd;
}

void mpf_ILLfct_free_bndinfo (mpf_bndinfo * binfo)
{
    mpf_EGlpNumClearVar ((binfo->pbound));
    mpf_EGlpNumClearVar ((binfo->cbound));
    ILL_IFFREE (binfo, mpf_bndinfo);
    return;
}

static int mpf_compute_zA1 (mpf_lpinfo * lp,
      mpf_svector * z,
      mpf_svector * zA,
      mpf_t ztoler),
/*
  compute_zA2 (mpf_lpinfo * lp,
                             mpf_svector * z,
                             mpf_svector * zA,
                             const mpf_t* ztoler), */
  mpf_compute_zA3 (mpf_lpinfo * lp,
      mpf_svector * z,
      mpf_svector * zA,
      mpf_t ztoler),
  mpf_expand_var_bounds (mpf_lpinfo * lp,
      mpf_t ftol,
      int *chgb),
  mpf_expand_var_coefs (mpf_lpinfo * lp,
      mpf_t ftol,
      int *chgc);

static void mpf_update_piv_values (mpf_count_struct * c,
      int phase,
      mpf_t piv),
/* copy_vectors (mpf_svector * a, mpf_svector * b), */
  mpf_add_vectors (mpf_lpinfo * lp,
      mpf_svector * a,
      mpf_svector * b,
      mpf_svector * c,
      mpf_t t);

static double mpf_my_rand (int bound,
      ILLrandstate * r);


void mpf_ILLfct_load_workvector (mpf_lpinfo * lp,
      mpf_svector * s)
{
    int i;

    for (i = 0; i < s->nzcnt; i++) {
	lp->work.indx[i] = s->indx[i];
	mpf_EGlpNumCopy (lp->work.coef[s->indx[i]], s->coef[i]);
    }
    lp->work.nzcnt = s->nzcnt;
}

void mpf_ILLfct_zero_workvector (mpf_lpinfo * lp)
{
    int i;

    for (i = 0; i < lp->work.nzcnt; i++)
	mpf_EGlpNumZero (lp->work.coef[lp->work.indx[i]]);
    lp->work.nzcnt = 0;
}

void mpf_ILLfct_set_variable_type (mpf_lpinfo * lp)
{
    int j;

    for (j = 0; j < lp->ncols; j++) {

	if (lp->matcnt[j] == 1 && lp->O->rowmap[lp->matind[lp->matbeg[j]]] == j)
	    lp->vclass[j] = CLASS_LOGICAL;
	else
	    lp->vclass[j] = CLASS_STRUCT;
	switch ((mpf_EGlpNumIsEqqual (lp->uz[j], mpf_INFTY) ? 1U : 0U) |
	    (mpf_EGlpNumIsEqqual (lp->lz[j], mpf_NINFTY) ? 2U : 0U)) {
	case 0:
	    if (mpf_EGlpNumIsLess (lp->lz[j], lp->uz[j]))
		lp->vtype[j] = VBOUNDED;
	    else if (mpf_EGlpNumIsEqqual (lp->lz[j], mpf_zeroLpNum) &&
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

void mpf_ILLfct_compute_pobj (mpf_lpinfo * lp)
{
    int i, j;
    int col;
    mpf_t sum;
    mpf_EGlpNumInitVar (sum);
    mpf_EGlpNumZero (sum);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumAddInnProdTo (sum, lp->cz[lp->baz[i]], lp->xbz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	if (lp->vstat[col] == STAT_UPPER)
	    mpf_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER)
	    mpf_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->lz[col]);
    }
    mpf_EGlpNumCopy (lp->pobjval, sum);
    mpf_EGlpNumCopy (lp->objval, sum);
    mpf_EGlpNumClearVar (sum);
}

void mpf_ILLfct_compute_dobj (mpf_lpinfo * lp)
{
    int i, j;
    int col;
    mpf_t sum;
    mpf_EGlpNumInitVar (sum);
    mpf_EGlpNumZero (sum);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumAddInnProdTo (sum, lp->piz[i], lp->bz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	if (lp->vstat[col] == STAT_UPPER)
	    mpf_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER)
	    mpf_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->lz[col]);
    }
    mpf_EGlpNumCopy (lp->dobjval, sum);
    mpf_EGlpNumCopy (lp->objval, sum);
    mpf_EGlpNumClearVar (sum);
}

void mpf_ILLfct_compute_xbz (mpf_lpinfo * lp)
{
    int i, j, r;
    int col, mcnt, mbeg;
    mpf_svector *srhs = &(lp->srhs);
    mpf_svector *ssoln = &(lp->ssoln);
    mpf_t xval;
    mpf_EGlpNumInitVar (xval);

    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumZero (lp->xbz[i]);
	mpf_EGlpNumCopy (srhs->coef[i], lp->bz[i]);
    }
    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];
	mpf_EGlpNumZero (xval);
	if (lp->vstat[col] == STAT_UPPER && mpf_EGlpNumIsNeqqZero (lp->uz[col]))
	    mpf_EGlpNumCopy (xval, lp->uz[col]);
	else if (lp->vstat[col] == STAT_LOWER && mpf_EGlpNumIsNeqqZero (lp->lz[col]))
	    mpf_EGlpNumCopy (xval, lp->lz[col]);

	if (mpf_EGlpNumIsNeqqZero (xval)) {
	    mcnt = lp->matcnt[col];
	    mbeg = lp->matbeg[col];
	    for (i = 0; i < mcnt; i++)
		mpf_EGlpNumSubInnProdTo (srhs->coef[lp->matind[mbeg + i]], xval,
		    lp->matval[mbeg + i]);
	}
    }
    for (i = 0, r = 0; i < lp->nrows; i++)
	if (mpf_EGlpNumIsNeqqZero (srhs->coef[i])) {
	    mpf_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
	    srhs->indx[r] = i;
	    r++;
	}
    srhs->nzcnt = r;

    mpf_ILLbasis_column_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpf_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
    mpf_EGlpNumClearVar (xval);
}

void mpf_ILLfct_compute_piz (mpf_lpinfo * lp)
{
    int i, r;
    mpf_svector *srhs = &(lp->srhs);
    mpf_svector *ssoln = &(lp->ssoln);

    for (i = 0, r = 0; i < lp->nrows; i++) {
	mpf_EGlpNumZero (lp->piz[i]);
	if (mpf_EGlpNumIsNeqqZero (lp->cz[lp->baz[i]])) {
	    srhs->indx[r] = i;
	    mpf_EGlpNumCopy (srhs->coef[r], lp->cz[lp->baz[i]]);
	    r++;
	}
    }
    srhs->nzcnt = r;

    mpf_ILLbasis_row_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpf_EGlpNumCopy (lp->piz[ssoln->indx[i]], ssoln->coef[i]);
}

void mpf_ILLfct_compute_dz (mpf_lpinfo * lp)
{
    int i, j;
    int col;
    int mcnt, mbeg;
    mpf_t sum;
    mpf_EGlpNumInitVar (sum);

    for (j = 0; j < lp->nnbasic; j++) {
	mpf_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpf_EGlpNumCopyDiff (lp->dz[j], lp->cz[col], sum);
    }
    mpf_EGlpNumClearVar (sum);
}

void mpf_ILLfct_compute_phaseI_xbz (mpf_lpinfo * lp)
{
    int i, j, r;
    int col, mcnt, mbeg;
    mpf_svector *srhs = &(lp->srhs);
    mpf_svector *ssoln = &(lp->ssoln);

    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumZero (lp->xbz[i]);
	mpf_EGlpNumZero (srhs->coef[i]);
    }
    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];

	if (lp->dfeas[j]) {
	    mcnt = lp->matcnt[col];
	    mbeg = lp->matbeg[col];
	    if (lp->dfeas[j] == -1)
		for (i = 0; i < mcnt; i++)
		    mpf_EGlpNumSubTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
	    else
		for (i = 0; i < mcnt; i++)
		    mpf_EGlpNumAddTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
	}
    }
    for (i = 0, r = 0; i < lp->nrows; i++)
	if (mpf_EGlpNumIsNeqqZero (srhs->coef[i])) {
	    mpf_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
	    srhs->indx[r] = i;
	    r++;
	}
    srhs->nzcnt = r;

    mpf_ILLbasis_column_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpf_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
}

void mpf_ILLfct_compute_phaseI_piz (mpf_lpinfo * lp)
{
    int i, r;
    mpf_svector *srhs = &(lp->srhs);
    mpf_svector *ssoln = &(lp->ssoln);

    for (i = 0, r = 0; i < lp->nrows; i++) {
	mpf_EGlpNumZero (lp->pIpiz[i]);
	if (lp->bfeas[i] != 0) {
	    srhs->indx[r] = i;
	    mpf_EGlpNumSet (srhs->coef[r], (double) lp->bfeas[i]);
	    r++;
	}
    }
    srhs->nzcnt = r;

    mpf_ILLbasis_row_solve (lp, srhs, ssoln);
    for (i = 0; i < ssoln->nzcnt; i++)
	mpf_EGlpNumCopy (lp->pIpiz[ssoln->indx[i]], ssoln->coef[i]);
    mpf_ILLfct_update_counts (lp, CNT_P1PINZ, ssoln->nzcnt, mpf_zeroLpNum);
}

void mpf_ILLfct_compute_phaseI_dz (mpf_lpinfo * lp)
{
    int i, j;
    int col;
    int mcnt, mbeg;
    mpf_t sum;
    mpf_EGlpNumInitVar (sum);
    ILL_IFTRACE ("%s\n", __func__);

    for (j = 0; j < lp->nnbasic; j++) {
	mpf_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpf_EGlpNumCopyNeg (lp->pIdz[j], sum);
	ILL_IFTRACE ("%d:%d:%lf:%la\n", j, col, mpf_EGlpNumToLf (sum),
	    mpf_EGlpNumToLf (sum));
    }
    mpf_EGlpNumClearVar (sum);
}

void mpf_ILLfct_compute_yz (mpf_lpinfo * lp,
      mpf_svector * yz,
      mpf_svector * updz,
      int col)
{
    mpf_svector a;

    a.nzcnt = lp->matcnt[col];
    a.indx = &(lp->matind[lp->matbeg[col]]);
    a.coef = &(lp->matval[lp->matbeg[col]]);

    mpf_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, mpf_PIVZ_TOLER);
    if (updz)
	mpf_ILLbasis_column_solve_update (lp, &a, updz, yz);
    else
	mpf_ILLbasis_column_solve (lp, &a, yz);
    mpf_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, mpf_SZERO_TOLER);
}

void mpf_ILLfct_compute_zz (mpf_lpinfo * lp,
      mpf_svector * zz,
      int row)
{
    mpf_ILLfct_compute_binvrow (lp, zz, row, mpf_PIVZ_TOLER);
}

void mpf_ILLfct_compute_binvrow (mpf_lpinfo * lp,
      mpf_svector * zz,
      int row,
      mpf_t ztoler)
{
    mpf_svector a;
    mpf_t e;
    mpf_EGlpNumInitVar (e);
    mpf_EGlpNumOne (e);

    a.nzcnt = 1;
    a.coef = &e;
    a.indx = &row;

    if (mpf_EGlpNumIsLess (mpf_zeroLpNum, ztoler))
	mpf_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, ztoler);
    mpf_ILLbasis_row_solve (lp, &a, zz);
    if (mpf_EGlpNumIsLess (mpf_zeroLpNum, ztoler))
	mpf_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, mpf_SZERO_TOLER);
    mpf_EGlpNumClearVar (e);
}

void mpf_ILLfct_compute_psteep_upv (mpf_lpinfo * lp,
      mpf_svector * swz)
{
    mpf_ILLbasis_row_solve (lp, &(lp->yjz), swz);
}

void mpf_ILLfct_compute_dsteep_upv (mpf_lpinfo * lp,
      mpf_svector * swz)
{
    mpf_ILLbasis_column_solve (lp, &(lp->zz), swz);
}

static int mpf_compute_zA1 (mpf_lpinfo * lp,
      mpf_svector * z,
      mpf_svector * zA,
      mpf_t ztoler)
{
    int rval = 0;
    int i, j, nz = 0;
    int col, mcnt, mbeg;
    mpf_t sum;
    mpf_t *v = 0;
    mpf_EGlpNumInitVar (sum);
    v = mpf_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumZero (v[i]);
    for (i = 0; i < z->nzcnt; i++)
	mpf_EGlpNumCopy (v[z->indx[i]], z->coef[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	mpf_EGlpNumZero (sum);
	col = lp->nbaz[j];
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (sum, v[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

	if (mpf_EGlpNumIsNeqZero (sum, ztoler)) {
	    mpf_EGlpNumCopy (zA->coef[nz], sum);
	    zA->indx[nz] = j;
	    nz++;
	}
    }
    zA->nzcnt = nz;

    mpf_EGlpNumClearVar (sum);
    mpf_EGlpNumFreeArray (v);
    ILL_RETURN (rval, "mpf_compute_zA1");
}


static int mpf_compute_zA3 (mpf_lpinfo * lp,
      mpf_svector * z,
      mpf_svector * zA,
      mpf_t ztoler)
{
    int rval = 0;
    int i, j, k, ix;
    int nz = 0;
    int row, col;
    int rcnt, rbeg;
    mpf_t val;
    mpf_EGlpNumInitVar (val);
    k = 0;
    for (i = 0; i < z->nzcnt; i++) {
	row = z->indx[i];
	mpf_EGlpNumCopy (val, z->coef[i]);
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
		mpf_EGlpNumAddInnProdTo (lp->work.coef[ix], val, lp->rowval[rbeg + j]);
	    }
	}
    }
    for (j = 0; j < k; j++) {
	ix = lp->work.indx[j];
	mpf_EGlpNumCopy (val, lp->work.coef[ix]);
	mpf_EGlpNumZero (lp->work.coef[ix]);
	lp->iwork[ix] = 0;
	if (mpf_EGlpNumIsNeqZero (val, ztoler)) {
	    mpf_EGlpNumCopy (zA->coef[nz], val);
	    zA->indx[nz] = ix;
	    nz++;
	}
    }
    zA->nzcnt = nz;
    mpf_EGlpNumClearVar (val);
    ILL_RETURN (rval, "mpf_compute_zA3");
}

int mpf_ILLfct_compute_zA (mpf_lpinfo * lp,
      mpf_svector * z,
      mpf_svector * zA)
{
    if (z->nzcnt < lp->nrows / 2)
	return mpf_compute_zA3 (lp, z, zA, mpf_PIVZ_TOLER);
    else
	return mpf_compute_zA1 (lp, z, zA, mpf_PIVZ_TOLER);
}

/* compute v^T A */
void mpf_ILLfct_compute_vA (mpf_lpinfo * lp,
      mpf_svector * v,
      mpf_t * vA)
{
    int i, j;
    int row, col;
    int rcnt, rbeg;
    mpf_t val;
    mpf_EGlpNumInitVar (val);

    for (j = 0; j < lp->ncols; j++)
	mpf_EGlpNumZero (vA[j]);

    for (i = 0; i < v->nzcnt; i++) {
	row = v->indx[i];
	mpf_EGlpNumCopy (val, v->coef[i]);
	rcnt = lp->rowcnt[row];
	rbeg = lp->rowbeg[row];
	for (j = 0; j < rcnt; j++) {
	    col = lp->rowind[rbeg + j];
	    mpf_EGlpNumAddInnProdTo (vA[col], val, lp->rowval[rbeg + j]);
	}
    }

    for (j = 0; j < lp->ncols; j++)
	if (mpf_EGlpNumIsEqual (vA[j], mpf_zeroLpNum, mpf_SZERO_TOLER))
	    mpf_EGlpNumZero (vA[j]);

    mpf_EGlpNumClearVar (val);
    return;
}

/* update information */

/*
1) lvstat - new status of leaving var.
*/
void mpf_ILLfct_update_basis_info (mpf_lpinfo * lp,
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

void mpf_ILLfct_update_xz (mpf_lpinfo * lp,
      mpf_t tz,
      int eindex,
      int lindex)
{
    int i, evar, estat;
    ILL_IFTRACE ("%s:%la:%d:%d:%d\n", __func__, mpf_EGlpNumToLf (tz), eindex,
	lindex, lp->yjz.nzcnt);

    if (mpf_EGlpNumIsNeqqZero (tz))
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    mpf_EGlpNumSubInnProdTo (lp->xbz[lp->yjz.indx[i]], tz, lp->yjz.coef[i]);

    if (lindex >= 0) {		/* variable leaves basis */
	evar = lp->nbaz[eindex];
	estat = lp->vstat[evar];
	if (estat == STAT_LOWER)
	    mpf_EGlpNumCopySum (lp->xbz[lindex], lp->lz[evar], tz);
	else if (estat == STAT_UPPER)
	    mpf_EGlpNumCopySum (lp->xbz[lindex], lp->uz[evar], tz);
	else if (estat == STAT_ZERO)
	    mpf_EGlpNumCopy (lp->xbz[lindex], tz);
    }
}

void mpf_ILLfct_update_piz (mpf_lpinfo * lp,
      mpf_t alpha)
{
    int i;

    for (i = 0; i < lp->zz.nzcnt; i++)
	mpf_EGlpNumAddInnProdTo (lp->piz[lp->zz.indx[i]], alpha, lp->zz.coef[i]);
}

void mpf_ILLfct_update_pIpiz (mpf_lpinfo * lp,
      mpf_svector * z,
      mpf_t alpha)
{
    int i;
    if (mpf_EGlpNumIsEqqual (alpha, mpf_zeroLpNum))
	return;
    if (mpf_EGlpNumIsEqqual (alpha, mpf_oneLpNum)) {
	for (i = 0; i < z->nzcnt; i++)
	    mpf_EGlpNumAddTo (lp->pIpiz[z->indx[i]], z->coef[i]);
    } else {
	for (i = 0; i < z->nzcnt; i++)
	    mpf_EGlpNumAddInnProdTo (lp->pIpiz[z->indx[i]], alpha, z->coef[i]);
    }
}

void mpf_ILLfct_update_dz (mpf_lpinfo * lp,
      int eindex,
      mpf_t alpha)
{
    int i;

    for (i = 0; i < lp->zA.nzcnt; i++)
	mpf_EGlpNumSubInnProdTo (lp->dz[lp->zA.indx[i]], alpha, lp->zA.coef[i]);
    mpf_EGlpNumCopyNeg (lp->dz[eindex], alpha);
}

void mpf_ILLfct_update_pIdz (mpf_lpinfo * lp,
      mpf_svector * zA,
      int eindex,
      mpf_t alpha)
{
    int i;
    if (mpf_EGlpNumIsEqqual (alpha, mpf_zeroLpNum))
	return;

    if (mpf_EGlpNumIsEqqual (alpha, mpf_oneLpNum)) {
	for (i = 0; i < zA->nzcnt; i++)
	    mpf_EGlpNumSubTo (lp->pIdz[zA->indx[i]], zA->coef[i]);
    } else {
	for (i = 0; i < zA->nzcnt; i++)
	    mpf_EGlpNumSubInnProdTo (lp->pIdz[zA->indx[i]], alpha, zA->coef[i]);
    }
    if (eindex > -1)
	mpf_EGlpNumCopyNeg (lp->pIdz[eindex], alpha);
}

/* bound and coef shift routines */

/* scale bound in mpf_my_rand to get more random digits, unless bound is
   large */
static double mpf_my_rand (int bound,
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

static int mpf_expand_var_bounds (mpf_lpinfo * lp,
      mpf_t ftol,
      int *chgb)
{
    int rval = 0;
    int i, col, nchg = 0;
    mpf_t newb, cftol;
    mpf_t *x, *l, *u;
    ILLrandstate r;
    mpf_EGlpNumInitVar (newb);
    mpf_EGlpNumInitVar (cftol);
    mpf_EGlpNumCopyAbs (cftol, ftol);
    mpf_EGlpNumDivUiTo (cftol, 10);

    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
	    continue;
	x = &(lp->xbz[i]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	/* we use newb as temporal variable outside the if's scope */
	mpf_EGlpNumCopyDiff (newb, *x, ftol);
	if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY) && mpf_EGlpNumIsLess (newb, *l)) {
	    mpf_EGlpNumSet (newb, -1.0 * (mpf_my_rand (50, &(lp->rstate)) + 1.0));
	    mpf_EGlpNumMultTo (newb, cftol);
	    if (mpf_EGlpNumIsLess (*x, *l))
		mpf_EGlpNumAddTo (newb, *x);
	    else
		mpf_EGlpNumAddTo (newb, *l);
	    rval = mpf_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
	mpf_EGlpNumCopySum (newb, *x, ftol);
	if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY) && mpf_EGlpNumIsLess (*u, newb)) {
	    mpf_EGlpNumSet (newb, mpf_my_rand (50, &(lp->rstate)) + 1.0);
	    mpf_EGlpNumMultTo (newb, cftol);
	    if (mpf_EGlpNumIsLess (*x, *u))
		mpf_EGlpNumAddTo (newb, *u);
	    else
		mpf_EGlpNumAddTo (newb, *x);
	    rval = mpf_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
    }
    *chgb = nchg;

CLEANUP:
    mpf_EGlpNumClearVar (newb);
    mpf_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "mpf_expand_var_bounds");
}

static int mpf_expand_phaseI_bounds (mpf_lpinfo * lp,
      int *chgb)
{
    int rval = 0;
    int i, col, nchg = 0;
    mpf_t newb, cftol;
    mpf_t *u, *l, *x;
    ILLrandstate r;
    mpf_EGlpNumInitVar (newb);
    mpf_EGlpNumInitVar (cftol);
    mpf_EGlpNumCopyAbs (cftol, lp->tol->ip_tol);
    mpf_EGlpNumDivUiTo (cftol, 10);
    ILLutil_sprand (1, &r);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
	    continue;
	x = &(lp->xbz[i]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY) && mpf_EGlpNumIsEqual (*x, *l, cftol)) {
	    mpf_EGlpNumSet (newb, mpf_my_rand (50, &(lp->rstate)) + 1.0);
	    mpf_EGlpNumMultTo (newb, cftol);
	    mpf_EGlpNumSign (newb);
	    mpf_EGlpNumAddTo (newb, *l);
	    rval = mpf_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
	if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY) && mpf_EGlpNumIsEqual (*x, *u, cftol)) {
	    mpf_EGlpNumSet (newb, mpf_my_rand (50, &(lp->rstate)) + 1.0);
	    mpf_EGlpNumMultTo (newb, cftol);
	    mpf_EGlpNumAddTo (newb, *u);
	    rval = mpf_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	}
    }
    *chgb = nchg;

CLEANUP:
    mpf_EGlpNumClearVar (newb);
    mpf_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "mpf_expand_phaseI_bounds");
}

int mpf_ILLfct_adjust_viol_bounds (mpf_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;
    mpf_t tol;
    mpf_EGlpNumInitVar (tol);
    mpf_EGlpNumCopyNeg (tol, lp->tol->pfeas_tol);
    rval = mpf_expand_var_bounds (lp, tol, &chgb);
#if mpf_FCT_DEBUG > 0
    if (rval == 0)
	printf ("adjusting %d bounds\n", chgb);
#endif
    mpf_EGlpNumClearVar (tol);
    ILL_RETURN (rval, "mpf_ILLfct_adjust_viol_bounds");
}

int mpf_ILLfct_perturb_bounds (mpf_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;

    rval = mpf_expand_var_bounds (lp, lp->tol->ip_tol, &chgb);
#if mpf_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d bounds\n", chgb);
#endif
    ILL_RETURN (rval, "mpf_ILLfct_perturb_bounds");
}

int mpf_ILLfct_perturb_phaseI_bounds (mpf_lpinfo * lp)
{
    int rval = 0;
    int chgb = 0;

    rval = mpf_expand_phaseI_bounds (lp, &chgb);
#if mpf_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d phase I bounds\n", chgb);
#endif
    ILL_RETURN (rval, "mpf_ILLfct_perturb_phaseI_bounds");
}

int mpf_ILLfct_bound_shift (mpf_lpinfo * lp,
      int col,
      int bndtype,
      mpf_t newbnd)
{
    int rval = 0;
    mpf_bndinfo *nbnd = 0;
    ILL_IFTRACE ("\n%s:%d:%d:%la", __func__, col, bndtype, mpf_EGlpNumToLf (newbnd));
    nbnd = mpf_ILLfct_new_bndinfo ();

    nbnd->varnum = col;
    nbnd->btype = bndtype;
    if (bndtype == BOUND_LOWER) {
	mpf_EGlpNumCopy (nbnd->pbound, lp->lz[col]);
	mpf_EGlpNumCopy (nbnd->cbound, newbnd);
	mpf_EGlpNumCopy (lp->lz[col], newbnd);
    } else {
	mpf_EGlpNumCopy (nbnd->pbound, lp->uz[col]);
	mpf_EGlpNumCopy (nbnd->cbound, newbnd);
	mpf_EGlpNumCopy (lp->uz[col], newbnd);
    }
    ILL_IFTRACE (":%la", mpf_EGlpNumToLf (nbnd->pbound));
    if (lp->vtype[col] == VFIXED || lp->vtype[col] == VARTIFICIAL) {
	/* printf ("changing f/a bound\n"); */
	if (mpf_EGlpNumIsLess (lp->lz[col], lp->uz[col]))
	    lp->vtype[col] = VBOUNDED;
    }
    nbnd->next = lp->bchanges;
    lp->bchanges = nbnd;
    lp->nbchange++;

    /* CLEANUP: */
    if (rval)
	mpf_ILLfct_free_bndinfo (nbnd);
    ILL_IFTRACE ("\n");
    ILL_RETURN (rval, "mpf_ILLfct_bound_shift");
}

void mpf_ILLfct_unroll_bound_change (mpf_lpinfo * lp)
{
    int col;
    int changex = 0;
    mpf_bndinfo *bptr = lp->bchanges;
    mpf_bndinfo *nptr = 0;
    ILL_IFTRACE ("%s:", __func__);

    while (lp->nbchange != 0) {
	col = bptr->varnum;
	ILL_IFTRACE (":%d", col);

	if (bptr->btype == BOUND_UPPER)
	    mpf_EGlpNumCopy (lp->uz[col], bptr->pbound);
	else
	    mpf_EGlpNumCopy (lp->lz[col], bptr->pbound);

	if (lp->vtype[col] == VBOUNDED) {
	    if (mpf_EGlpNumIsEqqual (lp->lz[col], lp->uz[col]))
		lp->vtype[col] = (mpf_EGlpNumIsEqqual (lp->lz[col], mpf_zeroLpNum)) ?
		    VARTIFICIAL : VFIXED;
	}
	if (lp->vstat[col] != STAT_BASIC) {
	    if ((bptr->btype == BOUND_UPPER && lp->vstat[col] == STAT_UPPER) ||
		(bptr->btype == BOUND_LOWER && lp->vstat[col] == STAT_LOWER))
		changex++;
	}
	nptr = bptr->next;
	mpf_EGlpNumClearVar ((bptr->cbound));
	mpf_EGlpNumClearVar ((bptr->pbound));
	ILL_IFFREE (bptr, mpf_bndinfo);
	bptr = nptr;
	lp->nbchange--;
    }
    lp->bchanges = bptr;
    ILL_IFTRACE ("\n");
    if (changex)
	mpf_ILLfct_compute_xbz (lp);
}

static int mpf_expand_var_coefs (mpf_lpinfo * lp,
      mpf_t ftol,
      int *chgc)
{
    int rval = 0;
    int i, col, vs, vt;
    int nchg = 0;
    mpf_t newc, cftol, mftol[1];
    mpf_t *c, *dj;
    ILLrandstate r;
    mpf_EGlpNumInitVar (newc);
    mpf_EGlpNumInitVar (cftol);
    mpf_EGlpNumInitVar (mftol[0]);
    mpf_EGlpNumCopyAbs (cftol, ftol);
    mpf_EGlpNumDivUiTo (cftol, 10);
    mpf_EGlpNumCopyNeg (mftol[0], ftol);
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
	    mpf_EGlpNumCopyDiff (newc, *c, *dj);
	    rval = mpf_ILLfct_coef_shift (lp, col, newc);
	    ILL_CLEANUP_IF (rval);
	    nchg++;
	    break;
	case STAT_LOWER:
	    if (mpf_EGlpNumIsLess (*dj, ftol)) {
		mpf_EGlpNumSet (newc, mpf_my_rand (50, &(lp->rstate)) + 1.0);
		mpf_EGlpNumMultTo (newc, cftol);
		mpf_EGlpNumAddTo (newc, *c);
		if (mpf_EGlpNumIsLess (*dj, mpf_zeroLpNum))
		    mpf_EGlpNumSubTo (newc, *dj);
		rval = mpf_ILLfct_coef_shift (lp, col, newc);
		ILL_CLEANUP_IF (rval);
		nchg++;
	    }
	    break;
	case STAT_UPPER:
	    if (mpf_EGlpNumIsLess (mftol[0], *dj)) {
		mpf_EGlpNumSet (newc, mpf_my_rand (50, &(lp->rstate)) + 1.0);
		mpf_EGlpNumMultTo (newc, cftol);
		mpf_EGlpNumSign (newc);
		mpf_EGlpNumAddTo (newc, *c);
		if (mpf_EGlpNumIsLess (mpf_zeroLpNum, *dj))
		    mpf_EGlpNumSubTo (newc, *dj);
		rval = mpf_ILLfct_coef_shift (lp, col, newc);
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
    mpf_EGlpNumClearVar (mftol[0]);
    mpf_EGlpNumClearVar (newc);
    mpf_EGlpNumClearVar (cftol);
    ILL_RETURN (rval, "mpf_expand_var_coefs");
}

int mpf_ILLfct_adjust_viol_coefs (mpf_lpinfo * lp)
{
    int rval = 0;
    int chgc = 0;
    mpf_t tol;
    mpf_EGlpNumInitVar (tol);
    mpf_EGlpNumCopyNeg (tol, lp->tol->dfeas_tol);

    rval = mpf_expand_var_coefs (lp, tol, &chgc);
#if mpf_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d coefs\n", chgc);
#endif
    mpf_EGlpNumClearVar (tol);
    ILL_RETURN (rval, "mpf_ILLfct_adjust_viol_coefs");
}

int mpf_ILLfct_perturb_coefs (mpf_lpinfo * lp)
{
    int rval = 0;
    int chgc = 0;

    rval = mpf_expand_var_coefs (lp, lp->tol->id_tol, &chgc);
#if mpf_FCT_DEBUG > 0
    if (rval == 0)
	printf ("perturbing %d coefs\n", chgc);
#endif
    ILL_RETURN (rval, "mpf_ILLfct_perturb_coefs");
}

int mpf_ILLfct_coef_shift (mpf_lpinfo * lp,
      int col,
      mpf_t newcoef)
{
    int rval = 0;
    mpf_coefinfo *ncoef = 0;

    ILL_SAFE_MALLOC (ncoef, 1, mpf_coefinfo);
    mpf_EGlpNumInitVar ((ncoef->pcoef));
    mpf_EGlpNumInitVar ((ncoef->ccoef));

    ncoef->varnum = col;
    mpf_EGlpNumCopy (ncoef->pcoef, lp->cz[col]);
    mpf_EGlpNumCopy (ncoef->ccoef, newcoef);
    mpf_EGlpNumCopy (lp->cz[col], newcoef);
    ncoef->next = lp->cchanges;
    lp->cchanges = ncoef;
    mpf_EGlpNumAddTo (lp->dz[lp->vindex[col]], ncoef->ccoef);
    mpf_EGlpNumSubTo (lp->dz[lp->vindex[col]], ncoef->pcoef);
    lp->ncchange++;

CLEANUP:
    if (rval) {
	mpf_EGlpNumClearVar ((ncoef->pcoef));
	mpf_EGlpNumClearVar ((ncoef->ccoef));
	ILL_IFFREE (ncoef, mpf_coefinfo);
    }
    ILL_RETURN (rval, "mpf_ILLfct_coef_shift");
}

void mpf_ILLfct_unroll_coef_change (mpf_lpinfo * lp)
{
    int bascoef = 0;
    mpf_coefinfo *cptr = (mpf_coefinfo *) lp->cchanges;
    mpf_coefinfo *nptr = 0;

    while (lp->ncchange != 0) {
	mpf_EGlpNumCopy (lp->cz[cptr->varnum], cptr->pcoef);
	if (lp->vstat[cptr->varnum] != STAT_BASIC) {
	    mpf_EGlpNumAddTo (lp->dz[lp->vindex[cptr->varnum]], cptr->pcoef);
	    mpf_EGlpNumSubTo (lp->dz[lp->vindex[cptr->varnum]], cptr->ccoef);
	} else
	    bascoef++;

	nptr = cptr->next;
	mpf_EGlpNumClearVar ((cptr->pcoef));
	mpf_EGlpNumClearVar ((cptr->ccoef));
	ILL_IFFREE (cptr, mpf_coefinfo);
	cptr = nptr;
	lp->ncchange--;
    }
    lp->cchanges = cptr;
    if (bascoef) {
	mpf_ILLfct_compute_piz (lp);
	mpf_ILLfct_compute_dz (lp);
    }
}

/* feasibility routines */
void mpf_ILLfct_check_pfeasible (mpf_lpinfo * lp,
      mpf_feas_info * fs,
      mpf_t ftol)
{
    int i, col;
    mpf_t infeas, err1, err2;
    mpf_EGlpNumInitVar (infeas);
    mpf_EGlpNumInitVar (err1);
    mpf_EGlpNumInitVar (err2);
    mpf_EGlpNumZero (infeas);
    fs->pstatus = PRIMAL_FEASIBLE;
    mpf_EGlpNumZero (fs->totinfeas);
    ILL_IFTRACE ("%s:tol %la\n", __func__, mpf_EGlpNumToLf (ftol));

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];
	mpf_EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
	mpf_EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
	if (mpf_EGlpNumIsLess (ftol, err1)
	    && mpf_EGlpNumIsNeq (lp->uz[col], mpf_INFTY, mpf_oneLpNum)) {
	    mpf_EGlpNumAddTo (infeas, err1);
	    WARNING (mpf_EGlpNumIsLess (mpf_INFTY, err1),
		"This is imposible lu = %15lg xbz = %15lg" " mpf_INFTY = %15lg",
		mpf_EGlpNumToLf (lp->uz[col]), mpf_EGlpNumToLf (lp->xbz[i]),
		mpf_EGlpNumToLf (mpf_INFTY));
	    lp->bfeas[i] = 1;
	} else if (mpf_EGlpNumIsLess (ftol, err2)
	    && mpf_EGlpNumIsNeq (lp->lz[col], mpf_NINFTY, mpf_oneLpNum)) {
	    mpf_EGlpNumAddTo (infeas, err2);
	    WARNING (mpf_EGlpNumIsLess (mpf_INFTY, err2),
		"This is imposible lz = %15lg xbz = %15lg" " mpf_NINFTY = %15lg",
		mpf_EGlpNumToLf (lp->lz[col]), mpf_EGlpNumToLf (lp->xbz[i]),
		mpf_EGlpNumToLf (mpf_NINFTY));
	    lp->bfeas[i] = -1;
	} else
	    lp->bfeas[i] = 0;
    }
    if (mpf_EGlpNumIsNeqqZero (infeas)) {
	fs->pstatus = PRIMAL_INFEASIBLE;
	mpf_EGlpNumCopy (fs->totinfeas, infeas);
	ILL_IFTRACE ("%s:inf %la\n", __func__, mpf_EGlpNumToLf (infeas));
	if (mpf_EGlpNumIsLess (fs->totinfeas, mpf_zeroLpNum)) {
	    printf ("Negative infeasibility, Imposible! %lf %la\n",
		mpf_EGlpNumToLf (infeas), mpf_EGlpNumToLf (infeas));
	}
    }
    mpf_EGlpNumCopy (lp->pinfeas, infeas);
    mpf_EGlpNumClearVar (infeas);
    mpf_EGlpNumClearVar (err1);
    mpf_EGlpNumClearVar (err2);
}

/* feasibility routines */
void mpf_ILLfct_check_pIpfeasible (mpf_lpinfo * lp,
      mpf_feas_info * fs,
      mpf_t ftol)
{
    int i, col;
    int ninf = 0;

    fs->pstatus = PRIMAL_FEASIBLE;
    mpf_EGlpNumZero (fs->totinfeas);

    for (i = 0; i < lp->nrows; i++) {
	if (mpf_EGlpNumIsEqual (lp->xbz[i], mpf_zeroLpNum, ftol))
	    continue;
	col = lp->baz[i];
	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, lp->xbz[i]) &&
	    mpf_EGlpNumIsNeqq (lp->uz[col], mpf_INFTY)) {
	    ninf++;
	} else if (mpf_EGlpNumIsLess (lp->xbz[i], mpf_zeroLpNum) &&
	    mpf_EGlpNumIsNeqq (lp->lz[col], mpf_NINFTY)) {
	    ninf++;
	}
    }
    if (ninf != 0)
	fs->pstatus = PRIMAL_INFEASIBLE;
}

void mpf_ILLfct_check_dfeasible (mpf_lpinfo * lp,
      mpf_feas_info * fs,
      mpf_t ftol)
{
    int j, col;
    mpf_t infeas;
    mpf_EGlpNumInitVar (infeas);
    mpf_EGlpNumZero (infeas);
    fs->dstatus = DUAL_FEASIBLE;
    mpf_EGlpNumZero (fs->totinfeas);

    for (j = 0; j < lp->nnbasic; j++) {
	lp->dfeas[j] = 0;
	if (mpf_EGlpNumIsEqual (lp->dz[j], mpf_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	if (mpf_EGlpNumIsLess (lp->dz[j], mpf_zeroLpNum) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)) {
	    mpf_EGlpNumSubTo (infeas, lp->dz[j]);
	    lp->dfeas[j] = -1;
	} else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, lp->dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO)) {
	    mpf_EGlpNumAddTo (infeas, lp->dz[j]);
	    lp->dfeas[j] = 1;
	}
    }

    if (mpf_EGlpNumIsNeqqZero (infeas)) {
	mpf_EGlpNumCopy (fs->totinfeas, infeas);
	fs->dstatus = DUAL_INFEASIBLE;
	ILL_IFTRACE ("%s:inf %la\n", __func__, mpf_EGlpNumToLf (infeas));
	if (mpf_EGlpNumIsLess (fs->totinfeas, mpf_zeroLpNum)) {
	    printf ("Negative infeasibility, Imposible! %lf %la\n",
		mpf_EGlpNumToLf (infeas), mpf_EGlpNumToLf (infeas));
	}
    }
    mpf_EGlpNumCopy (lp->dinfeas, infeas);
    mpf_EGlpNumClearVar (infeas);
}

void mpf_ILLfct_check_pIdfeasible (mpf_lpinfo * lp,
      mpf_feas_info * fs,
      mpf_t ftol)
{
    int j, col;
    int ninf = 0;
    mpf_t *dz = lp->pIdz;

    fs->dstatus = DUAL_FEASIBLE;

    for (j = 0; j < lp->nnbasic; j++) {
	if (mpf_EGlpNumIsEqual (dz[j], mpf_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;

	if (mpf_EGlpNumIsLess (dz[j], mpf_zeroLpNum) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
	    ninf++;
	else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, dz[j]) &&
	    (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
	    ninf++;
    }

    if (ninf != 0)
	fs->dstatus = DUAL_INFEASIBLE;
}

void mpf_ILLfct_dual_adjust (mpf_lpinfo * lp,
      mpf_t ftol)
{
    int j, col;

    for (j = 0; j < lp->nnbasic; j++) {
	if (mpf_EGlpNumIsEqual (lp->dz[j], mpf_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (mpf_EGlpNumIsLess (lp->dz[j], mpf_zeroLpNum) &&
	    mpf_EGlpNumIsNeqq (lp->uz[col], mpf_INFTY))
	    lp->vstat[col] = STAT_UPPER;
	else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, lp->dz[j]) &&
	    mpf_EGlpNumIsNeqq (lp->lz[col], mpf_NINFTY))
	    lp->vstat[col] = STAT_LOWER;
    }
}

void mpf_ILLfct_dphaseI_simple_update (mpf_lpinfo * lp,
      mpf_t ftol)
{
    int j, col;

    for (j = 0; j < lp->nnbasic; j++) {
	if (mpf_EGlpNumIsEqual (lp->dz[j], mpf_zeroLpNum, ftol))
	    continue;
	col = lp->nbaz[j];
	if (mpf_EGlpNumIsLess (lp->dz[j], mpf_zeroLpNum) && lp->vtype[col] == VBOUNDED)
	    lp->vstat[col] = STAT_UPPER;
	else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, lp->dz[j]) && lp->vtype[col] == VBOUNDED)
	    lp->vstat[col] = STAT_LOWER;
    }
}

/* set status values */
void mpf_ILLfct_set_status_values (mpf_lpinfo * lp,
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

void mpf_ILLfct_init_counts (mpf_lpinfo * lp)
{
    int i;
    mpf_count_struct *c = lp->cnts;
#define mpf_C_VALUE(a) (1.0+(double)(a)/(PARAM_HEAP_RATIO*mpf_ILLutil_our_log2(a)))
    mpf_EGlpNumSet (c->y_ravg, mpf_C_VALUE (lp->nrows));
    mpf_EGlpNumSet (c->za_ravg, mpf_C_VALUE (lp->nnbasic));
    ILL_IFTRACE ("%s:%la\n", __func__, mpf_EGlpNumToLf (c->za_ravg));
#undef mpf_C_VALUE
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

static void mpf_update_piv_values (mpf_count_struct * c,
      int phase,
      mpf_t piv2)
{
    int i = 0;
    mpf_t v, piv;

    if (mpf_EGlpNumIsEqqual (piv2, mpf_zeroLpNum))
	return;
    mpf_EGlpNumInitVar (v);
    mpf_EGlpNumInitVar (piv);
    mpf_EGlpNumCopyAbs (piv, piv2);
    mpf_EGlpNumOne (v);
    while (mpf_EGlpNumIsLess (piv, v) && i < 9) {
	mpf_EGlpNumDivUiTo (v, 10);
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
    mpf_EGlpNumClearVar (v);
    mpf_EGlpNumClearVar (piv);
}

void mpf_ILLfct_update_counts (mpf_lpinfo * lp,
      int f,
      int upi,
      mpf_t upd)
{
    mpf_count_struct *c = lp->cnts;

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
	mpf_update_piv_values (c, PRIMAL_PHASEI, upd);
	break;
    case CNT_PIIPIV:
	mpf_update_piv_values (c, PRIMAL_PHASEII, upd);
	break;
    case CNT_DIPIV:
	mpf_update_piv_values (c, DUAL_PHASEI, upd);
	break;
    case CNT_DIIPIV:
	mpf_update_piv_values (c, DUAL_PHASEII, upd);
	break;
    case CNT_YRAVG:
	mpf_EGlpNumMultUiTo (c->y_ravg, c->tot_iter);
	mpf_EGlpNumAddUiTo (c->y_ravg, upi);
	mpf_EGlpNumDivUiTo (c->y_ravg, c->tot_iter + 1);
	break;
    case CNT_ZARAVG:
	ILL_IFTRACE ("%s:%d:%d:%d:%la:%la", __func__, f, c->tot_iter, upi,
	    mpf_EGlpNumToLf (upd), mpf_EGlpNumToLf (c->za_ravg));
	mpf_EGlpNumMultUiTo (c->za_ravg, c->tot_iter);
	mpf_EGlpNumAddUiTo (c->za_ravg, upi);
	mpf_EGlpNumDivUiTo (c->za_ravg, c->tot_iter + 1);
	ILL_IFTRACE (":%la\n", mpf_EGlpNumToLf (c->za_ravg));
	break;
    }
}

void mpf_ILLfct_print_counts (mpf_lpinfo * lp)
{
    int i, niter;
    mpf_count_struct *c = lp->cnts;

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
static void mpf_add_vectors (mpf_lpinfo * lp,
      mpf_svector * a,
      mpf_svector * b,
      mpf_svector * c,
      mpf_t t)
{
    int i, r, l;
    mpf_svector *w = &(lp->work);

    for (i = 0; i < b->nzcnt; i++) {
	r = b->indx[i];
	w->indx[i] = r;
	mpf_EGlpNumCopy (w->coef[r], t);
	mpf_EGlpNumMultTo (w->coef[r], b->coef[i]);
	lp->iwork[r] = 1;
    }
    l = b->nzcnt;

    for (i = 0; i < a->nzcnt; i++) {
	r = a->indx[i];
	if (lp->iwork[r] == 0)
	    w->indx[l++] = r;
	mpf_EGlpNumAddTo (w->coef[r], a->coef[i]);
    }
    for (i = 0; i < l; i++) {
	r = w->indx[i];
	c->indx[i] = r;
	mpf_EGlpNumCopy (c->coef[i], w->coef[r]);
	mpf_EGlpNumZero (w->coef[r]);
	lp->iwork[r] = 0;
    }
    w->nzcnt = 0;
    c->nzcnt = l;
}

void mpf_ILLfct_update_pfeas (mpf_lpinfo * lp,
      int lindex,
      mpf_svector * srhs)
{
    int i, k, r;
    int col, nz = 0;
    int cbnd, f;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    mpf_t *t = lp->upd.t;
    mpf_t tz, *dty, ntmp;
    mpf_t *l, *x, *u, *pftol = &(lp->tol->ip_tol);
    mpf_EGlpNumInitVar (tz);
    mpf_EGlpNumInitVar (ntmp);
    dty = &(lp->upd.dty);
    mpf_EGlpNumZero (*dty);
    mpf_EGlpNumCopyAbs (tz, lp->upd.tz);
    mpf_EGlpNumDivUiTo (tz, 100);
    mpf_EGlpNumAddTo (tz, lp->upd.tz);
    ILL_IFTRACE ("%s:%d", __func__, tctr);
    for (i = 0; i < tctr && mpf_EGlpNumIsLeq (t[perm[i]], tz); i++) {
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
		mpf_EGlpNumCopyDiff (ntmp, *l, *x);
		if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY) && mpf_EGlpNumIsLess (*pftol, ntmp))
		    f = -1;
		else {
		    mpf_EGlpNumCopyDiff (ntmp, *x, *u);
		    if (mpf_EGlpNumIsNeqq (*u, mpf_INFTY) && mpf_EGlpNumIsLess (*pftol, ntmp))
			f = 1;
		}

		ILL_IFTRACE (":%d:%d", f, lp->bfeas[r]);
		if (f != lp->bfeas[r]) {
		    srhs->indx[nz] = r;
		    mpf_EGlpNumSet (srhs->coef[nz], (double) (f - lp->bfeas[r]));
		    mpf_EGlpNumAddInnProdTo (*dty, srhs->coef[nz], lp->yjz.coef[k]);
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
    mpf_EGlpNumClearVar (tz);
    mpf_EGlpNumClearVar (ntmp);
}

void mpf_ILLfct_compute_ppIzz (mpf_lpinfo * lp,
      mpf_svector * srhs,
      mpf_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	ILL_IFTRACE ("%s:\n", __func__);
	mpf_ILLbasis_row_solve (lp, srhs, ssoln);
    }
}

void mpf_ILLfct_update_ppI_prices (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * srhs,
      mpf_svector * ssoln,
      int eindex,
      int lindex,
      mpf_t alpha)
{
    mpf_t ntmp;
    mpf_EGlpNumInitVar (ntmp);
    mpf_EGlpNumCopy (ntmp, alpha);
    ILL_IFTRACE ("%s:\n", __func__);
    if (lindex == -1) {
	if (srhs->nzcnt != 0) {
	    mpf_ILLfct_update_pIpiz (lp, ssoln, mpf_oneLpNum);
	    if (pinf->p_strategy == COMPLETE_PRICING) {
		mpf_ILLfct_compute_zA (lp, ssoln, &(lp->zA));
		mpf_ILLfct_update_pIdz (lp, &(lp->zA), -1, mpf_oneLpNum);
	    }
	} else {
	    if (pinf->p_strategy == COMPLETE_PRICING)
		mpf_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
	    else
		mpf_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
	    mpf_EGlpNumClearVar (ntmp);
	    return;
	}
    } else {
	if (srhs->nzcnt == 0) {
	    mpf_ILLfct_update_pIpiz (lp, &(lp->zz), ntmp);
	    if (pinf->p_strategy == COMPLETE_PRICING)
		mpf_ILLfct_update_pIdz (lp, &(lp->zA), eindex, ntmp);
	} else {
	    mpf_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
	    mpf_EGlpNumSubTo (ntmp, alpha);
	    mpf_EGlpNumSign (ntmp);
	    mpf_add_vectors (lp, ssoln, &(lp->zz), &(lp->zz), ntmp);
	    mpf_ILLfct_update_pIpiz (lp, &(lp->zz), mpf_oneLpNum);
	    if (pinf->p_strategy == COMPLETE_PRICING) {
		mpf_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
		mpf_ILLfct_update_pIdz (lp, &(lp->zA), eindex, mpf_oneLpNum);
	    }
	}
	mpf_EGlpNumSet (lp->pIdz[eindex], (double) (lp->upd.fs));
	mpf_EGlpNumAddTo (lp->pIdz[eindex], ntmp);
	mpf_EGlpNumSign (lp->pIdz[eindex]);
    }
    if (pinf->p_strategy == COMPLETE_PRICING) {
	mpf_ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
	    PRIMAL_PHASEI);
	if (eindex > -1)
	    mpf_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
	mpf_ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, mpf_zeroLpNum);
    } else
	mpf_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
    mpf_EGlpNumClearVar (ntmp);
    return;
}

void mpf_ILLfct_update_dfeas (mpf_lpinfo * lp,
      int eindex,
      mpf_svector * srhs)
{
    int i, j, k, c;
    int cbnd, col, nz = 0;
    int vs, vt, f;
    int delta;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    int mcnt, mbeg;
    mpf_t *t = lp->upd.t;
    mpf_t *w = lp->work.coef;
    mpf_t tz;
    mpf_t *dty = &(lp->upd.dty);
    mpf_t *dftol = &(lp->tol->id_tol);
    mpf_t dj;
    mpf_EGlpNumInitVar (dj);
    mpf_EGlpNumInitVar (tz);
    mpf_EGlpNumZero (*dty);
    mpf_EGlpNumCopy (tz, lp->upd.tz);
    mpf_EGlpNumMultUiTo (tz, 101);
    mpf_EGlpNumDivUiTo (tz, 100);

    for (j = 0; j < tctr && mpf_EGlpNumIsLeq (t[perm[j]], tz); j++) {
	k = ix[perm[j]] / 10;
	c = lp->zA.indx[k];

	if (lp->iwork[c] != 1) {
	    lp->iwork[c] = 1;
	    cbnd = ix[perm[j]] % 10;
	    col = lp->nbaz[c];
	    mpf_EGlpNumCopy (dj, lp->dz[c]);
	    vs = lp->vstat[col];
	    vt = lp->vtype[col];

	    if (cbnd == BSKIP) {
		if (mpf_EGlpNumIsEqual (dj, mpf_zeroLpNum, *dftol));
		else if (mpf_EGlpNumIsLess (dj, mpf_zeroLpNum) && vs == STAT_LOWER)
		    lp->vstat[col] = STAT_UPPER;
		else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, dj) && vs == STAT_UPPER)
		    lp->vstat[col] = STAT_LOWER;
	    } else if (c != eindex) {
		if (mpf_EGlpNumIsEqual (dj, mpf_zeroLpNum, *dftol))
		    f = 0;
		else if (mpf_EGlpNumIsLess (dj, mpf_zeroLpNum) &&
		    (vs == STAT_LOWER || vs == STAT_ZERO))
		    f = -1;
		else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, dj) &&
		    (vs == STAT_UPPER || vs == STAT_ZERO))
		    f = 1;
		else
		    f = 0;

		if (f != lp->dfeas[c]) {
		    delta = f - lp->dfeas[c];
		    mcnt = lp->matcnt[col];
		    mbeg = lp->matbeg[col];
		    mpf_EGlpNumSet (dj, (double) (delta));
		    for (i = 0; i < mcnt; i++)
			mpf_EGlpNumAddInnProdTo (w[lp->matind[mbeg + i]], dj,
			    lp->matval[mbeg + i]);
		    mpf_EGlpNumAddInnProdTo (*dty, dj, lp->zA.coef[k]);
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
	    if (mpf_EGlpNumIsNeqqZero (w[i])) {
		mpf_EGlpNumCopy (srhs->coef[nz], w[i]);
		srhs->indx[nz] = i;
		nz++;
		mpf_EGlpNumZero (w[i]);
	    }
    }
    srhs->nzcnt = nz;
    mpf_EGlpNumClearVar (dj);
    mpf_EGlpNumClearVar (tz);
}

void mpf_ILLfct_compute_dpIy (mpf_lpinfo * lp,
      mpf_svector * srhs,
      mpf_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	mpf_ILLbasis_column_solve (lp, srhs, ssoln);
    }
}

void mpf_ILLfct_update_dpI_prices (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * srhs,
      mpf_svector * ssoln,
      int lindex,
      mpf_t alpha)
{
    int i;
    mpf_t ntmp;
    mpf_EGlpNumInitVar (ntmp);
    mpf_EGlpNumZero (ntmp);

    if (srhs->nzcnt == 0) {
	mpf_ILLfct_update_xz (lp, alpha, -1, -1);
    } else {
	mpf_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
	mpf_EGlpNumAddTo (ntmp, alpha);
	mpf_EGlpNumSign (ntmp);
	mpf_add_vectors (lp, ssoln, &(lp->yjz), &(lp->yjz), ntmp);
	mpf_EGlpNumSign (ntmp);
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    mpf_EGlpNumAddTo (lp->xbz[lp->yjz.indx[i]], lp->yjz.coef[i]);
    }
    mpf_EGlpNumSet (lp->xbz[lindex], ((double) (-lp->upd.fs)));
    mpf_EGlpNumAddTo (lp->xbz[lindex], ntmp);

    if (pinf->d_strategy == COMPLETE_PRICING) {
	mpf_ILLprice_compute_primal_inf (lp, pinf, lp->yjz.indx, lp->yjz.nzcnt,
	    DUAL_PHASEI);
	mpf_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEI);
	mpf_ILLfct_update_counts (lp, CNT_YRAVG, lp->yjz.nzcnt, mpf_zeroLpNum);
    } else
	mpf_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEI, ROW_PRICING);
    mpf_EGlpNumClearVar (ntmp);
}

void mpf_ILLfct_update_dIIfeas (mpf_lpinfo * lp,
      int eindex,
      mpf_svector * srhs)
{
    int j, k;
    int col, indx, vs;
    int *perm = lp->upd.perm;
    int *ix = lp->upd.ix;
    int tctr = lp->upd.tctr;
    mpf_t *zAj, *l, *u;
    mpf_t *dty = &(lp->upd.dty);
    mpf_t *t_max = &(lp->upd.tz);
    mpf_t *t = lp->upd.t;
    mpf_t delta;
    mpf_svector a;
    mpf_EGlpNumInitVar (delta);
    mpf_EGlpNumZero (delta);
    mpf_EGlpNumZero (*dty);

    srhs->nzcnt = 0;
    for (j = 0; j < tctr && mpf_EGlpNumIsLeq (t[perm[j]], *t_max); j++) {
	k = ix[perm[j]];
	indx = lp->zA.indx[k];

	if (indx != eindex) {
	    zAj = &(lp->zA.coef[k]);
	    col = lp->nbaz[indx];
	    l = &(lp->lz[col]);
	    u = &(lp->uz[col]);
	    vs = lp->vstat[col];
	    if (vs == STAT_UPPER)
		mpf_EGlpNumCopyDiff (delta, *l, *u);
	    else
		mpf_EGlpNumCopyDiff (delta, *u, *l);
	    mpf_EGlpNumAddInnProdTo (*dty, delta, *zAj);
	    lp->vstat[col] = (vs == STAT_UPPER) ? STAT_LOWER : STAT_UPPER;

	    a.nzcnt = lp->matcnt[col];
	    a.indx = &(lp->matind[lp->matbeg[col]]);
	    a.coef = &(lp->matval[lp->matbeg[col]]);
	    mpf_add_vectors (lp, srhs, &a, srhs, delta);
	}
    }
    mpf_EGlpNumClearVar (delta);
}

void mpf_ILLfct_compute_dpIIy (mpf_lpinfo * lp,
      mpf_svector * srhs,
      mpf_svector * ssoln)
{
    if (srhs->nzcnt != 0) {
	mpf_ILLbasis_column_solve (lp, srhs, ssoln);
    }
}

void mpf_ILLfct_update_dpII_prices (mpf_lpinfo * lp, mpf_price_info * pinf,
      mpf_svector * srhs, mpf_svector * ssoln, int lindex, mpf_t eval,
      mpf_t alpha)
{
    int i;
    mpf_svector *u;

    if (srhs->nzcnt == 0) {
	mpf_ILLfct_update_xz (lp, alpha, -1, -1);
	u = &(lp->yjz);
    } else {
	if (ssoln->nzcnt != 0)
	    for (i = 0; i < ssoln->nzcnt; i++)
		mpf_EGlpNumSubTo (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
	mpf_ILLfct_update_xz (lp, alpha, -1, -1);
	mpf_add_vectors (lp, ssoln, &(lp->yjz), ssoln, mpf_oneLpNum);
	u = ssoln;
    }
    mpf_EGlpNumCopySum (lp->xbz[lindex], eval, alpha);

    if (pinf->d_strategy == COMPLETE_PRICING) {
	mpf_ILLprice_compute_primal_inf (lp, pinf, u->indx, u->nzcnt, DUAL_PHASEII);
	mpf_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEII);
	mpf_ILLfct_update_counts (lp, CNT_YRAVG, u->nzcnt, mpf_zeroLpNum);
    } else
	mpf_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEII, ROW_PRICING);
}

int mpf_ILLfct_test_pivot (mpf_lpinfo * lp,
      int indx,
      int indxtype,
      mpf_t piv_val)
{
    int i;
    mpf_t pval, ntmp;
    mpf_EGlpNumInitVar (pval);
    mpf_EGlpNumInitVar (ntmp);
    mpf_EGlpNumZero (pval);

    if (indxtype == ROW_PIVOT) {
	for (i = 0; i < lp->yjz.nzcnt; i++)
	    if (lp->yjz.indx[i] == indx) {
		mpf_EGlpNumCopy (pval, lp->yjz.coef[i]);
		break;
	    }
    } else {
	for (i = 0; i < lp->zA.nzcnt; i++)
	    if (lp->zA.indx[i] == indx) {
		mpf_EGlpNumCopy (pval, lp->zA.coef[i]);
		break;
	    }
    }
    mpf_EGlpNumCopyDiff (ntmp, pval, piv_val);
    mpf_EGlpNumDivTo (ntmp, piv_val);
    if (mpf_EGlpNumIsLess (ntmp, mpf_zeroLpNum))
	mpf_EGlpNumSign (ntmp);
    if (mpf_EGlpNumIsLess (mpf_ALTPIV_TOLER, ntmp)) {
#if mpf_FCT_DEBUG > 1
	if (indxtype == ROW_PIVOT)
	    printf ("y_i = %.8f, z_j = %.8f %la %la\n", mpf_EGlpNumToLf (pval),
		mpf_EGlpNumToLf (piv_val), mpf_EGlpNumToLf (mpf_ALTPIV_TOLER),
		mpf_EGlpNumToLf (ntmp));
	else
	    printf ("z_j = %.8f, y_i = %.8f\n", mpf_EGlpNumToLf (pval),
		mpf_EGlpNumToLf (piv_val));
#endif
	mpf_EGlpNumClearVar (ntmp);
	mpf_EGlpNumClearVar (pval);
	return 1;
    }
    mpf_EGlpNumClearVar (pval);
    mpf_EGlpNumClearVar (ntmp);
    return 0;
}

#if mpf_FCT_DEBUG > 0

void mpf_fct_test_workvector (mpf_lpinfo * lp)
{
    int i, err = 0;
    for (i = 0; i < lp->ncols; i++) {
	if (mpf_EGlpNumIsNeqqZero (lp->work.coef[i])) {
	    err++;
	    mpf_EGlpNumZero (lp->work.coef[i]);
	}
	if (lp->iwork[i] != 0) {
	    err++;
	    lp->iwork[i] = 0;
	}
    }
    if (err)
	printf ("bad work vector, err=%d\n", err);
}

void mpf_fct_test_pfeasible (mpf_lpinfo * lp)
{
    int i, col;
    int err = 0;
    mpf_t *ftol = &(lp->tol->pfeas_tol);

    for (i = 0; i < lp->nrows; i++) {
	col = lp->baz[i];

	if (mpf_EGlpNumIsNeqq (lp->uz[col], mpf_INFTY)
	    && mpf_EGlpNumIsSumLess (*ftol, lp->uz[col], lp->xbz[i])) {
	    if (lp->bfeas[i] != 1) {
		err++;
		lp->bfeas[i] = 1;
	    }
	} else if (mpf_EGlpNumIsNeqq (lp->lz[col], mpf_NINFTY)
	    && mpf_EGlpNumIsSumLess (lp->xbz[i], *ftol, lp->lz[col])) {
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

void mpf_fct_test_dfeasible (mpf_lpinfo * lp)
{
    int j, col;
    int err = 0;
    mpf_t *ftol = &(lp->tol->dfeas_tol);
    mpf_t mftol[1];
    mpf_EGlpNumInitVar (mftol[0]);
    mpf_EGlpNumCopyNeg (mftol[0], *ftol);

    for (j = 0; j < lp->nnbasic; j++) {
	col = lp->nbaz[j];

	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	    continue;
	if (mpf_EGlpNumIsLess (lp->dz[j], mftol[0]) &&
	    (lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO)) {
	    if (lp->dfeas[j] != -1) {
		err++;
		lp->dfeas[j] = -1;
	    }
	}
	if (mpf_EGlpNumIsLess (*ftol, lp->dz[j]) &&
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

void mpf_fct_test_pI_x (mpf_lpinfo * lp,
      mpf_price_info * p)
{
    int i;
    int ern = 0;
    mpf_t *x;
    mpf_t err, diff;
    mpf_EGlpNumInitVar (err);
    mpf_EGlpNumInitVar (diff);
    mpf_EGlpNumZero (err);
    x = mpf_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumCopy (x[i], lp->xbz[i]);
    mpf_ILLfct_compute_phaseI_xbz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_PFEAS_TOLER, diff)) {
	    mpf_EGlpNumAddTo (err, diff);
	    ern++;
	    printf ("bad i = %d\n", i);
	}
    }
    if (mpf_EGlpNumIsNeqqZero (err))
	printf ("dI x err = %.7f, ern = %d\n", mpf_EGlpNumToLf (err), ern);
    mpf_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
    mpf_EGlpNumFreeArray (x);
    mpf_EGlpNumClearVar (diff);
    mpf_EGlpNumClearVar (err);
}

void mpf_fct_test_pII_x (mpf_lpinfo * lp,
      mpf_price_info * p)
{
    int i;
    int ern = 0;
    mpf_t *x;
    mpf_t err, diff;
    mpf_EGlpNumInitVar (err);
    mpf_EGlpNumInitVar (diff);
    mpf_EGlpNumZero (err);
    x = mpf_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumCopy (x[i], lp->xbz[i]);
    mpf_ILLfct_compute_xbz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_PFEAS_TOLER, diff)) {
	    mpf_EGlpNumAddTo (err, diff);
	    ern++;
	    printf ("bad i = %d\n", i);
	}
    }
    if (mpf_EGlpNumIsNeqqZero (err))
	printf ("dII x err = %.7f, ern = %d\n", mpf_EGlpNumToLf (err), ern);
    mpf_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
    mpf_EGlpNumFreeArray (x);
    mpf_EGlpNumClearVar (diff);
    mpf_EGlpNumClearVar (err);
}

void mpf_fct_test_pI_pi_dz (mpf_lpinfo * lp,
      mpf_price_info * p)
{
    int i;
    int ern = 0;
    mpf_t *pidz;
    mpf_t err, diff;
    mpf_EGlpNumInitVar (err);
    mpf_EGlpNumInitVar (diff);
    pidz = mpf_EGlpNumAllocArray (lp->ncols);
    mpf_EGlpNumZero (err);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumCopy (pidz[i], lp->pIpiz[i]);
    mpf_ILLfct_compute_phaseI_piz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumCopyDiff (diff, pidz[i], lp->pIpiz[i]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_DFEAS_TOLER, diff)) {
	    mpf_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpf_EGlpNumIsNeqqZero (err))
	printf ("pI pi err = %.7f, ern = %d\n", mpf_EGlpNumToLf (err), ern);

    mpf_EGlpNumZero (err);
    ern = 0;
    for (i = 0; i < lp->nnbasic; i++)
	mpf_EGlpNumCopy (pidz[i], lp->pIdz[i]);
    mpf_ILLfct_compute_phaseI_dz (lp);
    for (i = 0; i < lp->nnbasic; i++) {
	mpf_EGlpNumCopyDiff (diff, pidz[i], lp->pIdz[i]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_DFEAS_TOLER, diff)) {
	    mpf_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpf_EGlpNumIsNeqqZero (err))
	printf ("pI dz err = %.7f, ern = %d\n", mpf_EGlpNumToLf (err), ern);
    mpf_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
    mpf_EGlpNumClearVar (err);
    mpf_EGlpNumClearVar (diff);
    mpf_EGlpNumFreeArray (pidz);
}

void mpf_fct_test_pII_pi_dz (mpf_lpinfo * lp,
      mpf_price_info * p)
{
    int i;
    int ern = 0;
    mpf_t *pidz;
    mpf_t err, diff;
    mpf_EGlpNumInitVar (err);
    mpf_EGlpNumInitVar (diff);
    mpf_EGlpNumZero (err);
    pidz = mpf_EGlpNumAllocArray (lp->ncols);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumCopy (pidz[i], lp->piz[i]);
    mpf_ILLfct_compute_piz (lp);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumCopyDiff (diff, pidz[i], lp->piz[i]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_DFEAS_TOLER, diff)) {
	    mpf_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpf_EGlpNumIsNeqqZero (err))
	printf ("pII pi err = %.7f, ern = %d\n", mpf_EGlpNumToLf (err), ern);

    mpf_EGlpNumZero (err);
    ern = 0;
    for (i = 0; i < lp->nnbasic; i++)
	mpf_EGlpNumCopy (pidz[i], lp->dz[i]);
    mpf_ILLfct_compute_dz (lp);
    for (i = 0; i < lp->nnbasic; i++) {
	mpf_EGlpNumCopyDiff (diff, pidz[i], lp->dz[i]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_DFEAS_TOLER, diff)) {
	    mpf_EGlpNumAddTo (err, diff);
	    ern++;
	}
    }
    if (mpf_EGlpNumIsNeqqZero (err))
	printf ("pII dz err = %.7f, ern = %d\n", mpf_EGlpNumToLf (err), ern);
    /*
     * mpf_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
     */
    mpf_EGlpNumClearVar (err);
    mpf_EGlpNumClearVar (diff);
    mpf_EGlpNumFreeArray (pidz);
}

#endif
