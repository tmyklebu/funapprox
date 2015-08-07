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

/* RCS_INFO = "$RCSfile: mpf_basis.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "mpf_sortrus.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdefs.h"
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#include "mpf_basis.h"
#include "mpf_fct.h"
#include "mpf_lp.h"
#include "mpf_lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define mpf_DJZERO_TOLER mpf_PFEAS_TOLER
#define mpf_BASIS_STATS 0
#define mpf_BASIS_DEBUG 0

void mpf_ILLbasis_init_vardata (mpf_var_data * vd)
{
    memset (vd, 0, sizeof (mpf_var_data));
    mpf_EGlpNumInitVar (vd->cmax);
}

void mpf_ILLbasis_clear_vardata (mpf_var_data * vd)
{
    mpf_EGlpNumClearVar (vd->cmax);
    memset (vd, 0, sizeof (mpf_var_data));
}

static void mpf_get_var_info (mpf_lpinfo * lp,
      mpf_var_data * v);

static int mpf_init_slack_basis (mpf_lpinfo * lp,
      int *vstat,
      int *irow,
      int *rrow,
      int *unitcol,
      int *icol,
      int *rcol),
  mpf_get_initial_basis1 (mpf_lpinfo * lp,
      int *vstat),
  mpf_get_initial_basis2 (mpf_lpinfo * lp,
      int *vstat),
  mpf_set_basis_indices (mpf_lpinfo * lp,
      int *vstat),
  mpf_choose_basis (int algorithm,
      mpf_t pinf1,
      mpf_t dinf1,
      mpf_t pinf2,
      mpf_t dinf2);

void mpf_ILLbasis_init_basisinfo (mpf_lpinfo * lp)
{
    lp->baz = 0;
    lp->nbaz = 0;
    lp->vstat = 0;
    lp->vindex = 0;
    lp->f = 0;
}

void mpf_ILLbasis_free_basisinfo (mpf_lpinfo * lp)
{
    ILL_IFFREE (lp->baz, int);
    ILL_IFFREE (lp->nbaz, int);
    ILL_IFFREE (lp->vstat, int);
    ILL_IFFREE (lp->vindex, int);
    if (lp->f) {
	mpf_ILLfactor_free_factor_work (lp->f);
	mpf_EGlpNumClearVar (lp->f->fzero_tol);
	mpf_EGlpNumClearVar (lp->f->szero_tol);
	mpf_EGlpNumClearVar (lp->f->partial_tol);
	mpf_EGlpNumClearVar (lp->f->maxelem_orig);
	mpf_EGlpNumClearVar (lp->f->maxelem_factor);
	mpf_EGlpNumClearVar (lp->f->maxelem_cur);
	mpf_EGlpNumClearVar (lp->f->partial_cur);
	ILL_IFFREE (lp->f, mpf_factor_work);
    }
}

int mpf_ILLbasis_build_basisinfo (mpf_lpinfo * lp)
{
    int rval = 0;

    ILL_SAFE_MALLOC (lp->baz, lp->nrows, int);
    ILL_SAFE_MALLOC (lp->nbaz, lp->nnbasic, int);
    ILL_SAFE_MALLOC (lp->vstat, lp->ncols, int);
    ILL_SAFE_MALLOC (lp->vindex, lp->ncols, int);

    lp->fbasisid = -1;

CLEANUP:
    if (rval)
	mpf_ILLbasis_free_basisinfo (lp);
    ILL_RETURN (rval, "mpf_ILLbasis_build_basisinfo");
}

int mpf_ILLbasis_load (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B)
{
    int rval = 0;
    char *cstat = B->cstat;
    char *rstat = B->rstat;
    int *structmap = lp->O->structmap;
    int *rowmap = lp->O->rowmap;
    mpf_t *rng = lp->O->rangeval;
    int i, j, ncols = lp->ncols, nrows = lp->nrows, nstruct = lp->O->nstruct;
    int basic = 0, nonbasic = 0;

    mpf_ILLbasis_free_basisinfo (lp);
    mpf_ILLbasis_init_basisinfo (lp);
    rval = mpf_ILLbasis_build_basisinfo (lp);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < nstruct; i++) {
	j = structmap[i];
	if (cstat[i] == QS_COL_BSTAT_BASIC) {
	    lp->vstat[j] = STAT_BASIC;
	    lp->baz[basic] = j;
	    lp->vindex[j] = basic;
	    basic++;
	} else {
	    lp->nbaz[nonbasic] = j;
	    lp->vindex[j] = nonbasic;
	    nonbasic++;
	    switch (cstat[i]) {
	    case QS_COL_BSTAT_LOWER:
		lp->vstat[j] = STAT_LOWER;
		break;
	    case QS_COL_BSTAT_UPPER:
		lp->vstat[j] = STAT_UPPER;
		break;
	    case QS_COL_BSTAT_FREE:
		lp->vstat[j] = STAT_ZERO;
		break;
	    default:
		fprintf (stderr, "unknown col basis stat 1: %c\n", cstat[i]);
		rval = 1;
		goto CLEANUP;
	    }
	}
    }

    for (i = 0; i < nrows; i++) {
	j = rowmap[i];
	if (rng && (mpf_EGlpNumIsNeqqZero (rng[i]))) {
	    if (rstat[i] == QS_ROW_BSTAT_BASIC) {
		lp->vstat[j] = STAT_BASIC;
		lp->baz[basic] = j;
		lp->vindex[j] = basic;
		basic++;
	    } else {
		lp->nbaz[nonbasic] = j;
		lp->vindex[j] = nonbasic;
		nonbasic++;
		switch (rstat[i]) {
		case QS_ROW_BSTAT_LOWER:
		    lp->vstat[j] = STAT_LOWER;
		    break;
		case QS_ROW_BSTAT_UPPER:
		    lp->vstat[j] = STAT_UPPER;
		    break;
		default:
		    fprintf (stderr, "unknown range basis stat 2\n");
		    rval = 1;
		    goto CLEANUP;
		}
	    }
	} else {
	    switch (rstat[i]) {
	    case QS_ROW_BSTAT_BASIC:
		lp->vstat[j] = STAT_BASIC;
		lp->baz[basic] = j;
		lp->vindex[j] = basic;
		basic++;
		break;
	    case QS_ROW_BSTAT_LOWER:
		lp->vstat[j] = STAT_LOWER;
		lp->nbaz[nonbasic] = j;
		lp->vindex[j] = nonbasic;
		nonbasic++;
		break;
	    default:
		fprintf (stderr, "unknown row basis stat 3\n");
		rval = 1;
		goto CLEANUP;
	    }
	}
    }

    if (basic + nonbasic != ncols) {
	fprintf (stderr, "error in counts in ILLopt_load_basis\n");
	rval = 1;
	goto CLEANUP;
    }
    if (lp->fbasisid != 0)
	lp->basisid = 0;
    else
	lp->basisid = 1;

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLbasis_load");
}

int mpf_ILLbasis_tableau_row (mpf_lpinfo * lp,
      int row,
      mpf_t * brow,
      mpf_t * trow,
      mpf_t * rhs,
      int strict)
{
    int rval = 0;
    int i;
    int singular = 0;
    int indx;
    mpf_t coef;
    mpf_t sum;
    mpf_svector z, zA;
    mpf_EGlpNumInitVar (coef);
    mpf_EGlpNumInitVar (sum);
    mpf_EGlpNumZero (sum);

    mpf_ILLsvector_init (&z);
    mpf_ILLsvector_init (&zA);

    if (lp->basisid == -1) {
	fprintf (stderr, "mpf_ILLbasis_tableau_row: no basis\n");
	rval = E_GENERAL_ERROR;
	ILL_CLEANUP;
    }
    if (lp->fbasisid != lp->basisid) {	/* Needs to be changed */
	rval = mpf_ILLbasis_factor (lp, &singular);
	ILL_CLEANUP_IF (rval);
	if (singular) {
	    rval = E_BASIS_SINGULAR;
	    ILL_CLEANUP;
	}
    }
    if (brow == NULL) {
	fprintf (stderr, "No array for basis inverse row\n");
	rval = E_GENERAL_ERROR;
	ILL_CLEANUP;
    }
    rval = mpf_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);
    mpf_ILLfct_compute_zz (lp, &z, row);

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumZero (brow[i]);
    for (i = 0; i < z.nzcnt; i++) {
	indx = z.indx[i];
	mpf_EGlpNumCopy (coef, z.coef[i]);
	mpf_EGlpNumCopy (brow[indx], coef);
	mpf_EGlpNumAddInnProdTo (sum, coef, lp->bz[indx]);
    }

    if (rhs != NULL)
	mpf_EGlpNumCopy (*rhs, sum);
    if (trow != NULL) {
	if (!strict) {
	    rval = mpf_ILLsvector_alloc (&zA, lp->ncols);
	    if (rval)
		ILL_CLEANUP;
	    ILL_IFTRACE ("%s:\n", __func__);
	    rval = mpf_ILLfct_compute_zA (lp, &z, &zA);
	    ILL_CLEANUP_IF (rval);

	    for (i = 0; i < lp->ncols; i++)
		mpf_EGlpNumZero (trow[i]);
	    for (i = 0; i < zA.nzcnt; i++)
		mpf_EGlpNumCopy (trow[lp->nbaz[zA.indx[i]]], zA.coef[i]);
	    mpf_EGlpNumOne (trow[lp->baz[row]]);
	} else {
	    mpf_ILLfct_compute_vA (lp, &z, trow);
	}
    }
#if mpf_BASIS_DEBUG > 0
    if (rhs != NULL && trow != NULL) {
	mpf_t *tr = NULL;
	mpf_EGlpNumZero (sum);
	if (strict)
	    tr = trow;
	else {
	    tr = mpf_EGlpNumAllocArray (lp->ncols);
	    mpf_ILLfct_compute_vA (lp, &z, tr);
	}
	for (i = 0; i < lp->nrows; i++)
	    if (mpf_EGlpNumIsLess (mpf_zeroLpNum, tr[lp->baz[i]]))
		mpf_EGlpNumAddTo (sum, tr[lp->baz[i]]);
	    else
		mpf_EGlpNumSubTo (sum, tr[lp->baz[i]]);
	mpf_EGlpNumCopy (coef, mpf_oneLpNum);
	mpf_EGlpNumSubTo (coef, sum);
	if (mpf_EGlpNumIsLess (coef, mpf_zeroLpNum))
	    mpf_EGlpNumSign (coef);
	if (mpf_EGlpNumIsLess (mpf_PIVZ_TOLER, coef))
	    fprintf (stderr, "tableau: bas computed = %.12f\n", mpf_EGlpNumToLf (sum));
	if (!strict)
	    mpf_EGlpNumFreeArray (tr);
#if mpf_BASIS_DEBUG > 1
	mpf_EGlpNumZero (sum);
	for (i = 0; i < lp->ncols; i++) {
	    if (lp->vstat[i] == STAT_BASIC)
		mpf_EGlpNumAddInnProdTo (sum, lp->xbz[lp->vindex[i]], trow[i]);
	    else if (lp->vstat[i] == STAT_UPPER)
		mpf_EGlpNumAddInnProdTo (sum, lp->uz[i], trow[i]);
	    else if (lp->vstat[i] == STAT_LOWER)
		mpf_EGlpNumAddInnProdTo (sum, lp->lz[i], trow[i]);
	}
	mpf_EGlpNumSet (coef, 1e-10);
	if (mpf_EGlpNumIsNeq (sum, *rhs, coef))
	    fprintf (stderr, "tableau rhs = %.9f, computed = %.9f\n",
		mpf_EGlpNumToLf (*rhs), mpf_EGlpNumToLf (sum));
#endif
    }
#endif

CLEANUP:
    mpf_ILLsvector_free (&z);
    mpf_ILLsvector_free (&zA);
    mpf_EGlpNumClearVar (coef);
    mpf_EGlpNumClearVar (sum);
    return rval;
}

static void mpf_get_var_info (mpf_lpinfo * lp,
      mpf_var_data * v)
{
    int i = 0;

    v->nartif = 0;
    v->nslacks = 0;
    v->nfree = 0;
    v->nbndone = 0;
    v->nbounded = 0;
    v->nfixed = 0;
    mpf_EGlpNumCopy (v->cmax, mpf_NINFTY);

    for (i = 0; i < lp->ncols; i++) {
	switch (lp->vtype[i]) {
	case VARTIFICIAL:
	    v->nartif++;
	    break;
	case VFREE:
	    v->nfree++;
	    break;
	case VLOWER:
	case VUPPER:
	    if (lp->vclass[i] == CLASS_LOGICAL)
		v->nslacks++;
	    else
		v->nbndone++;
	    break;

	case VFIXED:
	    v->nfixed++;
	case VBOUNDED:
	    if (lp->vclass[i] == CLASS_LOGICAL)
		v->nslacks++;
	    else
		v->nbounded++;
	    break;
	}
	mpf_EGlpNumSetToMaxAbs (v->cmax, lp->cz[i]);
    }

#if mpf_BASIS_STATS > 0
    printf ("cols = %d, acols = %d, total  = %d, nrows = %d, nlog = %d\n",
	lp->ncols, lp->ncols - lp->nrows,
	v->nartif + v->nfree + v->nslacks + v->nbndone + v->nbounded,
	lp->nrows, v->nartif + v->nslacks);
#endif
}

static int mpf_init_slack_basis (mpf_lpinfo * lp,
      int *vstat,
      int *irow,
      int *rrow,
      int *unitcol,
      int *icol,
      int *rcol)
{
    int j, r, vt;
    int nslacks = 0;

    for (j = 0; j < lp->ncols; j++) {
	r = lp->matind[lp->matbeg[j]];
	vt = lp->vtype[j];

	if ((vt == VUPPER || vt == VLOWER || vt == VBOUNDED || vt == VFIXED) &&
	    lp->vclass[j] == CLASS_LOGICAL) {

	    vstat[j] = STAT_BASIC;
	    irow[r] = 1;
	    rrow[r] = 1;
	    unitcol[r] = j;
	    if (icol != NULL) {
		icol[j] = 1;
		rcol[j] = 1;
	    }
	    nslacks++;
	} else if (vt == VARTIFICIAL) {
	    unitcol[r] = j;
	    vstat[j] = STAT_UPPER;
	} else if (vt == VFREE)
	    vstat[j] = STAT_ZERO;
	else if (vt == VFIXED || vt == VUPPER)
	    vstat[j] = STAT_UPPER;
	else if (vt == VLOWER)
	    vstat[j] = STAT_LOWER;
	else if (vt == VBOUNDED) {
	    if (fabs (mpf_EGlpNumToLf (lp->lz[j])) < fabs (mpf_EGlpNumToLf (lp->uz[j])))
		vstat[j] = STAT_LOWER;
	    else
		vstat[j] = STAT_UPPER;
	}
    }
    return nslacks;
}

static int mpf_primal_col_select (mpf_lpinfo * lp,
      int *vstat,
      int *irow,
      int *rrow,
      int *unitcol,
      mpf_t * v,
      int *perm,
      int *porder,
      int nbelem,
      int pcols)
{
    int i, j, k, tr, r = 0;
    int mcnt, mbeg;
    int *matbeg = lp->matbeg;
    int *matcnt = lp->matcnt;
    int *matind = lp->matind;
    mpf_t *matval = lp->matval;
    mpf_t alpha, val, maxelem;
    mpf_EGlpNumInitVar (alpha);
    mpf_EGlpNumInitVar (val);
    mpf_EGlpNumInitVar (maxelem);

    for (k = 0; k < pcols; k++) {
	j = porder[perm[k]];
	mcnt = matcnt[j];
	mbeg = matbeg[j];

	mpf_EGlpNumCopy (alpha, mpf_NINFTY);
	mpf_EGlpNumCopy (maxelem, mpf_NINFTY);

	for (i = 0; i < mcnt; i++) {
	    mpf_EGlpNumCopyAbs (val, matval[mbeg + i]);
	    if (mpf_EGlpNumIsLess (maxelem, val))
		mpf_EGlpNumCopy (maxelem, val);
	    if (rrow[matind[mbeg + i]] == 0 && mpf_EGlpNumIsLess (alpha, val)) {
		mpf_EGlpNumCopy (alpha, val);
		r = matind[mbeg + i];
	    }
	}
	mpf_EGlpNumCopy (val, maxelem);
	mpf_EGlpNumMultTo (val, mpf_PARAM_IBASIS_RPIVOT);
	if (mpf_EGlpNumIsLess (val, alpha)) {
	    vstat[j] = STAT_BASIC;
	    nbelem++;
	    irow[r] = 1;
	    mpf_EGlpNumCopy (v[r], alpha);
	    for (i = 0; i < mcnt; i++)
		if (mpf_EGlpNumIsNeqqZero (matval[mbeg + i]))
		    rrow[matind[mbeg + i]]++;
	} else {
	    mpf_EGlpNumCopy (alpha, mpf_NINFTY);
	    for (i = 0; i < mcnt; i++) {
		tr = matind[mbeg + i];
		mpf_EGlpNumCopyAbs (val, matval[mbeg + i]);
		mpf_EGlpNumDivTo (val, mpf_PARAM_IBASIS_RTRIANG);
		if (mpf_EGlpNumIsNeqq (v[tr], mpf_INFTY) && mpf_EGlpNumIsLess (v[tr], val)) {
		    mpf_EGlpNumZero (alpha);
		    break;
		}
		mpf_EGlpNumCopyAbs (val, matval[mbeg + i]);
		if (irow[tr] == 0 && mpf_EGlpNumIsLess (alpha, val)) {
		    mpf_EGlpNumCopy (alpha, val);
		    r = tr;
		}
	    }
	    if (mpf_EGlpNumIsNeqqZero (alpha) && mpf_EGlpNumIsNeqq (alpha, mpf_NINFTY)) {
		vstat[j] = STAT_BASIC;
		nbelem++;
		irow[r] = 1;
		mpf_EGlpNumCopy (v[r], alpha);
		for (i = 0; i < mcnt; i++)
		    if (mpf_EGlpNumIsNeqqZero (matval[mbeg + i]))
			rrow[matind[mbeg + i]]++;
	    }
	}
    }
#if mpf_BASIS_STATS > 0
    printf ("nartifs = %d\n", lp->nrows - nbelem);
#endif

    if (nbelem < lp->nrows) {
	for (i = 0; i < lp->nrows; i++) {
	    if (irow[i] == 0) {
		if (unitcol[i] != -1) {
		    vstat[unitcol[i]] = STAT_BASIC;
		    nbelem++;
		} else {
		    fprintf (stderr, "Error: Not enough artificials\n");
		    return -1;
		}
	    }
	}
    }
    mpf_EGlpNumClearVar (alpha);
    mpf_EGlpNumClearVar (val);
    mpf_EGlpNumClearVar (maxelem);
    return nbelem;
}

/* This is an implementation of the initial basis procedure in: "Implementing
   the simplex method: the initial basis", by Bob Bixby. Goals: choose
   initial variables to go into basis which satisfy: 1) vars are slacks, 2)
   vars have freedom to move 3) initial submatrix is nonsingular, 4) low
   objective function contribution. */
static int mpf_get_initial_basis1 (mpf_lpinfo * lp,
      int *vstat)
{
    int rval = 0;
    int i, j, tot1 = 0, tot2 = 0;
    int nbelem = 0, nslacks = 0;
    int tfree = 0, tbndone = 0;
    int tbounded = 0;
    int *irow = NULL, *rrow = NULL;
    int *perm = NULL, *porder = NULL;
    int *unitcol = NULL;
    mpf_t cmax;
    mpf_t *v = NULL;
    mpf_t *qpenalty = NULL;
    mpf_var_data vd;
    mpf_ILLbasis_init_vardata (&vd);
    mpf_EGlpNumInitVar (cmax);

    mpf_get_var_info (lp, &vd);
    if (mpf_EGlpNumIsEqqual (vd.cmax, mpf_zeroLpNum))
	mpf_EGlpNumOne (cmax);
    else {
	mpf_EGlpNumCopy (cmax, vd.cmax);
	mpf_EGlpNumMultUiTo (cmax, 1000);
    }

    ILL_SAFE_MALLOC (irow, lp->nrows, int);
    ILL_SAFE_MALLOC (rrow, lp->nrows, int);
    v = mpf_EGlpNumAllocArray (lp->nrows);
    ILL_SAFE_MALLOC (unitcol, lp->nrows, int);

    for (i = 0; i < lp->nrows; i++) {
	unitcol[i] = -1;
	mpf_EGlpNumCopy (v[i], mpf_INFTY);
	irow[i] = 0;
	rrow[i] = 0;
    }

    nslacks = mpf_init_slack_basis (lp, vstat, irow, rrow, unitcol, NULL, NULL);
    if (nslacks != vd.nslacks) {
	printf ("complain: incorrect basis info(slacks)\n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
    if (nslacks == lp->nrows)
	ILL_CLEANUP;
    nbelem = nslacks;
    if (nbelem < lp->nrows) {
	for (i = 0; i < lp->nrows; i++) {
	    if (irow[i] == 0) {
		if (unitcol[i] != -1) {
		    vstat[unitcol[i]] = STAT_BASIC;
		    nbelem++;
		} else {
		    fprintf (stderr, "Error: Not enough artificials\n");
		    return -1;
		}
	    }
	}
    }
    ILL_CLEANUP;

    tot1 = vd.nfree + vd.nbndone;
    tot2 = vd.nfree + vd.nbndone + vd.nbounded;
    ILL_SAFE_MALLOC (perm, tot2, int);
    ILL_SAFE_MALLOC (porder, tot2, int);
    qpenalty = mpf_EGlpNumAllocArray (tot2);

    for (j = 0; j < lp->ncols; j++) {
	if (vstat[j] == STAT_BASIC)
	    continue;

	switch (lp->vtype[j]) {
	case VFREE:
	    porder[tfree] = j;
	    perm[tfree] = tfree;
	    mpf_EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
	    tfree++;
	    break;

	case VLOWER:
	case VUPPER:
	    porder[vd.nfree + tbndone] = j;
	    perm[vd.nfree + tbndone] = tbndone;
	    mpf_EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
	    if (lp->vtype[j] == VLOWER)
		mpf_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
	    else
		mpf_EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
	    tbndone++;
	    break;

	case VFIXED:
	case VBOUNDED:
	    porder[tot1 + tbounded] = j;
	    perm[tot1 + tbounded] = tbounded;
	    mpf_EGlpNumCopyFrac (qpenalty[tot1 + tbndone], lp->cz[j], cmax);
	    mpf_EGlpNumAddTo (qpenalty[tot1 + tbndone], lp->lz[j]);
	    mpf_EGlpNumSubTo (qpenalty[tot1 + tbndone], lp->uz[j]);
	    tbounded++;
	    break;
	}
    }
    if (tfree != vd.nfree || tbndone != vd.nbndone || tbounded != vd.nbounded) {
	printf ("complain: incorrect basis info \n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
    mpf_ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, vd.nfree);
    mpf_ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
	vd.nbndone);
    mpf_ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, vd.nbounded);

    for (i = 0; i < vd.nbndone; i++)
	perm[vd.nfree + i] += vd.nfree;
    for (i = 0; i < vd.nbounded; i++)
	perm[tot1 + i] += tot1;

    nbelem =
	mpf_primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
	tot2);
    if (nbelem != lp->nrows) {
	printf ("complain: incorrect final basis size\n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
CLEANUP:
    mpf_EGlpNumClearVar (cmax);
    if (rval)
	mpf_ILLbasis_free_basisinfo (lp);
    ILL_IFFREE (irow, int);
    ILL_IFFREE (rrow, int);
    mpf_EGlpNumFreeArray (v);
    ILL_IFFREE (perm, int);
    ILL_IFFREE (porder, int);
    ILL_IFFREE (unitcol, int);
    mpf_EGlpNumFreeArray (qpenalty);
    mpf_ILLbasis_clear_vardata (&vd);
    ILL_RETURN (rval, "mpf_ILLbasis_get_initial");
}

static int mpf_get_initial_basis2 (mpf_lpinfo * lp,
      int *vstat)
{
    int rval = 0;
    int i, j, k, tot1, tot2;
    int rbeg, rcnt, mcnt;
    int nbelem = 0, nslacks = 0;
    int tfree = 0, tbndone = 0;
    int tbounded = 0;
    int *irow = NULL, *rrow = NULL;
    int *perm = NULL, *porder = NULL;
    int *unitcol = NULL;
    mpf_t *v = NULL;
    mpf_t *qpenalty = NULL;
    int col = 0, s_i = 0, selc = 0;
    int *icol = NULL, *rcol = NULL;
    int *plen = NULL;
    mpf_t *dj = NULL;
    mpf_var_data vd;
    mpf_t seldj;
    mpf_t selv;
    mpf_t c_dj;
    mpf_t cmax;
    mpf_EGlpNumInitVar (seldj);
    mpf_EGlpNumInitVar (selv);
    mpf_EGlpNumInitVar (c_dj);
    mpf_EGlpNumInitVar (cmax);
    mpf_EGlpNumZero (c_dj);
    mpf_EGlpNumZero (selv);
    mpf_EGlpNumZero (seldj);
    mpf_ILLbasis_init_vardata (&vd);

    mpf_get_var_info (lp, &vd);

    ILL_SAFE_MALLOC (irow, lp->nrows, int);
    ILL_SAFE_MALLOC (rrow, lp->nrows, int);
    v = mpf_EGlpNumAllocArray (lp->nrows);
    ILL_SAFE_MALLOC (unitcol, lp->nrows, int);
    ILL_SAFE_MALLOC (icol, lp->ncols, int);
    ILL_SAFE_MALLOC (rcol, lp->ncols, int);
    dj = mpf_EGlpNumAllocArray (lp->ncols);

    for (i = 0; i < lp->nrows; i++) {
	unitcol[i] = -1;
	mpf_EGlpNumCopy (v[i], mpf_INFTY);
	irow[i] = 0;
	rrow[i] = 0;
    }
    /* assign all d_j */
    for (i = 0; i < lp->ncols; i++) {
	icol[i] = 0;
	rcol[i] = 0;
	mpf_EGlpNumCopy (dj[i], lp->cz[i]);
    }

    nslacks = mpf_init_slack_basis (lp, vstat, irow, rrow, unitcol, icol, rcol);
    if (nslacks != vd.nslacks) {
	printf ("complain: incorrect basis info\n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
    if (nslacks == lp->nrows)
	ILL_CLEANUP;
    nbelem = nslacks;

    /* allocate maximum required space for perm etc. */
    ILL_SAFE_MALLOC (perm, lp->ncols, int);
    ILL_SAFE_MALLOC (porder, lp->ncols, int);
    ILL_SAFE_MALLOC (plen, lp->nrows, int);
    qpenalty = mpf_EGlpNumAllocArray (lp->ncols);

    /* find all unit rows and record lengths */
    for (i = 0; i < lp->nrows; i++) {
	if (irow[i] != 1) {
	    rbeg = lp->rowbeg[i];
	    rcnt = lp->rowcnt[i];
	    for (j = 0; j < rcnt; j++) {
		mpf_EGlpNumCopyAbs (cmax, lp->rowval[rbeg + j]);
		if (mpf_EGlpNumIsNeqq (cmax, mpf_oneLpNum))
		    break;
	    }
	    if (j == rcnt) {
		perm[s_i] = s_i;
		porder[s_i] = i;
		plen[s_i] = rcnt;
		s_i++;
	    }
	}
    }

    /* sort all unit rows */
    mpf_ILLutil_int_perm_quicksort (perm, plen, s_i);

    /* now go through the unit rows */
    for (k = 0; k < s_i; k++) {
	i = porder[perm[k]];
	rbeg = lp->rowbeg[i];
	rcnt = lp->rowcnt[i];
	selc = -1;
	mpf_EGlpNumCopy (seldj, mpf_INFTY);
	mpf_EGlpNumZero (selv);

	/* for every row s_i, compute min {d_j : d_j <0 , j is u or l or fr} */
	for (j = 0; j < rcnt; j++) {
	    col = lp->rowind[rbeg + j];
	    if (rcol[col] == 1)
		break;
	    if (mpf_EGlpNumIsLess (dj[col], mpf_zeroLpNum)) {
		if (mpf_EGlpNumIsLess (dj[col], seldj)) {
		    selc = col;
		    mpf_EGlpNumCopy (seldj, dj[col]);
		    mpf_EGlpNumCopy (selv, lp->rowval[rbeg + j]);
		}
	    }
	}
	/* select pivot element and update all d_j's */
	if (selc != -1) {
	    nbelem++;
	    irow[i] = 1;
	    rrow[i] = 1;
	    icol[selc] = 1;
	    mpf_EGlpNumCopyFrac (c_dj, dj[selc], selv);
	    vstat[selc] = STAT_BASIC;
	    for (j = 0; j < rcnt; j++) {
		col = lp->rowind[rbeg + j];
		mpf_EGlpNumSubInnProdTo (dj[col], lp->rowval[rbeg + j], c_dj);
		rcol[col] = 1;
	    }
	}
    }
#if mpf_BASIS_STATS > 0
    printf ("unit rows = %d\n", s_i);
    printf ("nslacks %d, unit rows selected = %d\n", nslacks, nbelem - nslacks);
#endif
    /* now go through remaining cols with dj = 0 */
    tot1 = vd.nfree + vd.nbndone;

    if (mpf_EGlpNumIsEqqual (vd.cmax, mpf_zeroLpNum))
	mpf_EGlpNumOne (cmax);
    else {
	mpf_EGlpNumCopy (cmax, vd.cmax);
	mpf_EGlpNumMultUiTo (cmax, 1000);
    }
    for (j = 0; j < lp->ncols; j++) {
	if (vstat[j] == STAT_BASIC)
	    continue;
	if (icol[j] == 1 || mpf_EGlpNumIsNeqZero (dj[j], mpf_BD_TOLER))
	    continue;
	mcnt = lp->matcnt[j];

	mpf_EGlpNumSet (c_dj, (double) mcnt);
	switch (lp->vtype[j]) {
	case VFREE:
	    porder[tfree] = j;
	    perm[tfree] = tfree;
	    mpf_EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
	    mpf_EGlpNumAddTo (qpenalty[tfree], c_dj);
	    tfree++;
	    break;

	case VLOWER:
	case VUPPER:
	    porder[vd.nfree + tbndone] = j;
	    perm[vd.nfree + tbndone] = tbndone;
	    mpf_EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
	    mpf_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], c_dj);
	    if (lp->vtype[j] == VLOWER)
		mpf_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
	    else
		mpf_EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
	    tbndone++;
	    break;

	case VFIXED:
	case VBOUNDED:
	    porder[tot1 + tbounded] = j;
	    perm[tot1 + tbounded] = tbounded;
	    mpf_EGlpNumCopyFrac (qpenalty[tot1 + tbounded], lp->cz[j], cmax);
	    mpf_EGlpNumAddTo (qpenalty[tot1 + tbounded], lp->lz[j]);
	    mpf_EGlpNumSubTo (qpenalty[tot1 + tbounded], lp->uz[j]);
	    mpf_EGlpNumAddTo (qpenalty[tot1 + tbounded], c_dj);
	    tbounded++;
	    break;
	}
    }
#if mpf_BASIS_STATS > 0
    printf ("bfree %d, bone %d, bbnd %d\n", tfree, tbndone, tbounded);
#endif

    mpf_ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, tfree);
    mpf_ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
	tbndone);
    mpf_ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, tbounded);

    tot2 = tfree + tbndone;
    for (i = 0; i < tbndone; i++) {
	perm[tfree + i] = perm[vd.nfree + i] + tfree;
	porder[tfree + i] = porder[vd.nfree + i];
    }
    for (i = 0; i < tbounded; i++) {
	perm[tot2 + i] = perm[tot1 + i] + tot2;
	porder[tot2 + i] = porder[tot1 + i];
    }
    tot2 += tbounded;

    nbelem =
	mpf_primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
	tot2);
    if (nbelem != lp->nrows) {
	printf ("complain: incorrect final basis size\n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
CLEANUP:
    if (rval)
	mpf_ILLbasis_free_basisinfo (lp);

    ILL_IFFREE (irow, int);
    ILL_IFFREE (rrow, int);
    mpf_EGlpNumFreeArray (v);
    ILL_IFFREE (unitcol, int);
    ILL_IFFREE (icol, int);
    ILL_IFFREE (rcol, int);
    mpf_EGlpNumFreeArray (dj);
    ILL_IFFREE (perm, int);
    ILL_IFFREE (porder, int);
    ILL_IFFREE (plen, int);
    mpf_EGlpNumFreeArray (qpenalty);
    mpf_EGlpNumClearVar (seldj);
    mpf_EGlpNumClearVar (selv);
    mpf_EGlpNumClearVar (c_dj);
    mpf_EGlpNumClearVar (cmax);
    mpf_ILLbasis_clear_vardata (&vd);
    ILL_RETURN (rval, "mpf_ILLbasis_get_initial");
}

static int mpf_set_basis_indices (mpf_lpinfo * lp,
      int *vstat)
{
    int i, b = 0, nb = 0;
    int vs;

    for (i = 0; i < lp->ncols; i++) {
	vs = vstat[i];
	lp->vstat[i] = vs;

	if (vs == STAT_BASIC) {
	    lp->baz[b] = i;
	    lp->vindex[i] = b;
	    b++;
	} else if (vs == STAT_UPPER || vs == STAT_LOWER || vs == STAT_ZERO) {
	    lp->nbaz[nb] = i;
	    lp->vindex[i] = nb;
	    nb++;
	} else {
	    fprintf (stderr, "Error in basis creation\n");
	    return E_SIMPLEX_ERROR;
	}
    }
    if (b != lp->nrows) {
	fprintf (stderr, "Error 2 in basis creation\n");
	return E_SIMPLEX_ERROR;
    } else if (nb != lp->nnbasic) {
	fprintf (stderr, "Error 3 in basis creation\n");
	return E_SIMPLEX_ERROR;
    }
    return 0;
}

int mpf_ILLbasis_get_initial (mpf_lpinfo * lp,
      int algorithm)
{
    int rval = 0;
    int *vstat = NULL;

    mpf_ILLbasis_free_basisinfo (lp);
    mpf_ILLbasis_init_basisinfo (lp);
    rval = mpf_ILLbasis_build_basisinfo (lp);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (vstat, lp->ncols, int);

    if (algorithm == PRIMAL_SIMPLEX)
	rval = mpf_get_initial_basis1 (lp, vstat);
    else
	rval = mpf_get_initial_basis2 (lp, vstat);

    if (rval == E_SIMPLEX_ERROR) {
	FILE *f = fopen ("bad.lp", "w");
	int tval = mpf_ILLwrite_lp_file (lp->O, f, NULL);
	if (tval) {
	    fprintf (stderr, "Error writing bad lp\n");
	}
	if (f != NULL)
	    fclose (f);
    }
    ILL_CLEANUP_IF (rval);

    rval = mpf_set_basis_indices (lp, vstat);
    lp->basisid = 0;

CLEANUP:
    ILL_IFFREE (vstat, int);
    ILL_RETURN (rval, "mpf_ILLbasis_get_initial");
}

static int mpf_choose_basis (int algorithm,
      mpf_t pinf1,
      mpf_t dinf1,
      mpf_t pinf2,
      mpf_t dinf2)
{
    /* We changed the constant definitions outside here, the actual numbers
       are asigned in mpf_lpdata.c. the values are as follows: mpf_CB_EPS =
       0.001; mpf_CB_PRI_RLIMIT = 0.25; mpf_CB_INF_RATIO = 10.0; */
    int choice = 1;
    mpf_t rp, rd;
    if (algorithm == PRIMAL_SIMPLEX) {
	mpf_EGlpNumInitVar (rp);
	mpf_EGlpNumInitVar (rd);
	mpf_EGlpNumCopyDiff (rp, pinf1, pinf2);
	mpf_EGlpNumCopyDiff (rd, dinf1, dinf2);
	if (mpf_EGlpNumIsLeq (rp, mpf_CB_EPS) && mpf_EGlpNumIsLeq (rd, mpf_CB_EPS))
	    choice = 1;
	else {
	    mpf_EGlpNumSign (rp);
	    mpf_EGlpNumSign (rd);
	    if (mpf_EGlpNumIsLeq (rp, mpf_CB_EPS) && mpf_EGlpNumIsLeq (rd, mpf_CB_EPS))
		choice = 2;
	    else if (mpf_EGlpNumIsLess (pinf1, pinf2) && mpf_EGlpNumIsLess (dinf2, dinf1)) {
		choice = 1;
		mpf_EGlpNumCopyFrac (rp, pinf1, pinf2);
		mpf_EGlpNumCopyFrac (rd, dinf2, dinf1);
		mpf_EGlpNumMultTo (rd, mpf_CB_INF_RATIO);
		if (mpf_EGlpNumIsLess (mpf_CB_PRI_RLIMIT, rp) && (mpf_EGlpNumIsLess (rd, rp)))
		    choice = 2;
	    } else if (mpf_EGlpNumIsLess (pinf2, pinf1) && mpf_EGlpNumIsLess (dinf1, dinf2)) {
		choice = 2;
		mpf_EGlpNumCopyFrac (rp, pinf2, pinf1);
		mpf_EGlpNumCopyFrac (rd, dinf1, dinf2);
		mpf_EGlpNumMultTo (rd, mpf_CB_INF_RATIO);
		if (mpf_EGlpNumIsLess (mpf_CB_PRI_RLIMIT, rp) && mpf_EGlpNumIsLess (rd, rp))
		    choice = 1;
	    } else
		choice = 1;
	}
	mpf_EGlpNumClearVar (rp);
	mpf_EGlpNumClearVar (rd);
    }
    ILL_IFTRACE ("%s:%d\n", __func__, choice);
    return choice;
}

int mpf_ILLbasis_get_cinitial (mpf_lpinfo * lp,
      int algorithm)
{
    int rval = 0;
    int *vstat1 = NULL;
    int *vstat2 = NULL;
    int singular;
    int choice = 0;
#if mpf_BASIS_STATS > 0
    int i, nz1 = 0, nz2 = 0;
#endif
    mpf_t pinf1, pinf2, dinf1, dinf2;
    mpf_feas_info fi;
    mpf_EGlpNumInitVar (pinf1);
    mpf_EGlpNumInitVar (pinf2);
    mpf_EGlpNumInitVar (dinf1);
    mpf_EGlpNumInitVar (dinf2);
    mpf_EGlpNumInitVar (fi.totinfeas);

    mpf_ILLbasis_free_basisinfo (lp);
    mpf_ILLbasis_init_basisinfo (lp);
    rval = mpf_ILLbasis_build_basisinfo (lp);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (vstat1, lp->ncols, int);
    ILL_SAFE_MALLOC (vstat2, lp->ncols, int);

    if (algorithm != PRIMAL_SIMPLEX) {
	rval = mpf_get_initial_basis2 (lp, vstat2);
	ILL_CLEANUP_IF (rval);
	rval = mpf_set_basis_indices (lp, vstat2);
	lp->basisid = 0;
	ILL_CLEANUP;
    }
    rval = mpf_get_initial_basis1 (lp, vstat1);
    ILL_CLEANUP_IF (rval);
    rval = mpf_get_initial_basis2 (lp, vstat2);
    ILL_CLEANUP_IF (rval);
    lp->basisid = 0;

    /* handle first basis */
    rval = mpf_set_basis_indices (lp, vstat1);
    ILL_CLEANUP_IF (rval);
#if mpf_BASIS_STATS > 0
    for (i = 0; i < lp->nrows; i++)
	nz1 += lp->matcnt[lp->baz[i]];
#endif
    rval = mpf_ILLbasis_factor (lp, &singular);
    ILL_CLEANUP_IF (rval);

    mpf_ILLfct_compute_piz (lp);
    mpf_ILLfct_compute_dz (lp);
    mpf_ILLfct_dual_adjust (lp, mpf_zeroLpNum);
    mpf_ILLfct_compute_xbz (lp);

    mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
    mpf_EGlpNumCopy (pinf1, lp->pinfeas);
    mpf_EGlpNumCopy (dinf1, lp->dinfeas);
    /*
     * mpf_ILLfct_compute_pobj (lp);  obj1p = lp->objval;
     * mpf_ILLfct_compute_dobj (lp);  obj1d = lp->objval;
     */

    /* handle second basis */
    rval = mpf_set_basis_indices (lp, vstat2);
    ILL_CLEANUP_IF (rval);
#if mpf_BASIS_STATS > 0
    for (i = 0; i < lp->nrows; i++)
	nz2 += lp->matcnt[lp->baz[i]];
#endif
    rval = mpf_ILLbasis_factor (lp, &singular);
    ILL_CLEANUP_IF (rval);

    mpf_ILLfct_compute_piz (lp);
    mpf_ILLfct_compute_dz (lp);
    mpf_ILLfct_dual_adjust (lp, mpf_zeroLpNum);
    mpf_ILLfct_compute_xbz (lp);

    mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
    mpf_EGlpNumCopy (pinf2, lp->pinfeas);
    mpf_EGlpNumCopy (dinf2, lp->dinfeas);

#if mpf_BASIS_STATS > 0
    printf ("b1: nz %d pinf %.2f dinf %.2f\n", nz1, mpf_EGlpNumToLf (pinf1),
	mpf_EGlpNumToLf (dinf1));
    printf ("b2: nz %d pinf %.2f dinf %.2f\n", nz2, mpf_EGlpNumToLf (pinf2),
	mpf_EGlpNumToLf (dinf2));
#endif
    choice = mpf_choose_basis (algorithm, pinf1, dinf1, pinf2, dinf2);
    if (choice == 1) {
	lp->fbasisid = -1;
	rval = mpf_set_basis_indices (lp, vstat1);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:
    if (rval == E_SIMPLEX_ERROR) {
	FILE *fil = fopen ("bad.lp", "w");
	int tval = mpf_ILLwrite_lp_file (lp->O, fil, NULL);
	if (tval) {
	    fprintf (stderr, "Error writing bad lp\n");
	}
	if (fil != NULL)
	    fclose (fil);
    }
    ILL_IFFREE (vstat1, int);
    ILL_IFFREE (vstat2, int);
    mpf_EGlpNumClearVar (pinf1);
    mpf_EGlpNumClearVar (pinf2);
    mpf_EGlpNumClearVar (dinf1);
    mpf_EGlpNumClearVar (dinf2);
    mpf_EGlpNumClearVar (fi.totinfeas);
    ILL_RETURN (rval, "mpf_ILLbasis_get_initial");
}

int mpf_ILLbasis_factor (mpf_lpinfo * lp,
      int *singular)
{
    int rval = 0;
    int i;
    int eindex;
    int lindex;
    int ltype;
    int lvstat;
    int nsing = 0;
    int *singr = 0;
    int *singc = 0;

    *singular = 0;
    do {
	if (lp->f) {
	    mpf_ILLfactor_free_factor_work (lp->f);
	} else {
	    ILL_SAFE_MALLOC (lp->f, 1, mpf_factor_work);
	    mpf_EGlpNumInitVar (lp->f->fzero_tol);
	    mpf_EGlpNumInitVar (lp->f->szero_tol);
	    mpf_EGlpNumInitVar (lp->f->partial_tol);
	    mpf_EGlpNumInitVar (lp->f->maxelem_orig);
	    mpf_EGlpNumInitVar (lp->f->maxelem_factor);
	    mpf_EGlpNumInitVar (lp->f->maxelem_cur);
	    mpf_EGlpNumInitVar (lp->f->partial_cur);
	    mpf_ILLfactor_init_factor_work (lp->f);
	}
	rval = mpf_ILLfactor_create_factor_work (lp->f, lp->nrows);
	ILL_CLEANUP_IF (rval);

	rval = mpf_ILLfactor (lp->f, lp->baz, lp->matbeg, lp->matcnt,
	    lp->matind, lp->matval, &nsing, &singr, &singc);
	ILL_CLEANUP_IF (rval);

	if (nsing != 0) {
	    *singular = 1;
	    for (i = 0; i < nsing; i++) {
		eindex = lp->vindex[lp->O->rowmap[singr[i]]];
		lindex = singc[i];
		ltype = lp->vtype[lp->baz[lindex]];

		if (ltype == VBOUNDED || ltype == VLOWER || ltype == VARTIFICIAL)
		    lvstat = STAT_LOWER;
		else if (ltype == VUPPER)
		    lvstat = STAT_UPPER;
		else
		    lvstat = STAT_ZERO;

		mpf_ILLfct_update_basis_info (lp, eindex, lindex, lvstat);
		lp->basisid++;
	    }
	    ILL_IFFREE (singr, int);
	    ILL_IFFREE (singc, int);
	}
    } while (nsing != 0);

    lp->fbasisid = lp->basisid;

CLEANUP:
    ILL_IFFREE (singr, int);
    ILL_IFFREE (singc, int);
    ILL_RETURN (rval, "mpf_ILLbasis_factor");
}

int mpf_ILLbasis_refactor (mpf_lpinfo * lp)
{
    int sing = 0;
    int rval = 0;

    rval = mpf_ILLbasis_factor (lp, &sing);
    if (sing) {
	fprintf (stderr, "Singular basis in mpf_ILLbasis_refactor()\n");
	rval = -1;
    }
    ILL_RETURN (rval, "mpf_ILLbasis_refactor");
}

void mpf_ILLbasis_column_solve (mpf_lpinfo * lp,
      mpf_svector * rhs,
      mpf_svector * soln)
{
    mpf_ILLfactor_ftran (lp->f, rhs, soln);
}

void mpf_ILLbasis_column_solve_update (mpf_lpinfo * lp,
      mpf_svector * rhs,
      mpf_svector * upd,
      mpf_svector * soln)
{
    mpf_ILLfactor_ftran_update (lp->f, rhs, upd, soln);
}

void mpf_ILLbasis_row_solve (mpf_lpinfo * lp,
      mpf_svector * rhs,
      mpf_svector * soln)
{
    mpf_ILLfactor_btran (lp->f, rhs, soln);
}

int mpf_ILLbasis_update (mpf_lpinfo * lp,
      mpf_svector * y,
      int lindex,
      int *refactor,
      int *singular)
{
#if 0				/* To always refactor, change 0 to 1 */
    *refactor = 1;
    return mpf_ILLbasis_factor (lp, singular);
#else

    int rval = 0;

    *refactor = 0;
    rval = mpf_ILLfactor_update (lp->f, y, lindex, refactor);
    if (rval == E_FACTOR_BLOWUP || rval == E_UPDATE_SINGULAR_ROW
	|| rval == E_UPDATE_SINGULAR_COL) {
	/* Bico - comment out for dist fprintf(stderr, "Warning: numerically
	   bad basis in mpf_ILLfactor_update\n"); */
	*refactor = 1;
	rval = 0;
    }
    if (rval == E_UPDATE_NOSPACE) {
	*refactor = 1;
	rval = 0;
    }
    if (*refactor)
	rval = mpf_ILLbasis_factor (lp, singular);

    if (rval) {
	FILE *eout = 0;
	int tval;

	printf ("write bad lp to factor.lp\n");
	fflush (stdout);
	eout = fopen ("factor.lp", "w");
	if (!eout) {
	    fprintf (stderr, "could not open file to write bad factor lp\n");
	} else {
	    tval = mpf_ILLwrite_lp_file (lp->O, eout, NULL);
	    if (tval) {
		fprintf (stderr, "error while writing bad factor lp\n");
	    }
	    fclose (eout);
	}

	printf ("write bad basis to factor.bas\n");
	fflush (stdout);
	tval = mpf_ILLlib_writebasis (lp, 0, "factor.bas");
	if (tval) {
	    fprintf (stderr, "error while writing factor basis\n");
	}
    }
    ILL_RETURN (rval, "mpf_ILLbasis_update");
#endif
}
