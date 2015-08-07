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

/* RCS_INFO = "$RCSfile: mpq_basis.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "mpq_sortrus.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdefs.h"
#include "mpq_qstruct.h"
#include "mpq_qsopt.h"
#include "mpq_basis.h"
#include "mpq_fct.h"
#include "mpq_lp.h"
#include "mpq_lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define mpq_DJZERO_TOLER mpq_PFEAS_TOLER
#define mpq_BASIS_STATS 0
#define mpq_BASIS_DEBUG 0

void mpq_ILLbasis_init_vardata (mpq_var_data * vd)
{
    memset (vd, 0, sizeof (mpq_var_data));
    mpq_EGlpNumInitVar (vd->cmax);
}

void mpq_ILLbasis_clear_vardata (mpq_var_data * vd)
{
    mpq_EGlpNumClearVar (vd->cmax);
    memset (vd, 0, sizeof (mpq_var_data));
}

static void mpq_get_var_info (mpq_lpinfo * lp,
      mpq_var_data * v);

static int mpq_init_slack_basis (mpq_lpinfo * lp,
      int *vstat,
      int *irow,
      int *rrow,
      int *unitcol,
      int *icol,
      int *rcol),
  mpq_get_initial_basis1 (mpq_lpinfo * lp,
      int *vstat),
  mpq_get_initial_basis2 (mpq_lpinfo * lp,
      int *vstat),
  mpq_set_basis_indices (mpq_lpinfo * lp,
      int *vstat),
  mpq_choose_basis (int algorithm,
      mpq_t pinf1,
      mpq_t dinf1,
      mpq_t pinf2,
      mpq_t dinf2);

void mpq_ILLbasis_init_basisinfo (mpq_lpinfo * lp)
{
    lp->baz = 0;
    lp->nbaz = 0;
    lp->vstat = 0;
    lp->vindex = 0;
    lp->f = 0;
}

void mpq_ILLbasis_free_basisinfo (mpq_lpinfo * lp)
{
    ILL_IFFREE (lp->baz, int);
    ILL_IFFREE (lp->nbaz, int);
    ILL_IFFREE (lp->vstat, int);
    ILL_IFFREE (lp->vindex, int);
    if (lp->f) {
	mpq_ILLfactor_free_factor_work (lp->f);
	mpq_EGlpNumClearVar (lp->f->fzero_tol);
	mpq_EGlpNumClearVar (lp->f->szero_tol);
	mpq_EGlpNumClearVar (lp->f->partial_tol);
	mpq_EGlpNumClearVar (lp->f->maxelem_orig);
	mpq_EGlpNumClearVar (lp->f->maxelem_factor);
	mpq_EGlpNumClearVar (lp->f->maxelem_cur);
	mpq_EGlpNumClearVar (lp->f->partial_cur);
	ILL_IFFREE (lp->f, mpq_factor_work);
    }
}

int mpq_ILLbasis_build_basisinfo (mpq_lpinfo * lp)
{
    int rval = 0;

    ILL_SAFE_MALLOC (lp->baz, lp->nrows, int);
    ILL_SAFE_MALLOC (lp->nbaz, lp->nnbasic, int);
    ILL_SAFE_MALLOC (lp->vstat, lp->ncols, int);
    ILL_SAFE_MALLOC (lp->vindex, lp->ncols, int);

    lp->fbasisid = -1;

CLEANUP:
    if (rval)
	mpq_ILLbasis_free_basisinfo (lp);
    ILL_RETURN (rval, "mpq_ILLbasis_build_basisinfo");
}

int mpq_ILLbasis_load (mpq_lpinfo * lp,
      mpq_ILLlp_basis * B)
{
    int rval = 0;
    char *cstat = B->cstat;
    char *rstat = B->rstat;
    int *structmap = lp->O->structmap;
    int *rowmap = lp->O->rowmap;
    mpq_t *rng = lp->O->rangeval;
    int i, j, ncols = lp->ncols, nrows = lp->nrows, nstruct = lp->O->nstruct;
    int basic = 0, nonbasic = 0;

    mpq_ILLbasis_free_basisinfo (lp);
    mpq_ILLbasis_init_basisinfo (lp);
    rval = mpq_ILLbasis_build_basisinfo (lp);
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
	if (rng && (mpq_EGlpNumIsNeqqZero (rng[i]))) {
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

    ILL_RETURN (rval, "mpq_ILLbasis_load");
}

int mpq_ILLbasis_tableau_row (mpq_lpinfo * lp,
      int row,
      mpq_t * brow,
      mpq_t * trow,
      mpq_t * rhs,
      int strict)
{
    int rval = 0;
    int i;
    int singular = 0;
    int indx;
    mpq_t coef;
    mpq_t sum;
    mpq_svector z, zA;
    mpq_EGlpNumInitVar (coef);
    mpq_EGlpNumInitVar (sum);
    mpq_EGlpNumZero (sum);

    mpq_ILLsvector_init (&z);
    mpq_ILLsvector_init (&zA);

    if (lp->basisid == -1) {
	fprintf (stderr, "mpq_ILLbasis_tableau_row: no basis\n");
	rval = E_GENERAL_ERROR;
	ILL_CLEANUP;
    }
    if (lp->fbasisid != lp->basisid) {	/* Needs to be changed */
	rval = mpq_ILLbasis_factor (lp, &singular);
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
    rval = mpq_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);
    mpq_ILLfct_compute_zz (lp, &z, row);

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumZero (brow[i]);
    for (i = 0; i < z.nzcnt; i++) {
	indx = z.indx[i];
	mpq_EGlpNumCopy (coef, z.coef[i]);
	mpq_EGlpNumCopy (brow[indx], coef);
	mpq_EGlpNumAddInnProdTo (sum, coef, lp->bz[indx]);
    }

    if (rhs != NULL)
	mpq_EGlpNumCopy (*rhs, sum);
    if (trow != NULL) {
	if (!strict) {
	    rval = mpq_ILLsvector_alloc (&zA, lp->ncols);
	    if (rval)
		ILL_CLEANUP;
	    ILL_IFTRACE ("%s:\n", __func__);
	    rval = mpq_ILLfct_compute_zA (lp, &z, &zA);
	    ILL_CLEANUP_IF (rval);

	    for (i = 0; i < lp->ncols; i++)
		mpq_EGlpNumZero (trow[i]);
	    for (i = 0; i < zA.nzcnt; i++)
		mpq_EGlpNumCopy (trow[lp->nbaz[zA.indx[i]]], zA.coef[i]);
	    mpq_EGlpNumOne (trow[lp->baz[row]]);
	} else {
	    mpq_ILLfct_compute_vA (lp, &z, trow);
	}
    }
#if mpq_BASIS_DEBUG > 0
    if (rhs != NULL && trow != NULL) {
	mpq_t *tr = NULL;
	mpq_EGlpNumZero (sum);
	if (strict)
	    tr = trow;
	else {
	    tr = mpq_EGlpNumAllocArray (lp->ncols);
	    mpq_ILLfct_compute_vA (lp, &z, tr);
	}
	for (i = 0; i < lp->nrows; i++)
	    if (mpq_EGlpNumIsLess (mpq_zeroLpNum, tr[lp->baz[i]]))
		mpq_EGlpNumAddTo (sum, tr[lp->baz[i]]);
	    else
		mpq_EGlpNumSubTo (sum, tr[lp->baz[i]]);
	mpq_EGlpNumCopy (coef, mpq_oneLpNum);
	mpq_EGlpNumSubTo (coef, sum);
	if (mpq_EGlpNumIsLess (coef, mpq_zeroLpNum))
	    mpq_EGlpNumSign (coef);
	if (mpq_EGlpNumIsLess (mpq_PIVZ_TOLER, coef))
	    fprintf (stderr, "tableau: bas computed = %.12f\n", mpq_EGlpNumToLf (sum));
	if (!strict)
	    mpq_EGlpNumFreeArray (tr);
#if mpq_BASIS_DEBUG > 1
	mpq_EGlpNumZero (sum);
	for (i = 0; i < lp->ncols; i++) {
	    if (lp->vstat[i] == STAT_BASIC)
		mpq_EGlpNumAddInnProdTo (sum, lp->xbz[lp->vindex[i]], trow[i]);
	    else if (lp->vstat[i] == STAT_UPPER)
		mpq_EGlpNumAddInnProdTo (sum, lp->uz[i], trow[i]);
	    else if (lp->vstat[i] == STAT_LOWER)
		mpq_EGlpNumAddInnProdTo (sum, lp->lz[i], trow[i]);
	}
	mpq_EGlpNumSet (coef, 1e-10);
	if (mpq_EGlpNumIsNeq (sum, *rhs, coef))
	    fprintf (stderr, "tableau rhs = %.9f, computed = %.9f\n",
		mpq_EGlpNumToLf (*rhs), mpq_EGlpNumToLf (sum));
#endif
    }
#endif

CLEANUP:
    mpq_ILLsvector_free (&z);
    mpq_ILLsvector_free (&zA);
    mpq_EGlpNumClearVar (coef);
    mpq_EGlpNumClearVar (sum);
    return rval;
}

static void mpq_get_var_info (mpq_lpinfo * lp,
      mpq_var_data * v)
{
    int i = 0;

    v->nartif = 0;
    v->nslacks = 0;
    v->nfree = 0;
    v->nbndone = 0;
    v->nbounded = 0;
    v->nfixed = 0;
    mpq_EGlpNumCopy (v->cmax, mpq_NINFTY);

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
	mpq_EGlpNumSetToMaxAbs (v->cmax, lp->cz[i]);
    }

#if mpq_BASIS_STATS > 0
    printf ("cols = %d, acols = %d, total  = %d, nrows = %d, nlog = %d\n",
	lp->ncols, lp->ncols - lp->nrows,
	v->nartif + v->nfree + v->nslacks + v->nbndone + v->nbounded,
	lp->nrows, v->nartif + v->nslacks);
#endif
}

static int mpq_init_slack_basis (mpq_lpinfo * lp,
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
	    if (fabs (mpq_EGlpNumToLf (lp->lz[j])) < fabs (mpq_EGlpNumToLf (lp->uz[j])))
		vstat[j] = STAT_LOWER;
	    else
		vstat[j] = STAT_UPPER;
	}
    }
    return nslacks;
}

static int mpq_primal_col_select (mpq_lpinfo * lp,
      int *vstat,
      int *irow,
      int *rrow,
      int *unitcol,
      mpq_t * v,
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
    mpq_t *matval = lp->matval;
    mpq_t alpha, val, maxelem;
    mpq_EGlpNumInitVar (alpha);
    mpq_EGlpNumInitVar (val);
    mpq_EGlpNumInitVar (maxelem);

    for (k = 0; k < pcols; k++) {
	j = porder[perm[k]];
	mcnt = matcnt[j];
	mbeg = matbeg[j];

	mpq_EGlpNumCopy (alpha, mpq_NINFTY);
	mpq_EGlpNumCopy (maxelem, mpq_NINFTY);

	for (i = 0; i < mcnt; i++) {
	    mpq_EGlpNumCopyAbs (val, matval[mbeg + i]);
	    if (mpq_EGlpNumIsLess (maxelem, val))
		mpq_EGlpNumCopy (maxelem, val);
	    if (rrow[matind[mbeg + i]] == 0 && mpq_EGlpNumIsLess (alpha, val)) {
		mpq_EGlpNumCopy (alpha, val);
		r = matind[mbeg + i];
	    }
	}
	mpq_EGlpNumCopy (val, maxelem);
	mpq_EGlpNumMultTo (val, mpq_PARAM_IBASIS_RPIVOT);
	if (mpq_EGlpNumIsLess (val, alpha)) {
	    vstat[j] = STAT_BASIC;
	    nbelem++;
	    irow[r] = 1;
	    mpq_EGlpNumCopy (v[r], alpha);
	    for (i = 0; i < mcnt; i++)
		if (mpq_EGlpNumIsNeqqZero (matval[mbeg + i]))
		    rrow[matind[mbeg + i]]++;
	} else {
	    mpq_EGlpNumCopy (alpha, mpq_NINFTY);
	    for (i = 0; i < mcnt; i++) {
		tr = matind[mbeg + i];
		mpq_EGlpNumCopyAbs (val, matval[mbeg + i]);
		mpq_EGlpNumDivTo (val, mpq_PARAM_IBASIS_RTRIANG);
		if (mpq_EGlpNumIsNeqq (v[tr], mpq_INFTY) && mpq_EGlpNumIsLess (v[tr], val)) {
		    mpq_EGlpNumZero (alpha);
		    break;
		}
		mpq_EGlpNumCopyAbs (val, matval[mbeg + i]);
		if (irow[tr] == 0 && mpq_EGlpNumIsLess (alpha, val)) {
		    mpq_EGlpNumCopy (alpha, val);
		    r = tr;
		}
	    }
	    if (mpq_EGlpNumIsNeqqZero (alpha) && mpq_EGlpNumIsNeqq (alpha, mpq_NINFTY)) {
		vstat[j] = STAT_BASIC;
		nbelem++;
		irow[r] = 1;
		mpq_EGlpNumCopy (v[r], alpha);
		for (i = 0; i < mcnt; i++)
		    if (mpq_EGlpNumIsNeqqZero (matval[mbeg + i]))
			rrow[matind[mbeg + i]]++;
	    }
	}
    }
#if mpq_BASIS_STATS > 0
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
    mpq_EGlpNumClearVar (alpha);
    mpq_EGlpNumClearVar (val);
    mpq_EGlpNumClearVar (maxelem);
    return nbelem;
}

/* This is an implementation of the initial basis procedure in: "Implementing
   the simplex method: the initial basis", by Bob Bixby. Goals: choose
   initial variables to go into basis which satisfy: 1) vars are slacks, 2)
   vars have freedom to move 3) initial submatrix is nonsingular, 4) low
   objective function contribution. */
static int mpq_get_initial_basis1 (mpq_lpinfo * lp,
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
    mpq_t cmax;
    mpq_t *v = NULL;
    mpq_t *qpenalty = NULL;
    mpq_var_data vd;
    mpq_ILLbasis_init_vardata (&vd);
    mpq_EGlpNumInitVar (cmax);

    mpq_get_var_info (lp, &vd);
    if (mpq_EGlpNumIsEqqual (vd.cmax, mpq_zeroLpNum))
	mpq_EGlpNumOne (cmax);
    else {
	mpq_EGlpNumCopy (cmax, vd.cmax);
	mpq_EGlpNumMultUiTo (cmax, 1000);
    }

    ILL_SAFE_MALLOC (irow, lp->nrows, int);
    ILL_SAFE_MALLOC (rrow, lp->nrows, int);
    v = mpq_EGlpNumAllocArray (lp->nrows);
    ILL_SAFE_MALLOC (unitcol, lp->nrows, int);

    for (i = 0; i < lp->nrows; i++) {
	unitcol[i] = -1;
	mpq_EGlpNumCopy (v[i], mpq_INFTY);
	irow[i] = 0;
	rrow[i] = 0;
    }

    nslacks = mpq_init_slack_basis (lp, vstat, irow, rrow, unitcol, NULL, NULL);
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
    qpenalty = mpq_EGlpNumAllocArray (tot2);

    for (j = 0; j < lp->ncols; j++) {
	if (vstat[j] == STAT_BASIC)
	    continue;

	switch (lp->vtype[j]) {
	case VFREE:
	    porder[tfree] = j;
	    perm[tfree] = tfree;
	    mpq_EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
	    tfree++;
	    break;

	case VLOWER:
	case VUPPER:
	    porder[vd.nfree + tbndone] = j;
	    perm[vd.nfree + tbndone] = tbndone;
	    mpq_EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
	    if (lp->vtype[j] == VLOWER)
		mpq_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
	    else
		mpq_EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
	    tbndone++;
	    break;

	case VFIXED:
	case VBOUNDED:
	    porder[tot1 + tbounded] = j;
	    perm[tot1 + tbounded] = tbounded;
	    mpq_EGlpNumCopyFrac (qpenalty[tot1 + tbndone], lp->cz[j], cmax);
	    mpq_EGlpNumAddTo (qpenalty[tot1 + tbndone], lp->lz[j]);
	    mpq_EGlpNumSubTo (qpenalty[tot1 + tbndone], lp->uz[j]);
	    tbounded++;
	    break;
	}
    }
    if (tfree != vd.nfree || tbndone != vd.nbndone || tbounded != vd.nbounded) {
	printf ("complain: incorrect basis info \n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
    mpq_ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, vd.nfree);
    mpq_ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
	vd.nbndone);
    mpq_ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, vd.nbounded);

    for (i = 0; i < vd.nbndone; i++)
	perm[vd.nfree + i] += vd.nfree;
    for (i = 0; i < vd.nbounded; i++)
	perm[tot1 + i] += tot1;

    nbelem =
	mpq_primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
	tot2);
    if (nbelem != lp->nrows) {
	printf ("complain: incorrect final basis size\n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
CLEANUP:
    mpq_EGlpNumClearVar (cmax);
    if (rval)
	mpq_ILLbasis_free_basisinfo (lp);
    ILL_IFFREE (irow, int);
    ILL_IFFREE (rrow, int);
    mpq_EGlpNumFreeArray (v);
    ILL_IFFREE (perm, int);
    ILL_IFFREE (porder, int);
    ILL_IFFREE (unitcol, int);
    mpq_EGlpNumFreeArray (qpenalty);
    mpq_ILLbasis_clear_vardata (&vd);
    ILL_RETURN (rval, "mpq_ILLbasis_get_initial");
}

static int mpq_get_initial_basis2 (mpq_lpinfo * lp,
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
    mpq_t *v = NULL;
    mpq_t *qpenalty = NULL;
    int col = 0, s_i = 0, selc = 0;
    int *icol = NULL, *rcol = NULL;
    int *plen = NULL;
    mpq_t *dj = NULL;
    mpq_var_data vd;
    mpq_t seldj;
    mpq_t selv;
    mpq_t c_dj;
    mpq_t cmax;
    mpq_EGlpNumInitVar (seldj);
    mpq_EGlpNumInitVar (selv);
    mpq_EGlpNumInitVar (c_dj);
    mpq_EGlpNumInitVar (cmax);
    mpq_EGlpNumZero (c_dj);
    mpq_EGlpNumZero (selv);
    mpq_EGlpNumZero (seldj);
    mpq_ILLbasis_init_vardata (&vd);

    mpq_get_var_info (lp, &vd);

    ILL_SAFE_MALLOC (irow, lp->nrows, int);
    ILL_SAFE_MALLOC (rrow, lp->nrows, int);
    v = mpq_EGlpNumAllocArray (lp->nrows);
    ILL_SAFE_MALLOC (unitcol, lp->nrows, int);
    ILL_SAFE_MALLOC (icol, lp->ncols, int);
    ILL_SAFE_MALLOC (rcol, lp->ncols, int);
    dj = mpq_EGlpNumAllocArray (lp->ncols);

    for (i = 0; i < lp->nrows; i++) {
	unitcol[i] = -1;
	mpq_EGlpNumCopy (v[i], mpq_INFTY);
	irow[i] = 0;
	rrow[i] = 0;
    }
    /* assign all d_j */
    for (i = 0; i < lp->ncols; i++) {
	icol[i] = 0;
	rcol[i] = 0;
	mpq_EGlpNumCopy (dj[i], lp->cz[i]);
    }

    nslacks = mpq_init_slack_basis (lp, vstat, irow, rrow, unitcol, icol, rcol);
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
    qpenalty = mpq_EGlpNumAllocArray (lp->ncols);

    /* find all unit rows and record lengths */
    for (i = 0; i < lp->nrows; i++) {
	if (irow[i] != 1) {
	    rbeg = lp->rowbeg[i];
	    rcnt = lp->rowcnt[i];
	    for (j = 0; j < rcnt; j++) {
		mpq_EGlpNumCopyAbs (cmax, lp->rowval[rbeg + j]);
		if (mpq_EGlpNumIsNeqq (cmax, mpq_oneLpNum))
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
    mpq_ILLutil_int_perm_quicksort (perm, plen, s_i);

    /* now go through the unit rows */
    for (k = 0; k < s_i; k++) {
	i = porder[perm[k]];
	rbeg = lp->rowbeg[i];
	rcnt = lp->rowcnt[i];
	selc = -1;
	mpq_EGlpNumCopy (seldj, mpq_INFTY);
	mpq_EGlpNumZero (selv);

	/* for every row s_i, compute min {d_j : d_j <0 , j is u or l or fr} */
	for (j = 0; j < rcnt; j++) {
	    col = lp->rowind[rbeg + j];
	    if (rcol[col] == 1)
		break;
	    if (mpq_EGlpNumIsLess (dj[col], mpq_zeroLpNum)) {
		if (mpq_EGlpNumIsLess (dj[col], seldj)) {
		    selc = col;
		    mpq_EGlpNumCopy (seldj, dj[col]);
		    mpq_EGlpNumCopy (selv, lp->rowval[rbeg + j]);
		}
	    }
	}
	/* select pivot element and update all d_j's */
	if (selc != -1) {
	    nbelem++;
	    irow[i] = 1;
	    rrow[i] = 1;
	    icol[selc] = 1;
	    mpq_EGlpNumCopyFrac (c_dj, dj[selc], selv);
	    vstat[selc] = STAT_BASIC;
	    for (j = 0; j < rcnt; j++) {
		col = lp->rowind[rbeg + j];
		mpq_EGlpNumSubInnProdTo (dj[col], lp->rowval[rbeg + j], c_dj);
		rcol[col] = 1;
	    }
	}
    }
#if mpq_BASIS_STATS > 0
    printf ("unit rows = %d\n", s_i);
    printf ("nslacks %d, unit rows selected = %d\n", nslacks, nbelem - nslacks);
#endif
    /* now go through remaining cols with dj = 0 */
    tot1 = vd.nfree + vd.nbndone;

    if (mpq_EGlpNumIsEqqual (vd.cmax, mpq_zeroLpNum))
	mpq_EGlpNumOne (cmax);
    else {
	mpq_EGlpNumCopy (cmax, vd.cmax);
	mpq_EGlpNumMultUiTo (cmax, 1000);
    }
    for (j = 0; j < lp->ncols; j++) {
	if (vstat[j] == STAT_BASIC)
	    continue;
	if (icol[j] == 1 || mpq_EGlpNumIsNeqZero (dj[j], mpq_BD_TOLER))
	    continue;
	mcnt = lp->matcnt[j];

	mpq_EGlpNumSet (c_dj, (double) mcnt);
	switch (lp->vtype[j]) {
	case VFREE:
	    porder[tfree] = j;
	    perm[tfree] = tfree;
	    mpq_EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
	    mpq_EGlpNumAddTo (qpenalty[tfree], c_dj);
	    tfree++;
	    break;

	case VLOWER:
	case VUPPER:
	    porder[vd.nfree + tbndone] = j;
	    perm[vd.nfree + tbndone] = tbndone;
	    mpq_EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
	    mpq_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], c_dj);
	    if (lp->vtype[j] == VLOWER)
		mpq_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
	    else
		mpq_EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
	    tbndone++;
	    break;

	case VFIXED:
	case VBOUNDED:
	    porder[tot1 + tbounded] = j;
	    perm[tot1 + tbounded] = tbounded;
	    mpq_EGlpNumCopyFrac (qpenalty[tot1 + tbounded], lp->cz[j], cmax);
	    mpq_EGlpNumAddTo (qpenalty[tot1 + tbounded], lp->lz[j]);
	    mpq_EGlpNumSubTo (qpenalty[tot1 + tbounded], lp->uz[j]);
	    mpq_EGlpNumAddTo (qpenalty[tot1 + tbounded], c_dj);
	    tbounded++;
	    break;
	}
    }
#if mpq_BASIS_STATS > 0
    printf ("bfree %d, bone %d, bbnd %d\n", tfree, tbndone, tbounded);
#endif

    mpq_ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, tfree);
    mpq_ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
	tbndone);
    mpq_ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, tbounded);

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
	mpq_primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
	tot2);
    if (nbelem != lp->nrows) {
	printf ("complain: incorrect final basis size\n");
	rval = E_SIMPLEX_ERROR;
	ILL_CLEANUP;
    }
CLEANUP:
    if (rval)
	mpq_ILLbasis_free_basisinfo (lp);

    ILL_IFFREE (irow, int);
    ILL_IFFREE (rrow, int);
    mpq_EGlpNumFreeArray (v);
    ILL_IFFREE (unitcol, int);
    ILL_IFFREE (icol, int);
    ILL_IFFREE (rcol, int);
    mpq_EGlpNumFreeArray (dj);
    ILL_IFFREE (perm, int);
    ILL_IFFREE (porder, int);
    ILL_IFFREE (plen, int);
    mpq_EGlpNumFreeArray (qpenalty);
    mpq_EGlpNumClearVar (seldj);
    mpq_EGlpNumClearVar (selv);
    mpq_EGlpNumClearVar (c_dj);
    mpq_EGlpNumClearVar (cmax);
    mpq_ILLbasis_clear_vardata (&vd);
    ILL_RETURN (rval, "mpq_ILLbasis_get_initial");
}

static int mpq_set_basis_indices (mpq_lpinfo * lp,
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

int mpq_ILLbasis_get_initial (mpq_lpinfo * lp,
      int algorithm)
{
    int rval = 0;
    int *vstat = NULL;

    mpq_ILLbasis_free_basisinfo (lp);
    mpq_ILLbasis_init_basisinfo (lp);
    rval = mpq_ILLbasis_build_basisinfo (lp);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (vstat, lp->ncols, int);

    if (algorithm == PRIMAL_SIMPLEX)
	rval = mpq_get_initial_basis1 (lp, vstat);
    else
	rval = mpq_get_initial_basis2 (lp, vstat);

    if (rval == E_SIMPLEX_ERROR) {
	FILE *f = fopen ("bad.lp", "w");
	int tval = mpq_ILLwrite_lp_file (lp->O, f, NULL);
	if (tval) {
	    fprintf (stderr, "Error writing bad lp\n");
	}
	if (f != NULL)
	    fclose (f);
    }
    ILL_CLEANUP_IF (rval);

    rval = mpq_set_basis_indices (lp, vstat);
    lp->basisid = 0;

CLEANUP:
    ILL_IFFREE (vstat, int);
    ILL_RETURN (rval, "mpq_ILLbasis_get_initial");
}

static int mpq_choose_basis (int algorithm,
      mpq_t pinf1,
      mpq_t dinf1,
      mpq_t pinf2,
      mpq_t dinf2)
{
    /* We changed the constant definitions outside here, the actual numbers
       are asigned in mpq_lpdata.c. the values are as follows: mpq_CB_EPS =
       0.001; mpq_CB_PRI_RLIMIT = 0.25; mpq_CB_INF_RATIO = 10.0; */
    int choice = 1;
    mpq_t rp, rd;
    if (algorithm == PRIMAL_SIMPLEX) {
	mpq_EGlpNumInitVar (rp);
	mpq_EGlpNumInitVar (rd);
	mpq_EGlpNumCopyDiff (rp, pinf1, pinf2);
	mpq_EGlpNumCopyDiff (rd, dinf1, dinf2);
	if (mpq_EGlpNumIsLeq (rp, mpq_CB_EPS) && mpq_EGlpNumIsLeq (rd, mpq_CB_EPS))
	    choice = 1;
	else {
	    mpq_EGlpNumSign (rp);
	    mpq_EGlpNumSign (rd);
	    if (mpq_EGlpNumIsLeq (rp, mpq_CB_EPS) && mpq_EGlpNumIsLeq (rd, mpq_CB_EPS))
		choice = 2;
	    else if (mpq_EGlpNumIsLess (pinf1, pinf2) && mpq_EGlpNumIsLess (dinf2, dinf1)) {
		choice = 1;
		mpq_EGlpNumCopyFrac (rp, pinf1, pinf2);
		mpq_EGlpNumCopyFrac (rd, dinf2, dinf1);
		mpq_EGlpNumMultTo (rd, mpq_CB_INF_RATIO);
		if (mpq_EGlpNumIsLess (mpq_CB_PRI_RLIMIT, rp) && (mpq_EGlpNumIsLess (rd, rp)))
		    choice = 2;
	    } else if (mpq_EGlpNumIsLess (pinf2, pinf1) && mpq_EGlpNumIsLess (dinf1, dinf2)) {
		choice = 2;
		mpq_EGlpNumCopyFrac (rp, pinf2, pinf1);
		mpq_EGlpNumCopyFrac (rd, dinf1, dinf2);
		mpq_EGlpNumMultTo (rd, mpq_CB_INF_RATIO);
		if (mpq_EGlpNumIsLess (mpq_CB_PRI_RLIMIT, rp) && mpq_EGlpNumIsLess (rd, rp))
		    choice = 1;
	    } else
		choice = 1;
	}
	mpq_EGlpNumClearVar (rp);
	mpq_EGlpNumClearVar (rd);
    }
    ILL_IFTRACE ("%s:%d\n", __func__, choice);
    return choice;
}

int mpq_ILLbasis_get_cinitial (mpq_lpinfo * lp,
      int algorithm)
{
    int rval = 0;
    int *vstat1 = NULL;
    int *vstat2 = NULL;
    int singular;
    int choice = 0;
#if mpq_BASIS_STATS > 0
    int i, nz1 = 0, nz2 = 0;
#endif
    mpq_t pinf1, pinf2, dinf1, dinf2;
    mpq_feas_info fi;
    mpq_EGlpNumInitVar (pinf1);
    mpq_EGlpNumInitVar (pinf2);
    mpq_EGlpNumInitVar (dinf1);
    mpq_EGlpNumInitVar (dinf2);
    mpq_EGlpNumInitVar (fi.totinfeas);

    mpq_ILLbasis_free_basisinfo (lp);
    mpq_ILLbasis_init_basisinfo (lp);
    rval = mpq_ILLbasis_build_basisinfo (lp);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (vstat1, lp->ncols, int);
    ILL_SAFE_MALLOC (vstat2, lp->ncols, int);

    if (algorithm != PRIMAL_SIMPLEX) {
	rval = mpq_get_initial_basis2 (lp, vstat2);
	ILL_CLEANUP_IF (rval);
	rval = mpq_set_basis_indices (lp, vstat2);
	lp->basisid = 0;
	ILL_CLEANUP;
    }
    rval = mpq_get_initial_basis1 (lp, vstat1);
    ILL_CLEANUP_IF (rval);
    rval = mpq_get_initial_basis2 (lp, vstat2);
    ILL_CLEANUP_IF (rval);
    lp->basisid = 0;

    /* handle first basis */
    rval = mpq_set_basis_indices (lp, vstat1);
    ILL_CLEANUP_IF (rval);
#if mpq_BASIS_STATS > 0
    for (i = 0; i < lp->nrows; i++)
	nz1 += lp->matcnt[lp->baz[i]];
#endif
    rval = mpq_ILLbasis_factor (lp, &singular);
    ILL_CLEANUP_IF (rval);

    mpq_ILLfct_compute_piz (lp);
    mpq_ILLfct_compute_dz (lp);
    mpq_ILLfct_dual_adjust (lp);
    mpq_ILLfct_compute_xbz (lp);

    mpq_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
    mpq_ILLfct_check_dfeasible (lp, &fi);
    mpq_EGlpNumCopy (pinf1, lp->pinfeas);
    mpq_EGlpNumCopy (dinf1, lp->dinfeas);
    /*
     * mpq_ILLfct_compute_pobj (lp);  obj1p = lp->objval;
     * mpq_ILLfct_compute_dobj (lp);  obj1d = lp->objval;
     */

    /* handle second basis */
    rval = mpq_set_basis_indices (lp, vstat2);
    ILL_CLEANUP_IF (rval);
#if mpq_BASIS_STATS > 0
    for (i = 0; i < lp->nrows; i++)
	nz2 += lp->matcnt[lp->baz[i]];
#endif
    rval = mpq_ILLbasis_factor (lp, &singular);
    ILL_CLEANUP_IF (rval);

    mpq_ILLfct_compute_piz (lp);
    mpq_ILLfct_compute_dz (lp);
    mpq_ILLfct_dual_adjust (lp);
    mpq_ILLfct_compute_xbz (lp);

    mpq_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
    mpq_ILLfct_check_dfeasible (lp, &fi);
    mpq_EGlpNumCopy (pinf2, lp->pinfeas);
    mpq_EGlpNumCopy (dinf2, lp->dinfeas);

#if mpq_BASIS_STATS > 0
    printf ("b1: nz %d pinf %.2f dinf %.2f\n", nz1, mpq_EGlpNumToLf (pinf1),
	mpq_EGlpNumToLf (dinf1));
    printf ("b2: nz %d pinf %.2f dinf %.2f\n", nz2, mpq_EGlpNumToLf (pinf2),
	mpq_EGlpNumToLf (dinf2));
#endif
    choice = mpq_choose_basis (algorithm, pinf1, dinf1, pinf2, dinf2);
    if (choice == 1) {
	lp->fbasisid = -1;
	rval = mpq_set_basis_indices (lp, vstat1);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:
    if (rval == E_SIMPLEX_ERROR) {
	FILE *fil = fopen ("bad.lp", "w");
	int tval = mpq_ILLwrite_lp_file (lp->O, fil, NULL);
	if (tval) {
	    fprintf (stderr, "Error writing bad lp\n");
	}
	if (fil != NULL)
	    fclose (fil);
    }
    ILL_IFFREE (vstat1, int);
    ILL_IFFREE (vstat2, int);
    mpq_EGlpNumClearVar (pinf1);
    mpq_EGlpNumClearVar (pinf2);
    mpq_EGlpNumClearVar (dinf1);
    mpq_EGlpNumClearVar (dinf2);
    mpq_EGlpNumClearVar (fi.totinfeas);
    ILL_RETURN (rval, "mpq_ILLbasis_get_initial");
}

int mpq_ILLbasis_factor (mpq_lpinfo * lp,
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
	    mpq_ILLfactor_free_factor_work (lp->f);
	} else {
	    ILL_SAFE_MALLOC (lp->f, 1, mpq_factor_work);
	    mpq_EGlpNumInitVar (lp->f->fzero_tol);
	    mpq_EGlpNumInitVar (lp->f->szero_tol);
	    mpq_EGlpNumInitVar (lp->f->partial_tol);
	    mpq_EGlpNumInitVar (lp->f->maxelem_orig);
	    mpq_EGlpNumInitVar (lp->f->maxelem_factor);
	    mpq_EGlpNumInitVar (lp->f->maxelem_cur);
	    mpq_EGlpNumInitVar (lp->f->partial_cur);
	    mpq_ILLfactor_init_factor_work (lp->f);
	}
	rval = mpq_ILLfactor_create_factor_work (lp->f, lp->nrows);
	ILL_CLEANUP_IF (rval);

	rval = mpq_ILLfactor (lp->f, lp->baz, lp->matbeg, lp->matcnt,
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

		mpq_ILLfct_update_basis_info (lp, eindex, lindex, lvstat);
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
    ILL_RETURN (rval, "mpq_ILLbasis_factor");
}

int mpq_ILLbasis_refactor (mpq_lpinfo * lp)
{
    int sing = 0;
    int rval = 0;

    rval = mpq_ILLbasis_factor (lp, &sing);
    if (sing) {
	fprintf (stderr, "Singular basis in mpq_ILLbasis_refactor()\n");
	rval = -1;
    }
    ILL_RETURN (rval, "mpq_ILLbasis_refactor");
}

void mpq_ILLbasis_column_solve (mpq_lpinfo * lp,
      mpq_svector * rhs,
      mpq_svector * soln)
{
    mpq_ILLfactor_ftran (lp->f, rhs, soln);
}

void mpq_ILLbasis_column_solve_update (mpq_lpinfo * lp,
      mpq_svector * rhs,
      mpq_svector * upd,
      mpq_svector * soln)
{
    mpq_ILLfactor_ftran_update (lp->f, rhs, upd, soln);
}

void mpq_ILLbasis_row_solve (mpq_lpinfo * lp,
      mpq_svector * rhs,
      mpq_svector * soln)
{
    mpq_ILLfactor_btran (lp->f, rhs, soln);
}

int mpq_ILLbasis_update (mpq_lpinfo * lp,
      mpq_svector * y,
      int lindex,
      int *refactor,
      int *singular)
{
#if 0				/* To always refactor, change 0 to 1 */
    *refactor = 1;
    return mpq_ILLbasis_factor (lp, singular);
#else

    int rval = 0;

    *refactor = 0;
    rval = mpq_ILLfactor_update (lp->f, y, lindex, refactor);
    if (rval == E_FACTOR_BLOWUP || rval == E_UPDATE_SINGULAR_ROW
	|| rval == E_UPDATE_SINGULAR_COL) {
	/* Bico - comment out for dist fprintf(stderr, "Warning: numerically
	   bad basis in mpq_ILLfactor_update\n"); */
	*refactor = 1;
	rval = 0;
    }
    if (rval == E_UPDATE_NOSPACE) {
	*refactor = 1;
	rval = 0;
    }
    if (*refactor)
	rval = mpq_ILLbasis_factor (lp, singular);

    if (rval) {
	FILE *eout = 0;
	int tval;

	printf ("write bad lp to factor.lp\n");
	fflush (stdout);
	eout = fopen ("factor.lp", "w");
	if (!eout) {
	    fprintf (stderr, "could not open file to write bad factor lp\n");
	} else {
	    tval = mpq_ILLwrite_lp_file (lp->O, eout, NULL);
	    if (tval) {
		fprintf (stderr, "error while writing bad factor lp\n");
	    }
	    fclose (eout);
	}

	printf ("write bad basis to factor.bas\n");
	fflush (stdout);
	tval = mpq_ILLlib_writebasis (lp, 0, "factor.bas");
	if (tval) {
	    fprintf (stderr, "error while writing factor basis\n");
	}
    }
    ILL_RETURN (rval, "mpq_ILLbasis_update");
#endif
}
