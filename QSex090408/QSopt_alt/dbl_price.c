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

/* RCS_INFO = "$RCSfile: dbl_price.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;
#include "econfig.h"
#include "stddefs.h"
#include "dbl_qsopt.h"
#include "dbl_lpdefs.h"
#include "dbl_fct.h"
#include "dbl_price.h"
#include "dbl_basis.h"
#include "dbl_iqsutil.h"
#include "dbl_dstruct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define  dbl_MULTIP 1
#define  dbl_PRICE_DEBUG 0

static void dbl_update_d_scaleinf (dbl_price_info * const p,
      dbl_heap * const h,
      int const j,
      double inf,
      int const prule),
  dbl_update_p_scaleinf (dbl_price_info * const p,
      dbl_heap * const h,
      int const i,
      double inf,
      int const prule);

static void dbl_compute_dualI_inf (dbl_lpinfo * const lp,
      int const j,
      double *const inf),
  dbl_compute_dualII_inf (dbl_lpinfo * const lp,
      int const j,
      double *const inf),
  dbl_compute_primalI_inf (dbl_lpinfo * const lp,
      int const i,
      double *const inf),
  dbl_compute_primalII_inf (dbl_lpinfo * const lp,
      int const i,
      double *const inf);

void dbl_ILLprice_free_heap (dbl_price_info * const pinf)
{
    dbl_ILLheap_free (&(pinf->h));
}

int dbl_ILLprice_build_heap (dbl_price_info * const pinf,
      int const nkeys,
      double *keylist)
{
    dbl_ILLheap_init (&(pinf->h));
    dbl_EGlpNumSet (pinf->htrigger,
	1.0 +
	(double) nkeys / (PARAM_HEAP_RATIO * dbl_ILLutil_our_log2 (nkeys)));
    return dbl_ILLheap_build (&(pinf->h), nkeys, keylist);
}

int dbl_ILLprice_test_for_heap (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const nkeys,
      double *keylist,
      int const algo,
      int const upd)
{
    dbl_heap *const h = &(pinf->h);
    int rval = 0;
    double ravg;
    if (upd != 0) {
	dbl_EGlpNumInitVar (ravg);
	if (algo == PRIMAL_SIMPLEX)
	    dbl_EGlpNumCopy (ravg, lp->cnts->za_ravg);
	else
	    dbl_EGlpNumCopy (ravg, lp->cnts->y_ravg);
	if (dbl_EGlpNumIsLeq (ravg, pinf->htrigger))
	    pinf->hineff--;
	else {
	    dbl_EGlpNumDivUiTo (ravg, 2U);
	    if (dbl_EGlpNumIsLess (pinf->htrigger, ravg))
		pinf->hineff++;
	}
	dbl_EGlpNumClearVar (ravg);
    }
    if (h->hexist == 0 && pinf->hineff <= 0) {
	rval = dbl_ILLprice_build_heap (pinf, nkeys, keylist);
	ILL_CLEANUP_IF (rval);
    } else if (h->hexist != 0 && pinf->hineff >= PARAM_HEAP_UTRIGGER) {
	dbl_ILLprice_free_heap (pinf);
	/*
         * printf ("freeing dbl_heap ..\n");
         * printf ("iter = %d, ravg = %.2f, trigger = %.2f\n",
         * lp->cnts->tot_iter, ravg, pinf->htrigger);
         */
    }
CLEANUP:
    if (rval)
	dbl_ILLprice_free_heap (pinf);
    return rval;
}

void dbl_ILLprice_init_pricing_info (dbl_price_info * const pinf)
{
    pinf->p_strategy = -1;
    pinf->d_strategy = -1;
    pinf->pI_price = -1;
    pinf->pII_price = -1;
    pinf->dI_price = -1;
    pinf->dII_price = -1;
    pinf->cur_price = -1;
    pinf->p_scaleinf = 0;
    pinf->d_scaleinf = 0;
    pinf->pdinfo.norms = 0;
    pinf->pdinfo.refframe = 0;
    pinf->psinfo.norms = 0;
    pinf->ddinfo.norms = 0;
    pinf->ddinfo.refframe = 0;
    pinf->dsinfo.norms = 0;
    pinf->dmpinfo.gstart = pinf->pmpinfo.gstart = 0;
    pinf->dmpinfo.gshift = pinf->pmpinfo.gshift = 0;
    pinf->dmpinfo.gsize = pinf->pmpinfo.gsize = 0;
    pinf->dmpinfo.bucket = pinf->pmpinfo.bucket = 0;
    pinf->dmpinfo.perm = pinf->pmpinfo.perm = 0;
    pinf->dmpinfo.infeas = pinf->pmpinfo.infeas = 0;
    dbl_ILLheap_init (&(pinf->h));
    dbl_EGlpNumZero (pinf->htrigger);
    pinf->hineff = 0;
}

void dbl_ILLprice_free_pricing_info (dbl_price_info * const pinf)
{
    dbl_EGlpNumFreeArray (pinf->p_scaleinf);
    dbl_EGlpNumFreeArray (pinf->d_scaleinf);
    dbl_EGlpNumFreeArray (pinf->pdinfo.norms);
    ILL_IFFREE (pinf->pdinfo.refframe, int);
    dbl_EGlpNumFreeArray (pinf->psinfo.norms);
    dbl_EGlpNumFreeArray (pinf->ddinfo.norms);
    ILL_IFFREE (pinf->ddinfo.refframe, int);
    dbl_EGlpNumFreeArray (pinf->dsinfo.norms);
    dbl_ILLprice_free_mpartial_info (&(pinf->pmpinfo));
    dbl_ILLprice_free_mpartial_info (&(pinf->dmpinfo));
    dbl_ILLprice_free_heap (pinf);
}

int dbl_ILLprice_build_pricing_info (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const phase)
{
    int rval = 0;
    int p_price = -1;
    int d_price = -1;

    switch (phase) {
    case PRIMAL_PHASEI:
	p_price = pinf->pI_price;
	break;
    case PRIMAL_PHASEII:
	p_price = pinf->pII_price;
	break;
    case DUAL_PHASEI:
	d_price = pinf->dI_price;
	break;
    case DUAL_PHASEII:
	d_price = pinf->dII_price;
	break;
    }

    if (p_price != -1) {
	pinf->cur_price = p_price;

	if (p_price == QS_PRICE_PDANTZIG || p_price == QS_PRICE_PDEVEX ||
	    p_price == QS_PRICE_PSTEEP) {
	    pinf->p_strategy = COMPLETE_PRICING;
	    dbl_EGlpNumFreeArray (pinf->d_scaleinf);
	    pinf->d_scaleinf = dbl_EGlpNumAllocArray (lp->nnbasic);
	} else if (p_price == QS_PRICE_PMULTPARTIAL)
	    pinf->p_strategy = MULTI_PART_PRICING;

	switch (p_price) {
	case QS_PRICE_PDEVEX:
	    if (pinf->pdinfo.norms)
		return rval;
	    rval = dbl_ILLprice_build_pdevex_norms (lp, &(pinf->pdinfo), 0);
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_PSTEEP:
	    if (pinf->psinfo.norms)
		return rval;
	    rval = dbl_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_PMULTPARTIAL:
	    rval = dbl_ILLprice_build_mpartial_info (lp, pinf, COL_PRICING);
	    ILL_CLEANUP_IF (rval);
	    break;
	}
    } else if (d_price != -1) {
	pinf->cur_price = d_price;

	if (d_price == QS_PRICE_DDANTZIG || d_price == QS_PRICE_DSTEEP ||
	    d_price == QS_PRICE_DDEVEX) {
	    pinf->d_strategy = COMPLETE_PRICING;
	    dbl_EGlpNumFreeArray (pinf->p_scaleinf);
	    pinf->p_scaleinf = dbl_EGlpNumAllocArray (lp->nrows);
	} else if (d_price == QS_PRICE_DMULTPARTIAL)
	    pinf->d_strategy = MULTI_PART_PRICING;

	switch (d_price) {
	case QS_PRICE_DSTEEP:
	    if (pinf->dsinfo.norms)
		return rval;
	    rval = dbl_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_DMULTPARTIAL:
	    rval = dbl_ILLprice_build_mpartial_info (lp, pinf, ROW_PRICING);
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_DDEVEX:
	    if (pinf->ddinfo.norms)
		return rval;
	    rval = dbl_ILLprice_build_ddevex_norms (lp, &(pinf->ddinfo), 0);
	    ILL_CLEANUP_IF (rval);
	    break;
	}
    }
CLEANUP:
    if (rval)
	dbl_ILLprice_free_pricing_info (pinf);
    ILL_RETURN (rval, "dbl_ILLprice_build_pricing_info");
}

int dbl_ILLprice_update_pricing_info (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const phase,
      dbl_svector * const wz,
      int const eindex,
      int const lindex,
      double y)
{
    int rval = 0;
    int p_price = -1;
    int d_price = -1;

    switch (phase) {
    case PRIMAL_PHASEI:
	p_price = pinf->pI_price;
	break;
    case PRIMAL_PHASEII:
	p_price = pinf->pII_price;
	break;
    case DUAL_PHASEI:
	d_price = pinf->dI_price;
	break;
    case DUAL_PHASEII:
	d_price = pinf->dII_price;
	break;
    }

    if (p_price != -1) {
	if (p_price == QS_PRICE_PDEVEX) {
	    rval = dbl_ILLprice_update_pdevex_norms (lp, &(pinf->pdinfo), eindex, y);
	    ILL_CLEANUP_IF (rval);
	} else if (p_price == QS_PRICE_PSTEEP)
	    dbl_ILLprice_update_psteep_norms (lp, &(pinf->psinfo), wz, eindex, y);
    } else if (d_price != -1) {
	if (d_price == QS_PRICE_DSTEEP)
	    dbl_ILLprice_update_dsteep_norms (lp, &(pinf->dsinfo), wz, lindex, y);
	else if (d_price == QS_PRICE_DDEVEX) {
	    rval = dbl_ILLprice_update_ddevex_norms (lp, &(pinf->ddinfo), lindex, y);
	    ILL_CLEANUP_IF (rval);
	}
    }
CLEANUP:
    ILL_RETURN (rval, "dbl_ILLprice_update_pricing_info");
}

int dbl_ILLprice_get_price (dbl_price_info * const p,
      int const phase)
{
    int pri = -1;

    switch (phase) {
    case PRIMAL_PHASEI:
	return p->pI_price;
    case PRIMAL_PHASEII:
	return p->pII_price;
    case DUAL_PHASEI:
	return p->dI_price;
    case DUAL_PHASEII:
	return p->dII_price;
    }
    return pri;
}

void dbl_ILLprice_free_mpartial_info (dbl_mpart_info * p)
{
    ILL_IFFREE (p->gstart, int);
    ILL_IFFREE (p->gshift, int);
    ILL_IFFREE (p->gsize, int);
    ILL_IFFREE (p->bucket, int);
    dbl_EGlpNumFreeArray (p->infeas);
    ILL_IFFREE (p->perm, int);
}

int dbl_ILLprice_build_mpartial_info (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const pricetype)
{
    dbl_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int const nelems = (pricetype == COL_PRICING) ? lp->nnbasic : lp->nrows;
    int const extra = (nelems % 50) ? nelems - 50 * (nelems / 50) : 0;
    int i = 0;
    int rval = 0;
    p->k = 50;
    p->cgroup = 0;
    p->ngroups = nelems / p->k;
    if (extra != 0)
	p->ngroups++;

    ILL_SAFE_MALLOC (p->gstart, p->ngroups, int);
    ILL_SAFE_MALLOC (p->gshift, p->ngroups, int);
    ILL_SAFE_MALLOC (p->gsize, p->ngroups, int);
    ILL_SAFE_MALLOC (p->bucket, 2 * p->k, int);
    p->infeas = dbl_EGlpNumAllocArray (2 * p->k);
    ILL_SAFE_MALLOC (p->perm, 2 * p->k, int);

    p->bsize = 0;

    if (extra != 0) {
	p->gstart[0] = 0;
	p->gshift[0] = 1;
	p->gsize[0] = extra;
	for (i = 1; i < p->ngroups; i++) {
	    p->gstart[i] = extra + i - 1;
	    p->gshift[i] = p->ngroups - 1;
	    p->gsize[i] = p->k;
	}
    } else {
	for (i = 0; i < p->ngroups; i++) {
	    p->gstart[i] = i;
	    p->gshift[i] = p->ngroups;
	    p->gsize[i] = p->k;
	}
    }

CLEANUP:
    if (rval)
	dbl_ILLprice_free_mpartial_info (p);
    ILL_RETURN (rval, "dbl_ILLprice_build_mpartial_info");
}

void dbl_ILLprice_init_mpartial_price (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const phase,
      int const pricetype)
{
    dbl_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i;
    p->bsize = 0;
    i = p->cgroup;
    do {
	dbl_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
	i = (i + 1) % p->ngroups;
    } while (i != p->cgroup && p->bsize <= p->k);
    p->cgroup = i;
}

void dbl_ILLprice_update_mpartial_price (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const phase,
      int const pricetype)
{
    dbl_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i = 0;
    int csize = 0;
    double infeas;
    dbl_price_res pr;
    dbl_EGlpNumInitVar (pr.dinfeas);
    dbl_EGlpNumInitVar (pr.pinfeas);
    dbl_EGlpNumInitVar (infeas);

#ifdef dbl_MULTIP
    i = 0;
    while (i < p->bsize) {
	if (pricetype == COL_PRICING) {
	    dbl_ILLprice_column (lp, p->bucket[i], phase, &pr);
	    dbl_EGlpNumCopy (infeas, pr.dinfeas);
	} else {
	    dbl_ILLprice_row (lp, p->bucket[i], phase, &pr);
	    dbl_EGlpNumCopy (infeas, pr.pinfeas);
	}
	if (dbl_EGlpNumIsEqqual (infeas, dbl_zeroLpNum)) {
	    p->bucket[i] = p->bucket[p->bsize - 1];
	    p->bsize--;
	} else {
	    dbl_EGlpNumCopy (p->infeas[i], infeas);
	    i++;
	}
    }
    if (p->bsize > 0) {
	for (i = 0; i < p->bsize; i++)
	    p->perm[i] = i;
	dbl_EGutilPermSort ((size_t) (p->bsize), p->perm, (const double *const) p->infeas);

	csize = MIN (p->bsize, p->k);
	for (i = csize - 1; i >= 0; i--)
	    lp->iwork[p->bucket[p->perm[i]]] = 1;

	for (i = 0, csize = 0; i < p->bsize; i++)
	    if (lp->iwork[p->bucket[i]] == 1) {
		dbl_EGlpNumCopy (p->infeas[csize], p->infeas[i]);
		p->bucket[csize] = p->bucket[i];
		csize++;
	    }
	p->bsize = csize;
    }
#else
    p->bsize = 0;
#endif

    i = p->cgroup;
    do {
	dbl_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
	i = (i + 1) % p->ngroups;
    } while (i != p->cgroup && p->bsize <= p->k);
    p->cgroup = i;

#ifdef dbl_MULTIP
    for (i = 0; i < csize; i++)
	lp->iwork[p->bucket[i]] = 0;
#endif
    dbl_EGlpNumClearVar (infeas);
    dbl_EGlpNumClearVar (pr.pinfeas);
    dbl_EGlpNumClearVar (pr.dinfeas);
}

void dbl_ILLprice_delete_onempart_price (dbl_price_info * const pinf,
      int const indx, int const pricetype)
{
    dbl_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i = 0;

    for (i = 0; i < p->bsize; i++)
	if (p->bucket[i] == indx) {
	    p->bucket[i] = p->bucket[p->bsize - 1];
	    dbl_EGlpNumCopy (p->infeas[i], p->infeas[p->bsize - 1]);
	    p->bsize--;
	    break;
	}
}

void dbl_ILLprice_mpartial_group (dbl_lpinfo * const lp,
      dbl_mpart_info * const p,
      int const phase,
      int const g,
      int const pricetype)
{
    int i, ix;
    int gstart = p->gstart[g];
    int gsize = p->gsize[g];
    int gshift = p->gshift[g];
    double infeas;
    dbl_price_res pr;
    dbl_EGlpNumInitVar (pr.dinfeas);
    dbl_EGlpNumInitVar (pr.pinfeas);
    dbl_EGlpNumInitVar (infeas);

    for (i = 0, ix = gstart; i < gsize; i++, ix += gshift) {
#ifdef dbl_MULTIP
	if (lp->iwork[ix])
	    continue;
#endif
	if (pricetype == COL_PRICING) {
	    dbl_ILLprice_column (lp, ix, phase, &pr);
	    dbl_EGlpNumCopy (infeas, pr.dinfeas);
	} else {
	    dbl_ILLprice_row (lp, ix, phase, &pr);
	    dbl_EGlpNumCopy (infeas, pr.pinfeas);
	}
	if (dbl_EGlpNumIsNeqqZero (infeas)) {
	    dbl_EGlpNumCopy (p->infeas[p->bsize], infeas);
	    p->bucket[p->bsize] = ix;
	    p->bsize++;
	}
    }
    dbl_EGlpNumClearVar (infeas);
    dbl_EGlpNumClearVar (pr.dinfeas);
    dbl_EGlpNumClearVar (pr.pinfeas);
}

void dbl_ILLprice_column (dbl_lpinfo * const lp,
      int const ix,
      int const phase,
      dbl_price_res * const pr)
{
    int const col = lp->nbaz[ix];
    int const mcnt = lp->matcnt[col];
    int const mbeg = lp->matbeg[col];
    int i;
    double sum;
    dbl_EGlpNumZero (pr->dinfeas);
    if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	return;
    dbl_EGlpNumInitVar (sum);
    dbl_EGlpNumZero (sum);

    if (phase == PRIMAL_PHASEII) {
	for (i = 0; i < mcnt; i++)
	    dbl_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	dbl_EGlpNumCopyDiff (lp->dz[ix], lp->cz[col], sum);
	dbl_compute_dualII_inf (lp, ix, &(pr->dinfeas));
    } else {
	for (i = 0; i < mcnt; i++)
	    dbl_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	dbl_EGlpNumCopyNeg (lp->pIdz[ix], sum);
	dbl_compute_dualI_inf (lp, ix, &(pr->dinfeas));
    }
    dbl_EGlpNumClearVar (sum);
}

void dbl_ILLprice_row (dbl_lpinfo * const lp,
      int const ix,
      int const phase,
      dbl_price_res * const pr)
{
    if (phase == DUAL_PHASEII)
	dbl_compute_primalII_inf (lp, ix, &(pr->pinfeas));
    else
	dbl_compute_primalI_inf (lp, ix, &(pr->pinfeas));
}

int dbl_ILLprice_build_pdevex_norms (dbl_lpinfo * const lp,
      dbl_p_devex_info * const pdinfo,
      int const reinit)
{
    int j;
    int rval = 0;

    if (reinit == 0) {
	pdinfo->ninit = 0;
	pdinfo->norms = dbl_EGlpNumAllocArray (lp->nnbasic);
	ILL_SAFE_MALLOC (pdinfo->refframe, lp->ncols, int);
    }
    if (reinit != 0)
	pdinfo->ninit++;

    for (j = 0; j < lp->ncols; j++) {
	if (lp->vstat[j] == STAT_BASIC)
	    pdinfo->refframe[j] = 0;
	else {
	    dbl_EGlpNumOne (pdinfo->norms[lp->vindex[j]]);
	    pdinfo->refframe[j] = 1;
	}
    }

CLEANUP:
    if (rval) {
	dbl_EGlpNumFreeArray (pdinfo->norms);
	ILL_IFFREE (pdinfo->refframe, int);
    }
    ILL_RETURN (rval, "dbl_ILLprice_build_pdevex_norms");
}

int dbl_ILLprice_update_pdevex_norms (dbl_lpinfo * const lp,
      dbl_p_devex_info * const pdinfo,
      int const eindex,
      double yl)
{
    int i, j;
    double normj;
    double zAj;
    dbl_EGlpNumInitVar (normj);
    dbl_EGlpNumInitVar (zAj);
    dbl_EGlpNumZero (normj);

    for (i = 0; i < lp->yjz.nzcnt; i++)
	if (pdinfo->refframe[lp->baz[lp->yjz.indx[i]]])
	    dbl_EGlpNumAddInnProdTo (normj, lp->yjz.coef[i], lp->yjz.coef[i]);

    if (pdinfo->refframe[lp->nbaz[eindex]])
	dbl_EGlpNumAddTo (normj, dbl_oneLpNum);

    dbl_EGlpNumCopyFrac (zAj, normj, pdinfo->norms[eindex]);
    if (dbl_EGlpNumIsGreaDbl (zAj, 0x1p10) || dbl_EGlpNumIsLessDbl (zAj, 0x1p-10)) {
	dbl_EGlpNumClearVar (zAj);
	dbl_EGlpNumClearVar (normj);
	return dbl_ILLprice_build_pdevex_norms (lp, pdinfo, 1);
    }
    for (i = 0; i < lp->zA.nzcnt; i++) {
	j = lp->zA.indx[i];
	dbl_EGlpNumCopyFrac (zAj, lp->zA.coef[i], yl);
	dbl_EGlpNumMultTo (zAj, zAj);
	dbl_EGlpNumMultTo (zAj, normj);
	if (dbl_EGlpNumIsLess (pdinfo->norms[j], zAj))
	    dbl_EGlpNumCopy (pdinfo->norms[j], zAj);
    }
    dbl_EGlpNumDivTo (normj, yl);
    dbl_EGlpNumDivTo (normj, yl);
    if (dbl_EGlpNumIsLess (normj, dbl_oneLpNum))
	dbl_EGlpNumCopy (pdinfo->norms[eindex], dbl_oneLpNum);
    else
	dbl_EGlpNumCopy (pdinfo->norms[eindex], normj);
    dbl_EGlpNumClearVar (zAj);
    dbl_EGlpNumClearVar (normj);
    return 0;
}

int dbl_ILLprice_build_psteep_norms (dbl_lpinfo * const lp,
      dbl_p_steep_info * const psinfo)
{
    int j;
    int rval = 0;
    dbl_svector yz;

    dbl_ILLsvector_init (&yz);
    rval = dbl_ILLsvector_alloc (&yz, lp->nrows);
    ILL_CLEANUP_IF (rval);
    psinfo->norms = dbl_EGlpNumAllocArray (lp->nnbasic);

    for (j = 0; j < lp->nnbasic; j++) {
	rval = ILLstring_report (NULL, &lp->O->reporter);
	ILL_CLEANUP_IF (rval);
	dbl_ILLfct_compute_yz (lp, &yz, 0, lp->nbaz[j]);
	dbl_EGlpNumInnProd (psinfo->norms[j], yz.coef, yz.coef, (size_t) yz.nzcnt);
	dbl_EGlpNumAddTo (psinfo->norms[j], dbl_oneLpNum);
    }

CLEANUP:
    dbl_ILLsvector_free (&yz);
    if (rval)
	dbl_EGlpNumFreeArray (psinfo->norms);
    ILL_RETURN (rval, "dbl_ILLprice_build_psteep_norms");
}

void dbl_ILLprice_update_psteep_norms (dbl_lpinfo * const lp,
      dbl_p_steep_info * const psinfo,
      dbl_svector * const wz,
      int const eindex,
      double yl)
{
    int i, j, k;
    int mcnt, mbeg;
    double normj;
    double zAj, wAj;
    double *v = 0;
    dbl_EGlpNumInitVar (normj);
    dbl_EGlpNumInitVar (zAj);
    dbl_EGlpNumInitVar (wAj);
    dbl_EGlpNumInnProd (normj, lp->yjz.coef, lp->yjz.coef, (size_t) (lp->yjz.nzcnt));
    dbl_EGlpNumAddTo (normj, dbl_oneLpNum);

#if 0
    Bico - remove warnings for dist
	if (fabs ((normj - psinfo->norms[eindex]) / normj) > 1000.0 /* 0.01 */ ) {
	    printf ("warning: incorrect norm values\n");
	    printf ("anorm = %.6f, pnorm = %.6f\n", normj, psinfo->norms[eindex]);
	    fflush (stdout);
	}
#endif

    dbl_ILLfct_load_workvector (lp, wz);
    v = lp->work.coef;

    for (k = 0; k < lp->zA.nzcnt; k++) {
	j = lp->zA.indx[k];
	dbl_EGlpNumZero (wAj);
	mcnt = lp->matcnt[lp->nbaz[j]];
	mbeg = lp->matbeg[lp->nbaz[j]];
	for (i = 0; i < mcnt; i++)
	    dbl_EGlpNumAddInnProdTo (wAj, lp->matval[mbeg + i], v[lp->matind[mbeg + i]]);
	dbl_EGlpNumCopy (zAj, lp->zA.coef[k]);
	/* psinfo->norms[j] += (zAj * ((zAj * normj / yl) - (2.0 * wAj))) /
	   yl; */
	dbl_EGlpNumMultUiTo (wAj, 2U);
	dbl_EGlpNumMultTo (zAj, normj);
	dbl_EGlpNumDivTo (zAj, yl);
	dbl_EGlpNumSubTo (zAj, wAj);
	dbl_EGlpNumMultTo (zAj, lp->zA.coef[k]);
	dbl_EGlpNumDivTo (zAj, yl);
	dbl_EGlpNumAddTo (psinfo->norms[j], zAj);
	if (dbl_EGlpNumIsLess (psinfo->norms[j], dbl_oneLpNum))
	    dbl_EGlpNumOne (psinfo->norms[j]);
    }

    dbl_EGlpNumCopyFrac (psinfo->norms[eindex], normj, yl);
    dbl_EGlpNumDivTo (psinfo->norms[eindex], yl);
    if (dbl_EGlpNumIsLess (psinfo->norms[eindex], dbl_oneLpNum))
	dbl_EGlpNumOne (psinfo->norms[eindex]);

    dbl_ILLfct_zero_workvector (lp);
    dbl_EGlpNumClearVar (wAj);
    dbl_EGlpNumClearVar (zAj);
    dbl_EGlpNumClearVar (normj);
}

int dbl_ILLprice_build_ddevex_norms (dbl_lpinfo * const lp,
      dbl_d_devex_info * const ddinfo,
      int const reinit)
{
    int i;
    int rval = 0;

    if (reinit == 0) {
	ddinfo->ninit = 0;
	ddinfo->norms = dbl_EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (ddinfo->refframe, lp->ncols, int);
    }
    if (reinit != 0)
	ddinfo->ninit++;

    for (i = 0; i < lp->ncols; i++)
	ddinfo->refframe[i] = (lp->vstat[i] == STAT_BASIC) ? 1 : 0;

    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumOne (ddinfo->norms[i]);

CLEANUP:
    if (rval) {
	dbl_EGlpNumFreeArray (ddinfo->norms);
	ILL_IFFREE (ddinfo->refframe, int);
    }
    ILL_RETURN (rval, "dbl_ILLprice_build_ddevex_norms");
}

int dbl_ILLprice_update_ddevex_norms (dbl_lpinfo * const lp,
      dbl_d_devex_info * const ddinfo,
      int const lindex,
      double yl)
{
    int i, r;
    double normi;
    double yr;
    dbl_EGlpNumInitVar (normi);
    dbl_EGlpNumInitVar (yr);
    dbl_EGlpNumZero (normi);

    for (i = 0; i < lp->zA.nzcnt; i++)
	if (ddinfo->refframe[lp->nbaz[lp->zA.indx[i]]])
	    dbl_EGlpNumAddInnProdTo (normi, lp->zA.coef[i], lp->zA.coef[i]);

    if (ddinfo->refframe[lp->baz[lindex]])
	dbl_EGlpNumAddTo (normi, dbl_oneLpNum);

    dbl_EGlpNumCopyFrac (yr, normi, ddinfo->norms[lindex]);
    if (dbl_EGlpNumIsGreaDbl (yr, 0x1p10) || dbl_EGlpNumIsLessDbl (yr, 0x1p-10)) {
	dbl_EGlpNumClearVar (normi);
	dbl_EGlpNumClearVar (yr);
	return dbl_ILLprice_build_ddevex_norms (lp, ddinfo, 1);
    }
    dbl_EGlpNumDivTo (normi, yl);
    dbl_EGlpNumDivTo (normi, yl);

    for (i = 0; i < lp->yjz.nzcnt; i++) {
	r = lp->yjz.indx[i];
	dbl_EGlpNumCopy (yr, lp->yjz.coef[i]);
	dbl_EGlpNumMultTo (yr, yr);
	dbl_EGlpNumMultTo (yr, normi);
	if (dbl_EGlpNumIsLess (ddinfo->norms[r], yr))
	    dbl_EGlpNumCopy (ddinfo->norms[r], yr);
    }
    dbl_EGlpNumCopy (ddinfo->norms[lindex], normi);
    if (dbl_EGlpNumIsLess (ddinfo->norms[lindex], dbl_oneLpNum))
	dbl_EGlpNumOne (ddinfo->norms[lindex]);
    dbl_EGlpNumClearVar (normi);
    dbl_EGlpNumClearVar (yr);
    return 0;
}

int dbl_ILLprice_build_dsteep_norms (dbl_lpinfo * const lp,
      dbl_d_steep_info * const dsinfo)
{
    int i;
    int rval = 0;
    dbl_svector z;

    dbl_ILLsvector_init (&z);
    rval = dbl_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);
    dsinfo->norms = dbl_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++) {
	rval = ILLstring_report (NULL, &lp->O->reporter);
	ILL_CLEANUP_IF (rval);

	dbl_ILLfct_compute_zz (lp, &z, i);
	dbl_EGlpNumInnProd (dsinfo->norms[i], z.coef, z.coef, (size_t) z.nzcnt);
	if (dbl_EGlpNumIsLess (dsinfo->norms[i], dbl_PARAM_MIN_DNORM))
	    dbl_EGlpNumCopy (dsinfo->norms[i], dbl_PARAM_MIN_DNORM);
    }

CLEANUP:
    dbl_ILLsvector_free (&z);
    if (rval)
	dbl_EGlpNumFreeArray (dsinfo->norms);
    ILL_RETURN (rval, "dbl_ILLprice_build_dsteep_norms");
}

int dbl_ILLprice_get_dsteep_norms (dbl_lpinfo * const lp,
      int const count,
      int *const rowind,
      double *const norms)
{
    int i;
    int rval = 0;
    dbl_svector z;

    dbl_ILLsvector_init (&z);
    rval = dbl_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < count; i++) {
	dbl_ILLfct_compute_zz (lp, &z, rowind[i]);
	dbl_EGlpNumInnProd (norms[i], z.coef, z.coef, (size_t) z.nzcnt);
    }

CLEANUP:
    dbl_ILLsvector_free (&z);
    ILL_RETURN (rval, "dbl_ILLprice_get_dsteep_norms");
}

void dbl_ILLprice_update_dsteep_norms (dbl_lpinfo * const lp,
      dbl_d_steep_info * const dsinfo,
      dbl_svector * const wz,
      int const lindex,
      double yl)
{
    int i, k;
    double yij;
    double norml;
    double *v = 0;
    dbl_EGlpNumInitVar (norml);
    dbl_EGlpNumInitVar (yij);
    dbl_EGlpNumZero (norml);
    dbl_EGlpNumInnProd (norml, lp->zz.coef, lp->zz.coef, (size_t) (lp->zz.nzcnt));

#if 0
    Bico - remove warnings for dist
	if (fabs ((norml - dsinfo->norms[lindex]) / norml) > 1000.0 /* 0.01 */ ) {
	    printf ("warning: incorrect dnorm values\n");
	    printf ("anorm = %.6f, pnorm = %.6f\n", norml, dsinfo->norms[lindex]);
	    fflush (stdout);
	}
#endif

    dbl_ILLfct_load_workvector (lp, wz);
    v = lp->work.coef;

    for (k = 0; k < lp->yjz.nzcnt; k++) {
	i = lp->yjz.indx[k];
	dbl_EGlpNumCopy (yij, lp->yjz.coef[k]);
	dbl_EGlpNumMultTo (yij, norml);
	dbl_EGlpNumDivTo (yij, yl);
	dbl_EGlpNumSubTo (yij, v[i]);
	dbl_EGlpNumSubTo (yij, v[i]);
	dbl_EGlpNumMultTo (yij, lp->yjz.coef[k]);
	dbl_EGlpNumDivTo (yij, yl);
	dbl_EGlpNumAddTo (dsinfo->norms[i], yij);
	if (dbl_EGlpNumIsLess (dsinfo->norms[i], dbl_PARAM_MIN_DNORM))
	    dbl_EGlpNumCopy (dsinfo->norms[i], dbl_PARAM_MIN_DNORM);
    }
    dbl_EGlpNumCopyFrac (dsinfo->norms[lindex], norml, yl);
    dbl_EGlpNumDivTo (dsinfo->norms[lindex], yl);
    if (dbl_EGlpNumIsLess (dsinfo->norms[lindex], dbl_PARAM_MIN_DNORM))
	dbl_EGlpNumCopy (dsinfo->norms[lindex], dbl_PARAM_MIN_DNORM);

    dbl_ILLfct_zero_workvector (lp);
    dbl_EGlpNumClearVar (norml);
    dbl_EGlpNumClearVar (yij);
}

static void dbl_update_d_scaleinf (dbl_price_info * const p,
      dbl_heap * const h,
      int const j,
      double inf,
      int const prule)
{
    if (dbl_EGlpNumIsEqqual (inf, dbl_zeroLpNum)) {
	dbl_EGlpNumZero (p->d_scaleinf[j]);
	if (h->hexist != 0 && h->loc[j] != -1)
	    dbl_ILLheap_delete (h, j);
    } else {
	if (prule == QS_PRICE_PDANTZIG)
	    dbl_EGlpNumCopy (p->d_scaleinf[j], inf);
	else if (prule == QS_PRICE_PDEVEX)
	    dbl_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->pdinfo.norms[j]);
	else if (prule == QS_PRICE_PSTEEP)
	    dbl_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->psinfo.norms[j]);

	if (h->hexist != 0) {
	    if (h->loc[j] == -1)
		dbl_ILLheap_insert (h, j);
	    else
		dbl_ILLheap_modify (h, j);
	}
    }
}

static void dbl_compute_dualI_inf (dbl_lpinfo * const lp,
      const int j,
      double *const inf)
{
    int const col = lp->nbaz[j];
    int const vt = lp->vtype[col];
    int const vs = lp->vstat[col];
    dbl_EGlpNumZero (*inf);
    if (vt != VARTIFICIAL && vt != VFIXED) {
	if ((vs == STAT_LOWER || vs == STAT_ZERO) &&
	    dbl_EGlpNumIsSumLess (lp->pIdz[j], lp->tol->id_tol, dbl_zeroLpNum))
	    dbl_EGlpNumCopyAbs (*inf, lp->pIdz[j]);
	else if ((vs == STAT_UPPER || vs == STAT_ZERO) &&
	    dbl_EGlpNumIsLess (lp->tol->id_tol, lp->pIdz[j]))
	    dbl_EGlpNumCopy (*inf, lp->pIdz[j]);
    }
}

static void dbl_compute_dualII_inf (dbl_lpinfo * const lp,
      int const j,
      double *const inf)
{
    int const col = lp->nbaz[j];
    int const vt = lp->vtype[col];
    int const vs = lp->vstat[col];
    dbl_EGlpNumZero (*inf);
    if (vt != VARTIFICIAL && vt != VFIXED) {
	if ((vs == STAT_LOWER || vs == STAT_ZERO) &&
	    dbl_EGlpNumIsSumLess (lp->dz[j], lp->tol->dfeas_tol, dbl_zeroLpNum))
	    dbl_EGlpNumCopyAbs (*inf, lp->dz[j]);
	else if ((vs == STAT_UPPER || vs == STAT_ZERO) &&
	    dbl_EGlpNumIsLess (lp->tol->dfeas_tol, lp->dz[j]))
	    dbl_EGlpNumCopy (*inf, lp->dz[j]);
    }
}

void dbl_ILLprice_compute_dual_inf (dbl_lpinfo * const lp,
      dbl_price_info * const p,
      int *const ix,
      int const icnt,
      int const phase)
{
    int const price = (phase == PRIMAL_PHASEI) ? p->pI_price : p->pII_price;
    dbl_heap *const h = &(p->h);
    int i;
    double inf;
    dbl_EGlpNumInitVar (inf);
    dbl_EGlpNumZero (inf);

    if (phase == PRIMAL_PHASEI) {
	if (ix == NULL)
	    for (i = 0; i < lp->nnbasic; i++) {
		dbl_compute_dualI_inf (lp, i, &(inf));
		dbl_update_d_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		dbl_compute_dualI_inf (lp, ix[i], &(inf));
		dbl_update_d_scaleinf (p, h, ix[i], inf, price);
	    }
    } else if (phase == PRIMAL_PHASEII) {
	if (ix == NULL)
	    for (i = 0; i < lp->nnbasic; i++) {
		dbl_compute_dualII_inf (lp, i, &inf);
		dbl_update_d_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		dbl_compute_dualII_inf (lp, ix[i], &inf);
		dbl_update_d_scaleinf (p, h, ix[i], inf, price);
	    }
    }
    dbl_EGlpNumClearVar (inf);
}

void dbl_ILLprice_primal (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      dbl_price_res * const pr,
      int const phase)
{
    dbl_heap *const h = &(pinf->h);
    int j, vs;
    double d;
    dbl_EGlpNumInitVar (d);
    dbl_EGlpNumZero (d);
    pr->eindex = -1;

#if USEHEAP > 0
    dbl_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
	PRIMAL_SIMPLEX, 1);
#endif

    if (pinf->p_strategy == COMPLETE_PRICING) {
	if (h->hexist) {
	    pr->eindex = dbl_ILLheap_findmin (h);
	    if (pr->eindex != -1)
		dbl_ILLheap_delete (h, pr->eindex);
	} else {
	    for (j = 0; j < lp->nnbasic; j++) {
		if (dbl_EGlpNumIsLess (d, pinf->d_scaleinf[j])) {
		    dbl_EGlpNumCopy (d, pinf->d_scaleinf[j]);
		    pr->eindex = j;
		}
	    }
	}
    } else if (pinf->p_strategy == MULTI_PART_PRICING) {
	for (j = 0; j < pinf->pmpinfo.bsize; j++) {
	    if (dbl_EGlpNumIsLess (d, pinf->pmpinfo.infeas[j])) {
		dbl_EGlpNumCopy (d, pinf->pmpinfo.infeas[j]);
		pr->eindex = pinf->pmpinfo.bucket[j];
	    }
	}
    }
    if (pr->eindex < 0)
	pr->price_stat = PRICE_OPTIMAL;
    else {
	if (phase == PRIMAL_PHASEI)
	    dbl_EGlpNumCopy (d, lp->pIdz[pr->eindex]);
	else
	    dbl_EGlpNumCopy (d, lp->dz[pr->eindex]);
	vs = lp->vstat[lp->nbaz[pr->eindex]];

	pr->price_stat = PRICE_NONOPTIMAL;
	if (vs == STAT_UPPER || (vs == STAT_ZERO &&
		dbl_EGlpNumIsLess (lp->tol->dfeas_tol, d)))
	    pr->dir = VDECREASE;
	else
	    pr->dir = VINCREASE;
    }
    dbl_EGlpNumClearVar (d);
}

static void dbl_update_p_scaleinf (dbl_price_info * const p,
      dbl_heap * const h,
      int const i,
      double inf,
      int const prule)
{
    if (dbl_EGlpNumIsEqqual (inf, dbl_zeroLpNum)) {
	dbl_EGlpNumZero (p->p_scaleinf[i]);
	if (h->hexist != 0 && h->loc[i] != -1)
	    dbl_ILLheap_delete (h, i);
    } else {
	if (prule == QS_PRICE_DDANTZIG)
	    dbl_EGlpNumCopy (p->p_scaleinf[i], inf);
	else if (prule == QS_PRICE_DSTEEP)
	    dbl_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->dsinfo.norms[i]);
	else if (prule == QS_PRICE_DDEVEX)
	    dbl_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->ddinfo.norms[i]);

	if (h->hexist != 0) {
	    if (h->loc[i] == -1)
		dbl_ILLheap_insert (h, i);
	    else
		dbl_ILLheap_modify (h, i);
	}
    }
}

static void dbl_compute_primalI_inf (dbl_lpinfo * const lp,
      int const i,
      double *const inf)
{
    int const col = lp->baz[i];
    dbl_EGlpNumZero (*inf);
    if (dbl_EGlpNumIsLess (lp->tol->ip_tol, lp->xbz[i]) &&
	dbl_EGlpNumIsNeqq (lp->uz[col], dbl_INFTY))
	dbl_EGlpNumCopy (*inf, lp->xbz[i]);
    else if (dbl_EGlpNumIsNeqq (lp->lz[col], dbl_NINFTY) &&
	dbl_EGlpNumIsSumLess (lp->xbz[i], lp->tol->ip_tol, dbl_zeroLpNum))
	dbl_EGlpNumCopyAbs (*inf, lp->xbz[i]);
}

static void dbl_compute_primalII_inf (dbl_lpinfo * const lp,
      int const i,
      double *const inf)
{
    int const col = lp->baz[i];
    dbl_EGlpNumZero (*inf);
    if (dbl_EGlpNumIsNeqq (lp->uz[col], dbl_INFTY) &&
	dbl_EGlpNumIsSumLess (lp->uz[col], lp->tol->pfeas_tol, lp->xbz[i]))
	dbl_EGlpNumCopyDiff (*inf, lp->xbz[i], lp->uz[col]);
    else if (dbl_EGlpNumIsNeqq (lp->lz[col], dbl_NINFTY) &&
	dbl_EGlpNumIsSumLess (lp->xbz[i], lp->tol->pfeas_tol, lp->lz[col]))
	dbl_EGlpNumCopyDiff (*inf, lp->lz[col], lp->xbz[i]);
}

void dbl_ILLprice_compute_primal_inf (dbl_lpinfo * const lp,
      dbl_price_info * const p,
      int *const ix,
      int const icnt,
      int const phase)
{
    dbl_heap *const h = &(p->h);
    int const price = (phase == DUAL_PHASEI) ? p->dI_price : p->dII_price;
    int i;
    double inf;
    dbl_EGlpNumInitVar (inf);
    dbl_EGlpNumZero (inf);

    if (phase == DUAL_PHASEI) {
	if (ix == NULL)
	    for (i = 0; i < lp->nrows; i++) {
		dbl_compute_primalI_inf (lp, i, &inf);
		dbl_update_p_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		dbl_compute_primalI_inf (lp, ix[i], &inf);
		dbl_update_p_scaleinf (p, h, ix[i], inf, price);
	    }
    } else if (phase == DUAL_PHASEII) {
	if (ix == NULL)
	    for (i = 0; i < lp->nrows; i++) {
		dbl_compute_primalII_inf (lp, i, &inf);
		dbl_update_p_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		dbl_compute_primalII_inf (lp, ix[i], &inf);
		dbl_update_p_scaleinf (p, h, ix[i], inf, price);
	    }
    }
    dbl_EGlpNumClearVar (inf);
}

void dbl_ILLprice_dual (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      int const phase,
      dbl_price_res * const pr)
{
    dbl_heap *const h = &(pinf->h);
    int i;
    double d;
    dbl_EGlpNumInitVar (d);
    dbl_EGlpNumZero (d);
    pr->lindex = -1;

#if USEHEAP > 0
    dbl_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf, DUAL_SIMPLEX,
	1);
#endif

    if (pinf->d_strategy == COMPLETE_PRICING) {
	if (h->hexist) {
	    pr->lindex = dbl_ILLheap_findmin (h);
	    if (pr->lindex != -1)
		dbl_ILLheap_delete (h, pr->lindex);
	} else {
	    for (i = 0; i < lp->nrows; i++) {
		if (dbl_EGlpNumIsLess (d, pinf->p_scaleinf[i])) {
		    dbl_EGlpNumCopy (d, pinf->p_scaleinf[i]);
		    pr->lindex = i;
		}
	    }
	}
    } else if (pinf->d_strategy == MULTI_PART_PRICING) {
	for (i = 0; i < pinf->dmpinfo.bsize; i++) {
	    if (dbl_EGlpNumIsLess (d, pinf->dmpinfo.infeas[i])) {
		dbl_EGlpNumCopy (d, pinf->dmpinfo.infeas[i]);
		pr->lindex = pinf->dmpinfo.bucket[i];
	    }
	}
    }
    if (pr->lindex < 0)
	pr->price_stat = PRICE_OPTIMAL;
    else {
	pr->price_stat = NONOPTIMAL;

	if (dbl_EGlpNumIsNeqq (lp->uz[lp->baz[pr->lindex]], dbl_INFTY)) {
	    if (phase == DUAL_PHASEI)
		dbl_EGlpNumZero (d);
	    else
		dbl_EGlpNumCopy (d, lp->uz[lp->baz[pr->lindex]]);
	    if (dbl_EGlpNumIsSumLess (lp->tol->pfeas_tol, d, lp->xbz[pr->lindex]))
		pr->lvstat = STAT_UPPER;
	    else
		pr->lvstat = STAT_LOWER;
	} else
	    pr->lvstat = STAT_LOWER;
    }
    dbl_EGlpNumClearVar (d);
}

int dbl_ILLprice_get_rownorms (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      double *const rnorms)
{
    int rval = 0;
    int i;

    if (pinf->dsinfo.norms == NULL) {
	rval = dbl_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
	ILL_CLEANUP_IF (rval);
    }
    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumCopy (rnorms[i], pinf->dsinfo.norms[i]);

CLEANUP:
    if (rval)
	dbl_EGlpNumFreeArray (pinf->dsinfo.norms);
    return rval;
}

int dbl_ILLprice_get_colnorms (dbl_lpinfo * const lp,
      dbl_price_info * const pinf,
      double *const cnorms)
{
    int rval = 0;
    int i, j;

    if (pinf->psinfo.norms == NULL) {
	rval = dbl_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
	ILL_CLEANUP_IF (rval);
    }
    for (i = 0; i < lp->nrows; i++)
	dbl_EGlpNumZero (cnorms[lp->baz[i]]);
    for (j = 0; j < lp->nnbasic; j++)
	dbl_EGlpNumCopy (cnorms[lp->nbaz[j]], pinf->psinfo.norms[j]);

CLEANUP:
    if (rval)
	dbl_EGlpNumFreeArray (pinf->psinfo.norms);
    return rval;
}

int dbl_ILLprice_get_newnorms (dbl_lpinfo * const lp,
      int const nelems,
      double *const norms,
      int *const matcnt,
      int *const matbeg,
      int *const matind,
      double *const matval,
      int const option)
{
    int i, j;
    int rval = 0;
    dbl_svector a;
    dbl_svector y;

    dbl_ILLsvector_init (&y);
    rval = dbl_ILLsvector_alloc (&y, lp->nrows);
    ILL_CLEANUP_IF (rval);

    for (j = 0; j < nelems; j++) {
	a.nzcnt = matcnt[j];
	a.indx = &(matind[matbeg[j]]);
	a.coef = &(matval[matbeg[j]]);

	if (option == COLUMN_SOLVE)
	    dbl_ILLbasis_column_solve (lp, &a, &y);
	else
	    dbl_ILLbasis_row_solve (lp, &a, &y);

	dbl_EGlpNumOne (norms[j]);
	for (i = 0; i < y.nzcnt; i++)
	    dbl_EGlpNumAddInnProdTo (norms[j], y.coef[i], y.coef[i]);
    }

CLEANUP:
    dbl_ILLsvector_free (&y);
    ILL_RETURN (rval, "dbl_ILLprice_get_newnorms");
}

int dbl_ILLprice_get_new_rownorms (dbl_lpinfo * const lp,
      int const newrows,
      double *const rnorms,
      int *const rmatcnt,
      int *const rmatbeg,
      int *const rmatind,
      double *const rmatval)
{
    return dbl_ILLprice_get_newnorms (lp, newrows, rnorms, rmatcnt, rmatbeg, rmatind,
	rmatval, ROW_SOLVE);
}

int dbl_ILLprice_get_new_colnorms (dbl_lpinfo * const lp,
      int const newrows,
      double *const rnorms,
      int *const matcnt,
      int *const matbeg,
      int *const matind,
      double *const matval)
{
    return dbl_ILLprice_get_newnorms (lp, newrows, rnorms, matcnt, matbeg, matind,
	matval, COLUMN_SOLVE);
}

int dbl_ILLprice_load_rownorms (dbl_lpinfo * const lp,
      double *const rnorms,
      dbl_price_info * const pinf)
{
    int i;
    int rval = 0;

    dbl_EGlpNumFreeArray (pinf->dsinfo.norms);
    pinf->dsinfo.norms = dbl_EGlpNumAllocArray (lp->nrows);
    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumCopy (pinf->dsinfo.norms[i], rnorms[i]);
	if (dbl_EGlpNumIsLess (pinf->dsinfo.norms[i], dbl_PARAM_MIN_DNORM))
	    dbl_EGlpNumCopy (pinf->dsinfo.norms[i], dbl_PARAM_MIN_DNORM);
    }

    ILL_RETURN (rval, "dbl_ILLprice_load_rownorms");
}

int dbl_ILLprice_load_colnorms (dbl_lpinfo * const lp,
      double *const cnorms,
      dbl_price_info * const pinf)
{
    int j;
    int rval = 0;

    dbl_EGlpNumFreeArray (pinf->psinfo.norms);
    pinf->psinfo.norms = dbl_EGlpNumAllocArray (lp->nnbasic);
    for (j = 0; j < lp->nnbasic; j++) {
	dbl_EGlpNumCopy (pinf->psinfo.norms[j], cnorms[lp->nbaz[j]]);
	if (dbl_EGlpNumIsLess (pinf->psinfo.norms[j], dbl_oneLpNum))
	    dbl_EGlpNumOne (pinf->psinfo.norms[j]);
    }

    ILL_RETURN (rval, "dbl_ILLprice_load_colnorms");
}

#if dbl_PRICE_DEBUG > 0
void dbl_test_dsteep_norms (dbl_lpinfo * lp, dbl_price_info * p)
{
    int i, errn = 0;
    double err, diff;
    double *pn;
    dbl_EGlpNumInitVar (err);
    dbl_EGlpNumInitVar (diff);
    pn = dbl_EGlpNumAllocArray (lp->nrows);
    dbl_EGlpNumZero (err);

    dbl_ILLprice_get_dsteep_norms (lp, lp->yjz.nzcnt, lp->yjz.indx, pn);
    for (i = 0; i < lp->yjz.nzcnt; i++) {
	dbl_EGlpNumCopyDiff (diff, pn[i], p->dsinfo.norms[lp->yjz.indx[i]]);
	if (dbl_EGlpNumIsLess (diff, dbl_zeroLpNum))
	    dbl_EGlpNumSign (diff);
	if (dbl_EGlpNumIsLess (dbl_PFEAS_TOLER, diff)) {
	    errn++;
	    dbl_EGlpNumAddTo (err, diff);
	    dbl_EGlpNumCopy (p->dsinfo.norms[lp->yjz.indx[i]], pn[i]);
	}
    }
    if (errn)
	printf ("%d: dnorm errn = %d, err = %.6f\n", lp->cnts->tot_iter, errn,
	    dbl_EGlpNumToLf (err));
    dbl_EGlpNumFreeArray (pn);
    dbl_EGlpNumClearVar (diff);
    dbl_EGlpNumClearVar (err);
}
#endif
