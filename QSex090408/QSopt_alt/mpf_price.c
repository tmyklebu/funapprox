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

/* RCS_INFO = "$RCSfile: mpf_price.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;
#include "econfig.h"
#include "stddefs.h"
#include "mpf_qsopt.h"
#include "mpf_lpdefs.h"
#include "mpf_fct.h"
#include "mpf_price.h"
#include "mpf_basis.h"
#include "mpf_iqsutil.h"
#include "mpf_dstruct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define  mpf_MULTIP 1
#define  mpf_PRICE_DEBUG 0

static void mpf_update_d_scaleinf (mpf_price_info * const p,
      mpf_heap * const h,
      int const j,
      mpf_t inf,
      int const prule),
  mpf_update_p_scaleinf (mpf_price_info * const p,
      mpf_heap * const h,
      int const i,
      mpf_t inf,
      int const prule);

static void mpf_compute_dualI_inf (mpf_lpinfo * const lp,
      int const j,
      mpf_t * const inf),
  mpf_compute_dualII_inf (mpf_lpinfo * const lp,
      int const j,
      mpf_t * const inf),
  mpf_compute_primalI_inf (mpf_lpinfo * const lp,
      int const i,
      mpf_t * const inf),
  mpf_compute_primalII_inf (mpf_lpinfo * const lp,
      int const i,
      mpf_t * const inf);

void mpf_ILLprice_free_heap (mpf_price_info * const pinf)
{
    mpf_ILLheap_free (&(pinf->h));
}

int mpf_ILLprice_build_heap (mpf_price_info * const pinf,
      int const nkeys,
      mpf_t * keylist)
{
    mpf_ILLheap_init (&(pinf->h));
    mpf_EGlpNumSet (pinf->htrigger,
	1.0 +
	(double) nkeys / (PARAM_HEAP_RATIO * mpf_ILLutil_our_log2 (nkeys)));
    return mpf_ILLheap_build (&(pinf->h), nkeys, keylist);
}

int mpf_ILLprice_test_for_heap (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      int const nkeys,
      mpf_t * keylist,
      int const algo,
      int const upd)
{
    mpf_heap *const h = &(pinf->h);
    int rval = 0;
    mpf_t ravg;
    if (upd != 0) {
	mpf_EGlpNumInitVar (ravg);
	if (algo == PRIMAL_SIMPLEX)
	    mpf_EGlpNumCopy (ravg, lp->cnts->za_ravg);
	else
	    mpf_EGlpNumCopy (ravg, lp->cnts->y_ravg);
	if (mpf_EGlpNumIsLeq (ravg, pinf->htrigger))
	    pinf->hineff--;
	else {
	    mpf_EGlpNumDivUiTo (ravg, 2U);
	    if (mpf_EGlpNumIsLess (pinf->htrigger, ravg))
		pinf->hineff++;
	}
	mpf_EGlpNumClearVar (ravg);
    }
    if (h->hexist == 0 && pinf->hineff <= 0) {
	rval = mpf_ILLprice_build_heap (pinf, nkeys, keylist);
	ILL_CLEANUP_IF (rval);
    } else if (h->hexist != 0 && pinf->hineff >= PARAM_HEAP_UTRIGGER) {
	mpf_ILLprice_free_heap (pinf);
	/*
         * printf ("freeing mpf_heap ..\n");
         * printf ("iter = %d, ravg = %.2f, trigger = %.2f\n",
         * lp->cnts->tot_iter, ravg, pinf->htrigger);
         */
    }
CLEANUP:
    if (rval)
	mpf_ILLprice_free_heap (pinf);
    return rval;
}

void mpf_ILLprice_init_pricing_info (mpf_price_info * const pinf)
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
    mpf_ILLheap_init (&(pinf->h));
    mpf_EGlpNumZero (pinf->htrigger);
    pinf->hineff = 0;
}

void mpf_ILLprice_free_pricing_info (mpf_price_info * const pinf)
{
    mpf_EGlpNumFreeArray (pinf->p_scaleinf);
    mpf_EGlpNumFreeArray (pinf->d_scaleinf);
    mpf_EGlpNumFreeArray (pinf->pdinfo.norms);
    ILL_IFFREE (pinf->pdinfo.refframe, int);
    mpf_EGlpNumFreeArray (pinf->psinfo.norms);
    mpf_EGlpNumFreeArray (pinf->ddinfo.norms);
    ILL_IFFREE (pinf->ddinfo.refframe, int);
    mpf_EGlpNumFreeArray (pinf->dsinfo.norms);
    mpf_ILLprice_free_mpartial_info (&(pinf->pmpinfo));
    mpf_ILLprice_free_mpartial_info (&(pinf->dmpinfo));
    mpf_ILLprice_free_heap (pinf);
}

int mpf_ILLprice_build_pricing_info (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
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
	    mpf_EGlpNumFreeArray (pinf->d_scaleinf);
	    pinf->d_scaleinf = mpf_EGlpNumAllocArray (lp->nnbasic);
	} else if (p_price == QS_PRICE_PMULTPARTIAL)
	    pinf->p_strategy = MULTI_PART_PRICING;

	switch (p_price) {
	case QS_PRICE_PDEVEX:
	    if (pinf->pdinfo.norms)
		return rval;
	    rval = mpf_ILLprice_build_pdevex_norms (lp, &(pinf->pdinfo), 0);
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_PSTEEP:
	    if (pinf->psinfo.norms)
		return rval;
	    rval = mpf_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_PMULTPARTIAL:
	    rval = mpf_ILLprice_build_mpartial_info (lp, pinf, COL_PRICING);
	    ILL_CLEANUP_IF (rval);
	    break;
	}
    } else if (d_price != -1) {
	pinf->cur_price = d_price;

	if (d_price == QS_PRICE_DDANTZIG || d_price == QS_PRICE_DSTEEP ||
	    d_price == QS_PRICE_DDEVEX) {
	    pinf->d_strategy = COMPLETE_PRICING;
	    mpf_EGlpNumFreeArray (pinf->p_scaleinf);
	    pinf->p_scaleinf = mpf_EGlpNumAllocArray (lp->nrows);
	} else if (d_price == QS_PRICE_DMULTPARTIAL)
	    pinf->d_strategy = MULTI_PART_PRICING;

	switch (d_price) {
	case QS_PRICE_DSTEEP:
	    if (pinf->dsinfo.norms)
		return rval;
	    rval = mpf_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_DMULTPARTIAL:
	    rval = mpf_ILLprice_build_mpartial_info (lp, pinf, ROW_PRICING);
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_DDEVEX:
	    if (pinf->ddinfo.norms)
		return rval;
	    rval = mpf_ILLprice_build_ddevex_norms (lp, &(pinf->ddinfo), 0);
	    ILL_CLEANUP_IF (rval);
	    break;
	}
    }
CLEANUP:
    if (rval)
	mpf_ILLprice_free_pricing_info (pinf);
    ILL_RETURN (rval, "mpf_ILLprice_build_pricing_info");
}

int mpf_ILLprice_update_pricing_info (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      int const phase,
      mpf_svector * const wz,
      int const eindex,
      int const lindex,
      mpf_t y)
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
	    rval = mpf_ILLprice_update_pdevex_norms (lp, &(pinf->pdinfo), eindex, y);
	    ILL_CLEANUP_IF (rval);
	} else if (p_price == QS_PRICE_PSTEEP)
	    mpf_ILLprice_update_psteep_norms (lp, &(pinf->psinfo), wz, eindex, y);
    } else if (d_price != -1) {
	if (d_price == QS_PRICE_DSTEEP)
	    mpf_ILLprice_update_dsteep_norms (lp, &(pinf->dsinfo), wz, lindex, y);
	else if (d_price == QS_PRICE_DDEVEX) {
	    rval = mpf_ILLprice_update_ddevex_norms (lp, &(pinf->ddinfo), lindex, y);
	    ILL_CLEANUP_IF (rval);
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLprice_update_pricing_info");
}

int mpf_ILLprice_get_price (mpf_price_info * const p,
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

void mpf_ILLprice_free_mpartial_info (mpf_mpart_info * p)
{
    ILL_IFFREE (p->gstart, int);
    ILL_IFFREE (p->gshift, int);
    ILL_IFFREE (p->gsize, int);
    ILL_IFFREE (p->bucket, int);
    mpf_EGlpNumFreeArray (p->infeas);
    ILL_IFFREE (p->perm, int);
}

int mpf_ILLprice_build_mpartial_info (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      int const pricetype)
{
    mpf_mpart_info *const p =
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
    p->infeas = mpf_EGlpNumAllocArray (2 * p->k);
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
	mpf_ILLprice_free_mpartial_info (p);
    ILL_RETURN (rval, "mpf_ILLprice_build_mpartial_info");
}

void mpf_ILLprice_init_mpartial_price (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      int const phase,
      int const pricetype)
{
    mpf_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i;
    p->bsize = 0;
    i = p->cgroup;
    do {
	mpf_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
	i = (i + 1) % p->ngroups;
    } while (i != p->cgroup && p->bsize <= p->k);
    p->cgroup = i;
}

void mpf_ILLprice_update_mpartial_price (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      int const phase,
      int const pricetype)
{
    mpf_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i = 0;
    int csize = 0;
    mpf_t infeas;
    mpf_price_res pr;
    mpf_EGlpNumInitVar (pr.dinfeas);
    mpf_EGlpNumInitVar (pr.pinfeas);
    mpf_EGlpNumInitVar (infeas);

#ifdef mpf_MULTIP
    i = 0;
    while (i < p->bsize) {
	if (pricetype == COL_PRICING) {
	    mpf_ILLprice_column (lp, p->bucket[i], phase, &pr);
	    mpf_EGlpNumCopy (infeas, pr.dinfeas);
	} else {
	    mpf_ILLprice_row (lp, p->bucket[i], phase, &pr);
	    mpf_EGlpNumCopy (infeas, pr.pinfeas);
	}
	if (mpf_EGlpNumIsEqqual (infeas, mpf_zeroLpNum)) {
	    p->bucket[i] = p->bucket[p->bsize - 1];
	    p->bsize--;
	} else {
	    mpf_EGlpNumCopy (p->infeas[i], infeas);
	    i++;
	}
    }
    if (p->bsize > 0) {
	for (i = 0; i < p->bsize; i++)
	    p->perm[i] = i;
	mpf_EGutilPermSort ((size_t) (p->bsize), p->perm, (const mpf_t * const) p->infeas);

	csize = MIN (p->bsize, p->k);
	for (i = csize - 1; i >= 0; i--)
	    lp->iwork[p->bucket[p->perm[i]]] = 1;

	for (i = 0, csize = 0; i < p->bsize; i++)
	    if (lp->iwork[p->bucket[i]] == 1) {
		mpf_EGlpNumCopy (p->infeas[csize], p->infeas[i]);
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
	mpf_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
	i = (i + 1) % p->ngroups;
    } while (i != p->cgroup && p->bsize <= p->k);
    p->cgroup = i;

#ifdef mpf_MULTIP
    for (i = 0; i < csize; i++)
	lp->iwork[p->bucket[i]] = 0;
#endif
    mpf_EGlpNumClearVar (infeas);
    mpf_EGlpNumClearVar (pr.pinfeas);
    mpf_EGlpNumClearVar (pr.dinfeas);
}

void mpf_ILLprice_delete_onempart_price (mpf_price_info * const pinf,
      int const indx, int const pricetype)
{
    mpf_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i = 0;

    for (i = 0; i < p->bsize; i++)
	if (p->bucket[i] == indx) {
	    p->bucket[i] = p->bucket[p->bsize - 1];
	    mpf_EGlpNumCopy (p->infeas[i], p->infeas[p->bsize - 1]);
	    p->bsize--;
	    break;
	}
}

void mpf_ILLprice_mpartial_group (mpf_lpinfo * const lp,
      mpf_mpart_info * const p,
      int const phase,
      int const g,
      int const pricetype)
{
    int i, ix;
    int gstart = p->gstart[g];
    int gsize = p->gsize[g];
    int gshift = p->gshift[g];
    mpf_t infeas;
    mpf_price_res pr;
    mpf_EGlpNumInitVar (pr.dinfeas);
    mpf_EGlpNumInitVar (pr.pinfeas);
    mpf_EGlpNumInitVar (infeas);

    for (i = 0, ix = gstart; i < gsize; i++, ix += gshift) {
#ifdef mpf_MULTIP
	if (lp->iwork[ix])
	    continue;
#endif
	if (pricetype == COL_PRICING) {
	    mpf_ILLprice_column (lp, ix, phase, &pr);
	    mpf_EGlpNumCopy (infeas, pr.dinfeas);
	} else {
	    mpf_ILLprice_row (lp, ix, phase, &pr);
	    mpf_EGlpNumCopy (infeas, pr.pinfeas);
	}
	if (mpf_EGlpNumIsNeqqZero (infeas)) {
	    mpf_EGlpNumCopy (p->infeas[p->bsize], infeas);
	    p->bucket[p->bsize] = ix;
	    p->bsize++;
	}
    }
    mpf_EGlpNumClearVar (infeas);
    mpf_EGlpNumClearVar (pr.dinfeas);
    mpf_EGlpNumClearVar (pr.pinfeas);
}

void mpf_ILLprice_column (mpf_lpinfo * const lp,
      int const ix,
      int const phase,
      mpf_price_res * const pr)
{
    int const col = lp->nbaz[ix];
    int const mcnt = lp->matcnt[col];
    int const mbeg = lp->matbeg[col];
    int i;
    mpf_t sum;
    mpf_EGlpNumZero (pr->dinfeas);
    if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	return;
    mpf_EGlpNumInitVar (sum);
    mpf_EGlpNumZero (sum);

    if (phase == PRIMAL_PHASEII) {
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpf_EGlpNumCopyDiff (lp->dz[ix], lp->cz[col], sum);
	mpf_compute_dualII_inf (lp, ix, &(pr->dinfeas));
    } else {
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpf_EGlpNumCopyNeg (lp->pIdz[ix], sum);
	mpf_compute_dualI_inf (lp, ix, &(pr->dinfeas));
    }
    mpf_EGlpNumClearVar (sum);
}

void mpf_ILLprice_row (mpf_lpinfo * const lp,
      int const ix,
      int const phase,
      mpf_price_res * const pr)
{
    if (phase == DUAL_PHASEII)
	mpf_compute_primalII_inf (lp, ix, &(pr->pinfeas));
    else
	mpf_compute_primalI_inf (lp, ix, &(pr->pinfeas));
}

int mpf_ILLprice_build_pdevex_norms (mpf_lpinfo * const lp,
      mpf_p_devex_info * const pdinfo,
      int const reinit)
{
    int j;
    int rval = 0;

    if (reinit == 0) {
	pdinfo->ninit = 0;
	pdinfo->norms = mpf_EGlpNumAllocArray (lp->nnbasic);
	ILL_SAFE_MALLOC (pdinfo->refframe, lp->ncols, int);
    }
    if (reinit != 0)
	pdinfo->ninit++;

    for (j = 0; j < lp->ncols; j++) {
	if (lp->vstat[j] == STAT_BASIC)
	    pdinfo->refframe[j] = 0;
	else {
	    mpf_EGlpNumOne (pdinfo->norms[lp->vindex[j]]);
	    pdinfo->refframe[j] = 1;
	}
    }

CLEANUP:
    if (rval) {
	mpf_EGlpNumFreeArray (pdinfo->norms);
	ILL_IFFREE (pdinfo->refframe, int);
    }
    ILL_RETURN (rval, "mpf_ILLprice_build_pdevex_norms");
}

int mpf_ILLprice_update_pdevex_norms (mpf_lpinfo * const lp,
      mpf_p_devex_info * const pdinfo,
      int const eindex,
      mpf_t yl)
{
    int i, j;
    mpf_t normj;
    mpf_t zAj;
    mpf_EGlpNumInitVar (normj);
    mpf_EGlpNumInitVar (zAj);
    mpf_EGlpNumZero (normj);

    for (i = 0; i < lp->yjz.nzcnt; i++)
	if (pdinfo->refframe[lp->baz[lp->yjz.indx[i]]])
	    mpf_EGlpNumAddInnProdTo (normj, lp->yjz.coef[i], lp->yjz.coef[i]);

    if (pdinfo->refframe[lp->nbaz[eindex]])
	mpf_EGlpNumAddTo (normj, mpf_oneLpNum);

    mpf_EGlpNumCopyFrac (zAj, normj, pdinfo->norms[eindex]);
    if (mpf_EGlpNumIsGreaDbl (zAj, 0x1p10) || mpf_EGlpNumIsLessDbl (zAj, 0x1p-10)) {
	mpf_EGlpNumClearVar (zAj);
	mpf_EGlpNumClearVar (normj);
	return mpf_ILLprice_build_pdevex_norms (lp, pdinfo, 1);
    }
    for (i = 0; i < lp->zA.nzcnt; i++) {
	j = lp->zA.indx[i];
	mpf_EGlpNumCopyFrac (zAj, lp->zA.coef[i], yl);
	mpf_EGlpNumMultTo (zAj, zAj);
	mpf_EGlpNumMultTo (zAj, normj);
	if (mpf_EGlpNumIsLess (pdinfo->norms[j], zAj))
	    mpf_EGlpNumCopy (pdinfo->norms[j], zAj);
    }
    mpf_EGlpNumDivTo (normj, yl);
    mpf_EGlpNumDivTo (normj, yl);
    if (mpf_EGlpNumIsLess (normj, mpf_oneLpNum))
	mpf_EGlpNumCopy (pdinfo->norms[eindex], mpf_oneLpNum);
    else
	mpf_EGlpNumCopy (pdinfo->norms[eindex], normj);
    mpf_EGlpNumClearVar (zAj);
    mpf_EGlpNumClearVar (normj);
    return 0;
}

int mpf_ILLprice_build_psteep_norms (mpf_lpinfo * const lp,
      mpf_p_steep_info * const psinfo)
{
    int j;
    int rval = 0;
    mpf_svector yz;

    mpf_ILLsvector_init (&yz);
    rval = mpf_ILLsvector_alloc (&yz, lp->nrows);
    ILL_CLEANUP_IF (rval);
    psinfo->norms = mpf_EGlpNumAllocArray (lp->nnbasic);

    for (j = 0; j < lp->nnbasic; j++) {
	rval = ILLstring_report (NULL, &lp->O->reporter);
	ILL_CLEANUP_IF (rval);
	mpf_ILLfct_compute_yz (lp, &yz, 0, lp->nbaz[j]);
	mpf_EGlpNumInnProd (psinfo->norms[j], yz.coef, yz.coef, (size_t) yz.nzcnt);
	mpf_EGlpNumAddTo (psinfo->norms[j], mpf_oneLpNum);
    }

CLEANUP:
    mpf_ILLsvector_free (&yz);
    if (rval)
	mpf_EGlpNumFreeArray (psinfo->norms);
    ILL_RETURN (rval, "mpf_ILLprice_build_psteep_norms");
}

void mpf_ILLprice_update_psteep_norms (mpf_lpinfo * const lp,
      mpf_p_steep_info * const psinfo,
      mpf_svector * const wz,
      int const eindex,
      mpf_t yl)
{
    int i, j, k;
    int mcnt, mbeg;
    mpf_t normj;
    mpf_t zAj, wAj;
    mpf_t *v = 0;
    mpf_EGlpNumInitVar (normj);
    mpf_EGlpNumInitVar (zAj);
    mpf_EGlpNumInitVar (wAj);
    mpf_EGlpNumInnProd (normj, lp->yjz.coef, lp->yjz.coef, (size_t) (lp->yjz.nzcnt));
    mpf_EGlpNumAddTo (normj, mpf_oneLpNum);

#if 0
    Bico - remove warnings for dist
	if (fabs ((normj - psinfo->norms[eindex]) / normj) > 1000.0 /* 0.01 */ ) {
	    printf ("warning: incorrect norm values\n");
	    printf ("anorm = %.6f, pnorm = %.6f\n", normj, psinfo->norms[eindex]);
	    fflush (stdout);
	}
#endif

    mpf_ILLfct_load_workvector (lp, wz);
    v = lp->work.coef;

    for (k = 0; k < lp->zA.nzcnt; k++) {
	j = lp->zA.indx[k];
	mpf_EGlpNumZero (wAj);
	mcnt = lp->matcnt[lp->nbaz[j]];
	mbeg = lp->matbeg[lp->nbaz[j]];
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (wAj, lp->matval[mbeg + i], v[lp->matind[mbeg + i]]);
	mpf_EGlpNumCopy (zAj, lp->zA.coef[k]);
	/* psinfo->norms[j] += (zAj * ((zAj * normj / yl) - (2.0 * wAj))) /
	   yl; */
	mpf_EGlpNumMultUiTo (wAj, 2U);
	mpf_EGlpNumMultTo (zAj, normj);
	mpf_EGlpNumDivTo (zAj, yl);
	mpf_EGlpNumSubTo (zAj, wAj);
	mpf_EGlpNumMultTo (zAj, lp->zA.coef[k]);
	mpf_EGlpNumDivTo (zAj, yl);
	mpf_EGlpNumAddTo (psinfo->norms[j], zAj);
	if (mpf_EGlpNumIsLess (psinfo->norms[j], mpf_oneLpNum))
	    mpf_EGlpNumOne (psinfo->norms[j]);
    }

    mpf_EGlpNumCopyFrac (psinfo->norms[eindex], normj, yl);
    mpf_EGlpNumDivTo (psinfo->norms[eindex], yl);
    if (mpf_EGlpNumIsLess (psinfo->norms[eindex], mpf_oneLpNum))
	mpf_EGlpNumOne (psinfo->norms[eindex]);

    mpf_ILLfct_zero_workvector (lp);
    mpf_EGlpNumClearVar (wAj);
    mpf_EGlpNumClearVar (zAj);
    mpf_EGlpNumClearVar (normj);
}

int mpf_ILLprice_build_ddevex_norms (mpf_lpinfo * const lp,
      mpf_d_devex_info * const ddinfo,
      int const reinit)
{
    int i;
    int rval = 0;

    if (reinit == 0) {
	ddinfo->ninit = 0;
	ddinfo->norms = mpf_EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (ddinfo->refframe, lp->ncols, int);
    }
    if (reinit != 0)
	ddinfo->ninit++;

    for (i = 0; i < lp->ncols; i++)
	ddinfo->refframe[i] = (lp->vstat[i] == STAT_BASIC) ? 1 : 0;

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumOne (ddinfo->norms[i]);

CLEANUP:
    if (rval) {
	mpf_EGlpNumFreeArray (ddinfo->norms);
	ILL_IFFREE (ddinfo->refframe, int);
    }
    ILL_RETURN (rval, "mpf_ILLprice_build_ddevex_norms");
}

int mpf_ILLprice_update_ddevex_norms (mpf_lpinfo * const lp,
      mpf_d_devex_info * const ddinfo,
      int const lindex,
      mpf_t yl)
{
    int i, r;
    mpf_t normi;
    mpf_t yr;
    mpf_EGlpNumInitVar (normi);
    mpf_EGlpNumInitVar (yr);
    mpf_EGlpNumZero (normi);

    for (i = 0; i < lp->zA.nzcnt; i++)
	if (ddinfo->refframe[lp->nbaz[lp->zA.indx[i]]])
	    mpf_EGlpNumAddInnProdTo (normi, lp->zA.coef[i], lp->zA.coef[i]);

    if (ddinfo->refframe[lp->baz[lindex]])
	mpf_EGlpNumAddTo (normi, mpf_oneLpNum);

    mpf_EGlpNumCopyFrac (yr, normi, ddinfo->norms[lindex]);
    if (mpf_EGlpNumIsGreaDbl (yr, 0x1p10) || mpf_EGlpNumIsLessDbl (yr, 0x1p-10)) {
	mpf_EGlpNumClearVar (normi);
	mpf_EGlpNumClearVar (yr);
	return mpf_ILLprice_build_ddevex_norms (lp, ddinfo, 1);
    }
    mpf_EGlpNumDivTo (normi, yl);
    mpf_EGlpNumDivTo (normi, yl);

    for (i = 0; i < lp->yjz.nzcnt; i++) {
	r = lp->yjz.indx[i];
	mpf_EGlpNumCopy (yr, lp->yjz.coef[i]);
	mpf_EGlpNumMultTo (yr, yr);
	mpf_EGlpNumMultTo (yr, normi);
	if (mpf_EGlpNumIsLess (ddinfo->norms[r], yr))
	    mpf_EGlpNumCopy (ddinfo->norms[r], yr);
    }
    mpf_EGlpNumCopy (ddinfo->norms[lindex], normi);
    if (mpf_EGlpNumIsLess (ddinfo->norms[lindex], mpf_oneLpNum))
	mpf_EGlpNumOne (ddinfo->norms[lindex]);
    mpf_EGlpNumClearVar (normi);
    mpf_EGlpNumClearVar (yr);
    return 0;
}

int mpf_ILLprice_build_dsteep_norms (mpf_lpinfo * const lp,
      mpf_d_steep_info * const dsinfo)
{
    int i;
    int rval = 0;
    mpf_svector z;

    mpf_ILLsvector_init (&z);
    rval = mpf_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);
    dsinfo->norms = mpf_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++) {
	rval = ILLstring_report (NULL, &lp->O->reporter);
	ILL_CLEANUP_IF (rval);

	mpf_ILLfct_compute_zz (lp, &z, i);
	mpf_EGlpNumInnProd (dsinfo->norms[i], z.coef, z.coef, (size_t) z.nzcnt);
	if (mpf_EGlpNumIsLess (dsinfo->norms[i], mpf_PARAM_MIN_DNORM))
	    mpf_EGlpNumCopy (dsinfo->norms[i], mpf_PARAM_MIN_DNORM);
    }

CLEANUP:
    mpf_ILLsvector_free (&z);
    if (rval)
	mpf_EGlpNumFreeArray (dsinfo->norms);
    ILL_RETURN (rval, "mpf_ILLprice_build_dsteep_norms");
}

int mpf_ILLprice_get_dsteep_norms (mpf_lpinfo * const lp,
      int const count,
      int *const rowind,
      mpf_t * const norms)
{
    int i;
    int rval = 0;
    mpf_svector z;

    mpf_ILLsvector_init (&z);
    rval = mpf_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < count; i++) {
	mpf_ILLfct_compute_zz (lp, &z, rowind[i]);
	mpf_EGlpNumInnProd (norms[i], z.coef, z.coef, (size_t) z.nzcnt);
    }

CLEANUP:
    mpf_ILLsvector_free (&z);
    ILL_RETURN (rval, "mpf_ILLprice_get_dsteep_norms");
}

void mpf_ILLprice_update_dsteep_norms (mpf_lpinfo * const lp,
      mpf_d_steep_info * const dsinfo,
      mpf_svector * const wz,
      int const lindex,
      mpf_t yl)
{
    int i, k;
    mpf_t yij;
    mpf_t norml;
    mpf_t *v = 0;
    mpf_EGlpNumInitVar (norml);
    mpf_EGlpNumInitVar (yij);
    mpf_EGlpNumZero (norml);
    mpf_EGlpNumInnProd (norml, lp->zz.coef, lp->zz.coef, (size_t) (lp->zz.nzcnt));

#if 0
    Bico - remove warnings for dist
	if (fabs ((norml - dsinfo->norms[lindex]) / norml) > 1000.0 /* 0.01 */ ) {
	    printf ("warning: incorrect dnorm values\n");
	    printf ("anorm = %.6f, pnorm = %.6f\n", norml, dsinfo->norms[lindex]);
	    fflush (stdout);
	}
#endif

    mpf_ILLfct_load_workvector (lp, wz);
    v = lp->work.coef;

    for (k = 0; k < lp->yjz.nzcnt; k++) {
	i = lp->yjz.indx[k];
	mpf_EGlpNumCopy (yij, lp->yjz.coef[k]);
	mpf_EGlpNumMultTo (yij, norml);
	mpf_EGlpNumDivTo (yij, yl);
	mpf_EGlpNumSubTo (yij, v[i]);
	mpf_EGlpNumSubTo (yij, v[i]);
	mpf_EGlpNumMultTo (yij, lp->yjz.coef[k]);
	mpf_EGlpNumDivTo (yij, yl);
	mpf_EGlpNumAddTo (dsinfo->norms[i], yij);
	if (mpf_EGlpNumIsLess (dsinfo->norms[i], mpf_PARAM_MIN_DNORM))
	    mpf_EGlpNumCopy (dsinfo->norms[i], mpf_PARAM_MIN_DNORM);
    }
    mpf_EGlpNumCopyFrac (dsinfo->norms[lindex], norml, yl);
    mpf_EGlpNumDivTo (dsinfo->norms[lindex], yl);
    if (mpf_EGlpNumIsLess (dsinfo->norms[lindex], mpf_PARAM_MIN_DNORM))
	mpf_EGlpNumCopy (dsinfo->norms[lindex], mpf_PARAM_MIN_DNORM);

    mpf_ILLfct_zero_workvector (lp);
    mpf_EGlpNumClearVar (norml);
    mpf_EGlpNumClearVar (yij);
}

static void mpf_update_d_scaleinf (mpf_price_info * const p,
      mpf_heap * const h,
      int const j,
      mpf_t inf,
      int const prule)
{
    if (mpf_EGlpNumIsEqqual (inf, mpf_zeroLpNum)) {
	mpf_EGlpNumZero (p->d_scaleinf[j]);
	if (h->hexist != 0 && h->loc[j] != -1)
	    mpf_ILLheap_delete (h, j);
    } else {
	if (prule == QS_PRICE_PDANTZIG)
	    mpf_EGlpNumCopy (p->d_scaleinf[j], inf);
	else if (prule == QS_PRICE_PDEVEX)
	    mpf_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->pdinfo.norms[j]);
	else if (prule == QS_PRICE_PSTEEP)
	    mpf_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->psinfo.norms[j]);

	if (h->hexist != 0) {
	    if (h->loc[j] == -1)
		mpf_ILLheap_insert (h, j);
	    else
		mpf_ILLheap_modify (h, j);
	}
    }
}

static void mpf_compute_dualI_inf (mpf_lpinfo * const lp,
      const int j,
      mpf_t * const inf)
{
    int const col = lp->nbaz[j];
    int const vt = lp->vtype[col];
    int const vs = lp->vstat[col];
    mpf_EGlpNumZero (*inf);
    if (vt != VARTIFICIAL && vt != VFIXED) {
	if ((vs == STAT_LOWER || vs == STAT_ZERO) &&
	    mpf_EGlpNumIsSumLess (lp->pIdz[j], lp->tol->id_tol, mpf_zeroLpNum))
	    mpf_EGlpNumCopyAbs (*inf, lp->pIdz[j]);
	else if ((vs == STAT_UPPER || vs == STAT_ZERO) &&
	    mpf_EGlpNumIsLess (lp->tol->id_tol, lp->pIdz[j]))
	    mpf_EGlpNumCopy (*inf, lp->pIdz[j]);
    }
}

static void mpf_compute_dualII_inf (mpf_lpinfo * const lp,
      int const j,
      mpf_t * const inf)
{
    int const col = lp->nbaz[j];
    int const vt = lp->vtype[col];
    int const vs = lp->vstat[col];
    mpf_EGlpNumZero (*inf);
    if (vt != VARTIFICIAL && vt != VFIXED) {
	if ((vs == STAT_LOWER || vs == STAT_ZERO) &&
	    mpf_EGlpNumIsSumLess (lp->dz[j], lp->tol->dfeas_tol, mpf_zeroLpNum))
	    mpf_EGlpNumCopyAbs (*inf, lp->dz[j]);
	else if ((vs == STAT_UPPER || vs == STAT_ZERO) &&
	    mpf_EGlpNumIsLess (lp->tol->dfeas_tol, lp->dz[j]))
	    mpf_EGlpNumCopy (*inf, lp->dz[j]);
    }
}

void mpf_ILLprice_compute_dual_inf (mpf_lpinfo * const lp,
      mpf_price_info * const p,
      int *const ix,
      int const icnt,
      int const phase)
{
    int const price = (phase == PRIMAL_PHASEI) ? p->pI_price : p->pII_price;
    mpf_heap *const h = &(p->h);
    int i;
    mpf_t inf;
    mpf_EGlpNumInitVar (inf);
    mpf_EGlpNumZero (inf);

    if (phase == PRIMAL_PHASEI) {
	if (ix == NULL)
	    for (i = 0; i < lp->nnbasic; i++) {
		mpf_compute_dualI_inf (lp, i, &(inf));
		mpf_update_d_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpf_compute_dualI_inf (lp, ix[i], &(inf));
		mpf_update_d_scaleinf (p, h, ix[i], inf, price);
	    }
    } else if (phase == PRIMAL_PHASEII) {
	if (ix == NULL)
	    for (i = 0; i < lp->nnbasic; i++) {
		mpf_compute_dualII_inf (lp, i, &inf);
		mpf_update_d_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpf_compute_dualII_inf (lp, ix[i], &inf);
		mpf_update_d_scaleinf (p, h, ix[i], inf, price);
	    }
    }
    mpf_EGlpNumClearVar (inf);
}

void mpf_ILLprice_primal (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      mpf_price_res * const pr,
      int const phase)
{
    mpf_heap *const h = &(pinf->h);
    int j, vs;
    mpf_t d;
    mpf_EGlpNumInitVar (d);
    mpf_EGlpNumZero (d);
    pr->eindex = -1;

#if USEHEAP > 0
    mpf_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
	PRIMAL_SIMPLEX, 1);
#endif

    if (pinf->p_strategy == COMPLETE_PRICING) {
	if (h->hexist) {
	    pr->eindex = mpf_ILLheap_findmin (h);
	    if (pr->eindex != -1)
		mpf_ILLheap_delete (h, pr->eindex);
	} else {
	    for (j = 0; j < lp->nnbasic; j++) {
		if (mpf_EGlpNumIsLess (d, pinf->d_scaleinf[j])) {
		    mpf_EGlpNumCopy (d, pinf->d_scaleinf[j]);
		    pr->eindex = j;
		}
	    }
	}
    } else if (pinf->p_strategy == MULTI_PART_PRICING) {
	for (j = 0; j < pinf->pmpinfo.bsize; j++) {
	    if (mpf_EGlpNumIsLess (d, pinf->pmpinfo.infeas[j])) {
		mpf_EGlpNumCopy (d, pinf->pmpinfo.infeas[j]);
		pr->eindex = pinf->pmpinfo.bucket[j];
	    }
	}
    }
    if (pr->eindex < 0)
	pr->price_stat = PRICE_OPTIMAL;
    else {
	if (phase == PRIMAL_PHASEI)
	    mpf_EGlpNumCopy (d, lp->pIdz[pr->eindex]);
	else
	    mpf_EGlpNumCopy (d, lp->dz[pr->eindex]);
	vs = lp->vstat[lp->nbaz[pr->eindex]];

	pr->price_stat = PRICE_NONOPTIMAL;
	if (vs == STAT_UPPER || (vs == STAT_ZERO &&
		mpf_EGlpNumIsLess (lp->tol->dfeas_tol, d)))
	    pr->dir = VDECREASE;
	else
	    pr->dir = VINCREASE;
    }
    mpf_EGlpNumClearVar (d);
}

static void mpf_update_p_scaleinf (mpf_price_info * const p,
      mpf_heap * const h,
      int const i,
      mpf_t inf,
      int const prule)
{
    if (mpf_EGlpNumIsEqqual (inf, mpf_zeroLpNum)) {
	mpf_EGlpNumZero (p->p_scaleinf[i]);
	if (h->hexist != 0 && h->loc[i] != -1)
	    mpf_ILLheap_delete (h, i);
    } else {
	if (prule == QS_PRICE_DDANTZIG)
	    mpf_EGlpNumCopy (p->p_scaleinf[i], inf);
	else if (prule == QS_PRICE_DSTEEP)
	    mpf_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->dsinfo.norms[i]);
	else if (prule == QS_PRICE_DDEVEX)
	    mpf_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->ddinfo.norms[i]);

	if (h->hexist != 0) {
	    if (h->loc[i] == -1)
		mpf_ILLheap_insert (h, i);
	    else
		mpf_ILLheap_modify (h, i);
	}
    }
}

static void mpf_compute_primalI_inf (mpf_lpinfo * const lp,
      int const i,
      mpf_t * const inf)
{
    int const col = lp->baz[i];
    mpf_EGlpNumZero (*inf);
    if (mpf_EGlpNumIsLess (lp->tol->ip_tol, lp->xbz[i]) &&
	mpf_EGlpNumIsNeqq (lp->uz[col], mpf_INFTY))
	mpf_EGlpNumCopy (*inf, lp->xbz[i]);
    else if (mpf_EGlpNumIsNeqq (lp->lz[col], mpf_NINFTY) &&
	mpf_EGlpNumIsSumLess (lp->xbz[i], lp->tol->ip_tol, mpf_zeroLpNum))
	mpf_EGlpNumCopyAbs (*inf, lp->xbz[i]);
}

static void mpf_compute_primalII_inf (mpf_lpinfo * const lp,
      int const i,
      mpf_t * const inf)
{
    int const col = lp->baz[i];
    mpf_EGlpNumZero (*inf);
    if (mpf_EGlpNumIsNeqq (lp->uz[col], mpf_INFTY) &&
	mpf_EGlpNumIsSumLess (lp->uz[col], lp->tol->pfeas_tol, lp->xbz[i]))
	mpf_EGlpNumCopyDiff (*inf, lp->xbz[i], lp->uz[col]);
    else if (mpf_EGlpNumIsNeqq (lp->lz[col], mpf_NINFTY) &&
	mpf_EGlpNumIsSumLess (lp->xbz[i], lp->tol->pfeas_tol, lp->lz[col]))
	mpf_EGlpNumCopyDiff (*inf, lp->lz[col], lp->xbz[i]);
}

void mpf_ILLprice_compute_primal_inf (mpf_lpinfo * const lp,
      mpf_price_info * const p,
      int *const ix,
      int const icnt,
      int const phase)
{
    mpf_heap *const h = &(p->h);
    int const price = (phase == DUAL_PHASEI) ? p->dI_price : p->dII_price;
    int i;
    mpf_t inf;
    mpf_EGlpNumInitVar (inf);
    mpf_EGlpNumZero (inf);

    if (phase == DUAL_PHASEI) {
	if (ix == NULL)
	    for (i = 0; i < lp->nrows; i++) {
		mpf_compute_primalI_inf (lp, i, &inf);
		mpf_update_p_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpf_compute_primalI_inf (lp, ix[i], &inf);
		mpf_update_p_scaleinf (p, h, ix[i], inf, price);
	    }
    } else if (phase == DUAL_PHASEII) {
	if (ix == NULL)
	    for (i = 0; i < lp->nrows; i++) {
		mpf_compute_primalII_inf (lp, i, &inf);
		mpf_update_p_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpf_compute_primalII_inf (lp, ix[i], &inf);
		mpf_update_p_scaleinf (p, h, ix[i], inf, price);
	    }
    }
    mpf_EGlpNumClearVar (inf);
}

void mpf_ILLprice_dual (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      int const phase,
      mpf_price_res * const pr)
{
    mpf_heap *const h = &(pinf->h);
    int i;
    mpf_t d;
    mpf_EGlpNumInitVar (d);
    mpf_EGlpNumZero (d);
    pr->lindex = -1;

#if USEHEAP > 0
    mpf_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf, DUAL_SIMPLEX,
	1);
#endif

    if (pinf->d_strategy == COMPLETE_PRICING) {
	if (h->hexist) {
	    pr->lindex = mpf_ILLheap_findmin (h);
	    if (pr->lindex != -1)
		mpf_ILLheap_delete (h, pr->lindex);
	} else {
	    for (i = 0; i < lp->nrows; i++) {
		if (mpf_EGlpNumIsLess (d, pinf->p_scaleinf[i])) {
		    mpf_EGlpNumCopy (d, pinf->p_scaleinf[i]);
		    pr->lindex = i;
		}
	    }
	}
    } else if (pinf->d_strategy == MULTI_PART_PRICING) {
	for (i = 0; i < pinf->dmpinfo.bsize; i++) {
	    if (mpf_EGlpNumIsLess (d, pinf->dmpinfo.infeas[i])) {
		mpf_EGlpNumCopy (d, pinf->dmpinfo.infeas[i]);
		pr->lindex = pinf->dmpinfo.bucket[i];
	    }
	}
    }
    if (pr->lindex < 0)
	pr->price_stat = PRICE_OPTIMAL;
    else {
	pr->price_stat = NONOPTIMAL;

	if (mpf_EGlpNumIsNeqq (lp->uz[lp->baz[pr->lindex]], mpf_INFTY)) {
	    if (phase == DUAL_PHASEI)
		mpf_EGlpNumZero (d);
	    else
		mpf_EGlpNumCopy (d, lp->uz[lp->baz[pr->lindex]]);
	    if (mpf_EGlpNumIsSumLess (lp->tol->pfeas_tol, d, lp->xbz[pr->lindex]))
		pr->lvstat = STAT_UPPER;
	    else
		pr->lvstat = STAT_LOWER;
	} else
	    pr->lvstat = STAT_LOWER;
    }
    mpf_EGlpNumClearVar (d);
}

int mpf_ILLprice_get_rownorms (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      mpf_t * const rnorms)
{
    int rval = 0;
    int i;

    if (pinf->dsinfo.norms == NULL) {
	rval = mpf_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
	ILL_CLEANUP_IF (rval);
    }
    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumCopy (rnorms[i], pinf->dsinfo.norms[i]);

CLEANUP:
    if (rval)
	mpf_EGlpNumFreeArray (pinf->dsinfo.norms);
    return rval;
}

int mpf_ILLprice_get_colnorms (mpf_lpinfo * const lp,
      mpf_price_info * const pinf,
      mpf_t * const cnorms)
{
    int rval = 0;
    int i, j;

    if (pinf->psinfo.norms == NULL) {
	rval = mpf_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
	ILL_CLEANUP_IF (rval);
    }
    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumZero (cnorms[lp->baz[i]]);
    for (j = 0; j < lp->nnbasic; j++)
	mpf_EGlpNumCopy (cnorms[lp->nbaz[j]], pinf->psinfo.norms[j]);

CLEANUP:
    if (rval)
	mpf_EGlpNumFreeArray (pinf->psinfo.norms);
    return rval;
}

int mpf_ILLprice_get_newnorms (mpf_lpinfo * const lp,
      int const nelems,
      mpf_t * const norms,
      int *const matcnt,
      int *const matbeg,
      int *const matind,
      mpf_t * const matval,
      int const option)
{
    int i, j;
    int rval = 0;
    mpf_svector a;
    mpf_svector y;

    mpf_ILLsvector_init (&y);
    rval = mpf_ILLsvector_alloc (&y, lp->nrows);
    ILL_CLEANUP_IF (rval);

    for (j = 0; j < nelems; j++) {
	a.nzcnt = matcnt[j];
	a.indx = &(matind[matbeg[j]]);
	a.coef = &(matval[matbeg[j]]);

	if (option == COLUMN_SOLVE)
	    mpf_ILLbasis_column_solve (lp, &a, &y);
	else
	    mpf_ILLbasis_row_solve (lp, &a, &y);

	mpf_EGlpNumOne (norms[j]);
	for (i = 0; i < y.nzcnt; i++)
	    mpf_EGlpNumAddInnProdTo (norms[j], y.coef[i], y.coef[i]);
    }

CLEANUP:
    mpf_ILLsvector_free (&y);
    ILL_RETURN (rval, "mpf_ILLprice_get_newnorms");
}

int mpf_ILLprice_get_new_rownorms (mpf_lpinfo * const lp,
      int const newrows,
      mpf_t * const rnorms,
      int *const rmatcnt,
      int *const rmatbeg,
      int *const rmatind,
      mpf_t * const rmatval)
{
    return mpf_ILLprice_get_newnorms (lp, newrows, rnorms, rmatcnt, rmatbeg, rmatind,
	rmatval, ROW_SOLVE);
}

int mpf_ILLprice_get_new_colnorms (mpf_lpinfo * const lp,
      int const newrows,
      mpf_t * const rnorms,
      int *const matcnt,
      int *const matbeg,
      int *const matind,
      mpf_t * const matval)
{
    return mpf_ILLprice_get_newnorms (lp, newrows, rnorms, matcnt, matbeg, matind,
	matval, COLUMN_SOLVE);
}

int mpf_ILLprice_load_rownorms (mpf_lpinfo * const lp,
      mpf_t * const rnorms,
      mpf_price_info * const pinf)
{
    int i;
    int rval = 0;

    mpf_EGlpNumFreeArray (pinf->dsinfo.norms);
    pinf->dsinfo.norms = mpf_EGlpNumAllocArray (lp->nrows);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumCopy (pinf->dsinfo.norms[i], rnorms[i]);
	if (mpf_EGlpNumIsLess (pinf->dsinfo.norms[i], mpf_PARAM_MIN_DNORM))
	    mpf_EGlpNumCopy (pinf->dsinfo.norms[i], mpf_PARAM_MIN_DNORM);
    }

    ILL_RETURN (rval, "mpf_ILLprice_load_rownorms");
}

int mpf_ILLprice_load_colnorms (mpf_lpinfo * const lp,
      mpf_t * const cnorms,
      mpf_price_info * const pinf)
{
    int j;
    int rval = 0;

    mpf_EGlpNumFreeArray (pinf->psinfo.norms);
    pinf->psinfo.norms = mpf_EGlpNumAllocArray (lp->nnbasic);
    for (j = 0; j < lp->nnbasic; j++) {
	mpf_EGlpNumCopy (pinf->psinfo.norms[j], cnorms[lp->nbaz[j]]);
	if (mpf_EGlpNumIsLess (pinf->psinfo.norms[j], mpf_oneLpNum))
	    mpf_EGlpNumOne (pinf->psinfo.norms[j]);
    }

    ILL_RETURN (rval, "mpf_ILLprice_load_colnorms");
}

#if mpf_PRICE_DEBUG > 0
void mpf_test_dsteep_norms (mpf_lpinfo * lp,
      mpf_price_info * p)
{
    int i, errn = 0;
    mpf_t err, diff;
    mpf_t *pn;
    mpf_EGlpNumInitVar (err);
    mpf_EGlpNumInitVar (diff);
    pn = mpf_EGlpNumAllocArray (lp->nrows);
    mpf_EGlpNumZero (err);

    mpf_ILLprice_get_dsteep_norms (lp, lp->yjz.nzcnt, lp->yjz.indx, pn);
    for (i = 0; i < lp->yjz.nzcnt; i++) {
	mpf_EGlpNumCopyDiff (diff, pn[i], p->dsinfo.norms[lp->yjz.indx[i]]);
	if (mpf_EGlpNumIsLess (diff, mpf_zeroLpNum))
	    mpf_EGlpNumSign (diff);
	if (mpf_EGlpNumIsLess (mpf_PFEAS_TOLER, diff)) {
	    errn++;
	    mpf_EGlpNumAddTo (err, diff);
	    mpf_EGlpNumCopy (p->dsinfo.norms[lp->yjz.indx[i]], pn[i]);
	}
    }
    if (errn)
	printf ("%d: dnorm errn = %d, err = %.6f\n", lp->cnts->tot_iter, errn,
	    mpf_EGlpNumToLf (err));
    mpf_EGlpNumFreeArray (pn);
    mpf_EGlpNumClearVar (diff);
    mpf_EGlpNumClearVar (err);
}
#endif
