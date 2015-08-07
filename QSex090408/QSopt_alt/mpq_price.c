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

/* RCS_INFO = "$RCSfile: mpq_price.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;
#include "econfig.h"
#include "stddefs.h"
#include "mpq_qsopt.h"
#include "mpq_lpdefs.h"
#include "mpq_fct.h"
#include "mpq_price.h"
#include "mpq_basis.h"
#include "mpq_iqsutil.h"
#include "mpq_dstruct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#define  mpq_MULTIP 1
#define  mpq_PRICE_DEBUG 0

static void mpq_update_d_scaleinf (mpq_price_info * const p,
      mpq_heap * const h,
      int const j,
      mpq_t inf,
      int const prule),
  mpq_update_p_scaleinf (mpq_price_info * const p,
      mpq_heap * const h,
      int const i,
      mpq_t inf,
      int const prule);

static void mpq_compute_dualI_inf (mpq_lpinfo * const lp,
      int const j,
      mpq_t * const inf),
  mpq_compute_dualII_inf (mpq_lpinfo * const lp,
      int const j,
      mpq_t * const inf),
  mpq_compute_primalI_inf (mpq_lpinfo * const lp,
      int const i,
      mpq_t * const inf),
  mpq_compute_primalII_inf (mpq_lpinfo * const lp,
      int const i,
      mpq_t * const inf);

void mpq_ILLprice_free_heap (mpq_price_info * const pinf)
{
    mpq_ILLheap_free (&(pinf->h));
}

int mpq_ILLprice_build_heap (mpq_price_info * const pinf,
      int const nkeys,
      mpq_t * keylist)
{
    mpq_ILLheap_init (&(pinf->h));
    mpq_EGlpNumSet (pinf->htrigger,
	1.0 +
	(double) nkeys / (PARAM_HEAP_RATIO * mpq_ILLutil_our_log2 (nkeys)));
    return mpq_ILLheap_build (&(pinf->h), nkeys, keylist);
}

int mpq_ILLprice_test_for_heap (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      int const nkeys,
      mpq_t * keylist,
      int const algo,
      int const upd)
{
    mpq_heap *const h = &(pinf->h);
    int rval = 0;
    mpq_t ravg;
    if (upd != 0) {
	mpq_EGlpNumInitVar (ravg);
	if (algo == PRIMAL_SIMPLEX)
	    mpq_EGlpNumCopy (ravg, lp->cnts->za_ravg);
	else
	    mpq_EGlpNumCopy (ravg, lp->cnts->y_ravg);
	if (mpq_EGlpNumIsLeq (ravg, pinf->htrigger))
	    pinf->hineff--;
	else {
	    mpq_EGlpNumDivUiTo (ravg, 2U);
	    if (mpq_EGlpNumIsLess (pinf->htrigger, ravg))
		pinf->hineff++;
	}
	mpq_EGlpNumClearVar (ravg);
    }
    if (h->hexist == 0 && pinf->hineff <= 0) {
	rval = mpq_ILLprice_build_heap (pinf, nkeys, keylist);
	ILL_CLEANUP_IF (rval);
    } else if (h->hexist != 0 && pinf->hineff >= PARAM_HEAP_UTRIGGER) {
	mpq_ILLprice_free_heap (pinf);
	/*
         * printf ("freeing mpq_heap ..\n");
         * printf ("iter = %d, ravg = %.2f, trigger = %.2f\n",
         * lp->cnts->tot_iter, ravg, pinf->htrigger);
         */
    }
CLEANUP:
    if (rval)
	mpq_ILLprice_free_heap (pinf);
    return rval;
}

void mpq_ILLprice_init_pricing_info (mpq_price_info * const pinf)
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
    mpq_ILLheap_init (&(pinf->h));
    mpq_EGlpNumZero (pinf->htrigger);
    pinf->hineff = 0;
}

void mpq_ILLprice_free_pricing_info (mpq_price_info * const pinf)
{
    mpq_EGlpNumFreeArray (pinf->p_scaleinf);
    mpq_EGlpNumFreeArray (pinf->d_scaleinf);
    mpq_EGlpNumFreeArray (pinf->pdinfo.norms);
    ILL_IFFREE (pinf->pdinfo.refframe, int);
    mpq_EGlpNumFreeArray (pinf->psinfo.norms);
    mpq_EGlpNumFreeArray (pinf->ddinfo.norms);
    ILL_IFFREE (pinf->ddinfo.refframe, int);
    mpq_EGlpNumFreeArray (pinf->dsinfo.norms);
    mpq_ILLprice_free_mpartial_info (&(pinf->pmpinfo));
    mpq_ILLprice_free_mpartial_info (&(pinf->dmpinfo));
    mpq_ILLprice_free_heap (pinf);
}

int mpq_ILLprice_build_pricing_info (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
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
	    mpq_EGlpNumFreeArray (pinf->d_scaleinf);
	    pinf->d_scaleinf = mpq_EGlpNumAllocArray (lp->nnbasic);
	} else if (p_price == QS_PRICE_PMULTPARTIAL)
	    pinf->p_strategy = MULTI_PART_PRICING;

	switch (p_price) {
	case QS_PRICE_PDEVEX:
	    if (pinf->pdinfo.norms)
		return rval;
	    rval = mpq_ILLprice_build_pdevex_norms (lp, &(pinf->pdinfo), 0);
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_PSTEEP:
	    if (pinf->psinfo.norms)
		return rval;
	    rval = mpq_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_PMULTPARTIAL:
	    rval = mpq_ILLprice_build_mpartial_info (lp, pinf, COL_PRICING);
	    ILL_CLEANUP_IF (rval);
	    break;
	}
    } else if (d_price != -1) {
	pinf->cur_price = d_price;

	if (d_price == QS_PRICE_DDANTZIG || d_price == QS_PRICE_DSTEEP ||
	    d_price == QS_PRICE_DDEVEX) {
	    pinf->d_strategy = COMPLETE_PRICING;
	    mpq_EGlpNumFreeArray (pinf->p_scaleinf);
	    pinf->p_scaleinf = mpq_EGlpNumAllocArray (lp->nrows);
	} else if (d_price == QS_PRICE_DMULTPARTIAL)
	    pinf->d_strategy = MULTI_PART_PRICING;

	switch (d_price) {
	case QS_PRICE_DSTEEP:
	    if (pinf->dsinfo.norms)
		return rval;
	    rval = mpq_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_DMULTPARTIAL:
	    rval = mpq_ILLprice_build_mpartial_info (lp, pinf, ROW_PRICING);
	    ILL_CLEANUP_IF (rval);
	    break;
	case QS_PRICE_DDEVEX:
	    if (pinf->ddinfo.norms)
		return rval;
	    rval = mpq_ILLprice_build_ddevex_norms (lp, &(pinf->ddinfo), 0);
	    ILL_CLEANUP_IF (rval);
	    break;
	}
    }
CLEANUP:
    if (rval)
	mpq_ILLprice_free_pricing_info (pinf);
    ILL_RETURN (rval, "mpq_ILLprice_build_pricing_info");
}

int mpq_ILLprice_update_pricing_info (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      int const phase,
      mpq_svector * const wz,
      int const eindex,
      int const lindex,
      mpq_t y)
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
	    rval = mpq_ILLprice_update_pdevex_norms (lp, &(pinf->pdinfo), eindex, y);
	    ILL_CLEANUP_IF (rval);
	} else if (p_price == QS_PRICE_PSTEEP)
	    mpq_ILLprice_update_psteep_norms (lp, &(pinf->psinfo), wz, eindex, y);
    } else if (d_price != -1) {
	if (d_price == QS_PRICE_DSTEEP)
	    mpq_ILLprice_update_dsteep_norms (lp, &(pinf->dsinfo), wz, lindex, y);
	else if (d_price == QS_PRICE_DDEVEX) {
	    rval = mpq_ILLprice_update_ddevex_norms (lp, &(pinf->ddinfo), lindex, y);
	    ILL_CLEANUP_IF (rval);
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpq_ILLprice_update_pricing_info");
}

int mpq_ILLprice_get_price (mpq_price_info * const p,
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

void mpq_ILLprice_free_mpartial_info (mpq_mpart_info * p)
{
    ILL_IFFREE (p->gstart, int);
    ILL_IFFREE (p->gshift, int);
    ILL_IFFREE (p->gsize, int);
    ILL_IFFREE (p->bucket, int);
    mpq_EGlpNumFreeArray (p->infeas);
    ILL_IFFREE (p->perm, int);
}

int mpq_ILLprice_build_mpartial_info (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      int const pricetype)
{
    mpq_mpart_info *const p =
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
    p->infeas = mpq_EGlpNumAllocArray (2 * p->k);
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
	mpq_ILLprice_free_mpartial_info (p);
    ILL_RETURN (rval, "mpq_ILLprice_build_mpartial_info");
}

void mpq_ILLprice_init_mpartial_price (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      int const phase,
      int const pricetype)
{
    mpq_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i;
    p->bsize = 0;
    i = p->cgroup;
    do {
	mpq_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
	i = (i + 1) % p->ngroups;
    } while (i != p->cgroup && p->bsize <= p->k);
    p->cgroup = i;
}

void mpq_ILLprice_update_mpartial_price (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      int const phase,
      int const pricetype)
{
    mpq_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i = 0;
    int csize = 0;
    mpq_t infeas;
    mpq_price_res pr;
    mpq_EGlpNumInitVar (pr.dinfeas);
    mpq_EGlpNumInitVar (pr.pinfeas);
    mpq_EGlpNumInitVar (infeas);

#ifdef mpq_MULTIP
    i = 0;
    while (i < p->bsize) {
	if (pricetype == COL_PRICING) {
	    mpq_ILLprice_column (lp, p->bucket[i], phase, &pr);
	    mpq_EGlpNumCopy (infeas, pr.dinfeas);
	} else {
	    mpq_ILLprice_row (lp, p->bucket[i], phase, &pr);
	    mpq_EGlpNumCopy (infeas, pr.pinfeas);
	}
	if (mpq_EGlpNumIsEqqual (infeas, mpq_zeroLpNum)) {
	    p->bucket[i] = p->bucket[p->bsize - 1];
	    p->bsize--;
	} else {
	    mpq_EGlpNumCopy (p->infeas[i], infeas);
	    i++;
	}
    }
    if (p->bsize > 0) {
	for (i = 0; i < p->bsize; i++)
	    p->perm[i] = i;
	mpq_EGutilPermSort ((size_t) (p->bsize), p->perm, (const mpq_t * const) p->infeas);

	csize = MIN (p->bsize, p->k);
	for (i = csize - 1; i >= 0; i--)
	    lp->iwork[p->bucket[p->perm[i]]] = 1;

	for (i = 0, csize = 0; i < p->bsize; i++)
	    if (lp->iwork[p->bucket[i]] == 1) {
		mpq_EGlpNumCopy (p->infeas[csize], p->infeas[i]);
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
	mpq_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
	i = (i + 1) % p->ngroups;
    } while (i != p->cgroup && p->bsize <= p->k);
    p->cgroup = i;

#ifdef mpq_MULTIP
    for (i = 0; i < csize; i++)
	lp->iwork[p->bucket[i]] = 0;
#endif
    mpq_EGlpNumClearVar (infeas);
    mpq_EGlpNumClearVar (pr.pinfeas);
    mpq_EGlpNumClearVar (pr.dinfeas);
}

void mpq_ILLprice_delete_onempart_price (mpq_price_info * const pinf,
      int const indx, int const pricetype)
{
    mpq_mpart_info *const p =
    (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
    int i = 0;

    for (i = 0; i < p->bsize; i++)
	if (p->bucket[i] == indx) {
	    p->bucket[i] = p->bucket[p->bsize - 1];
	    mpq_EGlpNumCopy (p->infeas[i], p->infeas[p->bsize - 1]);
	    p->bsize--;
	    break;
	}
}

void mpq_ILLprice_mpartial_group (mpq_lpinfo * const lp,
      mpq_mpart_info * const p,
      int const phase,
      int const g,
      int const pricetype)
{
    int i, ix;
    int gstart = p->gstart[g];
    int gsize = p->gsize[g];
    int gshift = p->gshift[g];
    mpq_t infeas;
    mpq_price_res pr;
    mpq_EGlpNumInitVar (pr.dinfeas);
    mpq_EGlpNumInitVar (pr.pinfeas);
    mpq_EGlpNumInitVar (infeas);

    for (i = 0, ix = gstart; i < gsize; i++, ix += gshift) {
#ifdef mpq_MULTIP
	if (lp->iwork[ix])
	    continue;
#endif
	if (pricetype == COL_PRICING) {
	    mpq_ILLprice_column (lp, ix, phase, &pr);
	    mpq_EGlpNumCopy (infeas, pr.dinfeas);
	} else {
	    mpq_ILLprice_row (lp, ix, phase, &pr);
	    mpq_EGlpNumCopy (infeas, pr.pinfeas);
	}
	if (mpq_EGlpNumIsNeqqZero (infeas)) {
	    mpq_EGlpNumCopy (p->infeas[p->bsize], infeas);
	    p->bucket[p->bsize] = ix;
	    p->bsize++;
	}
    }
    mpq_EGlpNumClearVar (infeas);
    mpq_EGlpNumClearVar (pr.dinfeas);
    mpq_EGlpNumClearVar (pr.pinfeas);
}

void mpq_ILLprice_column (mpq_lpinfo * const lp,
      int const ix,
      int const phase,
      mpq_price_res * const pr)
{
    int const col = lp->nbaz[ix];
    int const mcnt = lp->matcnt[col];
    int const mbeg = lp->matbeg[col];
    int i;
    mpq_t sum;
    mpq_EGlpNumZero (pr->dinfeas);
    if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
	return;
    mpq_EGlpNumInitVar (sum);
    mpq_EGlpNumZero (sum);

    if (phase == PRIMAL_PHASEII) {
	for (i = 0; i < mcnt; i++)
	    mpq_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpq_EGlpNumCopyDiff (lp->dz[ix], lp->cz[col], sum);
	mpq_compute_dualII_inf (lp, ix, &(pr->dinfeas));
    } else {
	for (i = 0; i < mcnt; i++)
	    mpq_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
		lp->matval[mbeg + i]);
	mpq_EGlpNumCopyNeg (lp->pIdz[ix], sum);
	mpq_compute_dualI_inf (lp, ix, &(pr->dinfeas));
    }
    mpq_EGlpNumClearVar (sum);
}

void mpq_ILLprice_row (mpq_lpinfo * const lp,
      int const ix,
      int const phase,
      mpq_price_res * const pr)
{
    if (phase == DUAL_PHASEII)
	mpq_compute_primalII_inf (lp, ix, &(pr->pinfeas));
    else
	mpq_compute_primalI_inf (lp, ix, &(pr->pinfeas));
}

int mpq_ILLprice_build_pdevex_norms (mpq_lpinfo * const lp,
      mpq_p_devex_info * const pdinfo,
      int const reinit)
{
    int j;
    int rval = 0;

    if (reinit == 0) {
	pdinfo->ninit = 0;
	pdinfo->norms = mpq_EGlpNumAllocArray (lp->nnbasic);
	ILL_SAFE_MALLOC (pdinfo->refframe, lp->ncols, int);
    }
    if (reinit != 0)
	pdinfo->ninit++;

    for (j = 0; j < lp->ncols; j++) {
	if (lp->vstat[j] == STAT_BASIC)
	    pdinfo->refframe[j] = 0;
	else {
	    mpq_EGlpNumOne (pdinfo->norms[lp->vindex[j]]);
	    pdinfo->refframe[j] = 1;
	}
    }

CLEANUP:
    if (rval) {
	mpq_EGlpNumFreeArray (pdinfo->norms);
	ILL_IFFREE (pdinfo->refframe, int);
    }
    ILL_RETURN (rval, "mpq_ILLprice_build_pdevex_norms");
}

int mpq_ILLprice_update_pdevex_norms (mpq_lpinfo * const lp,
      mpq_p_devex_info * const pdinfo,
      int const eindex,
      mpq_t yl)
{
    int i, j;
    mpq_t normj;
    mpq_t zAj;
    mpq_EGlpNumInitVar (normj);
    mpq_EGlpNumInitVar (zAj);
    mpq_EGlpNumZero (normj);

    for (i = 0; i < lp->yjz.nzcnt; i++)
	if (pdinfo->refframe[lp->baz[lp->yjz.indx[i]]])
	    mpq_EGlpNumAddInnProdTo (normj, lp->yjz.coef[i], lp->yjz.coef[i]);

    if (pdinfo->refframe[lp->nbaz[eindex]])
	mpq_EGlpNumAddTo (normj, mpq_oneLpNum);

    mpq_EGlpNumCopyFrac (zAj, normj, pdinfo->norms[eindex]);
    if (mpq_EGlpNumIsGreaDbl (zAj, 0x1p10) || mpq_EGlpNumIsLessDbl (zAj, 0x1p-10)) {
	mpq_EGlpNumClearVar (zAj);
	mpq_EGlpNumClearVar (normj);
	return mpq_ILLprice_build_pdevex_norms (lp, pdinfo, 1);
    }
    for (i = 0; i < lp->zA.nzcnt; i++) {
	j = lp->zA.indx[i];
	mpq_EGlpNumCopyFrac (zAj, lp->zA.coef[i], yl);
	mpq_EGlpNumMultTo (zAj, zAj);
	mpq_EGlpNumMultTo (zAj, normj);
	if (mpq_EGlpNumIsLess (pdinfo->norms[j], zAj))
	    mpq_EGlpNumCopy (pdinfo->norms[j], zAj);
    }
    mpq_EGlpNumDivTo (normj, yl);
    mpq_EGlpNumDivTo (normj, yl);
    if (mpq_EGlpNumIsLess (normj, mpq_oneLpNum))
	mpq_EGlpNumCopy (pdinfo->norms[eindex], mpq_oneLpNum);
    else
	mpq_EGlpNumCopy (pdinfo->norms[eindex], normj);
    mpq_EGlpNumClearVar (zAj);
    mpq_EGlpNumClearVar (normj);
    return 0;
}

int mpq_ILLprice_build_psteep_norms (mpq_lpinfo * const lp,
      mpq_p_steep_info * const psinfo)
{
    int j;
    int rval = 0;
    mpq_svector yz;

    mpq_ILLsvector_init (&yz);
    rval = mpq_ILLsvector_alloc (&yz, lp->nrows);
    ILL_CLEANUP_IF (rval);
    psinfo->norms = mpq_EGlpNumAllocArray (lp->nnbasic);

    for (j = 0; j < lp->nnbasic; j++) {
	rval = ILLstring_report (NULL, &lp->O->reporter);
	ILL_CLEANUP_IF (rval);
	mpq_ILLfct_compute_yz (lp, &yz, 0, lp->nbaz[j]);
	mpq_EGlpNumInnProd (psinfo->norms[j], yz.coef, yz.coef, (size_t) yz.nzcnt);
	mpq_EGlpNumAddTo (psinfo->norms[j], mpq_oneLpNum);
    }

CLEANUP:
    mpq_ILLsvector_free (&yz);
    if (rval)
	mpq_EGlpNumFreeArray (psinfo->norms);
    ILL_RETURN (rval, "mpq_ILLprice_build_psteep_norms");
}

void mpq_ILLprice_update_psteep_norms (mpq_lpinfo * const lp,
      mpq_p_steep_info * const psinfo,
      mpq_svector * const wz,
      int const eindex,
      mpq_t yl)
{
    int i, j, k;
    int mcnt, mbeg;
    mpq_t normj;
    mpq_t zAj, wAj;
    mpq_t *v = 0;
    mpq_EGlpNumInitVar (normj);
    mpq_EGlpNumInitVar (zAj);
    mpq_EGlpNumInitVar (wAj);
    mpq_EGlpNumInnProd (normj, lp->yjz.coef, lp->yjz.coef, (size_t) (lp->yjz.nzcnt));
    mpq_EGlpNumAddTo (normj, mpq_oneLpNum);

#if 0
    Bico - remove warnings for dist
	if (fabs ((normj - psinfo->norms[eindex]) / normj) > 1000.0 /* 0.01 */ ) {
	    printf ("warning: incorrect norm values\n");
	    printf ("anorm = %.6f, pnorm = %.6f\n", normj, psinfo->norms[eindex]);
	    fflush (stdout);
	}
#endif

    mpq_ILLfct_load_workvector (lp, wz);
    v = lp->work.coef;

    for (k = 0; k < lp->zA.nzcnt; k++) {
	j = lp->zA.indx[k];
	mpq_EGlpNumZero (wAj);
	mcnt = lp->matcnt[lp->nbaz[j]];
	mbeg = lp->matbeg[lp->nbaz[j]];
	for (i = 0; i < mcnt; i++)
	    mpq_EGlpNumAddInnProdTo (wAj, lp->matval[mbeg + i], v[lp->matind[mbeg + i]]);
	mpq_EGlpNumCopy (zAj, lp->zA.coef[k]);
	/* psinfo->norms[j] += (zAj * ((zAj * normj / yl) - (2.0 * wAj))) /
	   yl; */
	mpq_EGlpNumMultUiTo (wAj, 2U);
	mpq_EGlpNumMultTo (zAj, normj);
	mpq_EGlpNumDivTo (zAj, yl);
	mpq_EGlpNumSubTo (zAj, wAj);
	mpq_EGlpNumMultTo (zAj, lp->zA.coef[k]);
	mpq_EGlpNumDivTo (zAj, yl);
	mpq_EGlpNumAddTo (psinfo->norms[j], zAj);
	if (mpq_EGlpNumIsLess (psinfo->norms[j], mpq_oneLpNum))
	    mpq_EGlpNumOne (psinfo->norms[j]);
    }

    mpq_EGlpNumCopyFrac (psinfo->norms[eindex], normj, yl);
    mpq_EGlpNumDivTo (psinfo->norms[eindex], yl);
    if (mpq_EGlpNumIsLess (psinfo->norms[eindex], mpq_oneLpNum))
	mpq_EGlpNumOne (psinfo->norms[eindex]);

    mpq_ILLfct_zero_workvector (lp);
    mpq_EGlpNumClearVar (wAj);
    mpq_EGlpNumClearVar (zAj);
    mpq_EGlpNumClearVar (normj);
}

int mpq_ILLprice_build_ddevex_norms (mpq_lpinfo * const lp,
      mpq_d_devex_info * const ddinfo,
      int const reinit)
{
    int i;
    int rval = 0;

    if (reinit == 0) {
	ddinfo->ninit = 0;
	ddinfo->norms = mpq_EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (ddinfo->refframe, lp->ncols, int);
    }
    if (reinit != 0)
	ddinfo->ninit++;

    for (i = 0; i < lp->ncols; i++)
	ddinfo->refframe[i] = (lp->vstat[i] == STAT_BASIC) ? 1 : 0;

    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumOne (ddinfo->norms[i]);

CLEANUP:
    if (rval) {
	mpq_EGlpNumFreeArray (ddinfo->norms);
	ILL_IFFREE (ddinfo->refframe, int);
    }
    ILL_RETURN (rval, "mpq_ILLprice_build_ddevex_norms");
}

int mpq_ILLprice_update_ddevex_norms (mpq_lpinfo * const lp,
      mpq_d_devex_info * const ddinfo,
      int const lindex,
      mpq_t yl)
{
    int i, r;
    mpq_t normi;
    mpq_t yr;
    mpq_EGlpNumInitVar (normi);
    mpq_EGlpNumInitVar (yr);
    mpq_EGlpNumZero (normi);

    for (i = 0; i < lp->zA.nzcnt; i++)
	if (ddinfo->refframe[lp->nbaz[lp->zA.indx[i]]])
	    mpq_EGlpNumAddInnProdTo (normi, lp->zA.coef[i], lp->zA.coef[i]);

    if (ddinfo->refframe[lp->baz[lindex]])
	mpq_EGlpNumAddTo (normi, mpq_oneLpNum);

    mpq_EGlpNumCopyFrac (yr, normi, ddinfo->norms[lindex]);
    if (mpq_EGlpNumIsGreaDbl (yr, 0x1p10) || mpq_EGlpNumIsLessDbl (yr, 0x1p-10)) {
	mpq_EGlpNumClearVar (normi);
	mpq_EGlpNumClearVar (yr);
	return mpq_ILLprice_build_ddevex_norms (lp, ddinfo, 1);
    }
    mpq_EGlpNumDivTo (normi, yl);
    mpq_EGlpNumDivTo (normi, yl);

    for (i = 0; i < lp->yjz.nzcnt; i++) {
	r = lp->yjz.indx[i];
	mpq_EGlpNumCopy (yr, lp->yjz.coef[i]);
	mpq_EGlpNumMultTo (yr, yr);
	mpq_EGlpNumMultTo (yr, normi);
	if (mpq_EGlpNumIsLess (ddinfo->norms[r], yr))
	    mpq_EGlpNumCopy (ddinfo->norms[r], yr);
    }
    mpq_EGlpNumCopy (ddinfo->norms[lindex], normi);
    if (mpq_EGlpNumIsLess (ddinfo->norms[lindex], mpq_oneLpNum))
	mpq_EGlpNumOne (ddinfo->norms[lindex]);
    mpq_EGlpNumClearVar (normi);
    mpq_EGlpNumClearVar (yr);
    return 0;
}

int mpq_ILLprice_build_dsteep_norms (mpq_lpinfo * const lp,
      mpq_d_steep_info * const dsinfo)
{
    int i;
    int rval = 0;
    mpq_svector z;

    mpq_ILLsvector_init (&z);
    rval = mpq_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);
    dsinfo->norms = mpq_EGlpNumAllocArray (lp->nrows);

    for (i = 0; i < lp->nrows; i++) {
	rval = ILLstring_report (NULL, &lp->O->reporter);
	ILL_CLEANUP_IF (rval);

	mpq_ILLfct_compute_zz (lp, &z, i);
	mpq_EGlpNumInnProd (dsinfo->norms[i], z.coef, z.coef, (size_t) z.nzcnt);
	if (mpq_EGlpNumIsLess (dsinfo->norms[i], mpq_PARAM_MIN_DNORM))
	    mpq_EGlpNumCopy (dsinfo->norms[i], mpq_PARAM_MIN_DNORM);
    }

CLEANUP:
    mpq_ILLsvector_free (&z);
    if (rval)
	mpq_EGlpNumFreeArray (dsinfo->norms);
    ILL_RETURN (rval, "mpq_ILLprice_build_dsteep_norms");
}

int mpq_ILLprice_get_dsteep_norms (mpq_lpinfo * const lp,
      int const count,
      int *const rowind,
      mpq_t * const norms)
{
    int i;
    int rval = 0;
    mpq_svector z;

    mpq_ILLsvector_init (&z);
    rval = mpq_ILLsvector_alloc (&z, lp->nrows);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < count; i++) {
	mpq_ILLfct_compute_zz (lp, &z, rowind[i]);
	mpq_EGlpNumInnProd (norms[i], z.coef, z.coef, (size_t) z.nzcnt);
    }

CLEANUP:
    mpq_ILLsvector_free (&z);
    ILL_RETURN (rval, "mpq_ILLprice_get_dsteep_norms");
}

void mpq_ILLprice_update_dsteep_norms (mpq_lpinfo * const lp,
      mpq_d_steep_info * const dsinfo,
      mpq_svector * const wz,
      int const lindex,
      mpq_t yl)
{
    int i, k;
    mpq_t yij;
    mpq_t norml;
    mpq_t *v = 0;
    mpq_EGlpNumInitVar (norml);
    mpq_EGlpNumInitVar (yij);
    mpq_EGlpNumZero (norml);
    mpq_EGlpNumInnProd (norml, lp->zz.coef, lp->zz.coef, (size_t) (lp->zz.nzcnt));

#if 0
    Bico - remove warnings for dist
	if (fabs ((norml - dsinfo->norms[lindex]) / norml) > 1000.0 /* 0.01 */ ) {
	    printf ("warning: incorrect dnorm values\n");
	    printf ("anorm = %.6f, pnorm = %.6f\n", norml, dsinfo->norms[lindex]);
	    fflush (stdout);
	}
#endif

    mpq_ILLfct_load_workvector (lp, wz);
    v = lp->work.coef;

    for (k = 0; k < lp->yjz.nzcnt; k++) {
	i = lp->yjz.indx[k];
	mpq_EGlpNumCopy (yij, lp->yjz.coef[k]);
	mpq_EGlpNumMultTo (yij, norml);
	mpq_EGlpNumDivTo (yij, yl);
	mpq_EGlpNumSubTo (yij, v[i]);
	mpq_EGlpNumSubTo (yij, v[i]);
	mpq_EGlpNumMultTo (yij, lp->yjz.coef[k]);
	mpq_EGlpNumDivTo (yij, yl);
	mpq_EGlpNumAddTo (dsinfo->norms[i], yij);
	if (mpq_EGlpNumIsLess (dsinfo->norms[i], mpq_PARAM_MIN_DNORM))
	    mpq_EGlpNumCopy (dsinfo->norms[i], mpq_PARAM_MIN_DNORM);
    }
    mpq_EGlpNumCopyFrac (dsinfo->norms[lindex], norml, yl);
    mpq_EGlpNumDivTo (dsinfo->norms[lindex], yl);
    if (mpq_EGlpNumIsLess (dsinfo->norms[lindex], mpq_PARAM_MIN_DNORM))
	mpq_EGlpNumCopy (dsinfo->norms[lindex], mpq_PARAM_MIN_DNORM);

    mpq_ILLfct_zero_workvector (lp);
    mpq_EGlpNumClearVar (norml);
    mpq_EGlpNumClearVar (yij);
}

static void mpq_update_d_scaleinf (mpq_price_info * const p,
      mpq_heap * const h,
      int const j,
      mpq_t inf,
      int const prule)
{
    if (mpq_EGlpNumIsEqqual (inf, mpq_zeroLpNum)) {
	mpq_EGlpNumZero (p->d_scaleinf[j]);
	if (h->hexist != 0 && h->loc[j] != -1)
	    mpq_ILLheap_delete (h, j);
    } else {
	if (prule == QS_PRICE_PDANTZIG)
	    mpq_EGlpNumCopy (p->d_scaleinf[j], inf);
	else if (prule == QS_PRICE_PDEVEX)
	    mpq_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->pdinfo.norms[j]);
	else if (prule == QS_PRICE_PSTEEP)
	    mpq_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->psinfo.norms[j]);

	if (h->hexist != 0) {
	    if (h->loc[j] == -1)
		mpq_ILLheap_insert (h, j);
	    else
		mpq_ILLheap_modify (h, j);
	}
    }
}

static void mpq_compute_dualI_inf (mpq_lpinfo * const lp,
      const int j,
      mpq_t * const inf)
{
    int const col = lp->nbaz[j];
    int const vt = lp->vtype[col];
    int const vs = lp->vstat[col];
    mpq_EGlpNumZero (*inf);
    if (vt != VARTIFICIAL && vt != VFIXED) {
	if ((vs == STAT_LOWER || vs == STAT_ZERO) &&
	    mpq_EGlpNumIsSumLess (lp->pIdz[j], lp->tol->id_tol, mpq_zeroLpNum))
	    mpq_EGlpNumCopyAbs (*inf, lp->pIdz[j]);
	else if ((vs == STAT_UPPER || vs == STAT_ZERO) &&
	    mpq_EGlpNumIsLess (lp->tol->id_tol, lp->pIdz[j]))
	    mpq_EGlpNumCopy (*inf, lp->pIdz[j]);
    }
}

static void mpq_compute_dualII_inf (mpq_lpinfo * const lp,
      int const j,
      mpq_t * const inf)
{
    int const col = lp->nbaz[j];
    int const vt = lp->vtype[col];
    int const vs = lp->vstat[col];
    mpq_EGlpNumZero (*inf);
    if (vt != VARTIFICIAL && vt != VFIXED) {
	if ((vs == STAT_LOWER || vs == STAT_ZERO) &&
	    mpq_EGlpNumIsSumLess (lp->dz[j], lp->tol->dfeas_tol, mpq_zeroLpNum))
	    mpq_EGlpNumCopyAbs (*inf, lp->dz[j]);
	else if ((vs == STAT_UPPER || vs == STAT_ZERO) &&
	    mpq_EGlpNumIsLess (lp->tol->dfeas_tol, lp->dz[j]))
	    mpq_EGlpNumCopy (*inf, lp->dz[j]);
    }
}

void mpq_ILLprice_compute_dual_inf (mpq_lpinfo * const lp,
      mpq_price_info * const p,
      int *const ix,
      int const icnt,
      int const phase)
{
    int const price = (phase == PRIMAL_PHASEI) ? p->pI_price : p->pII_price;
    mpq_heap *const h = &(p->h);
    int i;
    mpq_t inf;
    mpq_EGlpNumInitVar (inf);
    mpq_EGlpNumZero (inf);

    if (phase == PRIMAL_PHASEI) {
	if (ix == NULL)
	    for (i = 0; i < lp->nnbasic; i++) {
		mpq_compute_dualI_inf (lp, i, &(inf));
		mpq_update_d_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpq_compute_dualI_inf (lp, ix[i], &(inf));
		mpq_update_d_scaleinf (p, h, ix[i], inf, price);
	    }
    } else if (phase == PRIMAL_PHASEII) {
	if (ix == NULL)
	    for (i = 0; i < lp->nnbasic; i++) {
		mpq_compute_dualII_inf (lp, i, &inf);
		mpq_update_d_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpq_compute_dualII_inf (lp, ix[i], &inf);
		mpq_update_d_scaleinf (p, h, ix[i], inf, price);
	    }
    }
    mpq_EGlpNumClearVar (inf);
}

void mpq_ILLprice_primal (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      mpq_price_res * const pr,
      int const phase)
{
    mpq_heap *const h = &(pinf->h);
    int j, vs;
    mpq_t d;
    mpq_EGlpNumInitVar (d);
    mpq_EGlpNumZero (d);
    pr->eindex = -1;

#if USEHEAP > 0
    mpq_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
	PRIMAL_SIMPLEX, 1);
#endif

    if (pinf->p_strategy == COMPLETE_PRICING) {
	if (h->hexist) {
	    pr->eindex = mpq_ILLheap_findmin (h);
	    if (pr->eindex != -1)
		mpq_ILLheap_delete (h, pr->eindex);
	} else {
	    for (j = 0; j < lp->nnbasic; j++) {
		if (mpq_EGlpNumIsLess (d, pinf->d_scaleinf[j])) {
		    mpq_EGlpNumCopy (d, pinf->d_scaleinf[j]);
		    pr->eindex = j;
		}
	    }
	}
    } else if (pinf->p_strategy == MULTI_PART_PRICING) {
	for (j = 0; j < pinf->pmpinfo.bsize; j++) {
	    if (mpq_EGlpNumIsLess (d, pinf->pmpinfo.infeas[j])) {
		mpq_EGlpNumCopy (d, pinf->pmpinfo.infeas[j]);
		pr->eindex = pinf->pmpinfo.bucket[j];
	    }
	}
    }
    if (pr->eindex < 0)
	pr->price_stat = PRICE_OPTIMAL;
    else {
	if (phase == PRIMAL_PHASEI)
	    mpq_EGlpNumCopy (d, lp->pIdz[pr->eindex]);
	else
	    mpq_EGlpNumCopy (d, lp->dz[pr->eindex]);
	vs = lp->vstat[lp->nbaz[pr->eindex]];

	pr->price_stat = PRICE_NONOPTIMAL;
	if (vs == STAT_UPPER || (vs == STAT_ZERO &&
		mpq_EGlpNumIsLess (lp->tol->dfeas_tol, d)))
	    pr->dir = VDECREASE;
	else
	    pr->dir = VINCREASE;
    }
    mpq_EGlpNumClearVar (d);
}

static void mpq_update_p_scaleinf (mpq_price_info * const p,
      mpq_heap * const h,
      int const i,
      mpq_t inf,
      int const prule)
{
    if (mpq_EGlpNumIsEqqual (inf, mpq_zeroLpNum)) {
	mpq_EGlpNumZero (p->p_scaleinf[i]);
	if (h->hexist != 0 && h->loc[i] != -1)
	    mpq_ILLheap_delete (h, i);
    } else {
	if (prule == QS_PRICE_DDANTZIG)
	    mpq_EGlpNumCopy (p->p_scaleinf[i], inf);
	else if (prule == QS_PRICE_DSTEEP)
	    mpq_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->dsinfo.norms[i]);
	else if (prule == QS_PRICE_DDEVEX)
	    mpq_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->ddinfo.norms[i]);

	if (h->hexist != 0) {
	    if (h->loc[i] == -1)
		mpq_ILLheap_insert (h, i);
	    else
		mpq_ILLheap_modify (h, i);
	}
    }
}

static void mpq_compute_primalI_inf (mpq_lpinfo * const lp,
      int const i,
      mpq_t * const inf)
{
    int const col = lp->baz[i];
    mpq_EGlpNumZero (*inf);
    if (mpq_EGlpNumIsLess (lp->tol->ip_tol, lp->xbz[i]) &&
	mpq_EGlpNumIsNeqq (lp->uz[col], mpq_INFTY))
	mpq_EGlpNumCopy (*inf, lp->xbz[i]);
    else if (mpq_EGlpNumIsNeqq (lp->lz[col], mpq_NINFTY) &&
	mpq_EGlpNumIsSumLess (lp->xbz[i], lp->tol->ip_tol, mpq_zeroLpNum))
	mpq_EGlpNumCopyAbs (*inf, lp->xbz[i]);
}

static void mpq_compute_primalII_inf (mpq_lpinfo * const lp,
      int const i,
      mpq_t * const inf)
{
    int const col = lp->baz[i];
    mpq_EGlpNumZero (*inf);
    if (mpq_EGlpNumIsNeqq (lp->uz[col], mpq_INFTY) &&
	mpq_EGlpNumIsSumLess (lp->uz[col], lp->tol->pfeas_tol, lp->xbz[i]))
	mpq_EGlpNumCopyDiff (*inf, lp->xbz[i], lp->uz[col]);
    else if (mpq_EGlpNumIsNeqq (lp->lz[col], mpq_NINFTY) &&
	mpq_EGlpNumIsSumLess (lp->xbz[i], lp->tol->pfeas_tol, lp->lz[col]))
	mpq_EGlpNumCopyDiff (*inf, lp->lz[col], lp->xbz[i]);
}

void mpq_ILLprice_compute_primal_inf (mpq_lpinfo * const lp,
      mpq_price_info * const p,
      int *const ix,
      int const icnt,
      int const phase)
{
    mpq_heap *const h = &(p->h);
    int const price = (phase == DUAL_PHASEI) ? p->dI_price : p->dII_price;
    int i;
    mpq_t inf;
    mpq_EGlpNumInitVar (inf);
    mpq_EGlpNumZero (inf);

    if (phase == DUAL_PHASEI) {
	if (ix == NULL)
	    for (i = 0; i < lp->nrows; i++) {
		mpq_compute_primalI_inf (lp, i, &inf);
		mpq_update_p_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpq_compute_primalI_inf (lp, ix[i], &inf);
		mpq_update_p_scaleinf (p, h, ix[i], inf, price);
	    }
    } else if (phase == DUAL_PHASEII) {
	if (ix == NULL)
	    for (i = 0; i < lp->nrows; i++) {
		mpq_compute_primalII_inf (lp, i, &inf);
		mpq_update_p_scaleinf (p, h, i, inf, price);
	    }
	else
	    for (i = 0; i < icnt; i++) {
		mpq_compute_primalII_inf (lp, ix[i], &inf);
		mpq_update_p_scaleinf (p, h, ix[i], inf, price);
	    }
    }
    mpq_EGlpNumClearVar (inf);
}

void mpq_ILLprice_dual (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      int const phase,
      mpq_price_res * const pr)
{
    mpq_heap *const h = &(pinf->h);
    int i;
    mpq_t d;
    mpq_EGlpNumInitVar (d);
    mpq_EGlpNumZero (d);
    pr->lindex = -1;

#if USEHEAP > 0
    mpq_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf, DUAL_SIMPLEX,
	1);
#endif

    if (pinf->d_strategy == COMPLETE_PRICING) {
	if (h->hexist) {
	    pr->lindex = mpq_ILLheap_findmin (h);
	    if (pr->lindex != -1)
		mpq_ILLheap_delete (h, pr->lindex);
	} else {
	    for (i = 0; i < lp->nrows; i++) {
		if (mpq_EGlpNumIsLess (d, pinf->p_scaleinf[i])) {
		    mpq_EGlpNumCopy (d, pinf->p_scaleinf[i]);
		    pr->lindex = i;
		}
	    }
	}
    } else if (pinf->d_strategy == MULTI_PART_PRICING) {
	for (i = 0; i < pinf->dmpinfo.bsize; i++) {
	    if (mpq_EGlpNumIsLess (d, pinf->dmpinfo.infeas[i])) {
		mpq_EGlpNumCopy (d, pinf->dmpinfo.infeas[i]);
		pr->lindex = pinf->dmpinfo.bucket[i];
	    }
	}
    }
    if (pr->lindex < 0)
	pr->price_stat = PRICE_OPTIMAL;
    else {
	pr->price_stat = NONOPTIMAL;

	if (mpq_EGlpNumIsNeqq (lp->uz[lp->baz[pr->lindex]], mpq_INFTY)) {
	    if (phase == DUAL_PHASEI)
		mpq_EGlpNumZero (d);
	    else
		mpq_EGlpNumCopy (d, lp->uz[lp->baz[pr->lindex]]);
	    if (mpq_EGlpNumIsSumLess (lp->tol->pfeas_tol, d, lp->xbz[pr->lindex]))
		pr->lvstat = STAT_UPPER;
	    else
		pr->lvstat = STAT_LOWER;
	} else
	    pr->lvstat = STAT_LOWER;
    }
    mpq_EGlpNumClearVar (d);
}

int mpq_ILLprice_get_rownorms (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      mpq_t * const rnorms)
{
    int rval = 0;
    int i;

    if (pinf->dsinfo.norms == NULL) {
	rval = mpq_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
	ILL_CLEANUP_IF (rval);
    }
    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumCopy (rnorms[i], pinf->dsinfo.norms[i]);

CLEANUP:
    if (rval)
	mpq_EGlpNumFreeArray (pinf->dsinfo.norms);
    return rval;
}

int mpq_ILLprice_get_colnorms (mpq_lpinfo * const lp,
      mpq_price_info * const pinf,
      mpq_t * const cnorms)
{
    int rval = 0;
    int i, j;

    if (pinf->psinfo.norms == NULL) {
	rval = mpq_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
	ILL_CLEANUP_IF (rval);
    }
    for (i = 0; i < lp->nrows; i++)
	mpq_EGlpNumZero (cnorms[lp->baz[i]]);
    for (j = 0; j < lp->nnbasic; j++)
	mpq_EGlpNumCopy (cnorms[lp->nbaz[j]], pinf->psinfo.norms[j]);

CLEANUP:
    if (rval)
	mpq_EGlpNumFreeArray (pinf->psinfo.norms);
    return rval;
}

int mpq_ILLprice_get_newnorms (mpq_lpinfo * const lp,
      int const nelems,
      mpq_t * const norms,
      int *const matcnt,
      int *const matbeg,
      int *const matind,
      mpq_t * const matval,
      int const option)
{
    int i, j;
    int rval = 0;
    mpq_svector a;
    mpq_svector y;

    mpq_ILLsvector_init (&y);
    rval = mpq_ILLsvector_alloc (&y, lp->nrows);
    ILL_CLEANUP_IF (rval);

    for (j = 0; j < nelems; j++) {
	a.nzcnt = matcnt[j];
	a.indx = &(matind[matbeg[j]]);
	a.coef = &(matval[matbeg[j]]);

	if (option == COLUMN_SOLVE)
	    mpq_ILLbasis_column_solve (lp, &a, &y);
	else
	    mpq_ILLbasis_row_solve (lp, &a, &y);

	mpq_EGlpNumOne (norms[j]);
	for (i = 0; i < y.nzcnt; i++)
	    mpq_EGlpNumAddInnProdTo (norms[j], y.coef[i], y.coef[i]);
    }

CLEANUP:
    mpq_ILLsvector_free (&y);
    ILL_RETURN (rval, "mpq_ILLprice_get_newnorms");
}

int mpq_ILLprice_get_new_rownorms (mpq_lpinfo * const lp,
      int const newrows,
      mpq_t * const rnorms,
      int *const rmatcnt,
      int *const rmatbeg,
      int *const rmatind,
      mpq_t * const rmatval)
{
    return mpq_ILLprice_get_newnorms (lp, newrows, rnorms, rmatcnt, rmatbeg, rmatind,
	rmatval, ROW_SOLVE);
}

int mpq_ILLprice_get_new_colnorms (mpq_lpinfo * const lp,
      int const newrows,
      mpq_t * const rnorms,
      int *const matcnt,
      int *const matbeg,
      int *const matind,
      mpq_t * const matval)
{
    return mpq_ILLprice_get_newnorms (lp, newrows, rnorms, matcnt, matbeg, matind,
	matval, COLUMN_SOLVE);
}

int mpq_ILLprice_load_rownorms (mpq_lpinfo * const lp,
      mpq_t * const rnorms,
      mpq_price_info * const pinf)
{
    int i;
    int rval = 0;

    mpq_EGlpNumFreeArray (pinf->dsinfo.norms);
    pinf->dsinfo.norms = mpq_EGlpNumAllocArray (lp->nrows);
    for (i = 0; i < lp->nrows; i++) {
	mpq_EGlpNumCopy (pinf->dsinfo.norms[i], rnorms[i]);
	if (mpq_EGlpNumIsLess (pinf->dsinfo.norms[i], mpq_PARAM_MIN_DNORM))
	    mpq_EGlpNumCopy (pinf->dsinfo.norms[i], mpq_PARAM_MIN_DNORM);
    }

    ILL_RETURN (rval, "mpq_ILLprice_load_rownorms");
}

int mpq_ILLprice_load_colnorms (mpq_lpinfo * const lp,
      mpq_t * const cnorms,
      mpq_price_info * const pinf)
{
    int j;
    int rval = 0;

    mpq_EGlpNumFreeArray (pinf->psinfo.norms);
    pinf->psinfo.norms = mpq_EGlpNumAllocArray (lp->nnbasic);
    for (j = 0; j < lp->nnbasic; j++) {
	mpq_EGlpNumCopy (pinf->psinfo.norms[j], cnorms[lp->nbaz[j]]);
	if (mpq_EGlpNumIsLess (pinf->psinfo.norms[j], mpq_oneLpNum))
	    mpq_EGlpNumOne (pinf->psinfo.norms[j]);
    }

    ILL_RETURN (rval, "mpq_ILLprice_load_colnorms");
}

#if mpq_PRICE_DEBUG > 0
void mpq_test_dsteep_norms (mpq_lpinfo * lp,
      mpq_price_info * p)
{
    int i, errn = 0;
    mpq_t err, diff;
    mpq_t *pn;
    mpq_EGlpNumInitVar (err);
    mpq_EGlpNumInitVar (diff);
    pn = mpq_EGlpNumAllocArray (lp->nrows);
    mpq_EGlpNumZero (err);

    mpq_ILLprice_get_dsteep_norms (lp, lp->yjz.nzcnt, lp->yjz.indx, pn);
    for (i = 0; i < lp->yjz.nzcnt; i++) {
	mpq_EGlpNumCopyDiff (diff, pn[i], p->dsinfo.norms[lp->yjz.indx[i]]);
	if (mpq_EGlpNumIsLess (diff, mpq_zeroLpNum))
	    mpq_EGlpNumSign (diff);
	if (mpq_EGlpNumIsLess (mpq_PFEAS_TOLER, diff)) {
	    errn++;
	    mpq_EGlpNumAddTo (err, diff);
	    mpq_EGlpNumCopy (p->dsinfo.norms[lp->yjz.indx[i]], pn[i]);
	}
    }
    if (errn)
	printf ("%d: dnorm errn = %d, err = %.6f\n", lp->cnts->tot_iter, errn,
	    mpq_EGlpNumToLf (err));
    mpq_EGlpNumFreeArray (pn);
    mpq_EGlpNumClearVar (diff);
    mpq_EGlpNumClearVar (err);
}
#endif
