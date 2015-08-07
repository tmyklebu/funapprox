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

/* RCS_INFO = "$RCSfile: mpf_simplex.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

#include "basicdefs.h"
#include "config.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdata.h"
#include "mpf_lpdefs.h"

#include "stddefs.h"
#include "mpf_fct.h"
#include "mpf_ratio.h"
#include "mpf_price.h"
#include "mpf_basis.h"
#include "mpf_simplex.h"
#include "mpf_dstruct.h"
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#include "mpf_lib.h"		/* for mpf_ILLlib_writebasis */
#include "mpf_lp.h"		/* for mpf_ILLwrite_lp */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static void mpf_init_lp_status_info (mpf_lp_status_info * ls),
  mpf_init_simplex_tols (mpf_lpinfo * lp),
  mpf_monitor_iter (mpf_lpinfo * lp, mpf_iter_info * it, int cphase),
  mpf_get_current_stat (mpf_lp_status_info * p,
      int algorithm,
      int *bstat);

static int mpf_terminate_simplex (mpf_lpinfo * lp,
      int phase,
      mpf_iter_info * it),
  mpf_primal_phaseI_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it),
  mpf_primal_phaseII_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it),
  mpf_dual_phaseI_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it),
  mpf_dual_phaseII_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it),
  mpf_report_value (mpf_lpinfo * lp,
      mpf_iter_info * it,
      const char *value_name,
      mpf_t value);


void mpf_ILLsimplex_init_lpinfo (mpf_lpinfo * lp)
{
    mpf_ILLbasis_init_basisinfo (lp);
    mpf_init_internal_lpinfo (lp);
}

void mpf_ILLsimplex_free_lpinfo (mpf_lpinfo * lp)
{
    if (lp) {
	mpf_EGlpNumFreeArray (lp->lz);
	mpf_EGlpNumFreeArray (lp->uz);
	mpf_EGlpNumFreeArray (lp->cz);
	mpf_ILLbasis_free_basisinfo (lp);
	mpf_free_internal_lpinfo (lp);
    }
}

void mpf_ILLsimplex_load_lpinfo (mpf_ILLlpdata * qslp,
      mpf_lpinfo * lp)
{
    lp->basisid = -1;
    lp->maxiter = 10000000;
    lp->maxtime = 300000;
    /* lp->iterskip = 10; */
    lp->iterskip = 100;
    mpf_EGlpNumCopy (lp->objbound, mpf_INFTY);
    lp->O = qslp;
}

void mpf_ILLsimplex_set_bound (mpf_lpinfo * lp,
      const mpf_t * objbound,
      int sense)
{
    mpf_EGlpNumCopy (lp->objbound, *objbound);
    if (sense == mpf_ILL_MAX)
	mpf_EGlpNumSign (lp->objbound);
}

static void mpf_init_lp_status_info (mpf_lp_status_info * ls)
{
    ls->optimal = 0;
    ls->primal_feasible = 0;
    ls->primal_infeasible = 0;
    ls->primal_unbounded = 0;
    ls->dual_feasible = 0;
    ls->dual_infeasible = 0;
    ls->dual_unbounded = 0;
}

static void mpf_init_simplex_tols (mpf_lpinfo * lp)
{
#if VERBOSE_LEVEL <= DEBUG
    char *strtmp = 0;
#endif
    mpf_EGlpNumCopy (lp->tol->pfeas_tol, mpf_PFEAS_TOLER);
    mpf_EGlpNumCopy (lp->tol->dfeas_tol, mpf_DFEAS_TOLER);
    mpf_EGlpNumCopy (lp->tol->pivot_tol, mpf_PIVOT_TOLER);
    mpf_EGlpNumCopy (lp->tol->szero_tol, mpf_SZERO_TOLER);
    mpf_EGlpNumCopy (lp->tol->ip_tol, lp->tol->pfeas_tol);
    mpf_EGlpNumCopy (lp->tol->id_tol, lp->tol->dfeas_tol);
    if (mpf_EGlpNumIsNeqqZero (lp->tol->ip_tol)) {
#if VERBOSE_LEVEL <= DEBUG
	strtmp = mpf_EGlpNumGetStr (lp->tol->ip_tol);
	MESSAGE (VERBOSE_LEVEL, "ip_tol %lg %s", mpf_EGlpNumToLf (lp->tol->ip_tol), strtmp);
	EGfree (strtmp);
	strtmp = mpf_EGlpNumGetStr (mpf_epsLpNum);
	MESSAGE (VERBOSE_LEVEL, "eps %lg %s", mpf_EGlpNumToLf (mpf_epsLpNum), strtmp);
	EGfree (strtmp);
	strtmp = mpf_EGlpNumGetStr (mpf_PFEAS_TOLER);
	MESSAGE (VERBOSE_LEVEL, "mpf_PFEAS_TOLER %lg %s", mpf_EGlpNumToLf (mpf_PFEAS_TOLER), strtmp);
	EGfree (strtmp);
#endif
	mpf_EGlpNumDivUiTo (lp->tol->ip_tol, 2UL);
    }
    if (mpf_EGlpNumIsNeqqZero (lp->tol->id_tol)) {
#if VERBOSE_LEVEL <= DEBUG
	strtmp = mpf_EGlpNumGetStr (lp->tol->id_tol);
	MESSAGE (VERBOSE_LEVEL, "id_tol %lg %s", mpf_EGlpNumToLf (lp->tol->id_tol), strtmp);
	EGfree (strtmp);
#endif
	mpf_EGlpNumDivUiTo (lp->tol->id_tol, 2UL);
    }
}

void mpf_init_internal_lpinfo (mpf_lpinfo * lp)
{
    int rval = 0;
    lp->nrows = 0;
    lp->nnbasic = 0;
    lp->localrows = 0;
    lp->rowcnt = 0;
    lp->rowbeg = 0;
    lp->rowind = 0;
    lp->rowval = 0;
    lp->cz = 0;
    lp->lz = 0;
    lp->uz = 0;
    lp->xbz = 0;
    lp->piz = 0;
    lp->dz = 0;
    lp->pIxbz = 0;
    lp->pIpiz = 0;
    lp->pIdz = 0;
    lp->vtype = 0;
    lp->vclass = 0;
    lp->iwork = 0;
    lp->upd.perm = 0;
    lp->upd.ix = 0;
    lp->upd.t = 0;
    lp->bfeas = 0;
    lp->dfeas = 0;
    lp->tol = 0;
    lp->cnts = 0;
    lp->bchanges = 0;
    lp->cchanges = 0;
    mpf_ILLsvector_init (&(lp->zz));
    mpf_ILLsvector_init (&(lp->yjz));
    mpf_ILLsvector_init (&(lp->zA));
    mpf_ILLsvector_init (&(lp->work));
    mpf_ILLsvector_init (&(lp->srhs));
    mpf_ILLsvector_init (&(lp->ssoln));
    ILL_SAFE_MALLOC (lp->tol, 1, mpf_tol_struct);
    mpf_EGlpNumInitVar (lp->tol->pfeas_tol);
    mpf_EGlpNumInitVar (lp->tol->dfeas_tol);
    mpf_EGlpNumInitVar (lp->tol->pivot_tol);
    mpf_EGlpNumInitVar (lp->tol->szero_tol);
    mpf_EGlpNumInitVar (lp->tol->ip_tol);
    mpf_EGlpNumInitVar (lp->tol->id_tol);
    ILL_SAFE_MALLOC (lp->cnts, 1, mpf_count_struct);
    mpf_EGlpNumInitVar (lp->cnts->y_ravg);
    mpf_EGlpNumInitVar (lp->cnts->z_ravg);
    mpf_EGlpNumInitVar (lp->cnts->za_ravg);
CLEANUP:
    if (rval) {
	fprintf (stderr, "\nno memory, in %s, exit\n", __func__);
	exit (1);
    }
}

void mpf_free_internal_lpinfo (mpf_lpinfo * lp)
{
    mpf_bndinfo *binfo = 0;
    mpf_coefinfo *cinfo = 0;

    if (lp->localrows) {
	ILL_IFFREE (lp->rowcnt, int);
	ILL_IFFREE (lp->rowbeg, int);
	ILL_IFFREE (lp->rowind, int);
	mpf_EGlpNumFreeArray (lp->rowval);
	lp->localrows = 0;
    }
    mpf_EGlpNumFreeArray (lp->lz);
    mpf_EGlpNumFreeArray (lp->uz);
    mpf_EGlpNumFreeArray (lp->cz);
    mpf_EGlpNumFreeArray (lp->xbz);
    mpf_EGlpNumFreeArray (lp->piz);
    mpf_EGlpNumFreeArray (lp->pIpiz);
    mpf_EGlpNumFreeArray (lp->dz);
    mpf_EGlpNumFreeArray (lp->pIdz);
    mpf_EGlpNumFreeArray (lp->pIxbz);

    ILL_IFFREE (lp->vtype, int);
    ILL_IFFREE (lp->vclass, char);

    mpf_ILLsvector_free (&(lp->zz));
    mpf_ILLsvector_free (&(lp->yjz));
    mpf_ILLsvector_free (&(lp->zA));
    mpf_ILLsvector_free (&(lp->work));
    mpf_ILLsvector_free (&(lp->srhs));
    mpf_ILLsvector_free (&(lp->ssoln));
    ILL_IFFREE (lp->iwork, int);
    ILL_IFFREE (lp->upd.perm, int);
    ILL_IFFREE (lp->upd.ix, int);
    mpf_EGlpNumFreeArray (lp->upd.t);

    ILL_IFFREE (lp->bfeas, int);
    ILL_IFFREE (lp->dfeas, int);
    if (lp->tol) {
	mpf_EGlpNumClearVar (lp->tol->pfeas_tol);
	mpf_EGlpNumClearVar (lp->tol->dfeas_tol);
	mpf_EGlpNumClearVar (lp->tol->pivot_tol);
	mpf_EGlpNumClearVar (lp->tol->szero_tol);
	mpf_EGlpNumClearVar (lp->tol->ip_tol);
	mpf_EGlpNumClearVar (lp->tol->id_tol);
	ILL_IFFREE (lp->tol, mpf_tol_struct);
    }
    if (lp->cnts) {
	mpf_EGlpNumClearVar (lp->cnts->y_ravg);
	mpf_EGlpNumClearVar (lp->cnts->z_ravg);
	mpf_EGlpNumClearVar (lp->cnts->za_ravg);
	ILL_IFFREE (lp->cnts, mpf_count_struct);
    }
    while (lp->bchanges) {
	binfo = lp->bchanges;
	mpf_EGlpNumClearVar (binfo->pbound);
	mpf_EGlpNumClearVar (binfo->cbound);
	lp->bchanges = binfo->next;
	ILL_IFFREE (binfo, mpf_bndinfo);
    }

    while (lp->cchanges) {
	cinfo = lp->cchanges;
	mpf_EGlpNumClearVar (cinfo->pcoef);
	mpf_EGlpNumClearVar (cinfo->ccoef);
	lp->cchanges = cinfo->next;
	ILL_IFFREE (cinfo, mpf_coefinfo);
    }
}

int mpf_build_internal_lpinfo (mpf_lpinfo * lp)
{
    int rval = 0;
    int i, n;
    mpf_ILLlpdata *qslp = lp->O;
    mpf_ILLlp_sinfo *S = lp->O->sinfo;
    mpf_t *lower, *upper, *obj;
    mpf_ILLlp_rows lprows;
    mpf_ILLmatrix *A;

    mpf_init_lp_status_info (&(lp->probstat));
    mpf_init_lp_status_info (&(lp->basisstat));

    if (S != 0) {
	lp->nrows = S->nrows;
	lp->ncols = S->ncols;
	lp->bz = S->rhs;
	lower = S->lower;
	upper = S->upper;
	obj = S->obj;
	A = &(S->A);
    } else {
	lp->nrows = qslp->nrows;
	lp->ncols = qslp->ncols;
	lp->bz = qslp->rhs;
	lower = qslp->lower;
	upper = qslp->upper;
	obj = qslp->obj;
	A = &(qslp->A);
    }

    lp->matbeg = A->matbeg;
    lp->matcnt = A->matcnt;
    lp->matind = A->matind;
    lp->matval = A->matval;

    lp->nnbasic = lp->ncols - lp->nrows;

    lp->lz = mpf_EGlpNumAllocArray (lp->ncols);
    lp->cz = mpf_EGlpNumAllocArray (lp->ncols);
    lp->uz = mpf_EGlpNumAllocArray (lp->ncols);
    if (!lp->lz || !lp->uz || !lp->cz) {
	fprintf (stderr, "mpf_build_internal_lpinfo\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < lp->ncols; i++) {
	mpf_EGlpNumCopy (lp->lz[i], lower[i]);
	mpf_EGlpNumCopy (lp->uz[i], upper[i]);
	mpf_EGlpNumCopy (lp->cz[i], obj[i]);
	if (qslp->objsense == mpf_ILL_MAX) {
	    mpf_EGlpNumSign (lp->cz[i]);
	}
    }

    if (!lp->O->rA) {
	rval = mpf_ILLlp_rows_init (&lprows, lp->O, 1);
	ILL_CLEANUP_IF (rval);
	lp->rowbeg = lprows.rowbeg;
	lp->rowcnt = lprows.rowcnt;
	lp->rowind = lprows.rowind;
	lp->rowval = lprows.rowval;
	lp->localrows = 1;
    } else {
	/* row format exists, just use pointers */
	lp->rowbeg = lp->O->rA->rowbeg;
	lp->rowcnt = lp->O->rA->rowcnt;
	lp->rowind = lp->O->rA->rowind;
	lp->rowval = lp->O->rA->rowval;
	lp->localrows = 0;
    }

    lp->xbz = mpf_EGlpNumAllocArray (lp->nrows);
    lp->piz = mpf_EGlpNumAllocArray (lp->nrows);
    lp->dz = mpf_EGlpNumAllocArray (lp->nnbasic);
    lp->final_phase = -1;
    lp->infub_ix = -1;

    ILL_SAFE_MALLOC (lp->vtype, lp->ncols, int);
    ILL_SAFE_MALLOC (lp->vclass, lp->ncols, char);

    rval = mpf_ILLsvector_alloc (&(lp->zz), lp->nrows);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLsvector_alloc (&(lp->yjz), lp->nrows);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLsvector_alloc (&(lp->zA), lp->nnbasic);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLsvector_alloc (&(lp->work), lp->ncols);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLsvector_alloc (&(lp->srhs), lp->nrows);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLsvector_alloc (&(lp->ssoln), lp->nrows);
    ILL_CLEANUP_IF (rval);
    ILL_SAFE_MALLOC (lp->iwork, lp->ncols, int);
    for (i = 0; i < lp->ncols; i++) {
	lp->work.indx[i] = 0;
	mpf_EGlpNumZero (lp->work.coef[i]);
	lp->iwork[i] = 0;
    }
    n = lp->nrows > lp->ncols ? 2 * (lp->nrows) + 1 : 2 * (lp->ncols) + 1;
    lp->upd.t = mpf_EGlpNumAllocArray (n);
    ILL_SAFE_MALLOC (lp->upd.perm, n, int);
    ILL_SAFE_MALLOC (lp->upd.ix, n, int);


    ILL_SAFE_MALLOC (lp->bfeas, lp->nrows, int);
    ILL_SAFE_MALLOC (lp->dfeas, lp->nnbasic, int);

    mpf_init_simplex_tols (lp);
    mpf_ILLfct_init_counts (lp);

    lp->nbchange = 0;
    lp->ncchange = 0;

    lp->pIratio = RATIOTEST_HARRIS;
    lp->pIIratio = RATIOTEST_HARRIS;
    lp->dIratio = RATIOTEST_HARRIS;
    lp->dIIratio = RATIOTEST_HARRIS;
    lp->starttime = ILLutil_zeit ();
    ILLutil_sprand (1, &(lp->rstate));

CLEANUP:
    if (rval)
	mpf_free_internal_lpinfo (lp);
    ILL_RETURN (rval, "mpf_build_internal_lpinfo");
}

int mpf_ILLsimplex_retest_psolution (mpf_lpinfo * lp,
      mpf_price_info * p,
      int phase,
      mpf_feas_info * fi)
{
    int rval = 0;
    int fbid = lp->fbasisid;
    int bid = lp->basisid;
    mpf_t *ptol = &(lp->tol->pfeas_tol);
    mpf_t *dtol = &(lp->tol->dfeas_tol);
    mpf_t *iptol = &(lp->tol->ip_tol);
    mpf_t *idtol = &(lp->tol->id_tol);

    fi->pstatus = -1;
    fi->dstatus = -1;
    if (fbid < bid - PARAM_PRIMAL_REFACTORGAP) {
	rval = mpf_ILLbasis_refactor (lp);
	ILL_CLEANUP_IF (rval);
    }
    if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP)
	mpf_ILLfct_compute_xbz (lp);

    if (phase == PRIMAL_PHASEII) {
	if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP) {
	    mpf_ILLfct_compute_piz (lp);
	    mpf_ILLfct_compute_dz (lp);
	    if (p != NULL && p->p_strategy == COMPLETE_PRICING)
		mpf_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
	}
	mpf_ILLfct_compute_pobj (lp);
	mpf_ILLfct_check_pfeasible (lp, fi, *ptol);
	mpf_ILLfct_check_dfeasible (lp, fi, *dtol);
    } else if (phase == PRIMAL_PHASEI) {
	mpf_ILLfct_check_pfeasible (lp, fi, *iptol);
	if (fi->pstatus != PRIMAL_FEASIBLE) {
	    if (lp->pIpiz) {
		mpf_ILLfct_compute_phaseI_piz (lp);
		mpf_ILLfct_compute_phaseI_dz (lp);
		mpf_ILLfct_check_pIdfeasible (lp, fi, *idtol);
		if (p != NULL && p->p_strategy == COMPLETE_PRICING)
		    mpf_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
	    }
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLsimplex_retest_psolution");
}

int mpf_ILLsimplex_retest_dsolution (mpf_lpinfo * lp,
      mpf_price_info * p,
      int phase,
      mpf_feas_info * fi)
{
    int rval = 0;
    int fbid = lp->fbasisid;
    int bid = lp->basisid;
    mpf_t *ptol = &(lp->tol->pfeas_tol);
    mpf_t *dtol = &(lp->tol->dfeas_tol);
    mpf_t *iptol = &(lp->tol->ip_tol);
    mpf_t *idtol = &(lp->tol->id_tol);
    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */

    fi->pstatus = -1;
    fi->dstatus = -1;
    if (fbid < bid - PARAM_DUAL_REFACTORGAP) {
	/* ILL_IFTRACE("Refactor: %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	rval = mpf_ILLbasis_refactor (lp);
	ILL_CLEANUP_IF (rval);
    }
    if (fbid < bid - PARAM_DUAL_RESOLVEGAP) {
	/* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	mpf_ILLfct_compute_piz (lp);
	mpf_ILLfct_compute_dz (lp);
    }
    if (phase == DUAL_PHASEII) {
	/* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	if (fbid < bid - PARAM_DUAL_RESOLVEGAP) {
	    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    mpf_ILLfct_compute_xbz (lp);
	    ILL_CLEANUP_IF (rval);
	    if (p != NULL) {
		/* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
		if (p->d_strategy == COMPLETE_PRICING) {
		    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		       */
		    mpf_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
		} else {
		    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		       */
		    mpf_ILLprice_update_mpartial_price (lp, p, DUAL_PHASEII, ROW_PRICING);
		}
	    }
	}
	/* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	mpf_ILLfct_compute_dobj (lp);
	mpf_ILLfct_check_dfeasible (lp, fi, *dtol);
	mpf_ILLfct_check_pfeasible (lp, fi, *ptol);
    } else if (phase == DUAL_PHASEI) {
	mpf_ILLfct_check_dfeasible (lp, fi, *idtol);
	if (fi->dstatus != DUAL_FEASIBLE) {
	    mpf_ILLfct_compute_phaseI_xbz (lp);
	    mpf_ILLfct_check_pIpfeasible (lp, fi, *iptol);
	    if (p != NULL) {
		if (p->d_strategy == COMPLETE_PRICING)
		    mpf_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
		else
		    mpf_ILLprice_update_mpartial_price (lp, p, DUAL_PHASEI, ROW_PRICING);
	    }
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLsimplex_retest_dsolution");
}

int mpf_ILLsimplex_solution (mpf_lpinfo * lp,
      mpf_t * xz,
      mpf_t * piz,
      mpf_t * dz,
      mpf_t * objval)
{
    int i, j;
    int col;

    if (xz != NULL) {
	if (lp->basisstat.optimal == 0) {
	    ILL_RETURN (1, "mpf_ILLsimplex_solution");
	}
	for (i = 0; i < lp->nrows; i++)
	    mpf_EGlpNumCopy (xz[lp->baz[i]], lp->xbz[i]);
	for (j = 0; j < lp->nnbasic; j++) {
	    col = lp->nbaz[j];
	    if (lp->vstat[col] == STAT_UPPER)
		mpf_EGlpNumCopy (xz[col], lp->uz[col]);
	    else if (lp->vstat[col] == STAT_LOWER)
		mpf_EGlpNumCopy (xz[col], lp->lz[col]);
	    else
		mpf_EGlpNumZero (xz[col]);
	}
    }
    if (piz != NULL) {
	if (lp->basisstat.optimal == 0) {
	    ILL_RETURN (1, "mpf_ILLsimplex_solution");
	}
	for (i = 0; i < lp->nrows; i++)
	    mpf_EGlpNumCopy (piz[i], lp->piz[i]);
    }
    if (dz != NULL) {
	if (lp->basisstat.optimal == 0) {
	    ILL_RETURN (1, "mpf_ILLsimplex_solution");
	}
	for (i = 0; i < lp->nrows; i++)
	    mpf_EGlpNumZero (dz[lp->baz[i]]);
	for (j = 0; j < lp->nnbasic; j++)
	    mpf_EGlpNumCopy (dz[lp->nbaz[j]], lp->dz[j]);
    }
    if (objval != NULL)
	mpf_EGlpNumCopy (*objval, lp->objval);
    return 0;
}

int mpf_ILLsimplex_infcertificate (mpf_lpinfo * lp,
      mpf_t * pi)
{
    int i, col, nz;
    char *sense;
    mpf_t *x, *l, *u;
    mpf_lp_status_info *ls;

    if (pi == NULL)
	return 0;

    ls = &(lp->basisstat);
    if (ls->primal_infeasible == 0 && ls->dual_unbounded == 0) {
	ILL_RETURN (1, "mpf_ILLsimplex_infcertificate");
    }
    if (lp->final_phase == PRIMAL_PHASEI && lp->pIpiz != NULL) {
	for (i = 0; i < lp->nrows; i++)
	    mpf_EGlpNumCopy (pi[i], lp->pIpiz[i]);
    } else if (lp->final_phase == DUAL_PHASEII && lp->infub_ix != -1) {
	col = lp->baz[lp->infub_ix];
	x = &(lp->xbz[lp->infub_ix]);
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);

	for (i = 0; i < lp->nrows; i++)
	    mpf_EGlpNumZero (pi[i]);

	if (mpf_EGlpNumIsNeqq (*l, mpf_NINFTY) && mpf_EGlpNumIsLess (*x, *l)) {
	    for (i = 0, nz = lp->zz.nzcnt; i < nz; i++)
		mpf_EGlpNumCopyNeg (pi[lp->zz.indx[i]], lp->zz.coef[i]);
	} else {
	    for (i = 0, nz = lp->zz.nzcnt; i < nz; i++)
		mpf_EGlpNumCopy (pi[lp->zz.indx[i]], lp->zz.coef[i]);
	}
    } else {
	fprintf (stderr, "Invalid call to inf. certificate routine\n");
	ILL_RETURN (1, "mpf_ILLsimplex_infcertificate");
    }

    sense = lp->O->sense;
    for (i = 0; i < lp->nrows; i++) {
	if (sense[i] == 'G' && mpf_EGlpNumIsLess (pi[i], mpf_zeroLpNum))
	    mpf_EGlpNumZero (pi[i]);
	if (sense[i] == 'L' && mpf_EGlpNumIsLess (mpf_zeroLpNum, pi[i]))
	    mpf_EGlpNumZero (pi[i]);
    }
    return 0;
}

#if SIMPLEX_DEBUG > 1
static void mpf_test_cert (mpf_lpinfo * lp,
      mpf_t * pi)
{
    int i, j;
    int mcnt, mbeg;
    mpf_t fsum, sum;
    mpf_EGlpNumInitVar (fsum);
    mpf_EGlpNumInitVar (sum);
    mpf_EGlpNumZero (fsum);

    for (i = 0; i < lp->nrows; i++) {
	if (lp->O->sense[i] == 'G' && mpf_EGlpNumIsLess (pi[i], mpf_zeroLpNum))
	    printf ("compl \n");
	if (lp->O->sense[i] == 'L' && mpf_EGlpNumIsLess (mpf_zeroLpNum, pi[i]))
	    printf ("compll \n");
    }

    for (i = 0; i < lp->nrows; i++)
	mpf_EGlpNumAddInnProdTo (fsum, pi[i], lp->bz[i]);

    for (j = 0; j < lp->nnbasic; j++) {
	mpf_EGlpNumZero (sum);
	mcnt = lp->matcnt[j];
	mbeg = lp->matbeg[j];
	for (i = 0; i < mcnt; i++)
	    mpf_EGlpNumAddInnProdTo (sum, pi[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

	if (mpf_EGlpNumIsLess (mpf_PFEAS_TOLER, sum) &&
	    (lp->vtype[j] == VLOWER || lp->vtype[j] == VFREE))
	    printf ("compl2\n");
	else {
	    mpf_EGlpNumSign (sum);
	    if (mpf_EGlpNumIsLess (mpf_PFEAS_TOLER, sum) &&
		(lp->vtype[j] == VUPPER || lp->vtype[j] == VFREE))
		printf ("compl1\n");
	    mpf_EGlpNumSign (sum);
	}

	if (mpf_EGlpNumIsLess (sum, mpf_zeroLpNum)
	    && (lp->vtype[j] & (VFREE | VUPPER)) == 0)
	    mpf_EGlpNumSubInnProdTo (fsum, sum, lp->lz[j]);
	else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, sum)
	    && (lp->vtype[j] & (VFREE | VLOWER)) == 0)
	    mpf_EGlpNumSubInnProdTo (fsum, sum, lp->uz[j]);
    }
    printf ("fsum = %.8f\n", mpf_EGlpNumToLf (fsum));
    mpf_EGlpNumClearVar (fsum);
    mpf_EGlpNumClearVar (sum);
}
#endif

static void mpf_save_paraminfo (mpf_price_info * pinf,
      mpf_iter_info * it)
{
    mpf_param_info *pr = &(it->oldinfo);

    pr->origalgo = it->algorithm;
    pr->pphaseI = pinf->pI_price;
    pr->pphaseII = pinf->pII_price;
    pr->dphaseI = pinf->dI_price;
    pr->dphaseII = pinf->dII_price;
    pr->p_strategy = pinf->p_strategy;
    pr->d_strategy = pinf->d_strategy;
}

static void mpf_restore_paraminfo (mpf_iter_info * it,
      mpf_price_info * pinf)
{
    mpf_param_info *pr = &(it->oldinfo);

    it->algorithm = pr->origalgo;
    pinf->pI_price = pr->pphaseI;
    pinf->pII_price = pr->pphaseII;
    pinf->dI_price = pr->dphaseI;
    pinf->dII_price = pr->dphaseII;
    pinf->p_strategy = pr->p_strategy;
    pinf->d_strategy = pr->d_strategy;
}

int mpf_ILLsimplex (mpf_lpinfo * lp,
      int algorithm,
      mpf_ILLlp_basis * B,
      mpf_price_info * pinf,
      int *status,
      int sdisplay)
{
    int phase = -1;
    int singular = -1;
    int rval = 0;
    int new_price = -1;
    mpf_svector wz;
    mpf_svector updz;
    mpf_feas_info fi;
    mpf_iter_info it;
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumInitVar (it.prevobj);
    mpf_EGlpNumInitVar (it.objtol);

    it.newphase = -1;
    it.nextphase = -1;
    it.nextstep = -1;
    it.sdisplay = sdisplay;
    it.n_pivot_fail = 0;
    it.itercnt = 0;
    it.n_restart = 0;
    it.solstatus = ILL_LP_UNSOLVED;
    it.curtime = 0;
    it.rounds = 0;
    mpf_EGlpNumCopy (it.prevobj, mpf_INFTY);
    it.nosolve = 0;
    it.noprog = 0;
    mpf_EGlpNumCopy (it.objtol, mpf_OBJBND_TOLER);
    it.chkobj = PARAM_MAX_NOPROG;
    it.inner = 0;
    it.algorithm = algorithm;
    it.pricetype = -1;
    it.resumeid = -1;
    mpf_save_paraminfo (pinf, &it);

#if SIMPLEX_DEBUG > 0
    if (lp->O->nrows > 1000)
	it.sdisplay = 1;
#endif
    if (status)
	*status = QS_LP_UNSOLVED;

    mpf_free_internal_lpinfo (lp);
    mpf_init_internal_lpinfo (lp);
    rval = mpf_build_internal_lpinfo (lp);
    ILL_CLEANUP_IF (rval);

    mpf_ILLsvector_init (&wz);
    rval = mpf_ILLsvector_alloc (&wz, lp->nrows);
    ILL_CLEANUP_IF (rval);
    mpf_ILLsvector_init (&updz);
    rval = mpf_ILLsvector_alloc (&updz, lp->nrows);
    ILL_CLEANUP_IF (rval);

    if (it.sdisplay) {
	char buffer[256];
	int nonzero = 0;
	register int i = lp->ncols;
	while (i--)
	    nonzero += lp->matcnt[i];
	sprintf (buffer, "starting mpf_ILLsimplex on %s...\n", lp->O->probname);
	/* depending on LP's reporter string is printed to stdout or handed
	   to GUI */
	rval = rval || ILLstring_report (buffer, &lp->O->reporter);
	printf ("Problem has %d rows and %d cols and %d nonzeros\n", lp->nrows, lp->ncols, nonzero);
	fflush (stdout);
    }
    mpf_ILLfct_set_variable_type (lp);

    if (B != 0) {
	rval = mpf_ILLbasis_load (lp, B);
	ILL_CLEANUP_IF (rval);
	if (it.algorithm == DUAL_SIMPLEX) {
	    if (B->rownorms) {
		rval = mpf_ILLprice_load_rownorms (lp, B->rownorms, pinf);
		ILL_CLEANUP_IF (rval);
	    } else
		mpf_EGlpNumFreeArray (pinf->dsinfo.norms);
	} else if (it.algorithm == PRIMAL_SIMPLEX) {
	    if (B->colnorms) {
		rval = mpf_ILLprice_load_colnorms (lp, B->colnorms, pinf);
		ILL_CLEANUP_IF (rval);
	    } else
		mpf_EGlpNumFreeArray (pinf->psinfo.norms);
	} else if (it.algorithm != PRIMAL_OR_DUAL) {
	    fprintf (stderr, "Unknown algorithm %d in mpf_ILLsimplex\n", it.algorithm);
	    rval = 1;
	    ILL_CLEANUP;
	}
    } else if (lp->basisid == -1) {
	if (lp->nrows < 200 && lp->ncols < 400)
	    rval = mpf_ILLbasis_get_initial (lp, it.algorithm);
	else
	    rval = mpf_ILLbasis_get_cinitial (lp, it.algorithm);
	ILL_CLEANUP_IF (rval);
	mpf_ILLprice_free_pricing_info (pinf);
    }
    if (lp->fbasisid != lp->basisid) {
	rval = mpf_ILLbasis_factor (lp, &singular);
	ILL_CLEANUP_IF (rval);
	if (singular)
	    mpf_ILLprice_free_pricing_info (pinf);
    }
START:
#if 0
    if (it.resumeid == SIMPLEX_RESUME_UNSHIFT)
	fprintf (stderr, "Resuming Unshift\n");
    else if (it.resumeid == SIMPLEX_RESUME_SING)
	fprintf (stderr, "Resuming Singular\n");
    else if (it.resumeid == SIMPLEX_RESUME_NUMER)
	fprintf (stderr, "Resuming Numer\n");
    else if (it.resumeid != -1)
	fprintf (stderr, "Resuming for other reason... %d\n", it.resumeid);
#endif
    it.solstatus = ILL_LP_UNSOLVED;
    mpf_init_lp_status_info (&(lp->basisstat));

    mpf_ILLfct_compute_piz (lp);
    mpf_ILLfct_compute_dz (lp);
    if (it.algorithm == DUAL_SIMPLEX) {
	if (B != NULL || it.resumeid == SIMPLEX_RESUME_UNSHIFT)
	    mpf_ILLfct_dual_adjust (lp, lp->tol->dfeas_tol);
	else
	    mpf_ILLfct_dual_adjust (lp, mpf_zeroLpNum);
    }
    mpf_ILLfct_compute_xbz (lp);

    mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
    mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

    if (fi.dstatus == DUAL_FEASIBLE && mpf_EGlpNumIsNeqq (lp->objbound, mpf_INFTY)) {
	mpf_ILLfct_compute_dobj (lp);
	if (mpf_EGlpNumIsLess (lp->objbound, lp->dobjval)) {
	    it.solstatus = ILL_BND_REACHED;
	    printf ("solstatus = ILL_BND_REACHED = 5 %lf %lf\n",
		mpf_EGlpNumToLf (lp->objbound), mpf_EGlpNumToLf (lp->dobjval));
	    goto TERMINATE;
	}
    }
    if (fi.pstatus == PRIMAL_FEASIBLE && fi.dstatus == DUAL_FEASIBLE) {
	it.solstatus = ILL_LP_SOLVED;
	mpf_ILLfct_compute_pobj (lp);
	goto TERMINATE;
    }
    if (it.algorithm == PRIMAL_OR_DUAL) {
	if (fi.pstatus == PRIMAL_FEASIBLE)
	    it.algorithm = PRIMAL_SIMPLEX;
	else if (fi.dstatus == DUAL_FEASIBLE)
	    it.algorithm = DUAL_SIMPLEX;
	else if (mpf_EGlpNumToLf (lp->pinfeas) < 10 * mpf_EGlpNumToLf (lp->dinfeas))
	    it.algorithm = PRIMAL_SIMPLEX;
	else
	    it.algorithm = DUAL_SIMPLEX;
    }
    if (it.algorithm == PRIMAL_SIMPLEX) {
	if (fi.pstatus == PRIMAL_FEASIBLE)
	    phase = PRIMAL_PHASEII;
	else
	    phase = PRIMAL_PHASEI;
    } else if (it.algorithm == DUAL_SIMPLEX) {
	if (fi.dstatus == DUAL_FEASIBLE)
	    phase = DUAL_PHASEII;
	else
	    phase = DUAL_PHASEI;
    }
    rval = mpf_ILLprice_build_pricing_info (lp, pinf, phase);
    ILL_CLEANUP_IF (rval);

    it.newphase = SIMPLEX_PHASE_NEW;
    it.nextstep = SIMPLEX_CONTINUE;

    while (it.nextstep == SIMPLEX_CONTINUE) {
	/* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */

	if (phase == PRIMAL_PHASEI) {
	    rval = mpf_primal_phaseI_step (lp, pinf, &updz, &wz, &it);
	    ILL_CLEANUP_IF (rval);
	} else if (phase == PRIMAL_PHASEII) {
	    rval = mpf_primal_phaseII_step (lp, pinf, &updz, &wz, &it);
	    ILL_CLEANUP_IF (rval);
	} else if (phase == DUAL_PHASEI) {
	    rval = mpf_dual_phaseI_step (lp, pinf, &updz, &wz, &it);
	    ILL_CLEANUP_IF (rval);
	} else if (phase == DUAL_PHASEII) {
	    rval = mpf_dual_phaseII_step (lp, pinf, &updz, &wz, &it);
	    ILL_CLEANUP_IF (rval);
	}
	if (it.nextstep == SIMPLEX_RESUME) {
	    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    mpf_ILLprice_free_pricing_info (pinf);
	    if (it.resumeid == SIMPLEX_RESUME_UNSHIFT) {
		if (it.pricetype == QS_PRICE_PDEVEX) {
		    pinf->pI_price = QS_PRICE_PDEVEX;
		    pinf->pII_price = QS_PRICE_PDEVEX;
		} else if (it.pricetype == QS_PRICE_DDEVEX) {
		    pinf->dI_price = QS_PRICE_DDEVEX;
		    pinf->dII_price = QS_PRICE_DDEVEX;
		}
	    } else if (it.resumeid == SIMPLEX_RESUME_NUMER) {
		mpf_ILLfct_unroll_bound_change (lp);
		mpf_ILLfct_unroll_coef_change (lp);
		/* we are disabling re-do under this circunstances ! */
		rval = mpf_ILLbasis_get_initial (lp, it.algorithm);
		ILL_CLEANUP_IF (rval);
		rval = mpf_ILLbasis_factor (lp, &singular);
		ILL_CLEANUP_IF (rval);
	    }
	    it.pricetype = -1;
	    if (it.n_restart > SIMPLEX_MAX_RESTART) {
		it.solstatus = ILL_MAX_ITER;
		goto LIMIT_TERMINATE;
	    }
	    goto START;
	} else if (it.nextstep == SIMPLEX_CONTINUE) {
	    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    it.itercnt++;

	    if (it.nextphase != phase) {
		/* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
		it.newphase = SIMPLEX_PHASE_NEW;
		phase = it.nextphase;
		new_price = mpf_ILLprice_get_price (pinf, phase);

		if (pinf->cur_price != new_price) {
		    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		       */
		    mpf_ILLprice_free_pricing_info (pinf);
		    rval = mpf_ILLprice_build_pricing_info (lp, pinf, phase);
		    ILL_CLEANUP_IF (rval);
		}
	    }
	}
    }

#if SIMPLEX_DEBUG > 0
    mpf_ILLfct_print_counts (lp);
#endif
LIMIT_TERMINATE:
    rval = mpf_terminate_simplex (lp, phase, &it);
    ILL_CLEANUP_IF (rval);

TERMINATE:
    mpf_restore_paraminfo (&it, pinf);

    if (it.sdisplay) {
	printf ("completed mpf_ILLsimplex\n");
	printf ("%s: ", lp->O->probname);
	fflush (stdout);
    }
    if (status) {
	if (it.solstatus == ILL_MAX_ITER) {
	    *status = QS_LP_ITER_LIMIT;
	} else if (it.solstatus == ILL_MAX_TIME) {
	    *status = QS_LP_TIME_LIMIT;
	} else if (it.solstatus == ILL_LP_ABORTED) {
	    *status = QS_LP_ABORTED;
	} else if (it.solstatus == ILL_PPHASEI_ERROR ||
		it.solstatus == ILL_PPHASEII_ERROR ||
		it.solstatus == ILL_DPHASEI_ERROR ||
		it.solstatus == ILL_DPHASEII_ERROR ||
	    it.solstatus == ILL_LP_UNSOLVED) {
	    *status = QS_LP_UNSOLVED;
	} else if (it.solstatus == ILL_LP_SOLVED) {
	    if (lp->basisstat.optimal) {
		*status = QS_LP_OPTIMAL;
	    } else if (lp->basisstat.primal_infeasible || lp->basisstat.dual_unbounded) {
		*status = QS_LP_INFEASIBLE;
		if (it.sdisplay) {
		    if (lp->basisstat.primal_infeasible)
			fprintf (stdout, "Primal Infeasible\n");
		    else
			fprintf (stdout, "Dual Unbounded\n");
		}
	    } else if (lp->basisstat.primal_unbounded) {
		*status = QS_LP_UNBOUNDED;
	    }
	} else {
	    fprintf (stderr, "unknown solution status in mpf_ILLsimplex %d\n",
		it.solstatus);
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    }
#if SIMPLEX_DEBUG > 1
    {
	int rva = 0;
	mpf_t *pi = NULL;
	pi = mpf_EGlpNumAllocArray (lp->nrows);
	rva = mpf_ILLsimplex_infcertificate (lp, pi);
	printf ("rva = %d\n", rva);
	if (!rva) {
	    mpf_test_cert (lp, pi);
	}
	mpf_EGlpNumFreeArray (pi);
    }
#endif
    if (it.sdisplay) {
	int bstat = 0;

	printf ("time = %.3f, pI = %d, pII = %d, dI = %d, dII = %d, ",
	    ILLutil_zeit () - lp->starttime, lp->cnts->pI_iter,
	    lp->cnts->pII_iter, lp->cnts->dI_iter, lp->cnts->dII_iter);
	fflush (stdout);
	mpf_get_current_stat (&(lp->basisstat), it.algorithm, &bstat);
	switch (bstat) {
	case OPTIMAL:
	    printf ("opt = %f\n", mpf_EGlpNumToLf (lp->objval));
	    break;
	case PRIMAL_INFEASIBLE:
	    printf ("no primal soln\n");
	    break;
	case PRIMAL_UNBOUNDED:
	    printf ("primal unbounded\n");
	    break;
	case PRIMAL_FEASIBLE:
	    printf ("primal obj = %f\n", mpf_EGlpNumToLf (lp->pobjval));
	    break;
	case DUAL_INFEASIBLE:
	    printf ("no dual soln\n");
	    break;
	case DUAL_UNBOUNDED:
	    printf ("dual unbounded\n");
	    break;
	case DUAL_FEASIBLE:
	    printf ("dual obj = %f\n", mpf_EGlpNumToLf (lp->dobjval));
	    break;
	}
	fflush (stdout);

	if (it.sdisplay > 1) {
	    if (it.algorithm == PRIMAL_SIMPLEX && pinf->pI_price == QS_PRICE_PDEVEX)
		printf ("Devex norms initialised %d times\n", pinf->pdinfo.ninit);
	    fflush (stdout);
	}
    }
CLEANUP:
    mpf_ILLsvector_free (&wz);
    mpf_ILLsvector_free (&updz);
    mpf_EGlpNumClearVar (it.prevobj);
    mpf_EGlpNumClearVar (it.objtol);
    mpf_EGlpNumClearVar (fi.totinfeas);
    ILL_RETURN (rval, "mpf_ILLsimplex");
}

static int mpf_terminate_simplex (mpf_lpinfo * lp,
      int phase,
      mpf_iter_info * it)
{
    int rval = 0;
    int sphase;
    mpf_feas_info fi;
    mpf_EGlpNumInitVar (fi.totinfeas);

    if (it->solstatus != ILL_MAX_TIME && it->solstatus != ILL_MAX_ITER)
	ILL_CLEANUP;

    if (it->algorithm == PRIMAL_SIMPLEX) {
	if (lp->nbchange != 0) {
	    if (it->sdisplay > 1) {
		printf ("unrolling %d bound shifts\n", lp->nbchange);
		fflush (stdout);
	    }
	    mpf_ILLfct_unroll_bound_change (lp);
	}
	rval = mpf_ILLsimplex_retest_psolution (lp, NULL, phase, &fi);
	ILL_CLEANUP_IF (rval);

	sphase = (phase == PRIMAL_PHASEI) ? PHASEI : PHASEII;
	mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, sphase);
    } else if (it->algorithm == DUAL_SIMPLEX) {
	if (lp->ncchange != 0) {
	    if (it->sdisplay > 1) {
		printf ("unrolling %d coef shifts\n", lp->ncchange);
		fflush (stdout);
	    }
	    mpf_ILLfct_unroll_coef_change (lp);
	}
	rval = mpf_ILLsimplex_retest_dsolution (lp, NULL, phase, &fi);
	ILL_CLEANUP_IF (rval);

	sphase = (phase == DUAL_PHASEI) ? PHASEI : PHASEII;
	mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, sphase, PHASEII);
    }
CLEANUP:
    mpf_EGlpNumClearVar (fi.totinfeas);
    ILL_RETURN (rval, "mpf_terminate_simplex");
}

static int mpf_test_progress (mpf_t objval,
      mpf_t prevobj)
{
    mpf_t denom;
    mpf_EGlpNumInitVar (denom);
    mpf_EGlpNumCopyDiff (denom, objval, prevobj);
    if (mpf_EGlpNumIsNeqZero (objval, mpf_PROGRESS_ZERO))
	mpf_EGlpNumDivTo (denom, objval);
    if (mpf_EGlpNumIsEqual (denom, mpf_zeroLpNum, mpf_PROGRESS_THRESH)) {
	mpf_EGlpNumClearVar (denom);
	return 0;
    } else {
	mpf_EGlpNumClearVar (denom);
	return 1;
    }
}

static void mpf_monitor_iter (mpf_lpinfo * lp, mpf_iter_info * it, int phase)
{
    mpf_t print_val;
    double tottime = ILLutil_zeit () - lp->starttime;
    int curtime = mpf_ILLutil_our_floor (tottime);	/* MONIKA */
    char print_str[20];
    mpf_feas_info fi;
    int aborted = 0;
    mpf_EGlpNumInitVar (print_val);
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumZero (print_val);

    /* one of the following two time display mechanisms */
    switch (phase) {
    case PRIMAL_PHASEI:
	mpf_EGlpNumCopy (print_val, lp->pinfeas);
	strcpy (print_str, "primal infeas");
	if (mpf_EGlpNumIsLess (lp->pinfeas, mpf_zeroLpNum) &&
	    (mpf_EGlpNumIsNeqZero (lp->pinfeas, mpf_oneLpNum))) {
	    /* printf ("Negative Infeasibility! Imposible %lg %la, iter
	       %d\n", mpf_EGlpNumToLf (print_val), mpf_EGlpNumToLf
	       (print_val), it->itercnt); */
	    /* exit(1); */
	}
	break;
    case PRIMAL_PHASEII:
	mpf_EGlpNumCopy (print_val, lp->pobjval);
	strcpy (print_str, "primal objval");
	break;
    case DUAL_PHASEI:
	mpf_EGlpNumCopy (print_val, lp->dinfeas);
	strcpy (print_str, "dual infeas");
	break;
    case DUAL_PHASEII:
	mpf_EGlpNumCopy (print_val, lp->dobjval);
	strcpy (print_str, "dual objval");
	break;
    }

    aborted = mpf_report_value (lp, it, print_str, print_val);
    /* if (it->sdisplay && it->itercnt % lp->iterskip == 0) { printf ("(%d):
       %s = %f\n", it->itercnt, print_str, print_val); fflush (stdout); } */
    if (curtime != it->curtime) {
	it->curtime = curtime;
	/*
         * if (it->sdisplay){
         * printf ("time = %d.0, ", curtime);
         * printf ("(%d): %s = %f\n", it->itercnt, print_str, print_val);
         * fflush (stdout);
         * }
         */
    }
    if (phase == DUAL_PHASEII && mpf_EGlpNumIsNeqq (lp->objbound, mpf_INFTY)) {
	/* if (lp->dobjval > lp->objbound + it->objtol) */
	mpf_EGlpNumCopyDiff (print_val, lp->dobjval, lp->objbound);
	if (mpf_EGlpNumIsLess (it->objtol, print_val)) {
	    mpf_ILLfct_unroll_coef_change (lp);
	    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	    mpf_ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

	    if (fi.dstatus == DUAL_FEASIBLE) {
		mpf_ILLfct_compute_dobj (lp);
		if (mpf_EGlpNumIsLess (lp->objbound, lp->dobjval)) {
		    it->solstatus = ILL_BND_REACHED;
		    it->nextstep = SIMPLEX_TERMINATE;
		    /* if (it->sdisplay) */
		    {
			printf ("bound reached %lf %lf\n", mpf_EGlpNumToLf (lp->objbound),
			    mpf_EGlpNumToLf (lp->dobjval));
			fflush (stdout);
		    }
		} else
		    mpf_EGlpNumMultUiTo (it->objtol, 10);
	    } else {
		it->nextphase = DUAL_PHASEI;
		it->newphase = SIMPLEX_PHASE_NEW;
		mpf_EGlpNumMultUiTo (it->objtol, 5);
	    }
	}
    }
    if (it->itercnt >= lp->maxiter) {
	it->solstatus = ILL_MAX_ITER;
	it->nextstep = SIMPLEX_TERMINATE;
	if (it->sdisplay) {
	    printf ("iter limit reached\n");
	    fflush (stdout);
	}
	ILL_CLEANUP;
    } else if (tottime >= lp->maxtime) {
	it->solstatus = ILL_MAX_TIME;
	it->nextstep = SIMPLEX_TERMINATE;
	if (it->sdisplay) {
	    printf ("time limit reached\n");
	    fflush (stdout);
	}
	ILL_CLEANUP;
    } else if (aborted) {
	it->solstatus = ILL_LP_ABORTED;
	it->nextstep = SIMPLEX_TERMINATE;
	if (it->sdisplay) {
	    printf ("aborted\n");
	    fflush (stdout);
	}
	ILL_CLEANUP;
    }
    /*
     * if (it->rounds && it->inner){
     * it->inner --;
     * if (it->inner == 0){
     * printf ("restoring ..\n");
     * mpf_restore_paraminfo (it, p);
     * it->newphase   = SIMPLEX_PHASE_NEW;
     * it->nextstep   = SIMPLEX_RESUME;
     * it->resumeid   = SIMPLEX_RESUME_OUTER;
     * ILL_CLEANUP;
     * }
     * }
     */
    if (phase == DUAL_PHASEII) {
	if (it->noprog > it->chkobj) {
	    mpf_ILLfct_perturb_coefs (lp);
	    it->noprog = 0;
	    mpf_EGlpNumCopy (it->prevobj, lp->dobjval);
	}
    } else if (phase == PRIMAL_PHASEII) {
	if (it->noprog > it->chkobj) {
	    mpf_ILLfct_perturb_bounds (lp);
	    it->noprog = 0;
	    mpf_EGlpNumCopy (it->prevobj, lp->pobjval);
	}
    } else if (phase == PRIMAL_PHASEI) {
	if (it->noprog > it->chkobj) {
	    it->algorithm = DUAL_SIMPLEX;
	    it->nextstep = SIMPLEX_RESUME;
	    it->resumeid = SIMPLEX_RESUME_NUMER;
	    /* this is to force to exit in the case of bad basis */
	    mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	    mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	    it->n_restart++;
	    /* fprintf(stderr,"Resume Numerical
	       %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	}
    } else if (phase == DUAL_PHASEI) {
	if (it->noprog > it->chkobj) {
	    it->algorithm = PRIMAL_SIMPLEX;
	    it->nextstep = SIMPLEX_RESUME;
	    it->resumeid = SIMPLEX_RESUME_NUMER;
	    /* this is to force to exit in the case of bad basis */
	    mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	    mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	    it->n_restart++;
	    /* fprintf(stderr,"Resume Numerical
	       %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	}
    }
CLEANUP:
    mpf_EGlpNumClearVar (fi.totinfeas);
    mpf_EGlpNumClearVar (print_val);
    return;
}

static int mpf_primal_phaseI_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it)
{
    int rval = 0;
    int singular = 0;
    int refactor = 0;
    int cphase = PRIMAL_PHASEI;
    mpf_t alpha;
    mpf_feas_info fi;
    mpf_ratio_res rs;
    mpf_price_res pr;
    mpf_EGlpNumInitVar (alpha);
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumInitVar (pr.dinfeas);
    mpf_EGlpNumInitVar (pr.pinfeas);
    mpf_EGlpNumInitVar (rs.tz);
    mpf_EGlpNumInitVar (rs.lbound);
    mpf_EGlpNumInitVar (rs.ecoeff);
    mpf_EGlpNumInitVar (rs.pivotval);
    mpf_EGlpNumZero (alpha);

    mpf_ILLfct_update_counts (lp, CNT_PPHASE1ITER, 0, mpf_zeroLpNum);
    it->nextstep = SIMPLEX_CONTINUE;
    it->nextphase = PRIMAL_PHASEI;
    lp->final_phase = PRIMAL_PHASEI;
    it->nosolve++;

    if (it->newphase != 0) {
	mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
	if (it->newphase == SIMPLEX_PHASE_NEW) {
	    it->noprog = 0;
	    if (it->sdisplay) {
		printf ("starting primal phase I\n");
		fflush (stdout);
	    }
	}
	it->newphase = 0;
	it->nosolve = 0;
	mpf_EGlpNumCopy (it->prevobj, lp->pinfeas);
	lp->pIpiz = mpf_EGlpNumAllocArray (lp->nrows);
	lp->pIdz = mpf_EGlpNumAllocArray (lp->nnbasic);

	mpf_ILLfct_compute_phaseI_piz (lp);
	if (pinf->p_strategy == COMPLETE_PRICING) {
	    mpf_ILLfct_compute_phaseI_dz (lp);
#if USEHEAP > 0
	    mpf_ILLprice_free_heap (pinf);
#endif
	    mpf_ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEI);
#if USEHEAP > 0
	    rval = mpf_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
		PRIMAL_SIMPLEX, 0);
	    ILL_CLEANUP_IF (rval);
#endif
	} else if (pinf->p_strategy == MULTI_PART_PRICING)
	    mpf_ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
    }
    mpf_monitor_iter (lp, it, cphase);
    if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
	it->newphase != 0)
	ILL_CLEANUP;

    mpf_ILLprice_primal (lp, pinf, &pr, cphase);
    ILL_IFTRACE2 ("%s:after_price\n", __func__);

    if (pr.price_stat == PRICE_OPTIMAL) {
	if (it->sdisplay > 1) {
	    printf ("primal phase I seemingly done\n");
	    printf ("retesting soln\n");
	    fflush (stdout);
	}
	rval = mpf_ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);

	ILL_CLEANUP_IF (rval);
	mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEI);

	if (fi.pstatus == PRIMAL_FEASIBLE) {
	    it->nextphase = PRIMAL_PHASEII;
	} else if (fi.dstatus == DUAL_FEASIBLE) {
	    it->solstatus = ILL_LP_SOLVED;
	    it->nextstep = SIMPLEX_TERMINATE;
	}
	ILL_CLEANUP;
    }
    mpf_ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
    mpf_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, mpf_zeroLpNum);
    mpf_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, mpf_zeroLpNum);

    mpf_ILLratio_pI_test (lp, pr.eindex, pr.dir, &rs);
    /* ILL_IFTRACE(":%d",rs.lindex); */

    if (rs.ratio_stat == RATIO_FAILED) {
	/*
         * rval = E_SIMPLEX_ERROR;
         * it->solstatus = ILL_PPHASEI_ERROR;
         */
	/* ILL_IFTRACE("ratio_failed\n"); */
	it->algorithm = DUAL_SIMPLEX;
	it->nextstep = SIMPLEX_RESUME;
	it->resumeid = SIMPLEX_RESUME_NUMER;
	/* this is to force to exit in the case of bad basis */
	it->n_restart++;
	mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	/* fprintf(stderr,"Resume Numerical
	   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	ILL_CLEANUP;
    } else if (rs.ratio_stat == RATIO_NEGATIVE) {
	mpf_t itol;
	mpf_EGlpNumInitVar (itol);
	/* ILL_IFTRACE("ratio_negative\n"); */
	mpf_EGlpNumCopy (itol, lp->tol->ip_tol);
	mpf_EGlpNumZero (lp->tol->ip_tol);
	mpf_EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);
	if (!mpf_test_progress (lp->pinfeas, it->prevobj))
	    it->noprog++;
	else {
	    mpf_EGlpNumCopy (it->prevobj, lp->pinfeas);
	    it->noprog = 0;
	}
	mpf_ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
	mpf_EGlpNumCopy (lp->tol->ip_tol, itol);
	mpf_ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
	mpf_ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
	    rs.lindex, mpf_zeroLpNum);
	mpf_EGlpNumClearVar (itol);
    } else if (rs.ratio_stat == RATIO_NOBCHANGE) {
	/* ILL_IFTRACE("ratio_nobchange\n"); */
	mpf_EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);
	if (!mpf_test_progress (lp->pinfeas, it->prevobj))
	    it->noprog++;
	else {
	    mpf_EGlpNumCopy (it->prevobj, lp->pinfeas);
	    it->noprog = 0;
	}

	/* ILL_IFTRACE("%s:a\n",__func__); */
	mpf_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
	mpf_ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
	mpf_ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
	mpf_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
	mpf_fct_test_workvector (lp);
	mpf_fct_test_pfeasible (lp);
#endif
	mpf_ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
	    rs.lindex, mpf_zeroLpNum);
    } else if (rs.ratio_stat == RATIO_BCHANGE) {
	/* ILL_IFTRACE("ratio_bchange\n"); */
	mpf_EGlpNumCopyFrac (alpha, lp->pIdz[pr.eindex], rs.pivotval);
	mpf_EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);

	if (!mpf_test_progress (lp->pinfeas, it->prevobj)) {
	    if (lp->vtype[lp->nbaz[pr.eindex]] == VFREE ||
		lp->vtype[lp->baz[rs.lindex]] == VARTIFICIAL) {
		if (it->noprog > 0)
		    it->noprog--;
	    } else
		it->noprog++;
	} else {
	    mpf_EGlpNumCopy (it->prevobj, lp->pinfeas);
	    it->noprog = 0;
	}

	mpf_ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
	mpf_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, mpf_zeroLpNum);
	if (pinf->p_strategy == COMPLETE_PRICING) {
	    /* ILL_IFTRACE("%s:a\n",__func__); */
	    mpf_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	    mpf_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, mpf_zeroLpNum);

	    if (pinf->pI_price == QS_PRICE_PSTEEP) {
		mpf_ILLfct_compute_psteep_upv (lp, wz);
	    }
	}
	rval =
	    mpf_ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex,
	    rs.pivotval);
	ILL_CLEANUP_IF (rval);

	/* ILL_IFTRACE("%s:b:%d\n",__func__,rs.lindex); */
	mpf_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
	mpf_ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
	/* ILL_IFTRACE("%s:%d:%d\n",__func__,rs.lindex,lp->srhs.nzcnt); */
	mpf_ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
	mpf_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
	mpf_fct_test_workvector (lp);
	mpf_fct_test_pfeasible (lp);
#endif
	rval = mpf_ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
	ILL_CLEANUP_IF (rval);

	if (singular) {
	    it->nextstep = SIMPLEX_RESUME;
	    it->resumeid = SIMPLEX_RESUME_SING;
	    /* this is to force to exit in the case of bad basis */
	    it->n_restart++;
	    mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	    mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	    /* fprintf(stderr,"Resume Singular
	       %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    ILL_CLEANUP;
	}
	if (!refactor) {
	    mpf_ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
		rs.lindex, alpha);
	}
	if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE) {
	    mpf_ILLfct_compute_xbz (lp);
	    mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
	    mpf_ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);
	    if (fi.pstatus == PRIMAL_FEASIBLE)
		it->nextphase = PRIMAL_PHASEII;

	    it->newphase = SIMPLEX_PHASE_RECOMP;
	    ILL_CLEANUP;
	}
    }
#if DENSE_PI > 1
    mpf_fct_test_workvector (lp);
    fct_test_pi_dz (lp, pinf);
#endif

CLEANUP:
    if (it->nextphase != PRIMAL_PHASEI || it->nextstep == SIMPLEX_RESUME ||
	it->newphase != 0 || rval != 0) {
	mpf_EGlpNumFreeArray (lp->pIpiz);
	mpf_EGlpNumFreeArray (lp->pIdz);
    }
    mpf_EGlpNumClearVar (alpha);
    mpf_EGlpNumClearVar (fi.totinfeas);
    mpf_EGlpNumClearVar (pr.dinfeas);
    mpf_EGlpNumClearVar (pr.pinfeas);
    mpf_EGlpNumClearVar (rs.tz);
    mpf_EGlpNumClearVar (rs.lbound);
    mpf_EGlpNumClearVar (rs.ecoeff);
    mpf_EGlpNumClearVar (rs.pivotval);
    /* ILL_RETURN (rval, "mpf_primal_phaseI_step"); */
    return rval;
}

static int mpf_primal_phaseII_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it)
{
    int boundch;
    int rval = 0;
    int bndtype = 0;
    int singular = 0;
    int refactor = 0;
    int ratio_iter = 0;
    int cphase = PRIMAL_PHASEII;
    mpf_t lbound;
    mpf_t alpha;
    mpf_feas_info fi;
    mpf_ratio_res rs;
    mpf_price_res pr;
    mpf_EGlpNumInitVar (alpha);
    mpf_EGlpNumInitVar (lbound);
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumInitVar (pr.dinfeas);
    mpf_EGlpNumInitVar (pr.pinfeas);
    mpf_EGlpNumInitVar (rs.tz);
    mpf_EGlpNumInitVar (rs.lbound);
    mpf_EGlpNumInitVar (rs.ecoeff);
    mpf_EGlpNumInitVar (rs.pivotval);

    mpf_ILLfct_update_counts (lp, CNT_PPHASE2ITER, 0, mpf_zeroLpNum);
    it->nextstep = SIMPLEX_CONTINUE;
    it->nextphase = PRIMAL_PHASEII;
    lp->final_phase = PRIMAL_PHASEII;
    it->nosolve++;

    if (it->newphase != 0) {
	mpf_ILLfct_compute_pobj (lp);
	if (it->newphase == SIMPLEX_PHASE_NEW) {
	    it->noprog = 0;
	    if (it->sdisplay) {
		printf ("starting primal phase II\n");
		fflush (stdout);
	    }
	}
	it->newphase = 0;
	it->nosolve = 0;
	mpf_EGlpNumCopy (it->prevobj, lp->pobjval);
	mpf_ILLfct_compute_piz (lp);
	if (pinf->p_strategy == COMPLETE_PRICING) {
	    mpf_ILLfct_compute_dz (lp);
#if USEHEAP > 0
	    mpf_ILLprice_free_heap (pinf);
#endif
	    mpf_ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEII);
#if USEHEAP > 0
	    rval = mpf_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
		PRIMAL_SIMPLEX, 0);
	    ILL_CLEANUP_IF (rval);
#endif
	} else if (pinf->p_strategy == MULTI_PART_PRICING) {
	    mpf_ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
	}
    }
    mpf_monitor_iter (lp, it, cphase);
    if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
	it->newphase != 0)
	ILL_CLEANUP;

    mpf_ILLprice_primal (lp, pinf, &pr, cphase);

    if (pr.price_stat == PRICE_OPTIMAL) {
	/* ILL_IFTRACE("%s:PRICE_OPTIMAL\n",__func__); */
	if (lp->nbchange != 0) {
	    if (it->sdisplay > 1) {
		printf ("unrolling %d bound shifts\n", lp->nbchange);
		fflush (stdout);
	    }
	    mpf_ILLfct_unroll_bound_change (lp);
	    mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	    mpf_ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);

	     /* HHH */ mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	    /* HHH* printf ("primal (opt) infeas %.6f\n", lp->pinfeas);
	       fflush (stdout); HHH* printf ("dual (opt) infeas %.6f\n",
	       lp->dinfeas); fflush (stdout); */

	    if (fi.pstatus != PRIMAL_FEASIBLE) {
		it->algorithm = DUAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_UNSHIFT;
		it->pricetype = QS_PRICE_DDEVEX;
		/* this is to force to exit in the case of bad basis */
		/* fprintf(stderr,"Resume Unshift
		   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
		mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		it->n_restart++;
		ILL_CLEANUP;
		/*
                 * it->nextphase = PRIMAL_PHASEI;
                 * lp->tol->ip_tol /= 5.0;
                 * lp->tol->id_tol /= 5.0;
                 * ILL_CLEANUP;
                 */
	    }
	}
	if (it->sdisplay > 1) {
	    printf ("problem seemingly solved\n");
	    printf ("seemingly opt = %f\nretesting soln\n",
		mpf_EGlpNumToLf (lp->pobjval));
	    fflush (stdout);
	}
	rval = mpf_ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);
	ILL_CLEANUP_IF (rval);
	mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

	if (fi.pstatus == PRIMAL_INFEASIBLE) {
	    it->nextphase = PRIMAL_PHASEI;
	    mpf_EGlpNumDivUiTo (lp->tol->ip_tol, 5);
	    mpf_EGlpNumDivUiTo (lp->tol->id_tol, 5);
	    ILL_IFTRACE ("%s:PINF:%lg\n", __func__, mpf_EGlpNumToLf (lp->tol->ip_tol));
	} else if (fi.dstatus == DUAL_FEASIBLE) {
	    /* ILL_IFTRACE("%s:PFEAS_DFEAS\n",__func__); */
	    it->solstatus = ILL_LP_SOLVED;
	    mpf_EGlpNumCopy (lp->objval, lp->pobjval);
	    it->nextstep = SIMPLEX_TERMINATE;
	} else
	    ILL_IFTRACE ("%s:DINF:%la:%lf\n", __func__, mpf_EGlpNumToLf (lp->dinfeas),
		mpf_EGlpNumToLf (lp->dinfeas));
	ILL_CLEANUP;
    }
    mpf_ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
    mpf_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, mpf_zeroLpNum);
    mpf_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, mpf_zeroLpNum);
    ratio_iter = 0;
    do {
	mpf_ILLratio_pII_test (lp, pr.eindex, pr.dir, &rs);
	/* ILL_IFTRACE("all:%d",rs.lindex); */
	mpf_EGlpNumCopy (lbound, rs.lbound);
	boundch = rs.boundch;
	ratio_iter++;

	if (boundch) {
	    /*
             * if (ratio_iter > PARAM_PRATIOTESTS){
             * lbound = lp->xbz[rs.lindex];
             * boundch = 0;
             * }
             */
	    boundch = 0;
	    bndtype = (rs.lvstat == STAT_UPPER) ? BOUND_UPPER : BOUND_LOWER;
	    rval = mpf_ILLfct_bound_shift (lp, lp->baz[rs.lindex], bndtype, lbound);
	    ILL_CLEANUP_IF (rval);
	}
    } while (boundch);

    if (rs.ratio_stat == RATIO_FAILED) {
	/* ILL_IFTRACE(":%d",rs.lindex); */
	/*
         * rval = E_SIMPLEX_ERROR;
         * it->solstatus = ILL_PPHASEII_ERROR;
         */
	it->algorithm = DUAL_SIMPLEX;
	it->nextstep = SIMPLEX_RESUME;
	it->resumeid = SIMPLEX_RESUME_NUMER;
	/* this is to force to exit in the case of bad basis */
	it->n_restart++;
	mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	/* fprintf(stderr,"Resume Numerical
	   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	ILL_CLEANUP;
    } else if (rs.ratio_stat == RATIO_UNBOUNDED) {
	/* ILL_IFTRACE(":%d",rs.lindex); */
	if (lp->nbchange != 0) {
	    if (it->sdisplay > 1) {
		printf ("unrolling %d bound shifts\n", lp->nbchange);
		fflush (stdout);
	    }
	    mpf_ILLfct_unroll_bound_change (lp);
	}
	mpf_ILLfct_set_status_values (lp, PRIMAL_UNBOUNDED, -1, PHASEII, -1);
	it->solstatus = ILL_LP_SOLVED;
	it->nextstep = SIMPLEX_TERMINATE;
	ILL_CLEANUP;
    } else if (rs.ratio_stat == RATIO_NOBCHANGE) {
	/* ILL_IFTRACE(":%d",rs.lindex); */
	mpf_EGlpNumAddInnProdTo (lp->pobjval, rs.tz, lp->dz[pr.eindex]);
	mpf_EGlpNumCopy (lp->objval, lp->pobjval);
	if (!mpf_test_progress (lp->pobjval, it->prevobj))
	    it->noprog++;
	else {
	    mpf_EGlpNumCopy (it->prevobj, lp->pobjval);
	    it->noprog = 0;
	}

	/* ILL_IFTRACE("%s:c:%d\n",__func__,rs.lindex); */
	mpf_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
	mpf_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
	if (pinf->p_strategy == COMPLETE_PRICING)
	    mpf_ILLprice_compute_dual_inf (lp, pinf, &pr.eindex, 1, PRIMAL_PHASEII);
	else if (pinf->p_strategy == MULTI_PART_PRICING)
	    mpf_ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
    } else if (rs.ratio_stat == RATIO_BCHANGE) {
	mpf_EGlpNumCopyFrac (alpha, lp->dz[pr.eindex], rs.pivotval);
	mpf_EGlpNumAddInnProdTo (lp->pobjval, rs.tz, lp->dz[pr.eindex]);
	mpf_EGlpNumCopy (lp->objval, lp->pobjval);

	if (!mpf_test_progress (lp->pobjval, it->prevobj)) {
	    /* ILL_IFTRACE(":%d",rs.lindex); */
	    if (lp->vtype[lp->nbaz[pr.eindex]] == VFREE ||
		lp->vtype[lp->baz[rs.lindex]] == VARTIFICIAL) {
		if (it->noprog > 0)
		    it->noprog--;
	    } else
		it->noprog++;
	} else {
	    mpf_EGlpNumCopy (it->prevobj, lp->pobjval);
	    it->noprog = 0;
	}

	/* ILL_IFTRACE(":%d",rs.lindex); */
	mpf_ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
	mpf_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, mpf_zeroLpNum);
	if (pinf->p_strategy == COMPLETE_PRICING) {
	    /* ILL_IFTRACE("%s:b\n",__func__); */
	    mpf_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	    mpf_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, mpf_zeroLpNum);
	    if (pinf->pII_price == QS_PRICE_PSTEEP)
		mpf_ILLfct_compute_psteep_upv (lp, wz);
	}
	rval =
	    mpf_ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex,
	    rs.pivotval);
	ILL_CLEANUP_IF (rval);

	/* ILL_IFTRACE("%s:d:%d\n",__func__,rs.lindex); */
	mpf_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
	mpf_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
	rval = mpf_ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
	ILL_CLEANUP_IF (rval);

	if (singular) {
	    it->nextstep = SIMPLEX_RESUME;
	    it->resumeid = SIMPLEX_RESUME_SING;
	    /* this is to force to exit in the case of bad basis */
	    it->n_restart++;
	    /* fprintf(stderr,"Resume Singular
	       %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	    mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	    ILL_CLEANUP;
	}
	if (!refactor) {
	    mpf_ILLfct_update_piz (lp, alpha);

	    if (pinf->p_strategy == COMPLETE_PRICING) {
		mpf_ILLfct_update_dz (lp, pr.eindex, alpha);
		mpf_ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
		    PRIMAL_PHASEII);
		mpf_ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, mpf_zeroLpNum);
	    } else if (pinf->p_strategy == MULTI_PART_PRICING) {
		mpf_ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
	    }
	}
	if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE) {
	    mpf_ILLfct_compute_xbz (lp);
	    it->newphase = SIMPLEX_PHASE_RECOMP;
	}
    }
CLEANUP:
    mpf_EGlpNumClearVar (alpha);
    mpf_EGlpNumClearVar (lbound);
    mpf_EGlpNumClearVar (fi.totinfeas);
    mpf_EGlpNumClearVar (pr.dinfeas);
    mpf_EGlpNumClearVar (pr.pinfeas);
    mpf_EGlpNumClearVar (rs.tz);
    mpf_EGlpNumClearVar (rs.lbound);
    mpf_EGlpNumClearVar (rs.ecoeff);
    mpf_EGlpNumClearVar (rs.pivotval);
    /* ILL_RETURN (rval, "mpf_primal_phaseII_step"); */
    return rval;
}

static int mpf_dual_phaseI_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it)
{
    int rval = 0;
    int singular = 0;
    int refactor = 0;
    int cphase = DUAL_PHASEI;
    mpf_t alpha;
    mpf_t alpha1;
    mpf_feas_info fi;
    mpf_ratio_res rs;
    mpf_price_res pr;
    mpf_EGlpNumInitVar (alpha);
    mpf_EGlpNumInitVar (alpha1);
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumInitVar (pr.dinfeas);
    mpf_EGlpNumInitVar (pr.pinfeas);
    mpf_EGlpNumInitVar (rs.tz);
    mpf_EGlpNumInitVar (rs.lbound);
    mpf_EGlpNumInitVar (rs.ecoeff);
    mpf_EGlpNumInitVar (rs.pivotval);
    mpf_EGlpNumZero (alpha1);

    mpf_ILLfct_update_counts (lp, CNT_DPHASE1ITER, 0, mpf_zeroLpNum);
    it->nextstep = SIMPLEX_CONTINUE;
    it->nextphase = DUAL_PHASEI;
    lp->final_phase = DUAL_PHASEI;
    it->nosolve++;

    if (it->newphase != 0) {
	mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
	if (it->newphase == SIMPLEX_PHASE_NEW) {
	    it->noprog = 0;
	    if (it->sdisplay) {
		printf ("starting dual phase I\n");
		fflush (stdout);
	    }
	}
	it->newphase = 0;
	it->nosolve = 0;
	mpf_EGlpNumCopy (it->prevobj, lp->dinfeas);

	mpf_ILLfct_compute_phaseI_xbz (lp);
	if (pinf->d_strategy == COMPLETE_PRICING) {
#if USEHEAP > 0
	    mpf_ILLprice_free_heap (pinf);
#endif
	    mpf_ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEI);
#if USEHEAP > 0
	    rval = mpf_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
		DUAL_SIMPLEX, 0);
	    ILL_CLEANUP_IF (rval);
#endif
	} else if (pinf->d_strategy == MULTI_PART_PRICING) {
	    mpf_ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
	}
    }
    mpf_monitor_iter (lp, it, cphase);
    if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
	it->newphase != 0)
	ILL_CLEANUP;

    mpf_ILLprice_dual (lp, pinf, cphase, &pr);

    if (pr.price_stat == PRICE_OPTIMAL) {
	if (it->sdisplay > 1) {
	    printf ("dual phase I seemingly done\n");
	    printf ("retesting soln\n");
	    fflush (stdout);
	}
	rval = mpf_ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
	ILL_CLEANUP_IF (rval);
	mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEI, PHASEII);

	if (fi.dstatus == DUAL_FEASIBLE) {
	    it->nextphase = DUAL_PHASEII;
	} else if (fi.pstatus == PRIMAL_FEASIBLE) {
	    it->solstatus = ILL_LP_SOLVED;
	    it->nextstep = SIMPLEX_TERMINATE;
	}
	it->newphase = SIMPLEX_PHASE_NEW;
	ILL_CLEANUP;
    }
    mpf_ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
    /* ILL_IFTRACE("%s:c\n",__func__); */
    mpf_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
    mpf_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, mpf_zeroLpNum);
    mpf_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, mpf_zeroLpNum);

    mpf_ILLratio_dI_test (lp, pr.lindex, pr.lvstat, &rs);

    if (rs.ratio_stat == RATIO_FAILED) {
	/*
         * rval = E_SIMPLEX_ERROR;
         * it->solstatus = ILL_DPHASEI_ERROR;
         */
	it->algorithm = PRIMAL_SIMPLEX;
	it->nextstep = SIMPLEX_RESUME;
	it->resumeid = SIMPLEX_RESUME_NUMER;
	/* this is to force to exit in the case of bad basis */
	it->n_restart++;
	mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	/* fprintf(stderr,"Resume Numerical
	   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	ILL_CLEANUP;
    } else if (rs.ratio_stat == RATIO_BCHANGE) {
	mpf_ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[rs.eindex]);
	rval = mpf_ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
	if (rval) {
	    it->n_pivot_fail++;
	    if (it->n_pivot_fail > SIMPLEX_MAX_PIVOT_FAIL) {
		it->n_pivot_fail = 0;
		/* this is to force to exit in the case of bad basis */
		it->n_restart++;
		it->algorithm = PRIMAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_NUMER;
		mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		/* fprintf(stderr,"Resume Pivot
		   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
		rval = 0;
		ILL_CLEANUP;
	    }
	    rval = mpf_ILLbasis_factor (lp, &singular);
	    ILL_CLEANUP_IF (rval);
	    if (singular == 0)
		refactor = 1;
	    goto END;
	}
	mpf_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, mpf_zeroLpNum);
	mpf_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, mpf_zeroLpNum);

	if (pinf->dI_price == QS_PRICE_DSTEEP)
	    mpf_ILLfct_compute_dsteep_upv (lp, wz);
	rval =
	    mpf_ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex,
	    rs.pivotval);
	ILL_CLEANUP_IF (rval);

	mpf_EGlpNumSubTo (lp->dinfeas, lp->upd.c_obj);

	if (!mpf_test_progress (lp->dinfeas, it->prevobj)) {
	    if (lp->vtype[lp->baz[pr.lindex]] == VARTIFICIAL ||
		lp->vtype[lp->nbaz[rs.eindex]] == VFREE) {
		if (it->noprog > 0)
		    it->noprog--;
	    } else
		it->noprog++;
	} else {
	    mpf_EGlpNumCopy (it->prevobj, lp->dinfeas);
	    it->noprog = 0;
	}

	mpf_EGlpNumCopyFrac (alpha, lp->dz[rs.eindex], rs.pivotval);
	mpf_EGlpNumCopyFrac (alpha1, lp->xbz[pr.lindex], rs.pivotval);

	mpf_ILLfct_update_piz (lp, alpha);
	mpf_ILLfct_update_dz (lp, rs.eindex, alpha);
	mpf_ILLfct_update_dfeas (lp, rs.eindex, &(lp->srhs));
	mpf_ILLfct_compute_dpIy (lp, &(lp->srhs), &(lp->ssoln));
	mpf_ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);

#if DENSE_PI > 0
	mpf_fct_test_workvector (lp);
	mpf_fct_test_dfeasible (lp);
#endif
	rval = mpf_ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
	ILL_CLEANUP_IF (rval);

#if DENSE_NORM > 0
	mpf_test_dsteep_norms (lp, pinf);
#endif

	mpf_ILLfct_update_dpI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.lindex,
	    alpha1);

END:
	if (singular) {
	    it->nextstep = SIMPLEX_RESUME;
	    it->resumeid = SIMPLEX_RESUME_SING;
	    /* this is to force to exit in the case of bad basis */
	    it->n_restart++;
	    /* fprintf(stderr,"Resume Singular
	       %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	    mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	    ILL_CLEANUP;
	}
	if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE) {
	    mpf_ILLfct_compute_piz (lp);
	    mpf_ILLfct_compute_dz (lp);
	    mpf_ILLfct_dual_adjust (lp, mpf_zeroLpNum);
	    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
	    mpf_ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);
	    if (fi.dstatus == DUAL_FEASIBLE)
		it->nextphase = DUAL_PHASEII;

	    it->newphase = SIMPLEX_PHASE_RECOMP;
	    ILL_CLEANUP;
	}
    }
#if DENSE_PI > 1
    mpf_fct_test_workvector (lp);
    mpf_fct_test_pI_x (lp, pinf);
#endif

CLEANUP:
    mpf_EGlpNumClearVar (alpha);
    mpf_EGlpNumClearVar (alpha1);
    mpf_EGlpNumClearVar (fi.totinfeas);
    mpf_EGlpNumClearVar (pr.dinfeas);
    mpf_EGlpNumClearVar (pr.pinfeas);
    mpf_EGlpNumClearVar (rs.tz);
    mpf_EGlpNumClearVar (rs.lbound);
    mpf_EGlpNumClearVar (rs.ecoeff);
    mpf_EGlpNumClearVar (rs.pivotval);
    /* ILL_RETURN (rval, "mpf_dual_phaseI_step"); */
    return rval;
}

static int mpf_dual_phaseII_step (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_svector * updz,
      mpf_svector * wz,
      mpf_iter_info * it)
{
    int coeffch;
    int rval = 0;
    int singular = 0;
    int refactor = 0;
    int ratio_iter = 0;
    int cphase = DUAL_PHASEII;
    int lcol, ecol;
    int estat, newphase;
    mpf_t x_bi, v_l, eval;
    mpf_t ecoeff;
    mpf_t alpha;
    mpf_t alpha1;
    mpf_feas_info fi;
    mpf_ratio_res rs;
    mpf_price_res pr;
    mpf_EGlpNumInitVar (x_bi);
    mpf_EGlpNumInitVar (v_l);
    mpf_EGlpNumInitVar (eval);
    mpf_EGlpNumInitVar (ecoeff);
    mpf_EGlpNumInitVar (alpha);
    mpf_EGlpNumInitVar (alpha1);
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumInitVar (pr.dinfeas);
    mpf_EGlpNumInitVar (pr.pinfeas);
    mpf_EGlpNumInitVar (rs.tz);
    mpf_EGlpNumInitVar (rs.lbound);
    mpf_EGlpNumInitVar (rs.ecoeff);
    mpf_EGlpNumInitVar (rs.pivotval);
    mpf_EGlpNumZero (rs.ecoeff);
    mpf_EGlpNumZero (alpha1);

    mpf_ILLfct_update_counts (lp, CNT_DPHASE2ITER, 0, mpf_zeroLpNum);
    it->nextstep = SIMPLEX_CONTINUE;
    it->nextphase = DUAL_PHASEII;
    lp->final_phase = DUAL_PHASEII;
    newphase = it->newphase;
    it->nosolve++;

    if (it->newphase != 0) {
	mpf_ILLfct_compute_dobj (lp);
	if (it->newphase == SIMPLEX_PHASE_NEW) {
	    it->noprog = 0;
	    if (it->sdisplay) {
		printf ("starting dual phase II\n");
		fflush (stdout);
	    }
	}
	it->newphase = 0;
	it->nosolve = 0;
	mpf_EGlpNumCopy (it->prevobj, lp->dobjval);
	mpf_ILLfct_compute_xbz (lp);

	if (pinf->d_strategy == COMPLETE_PRICING) {
#if USEHEAP > 0
	    mpf_ILLprice_free_heap (pinf);
#endif
	    mpf_ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEII);
#if USEHEAP > 0
	    rval = mpf_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
		DUAL_SIMPLEX, 0);
	    ILL_CLEANUP_IF (rval);
#endif
	} else if (pinf->d_strategy == MULTI_PART_PRICING) {
	    mpf_ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
	}
    }
    mpf_monitor_iter (lp, it, cphase);
    if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
	it->newphase != 0)
	ILL_CLEANUP;

    mpf_ILLprice_dual (lp, pinf, cphase, &pr);

    if (pr.price_stat == PRICE_OPTIMAL) {
	if (lp->ncchange != 0) {
	    if (it->sdisplay > 1) {
		printf ("unrolling %d coef shifts\n", lp->ncchange);
		fflush (stdout);
	    }
	    mpf_ILLfct_unroll_coef_change (lp);
	    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	    mpf_ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

	     /* HHH */ mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	    /* HHH* printf ("dual (opt) infeas %.6f\n", lp->dinfeas); fflush
	       (stdout); HHH* printf ("primal (opt) infeas %.6f\n",
	       lp->pinfeas); fflush (stdout); */

	    if (fi.dstatus != DUAL_FEASIBLE) {
		it->algorithm = PRIMAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_UNSHIFT;
		it->pricetype = QS_PRICE_PDEVEX;
		/* this is to force to exit in the case of bad basis */
		it->n_restart++;
		mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		/* fprintf(stderr,"Resume Unshift
		   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
		ILL_CLEANUP;
		/*
                 * it->nextphase = DUAL_PHASEI;
                 * lp->tol->ip_tol /= 5.0;
                 * lp->tol->id_tol /= 5.0;
                 * ILL_CLEANUP;
                 */
	    }
	}
	if (it->sdisplay > 1) {
	    /* ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    printf ("problem seemingly solved\n");
	    printf ("seemingly dual opt = %f\n", mpf_EGlpNumToLf (lp->dobjval));
	    printf ("retesting soln\n");
	    fflush (stdout);
	}
	rval = mpf_ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
	ILL_CLEANUP_IF (rval);
	mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

	if (fi.dstatus == DUAL_INFEASIBLE) {
	    ILL_IFTRACE ("DUAL_INFEAS: %s\n", __func__);
	    it->nextphase = DUAL_PHASEI;
	    mpf_EGlpNumDivUiTo (lp->tol->ip_tol, 5);
	    mpf_EGlpNumDivUiTo (lp->tol->id_tol, 5);
	} else if (fi.pstatus == PRIMAL_FEASIBLE) {
	    ILL_IFTRACE ("PRIM_FEAS: %s\n", __func__);
	    mpf_EGlpNumCopy (lp->objval, lp->dobjval);
	    it->solstatus = ILL_LP_SOLVED;
	    it->nextstep = SIMPLEX_TERMINATE;
	} else
	    ILL_IFTRACE ("PRIM_INFEAS: %s\n", __func__);
	ILL_CLEANUP;
    }
    mpf_ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
    /* ILL_IFTRACE("%s:d\n",__func__); */
    mpf_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
    mpf_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, mpf_zeroLpNum);
    mpf_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, mpf_zeroLpNum);

    ratio_iter = 0;
    do {
	mpf_ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
	if (rs.ratio_stat == RATIO_NEGATIVE) {
	    if (it->sdisplay > 1) {
		printf ("adjust coefs to remove negative ratio tests\n");
		fflush (stdout);
	    }
	    mpf_ILLfct_adjust_viol_coefs (lp);
	    mpf_ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
	    if (rs.ratio_stat == RATIO_NEGATIVE) {
		printf ("internal error: bad ratio test\n");
		fflush (stdout);
		rs.ratio_stat = RATIO_FAILED;
		break;
	    }
	}
	coeffch = rs.coeffch;
	mpf_EGlpNumCopy (ecoeff, rs.ecoeff);
	ratio_iter++;

	if (coeffch) {
	    /*
             * if (ratio_iter > PARAM_DRATIOTESTS){
             * ecoeff = lp->cz[lp->nbaz[rs.eindex]] - lp->dz[rs.eindex];
             * coeffch = 0;
             * }
             */
	    coeffch = 0;
	    rval = mpf_ILLfct_coef_shift (lp, lp->nbaz[rs.eindex], ecoeff);
	    ILL_CLEANUP_IF (rval);
	}
	if (rs.ratio_stat == RATIO_BCHANGE)
	    if (lp->vstat[lp->nbaz[rs.eindex]] == STAT_ZERO)
		break;

    } while (coeffch);

    if (rs.ratio_stat == RATIO_FAILED) {
	/*
         * rval = E_SIMPLEX_ERROR;
         * it->solstatus = ILL_DPHASEII_ERROR;
         */
	it->algorithm = PRIMAL_SIMPLEX;
	it->nextstep = SIMPLEX_RESUME;
	it->resumeid = SIMPLEX_RESUME_NUMER;
	/* this is to force to exit in the case of bad basis */
	it->n_restart++;
	mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	/* fprintf(stderr,"Resume Numerical
	   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	ILL_CLEANUP;
    } else if (rs.ratio_stat == RATIO_UNBOUNDED) {
	lp->infub_ix = pr.lindex;
	if (lp->ncchange != 0) {
	    if (it->sdisplay > 1) {
		printf ("unrolling %d coef shifts\n", lp->ncchange);
		fflush (stdout);
	    }
	    mpf_ILLfct_unroll_coef_change (lp);
	}
	mpf_ILLfct_set_status_values (lp, -1, DUAL_UNBOUNDED, -1, PHASEII);
	it->solstatus = ILL_LP_SOLVED;
	it->nextstep = SIMPLEX_TERMINATE;
    } else if (rs.ratio_stat == RATIO_BCHANGE) {
	lcol = lp->baz[pr.lindex];
	ecol = lp->nbaz[rs.eindex];

	mpf_ILLfct_compute_yz (lp, &(lp->yjz), updz, ecol);
	mpf_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, mpf_zeroLpNum);
	mpf_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, mpf_zeroLpNum);
	rval = mpf_ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
	if (rval != 0) {
	    it->n_pivot_fail++;
	    if (it->n_pivot_fail > SIMPLEX_MAX_PIVOT_FAIL) {
		it->n_pivot_fail = 0;
		/* this is to force to exit in the case of bad basis */
		it->algorithm = PRIMAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_NUMER;
		it->n_restart++;
		mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		/* fprintf(stderr,"Resume Pivot
		   %s:%s:%d\n",__func__,__FILE__,__LINE__); */
		rval = 0;
		ILL_CLEANUP;
	    }
	    if (newphase == 0) {
		rval = mpf_ILLbasis_factor (lp, &singular);
		ILL_CLEANUP_IF (rval);
		if (singular == 0)
		    refactor = 1;
		goto END;
	    } else {
		if (it->sdisplay > 1) {
		    printf ("warning: bad step\n");
		    fflush (stdout);
		}
	    }
	}
	mpf_EGlpNumAddTo (lp->dobjval, lp->upd.c_obj);
	mpf_EGlpNumCopy (lp->objval, lp->dobjval);

	if (!mpf_test_progress (lp->dobjval, it->prevobj)) {
	    if (lp->vtype[lcol] == VARTIFICIAL || lp->vtype[ecol] == VFREE) {
		if (it->noprog > 0)
		    it->noprog--;
	    } else
		it->noprog++;
	} else {
	    mpf_EGlpNumCopy (it->prevobj, lp->dobjval);
	    it->noprog = 0;
	}

	if (pinf->dII_price == QS_PRICE_DSTEEP)
	    mpf_ILLfct_compute_dsteep_upv (lp, wz);
	rval =
	    mpf_ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex,
	    rs.pivotval);
	ILL_CLEANUP_IF (rval);

	mpf_EGlpNumCopy (x_bi, lp->xbz[pr.lindex]);
	if (pr.lvstat == STAT_LOWER)
	    mpf_EGlpNumCopy (v_l, lp->lz[lcol]);
	else
	    mpf_EGlpNumCopy (v_l, lp->uz[lcol]);
	mpf_EGlpNumCopy (alpha, rs.tz);
	if (pr.lvstat == STAT_LOWER)
	    mpf_EGlpNumSign (alpha);
	estat = lp->vstat[ecol];
	if (estat == STAT_LOWER)
	    mpf_EGlpNumCopy (eval, lp->lz[ecol]);
	else if (estat == STAT_ZERO)
	    mpf_EGlpNumZero (eval);
	else
	    mpf_EGlpNumCopy (eval, lp->uz[ecol]);

	mpf_ILLfct_update_piz (lp, alpha);
	mpf_ILLfct_update_dz (lp, rs.eindex, alpha);
	mpf_ILLfct_update_dIIfeas (lp, rs.eindex, &(lp->srhs));
	mpf_ILLfct_compute_dpIIy (lp, &(lp->srhs), &(lp->ssoln));
	mpf_EGlpNumCopyDiff (alpha1, x_bi, v_l);
	mpf_EGlpNumSubTo (alpha1, lp->upd.dty);
	mpf_EGlpNumDivTo (alpha1, rs.pivotval);
	mpf_ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);
	rval = mpf_ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
	ILL_CLEANUP_IF (rval);

	mpf_ILLfct_update_dpII_prices (lp, pinf, &(lp->srhs), &(lp->ssoln),
	    pr.lindex, eval, alpha1);

#if DENSE_NORM > 0
	mpf_test_dsteep_norms (lp, pinf);
#endif

END:
	if (singular) {
	    it->nextstep = SIMPLEX_RESUME;
	    it->resumeid = SIMPLEX_RESUME_SING;
	    /* this is to force to exit in the case of bad basis */
	    it->n_restart++;
	    mpf_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
	    mpf_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
	    /* fprintf(stderr,"Resume Singular
	       %s:%s:%d\n",__func__,__FILE__,__LINE__); */
	    ILL_CLEANUP;
	}
	if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE) {
	    mpf_ILLfct_compute_piz (lp);
	    mpf_ILLfct_compute_dz (lp);
	    mpf_ILLfct_dual_adjust (lp, mpf_zeroLpNum);
	    it->newphase = SIMPLEX_PHASE_RECOMP;
	}
    }
#if DENSE_PIIPI > 0
    mpf_fct_test_workvector (lp);
    if (!refactor) {
	mpf_fct_test_pII_x (lp, pinf);
	mpf_fct_test_pII_pi_dz (lp, pinf);
    }
#endif

CLEANUP:
    mpf_EGlpNumClearVar (x_bi);
    mpf_EGlpNumClearVar (v_l);
    mpf_EGlpNumClearVar (eval);
    mpf_EGlpNumClearVar (ecoeff);
    mpf_EGlpNumClearVar (alpha);
    mpf_EGlpNumClearVar (alpha1);
    mpf_EGlpNumClearVar (fi.totinfeas);
    mpf_EGlpNumClearVar (pr.dinfeas);
    mpf_EGlpNumClearVar (pr.pinfeas);
    mpf_EGlpNumClearVar (rs.tz);
    mpf_EGlpNumClearVar (rs.lbound);
    mpf_EGlpNumClearVar (rs.ecoeff);
    mpf_EGlpNumClearVar (rs.pivotval);
    /* ILL_RETURN (rval, "mpf_dual_phaseII_step"); */
    return rval;
}

static void mpf_get_current_stat (mpf_lp_status_info * p,
      int algorithm,
      int *bstat)
{
    if (p->optimal)
	*bstat = OPTIMAL;
    else if (algorithm == PRIMAL_SIMPLEX) {
	if (p->primal_feasible)
	    *bstat = PRIMAL_FEASIBLE;
	else if (p->primal_infeasible)
	    *bstat = PRIMAL_INFEASIBLE;
	else if (p->primal_unbounded)
	    *bstat = PRIMAL_UNBOUNDED;
	else
	    *bstat = NONOPTIMAL;
    } else if (algorithm == DUAL_SIMPLEX) {
	if (p->dual_feasible)
	    *bstat = DUAL_FEASIBLE;
	else if (p->dual_infeasible)
	    *bstat = DUAL_INFEASIBLE;
	else if (p->dual_unbounded)
	    *bstat = DUAL_UNBOUNDED;
	else
	    *bstat = NONOPTIMAL;
    }
}

int mpf_ILLsimplex_pivotin (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      int rcnt,
      int *rlist,
      int pivot_opt,
      int *basis_mod)
{
    int i, npiv = 0;
    int eindex;
    int rval = 0;
    int singular = 0;
    int refactor = 0;
    int *rowmap = lp->O->rowmap;
    int *clist = NULL;
    mpf_svector wz;
    mpf_svector updz;
    mpf_t alpha;
    mpf_ratio_res rs;
    mpf_feas_info fi;
    mpf_EGlpNumInitVar (alpha);
    mpf_EGlpNumInitVar (fi.totinfeas);
    mpf_EGlpNumInitVar (rs.tz);
    mpf_EGlpNumInitVar (rs.lbound);
    mpf_EGlpNumInitVar (rs.ecoeff);
    mpf_EGlpNumInitVar (rs.pivotval);
    mpf_EGlpNumZero (alpha);

    *basis_mod = 0;
    if (rcnt <= 0) {
	ILL_RETURN (rval, "mpf_ILLsimplex_pivotin");
    }
    if (pivot_opt == SIMPLEX_PIVOTINROW) {
	ILL_SAFE_MALLOC (clist, rcnt, int);
	for (i = 0; i < rcnt; i++)
	    clist[i] = rowmap[rlist[i]];
    } else
	clist = rlist;

    for (i = 0; i < rcnt; i++) {
	if (lp->vstat[clist[i]] != STAT_BASIC) {
	    *basis_mod = 1;
	    break;
	}
    }
    if (*basis_mod == 0) {
	if (pivot_opt == SIMPLEX_PIVOTINROW) {
	    ILL_IFFREE (clist, int);
	}
	ILL_RETURN (rval, "mpf_ILLsimplex_pivotin");
    }
    /* printf ("Forcing vars into basis in mpf_ILLsimplex_pivotin \n"); */
    mpf_ILLsvector_init (&wz);
    rval = mpf_ILLsvector_alloc (&wz, lp->nrows);
    ILL_CLEANUP_IF (rval);
    mpf_ILLsvector_init (&updz);
    rval = mpf_ILLsvector_alloc (&updz, lp->nrows);
    ILL_CLEANUP_IF (rval);

    mpf_EGlpNumCopy (lp->pobjval, lp->dobjval);
    for (i = 0; i < rcnt; i++) {
	if (lp->vstat[clist[i]] == STAT_BASIC)
	    continue;
	npiv++;

	eindex = lp->vindex[clist[i]];
	mpf_ILLfct_compute_yz (lp, &(lp->yjz), &updz, lp->nbaz[eindex]);
	mpf_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, mpf_zeroLpNum);
	mpf_ILLfct_update_counts (lp, CNT_UPNZ, updz.nzcnt, mpf_zeroLpNum);

	mpf_ILLratio_pivotin_test (lp, clist, rcnt, &rs);

	if (rs.ratio_stat == RATIO_UNBOUNDED || rs.ratio_stat == RATIO_FAILED) {
	    fprintf (stderr, "Pivot_in failed\n");
	    rval = E_SIMPLEX_ERROR;
	    ILL_CLEANUP;
	} else if (rs.ratio_stat == RATIO_BCHANGE) {
	    if (rs.lvstat == STAT_LOWER) {
		mpf_EGlpNumCopyDiff (alpha, lp->lz[lp->baz[rs.lindex]], lp->xbz[rs.lindex]);
		mpf_EGlpNumAddInnProdTo (lp->dobjval, rs.tz, alpha);
	    } else {
		mpf_EGlpNumCopyDiff (alpha, lp->xbz[rs.lindex], lp->uz[lp->baz[rs.lindex]]);
		mpf_EGlpNumAddInnProdTo (lp->dobjval, rs.tz, alpha);
	    }
	    mpf_EGlpNumCopyFrac (alpha, lp->dz[eindex], rs.pivotval);
	    mpf_EGlpNumCopy (lp->objval, lp->dobjval);

	    mpf_ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
	    /* ILL_IFTRACE("%s:e\n",__func__); */
	    mpf_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	    mpf_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, mpf_zeroLpNum);
	    mpf_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, mpf_zeroLpNum);

	    if (pinf->dsinfo.norms && pinf->dII_price == QS_PRICE_DSTEEP) {
		mpf_ILLfct_compute_dsteep_upv (lp, &wz);
		rval = mpf_ILLprice_update_pricing_info (lp, pinf, DUAL_PHASEII, &wz,
		    eindex, rs.lindex, rs.pivotval);
		ILL_CLEANUP_IF (rval);
	    } else if (pinf->psinfo.norms && pinf->pII_price == QS_PRICE_PSTEEP) {
		mpf_ILLfct_compute_psteep_upv (lp, &wz);
		rval = mpf_ILLprice_update_pricing_info (lp, pinf, PRIMAL_PHASEII, &wz,
		    eindex, rs.lindex, rs.pivotval);
		ILL_CLEANUP_IF (rval);
	    }
	    /* ILL_IFTRACE("%s:e\n",__func__); */
	    mpf_ILLfct_update_xz (lp, rs.tz, eindex, rs.lindex);
	    mpf_ILLfct_update_basis_info (lp, eindex, rs.lindex, rs.lvstat);
	    rval = mpf_ILLbasis_update (lp, &updz, rs.lindex, &refactor, &singular);
	    ILL_CLEANUP_IF (rval);

	    if (singular) {
		fprintf (stderr, "singular matrix in pivot_in\n");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	    }
	    if (!refactor) {
		mpf_ILLfct_update_piz (lp, alpha);
		mpf_ILLfct_update_dz (lp, eindex, alpha);
	    } else {
		mpf_ILLfct_compute_xbz (lp);
		mpf_ILLfct_compute_piz (lp);
		mpf_ILLfct_compute_dz (lp);
		mpf_ILLfct_compute_dobj (lp);
	    }
	}
    }
    /*
     * mpf_ILLfct_dphaseI_simple_update (lp, lp->tol->dfeas_tol);
     * mpf_ILLfct_compute_xbz (lp);
     * mpf_ILLfct_compute_dobj (lp);
     */

    mpf_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
    mpf_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
    mpf_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

CLEANUP:
    if (pivot_opt == SIMPLEX_PIVOTINROW)
	ILL_IFFREE (clist, int);
    mpf_ILLsvector_free (&wz);
    mpf_ILLsvector_free (&updz);
    mpf_EGlpNumClearVar (alpha);
    mpf_EGlpNumClearVar (fi.totinfeas);
    mpf_EGlpNumClearVar (rs.tz);
    mpf_EGlpNumClearVar (rs.lbound);
    mpf_EGlpNumClearVar (rs.ecoeff);
    mpf_EGlpNumClearVar (rs.pivotval);
    ILL_RETURN (rval, "mpf_ILLsimplex_pivotin");
}

static int mpf_report_value (mpf_lpinfo * lp,
      mpf_iter_info * it,
      const char *value_name,
      mpf_t value)
{
    int rval = 0;

    if (it->sdisplay && it->itercnt % lp->iterskip == 0) {
	char buffer[1024], *nstr;
	nstr = mpf_EGlpNumGetStr (value);
	snprintf (buffer, (size_t) 1023, "(%d): %s = %10.7lf %s\n", it->itercnt,
	    value_name, mpf_EGlpNumToLf (value), nstr);
	free (nstr);
	buffer[1022] = '\n';
	buffer[1023] = '\0';
	rval = ILLstring_report (buffer, &lp->O->reporter);
	fflush (stdout);
    } else {
	/* make sure ILLstring_report is called at least every 10 iterations */
	if (it->itercnt % (lp->iterskip / 10)) {
	    rval = ILLstring_report (NULL, &lp->O->reporter);
	}
    }
    if (rval != 0) {		/* ILLstring_report was called and failed,
				   which means we should abort */
	it->solstatus = QS_LP_ABORTED;
    }
    return rval;
}
