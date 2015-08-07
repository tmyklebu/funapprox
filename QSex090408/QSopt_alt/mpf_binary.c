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

/* RCS_INFO = "$RCSfile: mpf_binary.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* Simple MIP Code to test LP Solver                    */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int mpf_ILLmip_bfs (mpf_lpinfo *lp, double *val, double *x)                   */
/* */
/* NOTES                                                                   */
/* */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpf_priority.h"
#include "mpf_sortrus.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdata.h"
#include "mpf_lpdefs.h"
#include "mpf_simplex.h"
#include "mpf_binary.h"
#include "mpf_price.h"
#include "mpf_lib.h"
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* #define  mpf_ILL_INTTOL (0.000001) */
#define mpf_ILL_INTTOL mpf_PFEAS_TOLER

#define  mpf_STRONG_PIVOTS     (50)
#define  mpf_STRONG_CANDIDATES (10)

#define mpf_ILL_BRANCH_STRONG_WEIGHT (10)
#define mpf_ILL_BRANCH_STRONG_VAL(v0,v1)                                 \
    (((v0) < (v1) ? (mpf_ILL_BRANCH_STRONG_WEIGHT * (v0) + (v1))         \
                  : (mpf_ILL_BRANCH_STRONG_WEIGHT * (v1) + (v0)))        \
                    / (mpf_ILL_BRANCH_STRONG_WEIGHT + 1.0))

#define mpf_ILL_BRANCH_PENALTY_WEIGHT (2)
#define mpf_ILL_BRANCH_PENALTY_VAL(v0,v1,f)                              \
    (((v0)*(f) < (v1)*(1.0-(f)) ?                                    \
        (mpf_ILL_BRANCH_PENALTY_WEIGHT * (v0)*(f) + (v1)*(1.0-(f)))    \
      : (mpf_ILL_BRANCH_PENALTY_WEIGHT * (v1)*(1.0-(f)) + (v0)*(f)))    \
                    / (mpf_ILL_BRANCH_PENALTY_WEIGHT + 1.0))



#define mpf_FIRSTBRANCH  1
#define mpf_MIDDLEBRANCH 2
#define mpf_STRONGBRANCH 3
#define mpf_PENALTYBRANCH 4


typedef struct mpf_bbnode {
    struct mpf_bbnode *next;
    struct mpf_bbnode *prev;
    int id;
    int depth;
    int handle;
    mpf_t bound;
    char *cstat;
    char *rstat;
    mpf_t *rownorms;
    int rownorms_size;
    int bound_cnt;
    int *bound_indx;
    char *lu;
    mpf_t *bounds;
    int bounds_size;
}
  mpf_bbnode;

typedef struct mpf_mipinfo {
    int branching_rule;
    int watch;
    int depth;
    int totalnodes;
    int activenodes;
    int totalpivots;
    int lastpivots;
    int objsense;
    mpf_t objectivebound;
    mpf_t value;
    mpf_t *downpen;
    mpf_t *uppen;
    mpf_t *x;
    mpf_t *bestx;
    mpf_t *orig_lower;
    mpf_t *orig_upper;
    mpf_t *lower;
    mpf_t *upper;
    int nstruct;		/* size of all mpf_t arrays */
    mpf_lpinfo *lp;
    mpf_price_info *pinf;
    mpf_bbnode head_bbnode;
    mpf_ILLpriority *que;
    ILLptrworld ptrworld;
}
  mpf_mipinfo;


ILL_PTRWORLD_ROUTINES (mpf_bbnode, bbnodealloc, bbnode_bulkalloc, bbnodefree)
ILL_PTRWORLD_LISTFREE_ROUTINE (mpf_bbnode, bbnode_listfree, bbnodefree)
ILL_PTRWORLD_LEAKS_ROUTINE (mpf_bbnode, bbnode_check_leaks, depth, int)
static void mpf_cleanup_mip (mpf_mipinfo * minf),
  mpf_choose_initial_price (mpf_price_info * pinf),
  mpf_best_bbnode (mpf_mipinfo * minf,
      mpf_bbnode ** best),
  mpf_put_bbnode (mpf_mipinfo * minf,
      mpf_bbnode * b),
  mpf_remove_bbnode (mpf_bbnode * b),
  mpf_find_first_branch (mpf_lpinfo * lp,
      mpf_t * x,
      int *bvar),
  mpf_find_middle_branch (mpf_lpinfo * lp,
      mpf_t * x,
      int *bvar),
  mpf_check_integral (mpf_lpinfo * lp,
      mpf_t * x,
      int *yesno),
  mpf_copy_x (int nstruct,
      mpf_t * from_x,
      mpf_t * to_x),
  mpf_init_mipinfo (mpf_mipinfo * minf),
  mpf_free_mipinfo (mpf_mipinfo * minf),
  mpf_init_bbnode (mpf_bbnode * b),
  mpf_free_bbnode (mpf_bbnode * b);

static int mpf_startup_mip (mpf_mipinfo * minf,
      mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * lpval),
  mpf_run_bfs (mpf_mipinfo * minf),
  mpf_process_bfs_bbnode (mpf_mipinfo * minf,
      mpf_bbnode * b),
  mpf_child_work (mpf_mipinfo * minf,
      mpf_bbnode * active,
      int bvar,
      int bdir,
      mpf_t * cval,
      int *cp),
  mpf_fix_variables (mpf_lpinfo * lp,
      mpf_t * bestval,
      mpf_bbnode * b,
      mpf_t * wupper,
      mpf_t * wlower,
      int *hit),
  mpf_find_branch (mpf_mipinfo * minf,
      mpf_t * x,
      mpf_t * lpval,
      int *bvar),
  mpf_find_penalty_branch (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * x,
      mpf_t * downpen,
      mpf_t * uppen,
      mpf_t * lpval,
      int *bvar),
  mpf_find_strong_branch (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * x,
      int *bvar),
  mpf_plunge (mpf_mipinfo * minf),
  mpf_plunge_work (mpf_mipinfo * minf,
      int depth),
  mpf_round_variables (mpf_mipinfo * minf,
      int *count,
      mpf_t * tol);


static void mpf_choose_initial_price (mpf_price_info * pinf)
{
    pinf->pI_price = QS_PRICE_PSTEEP;
    pinf->pII_price = QS_PRICE_PSTEEP;
    pinf->dI_price = QS_PRICE_DSTEEP;
    pinf->dII_price = QS_PRICE_DSTEEP;
}

int mpf_ILLmip_bfs (mpf_lpinfo * lp,
      mpf_t * val,
      mpf_t * x)
{
    int tval, rval = 0;
    mpf_price_info pinf;
    mpf_mipinfo minf;
    mpf_bbnode *b;
    mpf_t lpval;
    double szeit = ILLutil_zeit ();
    mpf_EGlpNumInitVar (lpval);
    mpf_EGlpNumInitVar (pinf.htrigger);

    mpf_ILLprice_init_pricing_info (&pinf);
    mpf_init_mipinfo (&minf);

    if (!lp) {
	fprintf (stderr, "mpf_ILLmip_bfs called without an LP\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpf_startup_mip (&minf, lp, &pinf, &lpval);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (minf.que, 1, mpf_ILLpriority);
    rval = mpf_ILLutil_priority_init (minf.que, lp->O->nstruct + 1);
    ILL_CLEANUP_IF (rval);

    b = bbnodealloc (&minf.ptrworld);
    mpf_init_bbnode (b);
    b->depth = 0;
    b->id = minf.totalnodes++;
    mpf_EGlpNumCopy (b->bound, lpval);
    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
    rval = mpf_ILLlib_getbasis (lp, b->cstat, b->rstat);
    ILL_CLEANUP_IF (rval);

    if (pinf.dII_price == QS_PRICE_DSTEEP) {
	b->rownorms = mpf_EGlpNumAllocArray (lp->nrows);
	tval = mpf_ILLlib_getrownorms (lp, &pinf, b->rownorms);
	if (tval) {
	    printf ("Row norms not available\n");
	    fflush (stdout);
	    mpf_EGlpNumFreeArray (b->rownorms);
	}
    }
    rval = mpf_ILLutil_priority_insert (minf.que, (void *) b, &lpval, &(b->handle));
    ILL_CLEANUP_IF (rval);

    b->prev = &(minf.head_bbnode);
    b->next = 0;
    minf.head_bbnode.next = b;
    minf.activenodes++;

    minf.branching_rule = mpf_PENALTYBRANCH;

    rval = mpf_run_bfs (&minf);
    ILL_CLEANUP_IF (rval);

    printf ("Total Number of Nodes: %d\n", minf.totalnodes);
    printf ("Total Number of Pivots: %d\n", minf.totalpivots);
    printf ("BFS MIP Runing Time: %.2f seconds\n", ILLutil_zeit () - szeit);
    fflush (stdout);

    mpf_EGlpNumCopy (*val, minf.value);
    if (minf.objsense == mpf_ILL_MAX)
	mpf_EGlpNumSign (*val);

    if (x && mpf_EGlpNumIsNeqq (minf.value, mpf_ILL_MAXDOUBLE)) {
	mpf_copy_x (lp->O->nstruct, minf.bestx, x);
    }
CLEANUP:

    if (minf.que) {
	mpf_ILLutil_priority_free (minf.que);
	ILL_IFFREE (minf.que, mpf_ILLpriority);
    }
    mpf_cleanup_mip (&minf);
    mpf_free_mipinfo (&minf);
    mpf_ILLprice_free_pricing_info (&pinf);
    mpf_EGlpNumClearVar (lpval);
    mpf_EGlpNumClearVar (pinf.htrigger);
    ILL_RETURN (rval, "mpf_ILLmip_bfs");
}

static int mpf_startup_mip (mpf_mipinfo * minf,
      mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * lpval)
{
    int rval = 0;
    int i, col, status, intcount = 0;
    mpf_t val;
    mpf_ILLlpdata *qlp;
    mpf_EGlpNumInitVar (val);

    mpf_choose_initial_price (pinf);

    qlp = lp->O;

    rval = mpf_ILLlib_optimize (lp, 0, pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += mpf_ILLlib_iter (lp);

    rval = mpf_ILLlib_objval (lp, 0, &val);
    ILL_CLEANUP_IF (rval);

    printf ("LP Value: %.6f\n", mpf_EGlpNumToLf (val));
    fflush (stdout);
    if (lpval)
	mpf_EGlpNumCopy (*lpval, val);

    if (qlp->intmarker) {
	for (i = 0; i < qlp->nstruct; i++) {
	    if (qlp->intmarker[i]) {
		col = qlp->structmap[i];
		intcount++;
		if (mpf_EGlpNumIsEqqual (qlp->lower[col], mpf_ILL_MINDOUBLE)
		    || mpf_EGlpNumIsEqqual (qlp->upper[col], mpf_ILL_MAXDOUBLE)) {
		    printf ("Instance has unbounded integer variable\n");
		    fflush (stdout);
		    rval = 1;
		    goto CLEANUP;
		}
	    }
	}
    }
    if (intcount == 0) {
	printf ("No integer variables\n");
	fflush (stdout);
	rval = 1;
	goto CLEANUP;
    } else {
	printf ("%d integer variables\n", intcount);
	fflush (stdout);
    }

    if (qlp->sinfo) {		/* Free the presolve LP and work with orginal */
	mpf_ILLlp_sinfo_free (qlp->sinfo);
	ILL_IFFREE (qlp->sinfo, mpf_ILLlp_sinfo);
    }
    minf->lp = lp;
    minf->pinf = pinf;
    minf->objsense = qlp->objsense;
    if (qlp->objsense == mpf_ILL_MAX) {	/* MIP codes work with min */
	for (i = 0; i < lp->ncols; i++) {
	    mpf_EGlpNumCopyNeg (qlp->obj[i], qlp->obj[i]);
	}
	qlp->objsense = mpf_ILL_MIN;
    }
    minf->x = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->bestx = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->lower = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->upper = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->orig_lower = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->orig_upper = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->downpen = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->uppen = mpf_EGlpNumAllocArray (qlp->nstruct);
    minf->nstruct = qlp->nstruct;

    rval = mpf_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
	mpf_EGlpNumCopy (minf->lower[i], qlp->lower[i]);
	mpf_EGlpNumCopy (minf->upper[i], qlp->upper[i]);
	mpf_EGlpNumCopy (minf->orig_lower[i], qlp->lower[i]);
	mpf_EGlpNumCopy (minf->orig_upper[i], qlp->upper[i]);
	mpf_EGlpNumOne (minf->downpen[i]);
	mpf_EGlpNumOne (minf->uppen[i]);
	mpf_EGlpNumSign (minf->downpen[i]);
	mpf_EGlpNumSign (minf->uppen[i]);
    }


CLEANUP:

    mpf_EGlpNumClearVar (val);
    ILL_RETURN (rval, "mpf_startup_mip");
}

static void mpf_cleanup_mip (mpf_mipinfo * minf)
{
    int i;
    mpf_ILLlpdata *qslp;

    if (minf && minf->lp) {
	qslp = minf->lp->O;
	if (minf->objsense == mpf_ILL_MAX) {
	    for (i = 0; i < minf->lp->ncols; i++) {
		mpf_EGlpNumSign (qslp->obj[i]);
	    }
	    qslp->objsense = mpf_ILL_MIN;
	}
    }
}

static int mpf_run_bfs (mpf_mipinfo * minf)
{
    int rval = 0;
    mpf_bbnode *b;

    while (minf->head_bbnode.next) {
	mpf_best_bbnode (minf, &b);
	rval = mpf_process_bfs_bbnode (minf, b);
	ILL_CLEANUP_IF (rval);
	mpf_remove_bbnode (b);
	mpf_free_bbnode (b);
	bbnodefree (&minf->ptrworld, b);
	minf->activenodes--;
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_run_bfs");
}

static int mpf_process_bfs_bbnode (mpf_mipinfo * minf,
      mpf_bbnode * active)
{
    mpf_lpinfo *lp = minf->lp;
    mpf_ILLlp_basis B;
    int status, bvar = 0;
    int i, j, hit, dnp = 0, upp = 0;
    int nstruct = lp->O->nstruct;
    mpf_t t, lpval, dnval, upval;
    mpf_t *wupper = 0;
    mpf_t *wlower = 0;
    int rval = 0;
    mpf_EGlpNumInitVar (t);
    mpf_EGlpNumInitVar (lpval);
    mpf_EGlpNumInitVar (dnval);
    mpf_EGlpNumInitVar (upval);

    mpf_ILLlp_basis_init (&B);

    if (minf->watch > 1) {
	printf ("Node %4d: %.3f", active->id, mpf_EGlpNumToLf (active->bound));
	if (mpf_EGlpNumIsNeqq (minf->value, mpf_ILL_MAXDOUBLE))
	    printf (" %.3f", mpf_EGlpNumToLf (minf->value));
	else
	    printf ("  None");
	printf (", Active %d ", minf->activenodes);
	fflush (stdout);
    } else if (minf->watch == 1) {
	if (minf->lastpivots > 1000) {
	    minf->lastpivots = 0;
	    printf ("Pivots %d, Active Nodes %d, Bound %.3f, Soln ",
		minf->totalpivots, minf->activenodes,
		mpf_EGlpNumToLf (active->bound));
	    if (!mpf_EGlpNumIsLess (minf->value, mpf_ILL_MAXDOUBLE))
		printf ("%.3f", mpf_EGlpNumToLf (minf->value));
	    else
		printf ("None\n");
	}
    }
    if (mpf_EGlpNumIsLeq (minf->objectivebound, active->bound)) {
	if (minf->watch > 1) {
	    printf ("  Node can be purged\n");
	    fflush (stdout);
	}
	goto CLEANUP;
    }
    /* Set the LP bounds for the mpf_node. */

    wlower = mpf_EGlpNumAllocArray (nstruct);
    wupper = mpf_EGlpNumAllocArray (nstruct);

    for (i = 0; i < nstruct; i++) {
	mpf_EGlpNumCopy (wlower[i], minf->orig_lower[i]);
	mpf_EGlpNumCopy (wupper[i], minf->orig_upper[i]);
    }
    for (i = 0; i < active->bound_cnt; i++) {
	j = active->bound_indx[i];
	if (active->lu[i] == 'L')
	    mpf_EGlpNumCopy (wlower[j], active->bounds[i]);
	else
	    mpf_EGlpNumCopy (wupper[j], active->bounds[i]);
    }

    if (active->bound_cnt > 0) {
	rval = mpf_ILLlib_chgbnds (lp, active->bound_cnt, active->bound_indx,
	    active->lu, active->bounds);
	ILL_CLEANUP_IF (rval);
    }
    /* Solve the LP. */

    rval = mpf_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
	active->rstat);
    ILL_CLEANUP_IF (rval);
    if (active->rownorms) {
	B.rownorms = mpf_EGlpNumAllocArray (lp->nrows);
	for (i = 0; i < lp->nrows; i++) {
	    mpf_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
	}
    }
    rval = mpf_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += mpf_ILLlib_iter (lp);
    minf->lastpivots += mpf_ILLlib_iter (lp);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    }
    if (status == QS_LP_INFEASIBLE) {
	printf ("  Infeasible LP, should have been purged earlier\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    }
    if (active->depth < 0) {
	for (i = 0; i < nstruct; i++) {
	    mpf_EGlpNumCopy (minf->lower[i], wlower[i]);
	    mpf_EGlpNumCopy (minf->upper[i], wupper[i]);
	}
	rval = mpf_plunge (minf);
	ILL_CLEANUP_IF (rval);
    }
    /* Fix variables. */

    if (mpf_EGlpNumIsLess (minf->value, mpf_ILL_MAXDOUBLE)) {
	rval = mpf_fix_variables (lp, &(minf->value), active, wupper, wlower, &hit);
	ILL_CLEANUP_IF (rval);

	if (hit) {
	    rval = mpf_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
	    ILL_CLEANUP_IF (rval);

	    minf->totalpivots += mpf_ILLlib_iter (lp);
	    minf->lastpivots += mpf_ILLlib_iter (lp);

	    if (status == QS_LP_UNSOLVED) {
		printf ("Simplex did not solve the LP\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	    }
	    if (status == QS_LP_INFEASIBLE) {
		printf ("  Infeasible LP after fixing\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	    }
	}
    }
    /* Branch. */

    rval = mpf_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = mpf_ILLlib_objval (lp, 0, &lpval);
    ILL_CLEANUP_IF (rval);

    rval = mpf_find_branch (minf, minf->x, &lpval, &bvar);
    ILL_CLEANUP_IF (rval);

    if (bvar == -1) {
	printf ("Found integral solution: %f\n", mpf_EGlpNumToLf (lpval));
	if (mpf_EGlpNumIsLess (lpval, minf->value)) {
	    mpf_EGlpNumCopy (minf->value, lpval);
	    mpf_EGlpNumCopy (minf->objectivebound, lpval);
	    mpf_EGlpNumSubTo (minf->objectivebound, mpf_ILL_INTTOL);
	    mpf_copy_x (nstruct, minf->x, minf->bestx);
	}
    } else {
	/* Create down child */

	rval = mpf_child_work (minf, active, bvar, 'D', &dnval, &dnp);
	ILL_CLEANUP_IF (rval);

	/* Restore parent basis */

	rval = mpf_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
	    active->rstat);
	ILL_CLEANUP_IF (rval);
	if (active->rownorms) {
	    B.rownorms = mpf_EGlpNumAllocArray (lp->nrows);
	    for (i = 0; i < lp->nrows; i++) {
		mpf_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
	    }
	}
	/* Create up child */

	rval = mpf_child_work (minf, active, bvar, 'U', &upval, &upp);
	ILL_CLEANUP_IF (rval);

	if (minf->watch > 1) {
	    if (mpf_EGlpNumIsEqqual (dnval, mpf_ILL_MAXDOUBLE)) {
		printf ("DN->XXX");
	    } else {
		printf ("DN->%.3f%c", mpf_EGlpNumToLf (dnval), dnp ? 'X' : ' ');
	    }
	    if (mpf_EGlpNumIsEqqual (upval, mpf_ILL_MAXDOUBLE)) {
		printf ("UP->XXX\n");
	    } else {
		printf ("UP->%.3f%c\n", mpf_EGlpNumToLf (upval), upp ? 'X' : ' ');
	    }
	    fflush (stdout);
	}
    }

    /* Set the LP bounds back to original values */

    for (i = 0; i < active->bound_cnt; i++) {
	if (active->lu[i] == 'L')
	    mpf_EGlpNumCopy (t, minf->orig_lower[active->bound_indx[i]]);
	else
	    mpf_EGlpNumCopy (t, minf->orig_upper[active->bound_indx[i]]);

	rval = mpf_ILLlib_chgbnd (lp, active->bound_indx[i], active->lu[i], t);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    mpf_EGlpNumFreeArray (wlower);
    mpf_EGlpNumFreeArray (wupper);
    mpf_ILLlp_basis_free (&B);
    mpf_EGlpNumClearVar (t);
    mpf_EGlpNumClearVar (lpval);
    mpf_EGlpNumClearVar (dnval);
    mpf_EGlpNumClearVar (upval);
    ILL_RETURN (rval, "mpf_process_bfs_bbnode");
}

static int mpf_child_work (mpf_mipinfo * minf,
      mpf_bbnode * active,
      int bvar,
      int bdir,
      mpf_t * cval,
      int *cp)
{
    int tval, rval = 0;
    int i, status, intsol;
    mpf_t t, oldt, lpval;
    mpf_t *xi = &(minf->x[bvar]);
    mpf_lpinfo *lp = minf->lp;
    mpf_bbnode *b;
    mpf_EGlpNumInitVar (t);
    mpf_EGlpNumInitVar (lpval);
    mpf_EGlpNumInitVar (oldt);

    *cp = 0;

    if (bdir == 'D') {
	rval = mpf_ILLlib_getbnd (lp, bvar, 'U', &oldt);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumFloor (t, *xi);
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'U', t);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = mpf_ILLlib_getbnd (lp, bvar, 'L', &oldt);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumCeil (t, *xi);
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'L', t);
	ILL_CLEANUP_IF (rval);
    }

    rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += mpf_ILLlib_iter (lp);
    minf->lastpivots += mpf_ILLlib_iter (lp);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve Child LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    }
    if (status == QS_LP_INFEASIBLE) {
	mpf_EGlpNumCopy (*cval, mpf_ILL_MAXDOUBLE);
	*cp = 1;
    } else {
	rval = mpf_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumCopy (*cval, lpval);

	/* What about the x vector?  Bico - 020531 */

	mpf_check_integral (lp, minf->x, &intsol);
	if (intsol) {
	    if (mpf_EGlpNumIsLess (lpval, minf->value)) {
		printf ("Found integral solution: %f\n", mpf_EGlpNumToLf (lpval));
		mpf_EGlpNumCopy (minf->value, lpval);
		mpf_EGlpNumCopy (minf->objectivebound, lpval);
		mpf_EGlpNumSubTo (minf->objectivebound, mpf_ILL_INTTOL);
		mpf_copy_x (lp->O->nstruct, minf->x, minf->bestx);
	    }
	}
	if (mpf_EGlpNumIsLeq (minf->objectivebound, lpval)) {
	    *cp = 1;
	} else {
	    b = bbnodealloc (&minf->ptrworld);
	    mpf_init_bbnode (b);
	    b->depth = active->depth + 1;
	    b->id = minf->totalnodes;
	    mpf_EGlpNumCopy (b->bound, lpval);
	    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
	    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
	    rval = mpf_ILLlib_getbasis (lp, b->cstat, b->rstat);
	    ILL_CLEANUP_IF (rval);
	    if (minf->pinf->dII_price == QS_PRICE_DSTEEP) {
		b->rownorms = mpf_EGlpNumAllocArray (lp->nrows);
		tval = mpf_ILLlib_getrownorms (lp, minf->pinf, b->rownorms);
		if (tval) {
		    printf ("Row norms not available\n");
		    fflush (stdout);
		    printf ("A\n");
		    exit (1);
		    mpf_EGlpNumFreeArray (b->rownorms);
		}
	    }
	    ILL_SAFE_MALLOC (b->bound_indx, active->bound_cnt + 1, int);
	    ILL_SAFE_MALLOC (b->lu, active->bound_cnt + 1, char);
	    b->bounds = mpf_EGlpNumAllocArray (active->bound_cnt + 1);
	    for (i = 0; i < active->bound_cnt; i++) {
		b->bound_indx[i] = active->bound_indx[i];
		b->lu[i] = active->lu[i];
		mpf_EGlpNumCopy (b->bounds[i], active->bounds[i]);
	    }
	    b->bound_indx[active->bound_cnt] = bvar;
	    if (bdir == 'D')
		b->lu[active->bound_cnt] = 'U';
	    else
		b->lu[active->bound_cnt] = 'L';
	    mpf_EGlpNumCopy (b->bounds[active->bound_cnt], t);
	    b->bound_cnt = active->bound_cnt + 1;

	    rval = mpf_ILLutil_priority_insert (minf->que, (void *) b, &lpval,
		&(b->handle));
	    ILL_CLEANUP_IF (rval);

	    mpf_put_bbnode (minf, b);
	    minf->activenodes++;
	}
    }
    minf->totalnodes++;

    if (bdir == 'D') {
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'U', oldt);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'L', oldt);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    mpf_EGlpNumClearVar (t);
    mpf_EGlpNumClearVar (lpval);
    mpf_EGlpNumClearVar (oldt);
    return rval;
}

static int mpf_fix_variables (mpf_lpinfo * lp,
      mpf_t * bestval,
      mpf_bbnode * b,
      mpf_t * wupper,
      mpf_t * wlower,
      int *hit)
{
    int rval = 0;
    int i, nnew = 0;
    int nstruct = lp->O->nstruct;
    mpf_t delta, lpval;
    int *new_indx = 0;
    char *new_lu = 0;
    mpf_t *new_bounds = 0;
    mpf_t *dj = 0;
    mpf_EGlpNumInitVar (delta);
    mpf_EGlpNumInitVar (lpval);

    *hit = 0;

    if (mpf_EGlpNumIsLess (*bestval, mpf_ILL_MAXDOUBLE)) {
	rval = mpf_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);
	/* delta = bestval - lpval + mpf_ILL_INTTOL; */
	mpf_EGlpNumCopy (delta, *bestval);
	mpf_EGlpNumSubTo (delta, lpval);
	mpf_EGlpNumAddTo (delta, mpf_ILL_INTTOL);

	ILL_SAFE_MALLOC (new_indx, nstruct, int);
	ILL_SAFE_MALLOC (new_lu, nstruct, char);
	dj = mpf_EGlpNumAllocArray (nstruct);
	new_bounds = mpf_EGlpNumAllocArray (nstruct);

	rval = mpf_ILLlib_solution (lp, 0, 0, 0, 0, 0, dj);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < nstruct; i++) {
	    if (lp->O->intmarker[i]) {
		if (mpf_EGlpNumIsNeqq (wlower[i], wupper[i])) {
		    if (mpf_EGlpNumIsLess (delta, dj[i])) {
			mpf_EGlpNumSubTo (wupper[i], mpf_oneLpNum);
			rval = mpf_ILLlib_chgbnd (lp, i, 'U', wupper[i]);
			ILL_CLEANUP_IF (rval);
			new_indx[nnew] = i;
			new_lu[nnew] = 'U';
			mpf_EGlpNumCopy (new_bounds[nnew], wupper[i]);
			nnew++;
		    }
		    /* if (-dj[i] > delta) */
		    mpf_EGlpNumSign (delta);
		    if (mpf_EGlpNumIsLess (delta, dj[i])) {
			mpf_EGlpNumAddTo (wlower[i], mpf_oneLpNum);
			rval = mpf_ILLlib_chgbnd (lp, i, 'L', wlower[i]);
			ILL_CLEANUP_IF (rval);
			new_indx[nnew] = i;
			new_lu[nnew] = 'L';
			mpf_EGlpNumCopy (new_bounds[nnew], wlower[i]);
			nnew++;
		    }
		    mpf_EGlpNumSign (delta);
		}
	    }
	}

	if (nnew) {
	    b->bound_indx = EGrealloc (b->bound_indx, sizeof (int) * (b->bound_cnt + nnew));
	    /*
				    rval = ILLutil_reallocrus_count ((void **) &(b->bound_indx), b->bound_cnt + nnew, sizeof (int));
				    ILL_CLEANUP_IF (rval);
	    */
	    b->lu = EGrealloc (b->lu, sizeof (char) * (b->bound_cnt + nnew));
	    /*
				    rval = ILLutil_reallocrus_count ((void **) &(b->lu), b->bound_cnt + nnew, sizeof (char));
				    ILL_CLEANUP_IF (rval);
	    */
	    mpf_EGlpNumReallocArray (&(b->bounds), b->bound_cnt + nnew);
	    for (i = 0; i < nnew; i++) {
		b->bound_indx[b->bound_cnt + i] = new_indx[i];
		b->lu[b->bound_cnt + i] = new_lu[i];
		mpf_EGlpNumCopy (b->bounds[b->bound_cnt + i], new_bounds[i]);
	    }
	    b->bound_cnt += nnew;
	}
    }
    *hit = nnew;

CLEANUP:

    ILL_IFFREE (new_indx, int);
    ILL_IFFREE (new_lu, char);
    mpf_EGlpNumFreeArray (dj);
    mpf_EGlpNumFreeArray (new_bounds);
    mpf_EGlpNumClearVar (delta);
    mpf_EGlpNumClearVar (lpval);
    return rval;
}

static void mpf_best_bbnode (mpf_mipinfo * minf,
      mpf_bbnode ** best)
{
#if 0
    mpf_bbnode *b;
    double bestval = mpf_ILL_MAXDOUBLE;

    for (b = minf->head_bbnode.next; b; b = b->next) {
	if (b->bound < bestval) {
	    *best = b;
	    bestval = b->bound;
	}
    }
#endif

    mpf_t val;
    mpf_EGlpNumInitVar (val);
    mpf_ILLutil_priority_deletemin (minf->que, &val, (void **) best);
    mpf_EGlpNumClearVar (val);
}

static void mpf_put_bbnode (mpf_mipinfo * minf,
      mpf_bbnode * b)
{
    b->next = minf->head_bbnode.next;
    b->prev = &(minf->head_bbnode);
    if (b->next)
	b->next->prev = b;
    minf->head_bbnode.next = b;
}

static void mpf_remove_bbnode (mpf_bbnode * b)
{
    b->prev->next = b->next;
    if (b->next)
	b->next->prev = b->prev;
}

static int mpf_find_branch (mpf_mipinfo * minf,
      mpf_t * x,
      mpf_t * lpval,
      int *bvar)
{
    mpf_lpinfo *lp = minf->lp;
    int rval = 0;

    switch (minf->branching_rule) {
    case mpf_PENALTYBRANCH:
	rval = mpf_find_penalty_branch (lp, minf->pinf, x, minf->downpen,
	    minf->uppen, lpval, bvar);
	ILL_CLEANUP_IF (rval);
	break;
    case mpf_FIRSTBRANCH:
	mpf_find_first_branch (lp, x, bvar);
	break;
    case mpf_MIDDLEBRANCH:
	mpf_find_middle_branch (lp, x, bvar);
	break;
    case mpf_STRONGBRANCH:
	rval = mpf_find_strong_branch (lp, minf->pinf, x, bvar);
	ILL_CLEANUP_IF (rval);
	break;
    default:
	fprintf (stderr, "Unknown branching rule.\n");
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_find_branch");
}

static void mpf_find_first_branch (mpf_lpinfo * lp,
      mpf_t * x,
      int *bvar)
{
    int i, ibest = -1;
    mpf_ILLlpdata *qslp = lp->O;
    mpf_t t;
    mpf_EGlpNumInitVar (t);

    for (i = 0; i < qslp->nstruct; i++) {
	if (qslp->intmarker[i]) {
	    /* t = mpf_ILLutil_our_frac (x[i]); */
	    mpf_EGlpNumFloor (t, x[i]);
	    mpf_EGlpNumSubTo (t, x[i]);
	    mpf_EGlpNumSign (t);
	    if ((mpf_EGlpNumIsNeqZero (t, mpf_ILL_INTTOL)) &&
		(mpf_EGlpNumIsNeq (t, mpf_oneLpNum, mpf_ILL_INTTOL))) {
		ibest = i;
		break;
	    }
	}
    }
    *bvar = ibest;
    mpf_EGlpNumClearVar (t);
}

static void mpf_find_middle_branch (mpf_lpinfo * lp,
      mpf_t * x,
      int *bvar)
{
    int i, ibest = -1;
    mpf_t t, tbest;
    mpf_ILLlpdata *qlp = lp->O;
    mpf_EGlpNumInitVar (t);
    mpf_EGlpNumInitVar (tbest);
    mpf_EGlpNumSet (tbest, 0.5);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = mpf_ILLutil_our_frac (x[i]) - 0.5; if (t < 0.0) t = -t; */
	    mpf_EGlpNumFloor (t, x[i]);
	    mpf_EGlpNumMultUiTo (t, 2);
	    mpf_EGlpNumSubTo (t, mpf_oneLpNum);
	    mpf_EGlpNumDivUiTo (t, 2);
	    if (mpf_EGlpNumIsLess (t, mpf_zeroLpNum))
		mpf_EGlpNumSign (t);
	    /* if (t < tbest) */
	    if (mpf_EGlpNumIsLess (t, tbest)) {
		mpf_EGlpNumCopy (tbest, t);
		ibest = i;
	    }
	}
    }

    /* if (tbest < (0.5 - mpf_ILL_INTTOL)) */
    mpf_EGlpNumAddTo (tbest, mpf_ILL_INTTOL);
    if (mpf_EGlpNumIsLessDbl (tbest, 0.5)) {
	*bvar = ibest;
    } else {
	*bvar = -1;
    }
    mpf_EGlpNumClearVar (t);
    mpf_EGlpNumClearVar (tbest);
}

static int mpf_find_penalty_branch (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * x,
      mpf_t * downpen,
      mpf_t * uppen,
      mpf_t * lpval,
      int *bvar)
{
    int rval = 0;
    int i, k, ibest = -1, ncand = 0, nneed = 0;
    mpf_ILLlpdata *qslp = lp->O;
    int *candidatelist = 0;
    int *needlist = 0;
    mpf_t *fval = 0;
    mpf_t *xlist = 0;
    mpf_t *newdown = 0;
    mpf_t *newup = 0;
    mpf_t a, t, tbest;
    mpf_EGlpNumInitVar (a);
    mpf_EGlpNumInitVar (t);
    mpf_EGlpNumInitVar (tbest);
    mpf_EGlpNumCopy (tbest, mpf_ILL_MINDOUBLE);

    ILL_SAFE_MALLOC (candidatelist, qslp->nstruct, int);
    ILL_SAFE_MALLOC (needlist, qslp->nstruct, int);
    fval = mpf_EGlpNumAllocArray (qslp->nstruct);
    xlist = mpf_EGlpNumAllocArray (qslp->nstruct);
    for (i = 0; i < qslp->nstruct; i++) {
	if (qslp->intmarker[i]) {
	    /* fval[i] = x[i] - floor(x[i]); */
	    mpf_EGlpNumFloor (fval[i], x[i]);
	    mpf_EGlpNumSubTo (fval[i], x[i]);
	    mpf_EGlpNumSign (fval[i]);
	    if ((mpf_EGlpNumIsNeqZero (fval[i], mpf_ILL_INTTOL)) &&
		(mpf_EGlpNumIsNeq (fval[i], mpf_oneLpNum, mpf_ILL_INTTOL))) {
		candidatelist[ncand++] = i;
		/* if (downpen[i] == -1.0) */
		mpf_EGlpNumSign (downpen[i]);
		if (mpf_EGlpNumIsEqqual (downpen[i], mpf_oneLpNum)) {
		    mpf_EGlpNumCopy (xlist[nneed], x[i]);
		    needlist[nneed++] = i;
		}
		mpf_EGlpNumSign (downpen[i]);
	    }
	}
    }

    if (nneed > 0) {
	newdown = mpf_EGlpNumAllocArray (nneed);
	newup = mpf_EGlpNumAllocArray (nneed);
	rval = mpf_ILLlib_strongbranch (lp, pinf, needlist, nneed,
	    0, newdown, newup,
	    5 * mpf_STRONG_PIVOTS, mpf_ILL_MAXDOUBLE);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < nneed; i++) {
	    k = needlist[i];
	    /* uppen[k] = (newup[i] - lpval) / (1.0 - fval[k]); */
	    mpf_EGlpNumCopyDiff (uppen[k], newup[i], *lpval);
	    mpf_EGlpNumCopyDiff (downpen[k], mpf_oneLpNum, fval[k]);
	    mpf_EGlpNumDivTo (uppen[k], downpen[k]);
	    /* downpen[k] = (newdown[i] - lpval) / fval[k]; */
	    mpf_EGlpNumCopyDiffRatio (downpen[k], newdown[i], *lpval, fval[k]);

	}
    }
    for (i = 0; i < ncand; i++) {
	k = candidatelist[i];
	/* t = mpf_ILL_BRANCH_PENALTY_VAL (downpen[k], uppen[k], fval[k]); */
	mpf_EGlpNumCopy (t, downpen[k]);
	mpf_EGlpNumMultTo (t, fval[k]);
	mpf_EGlpNumCopyDiff (a, mpf_oneLpNum, fval[k]);
	mpf_EGlpNumMultTo (a, uppen[k]);
	if (mpf_EGlpNumIsLess (t, a)) {
	    mpf_EGlpNumMultUiTo (t, mpf_ILL_BRANCH_PENALTY_WEIGHT);
	    mpf_EGlpNumAddTo (t, a);
	} else {
	    mpf_EGlpNumMultUiTo (a, mpf_ILL_BRANCH_PENALTY_WEIGHT);
	    mpf_EGlpNumAddTo (t, a);
	}
	mpf_EGlpNumDivUiTo (t, mpf_ILL_BRANCH_PENALTY_WEIGHT + 1);

	if (mpf_EGlpNumIsLess (tbest, t)) {
	    mpf_EGlpNumCopy (tbest, t);
	    ibest = k;
	}
    }

    *bvar = ibest;

CLEANUP:

    mpf_EGlpNumClearVar (a);
    mpf_EGlpNumClearVar (t);
    mpf_EGlpNumClearVar (tbest);
    mpf_EGlpNumFreeArray (newdown);
    mpf_EGlpNumFreeArray (newup);
    mpf_EGlpNumFreeArray (fval);
    mpf_EGlpNumFreeArray (xlist);
    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (needlist, int);
    return rval;
}

static int mpf_find_strong_branch (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * x,
      int *bvar)
{
    int rval = 0;
    int i, ibest = -1, ncand = 0;
    int maxtrys = mpf_STRONG_CANDIDATES;
    mpf_t t, tbest;
    mpf_ILLlpdata *qlp = lp->O;
    int *candidatelist = 0;
    int *newlist = 0;
    int *perm = 0;
    mpf_t *tval = 0;
    mpf_t *xlist = 0;
    mpf_t *downpen = 0;
    mpf_t *uppen = 0;
    ILLrandstate rstate;
    mpf_EGlpNumInitVar (t);
    mpf_EGlpNumInitVar (tbest);
    mpf_EGlpNumCopy (tbest, mpf_ILL_MINDOUBLE);

    ILLutil_sprand (999, &rstate);
    ILL_SAFE_MALLOC (candidatelist, qlp->nstruct, int);
    tval = mpf_EGlpNumAllocArray (qlp->nstruct);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = mpf_ILLutil_our_frac (x[i]) - 0.5; if (t < 0.0) t = -t; */
	    mpf_EGlpNumFloor (t, x[i]);
	    mpf_EGlpNumSubTo (t, x[i]);
	    mpf_EGlpNumSign (t);
	    mpf_EGlpNumMultUiTo (t, 2);
	    mpf_EGlpNumSubTo (t, mpf_oneLpNum);
	    if (mpf_EGlpNumIsLess (t, mpf_zeroLpNum))
		mpf_EGlpNumSign (t);
	    /* if (t < (0.5 - mpf_ILL_INTTOL)) */
	    if (mpf_EGlpNumIsNeq (t, mpf_oneLpNum, mpf_ILL_INTTOL)) {
		candidatelist[ncand] = i;
		mpf_EGlpNumDivUiTo (t, 2);
		mpf_EGlpNumCopy (tval[ncand++], t);
	    }
	}
    }

    if (ncand > 0) {
	if (ncand > maxtrys) {
	    ILL_SAFE_MALLOC (perm, ncand, int);

	    for (i = 0; i < ncand; i++) {
		perm[i] = i;
	    }
	    mpf_ILLutil_EGlpNum_rselect (perm, 0, ncand - 1, maxtrys, tval, &rstate);

	    ILL_SAFE_MALLOC (newlist, maxtrys, int);

	    for (i = 0; i < maxtrys; i++) {
		newlist[i] = candidatelist[perm[i]];
	    }
	    ILL_IFFREE (candidatelist, int);
	    candidatelist = newlist;
	    newlist = 0;
	    ncand = maxtrys;
	}
	downpen = mpf_EGlpNumAllocArray (ncand);
	uppen = mpf_EGlpNumAllocArray (ncand);
	xlist = mpf_EGlpNumAllocArray (ncand);

	for (i = 0; i < ncand; i++) {
	    mpf_EGlpNumCopy (xlist[i], x[candidatelist[i]]);
	}

	rval = mpf_ILLlib_strongbranch (lp, pinf, candidatelist, ncand,
	    0, downpen, uppen, mpf_STRONG_PIVOTS,
	    mpf_ILL_MAXDOUBLE);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < ncand; i++) {
	    /* t = mpf_ILL_BRANCH_STRONG_VAL (downpen[i], uppen[i]); */
	    if (mpf_EGlpNumIsLess (downpen[i], uppen[i])) {
		mpf_EGlpNumCopy (t, downpen[i]);
		mpf_EGlpNumMultUiTo (t, mpf_ILL_BRANCH_STRONG_WEIGHT);
		mpf_EGlpNumAddTo (t, uppen[i]);
	    } else {
		mpf_EGlpNumCopy (t, uppen[i]);
		mpf_EGlpNumMultUiTo (t, mpf_ILL_BRANCH_STRONG_WEIGHT);
		mpf_EGlpNumAddTo (t, downpen[i]);
	    }
	    mpf_EGlpNumDivUiTo (t, mpf_ILL_BRANCH_STRONG_WEIGHT + 1);
	    if (mpf_EGlpNumIsLess (tbest, t)) {
		mpf_EGlpNumCopy (tbest, t);
		ibest = candidatelist[i];
	    }
	}
    }
    *bvar = ibest;


CLEANUP:

    mpf_EGlpNumClearVar (t);
    mpf_EGlpNumClearVar (tbest);
    mpf_EGlpNumFreeArray (tval);
    mpf_EGlpNumFreeArray (xlist);
    mpf_EGlpNumFreeArray (uppen);
    mpf_EGlpNumFreeArray (downpen);
    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (newlist, int);
    ILL_IFFREE (perm, int);

    ILL_RETURN (rval, "mpf_find_strong_branch");
}

static void mpf_check_integral (mpf_lpinfo * lp,
      mpf_t * x,
      int *yesno)
{
    int i;
    mpf_t t;
    mpf_ILLlpdata *qlp = lp->O;
    mpf_EGlpNumInitVar (t);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = mpf_ILLutil_our_frac (x[i]); */
	    mpf_EGlpNumFloor (t, x[i]);
	    mpf_EGlpNumSubTo (t, x[i]);
	    mpf_EGlpNumSign (t);
	    /* if (t > mpf_ILL_INTTOL && t < 1.0 - mpf_ILL_INTTOL) */
	    if ((mpf_EGlpNumIsNeqZero (t, mpf_ILL_INTTOL)) &&
		(mpf_EGlpNumIsNeq (t, mpf_oneLpNum, mpf_ILL_INTTOL))) {
		*yesno = 0;
		mpf_EGlpNumClearVar (t);
		return;
	    }
	}
    }

    *yesno = 1;
    mpf_EGlpNumClearVar (t);
}

static int mpf_plunge (mpf_mipinfo * minf)
{
    int rval = 0;
    int i, status;
    mpf_lpinfo *lp = minf->lp;
    mpf_ILLlpdata *qlp = minf->lp->O;
    mpf_t *oldlower = 0;
    mpf_t *oldupper = 0;

    if (minf->watch) {
	printf ("Plunging ...\n");
	fflush (stdout);
    }
    oldlower = mpf_EGlpNumAllocArray (qlp->nstruct);
    oldupper = mpf_EGlpNumAllocArray (qlp->nstruct);

    for (i = 0; i < qlp->nstruct; i++) {
	mpf_EGlpNumCopy (oldlower[i], minf->lower[i]);
	mpf_EGlpNumCopy (oldupper[i], minf->upper[i]);
    }

    rval = mpf_plunge_work (minf, 0);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
	rval = mpf_ILLlib_chgbnd (lp, i, 'L', oldlower[i]);
	ILL_CLEANUP_IF (rval);
	rval = mpf_ILLlib_chgbnd (lp, i, 'U', oldupper[i]);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumCopy (minf->lower[i], oldlower[i]);
	mpf_EGlpNumCopy (minf->upper[i], oldupper[i]);
    }

    rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    mpf_EGlpNumFreeArray (oldlower);
    mpf_EGlpNumFreeArray (oldupper);

    ILL_RETURN (rval, "mpf_plunge");
}

static int mpf_plunge_work (mpf_mipinfo * minf,
      int depth)
{
    int rval = 0;
    int bvar, status, count;
    mpf_t lpval, val0, val1, int_tol;
    mpf_lpinfo *lp = minf->lp;
    mpf_EGlpNumInitVar (lpval);
    mpf_EGlpNumInitVar (val0);
    mpf_EGlpNumInitVar (val1);
    mpf_EGlpNumInitVar (int_tol);
    mpf_EGlpNumSet (int_tol, 0.001);

    rval = mpf_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = mpf_round_variables (minf, &count, &int_tol /* 0.001 */ );
    ILL_CLEANUP_IF (rval);
    if (count) {
	rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	if (status != QS_LP_OPTIMAL) {
	    goto CLEANUP;
	}
	rval = mpf_ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);
    }
    mpf_find_middle_branch (lp, minf->x, &bvar);
    if (bvar == -1) {
	rval = mpf_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);

	if (mpf_EGlpNumIsLess (lpval, minf->value)) {
	    printf ("Plunge Integral Solution: %.6f (Depth: %d)\n",
		mpf_EGlpNumToLf (lpval), depth);
	    fflush (stdout);

	    mpf_EGlpNumCopy (minf->value, lpval);
	    mpf_EGlpNumCopyDiff (minf->objectivebound, lpval, mpf_ILL_INTTOL);
	    mpf_copy_x (lp->O->nstruct, minf->x, minf->bestx);
	}
	goto CLEANUP;
    }
    mpf_EGlpNumOne (minf->lower[bvar]);
    rval = mpf_ILLlib_chgbnd (lp, bvar, 'L', mpf_oneLpNum);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the mpf_plunge LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
	mpf_EGlpNumCopy (val1, mpf_ILL_MAXDOUBLE);
    } else if (status == QS_LP_OPTIMAL) {
	rval = mpf_ILLlib_objval (lp, 0, &val1);
	ILL_CLEANUP_IF (rval);
    } else {
	ILL_CLEANUP;
    }

    rval = mpf_ILLlib_chgbnd (lp, bvar, 'L', mpf_zeroLpNum);
    ILL_CLEANUP_IF (rval);
    mpf_EGlpNumZero (minf->lower[bvar]);

    mpf_EGlpNumZero (minf->upper[bvar]);
    rval = mpf_ILLlib_chgbnd (lp, bvar, 'U', mpf_zeroLpNum);
    ILL_CLEANUP_IF (rval);
    rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the mpf_plunge LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
	mpf_EGlpNumCopy (val0, mpf_ILL_MAXDOUBLE);
    } else if (status == QS_LP_OPTIMAL) {
	rval = mpf_ILLlib_objval (lp, 0, &val0);
	ILL_CLEANUP_IF (rval);
    } else {
	ILL_CLEANUP;
    }

    rval = mpf_ILLlib_chgbnd (lp, bvar, 'U', mpf_oneLpNum);
    ILL_CLEANUP_IF (rval);
    mpf_EGlpNumCopy (minf->upper[bvar], mpf_oneLpNum);

    if (mpf_EGlpNumIsEqqual (val0, mpf_ILL_MAXDOUBLE) &&
	mpf_EGlpNumIsEqqual (val1, mpf_ILL_MAXDOUBLE)) {
	ILL_CLEANUP;
    }
    if (mpf_EGlpNumIsLess (val0, val1)) {
	mpf_EGlpNumZero (minf->upper[bvar]);
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'U', mpf_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	rval = mpf_plunge_work (minf, depth + 1);
	ILL_CLEANUP_IF (rval);
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'U', mpf_oneLpNum);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumOne (minf->upper[bvar]);
    } else {
	mpf_EGlpNumOne (minf->lower[bvar]);
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'L', mpf_oneLpNum);
	ILL_CLEANUP_IF (rval);
	rval = mpf_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	rval = mpf_plunge_work (minf, depth + 1);
	ILL_CLEANUP_IF (rval);
	rval = mpf_ILLlib_chgbnd (lp, bvar, 'L', mpf_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumZero (minf->lower[bvar]);
    }

CLEANUP:

    mpf_EGlpNumClearVar (lpval);
    mpf_EGlpNumClearVar (val0);
    mpf_EGlpNumClearVar (val1);
    mpf_EGlpNumClearVar (int_tol);
    ILL_RETURN (rval, "mpf_plunge_work");
}

static int mpf_round_variables (mpf_mipinfo * minf,
      int *count,
      mpf_t * tol)
{
    int rval = 0;
    int i, hit = 0;
    mpf_lpinfo *lp = minf->lp;
    mpf_ILLlpdata *qlp = lp->O;

    *count = 0;

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    if (mpf_EGlpNumIsNeqq (minf->lower[i], minf->upper[i])) {
		if (mpf_EGlpNumIsLess (minf->x[i], *tol)) {
		    mpf_EGlpNumZero (minf->upper[i]);
		    rval = mpf_ILLlib_chgbnd (lp, i, 'U', mpf_zeroLpNum);
		    ILL_CLEANUP_IF (rval);
		    hit++;
		} else if (mpf_EGlpNumIsEqual (minf->x[i], mpf_oneLpNum, *tol)) {
		    mpf_EGlpNumOne (minf->lower[i]);
		    rval = mpf_ILLlib_chgbnd (lp, i, 'L', mpf_oneLpNum);
		    ILL_CLEANUP_IF (rval);
		    hit++;
		}
	    }
	}
    }
    *count = hit;

CLEANUP:

    ILL_RETURN (rval, "mpf_round_variables");
}

static void mpf_copy_x (int nstruct,
      mpf_t * from_x,
      mpf_t * to_x)
{
    int j;

    for (j = 0; j < nstruct; j++) {
	mpf_EGlpNumCopy (to_x[j], from_x[j]);
    }
}

static void mpf_init_mipinfo (mpf_mipinfo * minf)
{
    if (minf) {
	minf->depth = 0;
	minf->totalnodes = 0;
	minf->activenodes = 0;
	minf->totalpivots = 0;
	minf->lastpivots = 0;
	minf->downpen = 0;
	minf->uppen = 0;
	minf->x = 0;
	minf->bestx = 0;
	minf->lower = 0;
	minf->upper = 0;
	minf->lp = 0;
	minf->pinf = 0;
	minf->head_bbnode.prev = 0;
	minf->head_bbnode.next = 0;
	minf->que = 0;
	minf->branching_rule = /* mpf_MIDDLEBRANCH */ mpf_STRONGBRANCH;
	minf->watch = 1;
	mpf_EGlpNumInitVar (minf->objectivebound);
	mpf_EGlpNumInitVar (minf->value);
	mpf_EGlpNumCopy (minf->objectivebound, mpf_ILL_MAXDOUBLE);
	mpf_EGlpNumCopy (minf->value, mpf_ILL_MAXDOUBLE);
	ILLptrworld_init (&minf->ptrworld);
    }
}

static void mpf_free_mipinfo (mpf_mipinfo * minf)
{
    int total, onlist;

    if (minf) {
	mpf_EGlpNumFreeArray (minf->downpen);
	mpf_EGlpNumFreeArray (minf->uppen);
	mpf_EGlpNumFreeArray (minf->x);
	mpf_EGlpNumFreeArray (minf->bestx);
	mpf_EGlpNumFreeArray (minf->lower);
	mpf_EGlpNumFreeArray (minf->upper);
	bbnode_listfree (&minf->ptrworld, minf->head_bbnode.next);
	if (bbnode_check_leaks (&minf->ptrworld, &total, &onlist)) {
	    fprintf (stderr, "WARNING: %d outstanding bbnodes\n", total - onlist);
	}
	ILLptrworld_delete (&minf->ptrworld);
	mpf_EGlpNumClearVar ((minf->objectivebound));
	mpf_EGlpNumClearVar ((minf->value));
	memset (minf, 0, sizeof (mpf_mipinfo));
	/* mpf_init_mipinfo (minf); */
    }
}

static void mpf_init_bbnode (mpf_bbnode * b)
{
    if (b) {
	b->next = 0;
	b->prev = 0;
	b->id = 0;
	b->depth = 0;
	b->handle = 0;
	b->cstat = 0;
	b->rstat = 0;
	b->rownorms = 0;
	b->bound_cnt = 0;
	b->bound_indx = 0;
	b->lu = 0;
	b->bounds = 0;
	mpf_EGlpNumInitVar ((b->bound));
	mpf_EGlpNumCopy (b->bound, mpf_ILL_MINDOUBLE);
    }
}

static void mpf_free_bbnode (mpf_bbnode * b)
{
    if (b) {
	mpf_EGlpNumFreeArray (b->rownorms);
	mpf_EGlpNumFreeArray (b->bounds);
	ILL_IFFREE (b->cstat, char);
	ILL_IFFREE (b->rstat, char);
	ILL_IFFREE (b->bound_indx, int);
	ILL_IFFREE (b->lu, char);
	mpf_EGlpNumClearVar ((b->bound));
	memset (b, 0, sizeof (mpf_bbnode));
    }
}
