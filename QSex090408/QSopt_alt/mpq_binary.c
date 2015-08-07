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

/* RCS_INFO = "$RCSfile: mpq_binary.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* Simple MIP Code to test LP Solver                    */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int mpq_ILLmip_bfs (mpq_lpinfo *lp, double *val, double *x)                   */
/* */
/* NOTES                                                                   */
/* */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpq_priority.h"
#include "mpq_sortrus.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdata.h"
#include "mpq_lpdefs.h"
#include "mpq_simplex.h"
#include "mpq_binary.h"
#include "mpq_price.h"
#include "mpq_lib.h"
#include "mpq_qstruct.h"
#include "mpq_qsopt.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* #define  mpq_ILL_INTTOL (0.000001) */
#define mpq_ILL_INTTOL mpq_PFEAS_TOLER

#define  mpq_STRONG_PIVOTS     (50)
#define  mpq_STRONG_CANDIDATES (10)

#define mpq_ILL_BRANCH_STRONG_WEIGHT (10)
#define mpq_ILL_BRANCH_STRONG_VAL(v0,v1)                                 \
    (((v0) < (v1) ? (mpq_ILL_BRANCH_STRONG_WEIGHT * (v0) + (v1))         \
                  : (mpq_ILL_BRANCH_STRONG_WEIGHT * (v1) + (v0)))        \
                    / (mpq_ILL_BRANCH_STRONG_WEIGHT + 1.0))

#define mpq_ILL_BRANCH_PENALTY_WEIGHT (2)
#define mpq_ILL_BRANCH_PENALTY_VAL(v0,v1,f)                              \
    (((v0)*(f) < (v1)*(1.0-(f)) ?                                    \
        (mpq_ILL_BRANCH_PENALTY_WEIGHT * (v0)*(f) + (v1)*(1.0-(f)))    \
      : (mpq_ILL_BRANCH_PENALTY_WEIGHT * (v1)*(1.0-(f)) + (v0)*(f)))    \
                    / (mpq_ILL_BRANCH_PENALTY_WEIGHT + 1.0))



#define mpq_FIRSTBRANCH  1
#define mpq_MIDDLEBRANCH 2
#define mpq_STRONGBRANCH 3
#define mpq_PENALTYBRANCH 4


typedef struct mpq_bbnode {
    struct mpq_bbnode *next;
    struct mpq_bbnode *prev;
    int id;
    int depth;
    int handle;
    mpq_t bound;
    char *cstat;
    char *rstat;
    mpq_t *rownorms;
    int rownorms_size;
    int bound_cnt;
    int *bound_indx;
    char *lu;
    mpq_t *bounds;
    int bounds_size;
}
  mpq_bbnode;

typedef struct mpq_mipinfo {
    int branching_rule;
    int watch;
    int depth;
    int totalnodes;
    int activenodes;
    int totalpivots;
    int lastpivots;
    int objsense;
    mpq_t objectivebound;
    mpq_t value;
    mpq_t *downpen;
    mpq_t *uppen;
    mpq_t *x;
    mpq_t *bestx;
    mpq_t *orig_lower;
    mpq_t *orig_upper;
    mpq_t *lower;
    mpq_t *upper;
    int nstruct;		/* size of all mpq_t arrays */
    mpq_lpinfo *lp;
    mpq_price_info *pinf;
    mpq_bbnode head_bbnode;
    mpq_ILLpriority *que;
    ILLptrworld ptrworld;
}
  mpq_mipinfo;


ILL_PTRWORLD_ROUTINES (mpq_bbnode, bbnodealloc, bbnode_bulkalloc, bbnodefree)
ILL_PTRWORLD_LISTFREE_ROUTINE (mpq_bbnode, bbnode_listfree, bbnodefree)
ILL_PTRWORLD_LEAKS_ROUTINE (mpq_bbnode, bbnode_check_leaks, depth, int)
static void mpq_cleanup_mip (mpq_mipinfo * minf),
  mpq_choose_initial_price (mpq_price_info * pinf),
  mpq_best_bbnode (mpq_mipinfo * minf,
      mpq_bbnode ** best),
  mpq_put_bbnode (mpq_mipinfo * minf,
      mpq_bbnode * b),
  mpq_remove_bbnode (mpq_bbnode * b),
  mpq_find_first_branch (mpq_lpinfo * lp,
      mpq_t * x,
      int *bvar),
  mpq_find_middle_branch (mpq_lpinfo * lp,
      mpq_t * x,
      int *bvar),
  mpq_check_integral (mpq_lpinfo * lp,
      mpq_t * x,
      int *yesno),
  mpq_copy_x (int nstruct,
      mpq_t * from_x,
      mpq_t * to_x),
  mpq_init_mipinfo (mpq_mipinfo * minf),
  mpq_free_mipinfo (mpq_mipinfo * minf),
  mpq_init_bbnode (mpq_bbnode * b),
  mpq_free_bbnode (mpq_bbnode * b);

static int mpq_startup_mip (mpq_mipinfo * minf,
      mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_t * lpval),
  mpq_run_bfs (mpq_mipinfo * minf),
  mpq_process_bfs_bbnode (mpq_mipinfo * minf,
      mpq_bbnode * b),
  mpq_child_work (mpq_mipinfo * minf,
      mpq_bbnode * active,
      int bvar,
      int bdir,
      mpq_t * cval,
      int *cp),
  mpq_fix_variables (mpq_lpinfo * lp,
      mpq_t * bestval,
      mpq_bbnode * b,
      mpq_t * wupper,
      mpq_t * wlower,
      int *hit),
  mpq_find_branch (mpq_mipinfo * minf,
      mpq_t * x,
      mpq_t * lpval,
      int *bvar),
  mpq_find_penalty_branch (mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_t * x,
      mpq_t * downpen,
      mpq_t * uppen,
      mpq_t * lpval,
      int *bvar),
  mpq_find_strong_branch (mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_t * x,
      int *bvar),
  mpq_plunge (mpq_mipinfo * minf),
  mpq_plunge_work (mpq_mipinfo * minf,
      int depth),
  mpq_round_variables (mpq_mipinfo * minf,
      int *count,
      mpq_t * tol);


static void mpq_choose_initial_price (mpq_price_info * pinf)
{
    pinf->pI_price = QS_PRICE_PSTEEP;
    pinf->pII_price = QS_PRICE_PSTEEP;
    pinf->dI_price = QS_PRICE_DSTEEP;
    pinf->dII_price = QS_PRICE_DSTEEP;
}

int mpq_ILLmip_bfs (mpq_lpinfo * lp,
      mpq_t * val,
      mpq_t * x)
{
    int tval, rval = 0;
    mpq_price_info pinf;
    mpq_mipinfo minf;
    mpq_bbnode *b;
    mpq_t lpval;
    double szeit = ILLutil_zeit ();
    mpq_EGlpNumInitVar (lpval);
    mpq_EGlpNumInitVar (pinf.htrigger);

    mpq_ILLprice_init_pricing_info (&pinf);
    mpq_init_mipinfo (&minf);

    if (!lp) {
	fprintf (stderr, "mpq_ILLmip_bfs called without an LP\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_startup_mip (&minf, lp, &pinf, &lpval);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (minf.que, 1, mpq_ILLpriority);
    rval = mpq_ILLutil_priority_init (minf.que, lp->O->nstruct + 1);
    ILL_CLEANUP_IF (rval);

    b = bbnodealloc (&minf.ptrworld);
    mpq_init_bbnode (b);
    b->depth = 0;
    b->id = minf.totalnodes++;
    mpq_EGlpNumCopy (b->bound, lpval);
    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
    rval = mpq_ILLlib_getbasis (lp, b->cstat, b->rstat);
    ILL_CLEANUP_IF (rval);

    if (pinf.dII_price == QS_PRICE_DSTEEP) {
	b->rownorms = mpq_EGlpNumAllocArray (lp->nrows);
	tval = mpq_ILLlib_getrownorms (lp, &pinf, b->rownorms);
	if (tval) {
	    printf ("Row norms not available\n");
	    fflush (stdout);
	    mpq_EGlpNumFreeArray (b->rownorms);
	}
    }
    rval = mpq_ILLutil_priority_insert (minf.que, (void *) b, &lpval, &(b->handle));
    ILL_CLEANUP_IF (rval);

    b->prev = &(minf.head_bbnode);
    b->next = 0;
    minf.head_bbnode.next = b;
    minf.activenodes++;

    minf.branching_rule = mpq_PENALTYBRANCH;

    rval = mpq_run_bfs (&minf);
    ILL_CLEANUP_IF (rval);

    printf ("Total Number of Nodes: %d\n", minf.totalnodes);
    printf ("Total Number of Pivots: %d\n", minf.totalpivots);
    printf ("BFS MIP Runing Time: %.2f seconds\n", ILLutil_zeit () - szeit);
    fflush (stdout);

    mpq_EGlpNumCopy (*val, minf.value);
    if (minf.objsense == mpq_ILL_MAX)
	mpq_EGlpNumSign (*val);

    if (x && mpq_EGlpNumIsNeqq (minf.value, mpq_ILL_MAXDOUBLE)) {
	mpq_copy_x (lp->O->nstruct, minf.bestx, x);
    }
CLEANUP:

    if (minf.que) {
	mpq_ILLutil_priority_free (minf.que);
	ILL_IFFREE (minf.que, mpq_ILLpriority);
    }
    mpq_cleanup_mip (&minf);
    mpq_free_mipinfo (&minf);
    mpq_ILLprice_free_pricing_info (&pinf);
    mpq_EGlpNumClearVar (lpval);
    mpq_EGlpNumClearVar (pinf.htrigger);
    ILL_RETURN (rval, "mpq_ILLmip_bfs");
}

static int mpq_startup_mip (mpq_mipinfo * minf,
      mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_t * lpval)
{
    int rval = 0;
    int i, col, status, intcount = 0;
    mpq_t val;
    mpq_ILLlpdata *qlp;
    mpq_EGlpNumInitVar (val);

    mpq_choose_initial_price (pinf);

    qlp = lp->O;

    rval = mpq_ILLlib_optimize (lp, 0, pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += mpq_ILLlib_iter (lp);

    rval = mpq_ILLlib_objval (lp, 0, &val);
    ILL_CLEANUP_IF (rval);

    printf ("LP Value: %.6f\n", mpq_EGlpNumToLf (val));
    fflush (stdout);
    if (lpval)
	mpq_EGlpNumCopy (*lpval, val);

    if (qlp->intmarker) {
	for (i = 0; i < qlp->nstruct; i++) {
	    if (qlp->intmarker[i]) {
		col = qlp->structmap[i];
		intcount++;
		if (mpq_EGlpNumIsEqqual (qlp->lower[col], mpq_ILL_MINDOUBLE)
		    || mpq_EGlpNumIsEqqual (qlp->upper[col], mpq_ILL_MAXDOUBLE)) {
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
	mpq_ILLlp_sinfo_free (qlp->sinfo);
	ILL_IFFREE (qlp->sinfo, mpq_ILLlp_sinfo);
    }
    minf->lp = lp;
    minf->pinf = pinf;
    minf->objsense = qlp->objsense;
    if (qlp->objsense == mpq_ILL_MAX) {	/* MIP codes work with min */
	for (i = 0; i < lp->ncols; i++) {
	    mpq_EGlpNumCopyNeg (qlp->obj[i], qlp->obj[i]);
	}
	qlp->objsense = mpq_ILL_MIN;
    }
    minf->x = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->bestx = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->lower = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->upper = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->orig_lower = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->orig_upper = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->downpen = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->uppen = mpq_EGlpNumAllocArray (qlp->nstruct);
    minf->nstruct = qlp->nstruct;

    rval = mpq_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
	mpq_EGlpNumCopy (minf->lower[i], qlp->lower[i]);
	mpq_EGlpNumCopy (minf->upper[i], qlp->upper[i]);
	mpq_EGlpNumCopy (minf->orig_lower[i], qlp->lower[i]);
	mpq_EGlpNumCopy (minf->orig_upper[i], qlp->upper[i]);
	mpq_EGlpNumOne (minf->downpen[i]);
	mpq_EGlpNumOne (minf->uppen[i]);
	mpq_EGlpNumSign (minf->downpen[i]);
	mpq_EGlpNumSign (minf->uppen[i]);
    }


CLEANUP:

    mpq_EGlpNumClearVar (val);
    ILL_RETURN (rval, "mpq_startup_mip");
}

static void mpq_cleanup_mip (mpq_mipinfo * minf)
{
    int i;
    mpq_ILLlpdata *qslp;

    if (minf && minf->lp) {
	qslp = minf->lp->O;
	if (minf->objsense == mpq_ILL_MAX) {
	    for (i = 0; i < minf->lp->ncols; i++) {
		mpq_EGlpNumSign (qslp->obj[i]);
	    }
	    qslp->objsense = mpq_ILL_MIN;
	}
    }
}

static int mpq_run_bfs (mpq_mipinfo * minf)
{
    int rval = 0;
    mpq_bbnode *b;

    while (minf->head_bbnode.next) {
	mpq_best_bbnode (minf, &b);
	rval = mpq_process_bfs_bbnode (minf, b);
	ILL_CLEANUP_IF (rval);
	mpq_remove_bbnode (b);
	mpq_free_bbnode (b);
	bbnodefree (&minf->ptrworld, b);
	minf->activenodes--;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_run_bfs");
}

static int mpq_process_bfs_bbnode (mpq_mipinfo * minf,
      mpq_bbnode * active)
{
    mpq_lpinfo *lp = minf->lp;
    mpq_ILLlp_basis B;
    int status, bvar = 0;
    int i, j, hit, dnp = 0, upp = 0;
    int nstruct = lp->O->nstruct;
    mpq_t t, lpval, dnval, upval;
    mpq_t *wupper = 0;
    mpq_t *wlower = 0;
    int rval = 0;
    mpq_EGlpNumInitVar (t);
    mpq_EGlpNumInitVar (lpval);
    mpq_EGlpNumInitVar (dnval);
    mpq_EGlpNumInitVar (upval);

    mpq_ILLlp_basis_init (&B);

    if (minf->watch > 1) {
	printf ("Node %4d: %.3f", active->id, mpq_EGlpNumToLf (active->bound));
	if (mpq_EGlpNumIsNeqq (minf->value, mpq_ILL_MAXDOUBLE))
	    printf (" %.3f", mpq_EGlpNumToLf (minf->value));
	else
	    printf ("  None");
	printf (", Active %d ", minf->activenodes);
	fflush (stdout);
    } else if (minf->watch == 1) {
	if (minf->lastpivots > 1000) {
	    minf->lastpivots = 0;
	    printf ("Pivots %d, Active Nodes %d, Bound %.3f, Soln ",
		minf->totalpivots, minf->activenodes,
		mpq_EGlpNumToLf (active->bound));
	    if (!mpq_EGlpNumIsLess (minf->value, mpq_ILL_MAXDOUBLE))
		printf ("%.3f", mpq_EGlpNumToLf (minf->value));
	    else
		printf ("None\n");
	}
    }
    if (mpq_EGlpNumIsLeq (minf->objectivebound, active->bound)) {
	if (minf->watch > 1) {
	    printf ("  Node can be purged\n");
	    fflush (stdout);
	}
	goto CLEANUP;
    }
    /* Set the LP bounds for the mpq_node. */

    wlower = mpq_EGlpNumAllocArray (nstruct);
    wupper = mpq_EGlpNumAllocArray (nstruct);

    for (i = 0; i < nstruct; i++) {
	mpq_EGlpNumCopy (wlower[i], minf->orig_lower[i]);
	mpq_EGlpNumCopy (wupper[i], minf->orig_upper[i]);
    }
    for (i = 0; i < active->bound_cnt; i++) {
	j = active->bound_indx[i];
	if (active->lu[i] == 'L')
	    mpq_EGlpNumCopy (wlower[j], active->bounds[i]);
	else
	    mpq_EGlpNumCopy (wupper[j], active->bounds[i]);
    }

    if (active->bound_cnt > 0) {
	rval = mpq_ILLlib_chgbnds (lp, active->bound_cnt, active->bound_indx,
	    active->lu, active->bounds);
	ILL_CLEANUP_IF (rval);
    }
    /* Solve the LP. */

    rval = mpq_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
	active->rstat);
    ILL_CLEANUP_IF (rval);
    if (active->rownorms) {
	B.rownorms = mpq_EGlpNumAllocArray (lp->nrows);
	for (i = 0; i < lp->nrows; i++) {
	    mpq_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
	}
    }
    rval = mpq_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += mpq_ILLlib_iter (lp);
    minf->lastpivots += mpq_ILLlib_iter (lp);

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
	    mpq_EGlpNumCopy (minf->lower[i], wlower[i]);
	    mpq_EGlpNumCopy (minf->upper[i], wupper[i]);
	}
	rval = mpq_plunge (minf);
	ILL_CLEANUP_IF (rval);
    }
    /* Fix variables. */

    if (mpq_EGlpNumIsLess (minf->value, mpq_ILL_MAXDOUBLE)) {
	rval = mpq_fix_variables (lp, &(minf->value), active, wupper, wlower, &hit);
	ILL_CLEANUP_IF (rval);

	if (hit) {
	    rval = mpq_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
	    ILL_CLEANUP_IF (rval);

	    minf->totalpivots += mpq_ILLlib_iter (lp);
	    minf->lastpivots += mpq_ILLlib_iter (lp);

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

    rval = mpq_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_objval (lp, 0, &lpval);
    ILL_CLEANUP_IF (rval);

    rval = mpq_find_branch (minf, minf->x, &lpval, &bvar);
    ILL_CLEANUP_IF (rval);

    if (bvar == -1) {
	printf ("Found integral solution: %f\n", mpq_EGlpNumToLf (lpval));
	if (mpq_EGlpNumIsLess (lpval, minf->value)) {
	    mpq_EGlpNumCopy (minf->value, lpval);
	    mpq_EGlpNumCopy (minf->objectivebound, lpval);
	    mpq_EGlpNumSubTo (minf->objectivebound, mpq_ILL_INTTOL);
	    mpq_copy_x (nstruct, minf->x, minf->bestx);
	}
    } else {
	/* Create down child */

	rval = mpq_child_work (minf, active, bvar, 'D', &dnval, &dnp);
	ILL_CLEANUP_IF (rval);

	/* Restore parent basis */

	rval = mpq_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
	    active->rstat);
	ILL_CLEANUP_IF (rval);
	if (active->rownorms) {
	    B.rownorms = mpq_EGlpNumAllocArray (lp->nrows);
	    for (i = 0; i < lp->nrows; i++) {
		mpq_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
	    }
	}
	/* Create up child */

	rval = mpq_child_work (minf, active, bvar, 'U', &upval, &upp);
	ILL_CLEANUP_IF (rval);

	if (minf->watch > 1) {
	    if (mpq_EGlpNumIsEqqual (dnval, mpq_ILL_MAXDOUBLE)) {
		printf ("DN->XXX");
	    } else {
		printf ("DN->%.3f%c", mpq_EGlpNumToLf (dnval), dnp ? 'X' : ' ');
	    }
	    if (mpq_EGlpNumIsEqqual (upval, mpq_ILL_MAXDOUBLE)) {
		printf ("UP->XXX\n");
	    } else {
		printf ("UP->%.3f%c\n", mpq_EGlpNumToLf (upval), upp ? 'X' : ' ');
	    }
	    fflush (stdout);
	}
    }

    /* Set the LP bounds back to original values */

    for (i = 0; i < active->bound_cnt; i++) {
	if (active->lu[i] == 'L')
	    mpq_EGlpNumCopy (t, minf->orig_lower[active->bound_indx[i]]);
	else
	    mpq_EGlpNumCopy (t, minf->orig_upper[active->bound_indx[i]]);

	rval = mpq_ILLlib_chgbnd (lp, active->bound_indx[i], active->lu[i], t);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    mpq_EGlpNumFreeArray (wlower);
    mpq_EGlpNumFreeArray (wupper);
    mpq_ILLlp_basis_free (&B);
    mpq_EGlpNumClearVar (t);
    mpq_EGlpNumClearVar (lpval);
    mpq_EGlpNumClearVar (dnval);
    mpq_EGlpNumClearVar (upval);
    ILL_RETURN (rval, "mpq_process_bfs_bbnode");
}

static int mpq_child_work (mpq_mipinfo * minf,
      mpq_bbnode * active,
      int bvar,
      int bdir,
      mpq_t * cval,
      int *cp)
{
    int tval, rval = 0;
    int i, status, intsol;
    mpq_t t, oldt, lpval;
    mpq_t *xi = &(minf->x[bvar]);
    mpq_lpinfo *lp = minf->lp;
    mpq_bbnode *b;
    mpq_EGlpNumInitVar (t);
    mpq_EGlpNumInitVar (lpval);
    mpq_EGlpNumInitVar (oldt);

    *cp = 0;

    if (bdir == 'D') {
	rval = mpq_ILLlib_getbnd (lp, bvar, 'U', &oldt);
	ILL_CLEANUP_IF (rval);
	mpq_EGlpNumFloor (t, *xi);
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'U', t);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = mpq_ILLlib_getbnd (lp, bvar, 'L', &oldt);
	ILL_CLEANUP_IF (rval);
	mpq_EGlpNumCeil (t, *xi);
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'L', t);
	ILL_CLEANUP_IF (rval);
    }

    rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += mpq_ILLlib_iter (lp);
    minf->lastpivots += mpq_ILLlib_iter (lp);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve Child LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    }
    if (status == QS_LP_INFEASIBLE) {
	mpq_EGlpNumCopy (*cval, mpq_ILL_MAXDOUBLE);
	*cp = 1;
    } else {
	rval = mpq_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);
	mpq_EGlpNumCopy (*cval, lpval);

	/* What about the x vector?  Bico - 020531 */

	mpq_check_integral (lp, minf->x, &intsol);
	if (intsol) {
	    if (mpq_EGlpNumIsLess (lpval, minf->value)) {
		printf ("Found integral solution: %f\n", mpq_EGlpNumToLf (lpval));
		mpq_EGlpNumCopy (minf->value, lpval);
		mpq_EGlpNumCopy (minf->objectivebound, lpval);
		mpq_EGlpNumSubTo (minf->objectivebound, mpq_ILL_INTTOL);
		mpq_copy_x (lp->O->nstruct, minf->x, minf->bestx);
	    }
	}
	if (mpq_EGlpNumIsLeq (minf->objectivebound, lpval)) {
	    *cp = 1;
	} else {
	    b = bbnodealloc (&minf->ptrworld);
	    mpq_init_bbnode (b);
	    b->depth = active->depth + 1;
	    b->id = minf->totalnodes;
	    mpq_EGlpNumCopy (b->bound, lpval);
	    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
	    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
	    rval = mpq_ILLlib_getbasis (lp, b->cstat, b->rstat);
	    ILL_CLEANUP_IF (rval);
	    if (minf->pinf->dII_price == QS_PRICE_DSTEEP) {
		b->rownorms = mpq_EGlpNumAllocArray (lp->nrows);
		tval = mpq_ILLlib_getrownorms (lp, minf->pinf, b->rownorms);
		if (tval) {
		    printf ("Row norms not available\n");
		    fflush (stdout);
		    printf ("A\n");
		    exit (1);
		    mpq_EGlpNumFreeArray (b->rownorms);
		}
	    }
	    ILL_SAFE_MALLOC (b->bound_indx, active->bound_cnt + 1, int);
	    ILL_SAFE_MALLOC (b->lu, active->bound_cnt + 1, char);
	    b->bounds = mpq_EGlpNumAllocArray (active->bound_cnt + 1);
	    for (i = 0; i < active->bound_cnt; i++) {
		b->bound_indx[i] = active->bound_indx[i];
		b->lu[i] = active->lu[i];
		mpq_EGlpNumCopy (b->bounds[i], active->bounds[i]);
	    }
	    b->bound_indx[active->bound_cnt] = bvar;
	    if (bdir == 'D')
		b->lu[active->bound_cnt] = 'U';
	    else
		b->lu[active->bound_cnt] = 'L';
	    mpq_EGlpNumCopy (b->bounds[active->bound_cnt], t);
	    b->bound_cnt = active->bound_cnt + 1;

	    rval = mpq_ILLutil_priority_insert (minf->que, (void *) b, &lpval,
		&(b->handle));
	    ILL_CLEANUP_IF (rval);

	    mpq_put_bbnode (minf, b);
	    minf->activenodes++;
	}
    }
    minf->totalnodes++;

    if (bdir == 'D') {
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'U', oldt);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'L', oldt);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    mpq_EGlpNumClearVar (t);
    mpq_EGlpNumClearVar (lpval);
    mpq_EGlpNumClearVar (oldt);
    return rval;
}

static int mpq_fix_variables (mpq_lpinfo * lp,
      mpq_t * bestval,
      mpq_bbnode * b,
      mpq_t * wupper,
      mpq_t * wlower,
      int *hit)
{
    int rval = 0;
    int i, nnew = 0;
    int nstruct = lp->O->nstruct;
    mpq_t delta, lpval;
    int *new_indx = 0;
    char *new_lu = 0;
    mpq_t *new_bounds = 0;
    mpq_t *dj = 0;
    mpq_EGlpNumInitVar (delta);
    mpq_EGlpNumInitVar (lpval);

    *hit = 0;

    if (mpq_EGlpNumIsLess (*bestval, mpq_ILL_MAXDOUBLE)) {
	rval = mpq_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);
	/* delta = bestval - lpval + mpq_ILL_INTTOL; */
	mpq_EGlpNumCopy (delta, *bestval);
	mpq_EGlpNumSubTo (delta, lpval);
	mpq_EGlpNumAddTo (delta, mpq_ILL_INTTOL);

	ILL_SAFE_MALLOC (new_indx, nstruct, int);
	ILL_SAFE_MALLOC (new_lu, nstruct, char);
	dj = mpq_EGlpNumAllocArray (nstruct);
	new_bounds = mpq_EGlpNumAllocArray (nstruct);

	rval = mpq_ILLlib_solution (lp, 0, 0, 0, 0, 0, dj);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < nstruct; i++) {
	    if (lp->O->intmarker[i]) {
		if (mpq_EGlpNumIsNeqq (wlower[i], wupper[i])) {
		    if (mpq_EGlpNumIsLess (delta, dj[i])) {
			mpq_EGlpNumSubTo (wupper[i], mpq_oneLpNum);
			rval = mpq_ILLlib_chgbnd (lp, i, 'U', wupper[i]);
			ILL_CLEANUP_IF (rval);
			new_indx[nnew] = i;
			new_lu[nnew] = 'U';
			mpq_EGlpNumCopy (new_bounds[nnew], wupper[i]);
			nnew++;
		    }
		    /* if (-dj[i] > delta) */
		    mpq_EGlpNumSign (delta);
		    if (mpq_EGlpNumIsLess (delta, dj[i])) {
			mpq_EGlpNumAddTo (wlower[i], mpq_oneLpNum);
			rval = mpq_ILLlib_chgbnd (lp, i, 'L', wlower[i]);
			ILL_CLEANUP_IF (rval);
			new_indx[nnew] = i;
			new_lu[nnew] = 'L';
			mpq_EGlpNumCopy (new_bounds[nnew], wlower[i]);
			nnew++;
		    }
		    mpq_EGlpNumSign (delta);
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
	    mpq_EGlpNumReallocArray (&(b->bounds), b->bound_cnt + nnew);
	    for (i = 0; i < nnew; i++) {
		b->bound_indx[b->bound_cnt + i] = new_indx[i];
		b->lu[b->bound_cnt + i] = new_lu[i];
		mpq_EGlpNumCopy (b->bounds[b->bound_cnt + i], new_bounds[i]);
	    }
	    b->bound_cnt += nnew;
	}
    }
    *hit = nnew;

CLEANUP:

    ILL_IFFREE (new_indx, int);
    ILL_IFFREE (new_lu, char);
    mpq_EGlpNumFreeArray (dj);
    mpq_EGlpNumFreeArray (new_bounds);
    mpq_EGlpNumClearVar (delta);
    mpq_EGlpNumClearVar (lpval);
    return rval;
}

static void mpq_best_bbnode (mpq_mipinfo * minf,
      mpq_bbnode ** best)
{
#if 0
    mpq_bbnode *b;
    double bestval = mpq_ILL_MAXDOUBLE;

    for (b = minf->head_bbnode.next; b; b = b->next) {
	if (b->bound < bestval) {
	    *best = b;
	    bestval = b->bound;
	}
    }
#endif

    mpq_t val;
    mpq_EGlpNumInitVar (val);
    mpq_ILLutil_priority_deletemin (minf->que, &val, (void **) best);
    mpq_EGlpNumClearVar (val);
}

static void mpq_put_bbnode (mpq_mipinfo * minf,
      mpq_bbnode * b)
{
    b->next = minf->head_bbnode.next;
    b->prev = &(minf->head_bbnode);
    if (b->next)
	b->next->prev = b;
    minf->head_bbnode.next = b;
}

static void mpq_remove_bbnode (mpq_bbnode * b)
{
    b->prev->next = b->next;
    if (b->next)
	b->next->prev = b->prev;
}

static int mpq_find_branch (mpq_mipinfo * minf,
      mpq_t * x,
      mpq_t * lpval,
      int *bvar)
{
    mpq_lpinfo *lp = minf->lp;
    int rval = 0;

    switch (minf->branching_rule) {
    case mpq_PENALTYBRANCH:
	rval = mpq_find_penalty_branch (lp, minf->pinf, x, minf->downpen,
	    minf->uppen, lpval, bvar);
	ILL_CLEANUP_IF (rval);
	break;
    case mpq_FIRSTBRANCH:
	mpq_find_first_branch (lp, x, bvar);
	break;
    case mpq_MIDDLEBRANCH:
	mpq_find_middle_branch (lp, x, bvar);
	break;
    case mpq_STRONGBRANCH:
	rval = mpq_find_strong_branch (lp, minf->pinf, x, bvar);
	ILL_CLEANUP_IF (rval);
	break;
    default:
	fprintf (stderr, "Unknown branching rule.\n");
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_find_branch");
}

static void mpq_find_first_branch (mpq_lpinfo * lp,
      mpq_t * x,
      int *bvar)
{
    int i, ibest = -1;
    mpq_ILLlpdata *qslp = lp->O;
    mpq_t t;
    mpq_EGlpNumInitVar (t);

    for (i = 0; i < qslp->nstruct; i++) {
	if (qslp->intmarker[i]) {
	    /* t = mpq_ILLutil_our_frac (x[i]); */
	    mpq_EGlpNumFloor (t, x[i]);
	    mpq_EGlpNumSubTo (t, x[i]);
	    mpq_EGlpNumSign (t);
	    if ((mpq_EGlpNumIsNeqZero (t, mpq_ILL_INTTOL)) &&
		(mpq_EGlpNumIsNeq (t, mpq_oneLpNum, mpq_ILL_INTTOL))) {
		ibest = i;
		break;
	    }
	}
    }
    *bvar = ibest;
    mpq_EGlpNumClearVar (t);
}

static void mpq_find_middle_branch (mpq_lpinfo * lp,
      mpq_t * x,
      int *bvar)
{
    int i, ibest = -1;
    mpq_t t, tbest;
    mpq_ILLlpdata *qlp = lp->O;
    mpq_EGlpNumInitVar (t);
    mpq_EGlpNumInitVar (tbest);
    mpq_EGlpNumSet (tbest, 0.5);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = mpq_ILLutil_our_frac (x[i]) - 0.5; if (t < 0.0) t = -t; */
	    mpq_EGlpNumFloor (t, x[i]);
	    mpq_EGlpNumMultUiTo (t, 2);
	    mpq_EGlpNumSubTo (t, mpq_oneLpNum);
	    mpq_EGlpNumDivUiTo (t, 2);
	    if (mpq_EGlpNumIsLess (t, mpq_zeroLpNum))
		mpq_EGlpNumSign (t);
	    /* if (t < tbest) */
	    if (mpq_EGlpNumIsLess (t, tbest)) {
		mpq_EGlpNumCopy (tbest, t);
		ibest = i;
	    }
	}
    }

    /* if (tbest < (0.5 - mpq_ILL_INTTOL)) */
    mpq_EGlpNumAddTo (tbest, mpq_ILL_INTTOL);
    if (mpq_EGlpNumIsLessDbl (tbest, 0.5)) {
	*bvar = ibest;
    } else {
	*bvar = -1;
    }
    mpq_EGlpNumClearVar (t);
    mpq_EGlpNumClearVar (tbest);
}

static int mpq_find_penalty_branch (mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_t * x,
      mpq_t * downpen,
      mpq_t * uppen,
      mpq_t * lpval,
      int *bvar)
{
    int rval = 0;
    int i, k, ibest = -1, ncand = 0, nneed = 0;
    mpq_ILLlpdata *qslp = lp->O;
    int *candidatelist = 0;
    int *needlist = 0;
    mpq_t *fval = 0;
    mpq_t *xlist = 0;
    mpq_t *newdown = 0;
    mpq_t *newup = 0;
    mpq_t a, t, tbest;
    mpq_EGlpNumInitVar (a);
    mpq_EGlpNumInitVar (t);
    mpq_EGlpNumInitVar (tbest);
    mpq_EGlpNumCopy (tbest, mpq_ILL_MINDOUBLE);

    ILL_SAFE_MALLOC (candidatelist, qslp->nstruct, int);
    ILL_SAFE_MALLOC (needlist, qslp->nstruct, int);
    fval = mpq_EGlpNumAllocArray (qslp->nstruct);
    xlist = mpq_EGlpNumAllocArray (qslp->nstruct);
    for (i = 0; i < qslp->nstruct; i++) {
	if (qslp->intmarker[i]) {
	    /* fval[i] = x[i] - floor(x[i]); */
	    mpq_EGlpNumFloor (fval[i], x[i]);
	    mpq_EGlpNumSubTo (fval[i], x[i]);
	    mpq_EGlpNumSign (fval[i]);
	    if ((mpq_EGlpNumIsNeqZero (fval[i], mpq_ILL_INTTOL)) &&
		(mpq_EGlpNumIsNeq (fval[i], mpq_oneLpNum, mpq_ILL_INTTOL))) {
		candidatelist[ncand++] = i;
		/* if (downpen[i] == -1.0) */
		mpq_EGlpNumSign (downpen[i]);
		if (mpq_EGlpNumIsEqqual (downpen[i], mpq_oneLpNum)) {
		    mpq_EGlpNumCopy (xlist[nneed], x[i]);
		    needlist[nneed++] = i;
		}
		mpq_EGlpNumSign (downpen[i]);
	    }
	}
    }

    if (nneed > 0) {
	newdown = mpq_EGlpNumAllocArray (nneed);
	newup = mpq_EGlpNumAllocArray (nneed);
	rval = mpq_ILLlib_strongbranch (lp, pinf, needlist, nneed,
	    0, newdown, newup,
	    5 * mpq_STRONG_PIVOTS, mpq_ILL_MAXDOUBLE);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < nneed; i++) {
	    k = needlist[i];
	    /* uppen[k] = (newup[i] - lpval) / (1.0 - fval[k]); */
	    mpq_EGlpNumCopyDiff (uppen[k], newup[i], *lpval);
	    mpq_EGlpNumCopyDiff (downpen[k], mpq_oneLpNum, fval[k]);
	    mpq_EGlpNumDivTo (uppen[k], downpen[k]);
	    /* downpen[k] = (newdown[i] - lpval) / fval[k]; */
	    mpq_EGlpNumCopyDiffRatio (downpen[k], newdown[i], *lpval, fval[k]);

	}
    }
    for (i = 0; i < ncand; i++) {
	k = candidatelist[i];
	/* t = mpq_ILL_BRANCH_PENALTY_VAL (downpen[k], uppen[k], fval[k]); */
	mpq_EGlpNumCopy (t, downpen[k]);
	mpq_EGlpNumMultTo (t, fval[k]);
	mpq_EGlpNumCopyDiff (a, mpq_oneLpNum, fval[k]);
	mpq_EGlpNumMultTo (a, uppen[k]);
	if (mpq_EGlpNumIsLess (t, a)) {
	    mpq_EGlpNumMultUiTo (t, mpq_ILL_BRANCH_PENALTY_WEIGHT);
	    mpq_EGlpNumAddTo (t, a);
	} else {
	    mpq_EGlpNumMultUiTo (a, mpq_ILL_BRANCH_PENALTY_WEIGHT);
	    mpq_EGlpNumAddTo (t, a);
	}
	mpq_EGlpNumDivUiTo (t, mpq_ILL_BRANCH_PENALTY_WEIGHT + 1);

	if (mpq_EGlpNumIsLess (tbest, t)) {
	    mpq_EGlpNumCopy (tbest, t);
	    ibest = k;
	}
    }

    *bvar = ibest;

CLEANUP:

    mpq_EGlpNumClearVar (a);
    mpq_EGlpNumClearVar (t);
    mpq_EGlpNumClearVar (tbest);
    mpq_EGlpNumFreeArray (newdown);
    mpq_EGlpNumFreeArray (newup);
    mpq_EGlpNumFreeArray (fval);
    mpq_EGlpNumFreeArray (xlist);
    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (needlist, int);
    return rval;
}

static int mpq_find_strong_branch (mpq_lpinfo * lp,
      mpq_price_info * pinf,
      mpq_t * x,
      int *bvar)
{
    int rval = 0;
    int i, ibest = -1, ncand = 0;
    int maxtrys = mpq_STRONG_CANDIDATES;
    mpq_t t, tbest;
    mpq_ILLlpdata *qlp = lp->O;
    int *candidatelist = 0;
    int *newlist = 0;
    int *perm = 0;
    mpq_t *tval = 0;
    mpq_t *xlist = 0;
    mpq_t *downpen = 0;
    mpq_t *uppen = 0;
    ILLrandstate rstate;
    mpq_EGlpNumInitVar (t);
    mpq_EGlpNumInitVar (tbest);
    mpq_EGlpNumCopy (tbest, mpq_ILL_MINDOUBLE);

    ILLutil_sprand (999, &rstate);
    ILL_SAFE_MALLOC (candidatelist, qlp->nstruct, int);
    tval = mpq_EGlpNumAllocArray (qlp->nstruct);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = mpq_ILLutil_our_frac (x[i]) - 0.5; if (t < 0.0) t = -t; */
	    mpq_EGlpNumFloor (t, x[i]);
	    mpq_EGlpNumSubTo (t, x[i]);
	    mpq_EGlpNumSign (t);
	    mpq_EGlpNumMultUiTo (t, 2);
	    mpq_EGlpNumSubTo (t, mpq_oneLpNum);
	    if (mpq_EGlpNumIsLess (t, mpq_zeroLpNum))
		mpq_EGlpNumSign (t);
	    /* if (t < (0.5 - mpq_ILL_INTTOL)) */
	    if (mpq_EGlpNumIsNeq (t, mpq_oneLpNum, mpq_ILL_INTTOL)) {
		candidatelist[ncand] = i;
		mpq_EGlpNumDivUiTo (t, 2);
		mpq_EGlpNumCopy (tval[ncand++], t);
	    }
	}
    }

    if (ncand > 0) {
	if (ncand > maxtrys) {
	    ILL_SAFE_MALLOC (perm, ncand, int);

	    for (i = 0; i < ncand; i++) {
		perm[i] = i;
	    }
	    mpq_ILLutil_EGlpNum_rselect (perm, 0, ncand - 1, maxtrys, tval, &rstate);

	    ILL_SAFE_MALLOC (newlist, maxtrys, int);

	    for (i = 0; i < maxtrys; i++) {
		newlist[i] = candidatelist[perm[i]];
	    }
	    ILL_IFFREE (candidatelist, int);
	    candidatelist = newlist;
	    newlist = 0;
	    ncand = maxtrys;
	}
	downpen = mpq_EGlpNumAllocArray (ncand);
	uppen = mpq_EGlpNumAllocArray (ncand);
	xlist = mpq_EGlpNumAllocArray (ncand);

	for (i = 0; i < ncand; i++) {
	    mpq_EGlpNumCopy (xlist[i], x[candidatelist[i]]);
	}

	rval = mpq_ILLlib_strongbranch (lp, pinf, candidatelist, ncand,
	    0, downpen, uppen, mpq_STRONG_PIVOTS,
	    mpq_ILL_MAXDOUBLE);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < ncand; i++) {
	    /* t = mpq_ILL_BRANCH_STRONG_VAL (downpen[i], uppen[i]); */
	    if (mpq_EGlpNumIsLess (downpen[i], uppen[i])) {
		mpq_EGlpNumCopy (t, downpen[i]);
		mpq_EGlpNumMultUiTo (t, mpq_ILL_BRANCH_STRONG_WEIGHT);
		mpq_EGlpNumAddTo (t, uppen[i]);
	    } else {
		mpq_EGlpNumCopy (t, uppen[i]);
		mpq_EGlpNumMultUiTo (t, mpq_ILL_BRANCH_STRONG_WEIGHT);
		mpq_EGlpNumAddTo (t, downpen[i]);
	    }
	    mpq_EGlpNumDivUiTo (t, mpq_ILL_BRANCH_STRONG_WEIGHT + 1);
	    if (mpq_EGlpNumIsLess (tbest, t)) {
		mpq_EGlpNumCopy (tbest, t);
		ibest = candidatelist[i];
	    }
	}
    }
    *bvar = ibest;


CLEANUP:

    mpq_EGlpNumClearVar (t);
    mpq_EGlpNumClearVar (tbest);
    mpq_EGlpNumFreeArray (tval);
    mpq_EGlpNumFreeArray (xlist);
    mpq_EGlpNumFreeArray (uppen);
    mpq_EGlpNumFreeArray (downpen);
    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (newlist, int);
    ILL_IFFREE (perm, int);

    ILL_RETURN (rval, "mpq_find_strong_branch");
}

static void mpq_check_integral (mpq_lpinfo * lp,
      mpq_t * x,
      int *yesno)
{
    int i;
    mpq_t t;
    mpq_ILLlpdata *qlp = lp->O;
    mpq_EGlpNumInitVar (t);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = mpq_ILLutil_our_frac (x[i]); */
	    mpq_EGlpNumFloor (t, x[i]);
	    mpq_EGlpNumSubTo (t, x[i]);
	    mpq_EGlpNumSign (t);
	    /* if (t > mpq_ILL_INTTOL && t < 1.0 - mpq_ILL_INTTOL) */
	    if ((mpq_EGlpNumIsNeqZero (t, mpq_ILL_INTTOL)) &&
		(mpq_EGlpNumIsNeq (t, mpq_oneLpNum, mpq_ILL_INTTOL))) {
		*yesno = 0;
		mpq_EGlpNumClearVar (t);
		return;
	    }
	}
    }

    *yesno = 1;
    mpq_EGlpNumClearVar (t);
}

static int mpq_plunge (mpq_mipinfo * minf)
{
    int rval = 0;
    int i, status;
    mpq_lpinfo *lp = minf->lp;
    mpq_ILLlpdata *qlp = minf->lp->O;
    mpq_t *oldlower = 0;
    mpq_t *oldupper = 0;

    if (minf->watch) {
	printf ("Plunging ...\n");
	fflush (stdout);
    }
    oldlower = mpq_EGlpNumAllocArray (qlp->nstruct);
    oldupper = mpq_EGlpNumAllocArray (qlp->nstruct);

    for (i = 0; i < qlp->nstruct; i++) {
	mpq_EGlpNumCopy (oldlower[i], minf->lower[i]);
	mpq_EGlpNumCopy (oldupper[i], minf->upper[i]);
    }

    rval = mpq_plunge_work (minf, 0);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
	rval = mpq_ILLlib_chgbnd (lp, i, 'L', oldlower[i]);
	ILL_CLEANUP_IF (rval);
	rval = mpq_ILLlib_chgbnd (lp, i, 'U', oldupper[i]);
	ILL_CLEANUP_IF (rval);
	mpq_EGlpNumCopy (minf->lower[i], oldlower[i]);
	mpq_EGlpNumCopy (minf->upper[i], oldupper[i]);
    }

    rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    mpq_EGlpNumFreeArray (oldlower);
    mpq_EGlpNumFreeArray (oldupper);

    ILL_RETURN (rval, "mpq_plunge");
}

static int mpq_plunge_work (mpq_mipinfo * minf,
      int depth)
{
    int rval = 0;
    int bvar, status, count;
    mpq_t lpval, val0, val1, int_tol;
    mpq_lpinfo *lp = minf->lp;
    mpq_EGlpNumInitVar (lpval);
    mpq_EGlpNumInitVar (val0);
    mpq_EGlpNumInitVar (val1);
    mpq_EGlpNumInitVar (int_tol);
    mpq_EGlpNumSet (int_tol, 0.001);

    rval = mpq_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = mpq_round_variables (minf, &count, &int_tol /* 0.001 */ );
    ILL_CLEANUP_IF (rval);
    if (count) {
	rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	if (status != QS_LP_OPTIMAL) {
	    goto CLEANUP;
	}
	rval = mpq_ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);
    }
    mpq_find_middle_branch (lp, minf->x, &bvar);
    if (bvar == -1) {
	rval = mpq_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);

	if (mpq_EGlpNumIsLess (lpval, minf->value)) {
	    printf ("Plunge Integral Solution: %.6f (Depth: %d)\n",
		mpq_EGlpNumToLf (lpval), depth);
	    fflush (stdout);

	    mpq_EGlpNumCopy (minf->value, lpval);
	    mpq_EGlpNumCopyDiff (minf->objectivebound, lpval, mpq_ILL_INTTOL);
	    mpq_copy_x (lp->O->nstruct, minf->x, minf->bestx);
	}
	goto CLEANUP;
    }
    mpq_EGlpNumOne (minf->lower[bvar]);
    rval = mpq_ILLlib_chgbnd (lp, bvar, 'L', mpq_oneLpNum);
    ILL_CLEANUP_IF (rval);
    rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the mpq_plunge LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
	mpq_EGlpNumCopy (val1, mpq_ILL_MAXDOUBLE);
    } else if (status == QS_LP_OPTIMAL) {
	rval = mpq_ILLlib_objval (lp, 0, &val1);
	ILL_CLEANUP_IF (rval);
    } else {
	ILL_CLEANUP;
    }

    rval = mpq_ILLlib_chgbnd (lp, bvar, 'L', mpq_zeroLpNum);
    ILL_CLEANUP_IF (rval);
    mpq_EGlpNumZero (minf->lower[bvar]);

    mpq_EGlpNumZero (minf->upper[bvar]);
    rval = mpq_ILLlib_chgbnd (lp, bvar, 'U', mpq_zeroLpNum);
    ILL_CLEANUP_IF (rval);
    rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the mpq_plunge LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
	mpq_EGlpNumCopy (val0, mpq_ILL_MAXDOUBLE);
    } else if (status == QS_LP_OPTIMAL) {
	rval = mpq_ILLlib_objval (lp, 0, &val0);
	ILL_CLEANUP_IF (rval);
    } else {
	ILL_CLEANUP;
    }

    rval = mpq_ILLlib_chgbnd (lp, bvar, 'U', mpq_oneLpNum);
    ILL_CLEANUP_IF (rval);
    mpq_EGlpNumCopy (minf->upper[bvar], mpq_oneLpNum);

    if (mpq_EGlpNumIsEqqual (val0, mpq_ILL_MAXDOUBLE) &&
	mpq_EGlpNumIsEqqual (val1, mpq_ILL_MAXDOUBLE)) {
	ILL_CLEANUP;
    }
    if (mpq_EGlpNumIsLess (val0, val1)) {
	mpq_EGlpNumZero (minf->upper[bvar]);
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'U', mpq_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	rval = mpq_plunge_work (minf, depth + 1);
	ILL_CLEANUP_IF (rval);
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'U', mpq_oneLpNum);
	ILL_CLEANUP_IF (rval);
	mpq_EGlpNumOne (minf->upper[bvar]);
    } else {
	mpq_EGlpNumOne (minf->lower[bvar]);
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'L', mpq_oneLpNum);
	ILL_CLEANUP_IF (rval);
	rval = mpq_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	rval = mpq_plunge_work (minf, depth + 1);
	ILL_CLEANUP_IF (rval);
	rval = mpq_ILLlib_chgbnd (lp, bvar, 'L', mpq_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	mpq_EGlpNumZero (minf->lower[bvar]);
    }

CLEANUP:

    mpq_EGlpNumClearVar (lpval);
    mpq_EGlpNumClearVar (val0);
    mpq_EGlpNumClearVar (val1);
    mpq_EGlpNumClearVar (int_tol);
    ILL_RETURN (rval, "mpq_plunge_work");
}

static int mpq_round_variables (mpq_mipinfo * minf,
      int *count,
      mpq_t * tol)
{
    int rval = 0;
    int i, hit = 0;
    mpq_lpinfo *lp = minf->lp;
    mpq_ILLlpdata *qlp = lp->O;

    *count = 0;

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    if (mpq_EGlpNumIsNeqq (minf->lower[i], minf->upper[i])) {
		if (mpq_EGlpNumIsLess (minf->x[i], *tol)) {
		    mpq_EGlpNumZero (minf->upper[i]);
		    rval = mpq_ILLlib_chgbnd (lp, i, 'U', mpq_zeroLpNum);
		    ILL_CLEANUP_IF (rval);
		    hit++;
		} else if (mpq_EGlpNumIsEqual (minf->x[i], mpq_oneLpNum, *tol)) {
		    mpq_EGlpNumOne (minf->lower[i]);
		    rval = mpq_ILLlib_chgbnd (lp, i, 'L', mpq_oneLpNum);
		    ILL_CLEANUP_IF (rval);
		    hit++;
		}
	    }
	}
    }
    *count = hit;

CLEANUP:

    ILL_RETURN (rval, "mpq_round_variables");
}

static void mpq_copy_x (int nstruct,
      mpq_t * from_x,
      mpq_t * to_x)
{
    int j;

    for (j = 0; j < nstruct; j++) {
	mpq_EGlpNumCopy (to_x[j], from_x[j]);
    }
}

static void mpq_init_mipinfo (mpq_mipinfo * minf)
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
	minf->branching_rule = /* mpq_MIDDLEBRANCH */ mpq_STRONGBRANCH;
	minf->watch = 1;
	mpq_EGlpNumInitVar (minf->objectivebound);
	mpq_EGlpNumInitVar (minf->value);
	mpq_EGlpNumCopy (minf->objectivebound, mpq_ILL_MAXDOUBLE);
	mpq_EGlpNumCopy (minf->value, mpq_ILL_MAXDOUBLE);
	ILLptrworld_init (&minf->ptrworld);
    }
}

static void mpq_free_mipinfo (mpq_mipinfo * minf)
{
    int total, onlist;

    if (minf) {
	mpq_EGlpNumFreeArray (minf->downpen);
	mpq_EGlpNumFreeArray (minf->uppen);
	mpq_EGlpNumFreeArray (minf->x);
	mpq_EGlpNumFreeArray (minf->bestx);
	mpq_EGlpNumFreeArray (minf->lower);
	mpq_EGlpNumFreeArray (minf->upper);
	bbnode_listfree (&minf->ptrworld, minf->head_bbnode.next);
	if (bbnode_check_leaks (&minf->ptrworld, &total, &onlist)) {
	    fprintf (stderr, "WARNING: %d outstanding bbnodes\n", total - onlist);
	}
	ILLptrworld_delete (&minf->ptrworld);
	mpq_EGlpNumClearVar ((minf->objectivebound));
	mpq_EGlpNumClearVar ((minf->value));
	memset (minf, 0, sizeof (mpq_mipinfo));
	/* mpq_init_mipinfo (minf); */
    }
}

static void mpq_init_bbnode (mpq_bbnode * b)
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
	mpq_EGlpNumInitVar ((b->bound));
	mpq_EGlpNumCopy (b->bound, mpq_ILL_MINDOUBLE);
    }
}

static void mpq_free_bbnode (mpq_bbnode * b)
{
    if (b) {
	mpq_EGlpNumFreeArray (b->rownorms);
	mpq_EGlpNumFreeArray (b->bounds);
	ILL_IFFREE (b->cstat, char);
	ILL_IFFREE (b->rstat, char);
	ILL_IFFREE (b->bound_indx, int);
	ILL_IFFREE (b->lu, char);
	mpq_EGlpNumClearVar ((b->bound));
	memset (b, 0, sizeof (mpq_bbnode));
    }
}
