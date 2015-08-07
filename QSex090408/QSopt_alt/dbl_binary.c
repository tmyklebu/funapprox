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

/* RCS_INFO = "$RCSfile: dbl_binary.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* Simple MIP Code to test LP Solver                    */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int dbl_ILLmip_bfs (dbl_lpinfo *lp, double *val, double *x)                   */
/* */
/* NOTES                                                                   */
/* */
/* */
/****************************************************************************/

#include "econfig.h"
#include "dbl_priority.h"
#include "dbl_sortrus.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdata.h"
#include "dbl_lpdefs.h"
#include "dbl_simplex.h"
#include "dbl_binary.h"
#include "dbl_price.h"
#include "dbl_lib.h"
#include "dbl_qstruct.h"
#include "dbl_qsopt.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/* #define  dbl_ILL_INTTOL (0.000001) */
#define dbl_ILL_INTTOL dbl_PFEAS_TOLER

#define  dbl_STRONG_PIVOTS     (50)
#define  dbl_STRONG_CANDIDATES (10)

#define dbl_ILL_BRANCH_STRONG_WEIGHT (10)
#define dbl_ILL_BRANCH_STRONG_VAL(v0,v1)                                 \
    (((v0) < (v1) ? (dbl_ILL_BRANCH_STRONG_WEIGHT * (v0) + (v1))         \
                  : (dbl_ILL_BRANCH_STRONG_WEIGHT * (v1) + (v0)))        \
                    / (dbl_ILL_BRANCH_STRONG_WEIGHT + 1.0))

#define dbl_ILL_BRANCH_PENALTY_WEIGHT (2)
#define dbl_ILL_BRANCH_PENALTY_VAL(v0,v1,f)                              \
    (((v0)*(f) < (v1)*(1.0-(f)) ?                                    \
        (dbl_ILL_BRANCH_PENALTY_WEIGHT * (v0)*(f) + (v1)*(1.0-(f)))    \
      : (dbl_ILL_BRANCH_PENALTY_WEIGHT * (v1)*(1.0-(f)) + (v0)*(f)))    \
                    / (dbl_ILL_BRANCH_PENALTY_WEIGHT + 1.0))



#define dbl_FIRSTBRANCH  1
#define dbl_MIDDLEBRANCH 2
#define dbl_STRONGBRANCH 3
#define dbl_PENALTYBRANCH 4


typedef struct dbl_bbnode {
    struct dbl_bbnode *next;
    struct dbl_bbnode *prev;
    int id;
    int depth;
    int handle;
    double bound;
    char *cstat;
    char *rstat;
    double *rownorms;
    int rownorms_size;
    int bound_cnt;
    int *bound_indx;
    char *lu;
    double *bounds;
    int bounds_size;
}
  dbl_bbnode;

typedef struct dbl_mipinfo {
    int branching_rule;
    int watch;
    int depth;
    int totalnodes;
    int activenodes;
    int totalpivots;
    int lastpivots;
    int objsense;
    double objectivebound;
    double value;
    double *downpen;
    double *uppen;
    double *x;
    double *bestx;
    double *orig_lower;
    double *orig_upper;
    double *lower;
    double *upper;
    int nstruct;		/* size of all double arrays */
    dbl_lpinfo *lp;
    dbl_price_info *pinf;
    dbl_bbnode head_bbnode;
    dbl_ILLpriority *que;
    ILLptrworld ptrworld;
}
  dbl_mipinfo;


ILL_PTRWORLD_ROUTINES (dbl_bbnode, bbnodealloc, bbnode_bulkalloc, bbnodefree)
ILL_PTRWORLD_LISTFREE_ROUTINE (dbl_bbnode, bbnode_listfree, bbnodefree)
ILL_PTRWORLD_LEAKS_ROUTINE (dbl_bbnode, bbnode_check_leaks, depth, int)
static void dbl_cleanup_mip (dbl_mipinfo * minf),
  dbl_choose_initial_price (dbl_price_info * pinf),
  dbl_best_bbnode (dbl_mipinfo * minf,
      dbl_bbnode ** best),
  dbl_put_bbnode (dbl_mipinfo * minf,
      dbl_bbnode * b),
  dbl_remove_bbnode (dbl_bbnode * b),
  dbl_find_first_branch (dbl_lpinfo * lp,
      double *x,
      int *bvar),
  dbl_find_middle_branch (dbl_lpinfo * lp,
      double *x,
      int *bvar),
  dbl_check_integral (dbl_lpinfo * lp,
      double *x,
      int *yesno),
  dbl_copy_x (int nstruct,
      double *from_x,
      double *to_x),
  dbl_init_mipinfo (dbl_mipinfo * minf),
  dbl_free_mipinfo (dbl_mipinfo * minf),
  dbl_init_bbnode (dbl_bbnode * b),
  dbl_free_bbnode (dbl_bbnode * b);

static int dbl_startup_mip (dbl_mipinfo * minf,
      dbl_lpinfo * lp,
      dbl_price_info * pinf,
      double *lpval),
  dbl_run_bfs (dbl_mipinfo * minf),
  dbl_process_bfs_bbnode (dbl_mipinfo * minf,
      dbl_bbnode * b),
  dbl_child_work (dbl_mipinfo * minf,
      dbl_bbnode * active,
      int bvar,
      int bdir,
      double *cval,
      int *cp),
  dbl_fix_variables (dbl_lpinfo * lp,
      double *bestval,
      dbl_bbnode * b,
      double *wupper,
      double *wlower,
      int *hit),
  dbl_find_branch (dbl_mipinfo * minf,
      double *x,
      double *lpval,
      int *bvar),
  dbl_find_penalty_branch (dbl_lpinfo * lp,
      dbl_price_info * pinf,
      double *x,
      double *downpen,
      double *uppen,
      double *lpval,
      int *bvar),
  dbl_find_strong_branch (dbl_lpinfo * lp,
      dbl_price_info * pinf,
      double *x,
      int *bvar),
  dbl_plunge (dbl_mipinfo * minf),
  dbl_plunge_work (dbl_mipinfo * minf,
      int depth),
  dbl_round_variables (dbl_mipinfo * minf,
      int *count,
      double *tol);


static void dbl_choose_initial_price (dbl_price_info * pinf)
{
    pinf->pI_price = QS_PRICE_PSTEEP;
    pinf->pII_price = QS_PRICE_PSTEEP;
    pinf->dI_price = QS_PRICE_DSTEEP;
    pinf->dII_price = QS_PRICE_DSTEEP;
}

int dbl_ILLmip_bfs (dbl_lpinfo * lp,
      double *val,
      double *x)
{
    int tval, rval = 0;
    dbl_price_info pinf;
    dbl_mipinfo minf;
    dbl_bbnode *b;
    double lpval;
    double szeit = ILLutil_zeit ();
    dbl_EGlpNumInitVar (lpval);
    dbl_EGlpNumInitVar (pinf.htrigger);

    dbl_ILLprice_init_pricing_info (&pinf);
    dbl_init_mipinfo (&minf);

    if (!lp) {
	fprintf (stderr, "dbl_ILLmip_bfs called without an LP\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_startup_mip (&minf, lp, &pinf, &lpval);
    ILL_CLEANUP_IF (rval);

    ILL_SAFE_MALLOC (minf.que, 1, dbl_ILLpriority);
    rval = dbl_ILLutil_priority_init (minf.que, lp->O->nstruct + 1);
    ILL_CLEANUP_IF (rval);

    b = bbnodealloc (&minf.ptrworld);
    dbl_init_bbnode (b);
    b->depth = 0;
    b->id = minf.totalnodes++;
    dbl_EGlpNumCopy (b->bound, lpval);
    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
    rval = dbl_ILLlib_getbasis (lp, b->cstat, b->rstat);
    ILL_CLEANUP_IF (rval);

    if (pinf.dII_price == QS_PRICE_DSTEEP) {
	b->rownorms = dbl_EGlpNumAllocArray (lp->nrows);
	tval = dbl_ILLlib_getrownorms (lp, &pinf, b->rownorms);
	if (tval) {
	    printf ("Row norms not available\n");
	    fflush (stdout);
	    dbl_EGlpNumFreeArray (b->rownorms);
	}
    }
    rval = dbl_ILLutil_priority_insert (minf.que, (void *) b, &lpval, &(b->handle));
    ILL_CLEANUP_IF (rval);

    b->prev = &(minf.head_bbnode);
    b->next = 0;
    minf.head_bbnode.next = b;
    minf.activenodes++;

    minf.branching_rule = dbl_PENALTYBRANCH;

    rval = dbl_run_bfs (&minf);
    ILL_CLEANUP_IF (rval);

    printf ("Total Number of Nodes: %d\n", minf.totalnodes);
    printf ("Total Number of Pivots: %d\n", minf.totalpivots);
    printf ("BFS MIP Runing Time: %.2f seconds\n", ILLutil_zeit () - szeit);
    fflush (stdout);

    dbl_EGlpNumCopy (*val, minf.value);
    if (minf.objsense == dbl_ILL_MAX)
	dbl_EGlpNumSign (*val);

    if (x && dbl_EGlpNumIsNeqq (minf.value, dbl_ILL_MAXDOUBLE)) {
	dbl_copy_x (lp->O->nstruct, minf.bestx, x);
    }
CLEANUP:

    if (minf.que) {
	dbl_ILLutil_priority_free (minf.que);
	ILL_IFFREE (minf.que, dbl_ILLpriority);
    }
    dbl_cleanup_mip (&minf);
    dbl_free_mipinfo (&minf);
    dbl_ILLprice_free_pricing_info (&pinf);
    dbl_EGlpNumClearVar (lpval);
    dbl_EGlpNumClearVar (pinf.htrigger);
    ILL_RETURN (rval, "dbl_ILLmip_bfs");
}

static int dbl_startup_mip (dbl_mipinfo * minf,
      dbl_lpinfo * lp,
      dbl_price_info * pinf,
      double *lpval)
{
    int rval = 0;
    int i, col, status, intcount = 0;
    double val;
    dbl_ILLlpdata *qlp;
    dbl_EGlpNumInitVar (val);

    dbl_choose_initial_price (pinf);

    qlp = lp->O;

    rval = dbl_ILLlib_optimize (lp, 0, pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += dbl_ILLlib_iter (lp);

    rval = dbl_ILLlib_objval (lp, 0, &val);
    ILL_CLEANUP_IF (rval);

    printf ("LP Value: %.6f\n", dbl_EGlpNumToLf (val));
    fflush (stdout);
    if (lpval)
	dbl_EGlpNumCopy (*lpval, val);

    if (qlp->intmarker) {
	for (i = 0; i < qlp->nstruct; i++) {
	    if (qlp->intmarker[i]) {
		col = qlp->structmap[i];
		intcount++;
		if (dbl_EGlpNumIsEqqual (qlp->lower[col], dbl_ILL_MINDOUBLE)
		    || dbl_EGlpNumIsEqqual (qlp->upper[col], dbl_ILL_MAXDOUBLE)) {
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
	dbl_ILLlp_sinfo_free (qlp->sinfo);
	ILL_IFFREE (qlp->sinfo, dbl_ILLlp_sinfo);
    }
    minf->lp = lp;
    minf->pinf = pinf;
    minf->objsense = qlp->objsense;
    if (qlp->objsense == dbl_ILL_MAX) {	/* MIP codes work with min */
	for (i = 0; i < lp->ncols; i++) {
	    dbl_EGlpNumCopyNeg (qlp->obj[i], qlp->obj[i]);
	}
	qlp->objsense = dbl_ILL_MIN;
    }
    minf->x = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->bestx = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->lower = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->upper = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->orig_lower = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->orig_upper = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->downpen = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->uppen = dbl_EGlpNumAllocArray (qlp->nstruct);
    minf->nstruct = qlp->nstruct;

    rval = dbl_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
	dbl_EGlpNumCopy (minf->lower[i], qlp->lower[i]);
	dbl_EGlpNumCopy (minf->upper[i], qlp->upper[i]);
	dbl_EGlpNumCopy (minf->orig_lower[i], qlp->lower[i]);
	dbl_EGlpNumCopy (minf->orig_upper[i], qlp->upper[i]);
	dbl_EGlpNumOne (minf->downpen[i]);
	dbl_EGlpNumOne (minf->uppen[i]);
	dbl_EGlpNumSign (minf->downpen[i]);
	dbl_EGlpNumSign (minf->uppen[i]);
    }


CLEANUP:

    dbl_EGlpNumClearVar (val);
    ILL_RETURN (rval, "dbl_startup_mip");
}

static void dbl_cleanup_mip (dbl_mipinfo * minf)
{
    int i;
    dbl_ILLlpdata *qslp;

    if (minf && minf->lp) {
	qslp = minf->lp->O;
	if (minf->objsense == dbl_ILL_MAX) {
	    for (i = 0; i < minf->lp->ncols; i++) {
		dbl_EGlpNumSign (qslp->obj[i]);
	    }
	    qslp->objsense = dbl_ILL_MIN;
	}
    }
}

static int dbl_run_bfs (dbl_mipinfo * minf)
{
    int rval = 0;
    dbl_bbnode *b;

    while (minf->head_bbnode.next) {
	dbl_best_bbnode (minf, &b);
	rval = dbl_process_bfs_bbnode (minf, b);
	ILL_CLEANUP_IF (rval);
	dbl_remove_bbnode (b);
	dbl_free_bbnode (b);
	bbnodefree (&minf->ptrworld, b);
	minf->activenodes--;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_run_bfs");
}

static int dbl_process_bfs_bbnode (dbl_mipinfo * minf,
      dbl_bbnode * active)
{
    dbl_lpinfo *lp = minf->lp;
    dbl_ILLlp_basis B;
    int status, bvar = 0;
    int i, j, hit, dnp = 0, upp = 0;
    int nstruct = lp->O->nstruct;
    double t, lpval, dnval, upval;
    double *wupper = 0;
    double *wlower = 0;
    int rval = 0;
    dbl_EGlpNumInitVar (t);
    dbl_EGlpNumInitVar (lpval);
    dbl_EGlpNumInitVar (dnval);
    dbl_EGlpNumInitVar (upval);

    dbl_ILLlp_basis_init (&B);

    if (minf->watch > 1) {
	printf ("Node %4d: %.3f", active->id, dbl_EGlpNumToLf (active->bound));
	if (dbl_EGlpNumIsNeqq (minf->value, dbl_ILL_MAXDOUBLE))
	    printf (" %.3f", dbl_EGlpNumToLf (minf->value));
	else
	    printf ("  None");
	printf (", Active %d ", minf->activenodes);
	fflush (stdout);
    } else if (minf->watch == 1) {
	if (minf->lastpivots > 1000) {
	    minf->lastpivots = 0;
	    printf ("Pivots %d, Active Nodes %d, Bound %.3f, Soln ",
		minf->totalpivots, minf->activenodes,
		dbl_EGlpNumToLf (active->bound));
	    if (!dbl_EGlpNumIsLess (minf->value, dbl_ILL_MAXDOUBLE))
		printf ("%.3f", dbl_EGlpNumToLf (minf->value));
	    else
		printf ("None\n");
	}
    }
    if (dbl_EGlpNumIsLeq (minf->objectivebound, active->bound)) {
	if (minf->watch > 1) {
	    printf ("  Node can be purged\n");
	    fflush (stdout);
	}
	goto CLEANUP;
    }
    /* Set the LP bounds for the dbl_node. */

    wlower = dbl_EGlpNumAllocArray (nstruct);
    wupper = dbl_EGlpNumAllocArray (nstruct);

    for (i = 0; i < nstruct; i++) {
	dbl_EGlpNumCopy (wlower[i], minf->orig_lower[i]);
	dbl_EGlpNumCopy (wupper[i], minf->orig_upper[i]);
    }
    for (i = 0; i < active->bound_cnt; i++) {
	j = active->bound_indx[i];
	if (active->lu[i] == 'L')
	    dbl_EGlpNumCopy (wlower[j], active->bounds[i]);
	else
	    dbl_EGlpNumCopy (wupper[j], active->bounds[i]);
    }

    if (active->bound_cnt > 0) {
	rval = dbl_ILLlib_chgbnds (lp, active->bound_cnt, active->bound_indx,
	    active->lu, active->bounds);
	ILL_CLEANUP_IF (rval);
    }
    /* Solve the LP. */

    rval = dbl_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
	active->rstat);
    ILL_CLEANUP_IF (rval);
    if (active->rownorms) {
	B.rownorms = dbl_EGlpNumAllocArray (lp->nrows);
	for (i = 0; i < lp->nrows; i++) {
	    dbl_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
	}
    }
    rval = dbl_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += dbl_ILLlib_iter (lp);
    minf->lastpivots += dbl_ILLlib_iter (lp);

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
	    dbl_EGlpNumCopy (minf->lower[i], wlower[i]);
	    dbl_EGlpNumCopy (minf->upper[i], wupper[i]);
	}
	rval = dbl_plunge (minf);
	ILL_CLEANUP_IF (rval);
    }
    /* Fix variables. */

    if (dbl_EGlpNumIsLess (minf->value, dbl_ILL_MAXDOUBLE)) {
	rval = dbl_fix_variables (lp, &(minf->value), active, wupper, wlower, &hit);
	ILL_CLEANUP_IF (rval);

	if (hit) {
	    rval = dbl_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0);
	    ILL_CLEANUP_IF (rval);

	    minf->totalpivots += dbl_ILLlib_iter (lp);
	    minf->lastpivots += dbl_ILLlib_iter (lp);

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

    rval = dbl_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_objval (lp, 0, &lpval);
    ILL_CLEANUP_IF (rval);

    rval = dbl_find_branch (minf, minf->x, &lpval, &bvar);
    ILL_CLEANUP_IF (rval);

    if (bvar == -1) {
	printf ("Found integral solution: %f\n", dbl_EGlpNumToLf (lpval));
	if (dbl_EGlpNumIsLess (lpval, minf->value)) {
	    dbl_EGlpNumCopy (minf->value, lpval);
	    dbl_EGlpNumCopy (minf->objectivebound, lpval);
	    dbl_EGlpNumSubTo (minf->objectivebound, dbl_ILL_INTTOL);
	    dbl_copy_x (nstruct, minf->x, minf->bestx);
	}
    } else {
	/* Create down child */

	rval = dbl_child_work (minf, active, bvar, 'D', &dnval, &dnp);
	ILL_CLEANUP_IF (rval);

	/* Restore parent basis */

	rval = dbl_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
	    active->rstat);
	ILL_CLEANUP_IF (rval);
	if (active->rownorms) {
	    B.rownorms = dbl_EGlpNumAllocArray (lp->nrows);
	    for (i = 0; i < lp->nrows; i++) {
		dbl_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
	    }
	}
	/* Create up child */

	rval = dbl_child_work (minf, active, bvar, 'U', &upval, &upp);
	ILL_CLEANUP_IF (rval);

	if (minf->watch > 1) {
	    if (dbl_EGlpNumIsEqqual (dnval, dbl_ILL_MAXDOUBLE)) {
		printf ("DN->XXX");
	    } else {
		printf ("DN->%.3f%c", dbl_EGlpNumToLf (dnval), dnp ? 'X' : ' ');
	    }
	    if (dbl_EGlpNumIsEqqual (upval, dbl_ILL_MAXDOUBLE)) {
		printf ("UP->XXX\n");
	    } else {
		printf ("UP->%.3f%c\n", dbl_EGlpNumToLf (upval), upp ? 'X' : ' ');
	    }
	    fflush (stdout);
	}
    }

    /* Set the LP bounds back to original values */

    for (i = 0; i < active->bound_cnt; i++) {
	if (active->lu[i] == 'L')
	    dbl_EGlpNumCopy (t, minf->orig_lower[active->bound_indx[i]]);
	else
	    dbl_EGlpNumCopy (t, minf->orig_upper[active->bound_indx[i]]);

	rval = dbl_ILLlib_chgbnd (lp, active->bound_indx[i], active->lu[i], t);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    dbl_EGlpNumFreeArray (wlower);
    dbl_EGlpNumFreeArray (wupper);
    dbl_ILLlp_basis_free (&B);
    dbl_EGlpNumClearVar (t);
    dbl_EGlpNumClearVar (lpval);
    dbl_EGlpNumClearVar (dnval);
    dbl_EGlpNumClearVar (upval);
    ILL_RETURN (rval, "dbl_process_bfs_bbnode");
}

static int dbl_child_work (dbl_mipinfo * minf,
      dbl_bbnode * active,
      int bvar,
      int bdir,
      double *cval,
      int *cp)
{
    int tval, rval = 0;
    int i, status, intsol;
    double t, oldt, lpval;
    double *xi = &(minf->x[bvar]);
    dbl_lpinfo *lp = minf->lp;
    dbl_bbnode *b;
    dbl_EGlpNumInitVar (t);
    dbl_EGlpNumInitVar (lpval);
    dbl_EGlpNumInitVar (oldt);

    *cp = 0;

    if (bdir == 'D') {
	rval = dbl_ILLlib_getbnd (lp, bvar, 'U', &oldt);
	ILL_CLEANUP_IF (rval);
	dbl_EGlpNumFloor (t, *xi);
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'U', t);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = dbl_ILLlib_getbnd (lp, bvar, 'L', &oldt);
	ILL_CLEANUP_IF (rval);
	dbl_EGlpNumCeil (t, *xi);
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'L', t);
	ILL_CLEANUP_IF (rval);
    }

    rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    minf->totalpivots += dbl_ILLlib_iter (lp);
    minf->lastpivots += dbl_ILLlib_iter (lp);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve Child LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    }
    if (status == QS_LP_INFEASIBLE) {
	dbl_EGlpNumCopy (*cval, dbl_ILL_MAXDOUBLE);
	*cp = 1;
    } else {
	rval = dbl_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);
	dbl_EGlpNumCopy (*cval, lpval);

	/* What about the x vector?  Bico - 020531 */

	dbl_check_integral (lp, minf->x, &intsol);
	if (intsol) {
	    if (dbl_EGlpNumIsLess (lpval, minf->value)) {
		printf ("Found integral solution: %f\n", dbl_EGlpNumToLf (lpval));
		dbl_EGlpNumCopy (minf->value, lpval);
		dbl_EGlpNumCopy (minf->objectivebound, lpval);
		dbl_EGlpNumSubTo (minf->objectivebound, dbl_ILL_INTTOL);
		dbl_copy_x (lp->O->nstruct, minf->x, minf->bestx);
	    }
	}
	if (dbl_EGlpNumIsLeq (minf->objectivebound, lpval)) {
	    *cp = 1;
	} else {
	    b = bbnodealloc (&minf->ptrworld);
	    dbl_init_bbnode (b);
	    b->depth = active->depth + 1;
	    b->id = minf->totalnodes;
	    dbl_EGlpNumCopy (b->bound, lpval);
	    ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
	    ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);
	    rval = dbl_ILLlib_getbasis (lp, b->cstat, b->rstat);
	    ILL_CLEANUP_IF (rval);
	    if (minf->pinf->dII_price == QS_PRICE_DSTEEP) {
		b->rownorms = dbl_EGlpNumAllocArray (lp->nrows);
		tval = dbl_ILLlib_getrownorms (lp, minf->pinf, b->rownorms);
		if (tval) {
		    printf ("Row norms not available\n");
		    fflush (stdout);
		    printf ("A\n");
		    exit (1);
		    dbl_EGlpNumFreeArray (b->rownorms);
		}
	    }
	    ILL_SAFE_MALLOC (b->bound_indx, active->bound_cnt + 1, int);
	    ILL_SAFE_MALLOC (b->lu, active->bound_cnt + 1, char);
	    b->bounds = dbl_EGlpNumAllocArray (active->bound_cnt + 1);
	    for (i = 0; i < active->bound_cnt; i++) {
		b->bound_indx[i] = active->bound_indx[i];
		b->lu[i] = active->lu[i];
		dbl_EGlpNumCopy (b->bounds[i], active->bounds[i]);
	    }
	    b->bound_indx[active->bound_cnt] = bvar;
	    if (bdir == 'D')
		b->lu[active->bound_cnt] = 'U';
	    else
		b->lu[active->bound_cnt] = 'L';
	    dbl_EGlpNumCopy (b->bounds[active->bound_cnt], t);
	    b->bound_cnt = active->bound_cnt + 1;

	    rval = dbl_ILLutil_priority_insert (minf->que, (void *) b, &lpval,
		&(b->handle));
	    ILL_CLEANUP_IF (rval);

	    dbl_put_bbnode (minf, b);
	    minf->activenodes++;
	}
    }
    minf->totalnodes++;

    if (bdir == 'D') {
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'U', oldt);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'L', oldt);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    dbl_EGlpNumClearVar (t);
    dbl_EGlpNumClearVar (lpval);
    dbl_EGlpNumClearVar (oldt);
    return rval;
}

static int dbl_fix_variables (dbl_lpinfo * lp,
      double *bestval,
      dbl_bbnode * b,
      double *wupper,
      double *wlower,
      int *hit)
{
    int rval = 0;
    int i, nnew = 0;
    int nstruct = lp->O->nstruct;
    double delta, lpval;
    int *new_indx = 0;
    char *new_lu = 0;
    double *new_bounds = 0;
    double *dj = 0;
    dbl_EGlpNumInitVar (delta);
    dbl_EGlpNumInitVar (lpval);

    *hit = 0;

    if (dbl_EGlpNumIsLess (*bestval, dbl_ILL_MAXDOUBLE)) {
	rval = dbl_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);
	/* delta = bestval - lpval + dbl_ILL_INTTOL; */
	dbl_EGlpNumCopy (delta, *bestval);
	dbl_EGlpNumSubTo (delta, lpval);
	dbl_EGlpNumAddTo (delta, dbl_ILL_INTTOL);

	ILL_SAFE_MALLOC (new_indx, nstruct, int);
	ILL_SAFE_MALLOC (new_lu, nstruct, char);
	dj = dbl_EGlpNumAllocArray (nstruct);
	new_bounds = dbl_EGlpNumAllocArray (nstruct);

	rval = dbl_ILLlib_solution (lp, 0, 0, 0, 0, 0, dj);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < nstruct; i++) {
	    if (lp->O->intmarker[i]) {
		if (dbl_EGlpNumIsNeqq (wlower[i], wupper[i])) {
		    if (dbl_EGlpNumIsLess (delta, dj[i])) {
			dbl_EGlpNumSubTo (wupper[i], dbl_oneLpNum);
			rval = dbl_ILLlib_chgbnd (lp, i, 'U', wupper[i]);
			ILL_CLEANUP_IF (rval);
			new_indx[nnew] = i;
			new_lu[nnew] = 'U';
			dbl_EGlpNumCopy (new_bounds[nnew], wupper[i]);
			nnew++;
		    }
		    /* if (-dj[i] > delta) */
		    dbl_EGlpNumSign (delta);
		    if (dbl_EGlpNumIsLess (delta, dj[i])) {
			dbl_EGlpNumAddTo (wlower[i], dbl_oneLpNum);
			rval = dbl_ILLlib_chgbnd (lp, i, 'L', wlower[i]);
			ILL_CLEANUP_IF (rval);
			new_indx[nnew] = i;
			new_lu[nnew] = 'L';
			dbl_EGlpNumCopy (new_bounds[nnew], wlower[i]);
			nnew++;
		    }
		    dbl_EGlpNumSign (delta);
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
	    dbl_EGlpNumReallocArray (&(b->bounds), b->bound_cnt + nnew);
	    for (i = 0; i < nnew; i++) {
		b->bound_indx[b->bound_cnt + i] = new_indx[i];
		b->lu[b->bound_cnt + i] = new_lu[i];
		dbl_EGlpNumCopy (b->bounds[b->bound_cnt + i], new_bounds[i]);
	    }
	    b->bound_cnt += nnew;
	}
    }
    *hit = nnew;

CLEANUP:

    ILL_IFFREE (new_indx, int);
    ILL_IFFREE (new_lu, char);
    dbl_EGlpNumFreeArray (dj);
    dbl_EGlpNumFreeArray (new_bounds);
    dbl_EGlpNumClearVar (delta);
    dbl_EGlpNumClearVar (lpval);
    return rval;
}

static void dbl_best_bbnode (dbl_mipinfo * minf,
      dbl_bbnode ** best)
{
#if 0
    dbl_bbnode *b;
    double bestval = dbl_ILL_MAXDOUBLE;

    for (b = minf->head_bbnode.next; b; b = b->next) {
	if (b->bound < bestval) {
	    *best = b;
	    bestval = b->bound;
	}
    }
#endif

    double val;
    dbl_EGlpNumInitVar (val);
    dbl_ILLutil_priority_deletemin (minf->que, &val, (void **) best);
    dbl_EGlpNumClearVar (val);
}

static void dbl_put_bbnode (dbl_mipinfo * minf,
      dbl_bbnode * b)
{
    b->next = minf->head_bbnode.next;
    b->prev = &(minf->head_bbnode);
    if (b->next)
	b->next->prev = b;
    minf->head_bbnode.next = b;
}

static void dbl_remove_bbnode (dbl_bbnode * b)
{
    b->prev->next = b->next;
    if (b->next)
	b->next->prev = b->prev;
}

static int dbl_find_branch (dbl_mipinfo * minf,
      double *x,
      double *lpval,
      int *bvar)
{
    dbl_lpinfo *lp = minf->lp;
    int rval = 0;

    switch (minf->branching_rule) {
    case dbl_PENALTYBRANCH:
	rval = dbl_find_penalty_branch (lp, minf->pinf, x, minf->downpen,
	    minf->uppen, lpval, bvar);
	ILL_CLEANUP_IF (rval);
	break;
    case dbl_FIRSTBRANCH:
	dbl_find_first_branch (lp, x, bvar);
	break;
    case dbl_MIDDLEBRANCH:
	dbl_find_middle_branch (lp, x, bvar);
	break;
    case dbl_STRONGBRANCH:
	rval = dbl_find_strong_branch (lp, minf->pinf, x, bvar);
	ILL_CLEANUP_IF (rval);
	break;
    default:
	fprintf (stderr, "Unknown branching rule.\n");
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_find_branch");
}

static void dbl_find_first_branch (dbl_lpinfo * lp,
      double *x,
      int *bvar)
{
    int i, ibest = -1;
    dbl_ILLlpdata *qslp = lp->O;
    double t;
    dbl_EGlpNumInitVar (t);

    for (i = 0; i < qslp->nstruct; i++) {
	if (qslp->intmarker[i]) {
	    /* t = dbl_ILLutil_our_frac (x[i]); */
	    dbl_EGlpNumFloor (t, x[i]);
	    dbl_EGlpNumSubTo (t, x[i]);
	    dbl_EGlpNumSign (t);
	    if ((dbl_EGlpNumIsNeqZero (t, dbl_ILL_INTTOL)) &&
		(dbl_EGlpNumIsNeq (t, dbl_oneLpNum, dbl_ILL_INTTOL))) {
		ibest = i;
		break;
	    }
	}
    }
    *bvar = ibest;
    dbl_EGlpNumClearVar (t);
}

static void dbl_find_middle_branch (dbl_lpinfo * lp,
      double *x,
      int *bvar)
{
    int i, ibest = -1;
    double t, tbest;
    dbl_ILLlpdata *qlp = lp->O;
    dbl_EGlpNumInitVar (t);
    dbl_EGlpNumInitVar (tbest);
    dbl_EGlpNumSet (tbest, 0.5);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = dbl_ILLutil_our_frac (x[i]) - 0.5; if (t < 0.0) t = -t; */
	    dbl_EGlpNumFloor (t, x[i]);
	    dbl_EGlpNumMultUiTo (t, 2);
	    dbl_EGlpNumSubTo (t, dbl_oneLpNum);
	    dbl_EGlpNumDivUiTo (t, 2);
	    if (dbl_EGlpNumIsLess (t, dbl_zeroLpNum))
		dbl_EGlpNumSign (t);
	    /* if (t < tbest) */
	    if (dbl_EGlpNumIsLess (t, tbest)) {
		dbl_EGlpNumCopy (tbest, t);
		ibest = i;
	    }
	}
    }

    /* if (tbest < (0.5 - dbl_ILL_INTTOL)) */
    dbl_EGlpNumAddTo (tbest, dbl_ILL_INTTOL);
    if (dbl_EGlpNumIsLessDbl (tbest, 0.5)) {
	*bvar = ibest;
    } else {
	*bvar = -1;
    }
    dbl_EGlpNumClearVar (t);
    dbl_EGlpNumClearVar (tbest);
}

static int dbl_find_penalty_branch (dbl_lpinfo * lp,
      dbl_price_info * pinf,
      double *x,
      double *downpen,
      double *uppen,
      double *lpval,
      int *bvar)
{
    int rval = 0;
    int i, k, ibest = -1, ncand = 0, nneed = 0;
    dbl_ILLlpdata *qslp = lp->O;
    int *candidatelist = 0;
    int *needlist = 0;
    double *fval = 0;
    double *xlist = 0;
    double *newdown = 0;
    double *newup = 0;
    double a, t, tbest;
    dbl_EGlpNumInitVar (a);
    dbl_EGlpNumInitVar (t);
    dbl_EGlpNumInitVar (tbest);
    dbl_EGlpNumCopy (tbest, dbl_ILL_MINDOUBLE);

    ILL_SAFE_MALLOC (candidatelist, qslp->nstruct, int);
    ILL_SAFE_MALLOC (needlist, qslp->nstruct, int);
    fval = dbl_EGlpNumAllocArray (qslp->nstruct);
    xlist = dbl_EGlpNumAllocArray (qslp->nstruct);
    for (i = 0; i < qslp->nstruct; i++) {
	if (qslp->intmarker[i]) {
	    /* fval[i] = x[i] - floor(x[i]); */
	    dbl_EGlpNumFloor (fval[i], x[i]);
	    dbl_EGlpNumSubTo (fval[i], x[i]);
	    dbl_EGlpNumSign (fval[i]);
	    if ((dbl_EGlpNumIsNeqZero (fval[i], dbl_ILL_INTTOL)) &&
		(dbl_EGlpNumIsNeq (fval[i], dbl_oneLpNum, dbl_ILL_INTTOL))) {
		candidatelist[ncand++] = i;
		/* if (downpen[i] == -1.0) */
		dbl_EGlpNumSign (downpen[i]);
		if (dbl_EGlpNumIsEqqual (downpen[i], dbl_oneLpNum)) {
		    dbl_EGlpNumCopy (xlist[nneed], x[i]);
		    needlist[nneed++] = i;
		}
		dbl_EGlpNumSign (downpen[i]);
	    }
	}
    }

    if (nneed > 0) {
	newdown = dbl_EGlpNumAllocArray (nneed);
	newup = dbl_EGlpNumAllocArray (nneed);
	rval = dbl_ILLlib_strongbranch (lp, pinf, needlist, nneed,
	    0, newdown, newup,
	    5 * dbl_STRONG_PIVOTS, dbl_ILL_MAXDOUBLE);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < nneed; i++) {
	    k = needlist[i];
	    /* uppen[k] = (newup[i] - lpval) / (1.0 - fval[k]); */
	    dbl_EGlpNumCopyDiff (uppen[k], newup[i], *lpval);
	    dbl_EGlpNumCopyDiff (downpen[k], dbl_oneLpNum, fval[k]);
	    dbl_EGlpNumDivTo (uppen[k], downpen[k]);
	    /* downpen[k] = (newdown[i] - lpval) / fval[k]; */
	    dbl_EGlpNumCopyDiffRatio (downpen[k], newdown[i], *lpval, fval[k]);

	}
    }
    for (i = 0; i < ncand; i++) {
	k = candidatelist[i];
	/* t = dbl_ILL_BRANCH_PENALTY_VAL (downpen[k], uppen[k], fval[k]); */
	dbl_EGlpNumCopy (t, downpen[k]);
	dbl_EGlpNumMultTo (t, fval[k]);
	dbl_EGlpNumCopyDiff (a, dbl_oneLpNum, fval[k]);
	dbl_EGlpNumMultTo (a, uppen[k]);
	if (dbl_EGlpNumIsLess (t, a)) {
	    dbl_EGlpNumMultUiTo (t, dbl_ILL_BRANCH_PENALTY_WEIGHT);
	    dbl_EGlpNumAddTo (t, a);
	} else {
	    dbl_EGlpNumMultUiTo (a, dbl_ILL_BRANCH_PENALTY_WEIGHT);
	    dbl_EGlpNumAddTo (t, a);
	}
	dbl_EGlpNumDivUiTo (t, dbl_ILL_BRANCH_PENALTY_WEIGHT + 1);

	if (dbl_EGlpNumIsLess (tbest, t)) {
	    dbl_EGlpNumCopy (tbest, t);
	    ibest = k;
	}
    }

    *bvar = ibest;

CLEANUP:

    dbl_EGlpNumClearVar (a);
    dbl_EGlpNumClearVar (t);
    dbl_EGlpNumClearVar (tbest);
    dbl_EGlpNumFreeArray (newdown);
    dbl_EGlpNumFreeArray (newup);
    dbl_EGlpNumFreeArray (fval);
    dbl_EGlpNumFreeArray (xlist);
    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (needlist, int);
    return rval;
}

static int dbl_find_strong_branch (dbl_lpinfo * lp,
      dbl_price_info * pinf,
      double *x,
      int *bvar)
{
    int rval = 0;
    int i, ibest = -1, ncand = 0;
    int maxtrys = dbl_STRONG_CANDIDATES;
    double t, tbest;
    dbl_ILLlpdata *qlp = lp->O;
    int *candidatelist = 0;
    int *newlist = 0;
    int *perm = 0;
    double *tval = 0;
    double *xlist = 0;
    double *downpen = 0;
    double *uppen = 0;
    ILLrandstate rstate;
    dbl_EGlpNumInitVar (t);
    dbl_EGlpNumInitVar (tbest);
    dbl_EGlpNumCopy (tbest, dbl_ILL_MINDOUBLE);

    ILLutil_sprand (999, &rstate);
    ILL_SAFE_MALLOC (candidatelist, qlp->nstruct, int);
    tval = dbl_EGlpNumAllocArray (qlp->nstruct);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = dbl_ILLutil_our_frac (x[i]) - 0.5; if (t < 0.0) t = -t; */
	    dbl_EGlpNumFloor (t, x[i]);
	    dbl_EGlpNumSubTo (t, x[i]);
	    dbl_EGlpNumSign (t);
	    dbl_EGlpNumMultUiTo (t, 2);
	    dbl_EGlpNumSubTo (t, dbl_oneLpNum);
	    if (dbl_EGlpNumIsLess (t, dbl_zeroLpNum))
		dbl_EGlpNumSign (t);
	    /* if (t < (0.5 - dbl_ILL_INTTOL)) */
	    if (dbl_EGlpNumIsNeq (t, dbl_oneLpNum, dbl_ILL_INTTOL)) {
		candidatelist[ncand] = i;
		dbl_EGlpNumDivUiTo (t, 2);
		dbl_EGlpNumCopy (tval[ncand++], t);
	    }
	}
    }

    if (ncand > 0) {
	if (ncand > maxtrys) {
	    ILL_SAFE_MALLOC (perm, ncand, int);

	    for (i = 0; i < ncand; i++) {
		perm[i] = i;
	    }
	    dbl_ILLutil_EGlpNum_rselect (perm, 0, ncand - 1, maxtrys, tval, &rstate);

	    ILL_SAFE_MALLOC (newlist, maxtrys, int);

	    for (i = 0; i < maxtrys; i++) {
		newlist[i] = candidatelist[perm[i]];
	    }
	    ILL_IFFREE (candidatelist, int);
	    candidatelist = newlist;
	    newlist = 0;
	    ncand = maxtrys;
	}
	downpen = dbl_EGlpNumAllocArray (ncand);
	uppen = dbl_EGlpNumAllocArray (ncand);
	xlist = dbl_EGlpNumAllocArray (ncand);

	for (i = 0; i < ncand; i++) {
	    dbl_EGlpNumCopy (xlist[i], x[candidatelist[i]]);
	}

	rval = dbl_ILLlib_strongbranch (lp, pinf, candidatelist, ncand,
	    0, downpen, uppen, dbl_STRONG_PIVOTS,
	    dbl_ILL_MAXDOUBLE);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < ncand; i++) {
	    /* t = dbl_ILL_BRANCH_STRONG_VAL (downpen[i], uppen[i]); */
	    if (dbl_EGlpNumIsLess (downpen[i], uppen[i])) {
		dbl_EGlpNumCopy (t, downpen[i]);
		dbl_EGlpNumMultUiTo (t, dbl_ILL_BRANCH_STRONG_WEIGHT);
		dbl_EGlpNumAddTo (t, uppen[i]);
	    } else {
		dbl_EGlpNumCopy (t, uppen[i]);
		dbl_EGlpNumMultUiTo (t, dbl_ILL_BRANCH_STRONG_WEIGHT);
		dbl_EGlpNumAddTo (t, downpen[i]);
	    }
	    dbl_EGlpNumDivUiTo (t, dbl_ILL_BRANCH_STRONG_WEIGHT + 1);
	    if (dbl_EGlpNumIsLess (tbest, t)) {
		dbl_EGlpNumCopy (tbest, t);
		ibest = candidatelist[i];
	    }
	}
    }
    *bvar = ibest;


CLEANUP:

    dbl_EGlpNumClearVar (t);
    dbl_EGlpNumClearVar (tbest);
    dbl_EGlpNumFreeArray (tval);
    dbl_EGlpNumFreeArray (xlist);
    dbl_EGlpNumFreeArray (uppen);
    dbl_EGlpNumFreeArray (downpen);
    ILL_IFFREE (candidatelist, int);
    ILL_IFFREE (newlist, int);
    ILL_IFFREE (perm, int);

    ILL_RETURN (rval, "dbl_find_strong_branch");
}

static void dbl_check_integral (dbl_lpinfo * lp,
      double *x,
      int *yesno)
{
    int i;
    double t;
    dbl_ILLlpdata *qlp = lp->O;
    dbl_EGlpNumInitVar (t);

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    /* t = dbl_ILLutil_our_frac (x[i]); */
	    dbl_EGlpNumFloor (t, x[i]);
	    dbl_EGlpNumSubTo (t, x[i]);
	    dbl_EGlpNumSign (t);
	    /* if (t > dbl_ILL_INTTOL && t < 1.0 - dbl_ILL_INTTOL) */
	    if ((dbl_EGlpNumIsNeqZero (t, dbl_ILL_INTTOL)) &&
		(dbl_EGlpNumIsNeq (t, dbl_oneLpNum, dbl_ILL_INTTOL))) {
		*yesno = 0;
		dbl_EGlpNumClearVar (t);
		return;
	    }
	}
    }

    *yesno = 1;
    dbl_EGlpNumClearVar (t);
}

static int dbl_plunge (dbl_mipinfo * minf)
{
    int rval = 0;
    int i, status;
    dbl_lpinfo *lp = minf->lp;
    dbl_ILLlpdata *qlp = minf->lp->O;
    double *oldlower = 0;
    double *oldupper = 0;

    if (minf->watch) {
	printf ("Plunging ...\n");
	fflush (stdout);
    }
    oldlower = dbl_EGlpNumAllocArray (qlp->nstruct);
    oldupper = dbl_EGlpNumAllocArray (qlp->nstruct);

    for (i = 0; i < qlp->nstruct; i++) {
	dbl_EGlpNumCopy (oldlower[i], minf->lower[i]);
	dbl_EGlpNumCopy (oldupper[i], minf->upper[i]);
    }

    rval = dbl_plunge_work (minf, 0);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < qlp->nstruct; i++) {
	rval = dbl_ILLlib_chgbnd (lp, i, 'L', oldlower[i]);
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLlib_chgbnd (lp, i, 'U', oldupper[i]);
	ILL_CLEANUP_IF (rval);
	dbl_EGlpNumCopy (minf->lower[i], oldlower[i]);
	dbl_EGlpNumCopy (minf->upper[i], oldupper[i]);
    }

    rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    dbl_EGlpNumFreeArray (oldlower);
    dbl_EGlpNumFreeArray (oldupper);

    ILL_RETURN (rval, "dbl_plunge");
}

static int dbl_plunge_work (dbl_mipinfo * minf,
      int depth)
{
    int rval = 0;
    int bvar, status, count;
    double lpval, val0, val1, int_tol;
    dbl_lpinfo *lp = minf->lp;
    dbl_EGlpNumInitVar (lpval);
    dbl_EGlpNumInitVar (val0);
    dbl_EGlpNumInitVar (val1);
    dbl_EGlpNumInitVar (int_tol);
    dbl_EGlpNumSet (int_tol, 0.001);

    rval = dbl_ILLlib_get_x (lp, 0, minf->x);
    ILL_CLEANUP_IF (rval);

    rval = dbl_round_variables (minf, &count, &int_tol /* 0.001 */ );
    ILL_CLEANUP_IF (rval);
    if (count) {
	rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	if (status != QS_LP_OPTIMAL) {
	    goto CLEANUP;
	}
	rval = dbl_ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);
    }
    dbl_find_middle_branch (lp, minf->x, &bvar);
    if (bvar == -1) {
	rval = dbl_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);

	if (dbl_EGlpNumIsLess (lpval, minf->value)) {
	    printf ("Plunge Integral Solution: %.6f (Depth: %d)\n",
		dbl_EGlpNumToLf (lpval), depth);
	    fflush (stdout);

	    dbl_EGlpNumCopy (minf->value, lpval);
	    dbl_EGlpNumCopyDiff (minf->objectivebound, lpval, dbl_ILL_INTTOL);
	    dbl_copy_x (lp->O->nstruct, minf->x, minf->bestx);
	}
	goto CLEANUP;
    }
    dbl_EGlpNumOne (minf->lower[bvar]);
    rval = dbl_ILLlib_chgbnd (lp, bvar, 'L', dbl_oneLpNum);
    ILL_CLEANUP_IF (rval);
    rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the dbl_plunge LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
	dbl_EGlpNumCopy (val1, dbl_ILL_MAXDOUBLE);
    } else if (status == QS_LP_OPTIMAL) {
	rval = dbl_ILLlib_objval (lp, 0, &val1);
	ILL_CLEANUP_IF (rval);
    } else {
	ILL_CLEANUP;
    }

    rval = dbl_ILLlib_chgbnd (lp, bvar, 'L', dbl_zeroLpNum);
    ILL_CLEANUP_IF (rval);
    dbl_EGlpNumZero (minf->lower[bvar]);

    dbl_EGlpNumZero (minf->upper[bvar]);
    rval = dbl_ILLlib_chgbnd (lp, bvar, 'U', dbl_zeroLpNum);
    ILL_CLEANUP_IF (rval);
    rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

    if (status == QS_LP_UNSOLVED) {
	printf ("Simplex did not solve the dbl_plunge LP\n");
	fflush (stdout);
	rval = 1;
	ILL_CLEANUP;
    } else if (status == QS_LP_INFEASIBLE) {
	dbl_EGlpNumCopy (val0, dbl_ILL_MAXDOUBLE);
    } else if (status == QS_LP_OPTIMAL) {
	rval = dbl_ILLlib_objval (lp, 0, &val0);
	ILL_CLEANUP_IF (rval);
    } else {
	ILL_CLEANUP;
    }

    rval = dbl_ILLlib_chgbnd (lp, bvar, 'U', dbl_oneLpNum);
    ILL_CLEANUP_IF (rval);
    dbl_EGlpNumCopy (minf->upper[bvar], dbl_oneLpNum);

    if (dbl_EGlpNumIsEqqual (val0, dbl_ILL_MAXDOUBLE) &&
	dbl_EGlpNumIsEqqual (val1, dbl_ILL_MAXDOUBLE)) {
	ILL_CLEANUP;
    }
    if (dbl_EGlpNumIsLess (val0, val1)) {
	dbl_EGlpNumZero (minf->upper[bvar]);
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'U', dbl_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	rval = dbl_plunge_work (minf, depth + 1);
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'U', dbl_oneLpNum);
	ILL_CLEANUP_IF (rval);
	dbl_EGlpNumOne (minf->upper[bvar]);
    } else {
	dbl_EGlpNumOne (minf->lower[bvar]);
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'L', dbl_oneLpNum);
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
	rval = dbl_plunge_work (minf, depth + 1);
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLlib_chgbnd (lp, bvar, 'L', dbl_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	dbl_EGlpNumZero (minf->lower[bvar]);
    }

CLEANUP:

    dbl_EGlpNumClearVar (lpval);
    dbl_EGlpNumClearVar (val0);
    dbl_EGlpNumClearVar (val1);
    dbl_EGlpNumClearVar (int_tol);
    ILL_RETURN (rval, "dbl_plunge_work");
}

static int dbl_round_variables (dbl_mipinfo * minf,
      int *count,
      double *tol)
{
    int rval = 0;
    int i, hit = 0;
    dbl_lpinfo *lp = minf->lp;
    dbl_ILLlpdata *qlp = lp->O;

    *count = 0;

    for (i = 0; i < qlp->nstruct; i++) {
	if (qlp->intmarker[i]) {
	    if (dbl_EGlpNumIsNeqq (minf->lower[i], minf->upper[i])) {
		if (dbl_EGlpNumIsLess (minf->x[i], *tol)) {
		    dbl_EGlpNumZero (minf->upper[i]);
		    rval = dbl_ILLlib_chgbnd (lp, i, 'U', dbl_zeroLpNum);
		    ILL_CLEANUP_IF (rval);
		    hit++;
		} else if (dbl_EGlpNumIsEqual (minf->x[i], dbl_oneLpNum, *tol)) {
		    dbl_EGlpNumOne (minf->lower[i]);
		    rval = dbl_ILLlib_chgbnd (lp, i, 'L', dbl_oneLpNum);
		    ILL_CLEANUP_IF (rval);
		    hit++;
		}
	    }
	}
    }
    *count = hit;

CLEANUP:

    ILL_RETURN (rval, "dbl_round_variables");
}

static void dbl_copy_x (int nstruct,
      double *from_x,
      double *to_x)
{
    int j;

    for (j = 0; j < nstruct; j++) {
	dbl_EGlpNumCopy (to_x[j], from_x[j]);
    }
}

static void dbl_init_mipinfo (dbl_mipinfo * minf)
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
	minf->branching_rule = /* dbl_MIDDLEBRANCH */ dbl_STRONGBRANCH;
	minf->watch = 1;
	dbl_EGlpNumInitVar (minf->objectivebound);
	dbl_EGlpNumInitVar (minf->value);
	dbl_EGlpNumCopy (minf->objectivebound, dbl_ILL_MAXDOUBLE);
	dbl_EGlpNumCopy (minf->value, dbl_ILL_MAXDOUBLE);
	ILLptrworld_init (&minf->ptrworld);
    }
}

static void dbl_free_mipinfo (dbl_mipinfo * minf)
{
    int total, onlist;

    if (minf) {
	dbl_EGlpNumFreeArray (minf->downpen);
	dbl_EGlpNumFreeArray (minf->uppen);
	dbl_EGlpNumFreeArray (minf->x);
	dbl_EGlpNumFreeArray (minf->bestx);
	dbl_EGlpNumFreeArray (minf->lower);
	dbl_EGlpNumFreeArray (minf->upper);
	bbnode_listfree (&minf->ptrworld, minf->head_bbnode.next);
	if (bbnode_check_leaks (&minf->ptrworld, &total, &onlist)) {
	    fprintf (stderr, "WARNING: %d outstanding bbnodes\n", total - onlist);
	}
	ILLptrworld_delete (&minf->ptrworld);
	dbl_EGlpNumClearVar ((minf->objectivebound));
	dbl_EGlpNumClearVar ((minf->value));
	memset (minf, 0, sizeof (dbl_mipinfo));
	/* dbl_init_mipinfo (minf); */
    }
}

static void dbl_init_bbnode (dbl_bbnode * b)
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
	dbl_EGlpNumInitVar ((b->bound));
	dbl_EGlpNumCopy (b->bound, dbl_ILL_MINDOUBLE);
    }
}

static void dbl_free_bbnode (dbl_bbnode * b)
{
    if (b) {
	dbl_EGlpNumFreeArray (b->rownorms);
	dbl_EGlpNumFreeArray (b->bounds);
	ILL_IFFREE (b->cstat, char);
	ILL_IFFREE (b->rstat, char);
	ILL_IFFREE (b->bound_indx, int);
	ILL_IFFREE (b->lu, char);
	dbl_EGlpNumClearVar ((b->bound));
	memset (b, 0, sizeof (dbl_bbnode));
    }
}
