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

/* RCS_INFO = "$RCSfile: dbl_presolve.c,v $ $Revision: 1.2 $ $Date:
   2003/11/05 16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* Presolve Routine for Simplex Method                      */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int dbl_ILLlp_add_logicals (ILLlpata *lp)                                 */
/* int dbl_ILLlp_presolve (dbl_ILLlpdata *lp)                                    */
/* int dbl_ILLlp_scale (dbl_ILLlpdata *lp)                                       */
/* void dbl_ILLlp_sinfo_init (dbl_ILLlp_sinfo *sinfo)                            */
/* void dbl_ILLlp_sinfo_free (dbl_ILLlp_sinfo *sinfo)                            */
/* void dbl_ILLlp_predata_init (dbl_ILLlp_predata *pre)                          */
/* void dbl_ILLlp_predata_free (dbl_ILLlp_predata *pre)                          */
/* */
/* NOTES                                                                   */
/* */
/* presolve will assume that logicals have been added.                   */
/* */
/****************************************************************************/

#include "econfig.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdata.h"
#include "dbl_lpdefs.h"
/* extern double dbl_SZERO_TOLER; */

#define dbl_ILL_LP_STATUS_OK (0)
#define dbl_ILL_PRE_FEAS_TOL dbl_PFEAS_TOLER	/* (1e-6) */
#define dbl_ILL_PRE_ZERO_TOL dbl_PIVOT_TOLER	/* (1e-10) */

#define dbl_ILL_PRE_DELETE_EMPTY_ROW               (1)
#define dbl_ILL_PRE_DELETE_SINGLETON_ROW           (2)
#define dbl_ILL_PRE_DELETE_FIXED_VARIABLE          (3)
#define dbl_ILL_PRE_DELETE_FORCED_VARIABLE         (4)
#define dbl_ILL_PRE_DELETE_SINGLETON_VARIABLE      (5)
#define dbl_ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE (6)
#define dbl_ILL_PRE_DELETE_EMPTY_COLUMN            (7)

#define dbl_ILL_PRE_COL_STRUC                      (0)
#define dbl_ILL_PRE_COL_LOGICAL                    (1)

static int dbl_debug = 0;

typedef struct {
    int row;
    int col;
    char coltype;
    char mark;
    char del;
    double coef;
}
  dbl_edge;

typedef struct dbl_node {
    dbl_edge **adj;
    double obj;
    double lower;
    double upper;
    double rhs;
    int deg;
    char mark;
    char del;
    char coltype;
    char rowsense;
}
  dbl_node;

typedef struct dbl_intptr {
    int this;
    struct dbl_intptr *next;
}
  dbl_intptr;

typedef struct dbl_graph {
    dbl_edge *edgelist;
    struct dbl_node *rows;
    struct dbl_node *cols;
    int ecount;
    int nrows;
    int ncols;
    int nzcount;
    dbl_edge **adjspace;
    ILLptrworld intptrworld;
    int objsense;
}
  dbl_graph;

void dbl_ILLlp_sinfo_init (dbl_ILLlp_sinfo * sinfo),
  dbl_ILLlp_sinfo_free (dbl_ILLlp_sinfo * sinfo),
  dbl_ILLlp_predata_init (dbl_ILLlp_predata * pre),
  dbl_ILLlp_predata_free (dbl_ILLlp_predata * pre),
  dbl_ILLlp_preop_init (dbl_ILLlp_preop * op),
  dbl_ILLlp_preop_free (dbl_ILLlp_preop * op),
  dbl_ILLlp_preline_init (dbl_ILLlp_preline * line),
  dbl_ILLlp_preline_free (dbl_ILLlp_preline * line);

int dbl_ILLlp_sinfo_print (dbl_ILLlp_sinfo * s);

static void dbl_set_fixed_variable (dbl_graph * G,
      int j,
      double val),
  dbl_get_implied_rhs_bounds (dbl_graph * G,
      int i,
      double *lb,
      double *ub),
  dbl_get_implied_variable_bounds (dbl_graph * G,
      int j,
      dbl_edge * a_ij,
      double *lb,
      double *ub),
  dbl_dump_line (dbl_ILLlp_preline * line),
  dbl_init_graph (dbl_graph * G),
  dbl_free_graph (dbl_graph * G),
  dbl_dump_graph (dbl_graph * G);

static int dbl_simple_presolve (dbl_ILLlpdata * lp,
      dbl_ILLlp_predata * pre,
      dbl_ILLlp_sinfo * info,
      int pre_types,
      int *status),
  dbl_grab_lp_line (dbl_graph * G,
      int indx,
      dbl_ILLlp_preline * line,
      int row_or_col),
  dbl_grab_lp_info (dbl_graph * G,
      char **colnames,
      dbl_ILLlp_sinfo * info),
  dbl_fixed_variables (dbl_graph * G,
      dbl_ILLlp_predata * pre),
  dbl_empty_columns (dbl_graph * G,
      dbl_ILLlp_predata * pre),
  dbl_singleton_rows (dbl_graph * G,
      dbl_ILLlp_predata * pre,
      int *hit),
  dbl_forcing_constraints (dbl_graph * G,
      dbl_ILLlp_predata * pre,
      int *hit),
  dbl_singleton_columns (dbl_graph * G,
      dbl_ILLlp_predata * pre,
      int *hit),
  dbl_duplicate_rows (dbl_graph * G,
      int *hit),
  dbl_duplicate_cols (dbl_graph * G,
      int *hit),
  dbl_gather_dup_lists (int *s,
      int count,
      int *duptotal,
      int **dupcnt,
      int **dupind),
  dbl_get_next_preop (dbl_ILLlp_predata * pre,
      dbl_ILLlp_preop ** op),
  dbl_add_to_list (ILLptrworld * world,
      dbl_intptr ** list,
      int i),
  dbl_build_graph (dbl_ILLlpdata * lp,
      dbl_graph * G);


ILL_PTRWORLD_ROUTINES (dbl_intptr, intptralloc, intptr_bulkalloc, intptrfree)
ILL_PTRWORLD_LISTFREE_ROUTINE (dbl_intptr, intptr_listfree, intptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (dbl_intptr, intptr_check_leaks, this, int)
int dbl_ILLlp_add_logicals (dbl_ILLlpdata * lp)
{
    int rval = 0;
    int ncols, nrows, nzcount, i, aindex;
    char *sense;
    dbl_ILLmatrix *A;

    if (!lp) {
	fprintf (stderr, "dbl_ILLlp_add_logicals called with a NULL pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    printf ("dbl_ILLlp_add_logicals ...\n");
    fflush (stdout);

    A = &lp->A;
    sense = lp->sense;
    ncols = lp->ncols;
    nrows = lp->nrows;
    nzcount = lp->nzcount;

    if (nrows == 0)
	goto CLEANUP;
    dbl_EGlpNumReallocArray (&(lp->obj), lp->colsize + nrows);
    dbl_EGlpNumReallocArray (&(lp->upper), lp->colsize + nrows);
    dbl_EGlpNumReallocArray (&(lp->lower), lp->colsize + nrows);
    lp->colnames = EGrealloc (lp->colnames, sizeof (char *) * (lp->colsize + nrows));
    /*
	    rval = ILLutil_reallocrus_count ((void **) &(lp->colnames), lp->colsize + nrows, sizeof (char *));
	    ILL_CLEANUP_IF (rval);
    */
    memset (lp->colnames + ncols, 0, sizeof (char *) * nrows);

    ILL_SAFE_MALLOC (lp->rowmap, lp->rowsize, int);


    A->matcnt = EGrealloc (A->matcnt, sizeof (int) * (A->matcolsize + nrows));
    /*
	    rval = ILLutil_reallocrus_count ((void **) &(A->matcnt), A->matcolsize + nrows, sizeof (int));
	    ILL_CLEANUP_IF (rval);
    */

    A->matbeg = EGrealloc (A->matbeg, sizeof (int) * (A->matcolsize + nrows));
    /*
	    rval = ILLutil_reallocrus_count ((void **) &(A->matbeg), A->matcolsize + nrows, sizeof (int));
	    ILL_CLEANUP_IF (rval);
    */

    A->matind = EGrealloc (A->matind, sizeof (int) * (A->matsize + nrows));
    /*
	    rval = ILLutil_reallocrus_count ((void **) &(A->matind), A->matsize + nrows, sizeof (int));
	    ILL_CLEANUP_IF (rval);
    */
    dbl_EGlpNumReallocArray (&(A->matval), A->matsize + nrows);

    for (i = 0; i < nrows; i++) {
	A->matind[A->matsize + i] = -1;
    }

    aindex = A->matsize - A->matfree;

    for (i = 0; i < nrows; i++) {
	lp->rowmap[i] = ncols;
	dbl_EGlpNumZero (lp->obj[ncols]);
	A->matcnt[ncols] = 1;
	A->matbeg[ncols] = aindex;
	A->matind[aindex] = i;
	switch (sense[i]) {
	case 'E':		/* Arificial */
	    dbl_EGlpNumZero (lp->lower[ncols]);
	    dbl_EGlpNumZero (lp->upper[ncols]);
	    dbl_EGlpNumOne (A->matval[aindex]);
	    break;
	case 'G':		/* Surplus   */
	    dbl_EGlpNumZero (lp->lower[ncols]);
	    dbl_EGlpNumCopy (lp->upper[ncols], dbl_ILL_MAXDOUBLE);
	    dbl_EGlpNumOne (A->matval[aindex]);
	    dbl_EGlpNumSign (A->matval[aindex]);
	    break;
	case 'L':		/* Slack     */
	    dbl_EGlpNumZero (lp->lower[ncols]);
	    dbl_EGlpNumCopy (lp->upper[ncols], dbl_ILL_MAXDOUBLE);
	    dbl_EGlpNumOne (A->matval[aindex]);
	    break;
	case 'R':		/* Range     */
	    dbl_EGlpNumZero (lp->lower[ncols]);
	    dbl_EGlpNumCopy (lp->upper[ncols], lp->rangeval[i]);
	    dbl_EGlpNumOne (A->matval[aindex]);
	    dbl_EGlpNumSign (A->matval[aindex]);
	    break;
	default:
	    fprintf (stderr, "unknown sense %c in dbl_ILLlp_add_logicals\n", sense[i]);
	    rval = 1;
	    goto CLEANUP;
	}
	ncols++;
	nzcount++;
	aindex++;
    }

    lp->ncols = ncols;
    lp->nzcount = nzcount;
    A->matcols = ncols;

    lp->colsize += nrows;
    A->matsize += nrows;
    A->matcolsize += nrows;

    if (lp->rA) {
	dbl_ILLlp_rows_clear (lp->rA);
    } else {
	ILL_SAFE_MALLOC (lp->rA, 1, dbl_ILLlp_rows);
    }

    rval = dbl_ILLlp_rows_init (lp->rA, lp, 1);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_ILLlp_add_logicals");
}

int dbl_ILLlp_scale (dbl_ILLlpdata * lp)
{
    int rval = 0;
    int i, j, k, col, row, nstruct, start, stop;
    dbl_ILLmatrix *A;
    double rho;
    double *gama = 0;
    dbl_EGlpNumInitVar (rho);

    /* Columns - divide by largest absolute value */

    if (!lp) {
	ILL_ERROR (rval, "dbl_ILLlp_scale called with a NULL pointer");
    }
    if (lp->nrows == 0 || lp->ncols == 0)
	goto CLEANUP;

    A = &lp->A;
    nstruct = lp->nstruct;

    for (j = 0; j < nstruct; j++) {
	col = lp->structmap[j];
	dbl_EGlpNumZero (rho);

	start = A->matbeg[col];
	stop = start + A->matcnt[col];

	for (k = start; k < stop; k++) {
	    dbl_EGlpNumSetToMaxAbs (rho, A->matval[k]);
	}

	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, rho)) {
	    for (k = start; k < stop; k++) {
		dbl_EGlpNumDivTo (A->matval[k], rho);
	    }
	    dbl_EGlpNumDivTo (lp->obj[col], rho);
	    if (dbl_EGlpNumIsNeqq (lp->lower[col], dbl_ILL_MINDOUBLE))
		dbl_EGlpNumMultTo (lp->lower[col], rho);
	    if (dbl_EGlpNumIsNeqq (lp->upper[col], dbl_ILL_MAXDOUBLE))
		dbl_EGlpNumMultTo (lp->upper[col], rho);
	}
    }

    gama = dbl_EGlpNumAllocArray (lp->nrows);
    for (i = 0; i < lp->nrows; i++) {
	dbl_EGlpNumZero (gama[i]);
    }

    for (j = 0; j < nstruct; j++) {
	col = lp->structmap[j];
	start = A->matbeg[col];
	stop = start + A->matcnt[col];

	for (k = start; k < stop; k++) {
	    row = A->matind[k];
	    dbl_EGlpNumSetToMaxAbs (gama[row], A->matval[k]);
	}
    }

    for (j = 0; j < nstruct; j++) {
	col = lp->structmap[j];
	start = A->matbeg[col];
	stop = start + A->matcnt[col];

	for (k = start; k < stop; k++) {
	    row = A->matind[k];
	    if (dbl_EGlpNumIsLess (dbl_zeroLpNum, gama[row])) {
		dbl_EGlpNumDivTo (A->matval[k], gama[row]);
	    }
	}
    }

    for (i = 0; i < lp->nrows; i++) {
	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, gama[i])) {
	    dbl_EGlpNumDivTo (lp->rhs[i], gama[i]);
	    col = lp->rowmap[i];
	    if (dbl_EGlpNumIsNeqq (lp->upper[col], dbl_ILL_MAXDOUBLE)) {
		dbl_EGlpNumDivTo (lp->upper[col], gama[i]);	/* Ranged row */
	    }
	}
    }

    if (lp->rA) {		/* Need to clear the row version of data */
	dbl_ILLlp_rows_clear (lp->rA);
	ILL_IFFREE (lp->rA, dbl_ILLlp_rows);
    }
CLEANUP:

    dbl_EGlpNumClearVar (rho);
    dbl_EGlpNumFreeArray (gama);
    ILL_RETURN (rval, "dbl_ILLlp_scale");
}

int dbl_ILLlp_presolve (dbl_ILLlpdata * lp,
      int pre_types)
{
    int rval = 0;
    int status = dbl_ILL_LP_STATUS_OK;
    dbl_ILLlp_predata *pre = 0;
    dbl_ILLlp_sinfo *info = 0;

    if (!lp) {
	fprintf (stderr, "dbl_ILLlp_presolve called with a NULL pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    /*
        ILLlpdata_writelp (lp, 0);
        printf ("\n"); fflush (stdout);
    */

    ILL_SAFE_MALLOC (pre, 1, dbl_ILLlp_predata);
    dbl_ILLlp_predata_init (pre);

    ILL_SAFE_MALLOC (info, 1, dbl_ILLlp_sinfo);
    dbl_ILLlp_sinfo_init (info);

    rval = dbl_simple_presolve (lp, pre, info, pre_types, &status);
    ILL_CLEANUP_IF (rval);
    if (status != dbl_ILL_LP_STATUS_OK) {
	printf ("dbl_simple_presolve returned with bad status\n");
	rval = 1;
	goto CLEANUP;
    }
    /*
        rval = dbl_ILLlp_sinfo_print (info);
        ILL_CLEANUP_IF (rval);
    */

CLEANUP:

    if (rval) {
	if (pre) {
	    dbl_ILLlp_predata_free (pre);
	    ILL_IFFREE (pre, dbl_ILLlp_predata);
	}
	if (info) {
	    dbl_ILLlp_sinfo_free (info);
	    ILL_IFFREE (info, dbl_ILLlp_sinfo);
	}
    } else {
	lp->presolve = pre;
	lp->sinfo = info;
    }

    ILL_RETURN (rval, "dbl_ILLlp_presolve");
}


#if 0
int ILLlp_presolve_addrow (dbl_lpinfo * lp,
      int cnt,
      int *ind,
      double *val,
      double rhs)
{
    int rval = 0;
    dbl_ILLlpdata *qslp;
    dbl_ILLlp_sinfo *S;
    dbl_ILLmatrix *A;

    /* This will need to evolve into a function that handles the task */
    /* of working through the presolve data to determine the new LP   */
    /* created when a row is added to the original LP.                */

    /* The copies of the obj and bound used in the simplex code are   */
    /* also updated in this function.                                 */

    if (!lp) {
	fprintf (stderr, "ILLlp_presolve_addrow is called without an LP\n");
	rval = 1;
	goto CLEANUP;
    }
    if (lp->presolve != 0) {
	fprintf (stderr, "Not yet set up to handle addrows after presolve\n");
	rval = 1;
	goto CLEANUP;
    }
    qslp = lp->O;
    S = qslp->sinfo;
    A = S->A;


    rval = ILLlib_matrix_addrow (A, cnt, ind, val, rhs);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    ILL_RETURN (rval, "ILLlp_presolve_addrow");
}
#endif


static int dbl_simple_presolve (dbl_ILLlpdata * lp,
      dbl_ILLlp_predata * pre,
      dbl_ILLlp_sinfo * info,
      int pre_types,
      int *status)
{
    int rval = 0;
    int i, hit, newhit;
    dbl_graph G;

    if (status)
	*status = dbl_ILL_LP_STATUS_OK;
    dbl_init_graph (&G);

    if (!lp) {
	fprintf (stderr, "dbl_simple_presolve called with a NULL pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    printf ("Initial Rows = %d, Cols = %d, Nzcount = %d\n",
	lp->nrows, lp->ncols, lp->nzcount);
    fflush (stdout);

    rval = dbl_build_graph (lp, &G);
    ILL_CLEANUP_IF (rval);
    if (dbl_debug)
	dbl_dump_graph (&G);

    if (pre_types & dbl_ILL_PRE_FIXED) {
	rval = dbl_fixed_variables (&G, pre);
	ILL_CLEANUP_IF (rval);
    }
    do {
	hit = 0;
	if (pre_types & dbl_ILL_PRE_SINGLE_ROW) {
	    rval = dbl_singleton_rows (&G, pre, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & dbl_ILL_PRE_FORCING) {
	    rval = dbl_forcing_constraints (&G, pre, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & dbl_ILL_PRE_SINGLE_COL) {
	    rval = dbl_singleton_columns (&G, pre, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & dbl_ILL_PRE_DUPLICATE_ROW) {
	    rval = dbl_duplicate_rows (&G, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & dbl_ILL_PRE_DUPLICATE_COL) {
	    rval = dbl_duplicate_cols (&G, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	/*
	        {
	            int k, cnt = 0;
	            for (i = 0; i < G.ncols; i++) {
	                if (G.cols[i].del == 0) {
	                    for (k = 0; k < G.cols[i].deg; k++)  {
	                        if (G.cols[i].adj[k]->del == 0) {
	                            cnt++;
	                        }
	                    }
	                }
	            }
	            printf ("Current NZCOUNT = %d\n", cnt); fflush (stdout);
	        }
	*/
    } while (hit);

    if (dbl_ILL_PRE_EMPTY_COL) {
	rval = dbl_empty_columns (&G, pre);
	ILL_CLEANUP_IF (rval);
    }
    if (dbl_debug) {
	printf ("Operations\n");
	for (i = 0; i < pre->opcount; i++) {
	    switch (pre->oplist[i].ptype) {
	    case dbl_ILL_PRE_DELETE_EMPTY_ROW:
		printf ("Delete Empty Row: %d\n", pre->oplist[i].rowindex);
		fflush (stdout);
		break;
	    case dbl_ILL_PRE_DELETE_SINGLETON_ROW:
		printf ("Delete Singleton Row: %d (col %d)\n",
		    pre->oplist[i].rowindex, pre->oplist[i].colindex);
		fflush (stdout);
		dbl_dump_line (&pre->oplist[i].line);
		break;
	    case dbl_ILL_PRE_DELETE_FIXED_VARIABLE:
		printf ("Delete Fixed Variable: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		dbl_dump_line (&pre->oplist[i].line);
		break;
	    case dbl_ILL_PRE_DELETE_FORCED_VARIABLE:
		printf ("Delete Forced Variable: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		dbl_dump_line (&pre->oplist[i].line);
		break;
	    case dbl_ILL_PRE_DELETE_SINGLETON_VARIABLE:
		printf ("Delete Singleton Variable: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		dbl_dump_line (&pre->oplist[i].line);
		break;
	    case dbl_ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE:
		printf ("Delete Free Singleton Variable: %d\n",
		    pre->oplist[i].colindex);
		fflush (stdout);
		dbl_dump_line (&pre->oplist[i].line);
		break;
	    case dbl_ILL_PRE_DELETE_EMPTY_COLUMN:
		printf ("Delete Empty Column: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		dbl_dump_line (&pre->oplist[i].line);
		break;
	    default:
		fprintf (stderr, "unknon presolve operation\n");
		rval = 1;
		goto CLEANUP;
	    }
	}
	printf ("\n");
    }
    rval = dbl_grab_lp_info (&G, lp->colnames, info);
    ILL_CLEANUP_IF (rval);

    /*
        printf ("Final Rows = %d, Cols = %d, Nzcount = %d\n",
                   info->nrows, info->ncols, info->nzcount);
        fflush (stdout);
    */


CLEANUP:

    dbl_free_graph (&G);
    ILL_RETURN (rval, "dbl_simple_presolve");
}

static int dbl_grab_lp_line (dbl_graph * G,
      int indx,
      dbl_ILLlp_preline * line,
      int row_or_col)
{
    int rval = 0;
    int k, cnt;
    dbl_node *n;

    if (row_or_col == 0)
	n = &G->rows[indx];
    else
	n = &G->cols[indx];

    line->count = 0;

    for (k = 0; k < n->deg; k++) {
	if (n->adj[k]->del == 0) {
	    line->count++;
	}
    }

    if (line->count) {
	ILL_SAFE_MALLOC (line->ind, line->count, int);
	line->val = dbl_EGlpNumAllocArray (line->count);
	if (!line->ind || !line->val) {
	    fprintf (stderr, "out of memory in dbl_grab_lp_line\n");
	    rval = 1;
	    goto CLEANUP;
	}
	for (k = 0, cnt = 0; k < n->deg; k++) {
	    if (n->adj[k]->del == 0) {
		line->ind[cnt] = n->adj[k]->row;
		dbl_EGlpNumCopy (line->val[cnt], n->adj[k]->coef);
		cnt++;
	    }
	}
    }
    if (row_or_col == 0) {
	dbl_EGlpNumCopy (line->rhs, n->rhs);
    } else {
	dbl_EGlpNumCopy (line->obj, n->obj);
	dbl_EGlpNumCopy (line->lower, n->lower);
	dbl_EGlpNumCopy (line->upper, n->upper);
    }

    line->row_or_col = row_or_col;

CLEANUP:

    ILL_RETURN (rval, "dbl_grab_lp_line");
}

static void dbl_dump_line (dbl_ILLlp_preline * line)
{
    int k;

    printf (" ");
    if (line->row_or_col == 0) {
	for (k = 0; k < line->count; k++) {
	    printf (" C%d->%g", line->ind[k], dbl_EGlpNumToLf (line->val[k]));
	}
	printf (" RHS->%g\n", dbl_EGlpNumToLf (line->rhs));
    } else {
	for (k = 0; k < line->count; k++) {
	    printf (" R%d->%g", line->ind[k], dbl_EGlpNumToLf (line->val[k]));
	}
	printf (" Obj->%g  LB->%g  UB->%g\n", dbl_EGlpNumToLf (line->obj),
	    dbl_EGlpNumToLf (line->lower), dbl_EGlpNumToLf (line->upper));
    }
    fflush (stdout);
}

static int dbl_grab_lp_info (dbl_graph * G,
      char **colnames,
      dbl_ILLlp_sinfo * info)
{
    int rval = 0;
    int ncols = 0, nrows = 0, nzcount = 0;
    int i, j, k, cnt, len;
    dbl_node *grows = G->rows;
    dbl_node *gcols = G->cols;
    int *tdeg = 0;
    int *map = 0;
    char *buf = 0;
    dbl_ILLmatrix *A = &info->A;

    ILL_SAFE_MALLOC (tdeg, G->ncols, int);
    ILL_SAFE_MALLOC (map, G->nrows, int);
    if (!tdeg || !map) {
	fprintf (stderr, "out of memory in dbl_grab_lp_info\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < G->nrows; i++) {
	if (grows[i].del == 0) {
	    map[i] = nrows;
	    nrows++;
	}
    }

    for (j = 0; j < G->ncols; j++) {
	if (gcols[j].del == 0) {
	    tdeg[ncols] = 0;
	    for (k = 0; k < gcols[j].deg; k++) {
		if (gcols[j].adj[k]->del == 0) {
		    tdeg[ncols]++;
		    nzcount++;
		}
	    }
	    ncols++;
	}
    }

    info->ncols = ncols;
    info->nrows = nrows;
    info->nzcount = nzcount;

    info->rowsize = nrows;
    info->colsize = ncols;

    info->rhs = dbl_EGlpNumAllocArray (nrows);
    info->obj = dbl_EGlpNumAllocArray (ncols);
    info->upper = dbl_EGlpNumAllocArray (ncols);
    info->lower = dbl_EGlpNumAllocArray (ncols);
    A->matval = dbl_EGlpNumAllocArray (info->nzcount + 1);
    ILL_SAFE_MALLOC (A->matind, info->nzcount + 1, int);
    ILL_SAFE_MALLOC (A->matcnt, info->colsize, int);
    ILL_SAFE_MALLOC (A->matbeg, info->colsize, int);

    if (!info->rhs || !info->obj || !info->lower || !info->upper ||
	!A->matval || !A->matind || !A->matcnt || !A->matbeg) {
	fprintf (stderr, "out of memory in grab_lp\n");
	rval = 1;
	goto CLEANUP;
    }
    A->matind[info->nzcount] = -1;
    A->matsize = info->nzcount + 1;
    A->matcolsize = info->colsize;
    A->matfree = 1;
    A->matcols = ncols;
    A->matrows = nrows;


    nrows = 0;
    for (i = 0; i < G->nrows; i++) {
	if (grows[i].del == 0) {
	    dbl_EGlpNumCopy (info->rhs[nrows], grows[i].rhs);
	    nrows++;
	}
    }

    ncols = 0;
    cnt = 0;
    for (j = 0; j < G->ncols; j++) {
	if (gcols[j].del == 0) {
	    dbl_EGlpNumCopy (info->obj[ncols], gcols[j].obj);
	    dbl_EGlpNumCopy (info->lower[ncols], gcols[j].lower);
	    dbl_EGlpNumCopy (info->upper[ncols], gcols[j].upper);
	    A->matcnt[ncols] = tdeg[ncols];
	    A->matbeg[ncols] = cnt;
	    for (k = 0; k < gcols[j].deg; k++) {
		if (gcols[j].adj[k]->del == 0) {
		    dbl_EGlpNumCopy (A->matval[cnt], gcols[j].adj[k]->coef);
		    A->matind[cnt] = map[gcols[j].adj[k]->row];
		    cnt++;
		}
	    }
	    ncols++;
	}
    }

    if (colnames) {
	ILL_SAFE_MALLOC (info->colnames, info->colsize, char *);
	if (!info->colnames) {
	    fprintf (stderr, "out of memory in grab_lp\n");
	    rval = 1;
	    goto CLEANUP;
	}
	for (j = 0; j < info->colsize; j++) {
	    info->colnames[j] = 0;
	}

	ILL_SAFE_MALLOC (buf, ILL_namebufsize, char);
	if (!buf) {
	    fprintf (stderr, "out of memory in grab_lp\n");
	    rval = 1;
	    goto CLEANUP;
	}
	ncols = 0;
	for (j = 0; j < G->ncols; j++) {
	    if (gcols[j].del == 0) {
		if (gcols[j].coltype == dbl_ILL_PRE_COL_STRUC) {
		    len = strlen (colnames[j]) + 1;
		    ILL_SAFE_MALLOC (info->colnames[ncols], len, char);
		    if (!info->colnames[ncols]) {
			fprintf (stderr, "out of memory in grab_lp\n");
			rval = 1;
			goto CLEANUP;
		    }
		    strcpy (info->colnames[ncols], colnames[j]);
		} else {
		    for (k = 0; k < gcols[j].deg; k++) {
			if (gcols[j].adj[k]->del == 0) {
			    i = gcols[j].adj[k]->row;
			    break;
			}
		    }
		    if (k == gcols[j].deg) {
			fprintf (stderr, "problem with dbl_graph in grab_lp\n");
			rval = 1;
			goto CLEANUP;
		    }
		    sprintf (buf, "s%d", i);
		    len = strlen (buf) + 1;
		    ILL_SAFE_MALLOC (info->colnames[ncols], len, char);
		    if (!info->colnames[ncols]) {
			fprintf (stderr, "out of memory in grab_lp\n");
			rval = 1;
			goto CLEANUP;
		    }
		    strcpy (info->colnames[ncols], buf);
		}
		ncols++;
	    }
	}
    }
    /* dbl_ADD STRUCT VARIABLE STUFF */


CLEANUP:

    if (rval) {
	dbl_ILLlp_sinfo_free (info);
    }
    ILL_IFFREE (tdeg, int);
    ILL_IFFREE (map, int);
    ILL_IFFREE (buf, char);

    ILL_RETURN (rval, "dbl_grab_lp_info");
}

static int dbl_fixed_variables (dbl_graph * G,
      dbl_ILLlp_predata * pre)
{
    int rval = 0;
    int j;
    int ncols = G->ncols;
    dbl_node *cols = G->cols;
    dbl_ILLlp_preop *op = 0;

    for (j = 0; j < ncols; j++) {
	if (cols[j].del == 0) {
	    if (dbl_EGlpNumIsEqqual (cols[j].lower, cols[j].upper)) {
		rval = dbl_get_next_preop (pre, &op);
		ILL_CLEANUP_IF (rval);

		op->colindex = j;
		op->rowindex = -1;
		op->ptype = dbl_ILL_PRE_DELETE_FIXED_VARIABLE;

		rval = dbl_grab_lp_line (G, op->colindex, &op->line, 1);
		ILL_CLEANUP_IF (rval);
		pre->opcount++;

		dbl_set_fixed_variable (G, j, cols[j].lower);
	    }
	}
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_fixed_variables");
}

static int dbl_empty_columns (dbl_graph * G,
      dbl_ILLlp_predata * pre)
{
    int rval = 0;
    int j, k;
    int ncols = G->ncols;
    dbl_node *cols = G->cols;
    dbl_ILLlp_preop *op = 0;
    double objtmp;
    dbl_EGlpNumInitVar (objtmp);

    for (j = 0; j < ncols; j++) {
	if (cols[j].del == 0) {
	    for (k = 0; k < cols[j].deg; k++) {
		if (cols[j].adj[k]->del == 0)
		    break;
	    }
	    if (k == cols[j].deg) {
		rval = dbl_get_next_preop (pre, &op);
		ILL_CLEANUP_IF (rval);

		op->colindex = j;
		op->rowindex = -1;
		op->ptype = dbl_ILL_PRE_DELETE_EMPTY_COLUMN;

		rval = dbl_grab_lp_line (G, op->colindex, &op->line, 1);
		ILL_CLEANUP_IF (rval);
		pre->opcount++;
		dbl_EGlpNumCopy (objtmp, cols[j].obj);
		if (G->objsense < 0)
		    dbl_EGlpNumSign (objtmp);
		if (dbl_EGlpNumIsEqual (objtmp, dbl_zeroLpNum, dbl_ILL_PRE_FEAS_TOL)) {
		    dbl_set_fixed_variable (G, j, cols[j].lower);
		} else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, objtmp)) {
		    if (dbl_EGlpNumIsEqqual (cols[j].lower, dbl_ILL_MINDOUBLE)) {
			printf ("unbounded prob detected in dbl_empty_columns\n");
			printf ("col %d, obj %g\n", j, dbl_EGlpNumToLf (cols[j].obj));
			fflush (stdout);
			rval = 1;
			goto CLEANUP;
		    } else {
			dbl_set_fixed_variable (G, j, cols[j].lower);
		    }
		} else if (dbl_EGlpNumIsLess (objtmp, dbl_zeroLpNum)) {
		    if (dbl_EGlpNumIsEqqual (cols[j].upper, dbl_ILL_MAXDOUBLE)) {
			printf ("unbounded prob detected in dbl_empty_columns\n");
			printf ("col %d, obj %g\n", j, dbl_EGlpNumToLf (cols[j].obj));
			fflush (stdout);
			rval = 1;
			goto CLEANUP;
		    } else {
			dbl_set_fixed_variable (G, j, cols[j].upper);
		    }
		} else {
		    dbl_set_fixed_variable (G, j, cols[j].lower);
		}
	    }
	}
    }

CLEANUP:

    dbl_EGlpNumClearVar (objtmp);
    ILL_RETURN (rval, "dbl_empty_columns");
}

static int dbl_singleton_rows (dbl_graph * G,
      dbl_ILLlp_predata * pre,
      int *hit)
{
    int rval = 0;
    int rowindex, i, k, h;
    int nrows = G->nrows;
    dbl_node *rows = G->rows;
    dbl_node *cols = G->cols;
    dbl_node *r, *c;
    dbl_edge *pivot, *f;
    dbl_intptr *next, *list = 0;
    int *tdeg = 0;
    double val;
    dbl_ILLlp_preop *op = 0;
    dbl_EGlpNumInitVar (val);

    *hit = 0;
    if (G->nrows == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (tdeg, G->nrows, int);
    if (!tdeg) {
	fprintf (stderr, "out of memory in dbl_singleton_rows\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < nrows; i++) {
	if (rows[i].del == 0) {
	    tdeg[i] = 0;
	    for (k = 0; k < rows[i].deg; k++) {
		if (rows[i].adj[k]->del == 0) {
		    tdeg[i]++;
		}
	    }
	    if (tdeg[i] <= 1) {
		rval = dbl_add_to_list (&G->intptrworld, &list, i);
		ILL_CLEANUP_IF (rval);
	    }
	}
    }

    while (list) {
	(*hit)++;
	rowindex = list->this;
	next = list->next;
	intptrfree (&G->intptrworld, list);
	list = next;

	rval = dbl_get_next_preop (pre, &op);
	ILL_CLEANUP_IF (rval);

	r = &rows[rowindex];

	if (tdeg[rowindex] == 0) {
	    if (dbl_EGlpNumIsNeqZero (r->rhs, dbl_ILL_PRE_FEAS_TOL)) {
		printf ("infeasible row detected in singleton_row\n");
		printf ("empty row with rhs = %g\n", dbl_EGlpNumToLf (r->rhs));
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	    }
	    op->ptype = dbl_ILL_PRE_DELETE_EMPTY_ROW;
	    op->rowindex = rowindex;
	} else {
	    /* Find the "pivot" entry and colum */

	    for (k = 0; k < r->deg; k++) {
		if (r->adj[k]->del == 0)
		    break;
	    }
	    if (k == r->deg) {
		fprintf (stderr, "lost an dbl_edge in dbl_singleton_rows\n");
		rval = 1;
		goto CLEANUP;
	    }
	    pivot = r->adj[k];
	    c = &cols[pivot->col];

	    /* Store data from operation (incluing the col coefs) */

	    op->ptype = dbl_ILL_PRE_DELETE_SINGLETON_ROW;
	    op->rowindex = rowindex;
	    op->colindex = c - cols;
	    dbl_EGlpNumCopy (op->line.rhs, r->rhs);
	    rval = dbl_grab_lp_line (G, op->colindex, &op->line, 1);
	    ILL_CLEANUP_IF (rval);

	    /* Fix the x[c] to its rhs value */
	    /* val = r->rhs / pivot->coef; */
	    dbl_EGlpNumCopyFrac (val, r->rhs, pivot->coef);
	    /* if (val < c->lower - dbl_ILL_PRE_FEAS_TOL || val > c->upper +
	       dbl_ILL_PRE_FEAS_TOL) */
	    if (dbl_EGlpNumIsSumLess (val, dbl_ILL_PRE_FEAS_TOL, c->lower) ||
		dbl_EGlpNumIsSumLess (c->upper, dbl_ILL_PRE_FEAS_TOL, val)) {
		printf ("infeasible bounds detected in singleton_row %d\n", rowindex);
		printf ("lower->%g  upper->%g  val = %g\n",
		    dbl_EGlpNumToLf (c->lower), dbl_EGlpNumToLf (c->upper),
		    dbl_EGlpNumToLf (val));
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	    } else {
		dbl_EGlpNumCopy (c->lower, val);
		dbl_EGlpNumCopy (c->upper, val);
	    }

	    /* Delete x[c] from other rows (and adjust their rhs) */

	    c->del = 1;

	    for (h = 0; h < c->deg; h++) {
		f = c->adj[h];
		if (f->del == 0) {
		    /* rows[f->row].rhs -= (f->coef * c->lower); */
		    dbl_EGlpNumSubInnProdTo (rows[f->row].rhs, f->coef, c->lower);
		    tdeg[f->row]--;
		    if (tdeg[f->row] == 1) {
			if (f == pivot) {
			    fprintf (stderr, "bad pivot element\n");
			    rval = 1;
			    goto CLEANUP;
			}
			rval = dbl_add_to_list (&G->intptrworld, &list, f->row);
			ILL_CLEANUP_IF (rval);
		    }
		    f->del = 1;
		}
	    }
	}

	r->del = 1;
	pre->opcount++;
    }

CLEANUP:

    ILL_IFFREE (tdeg, int);
    intptr_listfree (&G->intptrworld, list);
    dbl_EGlpNumClearVar (val);
    ILL_RETURN (rval, "dbl_singleton_rows");
}

static int dbl_forcing_constraints (dbl_graph * G,
      dbl_ILLlp_predata * pre,
      int *hit)
{
    int rval = 0;
    int i, j, k, ts;
    dbl_node *rows = G->rows;
    dbl_node *cols = G->cols;
    dbl_edge *e;
    int nrows = G->nrows;
    double ub, lb;
    dbl_ILLlp_preop *op = 0;
    dbl_EGlpNumInitVar (ub);
    dbl_EGlpNumInitVar (lb);

    *hit = 0;

    for (i = 0; i < nrows; i++) {
	if (rows[i].del == 0) {
	    dbl_get_implied_rhs_bounds (G, i, &lb, &ub);
	    if (dbl_EGlpNumIsSumLess (rows[i].rhs, dbl_ILL_PRE_FEAS_TOL, lb) ||
		dbl_EGlpNumIsSumLess (ub, dbl_ILL_PRE_FEAS_TOL, rows[i].rhs)) {
		printf ("infeasible row detected in dbl_forcing_constraints\n");
		printf ("Row %d:  RHS->%g  LBnd->%g  UBnd->%g\n",
		    i, dbl_EGlpNumToLf (rows[i].rhs),
		    dbl_EGlpNumToLf (lb), dbl_EGlpNumToLf (ub));
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	    } else if (dbl_EGlpNumIsDiffLess (rows[i].rhs, dbl_ILL_PRE_FEAS_TOL, lb) ||
		dbl_EGlpNumIsDiffLess (ub, dbl_ILL_PRE_FEAS_TOL, rows[i].rhs)) {
		(*hit)++;
		ts = (dbl_EGlpNumIsDiffLess (rows[i].rhs, dbl_ILL_PRE_FEAS_TOL, lb) ? 0 : 1);
		for (k = 0; k < rows[i].deg; k++) {
		    e = rows[i].adj[k];
		    if (e->del == 0) {
			j = e->col;

			rval = dbl_get_next_preop (pre, &op);
			ILL_CLEANUP_IF (rval);

			op->colindex = j;
			op->rowindex = i;
			op->ptype = dbl_ILL_PRE_DELETE_FORCED_VARIABLE;

			rval = dbl_grab_lp_line (G, j, &op->line, 1);
			ILL_CLEANUP_IF (rval);
			pre->opcount++;

			if ((ts == 0 && dbl_EGlpNumIsLess (e->coef, dbl_zeroLpNum)) ||
			    (ts == 1 && dbl_EGlpNumIsLess (dbl_zeroLpNum, e->coef))) {
			    dbl_set_fixed_variable (G, j, cols[j].upper);
			} else {
			    dbl_set_fixed_variable (G, j, cols[j].lower);
			}
		    }
		}
	    }
	}
    }

CLEANUP:

    dbl_EGlpNumClearVar (ub);
    dbl_EGlpNumClearVar (lb);
    ILL_RETURN (rval, "dbl_forcing_constraints");
}

static int dbl_singleton_columns (dbl_graph * G,
      dbl_ILLlp_predata * pre,
      int *hit)
{
    int rval = 0;
    int ncols = G->ncols;
    int j, k, deg, rdeg, single = 0, irow;
    double lb, ub, b, eb;
    dbl_node *cols = G->cols;
    dbl_node *rows = G->rows;
    dbl_edge *b_edge;
    dbl_ILLlp_preop *op = 0;
    double newub, newlb;
    double a, c, l, u;
    dbl_EGlpNumInitVar (lb);
    dbl_EGlpNumInitVar (ub);
    dbl_EGlpNumInitVar (eb);
    dbl_EGlpNumInitVar (b);
    dbl_EGlpNumInitVar (newlb);
    dbl_EGlpNumInitVar (newub);
    dbl_EGlpNumInitVar (a);
    dbl_EGlpNumInitVar (c);
    dbl_EGlpNumInitVar (l);
    dbl_EGlpNumInitVar (u);

    *hit = 0;
    if (G->ncols == 0)
	goto CLEANUP;

    for (j = 0; j < ncols; j++) {
	if (cols[j].del == 0) {
	    deg = 0;
	    for (k = 0; k < cols[j].deg && deg <= 1; k++) {
		if (cols[j].adj[k]->del == 0) {
		    single = k;
		    deg++;
		}
	    }
	    if (deg == 1) {
		irow = cols[j].adj[single]->row;
		dbl_EGlpNumCopy (b, cols[j].adj[single]->coef);
		b_edge = cols[j].adj[single];

		dbl_get_implied_variable_bounds (G, j, b_edge, &lb, &ub);

		/* if (lb >= cols[j].lower && ub <= cols[j].upper) */
		if (dbl_EGlpNumIsLeq (cols[j].lower, lb) &&
		    dbl_EGlpNumIsLeq (ub, cols[j].upper)) {
		    dbl_edge *a_edge;

		    /* The jth variable can be substituted out of problem */
		    /* x = (c/b) - (a/b)y                           */


		    rval = dbl_get_next_preop (pre, &op);
		    ILL_CLEANUP_IF (rval);

		    op->colindex = j;
		    op->rowindex = irow;
		    op->ptype = dbl_ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE;

		    rval = dbl_grab_lp_line (G, irow, &op->line, 0);
		    ILL_CLEANUP_IF (rval);
		    pre->opcount++;

		    /* Adjust the objective function                      */
		    /* dy ==> (d - (e/b))ay   (e is obj coef of y)     */
		    /* eb = cols[j].obj / b; */
		    dbl_EGlpNumCopyFrac (eb, cols[j].obj, b);

		    for (k = 0; k < rows[irow].deg; k++) {
			a_edge = rows[irow].adj[k];
			if (a_edge->del == 0 && a_edge != b_edge) {
			    /* cols[a_edge->col].obj -= (eb * a_edge->coef); */
			    dbl_EGlpNumSubInnProdTo (cols[a_edge->col].obj, eb, a_edge->coef);
			}
		    }


		    /* Delete y from dbl_graph */

		    cols[j].del = 1;

		    /* Delete equation ay + bx = c */

		    rows[irow].del = 1;
		    for (k = 0; k < rows[irow].deg; k++) {
			rows[irow].adj[k]->del = 1;
		    }

		} else {
		    rdeg = 0;
		    for (k = 0; k < rows[irow].deg && rdeg <= 2; k++) {
			if (rows[irow].adj[k]->del == 0) {
			    rdeg++;
			}
		    }
		    if (rdeg == 2) {
			dbl_edge *a_edge = 0;
			int col2 = 0;
			dbl_EGlpNumCopy (newub, dbl_ILL_MAXDOUBLE);
			dbl_EGlpNumCopy (newlb, dbl_ILL_MINDOUBLE);
			dbl_EGlpNumZero (a);

			/* ay + bx = c                                */
			/* l <= x <= u                                */
			/* x - is column singleton                  */
			/* derive bounds on y and substitute out x  */

			dbl_EGlpNumCopy (c, rows[irow].rhs);
			dbl_EGlpNumCopy (l, cols[j].lower);
			dbl_EGlpNumCopy (u, cols[j].upper);

			/* Find the ay term */

			for (k = 0; k < rows[irow].deg; k++) {
			    if (rows[irow].adj[k]->del == 0 && rows[irow].adj[k]->col != j) {
				a_edge = rows[irow].adj[k];
				dbl_EGlpNumCopy (a, rows[irow].adj[k]->coef);
				col2 = rows[irow].adj[k]->col;
				break;
			    }
			}
			if (k == rows[irow].deg) {
			    fprintf (stderr, "dbl_graph error in singleton_col\n");
			    rval = 1;
			    goto CLEANUP;
			}
			/* Record the operation             */
			/* x is column j,  y is column col2 */

			rval = dbl_get_next_preop (pre, &op);
			ILL_CLEANUP_IF (rval);

			op->colindex = j;
			op->rowindex = irow;
			op->ptype = dbl_ILL_PRE_DELETE_SINGLETON_VARIABLE;

			rval = dbl_grab_lp_line (G, irow, &op->line, 0);
			ILL_CLEANUP_IF (rval);
			pre->opcount++;

			/* Adjust the bounds on y           */
			/* Using x = c/b - (a/b)y            */
			/* we use eb as temporal variable here */
			/* if (a / b > 0) */
			dbl_EGlpNumCopyFrac (eb, a, b);
			if (dbl_EGlpNumIsLess (dbl_zeroLpNum, eb)) {
			    /* if (l > -dbl_ILL_MAXDOUBLE) */
			    if (dbl_EGlpNumIsLess (dbl_ILL_MINDOUBLE, l)) {
				/* newub = (c / a) - (l * b) / a; */
				dbl_EGlpNumCopy (newub, c);
				dbl_EGlpNumSubInnProdTo (newub, l, b);
				dbl_EGlpNumDivTo (newub, a);
			    }
			    /* if (u < dbl_ILL_MAXDOUBLE) */
			    if (dbl_EGlpNumIsLess (u, dbl_ILL_MAXDOUBLE)) {
				/* newlb = (c / a) - (u * b) / a; */
				dbl_EGlpNumCopy (newlb, c);
				dbl_EGlpNumSubInnProdTo (newlb, u, b);
				dbl_EGlpNumDivTo (newlb, a);
			    }
			} else {
			    /* if (l > -dbl_ILL_MAXDOUBLE) */
			    if (dbl_EGlpNumIsLess (dbl_ILL_MINDOUBLE, l)) {
				/* newlb = (c / a) - (l * b) / a; */
				dbl_EGlpNumCopy (newlb, c);
				dbl_EGlpNumSubInnProdTo (newlb, l, b);
				dbl_EGlpNumDivTo (newlb, a);
			    }
			    /* if (u < dbl_ILL_MAXDOUBLE) */
			    if (dbl_EGlpNumIsLess (u, dbl_ILL_MAXDOUBLE)) {
				/* newub = (c / a) - (u * b) / a; */
				dbl_EGlpNumCopy (newub, c);
				dbl_EGlpNumSubInnProdTo (newub, u, b);
				dbl_EGlpNumDivTo (newub, a);
			    }
			}

			if (dbl_EGlpNumIsLess (cols[col2].lower, newlb))
			    dbl_EGlpNumCopy (cols[col2].lower, newlb);
			if (dbl_EGlpNumIsLess (newub, cols[col2].upper))
			    dbl_EGlpNumCopy (cols[col2].upper, newub);
			dbl_EGlpNumSubTo (cols[col2].obj, eb);

			/* Delete x (and the bx term) from dbl_graph */

			cols[j].del = 1;
			b_edge->del = 1;

			/* Delete equation ay + bx = c (and the ax term) */

			rows[irow].del = 1;
			a_edge->del = 1;
		    }
		}
	    }
	}
    }


CLEANUP:

    dbl_EGlpNumClearVar (lb);
    dbl_EGlpNumClearVar (ub);
    dbl_EGlpNumClearVar (eb);
    dbl_EGlpNumClearVar (b);
    dbl_EGlpNumClearVar (newlb);
    dbl_EGlpNumClearVar (newub);
    dbl_EGlpNumClearVar (a);
    dbl_EGlpNumClearVar (c);
    dbl_EGlpNumClearVar (l);
    dbl_EGlpNumClearVar (u);
    ILL_RETURN (rval, "dbl_singleton_columns");
}

static int dbl_duplicate_rows (dbl_graph * G,
      int *hit)
{
    int rval = 0;
    dbl_node *cols = G->cols;
    dbl_node *rows = G->rows;
    int ncols = G->ncols;
    int nrows = G->nrows;
    int *s = 0;
    double *f = 0;
    double szeit = ILLutil_zeit ();
    double q;
    int i, j, k, k2, ri, r0 = 0, n, nu = 0, got, t0, t = 1;
    dbl_node *c;
    dbl_EGlpNumInitVar (q);


    /* Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

    *hit = 0;
    if (nrows == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (s, nrows, int);
    f = dbl_EGlpNumAllocArray (nrows);

    for (i = 0; i < nrows; i++) {
	if (rows[i].del || rows[i].rowsense != 'E') {
	    s[i] = dbl_ILL_MAXINT;	/* dbl_ILL_MAXINT means no longer
					   eligible    */
	} else {
	    s[i] = 0;		/* 0 means eligible, >0 means in a group */
	    nu++;		/* Tracks the number of eligible rows    */
	}
    }

    for (j = 0; j < ncols; j++) {
	c = &cols[j];
	if (c->del)
	    continue;
	if (c->coltype != dbl_ILL_PRE_COL_STRUC)
	    continue;

	n = 0;
	t0 = t++;

	for (k = 0; k < c->deg; k++) {
	    if (c->adj[k]->del)
		continue;

	    ri = c->adj[k]->row;
	    if (s[ri] == 0) {
		s[ri] = t0;
		dbl_EGlpNumCopy (f[ri], c->adj[k]->coef);
		r0 = ri;
		n++;
	    } else if (s[ri] < t0) {
		got = 0;
		for (k2 = k + 1; k2 < c->deg; k2++) {
		    if (c->adj[k2]->del)
			continue;

		    i = c->adj[k2]->row;
		    if (s[i] == s[ri]) {
			/* q = (c->adj[k]->coef * (f[i])) / (f[ri] *
			   (c->adj[k2]->coef)); */
			dbl_EGlpNumCopy (q, c->adj[k]->coef);
			dbl_EGlpNumMultTo (q, f[i]);
			dbl_EGlpNumDivTo (q, f[ri]);
			dbl_EGlpNumDivTo (q, c->adj[k2]->coef);
			if (dbl_EGlpNumIsEqual (q, dbl_oneLpNum, dbl_ILL_PRE_ZERO_TOL)) {
			    s[ri] = t;
			    s[i] = t;
			    got++;
			}
		    }
		}
		if (got) {
		    t++;
		} else {
		    s[ri] = dbl_ILL_MAXINT;
		    if (--nu == 0)
			goto DONE;
		}
	    }
	}

	if (n == 1) {
	    s[r0] = dbl_ILL_MAXINT;
	    if (--nu == 0)
		goto DONE;
	}
    }

DONE:

    {
	int idup = 0;

	for (i = 0; i < nrows; i++) {
	    if (s[i] > 0 && s[i] < dbl_ILL_MAXINT) {
		printf ("Row %d: %d\n", i, s[i]);
		idup++;
	    }
	}
	printf ("Number of duplicate rows: %d\n", idup);
    }

    printf ("Time in dbl_duplicate_rows: %.2f (seconds)\n", ILLutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    ILL_IFFREE (s, int);
    dbl_EGlpNumFreeArray (f);
    dbl_EGlpNumClearVar (q);
    ILL_RETURN (rval, "dbl_duplicate_rows");
}

static int dbl_duplicate_cols (dbl_graph * G,
      int *hit)
{
    int rval = 0;
    dbl_node *cols = G->cols;
    dbl_node *rows = G->rows;
    int ncols = G->ncols;
    int nrows = G->nrows;
    int *s = 0;
    double *f = 0;
    double szeit = ILLutil_zeit ();
    double q;
    int i, j, k, k2, ci, c0 = 0, n, nu = 0, got, t0, t = 1;
    dbl_node *r;
    dbl_EGlpNumInitVar (q);


    /* Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

    *hit = 0;
    if (ncols == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (s, ncols, int);
    f = dbl_EGlpNumAllocArray (ncols);

    for (j = 0; j < ncols; j++) {
	if (cols[j].del || cols[j].coltype != dbl_ILL_PRE_COL_STRUC) {
	    s[j] = dbl_ILL_MAXINT;	/* dbl_ILL_MAXINT means no longer
					   eligible    */
	} else {
	    s[j] = 0;		/* 0 means eligible, >0 means in a group */
	    nu++;		/* Tracks the number of eligible rows    */
	}
    }

    for (i = 0; i < nrows; i++) {
	r = &rows[i];
	if (r->del)
	    continue;

	n = 0;
	t0 = t++;

	for (k = 0; k < r->deg; k++) {
	    if (r->adj[k]->del)
		continue;

	    ci = r->adj[k]->col;
	    if (s[ci] == 0) {
		s[ci] = t0;
		dbl_EGlpNumCopy (f[ci], r->adj[k]->coef);
		c0 = ci;
		n++;
	    } else if (s[ci] < t0) {
		got = 0;
		for (k2 = k + 1; k2 < r->deg; k2++) {
		    if (r->adj[k2]->del)
			continue;

		    j = r->adj[k2]->col;
		    if (s[j] == s[ci]) {
			/* q = (r->adj[k]->coef * (f[j])) / (f[ci] *
			   (r->adj[k2]->coef)); */
			dbl_EGlpNumCopy (q, r->adj[k]->coef);
			dbl_EGlpNumMultTo (q, f[j]);
			dbl_EGlpNumDivTo (q, f[ci]);
			dbl_EGlpNumDivTo (q, r->adj[k2]->coef);
			if (dbl_EGlpNumIsEqual (q, dbl_oneLpNum, dbl_ILL_PRE_ZERO_TOL)) {
			    s[ci] = t;
			    s[j] = t;
			    got++;
			}
		    }
		}
		if (got) {
		    t++;
		} else {
		    s[ci] = dbl_ILL_MAXINT;
		    if (--nu == 0)
			goto DONE;
		}
	    }
	}

	if (n == 1) {
	    s[c0] = dbl_ILL_MAXINT;
	    if (--nu == 0)
		goto DONE;
	}
    }

DONE:

    {
	int dcount;
	int *dcnt;
	int *dlist;

	rval = dbl_gather_dup_lists (s, ncols, &dcount, &dcnt, &dlist);
	ILL_CLEANUP_IF (rval);
    }

    printf ("Time in dbl_duplicate_cols: %.2f (seconds)\n", ILLutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    ILL_IFFREE (s, int);
    dbl_EGlpNumFreeArray (f);
    dbl_EGlpNumClearVar (q);
    ILL_RETURN (rval, "dbl_duplicate_cols");
}

static int dbl_gather_dup_lists ( /* dbl_graph *G, */ int *s,	/* double *f, */

      int count,
      int *duptotal,
      int **dupcnt,
      int **dupind)
{
    int rval = 0;
    int *cnt = 0;
    int *ind = 0;
    int *beg = 0;
    int i, smax = 0, ndup = 0, total = 0;

    *duptotal = 0;
    *dupcnt = 0;
    *dupind = 0;


    for (i = 0; i < count; i++) {
	if (s[i] < dbl_ILL_MAXINT && s[i] > smax)
	    smax = s[i];
    }
    if (smax == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (cnt, smax + 1, int);

    ILL_SAFE_MALLOC (ind, smax + 1, int);

    for (i = 0; i < smax + 1; i++) {
	cnt[i] = 0;
    }

    for (i = 0; i < count; i++) {
	if (s[i] < dbl_ILL_MAXINT) {
	    cnt[s[i]]++;
	}
    }

    if (cnt[0] > 0)
	printf ("%d Empty Lines\n", cnt[0]);

    printf ("Duplicate Classes:");
    fflush (stdout);
    for (i = 1; i < smax + 1; i++) {
	if (cnt[i] > 1) {
	    ndup++;
	    printf (" %d", cnt[i]);
	}
    }
    printf ("  Number %d\n", ndup);
    fflush (stdout);

    if (ndup == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (beg, ndup, int);

    for (i = 1, ndup = 0; i < smax + 1; i++) {
	if (cnt[i] > 1) {
	    beg[ndup] = total;
	    total += cnt[i];
	    ind[i] = ndup;
	    ndup++;
	}
    }

    if (total == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (*dupcnt, ndup, int);

    ILL_SAFE_MALLOC (*dupind, total, int);

    for (i = 0; i < ndup; i++) {
	(*dupcnt)[i] = 0;
    }

    for (i = 0; i < count; i++) {
	if (s[i] < dbl_ILL_MAXINT && s[i] > 0) {
	    if (cnt[s[i]] > 1) {
		(*dupind)[beg[ind[s[i]]] + (*dupcnt)[ind[s[i]]]] = i;
		(*dupcnt)[ind[s[i]]]++;
	    }
	}
    }

    for (i = 0; i < ndup; i++) {
	int j;
	for (j = beg[i]; j < beg[i] + (*dupcnt)[i]; j++) {
	    printf (" %d", (*dupind)[j]);
	}
	printf (" | ");
	fflush (stdout);
    }

    *duptotal = ndup;

CLEANUP:

    ILL_IFFREE (cnt, int);
    ILL_IFFREE (ind, int);
    ILL_IFFREE (beg, int);

    ILL_RETURN (rval, "dbl_gather_dup_lists");
}

static void dbl_set_fixed_variable (dbl_graph * G,
      int j,
      double val)
{
    int k;
    dbl_edge *e;

    G->cols[j].del = 1;
    for (k = 0; k < G->cols[j].deg; k++) {
	e = G->cols[j].adj[k];
	if (e->del == 0) {
	    /* G->rows[e->row].rhs -= (e->coef * val); */
	    dbl_EGlpNumSubInnProdTo (G->rows[e->row].rhs, e->coef, val);
	    e->del = 1;
	}
    }
}

static void dbl_get_implied_rhs_bounds (dbl_graph * G,
      int i,
      double *lb,
      double *ub)
{
    int k;
    double l, u;
    dbl_node *cols = G->cols;
    dbl_node *rows = G->rows;
    dbl_edge *e;
    dbl_EGlpNumInitVar (u);
    dbl_EGlpNumInitVar (l);

    dbl_EGlpNumZero (l);
    for (k = 0; k < rows[i].deg; k++) {
	e = rows[i].adj[k];
	if (e->del == 0) {
	    if (dbl_EGlpNumIsLess (e->coef, dbl_zeroLpNum)) {
		if (dbl_EGlpNumIsEqqual (cols[e->col].upper, dbl_ILL_MAXDOUBLE)) {
		    dbl_EGlpNumCopy (l, dbl_ILL_MINDOUBLE);
		    break;
		} else {
		    /* l += (e->coef * cols[e->col].upper); */
		    dbl_EGlpNumAddInnProdTo (l, e->coef, cols[e->col].upper);
		}
	    } else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, e->coef)) {
		if (dbl_EGlpNumIsEqqual (cols[e->col].lower, dbl_ILL_MINDOUBLE)) {
		    dbl_EGlpNumCopy (l, dbl_ILL_MINDOUBLE);
		    break;
		} else {
		    /* l += (e->coef * cols[e->col].lower); */
		    dbl_EGlpNumAddInnProdTo (l, e->coef, cols[e->col].lower);
		}
	    }
	}
    }

    dbl_EGlpNumZero (u);
    for (k = 0; k < rows[i].deg; k++) {
	e = rows[i].adj[k];
	if (e->del == 0) {
	    if (dbl_EGlpNumIsLess (e->coef, dbl_zeroLpNum)) {
		if (dbl_EGlpNumIsEqqual (cols[e->col].lower, dbl_ILL_MINDOUBLE)) {
		    dbl_EGlpNumCopy (u, dbl_ILL_MAXDOUBLE);
		} else {
		    /* u += (e->coef * cols[e->col].lower); */
		    dbl_EGlpNumAddInnProdTo (u, e->coef, cols[e->col].lower);
		}
	    } else if (dbl_EGlpNumIsLess (dbl_zeroLpNum, e->coef)) {
		if (dbl_EGlpNumIsEqqual (cols[e->col].upper, dbl_ILL_MAXDOUBLE)) {
		    dbl_EGlpNumCopy (u, dbl_ILL_MAXDOUBLE);
		} else {
		    /* u += (e->coef * cols[e->col].upper); */
		    dbl_EGlpNumAddInnProdTo (u, e->coef, cols[e->col].upper);
		}
	    }
	}
    }

    dbl_EGlpNumCopy (*lb, l);
    dbl_EGlpNumCopy (*ub, u);
    dbl_EGlpNumClearVar (u);
    dbl_EGlpNumClearVar (l);
}

static void dbl_get_implied_variable_bounds (dbl_graph * G,
      int j,
      dbl_edge * a_ij,
      double *lb,
      double *ub)
{
    int i = a_ij->row;
    double l, u;
    dbl_EGlpNumInitVar (u);
    dbl_EGlpNumInitVar (l);

    dbl_get_implied_rhs_bounds (G, i, &l, &u);
    dbl_EGlpNumCopy (*lb, dbl_ILL_MINDOUBLE);
    dbl_EGlpNumCopy (*ub, dbl_ILL_MAXDOUBLE);

    if (dbl_EGlpNumIsLess (dbl_ILL_PRE_FEAS_TOL, a_ij->coef)) {
	if (dbl_EGlpNumIsLess (u, dbl_ILL_MAXDOUBLE)) {
	    /**lb = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].upper;*/
	    dbl_EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, u, a_ij->coef);
	    dbl_EGlpNumAddTo (*lb, G->cols[j].upper);
	}
	if (dbl_EGlpNumIsLess (dbl_ILL_MINDOUBLE, l)) {
	    /**ub = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].lower;*/
	    dbl_EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, l, a_ij->coef);
	    dbl_EGlpNumAddTo (*ub, G->cols[j].lower);
	}
    } else if (dbl_EGlpNumIsLess (a_ij->coef, dbl_ILL_PRE_FEAS_TOL)) {
	if (dbl_EGlpNumIsLess (dbl_ILL_MINDOUBLE, l)) {
	    /**lb = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].upper;*/
	    dbl_EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, l, a_ij->coef);
	    dbl_EGlpNumAddTo (*lb, G->cols[j].upper);
	}
	if (dbl_EGlpNumIsLess (u, dbl_ILL_MAXDOUBLE)) {
	    /**ub = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].lower;*/
	    dbl_EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, u, a_ij->coef);
	    dbl_EGlpNumAddTo (*ub, G->cols[j].lower);
	}
    }
    dbl_EGlpNumClearVar (u);
    dbl_EGlpNumClearVar (l);
}

static int dbl_get_next_preop (dbl_ILLlp_predata * pre,
      dbl_ILLlp_preop ** op)
{
    int rval = 0;

    if (pre->opcount >= pre->opsize) {
	pre->opsize *= 1.3;
	pre->opsize += 1000;
	if (pre->opsize < pre->opcount + 1)
	    pre->opsize = pre->opcount + 1;
	pre->oplist = EGrealloc (pre->oplist, sizeof (dbl_ILLlp_preop) * pre->opsize);
	/*
			rval = ILLutil_reallocrus_scale ((void **) &pre->oplist, &pre->opsize, pre->opcount + 1, 1.3, sizeof (dbl_ILLlp_preop));
			ILL_CLEANUP_IF (rval);
	*/
    }
    *op = &pre->oplist[pre->opcount];
    dbl_ILLlp_preop_init (*op);

    /* CLEANUP: */

    ILL_RETURN (rval, "dbl_get_next_preop");
}

static int dbl_add_to_list (ILLptrworld * world,
      dbl_intptr ** list,
      int i)
{
    int rval = 0;
    dbl_intptr *ip;

    ip = intptralloc (world);
    if (!ip) {
	rval = 1;
	goto CLEANUP;
    }
    ip->this = i;
    ip->next = *list;
    *list = ip;

CLEANUP:

    ILL_RETURN (rval, "dbl_add_to_list");
}

static int dbl_build_graph (dbl_ILLlpdata * lp,
      dbl_graph * G)
{
    int rval = 0;
    int ncols = lp->ncols;
    int nrows = lp->nrows;
    int nzcount = lp->nzcount;
    int i, j, k, stop, count;
    dbl_edge *edgelist;
    dbl_node *rows, *cols;
    dbl_ILLmatrix *A = &lp->A;

    G->objsense = lp->objsense;

    ILL_SAFE_MALLOC (G->rows, nrows, dbl_node);
    if (!G->rows) {
	fprintf (stderr, "out of memory in dbl_build_graph\n");
	rval = 1;
	goto CLEANUP;
    }
    rows = G->rows;

    for (i = 0; i < nrows; i++) {
	rows[i].rowsense = lp->sense[i];
	rows[i].deg = 0;
    }

    ILL_SAFE_MALLOC (G->cols, ncols, dbl_node);
    ILL_SAFE_MALLOC (G->edgelist, nzcount, dbl_edge);
    for (i = nzcount; i--;)
	dbl_EGlpNumInitVar ((G->edgelist[i].coef));
    G->nzcount = nzcount;
    ILL_SAFE_MALLOC (G->adjspace, 2 * nzcount, dbl_edge *);

    if (!G->cols || !G->edgelist || !G->adjspace) {
	fprintf (stderr, "out of memory in dbl_build_graph\n");
	rval = 1;
	goto CLEANUP;
    }
    cols = G->cols;
    edgelist = G->edgelist;

    for (j = 0; j < ncols; j++) {
	stop = A->matbeg[j] + A->matcnt[j];
	for (k = A->matbeg[j]; k < stop; k++) {
	    rows[A->matind[k]].deg++;
	}
    }

    for (i = 0, count = 0; i < nrows; i++) {
	rows[i].adj = G->adjspace + count;
	count += rows[i].deg;
	rows[i].deg = 0;
    }

    for (j = 0; j < ncols; j++) {
	cols[j].adj = G->adjspace + count;
	count += A->matcnt[j];
	cols[j].deg = 0;
	cols[j].coltype = dbl_ILL_PRE_COL_STRUC;
    }
    for (i = 0; i < nrows; i++) {
	cols[lp->rowmap[i]].coltype = dbl_ILL_PRE_COL_LOGICAL;
    }

    for (j = 0, count = 0; j < ncols; j++) {
	dbl_EGlpNumCopy (cols[j].obj, lp->obj[j]);
	dbl_EGlpNumCopy (cols[j].lower, lp->lower[j]);
	dbl_EGlpNumCopy (cols[j].upper, lp->upper[j]);
	stop = A->matbeg[j] + A->matcnt[j];
	for (k = A->matbeg[j]; k < stop; k++) {
	    i = A->matind[k];
	    rows[i].adj[rows[i].deg++] = &(edgelist[count]);
	    cols[j].adj[cols[j].deg++] = &(edgelist[count]);
	    edgelist[count].row = i;
	    edgelist[count].col = j;
	    dbl_EGlpNumCopy (edgelist[count].coef, A->matval[k]);
	    edgelist[count].mark = 0;
	    edgelist[count].del = 0;
	    edgelist[count].coltype = cols[j].coltype;
	    count++;
	}
    }
    if (count != nzcount) {
	fprintf (stderr, "counts are off in dbl_build_graph\n");
	rval = 1;
	goto CLEANUP;
    }
    G->ecount = count;
    G->nrows = nrows;
    G->ncols = ncols;

    for (i = 0; i < G->nrows; i++) {
	G->rows[i].del = 0;
	dbl_EGlpNumCopy (G->rows[i].rhs, lp->rhs[i]);
    }
    for (j = 0; j < G->ncols; j++) {
	G->cols[j].del = 0;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_build_graph");
}

static void dbl_dump_graph (dbl_graph * G)
{
    int i, j, k;

    printf ("ecount = %d, nrows = %d, ncols = %d\n",
	G->ecount, G->nrows, G->ncols);
    fflush (stdout);

    for (i = 0; i < G->nrows; i++) {
	printf ("Row %d:", i);
	for (k = 0; k < G->rows[i].deg; k++) {
	    printf (" %d", G->rows[i].adj[k]->col);
	    if (G->rows[i].adj[k]->coltype == dbl_ILL_PRE_COL_LOGICAL)
		printf ("S");
	    printf ("(%g)", dbl_EGlpNumToLf (G->rows[i].adj[k]->coef));
	}
	printf ("  rhs: %g", dbl_EGlpNumToLf (G->rows[i].rhs));
	if (G->rows[i].del) {
	    printf (" (deleted)\n");
	} else {
	    printf ("\n");
	}
    }

    for (j = 0; j < G->ncols; j++) {
	if (G->cols[j].coltype == dbl_ILL_PRE_COL_LOGICAL) {
	    printf ("Slk %d:", j);
	} else {
	    printf ("Col %d:", j);
	}
	for (k = 0; k < G->cols[j].deg; k++) {
	    printf (" %d", G->cols[j].adj[k]->row);
	}
	printf ("  obj: %g  bnd: (%g, %g)", dbl_EGlpNumToLf (G->cols[j].obj),
	    dbl_EGlpNumToLf (G->cols[j].lower), dbl_EGlpNumToLf (G->cols[j].upper));
	if (G->cols[j].del) {
	    printf (" (deleted)\n");
	} else {
	    printf ("\n");
	}
    }
}

static void dbl_init_graph (dbl_graph * G)
{
    if (G) {
	G->edgelist = 0;
	G->rows = 0;
	G->cols = 0;
	G->ecount = 0;
	G->nrows = 0;
	G->ncols = 0;
	G->adjspace = 0;
	ILLptrworld_init (&G->intptrworld);
    }
}

static void dbl_free_graph (dbl_graph * G)
{
    register int i;
    if (G) {
	int total, onlist;
	for (i = G->nzcount; i--;)
	    dbl_EGlpNumClearVar ((G->edgelist[i].coef));
	ILL_IFFREE (G->edgelist, dbl_edge);
	ILL_IFFREE (G->rows, dbl_node);
	ILL_IFFREE (G->cols, dbl_node);
	ILL_IFFREE (G->adjspace, dbl_edge *);
	if (intptr_check_leaks (&G->intptrworld, &total, &onlist)) {
	    fprintf (stderr, "WARNING: %d outstanding intptrs\n", total - onlist);
	}
	ILLptrworld_delete (&G->intptrworld);
	dbl_init_graph (G);
    }
}

int dbl_ILLlp_sinfo_print (dbl_ILLlp_sinfo * s)
{
    int rval = 0;
    int i;
    dbl_ILLlpdata lp;
    char *sense = 0;

    dbl_ILLlpdata_init (&lp);

    lp.nrows = s->nrows;
    lp.ncols = s->ncols;
    lp.nzcount = s->nzcount;
    lp.objsense = s->objsense;
    lp.obj = s->obj;
    lp.rhs = s->rhs;
    lp.lower = s->lower;
    lp.upper = s->upper;
    lp.A.matval = s->A.matval;
    lp.A.matcnt = s->A.matcnt;
    lp.A.matbeg = s->A.matbeg;
    lp.A.matind = s->A.matind;
    lp.rownames = 0;
    lp.colnames = s->colnames;
    lp.objname = 0;
    lp.probname = 0;
    lp.intmarker = 0;

    ILL_SAFE_MALLOC (sense, s->nrows, char);
    if (!sense) {
	fprintf (stderr, "out of memory in dbl_ILLlp_sinfo_print\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < s->nrows; i++) {
	sense[i] = 'E';
    }
    lp.sense = sense;

    /*
        rval = ILLlpdata_writelp (&lp, 0);
        ILL_CLEANUP_IF (rval);
    */

CLEANUP:

    ILL_IFFREE (sense, char);
    ILL_RETURN (rval, "dbl_ILLlp_sinfo_print");
}

void dbl_ILLlp_sinfo_init (dbl_ILLlp_sinfo * sinfo)
{
    if (sinfo) {
	sinfo->ncols = 0;
	sinfo->nrows = 0;
	sinfo->nzcount = 0;
	sinfo->rowsize = 0;
	sinfo->colsize = 0;
	sinfo->obj = 0;
	sinfo->rhs = 0;
	sinfo->lower = 0;
	sinfo->upper = 0;
	sinfo->colnames = 0;
	sinfo->objsense = dbl_ILL_MIN;
	dbl_ILLmatrix_init (&sinfo->A);
    }
}

void dbl_ILLlp_sinfo_free (dbl_ILLlp_sinfo * sinfo)
{
    if (sinfo) {
	dbl_EGlpNumFreeArray (sinfo->obj);
	dbl_EGlpNumFreeArray (sinfo->lower);
	dbl_EGlpNumFreeArray (sinfo->upper);
	dbl_EGlpNumFreeArray (sinfo->rhs);
	dbl_ILLmatrix_free (&sinfo->A);
	if (sinfo->colnames) {
	    int i;
	    for (i = 0; i < sinfo->ncols; i++) {
		ILL_IFFREE (sinfo->colnames[i], char);
	    }
	    ILL_IFFREE (sinfo->colnames, char *);
	}
	dbl_ILLlp_sinfo_init (sinfo);
    }
}

void dbl_ILLlp_predata_init (dbl_ILLlp_predata * pre)
{
    if (pre) {
	pre->opcount = 0;
	pre->opsize = 0;
	pre->oplist = 0;
	pre->r_nrows = 0;
	pre->r_ncols = 0;
	pre->colmap = 0;
	pre->rowmap = 0;
	pre->colscale = 0;
	pre->rowscale = 0;
	pre->colfixval = 0;
	pre->rowfixval = 0;
    }
}

void dbl_ILLlp_predata_free (dbl_ILLlp_predata * pre)
{
    if (pre) {
	int i;

	for (i = 0; i < pre->opcount; i++) {
	    dbl_ILLlp_preop_free (&pre->oplist[i]);
	}
	ILL_IFFREE (pre->oplist, dbl_ILLlp_preop);
	ILL_IFFREE (pre->colmap, int);
	ILL_IFFREE (pre->rowmap, int);
	ILL_IFFREE (pre->colscale, double);
	ILL_IFFREE (pre->rowscale, double);
	ILL_IFFREE (pre->colfixval, double);
	ILL_IFFREE (pre->rowfixval, double);
	dbl_ILLlp_predata_init (pre);
    }
}

void dbl_ILLlp_preop_init (dbl_ILLlp_preop * op)
{
    if (op) {
	op->ptype = 0;
	op->rowindex = -1;
	op->colindex = -1;
	dbl_ILLlp_preline_init (&op->line);
    }
}

void dbl_ILLlp_preop_free (dbl_ILLlp_preop * op)
{
    if (op) {
	dbl_ILLlp_preline_free (&op->line);
	dbl_ILLlp_preop_init (op);
    }
}

void dbl_ILLlp_preline_init (dbl_ILLlp_preline * line)
{
    if (line) {
	dbl_EGlpNumInitVar (line->rhs);
	dbl_EGlpNumInitVar (line->obj);
	dbl_EGlpNumInitVar (line->upper);
	dbl_EGlpNumInitVar (line->lower);
	dbl_EGlpNumZero (line->rhs);
	dbl_EGlpNumZero (line->obj);
	dbl_EGlpNumZero (line->upper);
	dbl_EGlpNumZero (line->lower);
	line->count = 0;
	line->ind = 0;
	line->val = 0;
    }
}

void dbl_ILLlp_preline_free (dbl_ILLlp_preline * line)
{
    if (line) {
	dbl_EGlpNumClearVar (line->rhs);
	dbl_EGlpNumClearVar (line->obj);
	dbl_EGlpNumClearVar (line->upper);
	dbl_EGlpNumClearVar (line->lower);
	ILL_IFFREE (line->ind, int);
	dbl_EGlpNumFreeArray (line->val);
	/* dbl_ILLlp_preline_init (line); */
    }
}
