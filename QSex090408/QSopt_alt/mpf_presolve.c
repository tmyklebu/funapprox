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

/* RCS_INFO = "$RCSfile: mpf_presolve.c,v $ $Revision: 1.2 $ $Date:
   2003/11/05 16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* Presolve Routine for Simplex Method                      */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int mpf_ILLlp_add_logicals (ILLlpata *lp)                                 */
/* int mpf_ILLlp_presolve (mpf_ILLlpdata *lp)                                    */
/* int mpf_ILLlp_scale (mpf_ILLlpdata *lp)                                       */
/* void mpf_ILLlp_sinfo_init (mpf_ILLlp_sinfo *sinfo)                            */
/* void mpf_ILLlp_sinfo_free (mpf_ILLlp_sinfo *sinfo)                            */
/* void mpf_ILLlp_predata_init (mpf_ILLlp_predata *pre)                          */
/* void mpf_ILLlp_predata_free (mpf_ILLlp_predata *pre)                          */
/* */
/* NOTES                                                                   */
/* */
/* presolve will assume that logicals have been added.                   */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdata.h"
#include "mpf_lpdefs.h"
/* extern mpf_t mpf_SZERO_TOLER; */

#define mpf_ILL_LP_STATUS_OK (0)
#define mpf_ILL_PRE_FEAS_TOL mpf_PFEAS_TOLER	/* (1e-6) */
#define mpf_ILL_PRE_ZERO_TOL mpf_PIVOT_TOLER	/* (1e-10) */

#define mpf_ILL_PRE_DELETE_EMPTY_ROW               (1)
#define mpf_ILL_PRE_DELETE_SINGLETON_ROW           (2)
#define mpf_ILL_PRE_DELETE_FIXED_VARIABLE          (3)
#define mpf_ILL_PRE_DELETE_FORCED_VARIABLE         (4)
#define mpf_ILL_PRE_DELETE_SINGLETON_VARIABLE      (5)
#define mpf_ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE (6)
#define mpf_ILL_PRE_DELETE_EMPTY_COLUMN            (7)

#define mpf_ILL_PRE_COL_STRUC                      (0)
#define mpf_ILL_PRE_COL_LOGICAL                    (1)

static int mpf_debug = 0;

typedef struct {
    int row;
    int col;
    char coltype;
    char mark;
    char del;
    mpf_t coef;
}
  mpf_edge;

typedef struct mpf_node {
    mpf_edge **adj;
    mpf_t obj;
    mpf_t lower;
    mpf_t upper;
    mpf_t rhs;
    int deg;
    char mark;
    char del;
    char coltype;
    char rowsense;
}
  mpf_node;

typedef struct mpf_intptr {
    int this;
    struct mpf_intptr *next;
}
  mpf_intptr;

typedef struct mpf_graph {
    mpf_edge *edgelist;
    struct mpf_node *rows;
    struct mpf_node *cols;
    int ecount;
    int nrows;
    int ncols;
    int nzcount;
    mpf_edge **adjspace;
    ILLptrworld intptrworld;
    int objsense;
}
  mpf_graph;

void mpf_ILLlp_sinfo_init (mpf_ILLlp_sinfo * sinfo),
  mpf_ILLlp_sinfo_free (mpf_ILLlp_sinfo * sinfo),
  mpf_ILLlp_predata_init (mpf_ILLlp_predata * pre),
  mpf_ILLlp_predata_free (mpf_ILLlp_predata * pre),
  mpf_ILLlp_preop_init (mpf_ILLlp_preop * op),
  mpf_ILLlp_preop_free (mpf_ILLlp_preop * op),
  mpf_ILLlp_preline_init (mpf_ILLlp_preline * line),
  mpf_ILLlp_preline_free (mpf_ILLlp_preline * line);

int mpf_ILLlp_sinfo_print (mpf_ILLlp_sinfo * s);

static void mpf_set_fixed_variable (mpf_graph * G,
      int j,
      mpf_t val),
  mpf_get_implied_rhs_bounds (mpf_graph * G,
      int i,
      mpf_t * lb,
      mpf_t * ub),
  mpf_get_implied_variable_bounds (mpf_graph * G,
      int j,
      mpf_edge * a_ij,
      mpf_t * lb,
      mpf_t * ub),
  mpf_dump_line (mpf_ILLlp_preline * line),
  mpf_init_graph (mpf_graph * G),
  mpf_free_graph (mpf_graph * G),
  mpf_dump_graph (mpf_graph * G);

static int mpf_simple_presolve (mpf_ILLlpdata * lp,
      mpf_ILLlp_predata * pre,
      mpf_ILLlp_sinfo * info,
      int pre_types,
      int *status),
  mpf_grab_lp_line (mpf_graph * G,
      int indx,
      mpf_ILLlp_preline * line,
      int row_or_col),
  mpf_grab_lp_info (mpf_graph * G,
      char **colnames,
      mpf_ILLlp_sinfo * info),
  mpf_fixed_variables (mpf_graph * G,
      mpf_ILLlp_predata * pre),
  mpf_empty_columns (mpf_graph * G,
      mpf_ILLlp_predata * pre),
  mpf_singleton_rows (mpf_graph * G,
      mpf_ILLlp_predata * pre,
      int *hit),
  mpf_forcing_constraints (mpf_graph * G,
      mpf_ILLlp_predata * pre,
      int *hit),
  mpf_singleton_columns (mpf_graph * G,
      mpf_ILLlp_predata * pre,
      int *hit),
  mpf_duplicate_rows (mpf_graph * G,
      int *hit),
  mpf_duplicate_cols (mpf_graph * G,
      int *hit),
  mpf_gather_dup_lists (int *s,
      int count,
      int *duptotal,
      int **dupcnt,
      int **dupind),
  mpf_get_next_preop (mpf_ILLlp_predata * pre,
      mpf_ILLlp_preop ** op),
  mpf_add_to_list (ILLptrworld * world,
      mpf_intptr ** list,
      int i),
  mpf_build_graph (mpf_ILLlpdata * lp,
      mpf_graph * G);


ILL_PTRWORLD_ROUTINES (mpf_intptr, intptralloc, intptr_bulkalloc, intptrfree)
ILL_PTRWORLD_LISTFREE_ROUTINE (mpf_intptr, intptr_listfree, intptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (mpf_intptr, intptr_check_leaks, this, int)
int mpf_ILLlp_add_logicals (mpf_ILLlpdata * lp)
{
    int rval = 0;
    int ncols, nrows, nzcount, i, aindex;
    char *sense;
    mpf_ILLmatrix *A;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlp_add_logicals called with a NULL pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    printf ("mpf_ILLlp_add_logicals ...\n");
    fflush (stdout);

    A = &lp->A;
    sense = lp->sense;
    ncols = lp->ncols;
    nrows = lp->nrows;
    nzcount = lp->nzcount;

    if (nrows == 0)
	goto CLEANUP;
    mpf_EGlpNumReallocArray (&(lp->obj), lp->colsize + nrows);
    mpf_EGlpNumReallocArray (&(lp->upper), lp->colsize + nrows);
    mpf_EGlpNumReallocArray (&(lp->lower), lp->colsize + nrows);
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
    mpf_EGlpNumReallocArray (&(A->matval), A->matsize + nrows);

    for (i = 0; i < nrows; i++) {
	A->matind[A->matsize + i] = -1;
    }

    aindex = A->matsize - A->matfree;

    for (i = 0; i < nrows; i++) {
	lp->rowmap[i] = ncols;
	mpf_EGlpNumZero (lp->obj[ncols]);
	A->matcnt[ncols] = 1;
	A->matbeg[ncols] = aindex;
	A->matind[aindex] = i;
	switch (sense[i]) {
	case 'E':		/* Arificial */
	    mpf_EGlpNumZero (lp->lower[ncols]);
	    mpf_EGlpNumZero (lp->upper[ncols]);
	    mpf_EGlpNumOne (A->matval[aindex]);
	    break;
	case 'G':		/* Surplus   */
	    mpf_EGlpNumZero (lp->lower[ncols]);
	    mpf_EGlpNumCopy (lp->upper[ncols], mpf_ILL_MAXDOUBLE);
	    mpf_EGlpNumOne (A->matval[aindex]);
	    mpf_EGlpNumSign (A->matval[aindex]);
	    break;
	case 'L':		/* Slack     */
	    mpf_EGlpNumZero (lp->lower[ncols]);
	    mpf_EGlpNumCopy (lp->upper[ncols], mpf_ILL_MAXDOUBLE);
	    mpf_EGlpNumOne (A->matval[aindex]);
	    break;
	case 'R':		/* Range     */
	    mpf_EGlpNumZero (lp->lower[ncols]);
	    mpf_EGlpNumCopy (lp->upper[ncols], lp->rangeval[i]);
	    mpf_EGlpNumOne (A->matval[aindex]);
	    mpf_EGlpNumSign (A->matval[aindex]);
	    break;
	default:
	    fprintf (stderr, "unknown sense %c in mpf_ILLlp_add_logicals\n", sense[i]);
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
	mpf_ILLlp_rows_clear (lp->rA);
    } else {
	ILL_SAFE_MALLOC (lp->rA, 1, mpf_ILLlp_rows);
    }

    rval = mpf_ILLlp_rows_init (lp->rA, lp, 1);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlp_add_logicals");
}

int mpf_ILLlp_scale (mpf_ILLlpdata * lp)
{
    int rval = 0;
    int i, j, k, col, row, nstruct, start, stop;
    mpf_ILLmatrix *A;
    mpf_t rho;
    mpf_t *gama = 0;
    mpf_EGlpNumInitVar (rho);

    /* Columns - divide by largest absolute value */

    if (!lp) {
	ILL_ERROR (rval, "mpf_ILLlp_scale called with a NULL pointer");
    }
    if (lp->nrows == 0 || lp->ncols == 0)
	goto CLEANUP;

    A = &lp->A;
    nstruct = lp->nstruct;

    for (j = 0; j < nstruct; j++) {
	col = lp->structmap[j];
	mpf_EGlpNumZero (rho);

	start = A->matbeg[col];
	stop = start + A->matcnt[col];

	for (k = start; k < stop; k++) {
	    mpf_EGlpNumSetToMaxAbs (rho, A->matval[k]);
	}

	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, rho)) {
	    for (k = start; k < stop; k++) {
		mpf_EGlpNumDivTo (A->matval[k], rho);
	    }
	    mpf_EGlpNumDivTo (lp->obj[col], rho);
	    if (mpf_EGlpNumIsNeqq (lp->lower[col], mpf_ILL_MINDOUBLE))
		mpf_EGlpNumMultTo (lp->lower[col], rho);
	    if (mpf_EGlpNumIsNeqq (lp->upper[col], mpf_ILL_MAXDOUBLE))
		mpf_EGlpNumMultTo (lp->upper[col], rho);
	}
    }

    gama = mpf_EGlpNumAllocArray (lp->nrows);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumZero (gama[i]);
    }

    for (j = 0; j < nstruct; j++) {
	col = lp->structmap[j];
	start = A->matbeg[col];
	stop = start + A->matcnt[col];

	for (k = start; k < stop; k++) {
	    row = A->matind[k];
	    mpf_EGlpNumSetToMaxAbs (gama[row], A->matval[k]);
	}
    }

    for (j = 0; j < nstruct; j++) {
	col = lp->structmap[j];
	start = A->matbeg[col];
	stop = start + A->matcnt[col];

	for (k = start; k < stop; k++) {
	    row = A->matind[k];
	    if (mpf_EGlpNumIsLess (mpf_zeroLpNum, gama[row])) {
		mpf_EGlpNumDivTo (A->matval[k], gama[row]);
	    }
	}
    }

    for (i = 0; i < lp->nrows; i++) {
	if (mpf_EGlpNumIsLess (mpf_zeroLpNum, gama[i])) {
	    mpf_EGlpNumDivTo (lp->rhs[i], gama[i]);
	    col = lp->rowmap[i];
	    if (mpf_EGlpNumIsNeqq (lp->upper[col], mpf_ILL_MAXDOUBLE)) {
		mpf_EGlpNumDivTo (lp->upper[col], gama[i]);	/* Ranged row */
	    }
	}
    }

    if (lp->rA) {		/* Need to clear the row version of data */
	mpf_ILLlp_rows_clear (lp->rA);
	ILL_IFFREE (lp->rA, mpf_ILLlp_rows);
    }
CLEANUP:

    mpf_EGlpNumClearVar (rho);
    mpf_EGlpNumFreeArray (gama);
    ILL_RETURN (rval, "mpf_ILLlp_scale");
}

int mpf_ILLlp_presolve (mpf_ILLlpdata * lp,
      int pre_types)
{
    int rval = 0;
    int status = mpf_ILL_LP_STATUS_OK;
    mpf_ILLlp_predata *pre = 0;
    mpf_ILLlp_sinfo *info = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlp_presolve called with a NULL pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    /*
        ILLlpdata_writelp (lp, 0);
        printf ("\n"); fflush (stdout);
    */

    ILL_SAFE_MALLOC (pre, 1, mpf_ILLlp_predata);
    mpf_ILLlp_predata_init (pre);

    ILL_SAFE_MALLOC (info, 1, mpf_ILLlp_sinfo);
    mpf_ILLlp_sinfo_init (info);

    rval = mpf_simple_presolve (lp, pre, info, pre_types, &status);
    ILL_CLEANUP_IF (rval);
    if (status != mpf_ILL_LP_STATUS_OK) {
	printf ("mpf_simple_presolve returned with bad status\n");
	rval = 1;
	goto CLEANUP;
    }
    /*
        rval = mpf_ILLlp_sinfo_print (info);
        ILL_CLEANUP_IF (rval);
    */

CLEANUP:

    if (rval) {
	if (pre) {
	    mpf_ILLlp_predata_free (pre);
	    ILL_IFFREE (pre, mpf_ILLlp_predata);
	}
	if (info) {
	    mpf_ILLlp_sinfo_free (info);
	    ILL_IFFREE (info, mpf_ILLlp_sinfo);
	}
    } else {
	lp->presolve = pre;
	lp->sinfo = info;
    }

    ILL_RETURN (rval, "mpf_ILLlp_presolve");
}


#if 0
int ILLlp_presolve_addrow (mpf_lpinfo * lp,
      int cnt,
      int *ind,
      double *val,
      double rhs)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLlp_sinfo *S;
    mpf_ILLmatrix *A;

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


static int mpf_simple_presolve (mpf_ILLlpdata * lp,
      mpf_ILLlp_predata * pre,
      mpf_ILLlp_sinfo * info,
      int pre_types,
      int *status)
{
    int rval = 0;
    int i, hit, newhit;
    mpf_graph G;

    if (status)
	*status = mpf_ILL_LP_STATUS_OK;
    mpf_init_graph (&G);

    if (!lp) {
	fprintf (stderr, "mpf_simple_presolve called with a NULL pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    printf ("Initial Rows = %d, Cols = %d, Nzcount = %d\n",
	lp->nrows, lp->ncols, lp->nzcount);
    fflush (stdout);

    rval = mpf_build_graph (lp, &G);
    ILL_CLEANUP_IF (rval);
    if (mpf_debug)
	mpf_dump_graph (&G);

    if (pre_types & mpf_ILL_PRE_FIXED) {
	rval = mpf_fixed_variables (&G, pre);
	ILL_CLEANUP_IF (rval);
    }
    do {
	hit = 0;
	if (pre_types & mpf_ILL_PRE_SINGLE_ROW) {
	    rval = mpf_singleton_rows (&G, pre, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & mpf_ILL_PRE_FORCING) {
	    rval = mpf_forcing_constraints (&G, pre, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & mpf_ILL_PRE_SINGLE_COL) {
	    rval = mpf_singleton_columns (&G, pre, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & mpf_ILL_PRE_DUPLICATE_ROW) {
	    rval = mpf_duplicate_rows (&G, &newhit);
	    ILL_CLEANUP_IF (rval);
	    hit += newhit;
	}
	if (pre_types & mpf_ILL_PRE_DUPLICATE_COL) {
	    rval = mpf_duplicate_cols (&G, &newhit);
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

    if (mpf_ILL_PRE_EMPTY_COL) {
	rval = mpf_empty_columns (&G, pre);
	ILL_CLEANUP_IF (rval);
    }
    if (mpf_debug) {
	printf ("Operations\n");
	for (i = 0; i < pre->opcount; i++) {
	    switch (pre->oplist[i].ptype) {
	    case mpf_ILL_PRE_DELETE_EMPTY_ROW:
		printf ("Delete Empty Row: %d\n", pre->oplist[i].rowindex);
		fflush (stdout);
		break;
	    case mpf_ILL_PRE_DELETE_SINGLETON_ROW:
		printf ("Delete Singleton Row: %d (col %d)\n",
		    pre->oplist[i].rowindex, pre->oplist[i].colindex);
		fflush (stdout);
		mpf_dump_line (&pre->oplist[i].line);
		break;
	    case mpf_ILL_PRE_DELETE_FIXED_VARIABLE:
		printf ("Delete Fixed Variable: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		mpf_dump_line (&pre->oplist[i].line);
		break;
	    case mpf_ILL_PRE_DELETE_FORCED_VARIABLE:
		printf ("Delete Forced Variable: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		mpf_dump_line (&pre->oplist[i].line);
		break;
	    case mpf_ILL_PRE_DELETE_SINGLETON_VARIABLE:
		printf ("Delete Singleton Variable: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		mpf_dump_line (&pre->oplist[i].line);
		break;
	    case mpf_ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE:
		printf ("Delete Free Singleton Variable: %d\n",
		    pre->oplist[i].colindex);
		fflush (stdout);
		mpf_dump_line (&pre->oplist[i].line);
		break;
	    case mpf_ILL_PRE_DELETE_EMPTY_COLUMN:
		printf ("Delete Empty Column: %d\n", pre->oplist[i].colindex);
		fflush (stdout);
		mpf_dump_line (&pre->oplist[i].line);
		break;
	    default:
		fprintf (stderr, "unknon presolve operation\n");
		rval = 1;
		goto CLEANUP;
	    }
	}
	printf ("\n");
    }
    rval = mpf_grab_lp_info (&G, lp->colnames, info);
    ILL_CLEANUP_IF (rval);

    /*
        printf ("Final Rows = %d, Cols = %d, Nzcount = %d\n",
                   info->nrows, info->ncols, info->nzcount);
        fflush (stdout);
    */


CLEANUP:

    mpf_free_graph (&G);
    ILL_RETURN (rval, "mpf_simple_presolve");
}

static int mpf_grab_lp_line (mpf_graph * G,
      int indx,
      mpf_ILLlp_preline * line,
      int row_or_col)
{
    int rval = 0;
    int k, cnt;
    mpf_node *n;

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
	line->val = mpf_EGlpNumAllocArray (line->count);
	if (!line->ind || !line->val) {
	    fprintf (stderr, "out of memory in mpf_grab_lp_line\n");
	    rval = 1;
	    goto CLEANUP;
	}
	for (k = 0, cnt = 0; k < n->deg; k++) {
	    if (n->adj[k]->del == 0) {
		line->ind[cnt] = n->adj[k]->row;
		mpf_EGlpNumCopy (line->val[cnt], n->adj[k]->coef);
		cnt++;
	    }
	}
    }
    if (row_or_col == 0) {
	mpf_EGlpNumCopy (line->rhs, n->rhs);
    } else {
	mpf_EGlpNumCopy (line->obj, n->obj);
	mpf_EGlpNumCopy (line->lower, n->lower);
	mpf_EGlpNumCopy (line->upper, n->upper);
    }

    line->row_or_col = row_or_col;

CLEANUP:

    ILL_RETURN (rval, "mpf_grab_lp_line");
}

static void mpf_dump_line (mpf_ILLlp_preline * line)
{
    int k;

    printf (" ");
    if (line->row_or_col == 0) {
	for (k = 0; k < line->count; k++) {
	    printf (" C%d->%g", line->ind[k], mpf_EGlpNumToLf (line->val[k]));
	}
	printf (" RHS->%g\n", mpf_EGlpNumToLf (line->rhs));
    } else {
	for (k = 0; k < line->count; k++) {
	    printf (" R%d->%g", line->ind[k], mpf_EGlpNumToLf (line->val[k]));
	}
	printf (" Obj->%g  LB->%g  UB->%g\n", mpf_EGlpNumToLf (line->obj),
	    mpf_EGlpNumToLf (line->lower), mpf_EGlpNumToLf (line->upper));
    }
    fflush (stdout);
}

static int mpf_grab_lp_info (mpf_graph * G,
      char **colnames,
      mpf_ILLlp_sinfo * info)
{
    int rval = 0;
    int ncols = 0, nrows = 0, nzcount = 0;
    int i, j, k, cnt, len;
    mpf_node *grows = G->rows;
    mpf_node *gcols = G->cols;
    int *tdeg = 0;
    int *map = 0;
    char *buf = 0;
    mpf_ILLmatrix *A = &info->A;

    ILL_SAFE_MALLOC (tdeg, G->ncols, int);
    ILL_SAFE_MALLOC (map, G->nrows, int);
    if (!tdeg || !map) {
	fprintf (stderr, "out of memory in mpf_grab_lp_info\n");
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

    info->rhs = mpf_EGlpNumAllocArray (nrows);
    info->obj = mpf_EGlpNumAllocArray (ncols);
    info->upper = mpf_EGlpNumAllocArray (ncols);
    info->lower = mpf_EGlpNumAllocArray (ncols);
    A->matval = mpf_EGlpNumAllocArray (info->nzcount + 1);
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
	    mpf_EGlpNumCopy (info->rhs[nrows], grows[i].rhs);
	    nrows++;
	}
    }

    ncols = 0;
    cnt = 0;
    for (j = 0; j < G->ncols; j++) {
	if (gcols[j].del == 0) {
	    mpf_EGlpNumCopy (info->obj[ncols], gcols[j].obj);
	    mpf_EGlpNumCopy (info->lower[ncols], gcols[j].lower);
	    mpf_EGlpNumCopy (info->upper[ncols], gcols[j].upper);
	    A->matcnt[ncols] = tdeg[ncols];
	    A->matbeg[ncols] = cnt;
	    for (k = 0; k < gcols[j].deg; k++) {
		if (gcols[j].adj[k]->del == 0) {
		    mpf_EGlpNumCopy (A->matval[cnt], gcols[j].adj[k]->coef);
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
		if (gcols[j].coltype == mpf_ILL_PRE_COL_STRUC) {
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
			fprintf (stderr, "problem with mpf_graph in grab_lp\n");
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
    /* mpf_ADD STRUCT VARIABLE STUFF */


CLEANUP:

    if (rval) {
	mpf_ILLlp_sinfo_free (info);
    }
    ILL_IFFREE (tdeg, int);
    ILL_IFFREE (map, int);
    ILL_IFFREE (buf, char);

    ILL_RETURN (rval, "mpf_grab_lp_info");
}

static int mpf_fixed_variables (mpf_graph * G,
      mpf_ILLlp_predata * pre)
{
    int rval = 0;
    int j;
    int ncols = G->ncols;
    mpf_node *cols = G->cols;
    mpf_ILLlp_preop *op = 0;

    for (j = 0; j < ncols; j++) {
	if (cols[j].del == 0) {
	    if (mpf_EGlpNumIsEqqual (cols[j].lower, cols[j].upper)) {
		rval = mpf_get_next_preop (pre, &op);
		ILL_CLEANUP_IF (rval);

		op->colindex = j;
		op->rowindex = -1;
		op->ptype = mpf_ILL_PRE_DELETE_FIXED_VARIABLE;

		rval = mpf_grab_lp_line (G, op->colindex, &op->line, 1);
		ILL_CLEANUP_IF (rval);
		pre->opcount++;

		mpf_set_fixed_variable (G, j, cols[j].lower);
	    }
	}
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_fixed_variables");
}

static int mpf_empty_columns (mpf_graph * G,
      mpf_ILLlp_predata * pre)
{
    int rval = 0;
    int j, k;
    int ncols = G->ncols;
    mpf_node *cols = G->cols;
    mpf_ILLlp_preop *op = 0;
    mpf_t objtmp;
    mpf_EGlpNumInitVar (objtmp);

    for (j = 0; j < ncols; j++) {
	if (cols[j].del == 0) {
	    for (k = 0; k < cols[j].deg; k++) {
		if (cols[j].adj[k]->del == 0)
		    break;
	    }
	    if (k == cols[j].deg) {
		rval = mpf_get_next_preop (pre, &op);
		ILL_CLEANUP_IF (rval);

		op->colindex = j;
		op->rowindex = -1;
		op->ptype = mpf_ILL_PRE_DELETE_EMPTY_COLUMN;

		rval = mpf_grab_lp_line (G, op->colindex, &op->line, 1);
		ILL_CLEANUP_IF (rval);
		pre->opcount++;
		mpf_EGlpNumCopy (objtmp, cols[j].obj);
		if (G->objsense < 0)
		    mpf_EGlpNumSign (objtmp);
		if (mpf_EGlpNumIsEqual (objtmp, mpf_zeroLpNum, mpf_ILL_PRE_FEAS_TOL)) {
		    mpf_set_fixed_variable (G, j, cols[j].lower);
		} else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, objtmp)) {
		    if (mpf_EGlpNumIsEqqual (cols[j].lower, mpf_ILL_MINDOUBLE)) {
			printf ("unbounded prob detected in mpf_empty_columns\n");
			printf ("col %d, obj %g\n", j, mpf_EGlpNumToLf (cols[j].obj));
			fflush (stdout);
			rval = 1;
			goto CLEANUP;
		    } else {
			mpf_set_fixed_variable (G, j, cols[j].lower);
		    }
		} else if (mpf_EGlpNumIsLess (objtmp, mpf_zeroLpNum)) {
		    if (mpf_EGlpNumIsEqqual (cols[j].upper, mpf_ILL_MAXDOUBLE)) {
			printf ("unbounded prob detected in mpf_empty_columns\n");
			printf ("col %d, obj %g\n", j, mpf_EGlpNumToLf (cols[j].obj));
			fflush (stdout);
			rval = 1;
			goto CLEANUP;
		    } else {
			mpf_set_fixed_variable (G, j, cols[j].upper);
		    }
		} else {
		    mpf_set_fixed_variable (G, j, cols[j].lower);
		}
	    }
	}
    }

CLEANUP:

    mpf_EGlpNumClearVar (objtmp);
    ILL_RETURN (rval, "mpf_empty_columns");
}

static int mpf_singleton_rows (mpf_graph * G,
      mpf_ILLlp_predata * pre,
      int *hit)
{
    int rval = 0;
    int rowindex, i, k, h;
    int nrows = G->nrows;
    mpf_node *rows = G->rows;
    mpf_node *cols = G->cols;
    mpf_node *r, *c;
    mpf_edge *pivot, *f;
    mpf_intptr *next, *list = 0;
    int *tdeg = 0;
    mpf_t val;
    mpf_ILLlp_preop *op = 0;
    mpf_EGlpNumInitVar (val);

    *hit = 0;
    if (G->nrows == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (tdeg, G->nrows, int);
    if (!tdeg) {
	fprintf (stderr, "out of memory in mpf_singleton_rows\n");
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
		rval = mpf_add_to_list (&G->intptrworld, &list, i);
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

	rval = mpf_get_next_preop (pre, &op);
	ILL_CLEANUP_IF (rval);

	r = &rows[rowindex];

	if (tdeg[rowindex] == 0) {
	    if (mpf_EGlpNumIsNeqZero (r->rhs, mpf_ILL_PRE_FEAS_TOL)) {
		printf ("infeasible row detected in singleton_row\n");
		printf ("empty row with rhs = %g\n", mpf_EGlpNumToLf (r->rhs));
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	    }
	    op->ptype = mpf_ILL_PRE_DELETE_EMPTY_ROW;
	    op->rowindex = rowindex;
	} else {
	    /* Find the "pivot" entry and colum */

	    for (k = 0; k < r->deg; k++) {
		if (r->adj[k]->del == 0)
		    break;
	    }
	    if (k == r->deg) {
		fprintf (stderr, "lost an mpf_edge in mpf_singleton_rows\n");
		rval = 1;
		goto CLEANUP;
	    }
	    pivot = r->adj[k];
	    c = &cols[pivot->col];

	    /* Store data from operation (incluing the col coefs) */

	    op->ptype = mpf_ILL_PRE_DELETE_SINGLETON_ROW;
	    op->rowindex = rowindex;
	    op->colindex = c - cols;
	    mpf_EGlpNumCopy (op->line.rhs, r->rhs);
	    rval = mpf_grab_lp_line (G, op->colindex, &op->line, 1);
	    ILL_CLEANUP_IF (rval);

	    /* Fix the x[c] to its rhs value */
	    /* val = r->rhs / pivot->coef; */
	    mpf_EGlpNumCopyFrac (val, r->rhs, pivot->coef);
	    /* if (val < c->lower - mpf_ILL_PRE_FEAS_TOL || val > c->upper +
	       mpf_ILL_PRE_FEAS_TOL) */
	    if (mpf_EGlpNumIsSumLess (val, mpf_ILL_PRE_FEAS_TOL, c->lower) ||
		mpf_EGlpNumIsSumLess (c->upper, mpf_ILL_PRE_FEAS_TOL, val)) {
		printf ("infeasible bounds detected in singleton_row %d\n", rowindex);
		printf ("lower->%g  upper->%g  val = %g\n",
		    mpf_EGlpNumToLf (c->lower), mpf_EGlpNumToLf (c->upper),
		    mpf_EGlpNumToLf (val));
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	    } else {
		mpf_EGlpNumCopy (c->lower, val);
		mpf_EGlpNumCopy (c->upper, val);
	    }

	    /* Delete x[c] from other rows (and adjust their rhs) */

	    c->del = 1;

	    for (h = 0; h < c->deg; h++) {
		f = c->adj[h];
		if (f->del == 0) {
		    /* rows[f->row].rhs -= (f->coef * c->lower); */
		    mpf_EGlpNumSubInnProdTo (rows[f->row].rhs, f->coef, c->lower);
		    tdeg[f->row]--;
		    if (tdeg[f->row] == 1) {
			if (f == pivot) {
			    fprintf (stderr, "bad pivot element\n");
			    rval = 1;
			    goto CLEANUP;
			}
			rval = mpf_add_to_list (&G->intptrworld, &list, f->row);
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
    mpf_EGlpNumClearVar (val);
    ILL_RETURN (rval, "mpf_singleton_rows");
}

static int mpf_forcing_constraints (mpf_graph * G,
      mpf_ILLlp_predata * pre,
      int *hit)
{
    int rval = 0;
    int i, j, k, ts;
    mpf_node *rows = G->rows;
    mpf_node *cols = G->cols;
    mpf_edge *e;
    int nrows = G->nrows;
    mpf_t ub, lb;
    mpf_ILLlp_preop *op = 0;
    mpf_EGlpNumInitVar (ub);
    mpf_EGlpNumInitVar (lb);

    *hit = 0;

    for (i = 0; i < nrows; i++) {
	if (rows[i].del == 0) {
	    mpf_get_implied_rhs_bounds (G, i, &lb, &ub);
	    if (mpf_EGlpNumIsSumLess (rows[i].rhs, mpf_ILL_PRE_FEAS_TOL, lb) ||
		mpf_EGlpNumIsSumLess (ub, mpf_ILL_PRE_FEAS_TOL, rows[i].rhs)) {
		printf ("infeasible row detected in mpf_forcing_constraints\n");
		printf ("Row %d:  RHS->%g  LBnd->%g  UBnd->%g\n",
		    i, mpf_EGlpNumToLf (rows[i].rhs),
		    mpf_EGlpNumToLf (lb), mpf_EGlpNumToLf (ub));
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	    } else if (mpf_EGlpNumIsDiffLess (rows[i].rhs, mpf_ILL_PRE_FEAS_TOL, lb) ||
		mpf_EGlpNumIsDiffLess (ub, mpf_ILL_PRE_FEAS_TOL, rows[i].rhs)) {
		(*hit)++;
		ts = (mpf_EGlpNumIsDiffLess (rows[i].rhs, mpf_ILL_PRE_FEAS_TOL, lb) ? 0 : 1);
		for (k = 0; k < rows[i].deg; k++) {
		    e = rows[i].adj[k];
		    if (e->del == 0) {
			j = e->col;

			rval = mpf_get_next_preop (pre, &op);
			ILL_CLEANUP_IF (rval);

			op->colindex = j;
			op->rowindex = i;
			op->ptype = mpf_ILL_PRE_DELETE_FORCED_VARIABLE;

			rval = mpf_grab_lp_line (G, j, &op->line, 1);
			ILL_CLEANUP_IF (rval);
			pre->opcount++;

			if ((ts == 0 && mpf_EGlpNumIsLess (e->coef, mpf_zeroLpNum)) ||
			    (ts == 1 && mpf_EGlpNumIsLess (mpf_zeroLpNum, e->coef))) {
			    mpf_set_fixed_variable (G, j, cols[j].upper);
			} else {
			    mpf_set_fixed_variable (G, j, cols[j].lower);
			}
		    }
		}
	    }
	}
    }

CLEANUP:

    mpf_EGlpNumClearVar (ub);
    mpf_EGlpNumClearVar (lb);
    ILL_RETURN (rval, "mpf_forcing_constraints");
}

static int mpf_singleton_columns (mpf_graph * G,
      mpf_ILLlp_predata * pre,
      int *hit)
{
    int rval = 0;
    int ncols = G->ncols;
    int j, k, deg, rdeg, single = 0, irow;
    mpf_t lb, ub, b, eb;
    mpf_node *cols = G->cols;
    mpf_node *rows = G->rows;
    mpf_edge *b_edge;
    mpf_ILLlp_preop *op = 0;
    mpf_t newub, newlb;
    mpf_t a, c, l, u;
    mpf_EGlpNumInitVar (lb);
    mpf_EGlpNumInitVar (ub);
    mpf_EGlpNumInitVar (eb);
    mpf_EGlpNumInitVar (b);
    mpf_EGlpNumInitVar (newlb);
    mpf_EGlpNumInitVar (newub);
    mpf_EGlpNumInitVar (a);
    mpf_EGlpNumInitVar (c);
    mpf_EGlpNumInitVar (l);
    mpf_EGlpNumInitVar (u);

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
		mpf_EGlpNumCopy (b, cols[j].adj[single]->coef);
		b_edge = cols[j].adj[single];

		mpf_get_implied_variable_bounds (G, j, b_edge, &lb, &ub);

		/* if (lb >= cols[j].lower && ub <= cols[j].upper) */
		if (mpf_EGlpNumIsLeq (cols[j].lower, lb) &&
		    mpf_EGlpNumIsLeq (ub, cols[j].upper)) {
		    mpf_edge *a_edge;

		    /* The jth variable can be substituted out of problem */
		    /* x = (c/b) - (a/b)y                           */


		    rval = mpf_get_next_preop (pre, &op);
		    ILL_CLEANUP_IF (rval);

		    op->colindex = j;
		    op->rowindex = irow;
		    op->ptype = mpf_ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE;

		    rval = mpf_grab_lp_line (G, irow, &op->line, 0);
		    ILL_CLEANUP_IF (rval);
		    pre->opcount++;

		    /* Adjust the objective function                      */
		    /* dy ==> (d - (e/b))ay   (e is obj coef of y)     */
		    /* eb = cols[j].obj / b; */
		    mpf_EGlpNumCopyFrac (eb, cols[j].obj, b);

		    for (k = 0; k < rows[irow].deg; k++) {
			a_edge = rows[irow].adj[k];
			if (a_edge->del == 0 && a_edge != b_edge) {
			    /* cols[a_edge->col].obj -= (eb * a_edge->coef); */
			    mpf_EGlpNumSubInnProdTo (cols[a_edge->col].obj, eb, a_edge->coef);
			}
		    }


		    /* Delete y from mpf_graph */

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
			mpf_edge *a_edge = 0;
			int col2 = 0;
			mpf_EGlpNumCopy (newub, mpf_ILL_MAXDOUBLE);
			mpf_EGlpNumCopy (newlb, mpf_ILL_MINDOUBLE);
			mpf_EGlpNumZero (a);

			/* ay + bx = c                                */
			/* l <= x <= u                                */
			/* x - is column singleton                  */
			/* derive bounds on y and substitute out x  */

			mpf_EGlpNumCopy (c, rows[irow].rhs);
			mpf_EGlpNumCopy (l, cols[j].lower);
			mpf_EGlpNumCopy (u, cols[j].upper);

			/* Find the ay term */

			for (k = 0; k < rows[irow].deg; k++) {
			    if (rows[irow].adj[k]->del == 0 && rows[irow].adj[k]->col != j) {
				a_edge = rows[irow].adj[k];
				mpf_EGlpNumCopy (a, rows[irow].adj[k]->coef);
				col2 = rows[irow].adj[k]->col;
				break;
			    }
			}
			if (k == rows[irow].deg) {
			    fprintf (stderr, "mpf_graph error in singleton_col\n");
			    rval = 1;
			    goto CLEANUP;
			}
			/* Record the operation             */
			/* x is column j,  y is column col2 */

			rval = mpf_get_next_preop (pre, &op);
			ILL_CLEANUP_IF (rval);

			op->colindex = j;
			op->rowindex = irow;
			op->ptype = mpf_ILL_PRE_DELETE_SINGLETON_VARIABLE;

			rval = mpf_grab_lp_line (G, irow, &op->line, 0);
			ILL_CLEANUP_IF (rval);
			pre->opcount++;

			/* Adjust the bounds on y           */
			/* Using x = c/b - (a/b)y            */
			/* we use eb as temporal variable here */
			/* if (a / b > 0) */
			mpf_EGlpNumCopyFrac (eb, a, b);
			if (mpf_EGlpNumIsLess (mpf_zeroLpNum, eb)) {
			    /* if (l > -mpf_ILL_MAXDOUBLE) */
			    if (mpf_EGlpNumIsLess (mpf_ILL_MINDOUBLE, l)) {
				/* newub = (c / a) - (l * b) / a; */
				mpf_EGlpNumCopy (newub, c);
				mpf_EGlpNumSubInnProdTo (newub, l, b);
				mpf_EGlpNumDivTo (newub, a);
			    }
			    /* if (u < mpf_ILL_MAXDOUBLE) */
			    if (mpf_EGlpNumIsLess (u, mpf_ILL_MAXDOUBLE)) {
				/* newlb = (c / a) - (u * b) / a; */
				mpf_EGlpNumCopy (newlb, c);
				mpf_EGlpNumSubInnProdTo (newlb, u, b);
				mpf_EGlpNumDivTo (newlb, a);
			    }
			} else {
			    /* if (l > -mpf_ILL_MAXDOUBLE) */
			    if (mpf_EGlpNumIsLess (mpf_ILL_MINDOUBLE, l)) {
				/* newlb = (c / a) - (l * b) / a; */
				mpf_EGlpNumCopy (newlb, c);
				mpf_EGlpNumSubInnProdTo (newlb, l, b);
				mpf_EGlpNumDivTo (newlb, a);
			    }
			    /* if (u < mpf_ILL_MAXDOUBLE) */
			    if (mpf_EGlpNumIsLess (u, mpf_ILL_MAXDOUBLE)) {
				/* newub = (c / a) - (u * b) / a; */
				mpf_EGlpNumCopy (newub, c);
				mpf_EGlpNumSubInnProdTo (newub, u, b);
				mpf_EGlpNumDivTo (newub, a);
			    }
			}

			if (mpf_EGlpNumIsLess (cols[col2].lower, newlb))
			    mpf_EGlpNumCopy (cols[col2].lower, newlb);
			if (mpf_EGlpNumIsLess (newub, cols[col2].upper))
			    mpf_EGlpNumCopy (cols[col2].upper, newub);
			mpf_EGlpNumSubTo (cols[col2].obj, eb);

			/* Delete x (and the bx term) from mpf_graph */

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

    mpf_EGlpNumClearVar (lb);
    mpf_EGlpNumClearVar (ub);
    mpf_EGlpNumClearVar (eb);
    mpf_EGlpNumClearVar (b);
    mpf_EGlpNumClearVar (newlb);
    mpf_EGlpNumClearVar (newub);
    mpf_EGlpNumClearVar (a);
    mpf_EGlpNumClearVar (c);
    mpf_EGlpNumClearVar (l);
    mpf_EGlpNumClearVar (u);
    ILL_RETURN (rval, "mpf_singleton_columns");
}

static int mpf_duplicate_rows (mpf_graph * G,
      int *hit)
{
    int rval = 0;
    mpf_node *cols = G->cols;
    mpf_node *rows = G->rows;
    int ncols = G->ncols;
    int nrows = G->nrows;
    int *s = 0;
    mpf_t *f = 0;
    double szeit = ILLutil_zeit ();
    mpf_t q;
    int i, j, k, k2, ri, r0 = 0, n, nu = 0, got, t0, t = 1;
    mpf_node *c;
    mpf_EGlpNumInitVar (q);


    /* Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

    *hit = 0;
    if (nrows == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (s, nrows, int);
    f = mpf_EGlpNumAllocArray (nrows);

    for (i = 0; i < nrows; i++) {
	if (rows[i].del || rows[i].rowsense != 'E') {
	    s[i] = mpf_ILL_MAXINT;	/* mpf_ILL_MAXINT means no longer
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
	if (c->coltype != mpf_ILL_PRE_COL_STRUC)
	    continue;

	n = 0;
	t0 = t++;

	for (k = 0; k < c->deg; k++) {
	    if (c->adj[k]->del)
		continue;

	    ri = c->adj[k]->row;
	    if (s[ri] == 0) {
		s[ri] = t0;
		mpf_EGlpNumCopy (f[ri], c->adj[k]->coef);
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
			mpf_EGlpNumCopy (q, c->adj[k]->coef);
			mpf_EGlpNumMultTo (q, f[i]);
			mpf_EGlpNumDivTo (q, f[ri]);
			mpf_EGlpNumDivTo (q, c->adj[k2]->coef);
			if (mpf_EGlpNumIsEqual (q, mpf_oneLpNum, mpf_ILL_PRE_ZERO_TOL)) {
			    s[ri] = t;
			    s[i] = t;
			    got++;
			}
		    }
		}
		if (got) {
		    t++;
		} else {
		    s[ri] = mpf_ILL_MAXINT;
		    if (--nu == 0)
			goto DONE;
		}
	    }
	}

	if (n == 1) {
	    s[r0] = mpf_ILL_MAXINT;
	    if (--nu == 0)
		goto DONE;
	}
    }

DONE:

    {
	int idup = 0;

	for (i = 0; i < nrows; i++) {
	    if (s[i] > 0 && s[i] < mpf_ILL_MAXINT) {
		printf ("Row %d: %d\n", i, s[i]);
		idup++;
	    }
	}
	printf ("Number of duplicate rows: %d\n", idup);
    }

    printf ("Time in mpf_duplicate_rows: %.2f (seconds)\n", ILLutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    ILL_IFFREE (s, int);
    mpf_EGlpNumFreeArray (f);
    mpf_EGlpNumClearVar (q);
    ILL_RETURN (rval, "mpf_duplicate_rows");
}

static int mpf_duplicate_cols (mpf_graph * G,
      int *hit)
{
    int rval = 0;
    mpf_node *cols = G->cols;
    mpf_node *rows = G->rows;
    int ncols = G->ncols;
    int nrows = G->nrows;
    int *s = 0;
    mpf_t *f = 0;
    double szeit = ILLutil_zeit ();
    mpf_t q;
    int i, j, k, k2, ci, c0 = 0, n, nu = 0, got, t0, t = 1;
    mpf_node *r;
    mpf_EGlpNumInitVar (q);


    /* Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

    *hit = 0;
    if (ncols == 0)
	goto CLEANUP;

    ILL_SAFE_MALLOC (s, ncols, int);
    f = mpf_EGlpNumAllocArray (ncols);

    for (j = 0; j < ncols; j++) {
	if (cols[j].del || cols[j].coltype != mpf_ILL_PRE_COL_STRUC) {
	    s[j] = mpf_ILL_MAXINT;	/* mpf_ILL_MAXINT means no longer
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
		mpf_EGlpNumCopy (f[ci], r->adj[k]->coef);
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
			mpf_EGlpNumCopy (q, r->adj[k]->coef);
			mpf_EGlpNumMultTo (q, f[j]);
			mpf_EGlpNumDivTo (q, f[ci]);
			mpf_EGlpNumDivTo (q, r->adj[k2]->coef);
			if (mpf_EGlpNumIsEqual (q, mpf_oneLpNum, mpf_ILL_PRE_ZERO_TOL)) {
			    s[ci] = t;
			    s[j] = t;
			    got++;
			}
		    }
		}
		if (got) {
		    t++;
		} else {
		    s[ci] = mpf_ILL_MAXINT;
		    if (--nu == 0)
			goto DONE;
		}
	    }
	}

	if (n == 1) {
	    s[c0] = mpf_ILL_MAXINT;
	    if (--nu == 0)
		goto DONE;
	}
    }

DONE:

    {
	int dcount;
	int *dcnt;
	int *dlist;

	rval = mpf_gather_dup_lists (s, ncols, &dcount, &dcnt, &dlist);
	ILL_CLEANUP_IF (rval);
    }

    printf ("Time in mpf_duplicate_cols: %.2f (seconds)\n", ILLutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    ILL_IFFREE (s, int);
    mpf_EGlpNumFreeArray (f);
    mpf_EGlpNumClearVar (q);
    ILL_RETURN (rval, "mpf_duplicate_cols");
}

static int mpf_gather_dup_lists ( /* mpf_graph *G, */ int *s,	/* double *f, */

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
	if (s[i] < mpf_ILL_MAXINT && s[i] > smax)
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
	if (s[i] < mpf_ILL_MAXINT) {
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
	if (s[i] < mpf_ILL_MAXINT && s[i] > 0) {
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

    ILL_RETURN (rval, "mpf_gather_dup_lists");
}

static void mpf_set_fixed_variable (mpf_graph * G,
      int j,
      mpf_t val)
{
    int k;
    mpf_edge *e;

    G->cols[j].del = 1;
    for (k = 0; k < G->cols[j].deg; k++) {
	e = G->cols[j].adj[k];
	if (e->del == 0) {
	    /* G->rows[e->row].rhs -= (e->coef * val); */
	    mpf_EGlpNumSubInnProdTo (G->rows[e->row].rhs, e->coef, val);
	    e->del = 1;
	}
    }
}

static void mpf_get_implied_rhs_bounds (mpf_graph * G,
      int i,
      mpf_t * lb,
      mpf_t * ub)
{
    int k;
    mpf_t l, u;
    mpf_node *cols = G->cols;
    mpf_node *rows = G->rows;
    mpf_edge *e;
    mpf_EGlpNumInitVar (u);
    mpf_EGlpNumInitVar (l);

    mpf_EGlpNumZero (l);
    for (k = 0; k < rows[i].deg; k++) {
	e = rows[i].adj[k];
	if (e->del == 0) {
	    if (mpf_EGlpNumIsLess (e->coef, mpf_zeroLpNum)) {
		if (mpf_EGlpNumIsEqqual (cols[e->col].upper, mpf_ILL_MAXDOUBLE)) {
		    mpf_EGlpNumCopy (l, mpf_ILL_MINDOUBLE);
		    break;
		} else {
		    /* l += (e->coef * cols[e->col].upper); */
		    mpf_EGlpNumAddInnProdTo (l, e->coef, cols[e->col].upper);
		}
	    } else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, e->coef)) {
		if (mpf_EGlpNumIsEqqual (cols[e->col].lower, mpf_ILL_MINDOUBLE)) {
		    mpf_EGlpNumCopy (l, mpf_ILL_MINDOUBLE);
		    break;
		} else {
		    /* l += (e->coef * cols[e->col].lower); */
		    mpf_EGlpNumAddInnProdTo (l, e->coef, cols[e->col].lower);
		}
	    }
	}
    }

    mpf_EGlpNumZero (u);
    for (k = 0; k < rows[i].deg; k++) {
	e = rows[i].adj[k];
	if (e->del == 0) {
	    if (mpf_EGlpNumIsLess (e->coef, mpf_zeroLpNum)) {
		if (mpf_EGlpNumIsEqqual (cols[e->col].lower, mpf_ILL_MINDOUBLE)) {
		    mpf_EGlpNumCopy (u, mpf_ILL_MAXDOUBLE);
		} else {
		    /* u += (e->coef * cols[e->col].lower); */
		    mpf_EGlpNumAddInnProdTo (u, e->coef, cols[e->col].lower);
		}
	    } else if (mpf_EGlpNumIsLess (mpf_zeroLpNum, e->coef)) {
		if (mpf_EGlpNumIsEqqual (cols[e->col].upper, mpf_ILL_MAXDOUBLE)) {
		    mpf_EGlpNumCopy (u, mpf_ILL_MAXDOUBLE);
		} else {
		    /* u += (e->coef * cols[e->col].upper); */
		    mpf_EGlpNumAddInnProdTo (u, e->coef, cols[e->col].upper);
		}
	    }
	}
    }

    mpf_EGlpNumCopy (*lb, l);
    mpf_EGlpNumCopy (*ub, u);
    mpf_EGlpNumClearVar (u);
    mpf_EGlpNumClearVar (l);
}

static void mpf_get_implied_variable_bounds (mpf_graph * G,
      int j,
      mpf_edge * a_ij,
      mpf_t * lb,
      mpf_t * ub)
{
    int i = a_ij->row;
    mpf_t l, u;
    mpf_EGlpNumInitVar (u);
    mpf_EGlpNumInitVar (l);

    mpf_get_implied_rhs_bounds (G, i, &l, &u);
    mpf_EGlpNumCopy (*lb, mpf_ILL_MINDOUBLE);
    mpf_EGlpNumCopy (*ub, mpf_ILL_MAXDOUBLE);

    if (mpf_EGlpNumIsLess (mpf_ILL_PRE_FEAS_TOL, a_ij->coef)) {
	if (mpf_EGlpNumIsLess (u, mpf_ILL_MAXDOUBLE)) {
	    /**lb = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].upper;*/
	    mpf_EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, u, a_ij->coef);
	    mpf_EGlpNumAddTo (*lb, G->cols[j].upper);
	}
	if (mpf_EGlpNumIsLess (mpf_ILL_MINDOUBLE, l)) {
	    /**ub = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].lower;*/
	    mpf_EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, l, a_ij->coef);
	    mpf_EGlpNumAddTo (*ub, G->cols[j].lower);
	}
    } else if (mpf_EGlpNumIsLess (a_ij->coef, mpf_ILL_PRE_FEAS_TOL)) {
	if (mpf_EGlpNumIsLess (mpf_ILL_MINDOUBLE, l)) {
	    /**lb = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].upper;*/
	    mpf_EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, l, a_ij->coef);
	    mpf_EGlpNumAddTo (*lb, G->cols[j].upper);
	}
	if (mpf_EGlpNumIsLess (u, mpf_ILL_MAXDOUBLE)) {
	    /**ub = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].lower;*/
	    mpf_EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, u, a_ij->coef);
	    mpf_EGlpNumAddTo (*ub, G->cols[j].lower);
	}
    }
    mpf_EGlpNumClearVar (u);
    mpf_EGlpNumClearVar (l);
}

static int mpf_get_next_preop (mpf_ILLlp_predata * pre,
      mpf_ILLlp_preop ** op)
{
    int rval = 0;

    if (pre->opcount >= pre->opsize) {
	pre->opsize *= 1.3;
	pre->opsize += 1000;
	if (pre->opsize < pre->opcount + 1)
	    pre->opsize = pre->opcount + 1;
	pre->oplist = EGrealloc (pre->oplist, sizeof (mpf_ILLlp_preop) * pre->opsize);
	/*
			rval = ILLutil_reallocrus_scale ((void **) &pre->oplist, &pre->opsize, pre->opcount + 1, 1.3, sizeof (mpf_ILLlp_preop));
			ILL_CLEANUP_IF (rval);
	*/
    }
    *op = &pre->oplist[pre->opcount];
    mpf_ILLlp_preop_init (*op);

    /* CLEANUP: */

    ILL_RETURN (rval, "mpf_get_next_preop");
}

static int mpf_add_to_list (ILLptrworld * world,
      mpf_intptr ** list,
      int i)
{
    int rval = 0;
    mpf_intptr *ip;

    ip = intptralloc (world);
    if (!ip) {
	rval = 1;
	goto CLEANUP;
    }
    ip->this = i;
    ip->next = *list;
    *list = ip;

CLEANUP:

    ILL_RETURN (rval, "mpf_add_to_list");
}

static int mpf_build_graph (mpf_ILLlpdata * lp,
      mpf_graph * G)
{
    int rval = 0;
    int ncols = lp->ncols;
    int nrows = lp->nrows;
    int nzcount = lp->nzcount;
    int i, j, k, stop, count;
    mpf_edge *edgelist;
    mpf_node *rows, *cols;
    mpf_ILLmatrix *A = &lp->A;

    G->objsense = lp->objsense;

    ILL_SAFE_MALLOC (G->rows, nrows, mpf_node);
    if (!G->rows) {
	fprintf (stderr, "out of memory in mpf_build_graph\n");
	rval = 1;
	goto CLEANUP;
    }
    rows = G->rows;

    for (i = 0; i < nrows; i++) {
	rows[i].rowsense = lp->sense[i];
	rows[i].deg = 0;
    }

    ILL_SAFE_MALLOC (G->cols, ncols, mpf_node);
    ILL_SAFE_MALLOC (G->edgelist, nzcount, mpf_edge);
    for (i = nzcount; i--;)
	mpf_EGlpNumInitVar ((G->edgelist[i].coef));
    G->nzcount = nzcount;
    ILL_SAFE_MALLOC (G->adjspace, 2 * nzcount, mpf_edge *);

    if (!G->cols || !G->edgelist || !G->adjspace) {
	fprintf (stderr, "out of memory in mpf_build_graph\n");
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
	cols[j].coltype = mpf_ILL_PRE_COL_STRUC;
    }
    for (i = 0; i < nrows; i++) {
	cols[lp->rowmap[i]].coltype = mpf_ILL_PRE_COL_LOGICAL;
    }

    for (j = 0, count = 0; j < ncols; j++) {
	mpf_EGlpNumCopy (cols[j].obj, lp->obj[j]);
	mpf_EGlpNumCopy (cols[j].lower, lp->lower[j]);
	mpf_EGlpNumCopy (cols[j].upper, lp->upper[j]);
	stop = A->matbeg[j] + A->matcnt[j];
	for (k = A->matbeg[j]; k < stop; k++) {
	    i = A->matind[k];
	    rows[i].adj[rows[i].deg++] = &(edgelist[count]);
	    cols[j].adj[cols[j].deg++] = &(edgelist[count]);
	    edgelist[count].row = i;
	    edgelist[count].col = j;
	    mpf_EGlpNumCopy (edgelist[count].coef, A->matval[k]);
	    edgelist[count].mark = 0;
	    edgelist[count].del = 0;
	    edgelist[count].coltype = cols[j].coltype;
	    count++;
	}
    }
    if (count != nzcount) {
	fprintf (stderr, "counts are off in mpf_build_graph\n");
	rval = 1;
	goto CLEANUP;
    }
    G->ecount = count;
    G->nrows = nrows;
    G->ncols = ncols;

    for (i = 0; i < G->nrows; i++) {
	G->rows[i].del = 0;
	mpf_EGlpNumCopy (G->rows[i].rhs, lp->rhs[i]);
    }
    for (j = 0; j < G->ncols; j++) {
	G->cols[j].del = 0;
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_build_graph");
}

static void mpf_dump_graph (mpf_graph * G)
{
    int i, j, k;

    printf ("ecount = %d, nrows = %d, ncols = %d\n",
	G->ecount, G->nrows, G->ncols);
    fflush (stdout);

    for (i = 0; i < G->nrows; i++) {
	printf ("Row %d:", i);
	for (k = 0; k < G->rows[i].deg; k++) {
	    printf (" %d", G->rows[i].adj[k]->col);
	    if (G->rows[i].adj[k]->coltype == mpf_ILL_PRE_COL_LOGICAL)
		printf ("S");
	    printf ("(%g)", mpf_EGlpNumToLf (G->rows[i].adj[k]->coef));
	}
	printf ("  rhs: %g", mpf_EGlpNumToLf (G->rows[i].rhs));
	if (G->rows[i].del) {
	    printf (" (deleted)\n");
	} else {
	    printf ("\n");
	}
    }

    for (j = 0; j < G->ncols; j++) {
	if (G->cols[j].coltype == mpf_ILL_PRE_COL_LOGICAL) {
	    printf ("Slk %d:", j);
	} else {
	    printf ("Col %d:", j);
	}
	for (k = 0; k < G->cols[j].deg; k++) {
	    printf (" %d", G->cols[j].adj[k]->row);
	}
	printf ("  obj: %g  bnd: (%g, %g)", mpf_EGlpNumToLf (G->cols[j].obj),
	    mpf_EGlpNumToLf (G->cols[j].lower), mpf_EGlpNumToLf (G->cols[j].upper));
	if (G->cols[j].del) {
	    printf (" (deleted)\n");
	} else {
	    printf ("\n");
	}
    }
}

static void mpf_init_graph (mpf_graph * G)
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

static void mpf_free_graph (mpf_graph * G)
{
    register int i;
    if (G) {
	int total, onlist;
	for (i = G->nzcount; i--;)
	    mpf_EGlpNumClearVar ((G->edgelist[i].coef));
	ILL_IFFREE (G->edgelist, mpf_edge);
	ILL_IFFREE (G->rows, mpf_node);
	ILL_IFFREE (G->cols, mpf_node);
	ILL_IFFREE (G->adjspace, mpf_edge *);
	if (intptr_check_leaks (&G->intptrworld, &total, &onlist)) {
	    fprintf (stderr, "WARNING: %d outstanding intptrs\n", total - onlist);
	}
	ILLptrworld_delete (&G->intptrworld);
	mpf_init_graph (G);
    }
}

int mpf_ILLlp_sinfo_print (mpf_ILLlp_sinfo * s)
{
    int rval = 0;
    int i;
    mpf_ILLlpdata lp;
    char *sense = 0;

    mpf_ILLlpdata_init (&lp);

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
	fprintf (stderr, "out of memory in mpf_ILLlp_sinfo_print\n");
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
    ILL_RETURN (rval, "mpf_ILLlp_sinfo_print");
}

void mpf_ILLlp_sinfo_init (mpf_ILLlp_sinfo * sinfo)
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
	sinfo->objsense = mpf_ILL_MIN;
	mpf_ILLmatrix_init (&sinfo->A);
    }
}

void mpf_ILLlp_sinfo_free (mpf_ILLlp_sinfo * sinfo)
{
    if (sinfo) {
	mpf_EGlpNumFreeArray (sinfo->obj);
	mpf_EGlpNumFreeArray (sinfo->lower);
	mpf_EGlpNumFreeArray (sinfo->upper);
	mpf_EGlpNumFreeArray (sinfo->rhs);
	mpf_ILLmatrix_free (&sinfo->A);
	if (sinfo->colnames) {
	    int i;
	    for (i = 0; i < sinfo->ncols; i++) {
		ILL_IFFREE (sinfo->colnames[i], char);
	    }
	    ILL_IFFREE (sinfo->colnames, char *);
	}
	mpf_ILLlp_sinfo_init (sinfo);
    }
}

void mpf_ILLlp_predata_init (mpf_ILLlp_predata * pre)
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

void mpf_ILLlp_predata_free (mpf_ILLlp_predata * pre)
{
    if (pre) {
	int i;

	for (i = 0; i < pre->opcount; i++) {
	    mpf_ILLlp_preop_free (&pre->oplist[i]);
	}
	ILL_IFFREE (pre->oplist, mpf_ILLlp_preop);
	ILL_IFFREE (pre->colmap, int);
	ILL_IFFREE (pre->rowmap, int);
	ILL_IFFREE (pre->colscale, mpf_t);
	ILL_IFFREE (pre->rowscale, mpf_t);
	ILL_IFFREE (pre->colfixval, mpf_t);
	ILL_IFFREE (pre->rowfixval, mpf_t);
	mpf_ILLlp_predata_init (pre);
    }
}

void mpf_ILLlp_preop_init (mpf_ILLlp_preop * op)
{
    if (op) {
	op->ptype = 0;
	op->rowindex = -1;
	op->colindex = -1;
	mpf_ILLlp_preline_init (&op->line);
    }
}

void mpf_ILLlp_preop_free (mpf_ILLlp_preop * op)
{
    if (op) {
	mpf_ILLlp_preline_free (&op->line);
	mpf_ILLlp_preop_init (op);
    }
}

void mpf_ILLlp_preline_init (mpf_ILLlp_preline * line)
{
    if (line) {
	mpf_EGlpNumInitVar (line->rhs);
	mpf_EGlpNumInitVar (line->obj);
	mpf_EGlpNumInitVar (line->upper);
	mpf_EGlpNumInitVar (line->lower);
	mpf_EGlpNumZero (line->rhs);
	mpf_EGlpNumZero (line->obj);
	mpf_EGlpNumZero (line->upper);
	mpf_EGlpNumZero (line->lower);
	line->count = 0;
	line->ind = 0;
	line->val = 0;
    }
}

void mpf_ILLlp_preline_free (mpf_ILLlp_preline * line)
{
    if (line) {
	mpf_EGlpNumClearVar (line->rhs);
	mpf_EGlpNumClearVar (line->obj);
	mpf_EGlpNumClearVar (line->upper);
	mpf_EGlpNumClearVar (line->lower);
	ILL_IFFREE (line->ind, int);
	mpf_EGlpNumFreeArray (line->val);
	/* mpf_ILLlp_preline_init (line); */
    }
}
