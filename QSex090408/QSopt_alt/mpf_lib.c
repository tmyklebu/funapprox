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

/* RCS_INFO = "$RCSfile: mpf_lib.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* Interface Routines to Core LP Solver                       */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int mpf_ILLlib_optimize (mpf_lpinfo *lp, mpf_ILLlp_basis *B,
   mpf_price_info *pinf,    */
/* int algo, int *status, int simplex_display)                   */
/* int mpf_ILLlib_cache_solution (mpf_lpinfo *lp, mpf_ILLlp_cache *C)                */
/* int mpf_ILLlib_solution (mpf_lpinfo *lp, mpf_ILLlp_cache *C, double *val,         */
/* double *x, double *pi, double *slack, double *rc)             */
/* int mpf_ILLlib_get_x (mpf_lpinfo *lp, mpf_ILLlp_cache *C, double *x)              */
/* int mpf_ILLlib_get_slack (mpf_lpinfo *lp, mpf_ILLlp_cache *C, double
   slack)      */
/* int mpf_ILLlib_objval (mpf_lpinfo *lp, mpf_ILLlp_cache *C, double *val)           */
/* int mpf_ILLlib_newrow (mpf_lpinfo *lp, mpf_ILLlp_basis *B, double rhs,            */
/* char sense, double range, const char *name)                   */
/* -range can specify a rangeval for the row (if sense is not 'R',   */
/* then range is ignored); it should be 0 if no range is needed;    */
/* if sense is 'R' but no rangeval array exists for the LP, the     */
/* array will be allocated and initialized.                         */
/* int mpf_ILLlib_newrows (mpf_lpinfo *lp, mpf_ILLlp_basis *B, int num,
   double *rhs, */
/* char *sense, double *range, const char **names)               */
/* -range is an array specifying the rangevals for the rows; range   */
/* should be NULL if no rangevals are needed.                       */
/* int mpf_ILLlib_addrow (mpf_lpinfo *lp, mpf_ILLlp_basis *B, int cnt, int
   ind,     */
/* double *val, double rhs, char sense, double range,            */
/* const char *name)                                             */
/* int mpf_ILLlib_addrows (mpf_lpinfo *lp, mpf_ILLlp_basis *B, int num,              */
/* int *rmatcnt, int *rmatbeg, int *rmatind, double *rmatval,    */
/* double *rhs, char *sense, double *range, const char **names,  */
/* int *factorok)                                                */
/* int mpf_ILLlib_delrows (mpf_lpinfo *lp, mpf_ILLlp_basis *B,                       */
/* int num, int *dellist, int *basis_ok)                         */
/* int mpf_ILLlib_newcol (mpf_lpinfo *lp, mpf_ILLlp_basis *B,                        */
/* double obj, double lower, double upper, const char *name,     */
/* int factorok)                                                 */
/* int mpf_ILLlib_newcols (mpf_lpinfo *lp, mpf_ILLlp_basis *B,                       */
/* int num, double *obj, double *lower, double *upper,           */
/* const char **names, int factorok)                             */
/* int mpf_ILLlib_addcol (mpf_lpinfo *lp, mpf_ILLlp_basis *B,                        */
/* int cnt, int *ind, double *val, double obj, double lower,     */
/* double upper, const char *name, int factorok)                 */
/* int mpf_ILLlib_addcols (mpf_lpinfo *lp, mpf_ILLlp_basis *B,                       */
/* int num, int *cmatcnt, int *cmatbeg, int *cmatind,            */
/* double *cmatval, double *obj, double *lower, double *upper,   */
/* const char **names, int factorok)                             */
/* int mpf_ILLlib_delcols (mpf_lpinfo *lp, mpf_ILLlp_basis *B, int num, int
   dellist */
/* int *basis_ok)                                                */
/* int mpf_ILLlib_chgcoef (mpf_lpinfo *lp, int rowindex, int colindex,           */
/* double coef)                                                  */
/* int mpf_ILLlib_chgsense (mpf_lpinfo *lp, int num, int *rowlist, char
   sense)  */
/* int mpf_ILLlib_getrows (mpf_lpinfo *lp, int num, int *rowlist, int
   *rowcnt,  */
/* int **rowbeg, int **rowind, double **rowval, double **rhs,    */
/* char **sense, char ***names)                                  */
/* int mpf_ILLlib_getcols (mpf_lpinfo *lp, int num, int *collist, int
   *colcnt,  */
/* int **colbeg, int **colind, double **colval, double **obj,    */
/* double **lower, double **upper, char ***names)                */
/* int mpf_ILLlib_getobj (mpf_lpinfo *lp, double *obj)                           */
/* int mpf_ILLlib_chgobj (mpf_lpinfo *lp, int indx, double coef)                 */
/* int mpf_ILLlib_getrhs (mpf_lpinfo *lp, double *rhs)                           */
/* int mpf_ILLlib_chgrhs (mpf_lpinfo *lp, int indx, double coef)                 */
/* int mpf_ILLlib_getintflags (mpf_lpinfo *lp, int *intflags)                    */
/* int mpf_ILLlib_rownames (mpf_lpinfo *lp, char **rownames)                     */
/* int mpf_ILLlib_colnames (mpf_lpinfo *lp, char **colnames)                     */
/* int mpf_ILLlib_colindex (mpf_lpinfo *lp, char *name, int *colindex)           */
/* int mpf_ILLlib_rowindex (mpf_lpinfo *lp, char *name, int *rowindex)           */
/* int mpf_ILLlib_chgbnd  (mpf_lpinfo *lp, int indx, char lu, double bnd)        */
/* int mpf_ILLlib_chgbnds (mpf_lpinfo *lp, int cnt, int *indx, char *lu,         */
/* double *bnd)                                                  */
/* int mpf_ILLlib_getbnd (mpf_lpinfo *lp, int indx, char lu, double *bnd)        */
/* int mpf_ILLlib_getbnds (mpf_lpinfo *lp, double *lower, double *upper)         */
/* int mpf_ILLlib_strongbranch (mpf_lpinfo *lp, mpf_price_info *pinf,                */
/* int *candidatelist, int ncand, double *xlist, double *downpen,      */
/* double *uppen, int iterations, double objbound)                     */
/* int mpf_ILLlib_getbasis (mpf_lpinfo *lp, char *cstat, char *rstat)            */
/* int mpf_ILLlib_loadbasis (mpf_ILLlp_basis *B, int nstruct, int nrows,         */
/* char *cstat, char *rstat)                                           */
/* int mpf_ILLlib_readbasis (mpf_lpinfo *lp, mpf_ILLlp_basis *B, char
   mpf_fname)        */
/* int mpf_ILLlib_writebasis (mpf_lpinfo *lp, const char *mpf_fname)                 */
/* int mpf_ILLlib_getrownorms (mpf_lpinfo *lp, mpf_price_info *pinf,                 */
/* double *rownorms)                                             */
/* int mpf_ILLlib_loadrownorms (mpf_lpinfo *lp, mpf_price_info *pinf,                */
/* double *rownorms)                                             */
/* int mpf_ILLlib_recompute_rownorms (mpf_lpinfo *lp, mpf_price_info *pinf)          */
/* int mpf_ILLlib_print_x (FILE *fd, mpf_lpinfo *lp, mpf_ILLlp_cache *C,
   double *x,  */
/* int nonZerosOnly)                                             */
/* int mpf_ILLlib_print_x (mpf_lpinfo *lp, mpf_ILLlp_cache *C)                       */
/* int mpf_ILLlib_iter (mpf_lpinfo *lp)                                          */
/* */
/* NOTES                                                                   */
/* */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdata.h"
#include "mpf_lpdefs.h"
#include "mpf_simplex.h"
#include "mpf_price.h"
#include "mpf_basis.h"
#include "mpf_lib.h"
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#include "mpf_lp.h"
#include "mpf_mps.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static void mpf_check_pinf (mpf_price_info * pinf,
      int *it_exists);

static int mpf_matrix_addrow (mpf_ILLmatrix * A,
      int rowcnt,
      int *rowind,
      mpf_t * rowval),
  mpf_matrix_addrow_end (mpf_ILLmatrix * A,
      int row,
      int rowcnt,
      int *rowind,
      mpf_t * rowval),
  mpf_matrix_addcoef (mpf_lpinfo * lp,
      mpf_ILLmatrix * A,
      int row,
      int col,
      mpf_t val),
  mpf_matrix_addcol (mpf_ILLmatrix * A,
      int colcnt,
      int *colind,
      mpf_t * colval),
  mpf_delcols_work (mpf_lpinfo * lp,
      char *colmark),
  mpf_reset_colindex (mpf_lpinfo * lp),
  mpf_reset_rowindex (mpf_lpinfo * lp);

int mpf_ILLlib_optimize (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      mpf_price_info * pinf,
      int algo,
      int *status,
      int simplex_display)
{
    int rval = 0;
    int sol_status;

    if (status)
	*status = QS_LP_UNSOLVED;

    /* mpf_ILLprice_free_pricing_info (pinf); */
    /* Should be removed later */

    rval = mpf_ILLsimplex (lp, algo, B, pinf, &sol_status, simplex_display);
    ILL_CLEANUP_IF (rval);

    if (status)
	*status = sol_status;

CLEANUP:

    if (rval == E_SIMPLEX_ERROR) {
	FILE *eout = 0;
	int tval;

	printf ("write bad lp to error.lp\n");
	fflush (stdout);
	eout = fopen ("error.lp", "w");
	if (!eout) {
	    fprintf (stderr, "could not open file to write bad lp\n");
	} else {
	    tval = mpf_ILLwrite_lp (lp->O, NULL);
	    if (tval) {
		fprintf (stderr, "error while writing bad lp\n");
	    }
	    fclose (eout);
	}

	printf ("write bad basis to error.bas\n");
	fflush (stdout);
	tval = mpf_ILLlib_writebasis (lp, 0, "error.bas");
	if (tval) {
	    fprintf (stderr, "error while writing bad basis\n");
	}
    }
    ILL_RETURN (rval, "mpf_ILLlib_optimize");
}

int mpf_ILLlib_cache_solution (mpf_lpinfo * lp,
      mpf_ILLlp_cache * C)
{
    int rval = 0;

    if (C) {
	if (C->nstruct != lp->O->nstruct || C->nrows != lp->O->nrows) {
	    fprintf (stderr, "lp_cache does not match size of lp\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
	rval = mpf_ILLlib_solution (lp, 0, &(C->val), C->x, C->pi, C->slack, C->rc);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_cache_solution");
}

int mpf_ILLlib_solution (mpf_lpinfo * lp,
      mpf_ILLlp_cache * C,
      mpf_t * val,
      mpf_t * x,
      mpf_t * pi,
      mpf_t * slack,
      mpf_t * rc)
{
    int i, rval = 0;
    mpf_t *tempx = 0;
    mpf_t *temprc = 0;
    int ncols = lp->O->ncols;
    int nrows = lp->O->nrows;
    int nstruct = lp->O->nstruct;
    mpf_ILLlpdata *qslp = lp->O;

    if (C) {
	if (C->nrows != nrows || C->nstruct != nstruct) {
	    fprintf (stderr, "cache mismatch in mpf_ILLlib_solution\n");
	    rval = 0;
	    ILL_CLEANUP;
	}
	if (val) {
	    mpf_EGlpNumCopy (*val, C->val);
	}
	if (x) {
	    for (i = 0; i < nstruct; i++) {
		mpf_EGlpNumCopy (x[i], C->x[i]);
	    }
	}
	if (pi) {
	    for (i = 0; i < nrows; i++) {
		mpf_EGlpNumCopy (pi[i], C->pi[i]);
	    }
	}
	if (slack) {
	    for (i = 0; i < nrows; i++) {
		mpf_EGlpNumCopy (slack[i], C->slack[i]);
	    }
	}
	if (rc) {
	    for (i = 0; i < nstruct; i++) {
		mpf_EGlpNumCopy (rc[i], C->rc[i]);
	    }
	}
    } else {
	if (x || slack)
	    tempx = mpf_EGlpNumAllocArray (ncols);

	if (rc)
	    temprc = mpf_EGlpNumAllocArray (ncols);

	rval = mpf_ILLsimplex_solution (lp, tempx, pi, temprc, val);
	ILL_CLEANUP_IF (rval);

	if (x) {
	    for (i = 0; i < nstruct; i++) {
		mpf_EGlpNumCopy (x[i], tempx[qslp->structmap[i]]);
	    }
	}
	if (slack) {
	    for (i = 0; i < nrows; i++) {
		mpf_EGlpNumCopy (slack[i], tempx[qslp->rowmap[i]]);
	    }
	}
	if (rc) {
	    for (i = 0; i < nstruct; i++) {
		mpf_EGlpNumCopy (rc[i], temprc[qslp->structmap[i]]);
	    }
	}
	if (lp->O->objsense == mpf_ILL_MAX) {	/* Reverse signs for max prob */
	    if (val) {
		mpf_EGlpNumSign (*val);
	    }
	    if (pi) {
		for (i = 0; i < nrows; i++) {
		    mpf_EGlpNumSign (pi[i]);
		}
	    }
	    if (rc) {
		for (i = 0; i < nstruct; i++) {
		    mpf_EGlpNumSign (rc[i]);
		}
	    }
	}
    }

CLEANUP:

    mpf_EGlpNumFreeArray (tempx);
    mpf_EGlpNumFreeArray (temprc);
    ILL_RETURN (rval, "mpf_ILLlib_solution");
}

int mpf_ILLlib_get_x (mpf_lpinfo * lp,
      mpf_ILLlp_cache * C,
      mpf_t * x)
{
    int rval = 0;

    rval = mpf_ILLlib_solution (lp, C, 0, x, 0, 0, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_get_x");
}

int mpf_ILLlib_get_slack (mpf_lpinfo * lp,
      mpf_ILLlp_cache * C,
      mpf_t * slack)
{
    int rval = 0;

    rval = mpf_ILLlib_solution (lp, C, 0, 0, 0, slack, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_get_slack");
}


int mpf_ILLlib_objval (mpf_lpinfo * lp,
      mpf_ILLlp_cache * C,
      mpf_t * val)
{
    int rval = 0;

    if (lp->basisstat.optimal) {
	rval = mpf_ILLlib_solution (lp, C, val, 0, 0, 0, 0);
	ILL_CLEANUP_IF (rval);
    } else {
	mpf_EGlpNumCopy (*val, lp->dobjval);	/* Ask Sanjeeb */
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_objval");
}

int mpf_ILLlib_tableau (mpf_lpinfo * lp,
      int row,
      mpf_t * binv,
      mpf_t * tabrow)
{
    int rval = 0;
    int i;
    int ncols = lp->O->ncols;
    int nrows = lp->O->nrows;
    int nstruct = lp->O->nstruct;
    mpf_t *brow = 0;
    mpf_t *trow = 0;
    mpf_ILLlpdata *qslp = lp->O;

    if (row < 0 || row >= qslp->nrows) {
	fprintf (stderr, "mpf_ILLlib_tableau called with bad row: %d\n", row);
	rval = 1;
	ILL_CLEANUP;
    }
    brow = mpf_EGlpNumAllocArray (nrows);

    if (tabrow)
	trow = mpf_EGlpNumAllocArray (ncols);

    rval = mpf_ILLbasis_tableau_row (lp, row, brow, trow, 0, 0);
    ILL_CLEANUP_IF (rval);

    if (binv) {
	for (i = 0; i < nrows; i++) {
	    mpf_EGlpNumCopy (binv[i], brow[i]);
	}
    }
    if (tabrow) {
	for (i = 0; i < nstruct; i++) {
	    mpf_EGlpNumCopy (tabrow[i], trow[qslp->structmap[i]]);
	}
	for (i = 0; i < nrows; i++) {
	    mpf_EGlpNumCopy (tabrow[nstruct + i], trow[qslp->rowmap[i]]);
	}
    }
CLEANUP:

    mpf_EGlpNumFreeArray (brow);
    mpf_EGlpNumFreeArray (trow);
    ILL_RETURN (rval, "mpf_ILLlib_tableau");
}

int mpf_ILLlib_basis_order (mpf_lpinfo * lp,
      int *header)
{
    int rval = 0;
    int i, j;
    int ncols = lp->O->ncols;
    int nrows = lp->O->nrows;
    int nstruct = lp->O->nstruct;
    mpf_ILLlpdata *qslp = lp->O;
    int *invmap = 0;

    ILL_SAFE_MALLOC (invmap, ncols, int);
    for (j = 0; j < nstruct; j++) {
	invmap[qslp->structmap[j]] = j;
    }
    for (i = 0; i < nrows; i++) {
	invmap[qslp->rowmap[i]] = nstruct + i;
    }

    for (i = 0; i < nrows; i++) {
	header[i] = invmap[lp->baz[i]];
    }

CLEANUP:

    ILL_IFFREE (invmap, int);
    ILL_RETURN (rval, "mpf_ILLlib_basis_order");
}

int mpf_ILLlib_chgbnd (mpf_lpinfo * lp,
      int indx,
      int lu,
      mpf_t bnd)
{
    int rval = 0;
    int col;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_chgbnd called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (indx < 0 || indx > lp->O->nstruct) {
	fprintf (stderr, "mpf_ILLlib_chgbnd called with bad indx: %d\n", indx);
	rval = 1;
	ILL_CLEANUP;
    }
    if (lp->O->sinfo) {		/* Presolve LP is no longer valid, free the
				   data */
	mpf_ILLlp_sinfo_free (lp->O->sinfo);
	ILL_IFFREE (lp->O->sinfo, mpf_ILLlp_sinfo);
    }
    col = lp->O->structmap[indx];

    switch (lu) {
    case 'L':
	mpf_EGlpNumCopy (lp->O->lower[col], bnd);
	break;
    case 'U':
	mpf_EGlpNumCopy (lp->O->upper[col], bnd);
	break;
    case 'B':
	mpf_EGlpNumCopy (lp->O->lower[col], bnd);
	mpf_EGlpNumCopy (lp->O->upper[col], bnd);
	break;
    default:
	fprintf (stderr, "mpf_ILLlib_chgbnd called with lu: %c\n", lu);
	rval = 1;
	ILL_CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_chgbnd");
}

int mpf_ILLlib_chgbnds (mpf_lpinfo * lp,
      int cnt,
      int *indx,
      char *lu,
      mpf_t * bnd)
{
    int rval = 0;
    int i;

    for (i = 0; i < cnt; i++) {
	rval = mpf_ILLlib_chgbnd (lp, indx[i], lu[i], bnd[i]);
	if (rval)
	    ILL_CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_chgbnds");
}

int mpf_ILLlib_getbnd (mpf_lpinfo * lp,
      int indx,
      int lu,
      mpf_t * bnd)
{
    int rval = 0;
    int col;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getbnd called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (indx < 0 || indx > lp->O->nstruct) {
	fprintf (stderr, "mpf_ILLlib_getbnd called with bad indx: %d\n", indx);
	rval = 1;
	ILL_CLEANUP;
    }
    col = lp->O->structmap[indx];

    switch (lu) {
    case 'L':
	mpf_EGlpNumCopy (*bnd, lp->O->lower[col]);
	break;
    case 'U':
	mpf_EGlpNumCopy (*bnd, lp->O->upper[col]);
	break;
    default:
	fprintf (stderr, "mpf_ILLlib_getbnd called with lu: %c\n", lu);
	rval = 1;
	ILL_CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_getbnd");
}

int mpf_ILLlib_getbnds (mpf_lpinfo * lp,
      mpf_t * lower,
      mpf_t * upper)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;
    int nstruct;
    int j, col;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getbnd called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nstruct = qslp->nstruct;

    for (j = 0; j < nstruct; j++) {
	col = qslp->structmap[j];
	if (lower)
	    mpf_EGlpNumCopy (lower[j], qslp->lower[col]);
	if (upper)
	    mpf_EGlpNumCopy (upper[j], qslp->upper[col]);
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_getbnds");
}

int mpf_ILLlib_strongbranch (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      int *candidatelist,
      int ncand,
      mpf_t * xlist,
      mpf_t * downpen,
      mpf_t * uppen,
      int iterations,
      mpf_t objbound)
{
    int rval = 0;
    int i, k, status, have_norms;
    int olditer = lp->maxiter;
    int nstruct = lp->O->nstruct;
    int nrows = lp->O->nrows;
    mpf_t *myx = 0;
    mpf_t xi, t, oldbnd;
    mpf_price_info lpinf;
    mpf_ILLlp_basis B, origB;

    mpf_EGlpNumInitVar (lpinf.htrigger);
    mpf_EGlpNumInitVar (xi);
    mpf_EGlpNumInitVar (t);
    mpf_EGlpNumInitVar (oldbnd);
    mpf_EGlpNumZero (oldbnd);
    mpf_ILLlp_basis_init (&B);
    mpf_ILLlp_basis_init (&origB);
    mpf_ILLprice_init_pricing_info (&lpinf);
    lpinf.dI_price = QS_PRICE_DSTEEP;
    lpinf.dII_price = QS_PRICE_DSTEEP;

    if (xlist == 0) {
	myx = mpf_EGlpNumAllocArray (nstruct);
	rval = mpf_ILLlib_get_x (lp, 0, myx);
	ILL_CLEANUP_IF (rval);
    }
    rval = mpf_ILLlp_basis_alloc (&origB, nstruct, nrows);
    ILL_CLEANUP_IF (rval);

    rval = mpf_ILLlib_getbasis (lp, origB.cstat, origB.rstat);
    ILL_CLEANUP_IF (rval);

    mpf_check_pinf (pinf, &have_norms);
    if (have_norms == 0) {
	origB.rownorms = mpf_EGlpNumAllocArray (nrows);
	rval = mpf_ILLlib_getrownorms (lp, pinf, origB.rownorms);
	ILL_CLEANUP_IF (rval);
    } else {
	lp->basisid = -1;
	rval = mpf_ILLlib_optimize (lp, 0, &lpinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);
    }

    rval = mpf_ILLlp_basis_alloc (&B, nstruct, nrows);	/* Note: B and orgiB may */
    /* differ.               */
    ILL_CLEANUP_IF (rval);

    rval = mpf_ILLlib_getbasis (lp, B.cstat, B.rstat);
    ILL_CLEANUP_IF (rval);
    B.rownorms = mpf_EGlpNumAllocArray (nrows);

    if (have_norms == 0) {
	rval = mpf_ILLlib_getrownorms (lp, pinf, B.rownorms);
	ILL_CLEANUP_IF (rval);
    } else {
	rval = mpf_ILLlib_getrownorms (lp, &lpinf, B.rownorms);
	ILL_CLEANUP_IF (rval);
    }

    lp->maxiter = iterations;

    for (i = 0; i < ncand; i++) {
	k = candidatelist[i];
	rval = mpf_ILLlib_getbnd (lp, k, 'U', &oldbnd);
	ILL_CLEANUP_IF (rval);
	if (xlist)
	    mpf_EGlpNumCopy (xi, xlist[i]);
	else
	    mpf_EGlpNumCopy (xi, myx[k]);
	mpf_EGlpNumFloor (t, xi);
	if (mpf_EGlpNumIsLessDbl (t, 0.1) && mpf_EGlpNumIsGreaDbl (t, -0.1))
	    mpf_EGlpNumZero (t);

	rval = mpf_ILLlib_chgbnd (lp, k, 'U', t);
	ILL_CLEANUP_IF (rval);

	rval = mpf_ILLlib_optimize (lp, &B, &lpinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);

	mpf_EGlpNumCopy (downpen[i], lp->dobjval);
	rval = mpf_ILLlib_chgbnd (lp, k, 'U', oldbnd);
	ILL_CLEANUP_IF (rval);

	rval = mpf_ILLlib_getbnd (lp, k, 'L', &oldbnd);
	ILL_CLEANUP_IF (rval);
	mpf_EGlpNumCeil (t, xi);
	if (mpf_EGlpNumIsLessDbl (t, 1.1) && mpf_EGlpNumIsGreaDbl (t, 0.9))
	    mpf_EGlpNumOne (t);
	rval = mpf_ILLlib_chgbnd (lp, k, 'L', t);
	ILL_CLEANUP_IF (rval);

	rval = mpf_ILLlib_optimize (lp, &B, &lpinf, DUAL_SIMPLEX, &status, 0);
	ILL_CLEANUP_IF (rval);

	mpf_EGlpNumCopy (uppen[i], lp->dobjval);
	rval = mpf_ILLlib_chgbnd (lp, k, 'L', oldbnd);
	ILL_CLEANUP_IF (rval);
    }

    if (lp->O->objsense == mpf_ILL_MAX) {

    } else {
	for (i = 0; i < ncand; i++) {
	    if (mpf_EGlpNumIsLess (objbound, downpen[i]))
		mpf_EGlpNumCopy (downpen[i], objbound);
	    if (mpf_EGlpNumIsLess (objbound, uppen[i]))
		mpf_EGlpNumCopy (uppen[i], objbound);
	}
    }

    /* Restore the old optimal solution */

    lp->maxiter = olditer;
    rval = mpf_ILLlib_optimize (lp, &origB, pinf, DUAL_SIMPLEX, &status, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    mpf_EGlpNumClearVar (xi);
    mpf_EGlpNumClearVar (t);
    mpf_EGlpNumClearVar (oldbnd);
    lp->maxiter = olditer;
    mpf_ILLprice_free_pricing_info (&lpinf);
    mpf_ILLlp_basis_free (&B);
    mpf_ILLlp_basis_free (&origB);
    if (xlist == 0)
	mpf_EGlpNumFreeArray (myx);
    mpf_EGlpNumClearVar (lpinf.htrigger);
    ILL_RETURN (rval, "mpf_ILLlib_strongbranch");
}

#define mpf_EXTRA_ROWS (100)
#define mpf_EXTRA_COLS (100)
#define mpf_EXTRA_MAT  (1000)

int mpf_ILLlib_newrow (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      mpf_t rhs,
      int sense,
      mpf_t range,
      const char *name)
{
    int rval = 0;

    rval = mpf_ILLlib_addrow (lp, B, 0, 0, 0, rhs, sense, range, name);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_newrow");
}

int mpf_ILLlib_newrows (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int num,
      mpf_t * rhs,
      char *sense,
      mpf_t * range,
      const char **names)
{
    int rval = 0;
    int *rmatcnt = 0;
    int *rmatbeg = 0;
    int i;

    if (!num)
	ILL_CLEANUP;

    ILL_SAFE_MALLOC (rmatcnt, num, int);

    ILL_SAFE_MALLOC (rmatbeg, num, int);

    for (i = 0; i < num; i++) {
	rmatcnt[i] = 0;
	rmatbeg[i] = 0;
    }

    rval = mpf_ILLlib_addrows (lp, B, num, rmatcnt, rmatbeg, 0, 0, rhs, sense,
	range, names, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_IFFREE (rmatcnt, int);
    ILL_IFFREE (rmatbeg, int);

    ILL_RETURN (rval, "mpf_ILLlib_newrows");
}

int mpf_ILLlib_addrows (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int num,
      int *rmatcnt,
      int *rmatbeg,
      int *rmatind,
      mpf_t * rmatval,
      mpf_t * rhs,
      char *sense,
      mpf_t * range,
      const char **names,
      int *factorok)
{
    int rval = 0;
    int i, j, total, bsing;
    int *imap = 0;
    int *bbeg = 0;
    int *bcnt = 0;
    int *bindi = 0;
    int *rindi = 0;
    int *jstat = 0;
    mpf_t *bval = 0;
    mpf_t rng;
    int badfactor = 0;
    mpf_EGlpNumInitVar (rng);

    if (B == 0 || B->rownorms == 0) {
	if (factorok)
	    *factorok = 0;
    }
    if (B)
	mpf_EGlpNumFreeArray (B->colnorms);

    if (B && B->rownorms && factorok && *factorok == 1) {
	int *structmap = lp->O->structmap;

	lp->matbeg = lp->O->A.matbeg;
	lp->matcnt = lp->O->A.matcnt;
	lp->matind = lp->O->A.matind;
	lp->matval = lp->O->A.matval;

	lp->nrows = lp->O->nrows;
	lp->ncols = lp->O->ncols;
	if (B->rownorms_size < lp->O->nrows + num)
	    mpf_EGlpNumReallocArray (&(B->rownorms), lp->O->nrows + num);

	ILL_SAFE_MALLOC (bcnt, num, int);
	ILL_SAFE_MALLOC (bbeg, num, int);
	ILL_SAFE_MALLOC (imap, lp->O->nstruct, int);

	ILL_SAFE_MALLOC (jstat, lp->ncols, int);

	for (i = 0; i < lp->ncols; i++) {
	    jstat[i] = -1;
	}
	for (i = 0; i < lp->O->nstruct; i++) {
	    jstat[structmap[i]] = i;
	}

	for (i = 0; i < lp->O->nstruct; i++) {
	    imap[i] = -1;
	}
	for (i = 0; i < lp->nrows; i++) {
	    if (jstat[lp->baz[i]] != -1) {
		imap[jstat[lp->baz[i]]] = i;
	    }
	}

	for (i = 0, total = 0; i < num; i++) {
	    bcnt[i] = 0;
	    bbeg[i] = total;
	    for (j = 0; j < rmatcnt[i]; j++) {
		if (imap[rmatind[rmatbeg[i] + j]] != -1) {
		    bcnt[i]++;
		    total++;
		}
	    }
	}
	if (total) {
	    ILL_SAFE_MALLOC (bindi, total, int);
	    bval = mpf_EGlpNumAllocArray (total);
	}
	for (i = 0, total = 0; i < num; i++) {
	    for (j = 0; j < rmatcnt[i]; j++) {
		if (imap[rmatind[rmatbeg[i] + j]] != -1) {
		    mpf_EGlpNumCopy (bval[total], rmatval[rmatbeg[i] + j]);
		    bindi[total] = imap[rmatind[rmatbeg[i] + j]];
		    total++;
		}
	    }
	}

	rval = mpf_ILLprice_get_new_rownorms (lp, num, B->rownorms + lp->O->nrows,
	    bcnt, bbeg, bindi, bval);
	ILL_CLEANUP_IF (rval);

	ILL_IFFREE (bcnt, int);
	ILL_IFFREE (bbeg, int);
	ILL_IFFREE (bindi, int);
	mpf_EGlpNumFreeArray (bval);
	ILL_IFFREE (imap, int);

	badfactor = 1;
    }
    for (i = 0; i < num; i++) {
	if (range)
	    mpf_EGlpNumCopy (rng, range[i]);
	else
	    mpf_EGlpNumZero (rng);
	if (names) {
	    rval = mpf_ILLlib_addrow (lp, B, rmatcnt[i], rmatind + rmatbeg[i],
		rmatval + rmatbeg[i], rhs[i], sense[i], rng,
		names[i]);
	} else {
	    rval = mpf_ILLlib_addrow (lp, B, rmatcnt[i], rmatind + rmatbeg[i],
		rmatval + rmatbeg[i], rhs[i], sense[i], rng, 0);
	}
	ILL_CLEANUP_IF (rval);
    }


    if (B && B->rownorms && (factorok && *factorok == 0)) {
	lp->matbeg = lp->O->A.matbeg;
	lp->matcnt = lp->O->A.matcnt;
	lp->matind = lp->O->A.matind;
	lp->matval = lp->O->A.matval;
	lp->nrows = lp->O->nrows;
	lp->ncols = lp->O->ncols;
	lp->bz = lp->O->rhs;
	lp->nnbasic = lp->ncols - lp->nrows;

	rval = mpf_ILLbasis_load (lp, B);
	ILL_CLEANUP_IF (rval);

	if (lp->f)
	    mpf_ILLfactor_free_factor_work (lp->f);

	rval = mpf_ILLbasis_factor (lp, &bsing);
	ILL_CLEANUP_IF (rval);
	*factorok = 1;

	if (B->rownorms_size < lp->O->nrows)
	    mpf_EGlpNumReallocArray (&(B->rownorms), lp->O->nrows);

	ILL_SAFE_MALLOC (rindi, lp->O->nrows /* num */ , int);

	for (i = 0; i < num; i++) {
	    rindi[i] = lp->O->nrows - num + i;
	}

	rval = mpf_ILLprice_get_dsteep_norms (lp, num, rindi,
	    B->rownorms + lp->O->nrows - num);
	ILL_CLEANUP_IF (rval);
    }
    if (factorok != 0 && badfactor == 1) {
	*factorok = 0;
    }
CLEANUP:

    ILL_IFFREE (bcnt, int);
    ILL_IFFREE (bbeg, int);
    ILL_IFFREE (bindi, int);
    mpf_EGlpNumFreeArray (bval);
    ILL_IFFREE (imap, int);
    ILL_IFFREE (jstat, int);
    ILL_IFFREE (rindi, int);
    mpf_EGlpNumClearVar (rng);
    ILL_RETURN (rval, "mpf_ILLlib_addrows");
}

int mpf_ILLlib_addrow (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int cnt,
      int *ind,
      mpf_t * val,
      mpf_t rhs,
      int sense,
      mpf_t range,
      const char *name)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLmatrix *A;
    int i, nrows, ncols;
    char buf[ILL_namebufsize];
    int tind[1];
    mpf_t tval[1];
    int *tempind = 0;
    int pind, hit;
    mpf_EGlpNumInitVar (tval[0]);

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_addrow called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    A = &qslp->A;

    if (qslp->rA) {		/* After an addrow call, needs to be updated */
	mpf_ILLlp_rows_clear (qslp->rA);
	ILL_IFFREE (qslp->rA, mpf_ILLlp_rows);
    }
    if (qslp->sinfo) {		/* Presolve LP is no longer valid, free the
				   data */
	mpf_ILLlp_sinfo_free (qslp->sinfo);
	ILL_IFFREE (qslp->sinfo, mpf_ILLlp_sinfo);
    }
    nrows = qslp->nrows;
    ncols = qslp->ncols;

    /* If the row has a range, create the rangeval array if needed  */

    if (sense == 'R' && !(qslp->rangeval) && qslp->rowsize > 0) {
	qslp->rangeval = mpf_EGlpNumAllocArray (qslp->rowsize);
	for (i = 0; i < qslp->nrows; i++) {
	    mpf_EGlpNumZero (qslp->rangeval[i]);
	}
    }
    /* Add the row to the row structures */

    if (qslp->rowsize < nrows + 1) {
	mpf_EGlpNumReallocArray (&(qslp->rhs), qslp->rowsize + mpf_EXTRA_ROWS);
	qslp->sense = EGrealloc (qslp->sense,
	    sizeof (char) * (qslp->rowsize + mpf_EXTRA_ROWS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(qslp->sense), qslp->rowsize + mpf_EXTRA_ROWS, sizeof (char));
			ILL_CLEANUP_IF (rval);
	*/

	qslp->rowmap = EGrealloc (qslp->rowmap,
	    sizeof (int) * (qslp->rowsize + mpf_EXTRA_ROWS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(qslp->rowmap), qslp->rowsize + mpf_EXTRA_ROWS, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/

	if (qslp->rangeval || sense == 'R')
	    mpf_EGlpNumReallocArray (&(qslp->rangeval), qslp->rowsize + mpf_EXTRA_ROWS);

	qslp->rownames = EGrealloc (qslp->rownames,
	    sizeof (char *) * (qslp->rowsize + mpf_EXTRA_ROWS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(qslp->rownames), qslp->rowsize + mpf_EXTRA_ROWS, sizeof (char *));
			ILL_CLEANUP_IF (rval);
	*/
	qslp->rowsize += mpf_EXTRA_ROWS;
    }
    mpf_EGlpNumCopy (qslp->rhs[nrows], rhs);
    qslp->sense[nrows] = sense;
    qslp->rowmap[nrows] = ncols;/* this will be the new logical */
    if (qslp->rangeval) {
	if (sense == 'R')
	    mpf_EGlpNumCopy (qslp->rangeval[nrows], range);
	else
	    mpf_EGlpNumZero (qslp->rangeval[nrows]);
    }
    ILL_FAILtrue (qslp->rownames == NULL, "must always be non NULL");
    mpf_ILLlib_findName (qslp, 1 /* row */ , name, nrows, buf);
    mpf_ILL_UTIL_STR (qslp->rownames[nrows], buf);
    ILLsymboltab_register (&qslp->rowtab, buf, qslp->nrows, &pind, &hit);
    ILL_FAILfalse (hit == 0, "must be new");


    /* Add the logical variable to the column structures */

    if (qslp->colsize < ncols + 1) {
	mpf_EGlpNumReallocArray (&(qslp->lower), qslp->colsize + mpf_EXTRA_COLS);
	mpf_EGlpNumReallocArray (&(qslp->upper), qslp->colsize + mpf_EXTRA_COLS);
	mpf_EGlpNumReallocArray (&(qslp->obj), qslp->colsize + mpf_EXTRA_COLS);
	qslp->colsize += mpf_EXTRA_COLS;
    }
    mpf_EGlpNumZero (qslp->obj[ncols]);
    mpf_EGlpNumZero (qslp->lower[ncols]);
    if (sense == 'E') {
	mpf_EGlpNumZero (qslp->upper[ncols]);	/* Artificial */
    } else if (sense == 'R') {
	mpf_EGlpNumCopy (qslp->upper[ncols], range);	/* Range      */
    } else {
	mpf_EGlpNumCopy (qslp->upper[ncols], mpf_ILL_MAXDOUBLE);	/* Slack      */
    }

    /* Add new row and new logical col to matrix */

    /* Need to map the structural indices to their proper place */

    if (cnt) {
	ILL_SAFE_MALLOC (tempind, cnt, int);
	for (i = 0; i < cnt; i++) {
	    tempind[i] = qslp->structmap[ind[i]];
	}
    }
    rval = mpf_matrix_addrow (A, cnt, tempind, val);
    ILL_CLEANUP_IF (rval);

    tind[0] = nrows;
    mpf_EGlpNumOne (*tval);
    if (sense == 'G' || sense == 'R')
	mpf_EGlpNumSign (*tval);

    rval = mpf_matrix_addcol (A, 1, tind, tval);
    ILL_CLEANUP_IF (rval);

    if (B != 0) {
	B->rstat = EGrealloc (B->rstat, sizeof (char) * (nrows + 1));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(B->rstat), nrows + 1, sizeof (char));
			ILL_CLEANUP_IF (rval);
	*/
	B->rstat[nrows] = QS_ROW_BSTAT_BASIC;
    }
#if 0
    lp->basisid = -1;		/* To get optimizer to reload the basis */
#endif

    qslp->ncols++;
    qslp->nrows++;
    qslp->nzcount += (cnt + 1);

    if (B != 0) {
	B->nrows++;
    }
CLEANUP:
    ILL_IFFREE (tempind, int);
    mpf_EGlpNumClearVar (tval[0]);
    ILL_RETURN (rval, "mpf_ILLlib_addrow");
}

int mpf_ILLlib_delrows (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      mpf_ILLlp_cache * C,
      int num,
      int *dellist,
      int *basis_ok,
      int *cache_ok)
{
    int rval = 0;
    int i, j, k, nrows, ncols, nstruct, spot, dk, bok = 0, cok = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLmatrix *A;
    char *rowmark = 0;
    char *colmark = 0;
    int *newrowindex = 0;
    int *newcolindex = 0;
    int *ind, *beg, *cnt;
    mpf_t *val;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_delrows called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (num <= 0) {
	if (basis_ok)
	    *basis_ok = 1;
	if (cache_ok)
	    *cache_ok = 1;
	ILL_CLEANUP;
    }
    if (basis_ok)
	*basis_ok = 0;
    if (cache_ok)
	*cache_ok = 0;

    qslp = lp->O;
    A = &qslp->A;

    if (qslp->rA) {		/* After a delrow call, needs to be updated */
	mpf_ILLlp_rows_clear (qslp->rA);
	ILL_IFFREE (qslp->rA, mpf_ILLlp_rows);
    }
    nrows = A->matrows;
    ncols = A->matcols;
    ind = A->matind;
    beg = A->matbeg;
    cnt = A->matcnt;
    val = A->matval;
    nstruct = qslp->nstruct;

    ILL_SAFE_MALLOC (rowmark, nrows, char);

    for (i = 0; i < nrows; i++) {
	rowmark[i] = 0;
    }
    for (i = 0; i < num; i++) {
	rowmark[dellist[i]] = 1;
    }


    /* Try to update the basis */

    if (B) {
	bok = 1;
	cok = 1;
	for (i = 0; i < num; i++) {
	    j = dellist[i];
	    if (B->rstat[j] == QS_ROW_BSTAT_LOWER ||
		B->rstat[j] == QS_ROW_BSTAT_UPPER) {
		bok = 0;
		break;
	    }
	    if (C && mpf_EGlpNumIsLess (mpf_DFEAS_TOLER, C->pi[j])) {
		/*
		                printf ("XXXX: Postive pi (%f) at basic row\n", C->pi[j]);
		                fflush (stdout);
		*/
		cok = 0;
	    }
	}
	if (bok == 1) {
	    mpf_EGlpNumFreeArray (B->colnorms);
	    if (B->rownorms) {
		for (i = 0, k = 0; i < nstruct; i++) {
		    if (B->cstat[i] == QS_COL_BSTAT_BASIC)
			k++;
		}
		for (i = 0, j = k; i < nrows; i++) {
		    if (B->rstat[i] == QS_ROW_BSTAT_BASIC) {
			if (rowmark[i] == 0) {
			    mpf_EGlpNumCopy (B->rownorms[k++], B->rownorms[j]);
			}
			j++;
		    }
		}
		if (k != nrows - num) {
		    fprintf (stderr, "error in  mpf_ILLlib_delrows\n");
		    rval = 1;
		    ILL_CLEANUP;
		}
	    }
	    for (i = 0, j = 0; i < nrows; i++) {
		if (rowmark[i] == 0) {
		    B->rstat[j++] = B->rstat[i];
		}
	    }
	    B->nrows = j;

	    if (C && cok == 1) {
		for (i = 0, j = 0; i < nrows; i++) {
		    if (rowmark[i] == 0) {
			mpf_EGlpNumCopy (C->pi[j], C->pi[i]);
			mpf_EGlpNumCopy (C->slack[j++], C->slack[i]);
		    }
		}
		C->nrows = j;
		if (cache_ok)
		    *cache_ok = 1;
	    }
	    if (basis_ok)
		*basis_ok = 1;
	}
    }
    ILL_SAFE_MALLOC (newrowindex, nrows, int);


    /* Delete the marked rows */

    ILL_FAILtrue (qslp->rownames == NULL, "must always be non NULL");
    for (i = 0, j = 0; i < nrows; i++) {
	if (rowmark[i] == 0) {
	    if (i != j) {
		mpf_EGlpNumCopy (qslp->rhs[j], qslp->rhs[i]);
		qslp->sense[j] = qslp->sense[i];
		if (qslp->rangeval)
		    mpf_EGlpNumCopy (qslp->rangeval[j], qslp->rangeval[i]);
		if (qslp->rownames)
		    qslp->rownames[j] = qslp->rownames[i];
	    }
	    newrowindex[i] = j++;
	} else {
	    if (qslp->rownames) {
		rval = ILLsymboltab_delete (&qslp->rowtab, qslp->rownames[i]);
		ILL_CLEANUP_IF (rval);
		ILL_IFFREE (qslp->rownames[i], char);
	    }
	}
    }


    /* Delete the logicals */

    ILL_SAFE_MALLOC (colmark, ncols, char);

    for (i = 0; i < ncols; i++) {
	colmark[i] = 0;
    }
    for (i = 0; i < num; i++) {
	colmark[qslp->rowmap[dellist[i]]] = 1;
    }

    rval = mpf_delcols_work (lp, colmark);
    ILL_CLEANUP_IF (rval);

    A->matcols -= num;
    qslp->ncols -= num;


    /* Pack the rowmap  */

    for (i = 0, j = 0; i < nrows; i++) {
	if (rowmark[i] == 0) {
	    qslp->rowmap[j++] = qslp->rowmap[i];
	}
    }

    /* Remove the entries to deleted rows, and update the indices */

    for (i = 0; i < ncols - num; i++) {
	dk = 0;
	spot = beg[i];
	for (j = 0; j < cnt[i]; j++) {
	    if (rowmark[ind[beg[i] + j]] == 1) {
		dk++;
	    } else {
		mpf_EGlpNumCopy (val[spot], val[beg[i] + j]);
		ind[spot] = newrowindex[ind[beg[i] + j]];
		spot++;
	    }
	}
	for (; spot < beg[i] + cnt[i]; spot++) {
	    ind[spot] = -1;
	}

	cnt[i] -= dk;
	if (cnt[i] == 0) {
	    ind[beg[i]] = 1;	/* we always mark the empty cols */
	}
    }

    A->matrows -= num;
    qslp->nrows -= num;

#if 0
    lp->basisid = -1;		/* To get optimizer to reload the basis */
#endif

CLEANUP:

    ILL_IFFREE (rowmark, char);
    ILL_IFFREE (colmark, char);
    ILL_IFFREE (newcolindex, int);
    ILL_IFFREE (newrowindex, int);
    ILL_RETURN (rval, "mpf_ILLlib_delrows");
}

int mpf_ILLlib_delcols (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int num,
      int *dellist,
      int *basis_ok)
{
    int rval = 0;
    int i, j, bok, ncols;
    char *colmark = 0;
    mpf_ILLlpdata *qslp;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_delcols called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (basis_ok)
	*basis_ok = 0;

    if (num <= 0) {
	*basis_ok = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    ncols = qslp->A.matcols;

    if (qslp->rA) {		/* After a delcol call, needs to be updated */
	mpf_ILLlp_rows_clear (qslp->rA);
	ILL_IFFREE (qslp->rA, mpf_ILLlp_rows);
    }
    ILL_SAFE_MALLOC (colmark, ncols, char);

    for (i = 0; i < ncols; i++) {
	colmark[i] = 0;
    }
    for (i = 0; i < num; i++) {
	colmark[qslp->structmap[dellist[i]]] = 1;
    }

    if (B) {
	B->nstruct -= num;
	bok = 1;
	for (i = 0; i < num; i++) {
	    j = dellist[i];
	    if (B->cstat[j] == QS_COL_BSTAT_BASIC) {
		bok = 0;
		break;
	    }
	}
	if (bok == 1) {
	    mpf_EGlpNumFreeArray (B->colnorms);
	    for (i = 0, j = 0; i < qslp->nstruct; i++) {
		if (colmark[qslp->structmap[i]] == 0) {
		    B->cstat[j++] = B->cstat[i];
		}
	    }
	    if (basis_ok)
		*basis_ok = 1;
	}
    }
    rval = mpf_delcols_work (lp, colmark);
    ILL_CLEANUP_IF (rval);


    qslp->A.matcols -= num;
    qslp->ncols -= num;
    qslp->nstruct -= num;

#if 0
    lp->basisid = -1;		/* To get optimizer to reload the basis */
#endif

CLEANUP:

    ILL_IFFREE (colmark, char);
    ILL_RETURN (rval, "mpf_ILLlib_delcols");
}

static int mpf_delcols_work (mpf_lpinfo * lp,
      char *colmark)
{
    int rval = 0;
    int i, j, k, nrows, ncols;
    mpf_ILLlpdata *qslp;
    mpf_ILLmatrix *A;
    int *newcolindex = 0;
    int *ind, *beg, *cnt;

    /* Allows logicals to be deleted, to handle call from delcols. */

    qslp = lp->O;
    A = &qslp->A;
    nrows = A->matrows;
    ncols = A->matcols;
    ind = A->matind;
    beg = A->matbeg;
    cnt = A->matcnt;

    ILL_SAFE_MALLOC (newcolindex, ncols, int);

    /* Delete the columns */

    for (i = 0, j = 0; i < ncols; i++) {
	if (colmark[i] == 0) {
	    if (i != j) {
		beg[j] = beg[i];
		cnt[j] = cnt[i];
		mpf_EGlpNumCopy (qslp->obj[j], qslp->obj[i]);
		mpf_EGlpNumCopy (qslp->lower[j], qslp->lower[i]);
		mpf_EGlpNumCopy (qslp->upper[j], qslp->upper[i]);
	    }
	    newcolindex[i] = j++;
	} else {
	    for (k = 0; k < cnt[i]; k++) {
		ind[beg[i] + k] = -1;
	    }
	    newcolindex[i] = -1;
	}
    }

    /* Update the struct arrays */

    for (i = 0, j = 0; i < qslp->nstruct; i++) {
	k = qslp->structmap[i];
	if (colmark[k] == 0) {
	    qslp->structmap[j] = newcolindex[k];
	    qslp->colnames[j] = qslp->colnames[i];
	    if (qslp->intmarker)
		qslp->intmarker[j] = qslp->intmarker[i];
	    j++;
	} else {
	    rval = ILLsymboltab_delete (&qslp->coltab, qslp->colnames[i]);
	    ILL_CLEANUP_IF (rval);
	    ILL_IFFREE (qslp->colnames[i], char);
	}
    }

    /* Update the rowmap: note if logicals deleted, map will be -1 */

    for (i = 0; i < nrows; i++) {
	qslp->rowmap[i] = newcolindex[qslp->rowmap[i]];
    }

CLEANUP:

    ILL_IFFREE (newcolindex, int);
    ILL_RETURN (rval, "mpf_delcols_work");
}

int mpf_ILLlib_chgcoef (mpf_lpinfo * lp,
      int rowindex,
      int colindex,
      mpf_t coef)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLmatrix *A;
    int nrows, nstruct, j;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_chgcoef called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    A = &qslp->A;

    nrows = qslp->nrows;
    nstruct = qslp->nstruct;

    if (rowindex < 0 || rowindex >= nrows || colindex < 0 || colindex >= nstruct) {
	fprintf (stderr, "mpf_ILLlib_chgcoef called with out-of-range index\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (qslp->rA) {		/* After a chgcoef call, needs to be updated */
	mpf_ILLlp_rows_clear (qslp->rA);
	ILL_IFFREE (qslp->rA, mpf_ILLlp_rows);
    }
    if (qslp->sinfo) {		/* Presolve LP is no longer valid, free the
				   data */
	mpf_ILLlp_sinfo_free (qslp->sinfo);
	ILL_IFFREE (qslp->sinfo, mpf_ILLlp_sinfo);
    }
    j = qslp->structmap[colindex];

    rval = mpf_matrix_addcoef (lp, A, rowindex, j, coef);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_chgcoef");
}

int mpf_ILLlib_chgsense (mpf_lpinfo * lp,
      int num,
      int *rowlist,
      char *sense)
{
    int rval = 0;
    int i, j, k;
    mpf_ILLlpdata *qslp = lp->O;
    mpf_ILLmatrix *A = &(lp->O->A);

    for (i = 0; i < num; i++) {
	j = qslp->rowmap[rowlist[i]];
	if (A->matcnt[j] != 1) {
	    fprintf (stderr, "logical variable is not a singleton\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
	k = A->matbeg[j];
	switch (sense[i]) {
	case 'E':		/* Artificial */
	    qslp->sense[rowlist[i]] = 'E';
	    mpf_EGlpNumZero (qslp->lower[j]);
	    mpf_EGlpNumZero (qslp->upper[j]);
	    mpf_EGlpNumOne (A->matval[k]);
	    break;
	case 'G':		/* Surplus   */
	    qslp->sense[rowlist[i]] = 'G';
	    mpf_EGlpNumZero (qslp->lower[j]);
	    mpf_EGlpNumCopy (qslp->upper[j], mpf_ILL_MAXDOUBLE);
	    mpf_EGlpNumOne (A->matval[k]);
	    mpf_EGlpNumSign (A->matval[k]);
	    break;
	case 'L':		/* Slack     */
	    qslp->sense[rowlist[i]] = 'L';
	    mpf_EGlpNumZero (qslp->lower[j]);
	    mpf_EGlpNumCopy (qslp->upper[j], mpf_ILL_MAXDOUBLE);
	    mpf_EGlpNumOne (A->matval[k]);
	    break;
	default:
	    fprintf (stderr, "illegal sense %c in mpf_ILLlib_chgsense\n", sense[i]);
	    rval = 1;
	    ILL_CLEANUP;
	}
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_chgsense");
}

int mpf_ILLlib_newcol (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      mpf_t obj,
      mpf_t lower,
      mpf_t upper,
      const char *name,
      int factorok)
{
    int rval = 0;

    rval = mpf_ILLlib_addcol (lp, B, 0, 0, 0, obj, lower, upper, name, factorok);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_newcol");
}

int mpf_ILLlib_newcols (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int num,
      mpf_t * obj,
      mpf_t * lower,
      mpf_t * upper,
      const char **names,
      int factorok)
{
    int rval = 0;
    int *cmatcnt = 0;
    int *cmatbeg = 0;
    int i;

    ILL_SAFE_MALLOC (cmatcnt, num, int);

    ILL_SAFE_MALLOC (cmatbeg, num, int);

    for (i = 0; i < num; i++) {
	cmatcnt[i] = 0;
	cmatbeg[i] = 0;
    }

    rval = mpf_ILLlib_addcols (lp, B, num, cmatcnt, cmatbeg, 0,
	0, obj, lower, upper, names, factorok);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    ILL_IFFREE (cmatcnt, int);
    ILL_IFFREE (cmatbeg, int);

    ILL_RETURN (rval, "mpf_ILLlib_newcols");
}

int mpf_ILLlib_addcols (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int num,
      int *cmatcnt,
      int *cmatbeg,
      int *cmatind,
      mpf_t * cmatval,
      mpf_t * obj,
      mpf_t * lower,
      mpf_t * upper,
      const char **names,
      int factorok)
{
    int rval = 0;
    int i;

    for (i = 0; i < num; i++) {
	if (names) {
	    rval = mpf_ILLlib_addcol (lp, B, cmatcnt[i], cmatind + cmatbeg[i],
		cmatval + cmatbeg[i], obj[i], lower[i],
		upper[i], names[i], factorok);
	} else {
	    rval = mpf_ILLlib_addcol (lp, B, cmatcnt[i], cmatind + cmatbeg[i],
		cmatval + cmatbeg[i], obj[i], lower[i],
		upper[i], 0, factorok);
	}
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_addcols");
}

int mpf_ILLlib_addcol (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      int cnt,
      int *ind,
      mpf_t * val,
      mpf_t obj,
      mpf_t lower,
      mpf_t upper,
      const char *name,
      int factorok)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLmatrix *A;
    int ncols;
    char buf[ILL_namebufsize];
    int pind, hit;
    mpf_t l, u;
    mpf_EGlpNumInitVar (l);
    mpf_EGlpNumInitVar (u);

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_addcol called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    A = &qslp->A;
    ncols = qslp->ncols;

    if (qslp->rA) {		/* After an addcol call, needs to be updated */
	mpf_ILLlp_rows_clear (qslp->rA);
	ILL_IFFREE (qslp->rA, mpf_ILLlp_rows);
    }
    if (qslp->sinfo) {		/* Presolve LP is no longer valid, free the
				   data */
	mpf_ILLlp_sinfo_free (qslp->sinfo);
	ILL_IFFREE (qslp->sinfo, mpf_ILLlp_sinfo);
    }
    /* Add the new variable to the column structures */

    if (qslp->colsize < ncols + 1) {
	mpf_EGlpNumReallocArray (&(qslp->lower), qslp->colsize + mpf_EXTRA_COLS);
	mpf_EGlpNumReallocArray (&(qslp->upper), qslp->colsize + mpf_EXTRA_COLS);
	mpf_EGlpNumReallocArray (&(qslp->obj), qslp->colsize + mpf_EXTRA_COLS);
	qslp->colsize += mpf_EXTRA_COLS;
    }
    mpf_EGlpNumCopy (qslp->obj[ncols], obj);
    mpf_EGlpNumCopy (qslp->lower[ncols], lower);
    mpf_EGlpNumCopy (qslp->upper[ncols], upper);

    /* Add the variable to the structural arrays */

    if (qslp->structsize < qslp->nstruct + 1) {
	qslp->structmap = EGrealloc (qslp->structmap,
	    sizeof (int) * (qslp->structsize + mpf_EXTRA_COLS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(qslp->structmap), qslp->structsize + mpf_EXTRA_COLS, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/

	qslp->colnames = EGrealloc (qslp->colnames,
	    sizeof (char *) * (qslp->structsize + mpf_EXTRA_COLS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(qslp->colnames), qslp->structsize + mpf_EXTRA_COLS, sizeof (char *));
			ILL_CLEANUP_IF (rval);
	*/

	if (qslp->intmarker) {
	    qslp->intmarker = EGrealloc (qslp->intmarker,
		sizeof (char) * (qslp->structsize + mpf_EXTRA_COLS));
	    /*
				    rval = ILLutil_reallocrus_count ((void **) &(qslp->intmarker), qslp->structsize + mpf_EXTRA_COLS, sizeof (char));
				    ILL_CLEANUP_IF (rval);
	    */
	}
	qslp->structsize += mpf_EXTRA_COLS;
    }
    qslp->structmap[qslp->nstruct] = ncols;
    if (qslp->intmarker) {
	/* NOTE: If we want to add integer variables, this is the place. */
	qslp->intmarker[qslp->nstruct] = (char) 0;
    }
    ILL_FAILtrue (qslp->colnames == NULL, "must always be non NULL");
    mpf_ILLlib_findName (qslp, 0 /* isRow */ , name, qslp->nstruct, buf);
    ILLsymboltab_register (&qslp->coltab, buf, qslp->nstruct, &pind, &hit);
    ILL_FAILfalse ((pind == qslp->nstruct) && (hit == 0), "must be new");
    mpf_ILL_UTIL_STR (qslp->colnames[qslp->nstruct], buf);


    /* Add col to the matrix */

    rval = mpf_matrix_addcol (A, cnt, ind, val);
    ILL_CLEANUP_IF (rval);


    if (B) {
	B->cstat = EGrealloc (B->cstat, sizeof (char) * (qslp->nstruct + 1));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(B->cstat), qslp->nstruct + 1, sizeof (char));
			ILL_CLEANUP_IF (rval);
	*/
	if (mpf_EGlpNumIsEqual (lower, mpf_ILL_MINDOUBLE, mpf_oneLpNum) &&
	    mpf_EGlpNumIsEqual (upper, mpf_ILL_MAXDOUBLE, mpf_oneLpNum)) {
	    B->cstat[qslp->nstruct] = QS_COL_BSTAT_FREE;
	} else if (mpf_EGlpNumIsEqual (upper, mpf_ILL_MAXDOUBLE, mpf_oneLpNum)) {
	    B->cstat[qslp->nstruct] = QS_COL_BSTAT_LOWER;
	}
	/* else if (lower == mpf_ILL_MAXDOUBLE) */
	else if (mpf_EGlpNumIsEqual (lower, mpf_ILL_MAXDOUBLE, mpf_oneLpNum)) {
	    B->cstat[qslp->nstruct] = QS_COL_BSTAT_UPPER;
	} else {
	    /* l = fabs (lower); u = fabs (upper); */
	    mpf_EGlpNumCopyAbs (l, lower);
	    mpf_EGlpNumCopyAbs (u, upper);
	    if (mpf_EGlpNumIsLess (l, u)) {
		B->cstat[qslp->nstruct] = QS_COL_BSTAT_LOWER;
	    } else {
		B->cstat[qslp->nstruct] = QS_COL_BSTAT_UPPER;
	    }
	}

	/* UPDATE THE PINF mpf_PRIMAL NORMS */
	mpf_EGlpNumFreeArray (B->colnorms);
    }
    if (factorok == 0) {
#if 0
	lp->basisid = -1;	/* To get optimizer to reload the basis */
#endif
    } else {
	if (!lp->nbaz || !lp->vindex || !lp->vstat) {
	    fprintf (stderr, "ERROR: factorok set without a current basis\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
	lp->nbaz = EGrealloc (lp->nbaz, sizeof (int) * (qslp->nstruct + 1));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(lp->nbaz), qslp->nstruct + 1, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/

	lp->vindex = EGrealloc (lp->vindex, sizeof (int) * (qslp->ncols + 1));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(lp->vindex), qslp->ncols + 1, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/

	lp->vstat = EGrealloc (lp->vstat, sizeof (int) * (qslp->ncols + 1));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(lp->vstat), qslp->ncols + 1, sizeof (int));
	*/


	lp->nbaz[qslp->nstruct] = qslp->ncols;
	lp->vindex[qslp->ncols] = qslp->nstruct;

	if (mpf_EGlpNumIsEqual (lower, mpf_ILL_MINDOUBLE, mpf_oneLpNum) &&
	    mpf_EGlpNumIsEqual (upper, mpf_ILL_MAXDOUBLE, mpf_oneLpNum)) {
	    lp->vstat[qslp->ncols] = STAT_ZERO;
	} else if (mpf_EGlpNumIsEqual (upper, mpf_ILL_MAXDOUBLE, mpf_oneLpNum)) {
	    lp->vstat[qslp->ncols] = STAT_LOWER;
	}
	/* else if (lower == mpf_ILL_MAXDOUBLE) */
	else if (mpf_EGlpNumIsEqual (lower, mpf_ILL_MAXDOUBLE, mpf_oneLpNum)) {
	    lp->vstat[qslp->ncols] = STAT_UPPER;
	} else {
	    /* l = fabs (lower); u = fabs (upper); */
	    mpf_EGlpNumCopyAbs (l, lower);
	    mpf_EGlpNumCopyAbs (u, upper);
	    if (mpf_EGlpNumIsLess (l, u)) {
		lp->vstat[qslp->ncols] = STAT_LOWER;
	    } else {
		lp->vstat[qslp->ncols] = STAT_UPPER;
	    }
	}
    }


    qslp->ncols++;
    qslp->nstruct++;
    (qslp->nzcount) += cnt;

    if (B) {
	B->nstruct++;
    }
CLEANUP:
    mpf_EGlpNumClearVar (l);
    mpf_EGlpNumClearVar (u);
    ILL_RETURN (rval, "mpf_ILLlib_addcol");
}

static int mpf_matrix_addrow (mpf_ILLmatrix * A,
      int rowcnt,
      int *rowind,
      mpf_t * rowval)
{
    int rval = 0;
    int i, j, k, ind, memo, stop, delta = 0;

    /* matsize will be the length of the array.                   */
    /* matfree will keep track of the free space at end of array. */

    for (i = 0; i < rowcnt; i++) {
	if (rowind[i] >= A->matcols || rowind[i] < 0) {
	    fprintf (stderr, "illegal col index in mpf_matrix_addrow\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
    }

    for (i = 0; i < rowcnt; i++) {
	j = rowind[i];
	if (A->matcnt[j] > 0 &&
	    (A->matbeg[j] + A->matcnt[j] + 1 > A->matsize ||
		A->matind[A->matbeg[j] + A->matcnt[j]] != -1)) {
	    delta += (A->matcnt[j] + 2);	/* 1 for the new coef and 1
						   for */
	    /* an extra space               */
	}
    }

    if (delta < A->matfree) {
	for (i = 0; i < rowcnt; i++) {
	    j = rowind[i];
	    if (A->matcnt[j] == 0) {
		A->matind[A->matbeg[j]] = A->matrows;
		mpf_EGlpNumCopy (A->matval[A->matbeg[j]], rowval[i]);
		A->matcnt[j] = 1;
	    } else if (A->matind[A->matbeg[j] + A->matcnt[j]] == -1) {
		/* Since A->matfree is positive, we know that we are not */
		/* sitting at the end of the array.                      */
		A->matind[A->matbeg[j] + A->matcnt[j]] = A->matrows;
		mpf_EGlpNumCopy (A->matval[A->matbeg[j] + A->matcnt[j]], rowval[i]);
		if ((A->matbeg[j] + A->matcnt[j]) == (A->matsize - A->matfree)) {
		    A->matfree--;	/* at end of used space */
		}
		(A->matcnt[j])++;
	    } else {
		ind = A->matsize - A->matfree + 1;	/* leave space for -1 */
		memo = ind;
		stop = A->matbeg[j] + A->matcnt[j];
		for (k = A->matbeg[j]; k < stop; k++) {
		    if (ind >= A->matsize) {
			printf ("WHAT: %d, %d\n", A->matsize, ind);
			fflush (stdout);
			exit (1);
		    }
		    A->matind[ind] = A->matind[k];
		    mpf_EGlpNumCopy (A->matval[ind], A->matval[k]);
		    A->matind[k] = -1;
		    ind++;
		}
		A->matind[ind] = A->matrows;
		mpf_EGlpNumCopy (A->matval[ind], rowval[i]);
		A->matbeg[j] = memo;
		(A->matcnt[j])++;
		(A->matfree) -= (A->matcnt[j] + 1);
	    }
	}
    } else {
	rval = mpf_matrix_addrow_end (A, A->matrows, rowcnt, rowind, rowval);
	ILL_CLEANUP_IF (rval);
    }
    A->matrows++;

CLEANUP:

    ILL_RETURN (rval, "mpf_matrix_addrow");
}

static int mpf_matrix_addrow_end (mpf_ILLmatrix * A,
      int row,
      int rowcnt,
      int *rowind,
      mpf_t * rowval)
{
    int rval = 0;
    int i, j, k, start, stop, total;
    int *newbeg = 0;
    int *newind = 0;
    mpf_t *newval = 0;
    int ncols = A->matcols;

    if (A->matcolsize > 0) {
	ILL_SAFE_MALLOC (newbeg, A->matcolsize, int);
    }
    ILL_SAFE_MALLOC (newind, A->matsize + rowcnt + mpf_EXTRA_MAT, int);
    newval = mpf_EGlpNumAllocArray (A->matsize + rowcnt + mpf_EXTRA_MAT);

    A->matsize += (rowcnt + mpf_EXTRA_MAT);

    for (i = 0; i < rowcnt; i++) {
	A->matcnt[rowind[i]]++;
    }
    for (total = 0, j = 0; j < ncols; j++) {
	newbeg[j] = total;
	if (A->matcnt[j] > 0)
	    total += A->matcnt[j];
	else
	    total += 1;
    }
    for (i = 0; i < rowcnt; i++) {
	A->matcnt[rowind[i]]--;
    }
    for (j = total; j < A->matsize; j++) {
	newind[j] = -1;
    }
    A->matfree = A->matsize - total;

    for (j = 0; j < ncols; j++) {
	if (A->matcnt[j] > 0) {
	    stop = A->matbeg[j] + A->matcnt[j];
	    start = newbeg[j];
	    for (k = A->matbeg[j]; k < stop; k++) {
		newind[start] = A->matind[k];
		mpf_EGlpNumCopy (newval[start], A->matval[k]);
		start++;
	    }
	} else {
	    newind[newbeg[j]] = 1;
	}
    }
    for (i = 0; i < rowcnt; i++) {
	j = rowind[i];
	newind[newbeg[j] + A->matcnt[j]] = row;
	mpf_EGlpNumCopy (newval[newbeg[j] + A->matcnt[j]], rowval[i]);
	(A->matcnt[j])++;
    }

    ILL_IFFREE (A->matbeg, int);
    ILL_IFFREE (A->matind, int);
    mpf_EGlpNumFreeArray (A->matval);

    A->matbeg = newbeg;
    A->matind = newind;
    A->matval = newval;

CLEANUP:

    if (rval) {
	ILL_IFFREE (newbeg, int);
	ILL_IFFREE (newind, int);
	mpf_EGlpNumFreeArray (newval);
    }
    ILL_RETURN (rval, "mpf_matrix_addrow_end");
}

static int mpf_matrix_addcoef (mpf_lpinfo * lp,
      mpf_ILLmatrix * A,
      int row,
      int col,
      mpf_t val)
{
    int i, k, delta, ind, stop, memo;
    int tind[1];
    mpf_t tval[1];
    int rval = 0;
    mpf_EGlpNumInitVar (tval[0]);
    mpf_EGlpNumCopy (tval[0], val);

    if (row >= A->matrows || row < 0) {
	fprintf (stderr, "illegal row index in mpf_matrix_addcoef\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (col >= A->matcols || col < 0) {
	fprintf (stderr, "illegal col index in mpf_matrix_addcoef\n");
	rval = 1;
	ILL_CLEANUP;
    }
    for (i = A->matbeg[col]; i < A->matbeg[col] + A->matcnt[col]; i++) {
	if (A->matind[i] == row) {
	    mpf_EGlpNumCopy (A->matval[i], val);
	    ILL_CLEANUP;
	}
    }

    /* The coef is new, we need to add it to A */

    lp->O->nzcount++;
    delta = A->matcnt[col] + 2;

    if (A->matcnt[col] == 0) {
	/* First entry, always a free space */
	A->matind[A->matbeg[col]] = row;
	mpf_EGlpNumCopy (A->matval[A->matbeg[col]], val);
	A->matcnt[col] = 1;
    } else if (A->matbeg[col] + A->matcnt[col] < A->matsize &&
	A->matind[A->matbeg[col] + A->matcnt[col]] == -1) {
	/* Free space in the column */
	A->matind[A->matbeg[col] + A->matcnt[col]] = row;
	mpf_EGlpNumCopy (A->matval[A->matbeg[col] + A->matcnt[col]], val);
	if ((A->matbeg[col] + A->matcnt[col]) == (A->matsize - A->matfree)) {
	    A->matfree--;
	}
	(A->matcnt[col])++;
    } else if (A->matfree > delta) {
	/* Enough space to move column to end of array */
	ind = A->matsize - A->matfree + 1;
	memo = ind;
	stop = A->matbeg[col] + A->matcnt[col];
	for (k = A->matbeg[col]; k < stop; k++) {
	    A->matind[ind] = A->matind[k];
	    mpf_EGlpNumCopy (A->matval[ind], A->matval[k]);
	    A->matind[k] = -1;
	    ind++;
	}
	A->matind[ind] = row;
	mpf_EGlpNumCopy (A->matval[ind], val);

	A->matbeg[col] = memo;
	(A->matcnt[col])++;
	(A->matfree) -= (A->matcnt[col] + 1);
    } else {
	/* Need to malloc space to move column to end of array */

	tind[0] = col;

	rval = mpf_matrix_addrow_end (A, row, 1, tind, tval);
	ILL_CLEANUP_IF (rval);
    }

CLEANUP:

    mpf_EGlpNumClearVar (tval[0]);
    ILL_RETURN (rval, "mpf_matrix_addcoef");
}

static int mpf_matrix_addcol (mpf_ILLmatrix * A,
      int colcnt,
      int *colind,
      mpf_t * colval)
{
    int rval = 0;
    int i, ind;

    for (i = 0; i < colcnt; i++) {
	if (colind[i] >= A->matrows || colind[i] < 0) {
	    fprintf (stderr, "illegal row index in mpf_matrix_addcol\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
    }

    if (A->matcolsize < A->matcols + 1) {
	A->matbeg = EGrealloc (A->matbeg, sizeof (int) * (A->matcolsize + mpf_EXTRA_COLS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(A->matbeg), A->matcolsize + mpf_EXTRA_COLS, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/

	A->matcnt = EGrealloc (A->matcnt, sizeof (int) * (A->matcolsize + mpf_EXTRA_COLS));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(A->matcnt), A->matcolsize + mpf_EXTRA_COLS, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/

	(A->matcolsize) += mpf_EXTRA_COLS;
    }
    if (A->matfree < colcnt + 1) {
	A->matind = EGrealloc (A->matind,
	    sizeof (int) * (A->matsize + colcnt + mpf_EXTRA_MAT + 1));
	/*
			rval = ILLutil_reallocrus_count ((void **) &(A->matind), A->matsize + colcnt + mpf_EXTRA_MAT + 1, sizeof (int));
			ILL_CLEANUP_IF (rval);
	*/
	mpf_EGlpNumReallocArray (&(A->matval), A->matsize + colcnt + mpf_EXTRA_MAT + 1);

	for (i = 0; i < colcnt + mpf_EXTRA_MAT + 1; i++) {
	    A->matind[A->matsize + i] = -1;
	}
	A->matsize += (colcnt + mpf_EXTRA_MAT + 1);
	A->matfree += (colcnt + mpf_EXTRA_MAT + 1);
    }
    ind = A->matsize - A->matfree;
    A->matbeg[A->matcols] = ind;
    A->matcnt[A->matcols] = colcnt;
    if (colcnt == 0) {
	A->matind[ind] = 1;	/* Dummy value to stop columns from stealing */
	/* this space in addrows.                    */
	A->matfree -= 1;
    } else {
	for (i = 0; i < colcnt; i++) {
	    mpf_EGlpNumCopy (A->matval[ind], colval[i]);
	    A->matind[ind] = colind[i];
	    ind++;
	}
	A->matfree -= colcnt;
    }
    A->matcols++;

CLEANUP:

    ILL_RETURN (rval, "mpf_matrix_addcol");
}

int mpf_ILLlib_getrows (mpf_lpinfo * lp,
      int num,
      int *rowlist,
      int **rowcnt,
      int **rowbeg,
      int **rowind,
      mpf_t ** rowval,
      mpf_t ** rhs,
      char **sense,
      char ***names)
{
    int rval = 0;
    int *allbeg = 0;
    int *allcnt = 0;
    int *allind = 0;
    mpf_t *allval = 0;
    int i, row, k, start, stop, len, tcnt, cnt = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLlp_rows lprows;

    if (rowcnt)
	*rowcnt = 0;
    if (rowbeg)
	*rowbeg = 0;
    if (rowind)
	*rowind = 0;
    if (rowval)
	*rowval = 0;
    if (rhs)
	*rhs = 0;
    if (sense)
	*sense = 0;
    if (names)
	*names = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getrows called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (!num)
	ILL_CLEANUP;

    qslp = lp->O;

    rval = mpf_ILLlp_rows_init (&lprows, qslp, 0);
    ILL_CLEANUP_IF (rval);
    allbeg = lprows.rowbeg;
    allcnt = lprows.rowcnt;
    allind = lprows.rowind;
    allval = lprows.rowval;

    for (i = 0; i < num; i++) {
	cnt += allcnt[rowlist[i]];
    }

    if (rowcnt) {
	ILL_SAFE_MALLOC (*rowcnt, num, int);
	for (i = 0; i < num; i++) {
	    (*rowcnt)[i] = allcnt[rowlist[i]];
	}
    }
    if (rowbeg) {
	ILL_SAFE_MALLOC (*rowbeg, num, int);
	tcnt = 0;
	for (i = 0; i < num; i++) {
	    (*rowbeg)[i] = tcnt;
	    tcnt += allcnt[rowlist[i]];
	}
    }
    if (cnt && rowind) {
	ILL_SAFE_MALLOC (*rowind, cnt, int);
	tcnt = 0;
	for (i = 0; i < num; i++) {
	    row = rowlist[i];
	    start = allbeg[row];
	    stop = start + allcnt[row];
	    for (k = start; k < stop; k++) {
		(*rowind)[tcnt++] = allind[k];
	    }
	}
    }
    if (cnt && rowval) {
	*rowval = mpf_EGlpNumAllocArray (cnt);
	tcnt = 0;
	for (i = 0; i < num; i++) {
	    row = rowlist[i];
	    start = allbeg[row];
	    stop = start + allcnt[row];
	    for (k = start; k < stop; k++) {
		mpf_EGlpNumCopy ((*rowval)[tcnt++], allval[k]);
	    }
	}
    }
    if (rhs) {
	*rhs = mpf_EGlpNumAllocArray (num);
	for (i = 0; i < num; i++) {
	    mpf_EGlpNumCopy ((*rhs)[i], qslp->rhs[rowlist[i]]);
	}
    }
    if (sense) {
	ILL_SAFE_MALLOC (*sense, num, char);
	for (i = 0; i < num; i++) {
	    (*sense)[i] = qslp->sense[rowlist[i]];
	}
    }
    if (names) {
	if (qslp->rownames == 0) {
	    fprintf (stderr, "LP does not have row names\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
	ILL_SAFE_MALLOC (*names, num, char *);
	for (i = 0; i < num; i++) {
	    (*names)[i] = 0;
	}
	for (i = 0; i < num; i++) {
	    len = strlen (qslp->rownames[rowlist[i]]) + 1;
	    ILL_SAFE_MALLOC ((*names)[i], len, char);
	    strcpy ((*names)[i], qslp->rownames[rowlist[i]]);
	}
    }
CLEANUP:

    ILL_IFFREE (allbeg, int);
    ILL_IFFREE (allcnt, int);
    ILL_IFFREE (allind, int);
    mpf_EGlpNumFreeArray (allval);

    if (rval) {
	if (rowcnt)
	    ILL_IFFREE (*rowcnt, int);
	if (rowbeg)
	    ILL_IFFREE (*rowbeg, int);
	if (rowind)
	    ILL_IFFREE (*rowind, int);
	if (rowval)
	    mpf_EGlpNumFreeArray (*rowval);
	if (rhs)
	    mpf_EGlpNumFreeArray (*rhs);
	if (sense)
	    ILL_IFFREE (*sense, char);
	if (names && (*names)) {
	    for (i = 0; i < num; i++) {
		ILL_IFFREE ((*names)[i], char);
	    }
	    ILL_IFFREE (*names, char *);
	}
    }
    ILL_RETURN (rval, "mpf_ILLlib_getrows");
}

int mpf_ILLlib_getcols (mpf_lpinfo * lp,
      int num,
      int *collist,
      int **colcnt,
      int **colbeg,
      int **colind,
      mpf_t ** colval,
      mpf_t ** obj,
      mpf_t ** lower,
      mpf_t ** upper,
      char ***names)
{
    int rval = 0;
    int i, col, k, start, stop, len, tcnt, cnt = 0;
    mpf_ILLlpdata *qslp;
    mpf_ILLmatrix *A;
    int *tlist = 0;

    if (colcnt)
	*colcnt = 0;
    if (colbeg)
	*colbeg = 0;
    if (colind)
	*colind = 0;
    if (colval)
	*colval = 0;
    if (lower)
	*lower = 0;
    if (upper)
	*upper = 0;
    if (obj)
	*obj = 0;
    if (names)
	*names = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getcols called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (!num)
	ILL_CLEANUP;

    qslp = lp->O;
    A = &(qslp->A);

    ILL_SAFE_MALLOC (tlist, num, int);
    for (i = 0; i < num; i++) {
	tlist[i] = qslp->structmap[collist[i]];
    }

    for (i = 0; i < num; i++) {
	cnt += A->matcnt[tlist[i]];
    }

    if (colcnt) {
	ILL_SAFE_MALLOC (*colcnt, num, int);
	for (i = 0; i < num; i++) {
	    (*colcnt)[i] = A->matcnt[tlist[i]];
	}
    }
    if (colbeg) {
	ILL_SAFE_MALLOC (*colbeg, num, int);
	tcnt = 0;
	for (i = 0; i < num; i++) {
	    (*colbeg)[i] = tcnt;
	    tcnt += A->matcnt[tlist[i]];
	}
    }
    if (cnt && colind) {
	ILL_SAFE_MALLOC (*colind, cnt, int);
	tcnt = 0;
	for (i = 0; i < num; i++) {
	    col = tlist[i];
	    start = A->matbeg[col];
	    stop = start + A->matcnt[col];
	    for (k = start; k < stop; k++) {
		(*colind)[tcnt++] = A->matind[k];
	    }
	}
    }
    if (cnt && colval) {
	*colval = mpf_EGlpNumAllocArray (cnt);
	tcnt = 0;
	for (i = 0; i < num; i++) {
	    col = tlist[i];
	    start = A->matbeg[col];
	    stop = start + A->matcnt[col];
	    for (k = start; k < stop; k++) {
		mpf_EGlpNumCopy ((*colval)[tcnt++], A->matval[k]);
	    }
	}
    }
    if (obj) {
	*obj = mpf_EGlpNumAllocArray (num);
	for (i = 0; i < num; i++) {
	    mpf_EGlpNumCopy ((*obj)[i], qslp->obj[tlist[i]]);
	}
    }
    if (lower) {
	*lower = mpf_EGlpNumAllocArray (num);
	for (i = 0; i < num; i++) {
	    mpf_EGlpNumCopy ((*lower)[i], qslp->lower[tlist[i]]);
	}
    }
    if (upper) {
	*upper = mpf_EGlpNumAllocArray (num);
	for (i = 0; i < num; i++) {
	    mpf_EGlpNumCopy ((*upper)[i], qslp->upper[tlist[i]]);
	}
    }
    if (names) {
	if (qslp->colnames == 0) {
	    fprintf (stderr, "LP does not have col names\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
	ILL_SAFE_MALLOC (*names, num, char *);
	for (i = 0; i < num; i++) {
	    (*names)[i] = 0;
	}
	for (i = 0; i < num; i++) {
	    len = strlen (qslp->colnames[collist[i]]) + 1;
	    ILL_SAFE_MALLOC ((*names)[i], len, char);
	    strcpy ((*names)[i], qslp->colnames[collist[i]]);
	}
    }
CLEANUP:

    if (rval) {
	if (colcnt)
	    ILL_IFFREE (*colcnt, int);
	if (colbeg)
	    ILL_IFFREE (*colbeg, int);
	if (colind)
	    ILL_IFFREE (*colind, int);
	if (colval)
	    mpf_EGlpNumFreeArray (*colval);
	if (obj)
	    mpf_EGlpNumFreeArray (*obj);
	if (lower)
	    mpf_EGlpNumFreeArray (*lower);
	if (upper)
	    mpf_EGlpNumFreeArray (*upper);
	if (names && (*names)) {
	    for (i = 0; i < num; i++) {
		ILL_IFFREE ((*names)[i], char);
	    }
	    ILL_IFFREE (*names, char *);
	}
    }
    ILL_IFFREE (tlist, int);

    ILL_RETURN (rval, "mpf_ILLlib_getcols");
}

int mpf_ILLlib_getobj (mpf_lpinfo * lp,
      mpf_t * obj)
{
    mpf_ILLlpdata *qslp;
    int nstruct, j;
    int rval = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getobj called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nstruct = qslp->nstruct;

    for (j = 0; j < nstruct; j++) {
	mpf_EGlpNumCopy (obj[j], qslp->obj[j]);
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_getobj");
}

int mpf_ILLlib_chgobj (mpf_lpinfo * lp,
      int indx,
      mpf_t coef)
{
    int rval = 0;
    int col;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_chgobj called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (indx < 0 || indx >= lp->O->nstruct) {
	fprintf (stderr, "mpf_ILLlib_chgrhs called with bad indx: %d\n", indx);
	rval = 1;
	ILL_CLEANUP;
    }
    if (lp->O->sinfo) {		/* Presolve LP is no longer valid, free the
				   data */
	mpf_ILLlp_sinfo_free (lp->O->sinfo);
	ILL_IFFREE (lp->O->sinfo, mpf_ILLlp_sinfo);
    }
    col = lp->O->structmap[indx];
    mpf_EGlpNumCopy (lp->O->obj[col], coef);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_chgobj");
}

int mpf_ILLlib_getrhs (mpf_lpinfo * lp,
      mpf_t * rhs)
{
    mpf_ILLlpdata *qslp;
    int nrows, i;
    int rval = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getrhs called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nrows = qslp->nrows;

    for (i = 0; i < nrows; i++) {
	mpf_EGlpNumCopy (rhs[i], qslp->rhs[i]);
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_getrhs");
}

int mpf_ILLlib_chgrhs (mpf_lpinfo * lp,
      int indx,
      mpf_t coef)
{
    int rval = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_chgrhs called without an lp\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (indx < 0 || indx >= lp->O->nrows) {
	fprintf (stderr, "mpf_ILLlib_chgrhs called with bad indx: %d\n", indx);
	rval = 1;
	ILL_CLEANUP;
    }
    if (lp->O->sinfo) {		/* Presolve LP is no longer valid, free the
				   data */
	mpf_ILLlp_sinfo_free (lp->O->sinfo);
	ILL_IFFREE (lp->O->sinfo, mpf_ILLlp_sinfo);
    }
    mpf_EGlpNumCopy (lp->O->rhs[indx], coef);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_chgrhs");
}

int mpf_ILLlib_rownames (mpf_lpinfo * lp,
      char **rownames)
{
    mpf_ILLlpdata *qslp;
    int nrows, len, i, rcount = 0;
    int rval = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_rownames called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (!rownames) {
	fprintf (stderr, "mpf_ILLlib_rownames called with NULL rownames\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nrows = qslp->nrows;

    if (qslp->rownames == 0) {
	fprintf (stderr, "LP does not have rownames assigned\n");
	rval = 1;
	ILL_CLEANUP;
    }
    for (i = 0; i < nrows; i++) {
	len = strlen (qslp->rownames[i]) + 1;
	ILL_SAFE_MALLOC (rownames[i], len, char);
	strcpy (rownames[i], qslp->rownames[i]);
	rcount++;
    }

CLEANUP:

    if (rval) {
	for (i = 0; i < rcount; i++) {
	    ILL_IFFREE (rownames[i], char);
	}
    }
    ILL_RETURN (rval, "mpf_ILLlib_rownames");
}

int mpf_ILLlib_getintflags (mpf_lpinfo * lp,
      int *intflags)
{
    int j, nstruct, rval = 0;
    mpf_ILLlpdata *qslp;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getintflags called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nstruct = qslp->nstruct;

    if (qslp->intmarker == 0) {
	for (j = 0; j < nstruct; j++) {
	    intflags[j] = 0;
	}
    } else {
	for (j = 0; j < nstruct; j++) {
	    if (qslp->intmarker[j]) {
		intflags[j] = 1;
	    } else {
		intflags[j] = 0;
	    }
	}
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_getintflags");
}

int mpf_ILLlib_colnames (mpf_lpinfo * lp,
      char **colnames)
{
    mpf_ILLlpdata *qslp;
    int nstruct, len, i, ccount = 0;
    int rval = 0;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_colnames called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (!colnames) {
	fprintf (stderr, "mpf_ILLlib_colnames called with NULL colnames\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nstruct = qslp->nstruct;

    if (qslp->colnames == 0) {
	fprintf (stderr, "LP does not have colnames assigned\n");
	rval = 1;
	ILL_CLEANUP;
    }
    for (i = 0; i < nstruct; i++) {
	len = strlen (qslp->colnames[i]) + 1;
	ILL_SAFE_MALLOC (colnames[i], len, char);
	strcpy (colnames[i], qslp->colnames[i]);
	ccount++;
    }


CLEANUP:

    if (rval) {
	for (i = 0; i < ccount; i++) {
	    ILL_IFFREE (colnames[i], char);
	}
    }
    ILL_RETURN (rval, "mpf_ILLlib_colnames");
}

static int mpf_reset_colindex (mpf_lpinfo * lp)
{
    int rval = 0;
    int test;
    mpf_ILLlpdata *qslp = lp->O;

    test = ILLsymboltab_index_ok (&qslp->coltab);
    if (!test) {
	rval = ILLsymboltab_index_reset (&qslp->coltab, qslp->nstruct,
	    qslp->colnames);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_RETURN (rval, "mpf_reset_colindex");
}

static int mpf_reset_rowindex (mpf_lpinfo * lp)
{
    int rval = 0;
    int test;
    mpf_ILLlpdata *qslp = lp->O;

    test = ILLsymboltab_index_ok (&qslp->rowtab);
    if (!test) {
	rval = ILLsymboltab_index_reset (&qslp->rowtab, qslp->nrows,
	    qslp->rownames);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_RETURN (rval, "mpf_reset_rowindex");
}

int mpf_ILLlib_colindex (mpf_lpinfo * lp,
      const char *name,
      int *colindex)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;

    *colindex = -1;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_colindex called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;

    rval = mpf_reset_colindex (lp);
    ILL_CLEANUP_IF (rval);

    rval = ILLsymboltab_getindex (&qslp->coltab, name, colindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_colindex");
}

int mpf_ILLlib_rowindex (mpf_lpinfo * lp,
      const char *name,
      int *rowindex)
{
    int rval = 0;
    mpf_ILLlpdata *qslp;

    *rowindex = -1;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_rowindex called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;

    rval = mpf_reset_rowindex (lp);
    ILL_CLEANUP_IF (rval);

    rval = ILLsymboltab_getindex (&qslp->rowtab, name, rowindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_rowindex");
}

int mpf_ILLlib_getbasis (mpf_lpinfo * lp,
      char *cstat,
      char *rstat)
{
    int rval = 0;
    int i, j, nrows;
    mpf_ILLlpdata *qslp;

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_getbasis called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (lp->basisid == -1) {
	fprintf (stderr, "mpf_ILLlib_getbasis called with modifed LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    nrows = lp->nrows;
    qslp = lp->O;

    for (i = 0; i < qslp->nstruct; i++) {
	j = qslp->structmap[i];
	switch (lp->vstat[j]) {
	case STAT_BASIC:
	    cstat[i] = QS_COL_BSTAT_BASIC;
	    break;
	case STAT_LOWER:
	    cstat[i] = QS_COL_BSTAT_LOWER;
	    break;
	case STAT_UPPER:
	    cstat[i] = QS_COL_BSTAT_UPPER;
	    break;
	case STAT_ZERO:
	    cstat[i] = QS_COL_BSTAT_FREE;
	    break;
	default:
	    fprintf (stderr, "unknown vstat in mpf_ILLlib_getbasis: %d\n", lp->vstat[j]);
	    rval = 1;
	    ILL_CLEANUP;
	}
    }

    for (i = 0; i < nrows; i++) {
	j = qslp->rowmap[i];
	if (qslp->rangeval && mpf_EGlpNumIsNeqqZero (qslp->rangeval[i])) {
	    switch (lp->vstat[j]) {
	    case STAT_BASIC:
		rstat[i] = QS_ROW_BSTAT_BASIC;
		break;
	    case STAT_LOWER:
		rstat[i] = QS_ROW_BSTAT_LOWER;
		break;
	    case STAT_UPPER:
		rstat[i] = QS_ROW_BSTAT_UPPER;
		break;
	    default:
		fprintf (stderr, "unknown vstat in mpf_ILLlib_getbasis 2\n");
		rval = 1;
		ILL_CLEANUP;
	    }
	} else {
	    switch (lp->vstat[j]) {
	    case STAT_BASIC:
		rstat[i] = QS_ROW_BSTAT_BASIC;
		break;
	    case STAT_UPPER:
	    case STAT_LOWER:
		rstat[i] = QS_ROW_BSTAT_LOWER;
		break;
	    default:
		fprintf (stderr, "unknown vstat in mpf_ILLlib_getbasis 3: %d, %d\n",
		    i, lp->vstat[j]);
		rval = 1;
		ILL_CLEANUP;
	    }
	}
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_getbasis");
}

int mpf_ILLlib_loadbasis (mpf_ILLlp_basis * B,
      int nstruct,
      int nrows,
      char *cstat,
      char *rstat)
{
    int i;
    int rval = 0;

    mpf_ILLlp_basis_init (B);

    if (!cstat || !rstat) {
	rval = 1;
	ILL_CLEANUP_IF (rval);
    }
    rval = mpf_ILLlp_basis_alloc (B, nstruct, nrows);
    ILL_CLEANUP_IF (rval);

    for (i = 0; i < nstruct; i++) {
	B->cstat[i] = cstat[i];
    }
    for (i = 0; i < nrows; i++) {
	B->rstat[i] = rstat[i];
    }

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_loadbasis");
}

#define mpf_READ_BASIS_XL 0
#define mpf_READ_BASIS_XU 1
#define mpf_READ_BASIS_LL 2
#define mpf_READ_BASIS_UL 3

int mpf_ILLlib_readbasis (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      const char *mpf_fname)
{
    int rval = 0;
    mpf_ILLlpdata *qslp = lp->O;
    int nstruct = qslp->nstruct;
    int nrows = qslp->nrows;
    int i, j, end = 0, sec, havename = 0;
    int rowtype, row, col;
    char *bname = 0;
    FILE *file_in = 0;
    mpf_ILLread_mps_state state;
    mpf_qsline_reader *in = NULL;

    mpf_ILLlp_basis_init (B);

    ILL_SAFE_MALLOC (B->cstat, qslp->nstruct, char);
    ILL_SAFE_MALLOC (B->rstat, qslp->nrows, char);
    B->nstruct = nstruct;
    B->nrows = nrows;

    for (j = 0; j < nstruct; j++) {
	B->cstat[j] = QS_COL_BSTAT_LOWER;
    }
    for (i = 0; i < nrows; i++) {
	B->rstat[i] = QS_ROW_BSTAT_BASIC;
    }

    file_in = fopen (mpf_fname, "r");
    if (file_in == 0) {
	fprintf (stderr, "unable to open %s for reading\n", mpf_fname);
	rval = 1;
	ILL_CLEANUP;
    }
    in = mpf_ILLline_reader_new ((mpf_qsread_line_fct) fgets, file_in);
    rval = mpf_ILLmps_state_init (&state, in, mpf_fname);
    ILL_CLEANUP_IF (rval);

    while (mpf_ILLmps_next_line (&state) == 0) {
	if (mpf_ILLmps_empty_key (&state)) {

	    /* Get the XL XU LL UL line */

	    if (!havename) {
		rval = mpf_ILLmps_error (&state, "BASIS data before NAME\n");
		ILL_CLEANUP;
	    }
	    if (!strcmp (state.field, "XL")) {
		rowtype = mpf_READ_BASIS_XL;
	    } else if (!strcmp (state.field, "XU")) {
		rowtype = mpf_READ_BASIS_XU;
	    } else if (!strcmp (state.field, "LL")) {
		rowtype = mpf_READ_BASIS_LL;
	    } else if (!strcmp (state.field, "UL")) {
		rowtype = mpf_READ_BASIS_UL;
	    } else {
		rval = mpf_ILLmps_error (&state, "BASIS \"%s\" is invalid\n", state.field);
		ILL_CLEANUP;
	    }

	    if (mpf_ILLmps_next_field (&state) == 0) {

		rval = mpf_ILLlib_colindex (lp, (const char *) state.field, &col);
		ILL_CLEANUP_IF (rval);
		if (col == -1) {
		    rval = mpf_ILLmps_error (&state, "BASIS col not in LP\n");
		    ILL_CLEANUP;
		}
		if (rowtype == mpf_READ_BASIS_XL || rowtype == mpf_READ_BASIS_XU) {
		    if (mpf_ILLmps_next_field (&state) == 0) {
			rval = mpf_ILLlib_rowindex (lp, (const char *) state.field, &row);
			ILL_CLEANUP_IF (rval);
			if (row == -1) {
			    rval = mpf_ILLmps_error (&state, "BASIS row not in LP\n");
			    ILL_CLEANUP;
			}
			if (rowtype == mpf_READ_BASIS_XL) {
			    B->cstat[col] = QS_COL_BSTAT_BASIC;
			    B->rstat[row] = QS_ROW_BSTAT_LOWER;

			} else {
			    B->cstat[col] = QS_COL_BSTAT_BASIC;
			    B->rstat[row] = QS_ROW_BSTAT_UPPER;
			}
		    } else {
			rval = mpf_ILLmps_error (&state, "BASIS line needs row and column\n");
			ILL_CLEANUP;
		    }
		} else {
		    if (rowtype == mpf_READ_BASIS_LL) {
			B->cstat[col] = QS_COL_BSTAT_LOWER;
		    } else {
			B->cstat[col] = QS_COL_BSTAT_UPPER;
		    }
		}
	    } else {
		rval = mpf_ILLmps_error (&state, "BASIS line has no row/column\n");
		ILL_CLEANUP;
	    }
	} else {
	    /* found a section indicator in col 1 */
	    if (!strcmp (state.key, mpf_ILLmps_section_name[ILL_MPS_ENDATA])) {
		end = 1;
		break;		/* done reading */
	    }
	    sec = mpf_ILLutil_index (mpf_ILLmps_section_name, state.key);
	    if (sec < 0 || sec != ILL_MPS_NAME) {
		rval = mpf_ILLmps_error (&state, "BASIS \"%s\" is not a key\n", state.key);
		ILL_CLEANUP;
	    }
	    if (havename) {
		rval = mpf_ILLmps_error (&state, "BASIS two name sections\n");
		ILL_CLEANUP;
	    }
	    havename = 1;

	    if (mpf_ILLmps_empty_field (&state)) {
		mpf_ILLmps_warn (&state, "BASIS blank NAME.");
	    } else {
		mpf_ILL_UTIL_STR (bname, state.field);
		printf ("Basis Name: %s\n", bname);
		fflush (stdout);
		if (strcmp (bname, qslp->probname)) {
		    mpf_ILLmps_warn (&state, "BASIS name does not match LP.");
		}
	    }
	}
    }

    if (!end) {
	mpf_ILLmps_warn (&state, "Missing ENDATA in basis file.");
    }
    if (!mpf_ILLmps_next_line (&state)) {
	mpf_ILLmps_warn (&state, "Ignoring text after ENDATA.");
    }
    if (!havename) {
	rval = mpf_ILLmps_error (&state, "BASIS no name section\n");
	ILL_CLEANUP;
    }
    /* Correct the free variables */

    for (j = 0; j < nstruct; j++) {
	col = lp->O->structmap[j];
	if (mpf_EGlpNumIsEqqual (qslp->lower[col], mpf_ILL_MINDOUBLE) &&
	    mpf_EGlpNumIsEqqual (qslp->upper[col], mpf_ILL_MAXDOUBLE) &&
	    B->cstat[j] == QS_COL_BSTAT_LOWER) {
	    B->cstat[j] = QS_COL_BSTAT_FREE;
	}
    }

CLEANUP:

    if (file_in)
	fclose (file_in);
    mpf_ILLline_reader_free (in);

    if (rval) {
	mpf_ILLlp_basis_free (B);
    }
    ILL_IFFREE (bname, char);

    ILL_RETURN (rval, "mpf_ILLlib_readbasis");
}

int mpf_ILLlib_writebasis (mpf_lpinfo * lp,
      mpf_ILLlp_basis * B,
      const char *mpf_fname)
{
    int rval = 0;
    FILE *out = 0;
    char *cstat = 0;
    char *rstat = 0;
    mpf_ILLlpdata *qslp;
    int i, j, nstruct, nrows;

    /* NOTE: non-basic free variables are encoded as non-basic at lower */

    if (!lp) {
	fprintf (stderr, "mpf_ILLlib_writebasis called without an LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    if (!B && lp->basisid == -1) {
	fprintf (stderr, "mpf_ILLlib_writebasis called with unsolved LP\n");
	rval = 1;
	ILL_CLEANUP;
    }
    qslp = lp->O;
    nstruct = qslp->nstruct;
    nrows = qslp->nrows;

    out = fopen (mpf_fname, "w");
    if (out == 0) {
	fprintf (stderr, "unable to open %s for writing\n", mpf_fname);
	rval = 1;
	ILL_CLEANUP;
    }
    if (B) {
	cstat = B->cstat;
	rstat = B->rstat;
    } else {
	ILL_SAFE_MALLOC (cstat, nstruct, char);
	ILL_SAFE_MALLOC (rstat, nrows, char);

	rval = mpf_ILLlib_getbasis (lp, cstat, rstat);
	ILL_CLEANUP_IF (rval);
    }

    fprintf (out, "NAME    %s\n", qslp->probname);

    /* Pick out the non-basic rows and find a matching basic column */

    i = 0;
    j = 0;
    do {
	while (i < nrows && rstat[i] == QS_ROW_BSTAT_BASIC) {
	    i++;
	}
	if (i < nrows) {
	    while (j < nstruct && cstat[j] != QS_COL_BSTAT_BASIC) {
		j++;
	    }
	    if (j == nstruct) {
		/* No basic column to match the non-basic row */
		fprintf (stderr, "No basic column to match non-basic row %d\n", i);
		rval = 1;
		goto CLEANUP;
	    }
	    if (rstat[i] == QS_ROW_BSTAT_LOWER) {
		fprintf (out, " XL %s %s\n", qslp->colnames[j], qslp->rownames[i]);
	    } else {
		fprintf (out, " XU %s %s\n", qslp->colnames[j], qslp->rownames[i]);
	    }
	    i++;
	    j++;
	}
    } while (i < nrows);

    /* Now go through and output the non-basic cols at upper bound */

    for (j = 0; j < nstruct; j++) {
	if (cstat[j] == QS_COL_BSTAT_UPPER) {
	    fprintf (out, " UL %s\n", qslp->colnames[j]);
	}
    }

    fprintf (out, "ENDATA\n");

CLEANUP:

    if (out)
	fclose (out);
    if (!B) {
	ILL_IFFREE (cstat, char);
	ILL_IFFREE (rstat, char);
    }
    ILL_RETURN (rval, "mpf_ILLlib_writebasis");
}

int mpf_ILLlib_getrownorms (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * rownorms)
{
    int rval = 0;
    int i, j, basic = 0;
    mpf_ILLlpdata *qslp = lp->O;
    int nstruct = lp->O->nstruct;
    int nrows = lp->O->nrows;

    mpf_check_pinf (pinf, &rval);
    if (rval) {
	/*
	        fprintf (stderr, "dual steepest mpf_edge norms not available\n");
	*/
	ILL_CLEANUP;
    }
    for (i = 0; i < nstruct; i++) {
	j = qslp->structmap[i];
	if (lp->vstat[j] == STAT_BASIC) {
	    mpf_EGlpNumCopy (rownorms[basic++], pinf->dsinfo.norms[lp->vindex[j]]);
	}
    }
    for (i = 0; i < nrows; i++) {
	j = qslp->rowmap[i];
	if (lp->vstat[j] == STAT_BASIC) {
	    mpf_EGlpNumCopy (rownorms[basic++], pinf->dsinfo.norms[lp->vindex[j]]);
	}
    }

    if (basic != nrows) {
	fprintf (stderr, "error in mpf_ILLlib_getrownorms\n");
	rval = 1;
	ILL_CLEANUP;
    }
CLEANUP:

    /*
        ILL_RETURN (rval, "mpf_ILLlib_getrownorms");
    */
    return rval;		/* Don't want error message */
}

int mpf_ILLlib_loadrownorms (mpf_lpinfo * lp,
      mpf_price_info * pinf,
      mpf_t * rownorms)
{
    int rval = 0;

    rval = mpf_ILLprice_load_rownorms (lp, rownorms, pinf);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_loadrownorms");
}

int mpf_ILLlib_recompute_rownorms (mpf_lpinfo * lp,
      mpf_price_info * pinf)
{
    int rval = 0;

    rval = mpf_ILLprice_build_pricing_info (lp, pinf, DUAL_PHASEII);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpf_ILLlib_recompute_rownorms");
}

int mpf_ILLlib_iter (mpf_lpinfo * lp)
{
    int iter = 0;

    if (lp && lp->cnts) {
	iter = lp->cnts->pI_iter + lp->cnts->pII_iter +
	    lp->cnts->dI_iter + lp->cnts->dII_iter;
    }
    return iter;
}

/* #define mpf_PRINT_TOL 0.000001 */
#define mpf_PRINT_TOL mpf_PFEAS_TOLER

int mpf_ILLlib_print_x (FILE * fd,
      mpf_lpinfo * lp,
      mpf_ILLlp_cache * C,
      mpf_t * x,
      int nonZerosOnly)
{
    int rval = 0;
    int j;
    mpf_ILLlpdata *qslp = lp->O;
    mpf_t *dx, *myx = 0;
    char *strtmp;

    /* If x is not specified, grab the LP solution */

    if (!x) {
	myx = mpf_EGlpNumAllocArray (lp->ncols);
	rval = mpf_ILLlib_get_x (lp, C, myx);
	ILL_CLEANUP_IF (rval);
	dx = myx;
    } else {
	dx = x;
    }

    fprintf (fd, "Solution Values\n");
    for (j = 0; j < qslp->nstruct; j++) {
	/* if (!nonZerosOnly || dx[j] > mpf_PRINT_TOL || dx[j] <
	   -mpf_PRINT_TOL) */
	if (!nonZerosOnly || mpf_EGlpNumIsNeqZero (dx[j], mpf_PRINT_TOL)) {
	    strtmp = mpf_EGlpNumGetStr (dx[j]);
	    ILL_FAILfalse (qslp->colnames[j] != NULL, "no NULL names PLEASE!");
	    fprintf (fd, "%s = %s\n", qslp->colnames[j], strtmp);
	    fflush (fd);
	    EGfree (strtmp);
	}
    }

CLEANUP:
    mpf_EGlpNumFreeArray (myx);
    ILL_RETURN (rval, "mpf_ILLlib_print_x");
}

int mpf_ILLlib_findName (mpf_ILLlpdata * qslp,
      int forRow,
      const char *name,
      int id,
      char buf[ILL_namebufsize])
{
    ILLsymboltab *tab;
    const char *mode;
    const char *p1, *p2;
    int sind, rval = 0;

    id++;
    tab = (forRow) ? &qslp->rowtab : &qslp->coltab;
    if (tab->tablesize == 0)
	ILLsymboltab_create (tab, 100);
    p1 = (forRow) ? "c" : "x";
    p2 = (forRow) ? "c_" : "x_";
    mode = (forRow) ? "row" : "column";
    if (name == 0) {
	ILLsymboltab_unique_name (tab, id, p1, buf);
	/*
	 * fprintf(stderr, "Generating %s name \"%s\".\n", mode, buf);
	 */
    } else {
	strcpy (buf, name);
    }
    if (!ILLsymboltab_lookup (tab, buf, &sind)) {
	rval = ILLsymboltab_uname (&qslp->rowtab, buf, p1, p2);
	if (name != NULL) {
	    fprintf (stderr, "Changing %s name \"%s\" to \"%s\".\n", mode, name, buf);
	}
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:
    ILL_RETURN (rval, "findName");
}

int mpf_ILLwrite_lp_file (mpf_ILLlpdata * lp,
      FILE * out,
      mpf_qserror_collector * c)
{
    int rval = 0;
    qsstring_reporter rep;
    ILLstring_reporter_copy (&rep, &lp->reporter);
    ILLstring_reporter_init (&lp->reporter, (qsreport_string_fct) fprintf, out);
    rval = mpf_ILLwrite_lp (lp, c);
    ILLstring_reporter_copy (&lp->reporter, &rep);
    return rval;
}

static void mpf_check_pinf (mpf_price_info * pinf,
      int *it_exists)
{
    if (!pinf || pinf->dI_price != QS_PRICE_DSTEEP ||
	pinf->dII_price != QS_PRICE_DSTEEP || pinf->dsinfo.norms == 0) {
	*it_exists = 1;
    } else {
	*it_exists = 0;
    }
}

#if 0
static int test_matrix (mpf_ILLmatrix * A)
{
    int rval = 0;
    int i, j, k;
    int ncols = A->matcols;
    int nrows = A->matrows;
    int matsize = A->matsize;
    int *mbeg = A->matbeg;
    int *mcnt = A->matcnt;
    int *mind = A->matind;
    int *tempi = 0;


    if (matsize == 0)
	ILL_CLEANUP;

    ILL_SAFE_MALLOC (tempi, matsize, int);
    for (i = 0; i < matsize; i++)
	tempi[i] = 0;

    for (i = 0; i < ncols; i++) {
	for (j = 0; j < mcnt[i]; j++) {
	    k = mind[mbeg[i] + j];
	    if (k < 0 || k >= nrows) {
		printf ("ERROR IN MATRIX: %d\n", k);
		printf ("ncols = %d, bad col = %d\n", ncols, i);
		printf ("bad cnt = %d, bad index = %d\n", mcnt[i], mbeg[i] + j);
		printf ("matcolsize = %d, matsize = %d\n", A->matcolsize, A->matsize);
		rval = 1;
		ILL_CLEANUP;
	    }
	    if (tempi[mbeg[i] + j] != 0) {
		printf ("ERROR: over written matrix\n");
		printf ("ncols = %d, bad col = %d\n", ncols, i);
		printf ("nrows = %d\n", nrows);
		printf ("bad cnt = %d, bad index = %d\n", mcnt[i], mbeg[i] + j);
		rval = 1;
		ILL_CLEANUP;
	    } else {
		tempi[mbeg[i] + j] = 1;
	    }
	}
    }

    for (i = A->matsize - A->matfree; i < A->matsize; i++) {
	if (tempi[i] != 0) {
	    printf ("ERROR: free space is being used\n");
	    rval = 1;
	    ILL_CLEANUP;
	}
    }

CLEANUP:

    ILL_IFFREE (tempi, int);
    ILL_RETURN (rval, "test_matrix");
}
#endif
