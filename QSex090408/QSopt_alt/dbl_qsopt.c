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

/* RCS_INFO = "$RCSfile: dbl_qsopt.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/*                                                                          */
/* User-level Functions                                                     */
/*                                                                          */
/* EXPORTED FUNCTIONS                                                       */
/*                                                                          */
/* int dbl_QSopt_primal (dbl_QSdata *p, int *status)                        */
/* int dbl_QSopt_dual (dbl_QSdata *p, int *status)                          */
/* dbl_QSdata *dbl_QScreate_prob (const char *name, int objsense)           */
/* dbl_QSdata *dbl_QSread_prob (const char *filename, const char *filetype) */
/* dbl_QSdata *dbl_QSload_prob (const char *probname, int ncols, int nrows, */
/*   int *cmatcnt, int *cmatbeg, int *cmatind, double *cmatval,             */
/*   int objsense, double *obj, double *rhs, char *sense,                   */
/*   double *lower, double *upper, const char **colnames,                   */
/*   const char **rownames)                                                 */
/* dbl_QSdata *dbl_QScopy_prob (dbl_QSdata *p, const char *newname)         */
/* int dbl_QSchange_objsense (dbl_QSdata *p, int newsense)                  */
/* int dbl_QSget_objsense (dbl_QSdata *p, int *objsense)                    */
/* int dbl_QSnew_col (dbl_QSdata *p, double obj, double lower, double upper,*/
/*   const char *name)                                                      */
/* int dbl_QSadd_cols (dbl_QSdata *p, int num, int *cmatcnt, int *cmatbeg,  */
/*   int *cmatind, double *cmatval, double *obj, double *lower,             */
/*   double *upper, const char **names)                                     */
/* int dbl_QSadd_col (dbl_QSdata *p, int cnt, int *cmatind, double *cmatval,     */
/* double obj, double lower, double upper, const char *name)         */
/* int dbl_QSnew_row (dbl_QSdata *p, double rhs, const char sense, char
   *name)   */
/* int dbl_QSadd_rows (dbl_QSdata *p, int num, int *rmatcnt, int *rmatbeg,       */
/* int *rmatind, double *rmatval, double *rhs, char *sense,          */
/* char **names)                                                     */
/* int dbl_QSadd_row (dbl_QSdata *p, int cnt, int *rmatind, double *rmatval,     */
/* double rhs, char sense, const char *name)                         */
/* int dbl_QSdelete_rows (dbl_QSdata *p, int num, int *dellist)                  */
/* int dbl_QSdelete_row (dbl_QSdata *p, int rowindex)                            */
/* int dbl_QSdelete_setrows (dbl_QSdata *p, int *flags)                          */
/* int dbl_QSdelete_cols (dbl_QSdata *p, int num, int *dellist)                  */
/* int dbl_QSdelete_col (dbl_QSdata *p, int colindex)                            */
/* int dbl_QSdelete_setcols (dbl_QSdata *p, int *flags)                          */
/* int dbl_QSdelete_named_column (dbl_QSdata *p, const char *colname)            */
/* int dbl_QSdelete_named_columns_list (dbl_QSdata *p, int num,                  */
/* const char **colnames)                                            */
/* int dbl_QSdelete_named_row (dbl_QSdata *p, const char *rowname)               */
/* int dbl_QSdelete_named_rows_list (dbl_QSdata *p, int num,                     */
/* const char **rownames)                                            */
/* int dbl_QSchange_senses (dbl_QSdata *p, int num, int *rowlist, char
   *sense)   */
/* int dbl_QSchange_sense (dbl_QSdata *p, int rowindex, char sense)              */
/* int dbl_QSchange_coef (dbl_QSdata *p, int rowindex, int colindex,             */
/* double coef)                                                      */
/* int dbl_QSchange_objcoef (dbl_QSdata *p, int indx, double coef)               */
/* int dbl_QSchange_rhscoef (dbl_QSdata *p, int indx, double coef)               */
/* int dbl_QSchange_bounds (dbl_QSdata *p, int num, int *collist, char *lu,      */
/* double *bounds)                                                   */
/* int dbl_QSchange_bound (dbl_QSdata *p, int indx, char lu, double bound)       */
/* int dbl_QSwrite_basis (dbl_QSdata *p, QSbasis *B, const char *filename)       */
/* QSbasis *dbl_QSget_basis (dbl_QSdata *p)                                      */
/* QSbasis *dbl_QSread_basis (dbl_QSdata *p, const char *filename)               */
/* int dbl_QSload_basis (dbl_QSdata *p, QSbasis *B)                              */
/* int dbl_QSread_and_load_basis (dbl_QSdata *p, const char *filename)           */
/* int dbl_QSload_basis_array (dbl_QSdata *p, char *cstat, char *rstat)          */
/* int dbl_QSload_basis_and_row_norms_array (dbl_QSdata *p, char *cstat,         */
/* char *rstat, double *rownorms)                                      */
/* int dbl_QSget_basis_array (dbl_QSdata *p, char *cstat, char *rstat)           */
/* int dbl_QSget_basis_and_row_norms_array (dbl_QSdata *p, char *cstat,          */
/* char *rstat, double *rownorms)                                      */
/* int dbl_QSget_binv_row (dbl_QSdata *p, int indx, double *binvrow)             */
/* int dbl_QSget_tableau_row (dbl_QSdata *p, int indx, double *tableaurow)       */
/* int dbl_QSget_basis_order (dbl_QSdata *p, int *basorder)                      */
/* int dbl_QSget_status (dbl_QSdata *p, int *status)                             */
/* int dbl_QSget_solution (dbl_QSdata *p, double *value, double *x,              */
/* double *pi, double *slack, double *rc),                           */
/* int dbl_QSget_objval (dbl_QSdata *p, double *value)                           */
/* int dbl_QSget_x_array (dbl_QSdata *p, double *x)                              */
/* int dbl_QSget_rc_array (dbl_QSdata *p, double *rc)                            */
/* int dbl_QSget_pi_array (dbl_QSdata *p, double *pi)                            */
/* int dbl_QSget_slack_array (dbl_QSdata *p, double *slack)                      */
/* int dbl_QSget_infeas_array (dbl_QSdata *p, double *pi)                        */
/* int dbl_QSget_named_x (dbl_QSdata *p, const char *colname, double *val)       */
/* int dbl_QSget_named_rc (dbl_QSdata *p, const char *colname, double *val)      */
/* int dbl_QSget_named_pi (dbl_QSdata *p, const char *rowname, double *val)      */
/* int dbl_QSget_named_slack (dbl_QSdata *p, const char *rowname, double
   *val)   */
/* int dbl_QSget_colcount (dbl_QSdata *p)                                        */
/* int dbl_QSget_rowcount (dbl_QSdata *p)                                        */
/* int dbl_QSget_nzcount (dbl_QSdata *p)                                         */
/* int dbl_QSget_obj (dbl_QSdata *p, double *obj),                               */
/* int dbl_QSget_rhs (dbl_QSdata *p, double *rhs)                                */
/* char* dbl_QSget_probname (dbl_QSdata *p)                                      */
/* char* dbl_QSget_objname (dbl_QSdata *p)                                       */
/* int dbl_QSget_columns (dbl_QSdata *p, int **colcnt, int **colbeg,             */
/* int **colind, double **colval, double **obj, double **lower,      */
/* double **upper, char ***names)                                    */
/* int dbl_QSget_columns_list (dbl_QSdata *p, int num, int *collist,             */
/* int **colcnt, int **colbeg, int **colind, double **colval,        */
/* double **obj, double **lower, double **upper, char ***names)      */
/* int dbl_QSget_rows (dbl_QSdata *p, int **rowcnt, int **rowbeg, int
   **rowind,  */
/* double **rowval, double **rhs, char **sense, char ***names)       */
/* int dbl_QSget_rows_list (dbl_QSdata *p, int num, int *rowlist, int
   **rowcnt,  */
/* int **rowbeg, int **rowind, double **rowval, double **rhs,        */
/* char **sense, char ***names)                                      */
/* int dbl_QSget_column_index (dbl_QSdata *p, const char *name, int
   *colindex)   */
/* int dbl_QSget_row_index (dbl_QSdata *p, const char *name, int *rowindex)      */
/* int dbl_QSget_rownames (dbl_QSdata *p, char **rownames)                       */
/* int dbl_QSget_colnames (dbl_QSdata *p, char **colnames)                       */
/* int dbl_QSget_bound (dbl_QSdata *p, int colindex, char lu, double *bound)     */
/* int dbl_QSget_bounds (dbl_QSdata *p, double *lower, double *upper)            */
/* int dbl_QSget_intcount (dbl_QSdata *p, int *count)                            */
/* int dbl_QSget_intflags (dbl_QSdata *p, int *intflags)                         */
/* int dbl_QScompute_row_norms (dbl_QSdata *p)                                   */
/* void dbl_QSfree_prob (dbl_QSdata *p)                                          */
/* void dbl_QSfree_basis (QSbasis *B)                                        */
/* int dbl_QSwrite_prob (dbl_QSdata *p, const char *filename,                    */
/* const char *filetype)                                             */
/* int dbl_QSwrite_prob_file (dbl_QSdata *p, FILE *file, const char
   *filetype)   */
/* int dbl_QSset_param (dbl_QSdata *p, int whichparam, int newvalue)             */
/* int QSset_param_double (dbl_QSdata *p, int whichparam, double newvalue)   */
/* int dbl_QSget_param (dbl_QSdata *p, int whichparam, int *value)               */
/* int QSget_param_double (dbl_QSdata *p, int whichparam, double *value)     */
/* int dbl_QStest_row_norms (dbl_QSdata *p)                                      */
/* int dbl_QSopt_strongbranch (dbl_QSdata *p, int ncand, int *candidatelist,     */
/* double *xlist, double *down_vals, double *up_vals,                */
/* int iterations, double objbound)                                  */
/* int dbl_QSopt_pivotin_row (dbl_QSdata *p, int rcnt, int *rlist)               */
/* int dbl_QSopt_pivotin_col (dbl_QSdata *p, int ccnt, int *clist)               */
/* void dbl_QSfree (void *ptr)                                               */
/* void dbl_QSstart (void)                                                   */
/* void dbl_QSend (void)                                                     */
/* char *dbl_QSversion (void))                                               */
/* */
/* dbl_NEW FUNCTIONS - Add to Docs                                           */
/* */
/* char *dbl_QSversion (void))                                               */
/* int dbl_QSget_objsense (dbl_QSdata *p, int *objsense)                         */
/* */
/****************************************************************************/

#include "econfig.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdata.h"
#include "dbl_lpdefs.h"
#include "dbl_simplex.h"
#include "dbl_price.h"
#include "dbl_qstruct.h"
#include "dbl_qsopt.h"
#include "dbl_lib.h"
#include "dbl_mps.h"
#include "dbl_lp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
void dbl_QSset_precision (const unsigned prec)
{
    EGlpNumSetPrecision (prec);
    dbl_ILLchange_precision ();
    /* change the numbers */
}

static void dbl_init_basis (QSbasis * B),
  dbl_free_cache (dbl_QSdata * p);

int dbl_grab_cache (dbl_QSdata * p,
      int status);
static int dbl_opt_work (dbl_QSdata * p,
      int *status,
      int primal_or_dual),
  dbl_qsbasis_to_illbasis (QSbasis * qB,
      dbl_ILLlp_basis * B),
  dbl_illbasis_to_qsbasis (dbl_ILLlp_basis * B,
      QSbasis * qB),
  dbl_grab_basis (dbl_QSdata * p),
  dbl_check_qsdata_pointer (dbl_QSdata * p);


dbl_QSLIB_INTERFACE int dbl_QSopt_primal (dbl_QSdata * p,
      int *status)
{
    int rval = 0;

    if (status)
	*status = QS_LP_UNSOLVED;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    /* If both the basis and the cache exist, then skip the optimization */

    if (!p->basis || !p->cache) {
	rval = dbl_opt_work (p, status, 0);
	ILL_CLEANUP_IF (rval);
    } else {
	if (status)
	    *status = p->cache->status;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSopt_primal");
}

dbl_QSLIB_INTERFACE int dbl_QSopt_dual (dbl_QSdata * p,
      int *status)
{
    int rval = 0;

    if (status)
	*status = QS_LP_UNSOLVED;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (!p->basis || !p->cache || !p->factorok) {
	rval = dbl_opt_work (p, status, 1);
	ILL_CLEANUP_IF (rval);
    } else {
	if (status)
	    *status = p->cache->status;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSopt_dual");
}

static int dbl_opt_work (dbl_QSdata * p,
      int *status,
      int primal_or_dual)
{
    int rval = 0;
    int rstatus = QS_LP_UNSOLVED;
    dbl_QSdata *p2 = 0;

    if (p->basis) {
	if (p->basis->nstruct != p->qslp->nstruct ||
	    p->basis->nrows != p->qslp->nrows) {
	    fprintf (stderr, "Size of basis does not match LP\n");
	    rval = 1;
	    goto CLEANUP;
	}
    }
    if (!p->basis && p->lp->basisid == -1 && p->simplex_scaling == 1) {
	/* Try scaling by copying the LP and solving */

	dbl_ILLprice_free_pricing_info (p->pricing);	/* Just to be sure  */
	p->factorok = 0;	/* that p is clean. */

	p2 = dbl_QScopy_prob (p, "scaled_lp");
	if (p2 == 0)
	    goto CLEANUP;

	rval = dbl_ILLlp_scale (p2->qslp);
	ILL_CLEANUP_IF (rval);

	if (primal_or_dual == 0) {
	    rval = dbl_ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
		PRIMAL_SIMPLEX, 0, p2->simplex_display);
	} else {
	    rval = dbl_ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
		DUAL_SIMPLEX, 0, p2->simplex_display);
	}
	ILL_CLEANUP_IF (rval);

	rval = dbl_grab_basis (p2);
	ILL_CLEANUP_IF (rval);

	if (p->basis) {
	    dbl_ILLlp_basis_free (p->basis);
	    ILL_IFFREE (p->basis, dbl_ILLlp_basis);
	}
	p->basis = p2->basis;
	p2->basis = 0;
	dbl_QSfree_prob (p2);
	p2 = 0;
    }
    if (primal_or_dual == 0) {
	if (p->factorok == 0) {
	    if (p->basis == 0)
		p->lp->basisid = -1;
	    rval = dbl_ILLlib_optimize (p->lp, p->basis, p->pricing, PRIMAL_SIMPLEX,
		&rstatus, p->simplex_display);
	} else {
	    dbl_ILLprice_free_pricing_info (p->pricing);
	    if (p->lp->basisid != -1)
		p->lp->fbasisid = p->lp->basisid;
	    rval = dbl_ILLlib_optimize (p->lp, 0, p->pricing,
		PRIMAL_SIMPLEX, &rstatus, p->simplex_display);
	}
    } else {
	if (p->factorok == 0) {
	    if (p->basis == 0)
		p->lp->basisid = -1;
	    rval = dbl_ILLlib_optimize (p->lp, p->basis, p->pricing, DUAL_SIMPLEX,
		&rstatus, p->simplex_display);
	} else {
	    /* The factorization and rownorms should be up-to-date */
	    if (p->lp->basisid != -1) {
		p->lp->fbasisid = p->lp->basisid;
	    } else {
		dbl_ILLprice_free_pricing_info (p->pricing);
	    }
	    rval = dbl_ILLlib_optimize (p->lp, 0, p->pricing,
		DUAL_SIMPLEX, &rstatus, p->simplex_display);
	}
    }
    ILL_CLEANUP_IF (rval);

    rval = dbl_grab_basis (p);
    ILL_CLEANUP_IF (rval);

    if (rstatus == QS_LP_OPTIMAL) {
	rval = dbl_grab_cache (p, rstatus);
	ILL_CLEANUP_IF (rval);
    } else {
	dbl_free_cache (p);
    }

    p->factorok = 1;

#if 0
    p->lp->basisid = -1;	/* This will cause the basis to be reloaded
				   at the */
    /* next optimization - it could be moved into the  */
    /* add/del routines if we want to cache the        */
    /* factored basis.                                 */
    -switched to having qs_simplex load a basis whenever it is passed;
    the trouble with keeping basisid == -1 is that the QS - level routines
    will not know if the current lp has been optimized (so, for example,
	getbasis will not work)
	.
#endif
	CLEANUP:p->qstatus = rstatus;
    if (status)
	*status = rstatus;

    if (p2)
	dbl_QSfree_prob (p2);
    ILL_RETURN (rval, "dbl_opt_work");
}

dbl_QSLIB_INTERFACE int dbl_QSopt_pivotin_row (dbl_QSdata * p,
      int rcnt,
      int *rlist)
{
    int basismod = 0;
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing == 0) {
	ILL_ERROR (rval, "pricing info not available in dbl_QSopt_pivotin_row\n");
    }
    rval = dbl_ILLsimplex_pivotin (p->lp, p->pricing, rcnt, rlist,
	SIMPLEX_PIVOTINROW, &basismod);
    ILL_CLEANUP_IF (rval);

    rval = dbl_grab_basis (p);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSopt_pivotin_row");
}

dbl_QSLIB_INTERFACE int dbl_QSopt_pivotin_col (dbl_QSdata * p,
      int ccnt,
      int *clist)
{
    int basismod = 0;
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing == 0) {
	ILL_ERROR (rval, "pricing info not available in QSopt_pivotin\n");
    }
    rval = dbl_ILLsimplex_pivotin (p->lp, p->pricing, ccnt, clist,
	SIMPLEX_PIVOTINCOL, &basismod);
    ILL_CLEANUP_IF (rval);

    rval = dbl_grab_basis (p);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSopt_pivotin_col");
}


dbl_QSLIB_INTERFACE int dbl_QSopt_strongbranch (dbl_QSdata * p,
      int ncand,
      int *candidatelist,
      double *xlist,
      double *down_vals,
      double *up_vals,
      int iterations,
      double objbound)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing == 0) {
	rval = 1;
	ILL_CLEANUP_IF (rval);
    }
    rval = dbl_ILLlib_strongbranch (p->lp, p->pricing, candidatelist, ncand,
	xlist, down_vals, up_vals, iterations, objbound);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;
    dbl_free_cache (p);
    p->qstatus = QS_LP_UNSOLVED;/* Was set to MODIFIED in dbl_free_cache () */

CLEANUP:

    ILL_RETURN (rval, "dbl_QSopt_strongbranch");
}

dbl_QSLIB_INTERFACE dbl_QSdata *dbl_QScreate_prob (const char *name,
      int objsense)
{
    int rval = 0;
    dbl_QSdata *p = 0;
    int len;

    ILL_SAFE_MALLOC (p, 1, dbl_QSdata);
    if (!p) {
	fprintf (stderr, "out of memory in dbl_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    p->qslp = 0;
    p->lp = 0;
    p->pricing = 0;
    p->basis = 0;
    p->cache = 0;
    p->qstatus = QS_LP_UNSOLVED;
    p->factorok = 0;

    p->simplex_display = 0;
    p->simplex_scaling = 1;

    ILL_SAFE_MALLOC (p->qslp, 1, dbl_ILLlpdata);
    if (!p->qslp) {
	fprintf (stderr, "out of memory in dbl_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    dbl_ILLlpdata_init (p->qslp);

    ILL_SAFE_MALLOC (p->lp, 1, dbl_lpinfo);
    if (!p->lp) {
	fprintf (stderr, "out of memory in dbl_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    dbl_EGlpNumInitVar (p->lp->objval);
    dbl_EGlpNumInitVar (p->lp->pobjval);
    dbl_EGlpNumInitVar (p->lp->dobjval);
    dbl_EGlpNumInitVar (p->lp->pinfeas);
    dbl_EGlpNumInitVar (p->lp->dinfeas);
    dbl_EGlpNumInitVar (p->lp->objbound);
    dbl_EGlpNumInitVar (p->lp->upd.piv);
    dbl_EGlpNumInitVar (p->lp->upd.dty);
    dbl_EGlpNumInitVar (p->lp->upd.c_obj);
    dbl_EGlpNumInitVar (p->lp->upd.tz);
    dbl_ILLsimplex_init_lpinfo (p->lp);
    dbl_ILLsimplex_load_lpinfo (p->qslp, p->lp);

    ILL_SAFE_MALLOC (p->pricing, 1, dbl_price_info);
    if (!p->pricing) {
	fprintf (stderr, "out of memory in dbl_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    dbl_EGlpNumInitVar (p->pricing->htrigger);
    dbl_ILLprice_init_pricing_info (p->pricing);
    p->pricing->pI_price = QS_DEFAULT_PRICE_PI;
    p->pricing->pII_price = QS_DEFAULT_PRICE_PII;
    p->pricing->dI_price = QS_DEFAULT_PRICE_DI;
    p->pricing->dII_price = QS_DEFAULT_PRICE_DII;

    if (name) {
	len = strlen (name) + 1;
	ILL_SAFE_MALLOC (p->name, len, char);
	strcpy (p->name, name);
    } else {
	ILL_SAFE_MALLOC (p->name, 7, char);
	sprintf (p->name, "noname");
    }

    len = strlen (p->name) + 1;
    ILL_SAFE_MALLOC (p->qslp->probname, len, char);
    strcpy (p->qslp->probname, p->name);

    if (objsense == QS_MAX) {
	p->qslp->objsense = QS_MAX;
    }
CLEANUP:

    if (rval) {
	dbl_QSfree_prob (p);
	p = 0;
    }
    return p;
}

dbl_QSLIB_INTERFACE dbl_QSdata *dbl_QSread_prob (const char *filename,
      const char *filetype)
{
    dbl_QSdata *p = 0;
    FILE *file = 0;
    dbl_QSline_reader reader;

    if ((file = fopen (filename, "r")) == 0) {
	perror (filename);
	fprintf (stderr, "Unable to open \"%s\" for input.\n", filename);
    }
    if (file == NULL)
	goto CLEANUP;

    reader = dbl_ILLline_reader_new ((dbl_qsread_line_fct) fgets, file);
    p = dbl_QSget_prob (reader, filename, filetype);
    dbl_QSline_reader_free (reader);	/* Bico - 040723 */

CLEANUP:
    if (file != NULL) {
	fclose (file);
    }
    return p;
}

dbl_QSLIB_INTERFACE dbl_QSdata *dbl_QSload_prob (const char *probname,
      int ncols,
      int nrows,
      int *cmatcnt,
      int *cmatbeg,
      int *cmatind,
      double *cmatval,
      int objsense,
      double *obj,
      double *rhs,
      char *sense,
      double *lower,
      double *upper,
      const char **colnames,
      const char **rownames)
{
    int rval = 0;
    dbl_QSdata *p = 0;

    p = dbl_QScreate_prob (probname, objsense);
    if (p == 0)
	goto CLEANUP;

    rval = dbl_ILLlib_newrows (p->lp, 0, nrows, rhs, sense, 0, rownames);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_addcols (p->lp, 0, ncols, cmatcnt, cmatbeg, cmatind,
	cmatval, obj, lower, upper, colnames, 0);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;

CLEANUP:

    if (rval) {
	dbl_QSfree_prob (p);
	p = 0;
    }
    return p;
}

dbl_QSLIB_INTERFACE dbl_QSdata *dbl_QScopy_prob (dbl_QSdata * p,
      const char *newname)
{
    int rval = 0;
    int j, col, beg, pindex, hit;
    dbl_QSdata *p2 = 0;
    char *coln;
    char buf[ILL_namebufsize];

    /* printf ("dbl_QScopy_prob ...\n"); fflush (stdout); */

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    p2 = dbl_QScreate_prob (newname, p->qslp->objsense);
    if (p2 == 0)
	goto CLEANUP;

    rval = dbl_ILLlib_newrows (p2->lp, 0, p->qslp->nrows,
	p->qslp->rhs, p->qslp->sense, p->qslp->rangeval,
	(const char **) p->qslp->rownames);
    ILL_CLEANUP_IF (rval);

    for (j = 0; j < p->qslp->nstruct; j++) {
	col = p->qslp->structmap[j];
	if (p->qslp->colnames)
	    coln = p->qslp->colnames[j];
	else
	    coln = 0;
	beg = p->qslp->A.matbeg[col];

	/* Monika: Note that Java will need to handle these arrays */
	/* without using the beg offset.  The easiest way  */
	/* may be to copy the arrays, as in the getcols()  */
	/* code in dbl_lib.c.                                  */

	rval = dbl_ILLlib_addcol (p2->lp, 0,
	    p->qslp->A.matcnt[col],
	    p->qslp->A.matind + beg, p->qslp->A.matval + beg,
	    p->qslp->obj[col], p->qslp->lower[col],
	    p->qslp->upper[col], coln, 0);
	ILL_CLEANUP_IF (rval);
    }

    p2->qslp->objsense = p->qslp->objsense;

    p2->factorok = 0;
    p2->simplex_display = p->simplex_display;
    p2->simplex_scaling = p->simplex_scaling;
    dbl_EGlpNumClearVar (p2->pricing->htrigger);
    *(p2->pricing) = *(p->pricing);
    /* I added this line because copying the dbl_heap (as a pointer) doesn't
       make any sense ! */
    dbl_ILLheap_init (&(p2->pricing->h));
    dbl_EGlpNumInitVar (p2->pricing->htrigger);
    dbl_EGlpNumCopy (p2->pricing->htrigger, p->pricing->htrigger);

    if (p->qslp->intmarker != 0) {
	ILL_SAFE_MALLOC (p2->qslp->intmarker, p->qslp->nstruct, char);
	for (j = 0; j < p->qslp->nstruct; j++) {
	    p2->qslp->intmarker[j] = p->qslp->intmarker[j];
	}
    }
    if (p->qslp->objname != 0) {
	dbl_ILL_UTIL_STR (p2->qslp->objname, p->qslp->objname);
    } else {
	strcpy (buf, "obj");
	rval = ILLsymboltab_uname (&p2->qslp->rowtab, buf, "", NULL);
	ILL_CLEANUP_IF (rval);
	dbl_ILL_UTIL_STR (p2->qslp->objname, buf);
    }
    rval = ILLsymboltab_register (&p2->qslp->rowtab, p2->qslp->objname,
	-1, &pindex, &hit);
    rval = rval || hit;
    ILL_CLEANUP_IF (rval);

    ILLstring_reporter_copy (&p2->qslp->reporter, &p->qslp->reporter);

CLEANUP:

    if (rval) {
	dbl_QSfree_prob (p2);
	p2 = 0;
    }
    return p2;
}

dbl_QSLIB_INTERFACE int dbl_QSchange_objsense (dbl_QSdata * p,
      int newsense)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (newsense != QS_MIN && newsense != QS_MAX) {
	fprintf (stderr, "Illegal objective sense %d\n", newsense);
	rval = 1;
	goto CLEANUP;
    }
    if (p->qslp->objsense != newsense) {
	p->qslp->objsense = newsense;
	dbl_free_cache (p);
    }
CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_objsense");
}

dbl_QSLIB_INTERFACE int dbl_QSget_objsense (dbl_QSdata * p,
      int *objsense)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (objsense)
	*objsense = p->qslp->objsense;

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_objsense");
}


dbl_QSLIB_INTERFACE int dbl_QSnew_col (dbl_QSdata * p,
      double obj,
      double lower,
      double upper,
      const char *name)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_newcol (p->lp, p->basis, obj, lower, upper, name, p->factorok);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSnew_col");
}

dbl_QSLIB_INTERFACE int dbl_QSadd_cols (dbl_QSdata * p,
      int num,
      int *cmatcnt,
      int *cmatbeg,
      int *cmatind,
      double *cmatval,
      double *obj,
      double *lower,
      double *upper,
      const char **names)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_addcols (p->lp, p->basis, num, cmatcnt, cmatbeg,
	cmatind, cmatval, obj, lower, upper, names,
	p->factorok);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSadd_cols");
}

dbl_QSLIB_INTERFACE int dbl_QSadd_col (dbl_QSdata * p,
      int cnt,
      int *cmatind,
      double *cmatval,
      double obj,
      double lower,
      double upper,
      const char *name)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_addcol (p->lp, p->basis, cnt, cmatind, cmatval,
	obj, lower, upper, name, p->factorok);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSadd_col");
}

dbl_QSLIB_INTERFACE int dbl_QSnew_row (dbl_QSdata * p,
      double rhs,
      int sense,
      const char *name)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_newrow (p->lp, p->basis, rhs, sense, dbl_zeroLpNum, name);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;
    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSnew_row");
}

dbl_QSLIB_INTERFACE int dbl_QSadd_rows (dbl_QSdata * p,
      int num,
      int *rmatcnt,
      int *rmatbeg,
      int *rmatind,
      double *rmatval,
      double *rhs,
      char *sense,
      const char **names)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_addrows (p->lp, p->basis, num, rmatcnt, rmatbeg,
	rmatind, rmatval, rhs, sense, 0,
	names, &(p->factorok));
    ILL_CLEANUP_IF (rval);

    if (p->factorok == 1 && p->basis->rownorms) {
	rval = dbl_ILLlib_loadrownorms (p->lp, p->pricing, p->basis->rownorms);
	ILL_CLEANUP_IF (rval);
	/* This really should go inside of dbl_ILLlib_addrows, once pinf is  */
	/* is moved into the lp struct.                                  */
    }
    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSadd_rows");
}

dbl_QSLIB_INTERFACE int dbl_QSadd_row (dbl_QSdata * p,
      int cnt,
      int *rmatind,
      double *rmatval,
      double *rhs,
      int sense,
      const char *name)
{
    int rval = 0;
    int vmatcnt[1];
    int vmatbeg[1];
    char vsense[1];
    const char *vnames[1];

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vmatcnt[0] = cnt;
    vmatbeg[0] = 0;
    vsense[0] = sense;
    vnames[0] = name;

    rval = dbl_QSadd_rows (p, 1, vmatcnt, vmatbeg, rmatind, rmatval, rhs, vsense,
	vnames);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSadd_row");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_rows (dbl_QSdata * p,
      int num,
      int *dellist)
{
    int rval = 0;
    int basis_ok = 0;
    int cache_ok = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_delrows (p->lp, p->basis, p->cache, num, dellist, &basis_ok,
	&cache_ok);
    ILL_CLEANUP_IF (rval);

    /* For now, just remove the basis - wait for pivotin */

    if (p->basis && !basis_ok) {
	dbl_ILLlp_basis_free (p->basis);
	ILL_IFFREE (p->basis, dbl_ILLlp_basis);
    }
    p->factorok = 0;

    if (!p->basis || !basis_ok || !cache_ok) {
	/* Note: If we only delete basic rows then cached soln is valid */
	dbl_free_cache (p);
    }
CLEANUP:

    ILL_RETURN (rval, "dbl_QSdelete_rows");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_row (dbl_QSdata * p,
      int rowindex)
{
    int rval = 0;
    int vdellist[1];

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = rowindex;

    rval = dbl_QSdelete_rows (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSdelete_row");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_setrows (dbl_QSdata * p,
      int *flags)
{
    int rval = 0;
    int j, num = 0;
    int *dellist = 0;
    int nrows;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = p->qslp->nrows;

    for (j = 0; j < nrows; j++) {
	if (flags[j] == 1)
	    num++;
    }

    if (num > 0) {
	ILL_SAFE_MALLOC (dellist, num, int);

	for (j = 0, num = 0; j < nrows; j++) {
	    if (flags[j] == 1) {
		dellist[num++] = j;
	    }
	}

	rval = dbl_QSdelete_rows (p, num, dellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (dellist, int);
    ILL_RETURN (rval, "dbl_QSdelete_setrows");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_named_row (dbl_QSdata * p,
      const char *rowname)
{
    int rval = 0;
    int i, vdellist[1];

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_QSget_row_index (p, rowname, &i);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = i;

    rval = dbl_QSdelete_rows (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSdelete_named_row");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_named_rows_list (dbl_QSdata * p,
      int num,
      const char **rownames)
{
    int rval = 0;
    int i, k;
    int *vdellist = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (num > 0) {
	ILL_SAFE_MALLOC (vdellist, num, int);
	for (k = 0; k < num; k++) {
	    rval = dbl_QSget_row_index (p, rownames[k], &i);
	    ILL_CLEANUP_IF (rval);
	    vdellist[k] = i;
	}

	rval = dbl_QSdelete_rows (p, num, vdellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (vdellist, int);
    ILL_RETURN (rval, "dbl_QSdelete_named_rows_list");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_cols (dbl_QSdata * p,
      int num,
      int *dellist)
{
    int rval = 0;
    int basis_ok;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_delcols (p->lp, p->basis, num, dellist, &basis_ok);
    ILL_CLEANUP_IF (rval);

    /* For now, just remove the basis - wait for pivotout */

    if (p->basis && !basis_ok) {
	dbl_ILLlp_basis_free (p->basis);
	ILL_IFFREE (p->basis, dbl_ILLlp_basis);
    }
    p->factorok = 0;
    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSdelete_cols");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_col (dbl_QSdata * p,
      int colindex)
{
    int rval = 0;
    int vdellist[1];

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = colindex;

    rval = dbl_QSdelete_cols (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSdelete_col");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_setcols (dbl_QSdata * p,
      int *flags)
{
    int rval = 0;
    int j, num = 0;
    int *dellist = 0;
    int ncols;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = p->qslp->nstruct;

    for (j = 0; j < ncols; j++) {
	if (flags[j] == 1)
	    num++;
    }

    if (num > 0) {
	ILL_SAFE_MALLOC (dellist, num, int);

	for (j = 0, num = 0; j < ncols; j++) {
	    if (flags[j] == 1) {
		dellist[num++] = j;
	    }
	}

	rval = dbl_QSdelete_cols (p, num, dellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (dellist, int);
    ILL_RETURN (rval, "dbl_QSdelete_setcols");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_named_column (dbl_QSdata * p,
      const char *colname)
{
    int rval = 0;
    int j, vdellist[1];

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_QSget_column_index (p, colname, &j);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = j;

    rval = dbl_QSdelete_cols (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSdelete_named_column");
}

dbl_QSLIB_INTERFACE int dbl_QSdelete_named_columns_list (dbl_QSdata * p,
      int num,
      const char **colnames)
{
    int rval = 0;
    int i, j;
    int *vdellist = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (num > 0) {
	ILL_SAFE_MALLOC (vdellist, num, int);
	for (i = 0; i < num; i++) {
	    rval = dbl_QSget_column_index (p, colnames[i], &j);
	    ILL_CLEANUP_IF (rval);
	    vdellist[i] = j;
	}

	rval = dbl_QSdelete_cols (p, num, vdellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (vdellist, int);
    ILL_RETURN (rval, "dbl_QSdelete_named_columns_list");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_senses (dbl_QSdata * p,
      int num,
      int *rowlist,
      char *sense)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_chgsense (p->lp, num, rowlist, sense);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_senses");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_sense (dbl_QSdata * p,
      int rowindex,
      int sense)
{
    int rval = 0;
    int vrowlist[1];
    char vsenselist[1];

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vrowlist[0] = rowindex;
    vsenselist[0] = sense;

    rval = dbl_QSchange_senses (p, 1, vrowlist, vsenselist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_sense");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_coef (dbl_QSdata * p,
      int rowindex,
      int colindex,
      double coef)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_chgcoef (p->lp, rowindex, colindex, coef);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_coef");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_objcoef (dbl_QSdata * p,
      int indx,
      double coef)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_chgobj (p->lp, indx, coef);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_objcoef");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_rhscoef (dbl_QSdata * p,
      int indx,
      double coef)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_chgrhs (p->lp, indx, coef);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_rhscoef");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_bounds (dbl_QSdata * p,
      int num,
      int *collist,
      char *lu,
      double *bounds)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_chgbnds (p->lp, num, collist, lu, bounds);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_bounds");
}

dbl_QSLIB_INTERFACE int dbl_QSchange_bound (dbl_QSdata * p,
      int indx,
      int lu,
      double bound)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_chgbnd (p->lp, indx, lu, bound);
    ILL_CLEANUP_IF (rval);

    dbl_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSchange_bound");
}

#if 0
/*
 * Bico - I removed this on 02.04.22.  I don't think we need to support
 * this type of interface (the loading via arrays can do the job)
 */
dbl_QSLIB_INTERFACE QSbasis *QScreate_basis (int nstruct,
      int nrows)
{
    int rval = 0;
    int i;
    QSbasis *B = 0;

    ILL_SAFE_MALLOC (B, 1, QSbasis);

    B->nstruct = nstruct;
    B->nrows = nrows;
    B->cstat = 0;
    B->rstat = 0;

    if (nstruct) {
	ILL_SAFE_MALLOC (B->cstat, nstruct, char);
    }
    if (nrows) {
	ILL_SAFE_MALLOC (B->rstat, nrows, char);
    }
    for (i = 0; i < nstruct; i++)
	B->cstat[i] = 0;
    for (i = 0; i < nrows; i++)
	B->rstat[i] = 0;

CLEANUP:

    if (rval) {
	dbl_QSfree_basis (B);
	B = 0;
    }
    return B;
}
#endif

dbl_QSLIB_INTERFACE QSbasis *dbl_QSread_basis (dbl_QSdata * p,
      const char *filename)
{
    int rval = 0;
    QSbasis *qB = 0;
    dbl_ILLlp_basis B;

    dbl_ILLlp_basis_init (&B);

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ILL_NEW (qB, QSbasis);
    dbl_init_basis (qB);

    rval = dbl_ILLlib_readbasis (p->lp, &B, filename);
    ILL_CLEANUP_IF (rval);

    rval = dbl_illbasis_to_qsbasis (&B, qB);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    if (rval && qB) {
	dbl_QSfree_basis (qB);
	qB = 0;
    }
    dbl_ILLlp_basis_free (&B);

    return qB;
}

dbl_QSLIB_INTERFACE int dbl_QSload_basis (dbl_QSdata * p,
      QSbasis * B)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (B->nstruct != p->qslp->nstruct || B->nrows != p->qslp->nrows) {
	fprintf (stderr, "size of basis does not match lp\n");
	rval = 1;
	goto CLEANUP;
    }
    if (p->basis == 0) {
	ILL_SAFE_MALLOC (p->basis, 1, dbl_ILLlp_basis);
	dbl_ILLlp_basis_init (p->basis);
    } else {
	dbl_ILLlp_basis_free (p->basis);
    }

    rval = dbl_qsbasis_to_illbasis (B, p->basis);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;

CLEANUP:

    ILL_RETURN (rval, "dbl_QSload_basis");
}

dbl_QSLIB_INTERFACE int dbl_QSread_and_load_basis (dbl_QSdata * p,
      const char *filename)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->basis == 0) {
	ILL_SAFE_MALLOC (p->basis, 1, dbl_ILLlp_basis);
	dbl_ILLlp_basis_init (p->basis);
    } else {
	dbl_ILLlp_basis_free (p->basis);
    }

    rval = dbl_ILLlib_readbasis (p->lp, p->basis, filename);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    return rval;
}

dbl_QSLIB_INTERFACE int dbl_QSload_basis_array (dbl_QSdata * p,
      char *cstat,
      char *rstat)
{
    int rval = 0;
    int i;
    dbl_ILLlp_basis *B;
    dbl_ILLlpdata *qslp;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    qslp = p->qslp;

    if (qslp->nstruct > 0 && cstat == 0) {
	fprintf (stderr, "dbl_QSload_basis_array called without cstat\n");
	rval = 1;
	goto CLEANUP;
    }
    if (qslp->nrows > 0 && rstat == 0) {
	fprintf (stderr, "dbl_QSload_basis_array called without rstat\n");
	rval = 1;
	goto CLEANUP;
    }
    if (p->basis == 0) {
	ILL_SAFE_MALLOC (p->basis, 1, dbl_ILLlp_basis);
	dbl_ILLlp_basis_init (p->basis);
    } else {
	dbl_ILLlp_basis_free (p->basis);
    }

    B = p->basis;

    B->nstruct = qslp->nstruct;
    B->nrows = qslp->nrows;
    ILL_SAFE_MALLOC (B->cstat, qslp->nstruct, char);
    ILL_SAFE_MALLOC (B->rstat, qslp->nrows, char);

    for (i = 0; i < qslp->nstruct; i++) {
	B->cstat[i] = cstat[i];
    }

    for (i = 0; i < qslp->nrows; i++) {
	B->rstat[i] = rstat[i];
    }

    p->factorok = 0;

CLEANUP:

    ILL_RETURN (rval, "dbl_QSload_basis_array");
}

dbl_QSLIB_INTERFACE int dbl_QSload_basis_and_row_norms_array (dbl_QSdata * p,
      char *cstat,
      char *rstat,
      double *rownorms)
{
    int rval = 0;
    int i, nrows;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = p->qslp->nrows;

    rval = dbl_QSload_basis_array (p, cstat, rstat);
    ILL_CLEANUP_IF (rval);
    p->basis->rownorms = dbl_EGlpNumAllocArray (nrows);

    for (i = 0; i < nrows; i++) {
	dbl_EGlpNumCopy (p->basis->rownorms[i], rownorms[i]);
    }

    p->factorok = 0;

CLEANUP:

    ILL_RETURN (rval, "dbl_QSload_basis_and_row_norms_array");
}

dbl_QSLIB_INTERFACE int dbl_QSwrite_basis (dbl_QSdata * p,
      QSbasis * B,
      const char *filename)
{
    int rval = 0;
    dbl_ILLlp_basis iB, *basis = 0;

    dbl_ILLlp_basis_init (&iB);

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (B) {
	rval = dbl_qsbasis_to_illbasis (B, &iB);
	ILL_CLEANUP_IF (rval);
	basis = &iB;
    } else {
	if (p->basis == 0) {
	    fprintf (stderr, "no basis available in dbl_QSwrite_basis\n");
	    rval = 1;
	    goto CLEANUP;
	}
	basis = p->basis;
    }

    rval = dbl_ILLlib_writebasis (p->lp, basis, filename);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    dbl_ILLlp_basis_free (basis);
    ILL_RETURN (rval, "dbl_QSwrite_basis");
}

dbl_QSLIB_INTERFACE QSbasis *dbl_QSget_basis (dbl_QSdata * p)
{
    int rval = 0;
    QSbasis *B = 0;

    if (p->basis == 0) {
	fprintf (stderr, "no basis available in dbl_QSget_basis\n");
	rval = 1;
	goto CLEANUP;
    }
    ILL_SAFE_MALLOC (B, 1, QSbasis);
    dbl_init_basis (B);
    rval = dbl_illbasis_to_qsbasis (p->basis, B);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    if (rval) {
	dbl_QSfree_basis (B);
	B = 0;
    }
    return B;
}

dbl_QSLIB_INTERFACE int dbl_QSget_basis_array (dbl_QSdata * p,
      char *cstat,
      char *rstat)
{
    int rval = 0;
    int i;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->basis == 0) {
	fprintf (stderr, "no basis available in dbl_QSget_basis_array\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < p->basis->nstruct; i++)
	cstat[i] = p->basis->cstat[i];
    for (i = 0; i < p->basis->nrows; i++)
	rstat[i] = p->basis->rstat[i];

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_basis_array");
}

dbl_QSLIB_INTERFACE int dbl_QSget_basis_and_row_norms_array (dbl_QSdata * p,
      char *cstat,
      char *rstat,
      double *rownorms)
{
    int rval = 0;
    int i;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->basis == 0) {
	fprintf (stderr, "no basis available\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < p->basis->nstruct; i++)
	cstat[i] = p->basis->cstat[i];
    for (i = 0; i < p->basis->nrows; i++)
	rstat[i] = p->basis->rstat[i];

    if (p->basis->rownorms == 0) {
	fprintf (stderr, "no row norms available\n");
	rval = 1;
	goto CLEANUP;
    } else {
	for (i = 0; i < p->basis->nrows; i++) {
	    dbl_EGlpNumCopy (rownorms[i], p->basis->rownorms[i]);
	}
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_basis_and_row_norms_array");
}

static int dbl_illbasis_to_qsbasis (dbl_ILLlp_basis * B,
      QSbasis * qB)
{
    int rval = 0;
    int i;

    qB->nstruct = B->nstruct;
    qB->nrows = B->nrows;
    ILL_SAFE_MALLOC (qB->cstat, B->nstruct, char);
    ILL_SAFE_MALLOC (qB->rstat, B->nrows, char);

    for (i = 0; i < B->nstruct; i++) {
	qB->cstat[i] = B->cstat[i];
    }

    for (i = 0; i < B->nrows; i++) {
	qB->rstat[i] = B->rstat[i];
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_illbasis_to_qsbasis");
}

static int dbl_qsbasis_to_illbasis (QSbasis * qB,
      dbl_ILLlp_basis * B)
{
    int rval = 0;
    int i;

    B->nstruct = qB->nstruct;
    B->nrows = qB->nrows;
    ILL_SAFE_MALLOC (B->cstat, qB->nstruct, char);
    ILL_SAFE_MALLOC (B->rstat, qB->nrows, char);

    for (i = 0; i < qB->nstruct; i++) {
	B->cstat[i] = qB->cstat[i];
    }

    for (i = 0; i < qB->nrows; i++) {
	B->rstat[i] = qB->rstat[i];
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_qsbasis_to_illbasis");
}

static int dbl_grab_basis (dbl_QSdata * p)
{
    int rval = 0;
    dbl_ILLlp_basis *B = p->basis;
    int nstruct = p->qslp->nstruct;
    int nrows = p->qslp->nrows;

    if (!B) {
	ILL_SAFE_MALLOC (p->basis, 1, dbl_ILLlp_basis);
	dbl_ILLlp_basis_init (p->basis);
	B = p->basis;
    }
    if (nstruct != B->nstruct) {
	ILL_IFFREE (B->cstat, char);
	ILL_SAFE_MALLOC (B->cstat, nstruct, char);
	B->nstruct = nstruct;
    }
    if (nrows != B->nrows) {
	ILL_IFFREE (B->rstat, char);
	ILL_SAFE_MALLOC (B->rstat, nrows, char);
	B->nrows = nrows;
    }
    rval = dbl_ILLlib_getbasis (p->lp, B->cstat, B->rstat);
    ILL_CLEANUP_IF (rval);

    dbl_EGlpNumFreeArray (B->rownorms);
    dbl_EGlpNumFreeArray (B->colnorms);

    if (p->pricing->dII_price == QS_PRICE_DSTEEP) {
	B->rownorms = dbl_EGlpNumAllocArray (nrows);
	rval = dbl_ILLlib_getrownorms (p->lp, p->pricing, B->rownorms);
	if (rval) {
	    /*
	                fprintf (stderr, "no dbl_edge norms, continue anyway\n");
	    */
	    dbl_EGlpNumFreeArray (B->rownorms);
	    rval = 0;
	}
    }
CLEANUP:

    if (rval) {
	if (B) {
	    dbl_ILLlp_basis_free (B);
	    ILL_IFFREE (p->basis, dbl_ILLlp_basis);
	}
    }
    ILL_RETURN (rval, "dbl_grab_basis");
}

int dbl_grab_cache (dbl_QSdata * p,
      int status)
{
    int rval = 0;
    dbl_ILLlp_cache *C = p->cache;
    int nstruct = p->qslp->nstruct;
    int nrows = p->qslp->nrows;

    if (C == 0) {
	ILL_SAFE_MALLOC (p->cache, 1, dbl_ILLlp_cache);
	dbl_EGlpNumInitVar (p->cache->val);
	dbl_ILLlp_cache_init (p->cache);
	C = p->cache;
    }
    if (nstruct != C->nstruct || nrows != C->nrows) {
	dbl_ILLlp_cache_free (C);
	rval = dbl_ILLlp_cache_alloc (C, nstruct, nrows);
	ILL_CLEANUP_IF (rval);
    }
    rval = dbl_ILLlib_cache_solution (p->lp, C);
    ILL_CLEANUP_IF (rval);

    C->status = status;

CLEANUP:

    if (rval) {
	if (C) {
	    dbl_ILLlp_cache_free (C);
	    dbl_EGlpNumClearVar (p->cache->val);
	    ILL_IFFREE (p->cache, dbl_ILLlp_cache);
	}
    }
    ILL_RETURN (rval, "dbl_grab_cache");
}

void dbl_free_cache (dbl_QSdata * p)
{
    if (p->cache) {
	dbl_ILLlp_cache_free (p->cache);
	dbl_EGlpNumClearVar (p->cache->val);
	ILL_IFFREE (p->cache, dbl_ILLlp_cache);
    }
    p->qstatus = QS_LP_MODIFIED;
}

dbl_QSLIB_INTERFACE int dbl_QSget_binv_row (dbl_QSdata * p,
      int indx,
      double *binvrow)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "LP has not been optimized in dbl_QSget_binv_row\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_tableau (p->lp, indx, binvrow, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_binv_row");
}

dbl_QSLIB_INTERFACE int dbl_QSget_tableau_row (dbl_QSdata * p,
      int indx,
      double *tableaurow)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "LP has not been optimized in dbl_QSget_tableau_row\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_tableau (p->lp, indx, 0, tableaurow);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_tableau_row");
}

dbl_QSLIB_INTERFACE int dbl_QSget_basis_order (dbl_QSdata * p,
      int *basorder)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "LP has not been optimized in dbl_QSget_basis_order\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_basis_order (p->lp, basorder);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_basis_order");
}

dbl_QSLIB_INTERFACE int dbl_QScompute_row_norms (dbl_QSdata * p)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing->dII_price != QS_PRICE_DSTEEP) {
	fprintf (stderr, "not using dual steepest dbl_edge\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_recompute_rownorms (p->lp, p->pricing);
    ILL_CLEANUP_IF (rval);

    rval = dbl_grab_basis (p);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QScompute_row_norms");
}

dbl_QSLIB_INTERFACE void dbl_QSfree_prob (dbl_QSdata * p)
{
    if (p) {
	if (p->qslp) {
	    dbl_ILLlpdata_free (p->qslp);
	    ILL_IFFREE (p->qslp, dbl_ILLlpdata);
	}
	if (p->lp) {
	    dbl_ILLsimplex_free_lpinfo (p->lp);
	    dbl_EGlpNumClearVar (p->lp->objval);
	    dbl_EGlpNumClearVar (p->lp->pobjval);
	    dbl_EGlpNumClearVar (p->lp->dobjval);
	    dbl_EGlpNumClearVar (p->lp->pinfeas);
	    dbl_EGlpNumClearVar (p->lp->dinfeas);
	    dbl_EGlpNumClearVar (p->lp->objbound);
	    dbl_EGlpNumClearVar (p->lp->upd.piv);
	    dbl_EGlpNumClearVar (p->lp->upd.dty);
	    dbl_EGlpNumClearVar (p->lp->upd.c_obj);
	    dbl_EGlpNumClearVar (p->lp->upd.tz);
	    ILL_IFFREE (p->lp, dbl_lpinfo);
	}
	if (p->basis) {
	    dbl_ILLlp_basis_free (p->basis);
	    ILL_IFFREE (p->basis, dbl_ILLlp_basis);
	}
	if (p->cache) {
	    dbl_ILLlp_cache_free (p->cache);
	    dbl_EGlpNumClearVar (p->cache->val);
	    ILL_IFFREE (p->cache, dbl_ILLlp_cache);
	}
	if (p->pricing) {
	    dbl_EGlpNumClearVar (p->pricing->htrigger);
	    dbl_ILLprice_free_pricing_info (p->pricing);
	    ILL_IFFREE (p->pricing, dbl_price_info);
	}
	ILL_IFFREE (p->name, char);
	ILL_IFFREE (p, dbl_QSdata);
    }
}

dbl_QSLIB_INTERFACE void dbl_QSfree_basis (QSbasis * B)
{
    if (B) {
	ILL_IFFREE (B->rstat, char);
	ILL_IFFREE (B->cstat, char);
	ILL_IFFREE (B, QSbasis);
    }
}

static void dbl_init_basis (QSbasis * B)
{
    if (B) {
	B->nstruct = 0;
	B->nrows = 0;
	B->cstat = 0;
	B->rstat = 0;
    }
}

dbl_QSLIB_INTERFACE int dbl_QSget_status (dbl_QSdata * p,
      int *status)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (status)
	*status = p->qstatus;

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_status");
}

dbl_QSLIB_INTERFACE int dbl_QSget_solution (dbl_QSdata * p,
      double *value,
      double *x,
      double *pi,
      double *slack,
      double *rc)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_solution\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_solution (p->lp, p->cache, value, x, pi, slack, rc);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_solution");
}

dbl_QSLIB_INTERFACE int dbl_QSget_objval (dbl_QSdata * p,
      double *value)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    /* Want to get objval after limited number of pivots. */

    if (p->qstatus == QS_LP_MODIFIED) {
	fprintf (stderr, "QSmsg: LP has been modified since last solve.\n");
	rval = 1;
	ILL_CLEANUP;
    }
    rval = dbl_ILLlib_objval (p->lp, p->cache, value);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_objval");
}

dbl_QSLIB_INTERFACE int dbl_QSget_x_array (dbl_QSdata * p,
      double *x)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_x_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_get_x (p->lp, p->cache, x);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_x_array");
}

dbl_QSLIB_INTERFACE int dbl_QSget_slack_array (dbl_QSdata * p,
      double *slack)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_slack_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_get_slack (p->lp, p->cache, slack);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_slack_array");
}

dbl_QSLIB_INTERFACE int dbl_QSget_rc_array (dbl_QSdata * p,
      double *rc)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_rc_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_solution (p->lp, p->cache, 0, 0, 0, 0, rc);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_rc_array");
}

dbl_QSLIB_INTERFACE int dbl_QSget_pi_array (dbl_QSdata * p,
      double *pi)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_pi_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_ILLlib_solution (p->lp, p->cache, 0, 0, pi, 0, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_pi_array");
}

dbl_QSLIB_INTERFACE int dbl_QSget_infeas_array (dbl_QSdata * p,
      double *pi)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (pi == 0) {
	ILL_ERROR (rval, "QS_get_infeas_array called with NULL pi vector\n");
    }
    rval = dbl_ILLsimplex_infcertificate (p->lp, pi);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_infeas_array");
}

dbl_QSLIB_INTERFACE int dbl_QSget_named_x (dbl_QSdata * p,
      const char *colname,
      double *val)
{
    int rval = 0;
    int j;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_named_x\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_QSget_column_index (p, colname, &j);
    ILL_CLEANUP_IF (rval);

    if (j != -1) {
	dbl_EGlpNumCopy (*val, p->cache->x[j]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_named_x");
}

dbl_QSLIB_INTERFACE int dbl_QSget_named_rc (dbl_QSdata * p,
      const char *colname,
      double *val)
{
    int rval = 0;
    int j;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_named_rc\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_QSget_column_index (p, colname, &j);
    ILL_CLEANUP_IF (rval);

    if (j != -1) {
	dbl_EGlpNumCopy (*val, p->cache->rc[j]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_named_rc");
}

dbl_QSLIB_INTERFACE int dbl_QSget_named_pi (dbl_QSdata * p,
      const char *rowname,
      double *val)
{
    int rval = 0;
    int i;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_named_pi\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_QSget_row_index (p, rowname, &i);
    ILL_CLEANUP_IF (rval);

    if (i != -1) {
	dbl_EGlpNumCopy (*val, p->cache->pi[i]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_named_pi");
}

dbl_QSLIB_INTERFACE int dbl_QSget_named_slack (dbl_QSdata * p,
      const char *rowname,
      double *val)
{
    int rval = 0;
    int i;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in dbl_QSget_named_slack\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = dbl_QSget_row_index (p, rowname, &i);
    ILL_CLEANUP_IF (rval);

    if (i != -1) {
	dbl_EGlpNumCopy (*val, p->cache->slack[i]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_named_slack");
}

dbl_QSLIB_INTERFACE int dbl_QSget_colcount (dbl_QSdata * p)
{
    int cnt;

    if (dbl_check_qsdata_pointer (p))
	cnt = 0;
    else
	cnt = p->qslp->nstruct;

    return cnt;
}

dbl_QSLIB_INTERFACE int dbl_QSget_rowcount (dbl_QSdata * p)
{
    int cnt;

    if (dbl_check_qsdata_pointer (p))
	cnt = 0;
    else
	cnt = p->qslp->nrows;

    return cnt;
}

dbl_QSLIB_INTERFACE int dbl_QSget_nzcount (dbl_QSdata * p)
{
    int cnt;

    if (dbl_check_qsdata_pointer (p))
	cnt = 0;
    else
	cnt = p->qslp->nzcount - p->qslp->nrows;

    return cnt;
}

dbl_QSLIB_INTERFACE int dbl_QStest_row_norms (dbl_QSdata * p)
{
    int yesno;

    if (dbl_check_qsdata_pointer (p)) {
	yesno = 0;
    } else {
	if (p->basis && p->basis->rownorms) {
	    yesno = 1;
	} else {
	    yesno = 0;
	}
    }

    return yesno;
}

dbl_QSLIB_INTERFACE int dbl_QSget_obj (dbl_QSdata * p,
      double *obj)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_getobj (p->lp, obj);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_obj");
}

dbl_QSLIB_INTERFACE int dbl_QSget_rhs (dbl_QSdata * p,
      double *rhs)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_getrhs (p->lp, rhs);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_rhs");
}

dbl_QSLIB_INTERFACE int dbl_QSget_rows_list (dbl_QSdata * p,
      int num,
      int *rowlist,
      int **rowcnt,
      int **rowbeg,
      int **rowind,
      double **rowval,
      double **rhs,
      char **sense,
      char ***names)
{
    int rval = 0;
    int i, nrows;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = dbl_QSget_rowcount (p);
    for (i = 0; i < num; i++) {
	if (rowlist[i] < 0 || rowlist[i] >= nrows) {
	    fprintf (stderr, "entry %d in rowlist out of range\n", i);
	    rval = 1;
	    goto CLEANUP;
	}
    }

    rval = dbl_ILLlib_getrows (p->lp, num, rowlist, rowcnt, rowbeg, rowind,
	rowval, rhs, sense, names);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_rows_list");
}

dbl_QSLIB_INTERFACE int dbl_QSget_rows (dbl_QSdata * p,
      int **rowcnt,
      int **rowbeg,
      int **rowind,
      double **rowval,
      double **rhs,
      char **sense,
      char ***names)
{
    int rval = 0;
    int i, nrows;
    int *rowlist = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = dbl_QSget_rowcount (p);
    if (nrows > 0) {
	ILL_SAFE_MALLOC (rowlist, nrows, int);
	for (i = 0; i < nrows; i++) {
	    rowlist[i] = i;
	}
	rval = dbl_ILLlib_getrows (p->lp, nrows, rowlist, rowcnt, rowbeg, rowind,
	    rowval, rhs, sense, names);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (rowlist, int);
    ILL_RETURN (rval, "dbl_QSget_rows");
}

dbl_QSLIB_INTERFACE int dbl_QSget_columns_list (dbl_QSdata * p,
      int num,
      int *collist,
      int **colcnt,
      int **colbeg,
      int **colind,
      double **colval,
      double **obj,
      double **lower,
      double **upper,
      char ***names)
{
    int rval = 0;
    int j, ncols;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = dbl_QSget_colcount (p);
    for (j = 0; j < num; j++) {
	if (collist[j] < 0 || collist[j] >= ncols) {
	    fprintf (stderr, "entry %d in collist out of range\n", j);
	    rval = 1;
	    goto CLEANUP;
	}
    }

    rval = dbl_ILLlib_getcols (p->lp, num, collist, colcnt, colbeg, colind,
	colval, obj, lower, upper, names);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_columns_list");
}

dbl_QSLIB_INTERFACE int dbl_QSget_columns (dbl_QSdata * p,
      int **colcnt,
      int **colbeg,
      int **colind,
      double **colval,
      double **obj,
      double **lower,
      double **upper,
      char ***names)
{
    int rval = 0;
    int j, ncols;
    int *collist = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = dbl_QSget_colcount (p);
    if (ncols > 0) {
	ILL_SAFE_MALLOC (collist, ncols, int);
	for (j = 0; j < ncols; j++) {
	    collist[j] = j;
	}
	rval = dbl_ILLlib_getcols (p->lp, ncols, collist, colcnt, colbeg, colind,
	    colval, obj, lower, upper, names);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (collist, int);
    ILL_RETURN (rval, "dbl_QSget_columns");
}

dbl_QSLIB_INTERFACE char *dbl_QSget_probname (dbl_QSdata * p)
{
    int rval = 0;
    char *name = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    dbl_ILL_UTIL_STR (name, p->name);

CLEANUP:
    ILL_RETURN_PTR (name, "dbl_QSget_probname");
}

dbl_QSLIB_INTERFACE char *dbl_QSget_objname (dbl_QSdata * p)
{
    int rval = 0;
    char *name = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->qslp->objname != 0) {
	dbl_ILL_UTIL_STR (name, p->qslp->objname);
    }
CLEANUP:
    ILL_RETURN_PTR (name, "dbl_QSget_objname");
}

dbl_QSLIB_INTERFACE int dbl_QSget_rownames (dbl_QSdata * p,
      char **rownames)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_rownames (p->lp, rownames);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_rownames");
}

dbl_QSLIB_INTERFACE int dbl_QSget_colnames (dbl_QSdata * p,
      char **colnames)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_colnames (p->lp, colnames);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_colnames");
}

dbl_QSLIB_INTERFACE int dbl_QSget_bound (dbl_QSdata * p,
      int colindex,
      int lu,
      double *bound)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_getbnd (p->lp, colindex, lu, bound);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_bound");
}

dbl_QSLIB_INTERFACE int dbl_QSget_bounds (dbl_QSdata * p,
      double *lower,
      double *upper)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_getbnds (p->lp, lower, upper);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_bounds");
}

dbl_QSLIB_INTERFACE int dbl_QSget_intflags (dbl_QSdata * p,
      int *intflags)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (intflags == 0) {
	rval = 1;
	ILL_CLEANUP;
    }
    rval = dbl_ILLlib_getintflags (p->lp, intflags);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_intflags");
}

dbl_QSLIB_INTERFACE int dbl_QSget_intcount (dbl_QSdata * p,
      int *count)
{
    int j, ncols, cnt = 0, rval = 0;
    int *intflags = 0;

    *count = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = dbl_QSget_colcount (p);

    if (ncols > 0) {
	ILL_SAFE_MALLOC (intflags, ncols, int);

	rval = dbl_ILLlib_getintflags (p->lp, intflags);
	ILL_CLEANUP_IF (rval);

	for (j = 0; j < ncols; j++) {
	    if (intflags[j] > 0)
		cnt++;
	}
    }
CLEANUP:

    *count = cnt;
    ILL_IFFREE (intflags, int);
    ILL_RETURN (rval, "dbl_QSget_intcount");
}

dbl_QSLIB_INTERFACE int dbl_QSget_column_index (dbl_QSdata * p,
      const char *name,
      int *colindex)
{
    int rval = 0;

    *colindex = -1;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_colindex (p->lp, name, colindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_column_index");
}

dbl_QSLIB_INTERFACE int dbl_QSget_row_index (dbl_QSdata * p,
      const char *name,
      int *rowindex)
{
    int rval = 0;

    *rowindex = -1;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLlib_rowindex (p->lp, name, rowindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_row_index");
}

dbl_QSLIB_INTERFACE int dbl_QSwrite_prob (dbl_QSdata * p,
      const char *filename,
      const char *filetype)
{
    FILE *file = NULL;
    int rval = 0;

    if (filename == NULL) {
	file = stdout;
    } else {
	if ((file = fopen (filename, "w")) == 0) {
	    perror (filename);
	    fprintf (stderr, "Unable to open \"%s\" for output.\n", filename);
	}
    }
    ILL_CHECKnull (file, NULL);
    rval = dbl_QSwrite_prob_file (p, file, filetype);
    if (file && (file != stdout) && (file != stderr)) {
	fclose (file);
    }
CLEANUP:
    ILL_RETURN (rval, "dbl_QSwrite_prob_file");
}

dbl_QSLIB_INTERFACE int dbl_QSwrite_prob_file (dbl_QSdata * p,
      FILE * out,
      const char *filetype)
{
    int rval = 0;
    qsstring_reporter rep;

    ILLstring_reporter_copy (&rep, &p->qslp->reporter);
    ILLstring_reporter_init (&p->qslp->reporter,
	(qsreport_string_fct) fprintf, out);
    rval = dbl_QSreport_prob (p, filetype, NULL);
    ILLstring_reporter_copy (&p->qslp->reporter, &rep);
    ILL_RESULT (rval, "dbl_QSwrite_prob_file");
}

dbl_QSLIB_INTERFACE void dbl_QSfree (void *ptr)
{
    ILL_IFFREE (ptr, void);
}

dbl_QSLIB_INTERFACE int dbl_QSset_param (dbl_QSdata * p,
      int whichparam,
      int newvalue)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    switch (whichparam) {
    case QS_PARAM_PRIMAL_PRICING:
	if (newvalue == QS_PRICE_PDANTZIG ||
	    newvalue == QS_PRICE_PDEVEX ||
	    newvalue == QS_PRICE_PSTEEP || newvalue == QS_PRICE_PMULTPARTIAL) {
	    p->pricing->pI_price = newvalue;
	    p->pricing->pII_price = newvalue;
	} else {
	    fprintf (stderr, "illegal value for QS_PARAM_PRIMAL_PRICING\n");
	    rval = 1;
	    goto CLEANUP;
	}
	break;
    case QS_PARAM_DUAL_PRICING:
	if (newvalue == QS_PRICE_DDANTZIG ||
	    newvalue == QS_PRICE_DSTEEP ||
	    newvalue == QS_PRICE_DMULTPARTIAL || newvalue == QS_PRICE_DDEVEX) {
	    p->pricing->dI_price = newvalue;
	    p->pricing->dII_price = newvalue;
	} else {
	    fprintf (stderr, "illegal value for QS_PARAM_DUAL_PRICING\n");
	    rval = 1;
	    goto CLEANUP;
	}
	break;
    case QS_PARAM_SIMPLEX_DISPLAY:
	if (newvalue == 0 || (newvalue > 0 && newvalue < 4)) {
	    p->simplex_display = newvalue;
	} else {
	    fprintf (stderr, "illegal value for QS_PARAM_SIMPLEX_DISPLAY\n");
	    rval = 1;
	    goto CLEANUP;
	}
	break;
    case QS_PARAM_SIMPLEX_MAX_ITERATIONS:
	if (newvalue > 0) {
	    p->lp->maxiter = newvalue;
	} else {
	    fprintf (stderr, "illegal value for QS_PARAM_SIMPLEX_MAX_ITERATIONS\n");
	    rval = 1;
	    goto CLEANUP;
	}
	break;
    case QS_PARAM_SIMPLEX_SCALING:
	if (newvalue == 0 || newvalue == 1) {
	    p->simplex_scaling = newvalue;
	} else {
	    fprintf (stderr, "illegal value for QS_PARAM_SIMPLEX_SCALING\n");
	    rval = 1;
	    goto CLEANUP;
	}
	break;
    default:
	fprintf (stderr, "unknown parameter: %d\n", whichparam);
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSset_param");
}

dbl_QSLIB_INTERFACE int dbl_QSset_param_EGlpNum (dbl_QSdata * p,
      int whichparam,
      double newvalue)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    switch (whichparam) {
    case QS_PARAM_SIMPLEX_MAX_TIME:
	if (dbl_EGlpNumIsLess (dbl_zeroLpNum, newvalue)) {
	    p->lp->maxtime = dbl_EGlpNumToLf (newvalue);
	} else {
	    fprintf (stderr, "illegal value for QS_PARAM_SIMPLEX_MAX_TIME\n");
	    rval = 1;
	    goto CLEANUP;
	}
	break;
    default:
	fprintf (stderr, "unknown parameter: %d\n", whichparam);
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "QSset_param_double");
}

dbl_QSLIB_INTERFACE int dbl_QSget_param (dbl_QSdata * p,
      int whichparam,
      int *value)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (!value) {
	fprintf (stderr, "dbl_QSget_param call without a value pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    switch (whichparam) {
    case QS_PARAM_PRIMAL_PRICING:
	*value = p->pricing->pII_price;
	break;
    case QS_PARAM_DUAL_PRICING:
	*value = p->pricing->dII_price;
	break;
    case QS_PARAM_SIMPLEX_DISPLAY:
	*value = p->simplex_display;
	break;
    case QS_PARAM_SIMPLEX_MAX_ITERATIONS:
	*value = p->lp->maxiter;
	break;
    case QS_PARAM_SIMPLEX_SCALING:
	*value = p->simplex_scaling;
	break;
    default:
	fprintf (stderr, "unknown parameter: %d\n", whichparam);
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "dbl_QSget_param");
}

dbl_QSLIB_INTERFACE int dbl_QSget_param_EGlpNum (dbl_QSdata * p,
      int whichparam,
      double *value)
{
    int rval = 0;

    rval = dbl_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (!value) {
	fprintf (stderr, "QSget_param_double call without a value pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    switch (whichparam) {
    case QS_PARAM_SIMPLEX_MAX_TIME:
	dbl_EGlpNumSet (*value, p->lp->maxtime);
	break;
    default:
	fprintf (stderr, "unknown parameter: %d\n", whichparam);
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "QSget_param_double");
}

static int dbl_check_qsdata_pointer (dbl_QSdata * p)
{
    if (p == NULL) {
	fprintf (stderr, "NULL dbl_QSprob pointer\n");
	return 1;
    } else {
	return 0;
    }
}

static int dbl_formatIsMps (const char *filetype,
      int *isMps)
{
    int rval = 0;

    if (!strcasecmp (filetype, "MPS")) {
	*isMps = 1;
    } else if (!strcasecmp (filetype, "LP")) {
	*isMps = 0;
    } else {
	fprintf (stderr, "Unknown prob-file type: %s\n", filetype);
	rval = 1;
	ILL_CLEANUP;
    }
CLEANUP:
    return rval;
}


/****************************************************************************/
/* undocumentyed functions */

static void dbl_check_pointer (void *p,
      const char *fct,
      const char *param)
{
    if (p == NULL)
	fprintf (stderr, "NULL %s argument to %s\n", param, fct);
}

/* dbl_QSline_reader: used by mps/lp reader to get input lines by default
   input is read froma FILE* via fgets */
dbl_QSLIB_INTERFACE dbl_QSline_reader dbl_QSline_reader_new (void *fct,
      void *data_src)
{
    dbl_check_pointer (fct, "dbl_QSline_reader_new", "fct");
    dbl_check_pointer (data_src, "dbl_QSline_reader_new", "data_src");
    return dbl_ILLline_reader_new ((dbl_qsread_line_fct) fct, data_src);
}

dbl_QSLIB_INTERFACE void dbl_QSline_reader_set_error_collector (dbl_QSline_reader reader,
      dbl_QSerror_collector
      collector)
{
    dbl_check_pointer (reader, "dbl_QSline_reader_set_error_collector", "reader");
    dbl_check_pointer (collector, "dbl_QSline_reader_set_error_collector", "collector");
    reader->error_collector = collector;
}

dbl_QSLIB_INTERFACE void dbl_QSline_reader_free (dbl_QSline_reader reader)
{
    dbl_ILLline_reader_free (reader);
}

dbl_QSLIB_INTERFACE char *dbl_QSline_reader_get (dbl_QSline_reader reader,
      char *s,
      int size)
{
    dbl_check_pointer (reader, "dbl_QSline_reader_get", "reader");
    dbl_check_pointer (s, "dbl_QSline_reader_get", "s");
    return dbl_ILLline_reader_get (s, size, reader);
}


dbl_QSLIB_INTERFACE dbl_QSerror_collector dbl_QSerror_collector_new (void *fct,
      void *dest)
{
    dbl_check_pointer (fct, "dbl_QSerror_collector_new", "fct");
    return dbl_ILLerror_collector_new ((dbl_qsadd_error_fct) fct, dest);
}

dbl_QSLIB_INTERFACE
dbl_QSerror_collector dbl_QSerror_memory_collector_new (dbl_QSerror_memory mem)
{
    dbl_check_pointer (mem, "dbl_QSerror_memory_collector_new", "mem");
    return dbl_ILLerror_memory_collector_new (mem);
}

dbl_QSLIB_INTERFACE void dbl_QSerror_collector_free (dbl_QSerror_collector c)
{
    dbl_ILLerror_collector_free (c);
}

dbl_QSLIB_INTERFACE dbl_QSdata *dbl_QSget_prob (dbl_QSline_reader reader,
      const char *probname,
      const char *filetype)
{
    int isMps, rval = 0;
    dbl_QSdata *p = 0;

    if ((filetype != NULL) && !strcasecmp (filetype, "MPS")) {
	isMps = 1;
    } else if ((filetype != NULL) && !strcasecmp (filetype, "LP")) {
	isMps = 0;
    } else {
	fprintf (stderr, "Unknown prob-file type: %s\n",
	    (filetype != NULL) ? filetype : "NULL");
	rval = 1;
	ILL_CLEANUP;
    }

    p = dbl_ILLread (reader, probname, isMps);
    ILL_CHECKnull (p, NULL);

    ILL_FAILfalse (p->qslp != NULL, "If there's a p there must be a p-qslp");
    ILL_IFFREE (p->name, char);
    dbl_ILL_UTIL_STR (p->name, p->qslp->probname);
    dbl_ILLsimplex_load_lpinfo (p->qslp, p->lp);

CLEANUP:

    if (rval != 0) {
	dbl_QSfree_prob (p);
	p = 0;
    }
    return p;
}

dbl_QSLIB_INTERFACE int dbl_QSreport_prob (dbl_QSdata * p,
      const char *filetype,
      dbl_qserror_collector * c)
{
    int isMps, rval = 0;

    rval = dbl_formatIsMps (filetype, &isMps);
    ILL_CLEANUP_IF (rval);
    if (isMps) {
	rval = dbl_ILLwrite_mps (p->qslp, c);
    } else {
	rval = dbl_ILLwrite_lp (p->qslp, c);
    }
CLEANUP:
    ILL_RESULT (rval, "dbl_QSreport_prob");
}

dbl_QSLIB_INTERFACE char *dbl_QSversion (void)
{
    char *name = 0;
    name = EGsMalloc (char, 256);
    snprintf (name, (size_t) 255, "QSopt_ex 2.5.0 (build %s-%s)", __DATE__, __TIME__);
    return name;
}

/* QSstring_reporter: used by solver code to report dbl_feedback by default
   dbl_feedback is sent to stdout via fprintf */
dbl_QSLIB_INTERFACE void dbl_QSset_reporter (dbl_QSprob prob,
      int skip,
      void *fct,
      void *dest)
{
    int rval = 0;
    rval = dbl_check_qsdata_pointer (prob);
    if (rval != 0)
	return;

    dbl_check_pointer (fct, "dbl_QSset_reporter", "fct");

    ILL_FAILtrue (prob->lp == NULL, "dbl_QSprob internal error: prob->lp == NULL");
    ILLstring_reporter_init (&prob->qslp->reporter,
	(qsreport_string_fct) fct, dest);

    prob->lp->iterskip = skip;
CLEANUP:
    return;
}

dbl_QSLIB_INTERFACE const char *dbl_QSformat_error_type_string (int tp)
{
    const char *type = "Error";

    if (tp == QS_DATA_ERROR) {
	type = "Data Error";
    }
    if (tp == QS_DATA_WARN) {
	type = "Data Warning";
    }
    if (tp == QS_MPS_FORMAT_ERROR) {
	type = "MPS Error";
    }
    if (tp == QS_MPS_FORMAT_WARN) {
	type = "MPS Warning";
    }
    if (tp == QS_LP_FORMAT_ERROR) {
	type = "LP Error";
    }
    if (tp == QS_LP_FORMAT_WARN) {
	type = "LP Warning";
    }
    return type;
}

dbl_QSLIB_INTERFACE int dbl_QSerror_get_type (dbl_QSformat_error error)
{
    dbl_check_pointer (error, "dbl_QSerror_get_type", "error");
    return error->type;
}

dbl_QSLIB_INTERFACE const char *dbl_QSerror_get_desc (dbl_QSformat_error error)
{
    dbl_check_pointer (error, "dbl_QSerror_get_desc", "error");
    return error->desc;
}

dbl_QSLIB_INTERFACE int dbl_QSerror_get_line_number (dbl_QSformat_error error)
{
    dbl_check_pointer (error, "dbl_QSerror_get_line_number", "error");
    return error->lineNumber;
}

dbl_QSLIB_INTERFACE int dbl_QSerror_get_pos (dbl_QSformat_error error)
{
    dbl_check_pointer (error, "dbl_QSerror_get_pos", "error");
    return error->at;
}

dbl_QSLIB_INTERFACE const char *dbl_QSerror_get_line (dbl_QSformat_error error)
{
    dbl_check_pointer (error, "dbl_QSerror_get_line", "error");
    return error->theLine;
}

dbl_QSLIB_INTERFACE void dbl_QSerror_print (FILE * f,
      dbl_QSformat_error error)
{
    dbl_check_pointer (f, "dbl_QSerror_print", "f");
    if (error == NULL) {
	fprintf (stderr, "0\n");
    } else {
	dbl_ILLformat_error_print (f, error);
    }
}

dbl_QSLIB_INTERFACE dbl_QSerror_memory dbl_QSerror_memory_create (int takeErrorLines)
{
    return dbl_ILLerror_memory_create (takeErrorLines);
}

dbl_QSLIB_INTERFACE void dbl_QSerror_memory_free (dbl_QSerror_memory mem)
{
    if (mem != NULL) {
	dbl_ILLerror_memory_free (mem);
    }
}

dbl_QSLIB_INTERFACE int dbl_QSerror_memory_get_nerrors (dbl_QSerror_memory mem)
{
    dbl_check_pointer (mem, "dbl_QSerror_memory_get_nerrors", "mem");
    return mem->nerror;
}

dbl_QSLIB_INTERFACE int dbl_QSerror_memory_get_nof (dbl_QSerror_memory mem,
      int type)
{
    dbl_check_pointer (mem, "dbl_QSerror_memory_get_nerrors", "mem");
    if (0 <= type && type < QS_INPUT_NERROR) {
	return mem->has_error[type];
    } else {
	ILL_REPRT ("bad error type");
	return 0;
    }
}

dbl_QSLIB_INTERFACE dbl_QSformat_error
  dbl_QSerror_memory_get_last_error (dbl_QSerror_memory mem)
{
    dbl_check_pointer (mem, "QSerror_memory_get_last_errors", "mem");
    return mem->error_list;
}

dbl_QSLIB_INTERFACE dbl_QSformat_error dbl_QSerror_memory_get_prev_error (dbl_QSformat_error e)
{
    dbl_check_pointer (e, "QSerror_memory_get_prev_errors", "e");
    if (e != NULL)
	e = e->next;
    return e;
}
