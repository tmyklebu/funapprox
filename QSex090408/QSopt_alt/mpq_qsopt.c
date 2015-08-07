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

/* RCS_INFO = "$RCSfile: mpq_qsopt.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/* */
/* User-level Functions                                 */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int mpq_QSopt_primal (mpq_QSdata *p, int *status)                             */
/* int mpq_QSopt_dual (mpq_QSdata *p, int *status)                               */
/* mpq_QSdata *mpq_QScreate_prob (const char *name, int objsense)                */
/* mpq_QSdata *mpq_QSread_prob (const char *filename, const char *filetype)      */
/* mpq_QSdata *mpq_QSload_prob (const char *probname, int ncols, int nrows,      */
/* int *cmatcnt, int *cmatbeg, int *cmatind, double *cmatval,        */
/* int objsense, double *obj, double *rhs, char *sense,              */
/* double *lower, double *upper, const char **colnames,              */
/* const char **rownames)                                            */
/* mpq_QSdata *mpq_QScopy_prob (mpq_QSdata *p, const char *newname)                  */
/* int mpq_QSchange_objsense (mpq_QSdata *p, int newsense)                       */
/* int mpq_QSget_objsense (mpq_QSdata *p, int *objsense)                         */
/* int mpq_QSnew_col (mpq_QSdata *p, double obj, double lower, double upper,     */
/* const char *name)                                                 */
/* int mpq_QSadd_cols (mpq_QSdata *p, int num, int *cmatcnt, int *cmatbeg,       */
/* int *cmatind, double *cmatval, double *obj, double *lower,        */
/* double *upper, const char **names)                                */
/* int mpq_QSadd_col (mpq_QSdata *p, int cnt, int *cmatind, double *cmatval,     */
/* double obj, double lower, double upper, const char *name)         */
/* int mpq_QSnew_row (mpq_QSdata *p, double rhs, const char sense, char
   *name)   */
/* int mpq_QSadd_rows (mpq_QSdata *p, int num, int *rmatcnt, int *rmatbeg,       */
/* int *rmatind, double *rmatval, double *rhs, char *sense,          */
/* char **names)                                                     */
/* int mpq_QSadd_row (mpq_QSdata *p, int cnt, int *rmatind, double *rmatval,     */
/* double rhs, char sense, const char *name)                         */
/* int mpq_QSdelete_rows (mpq_QSdata *p, int num, int *dellist)                  */
/* int mpq_QSdelete_row (mpq_QSdata *p, int rowindex)                            */
/* int mpq_QSdelete_setrows (mpq_QSdata *p, int *flags)                          */
/* int mpq_QSdelete_cols (mpq_QSdata *p, int num, int *dellist)                  */
/* int mpq_QSdelete_col (mpq_QSdata *p, int colindex)                            */
/* int mpq_QSdelete_setcols (mpq_QSdata *p, int *flags)                          */
/* int mpq_QSdelete_named_column (mpq_QSdata *p, const char *colname)            */
/* int mpq_QSdelete_named_columns_list (mpq_QSdata *p, int num,                  */
/* const char **colnames)                                            */
/* int mpq_QSdelete_named_row (mpq_QSdata *p, const char *rowname)               */
/* int mpq_QSdelete_named_rows_list (mpq_QSdata *p, int num,                     */
/* const char **rownames)                                            */
/* int mpq_QSchange_senses (mpq_QSdata *p, int num, int *rowlist, char
   *sense)   */
/* int mpq_QSchange_sense (mpq_QSdata *p, int rowindex, char sense)              */
/* int mpq_QSchange_coef (mpq_QSdata *p, int rowindex, int colindex,             */
/* double coef)                                                      */
/* int mpq_QSchange_objcoef (mpq_QSdata *p, int indx, double coef)               */
/* int mpq_QSchange_rhscoef (mpq_QSdata *p, int indx, double coef)               */
/* int mpq_QSchange_bounds (mpq_QSdata *p, int num, int *collist, char *lu,      */
/* double *bounds)                                                   */
/* int mpq_QSchange_bound (mpq_QSdata *p, int indx, char lu, double bound)       */
/* int mpq_QSwrite_basis (mpq_QSdata *p, QSbasis *B, const char *filename)       */
/* QSbasis *mpq_QSget_basis (mpq_QSdata *p)                                      */
/* QSbasis *mpq_QSread_basis (mpq_QSdata *p, const char *filename)               */
/* int mpq_QSload_basis (mpq_QSdata *p, QSbasis *B)                              */
/* int mpq_QSread_and_load_basis (mpq_QSdata *p, const char *filename)           */
/* int mpq_QSload_basis_array (mpq_QSdata *p, char *cstat, char *rstat)          */
/* int mpq_QSload_basis_and_row_norms_array (mpq_QSdata *p, char *cstat,         */
/* char *rstat, double *rownorms)                                      */
/* int mpq_QSget_basis_array (mpq_QSdata *p, char *cstat, char *rstat)           */
/* int mpq_QSget_basis_and_row_norms_array (mpq_QSdata *p, char *cstat,          */
/* char *rstat, double *rownorms)                                      */
/* int mpq_QSget_binv_row (mpq_QSdata *p, int indx, double *binvrow)             */
/* int mpq_QSget_tableau_row (mpq_QSdata *p, int indx, double *tableaurow)       */
/* int mpq_QSget_basis_order (mpq_QSdata *p, int *basorder)                      */
/* int mpq_QSget_status (mpq_QSdata *p, int *status)                             */
/* int mpq_QSget_solution (mpq_QSdata *p, double *value, double *x,              */
/* double *pi, double *slack, double *rc),                           */
/* int mpq_QSget_objval (mpq_QSdata *p, double *value)                           */
/* int mpq_QSget_x_array (mpq_QSdata *p, double *x)                              */
/* int mpq_QSget_rc_array (mpq_QSdata *p, double *rc)                            */
/* int mpq_QSget_pi_array (mpq_QSdata *p, double *pi)                            */
/* int mpq_QSget_slack_array (mpq_QSdata *p, double *slack)                      */
/* int mpq_QSget_infeas_array (mpq_QSdata *p, double *pi)                        */
/* int mpq_QSget_named_x (mpq_QSdata *p, const char *colname, double *val)       */
/* int mpq_QSget_named_rc (mpq_QSdata *p, const char *colname, double *val)      */
/* int mpq_QSget_named_pi (mpq_QSdata *p, const char *rowname, double *val)      */
/* int mpq_QSget_named_slack (mpq_QSdata *p, const char *rowname, double
   *val)   */
/* int mpq_QSget_colcount (mpq_QSdata *p)                                        */
/* int mpq_QSget_rowcount (mpq_QSdata *p)                                        */
/* int mpq_QSget_nzcount (mpq_QSdata *p)                                         */
/* int mpq_QSget_obj (mpq_QSdata *p, double *obj),                               */
/* int mpq_QSget_rhs (mpq_QSdata *p, double *rhs)                                */
/* char* mpq_QSget_probname (mpq_QSdata *p)                                      */
/* char* mpq_QSget_objname (mpq_QSdata *p)                                       */
/* int mpq_QSget_columns (mpq_QSdata *p, int **colcnt, int **colbeg,             */
/* int **colind, double **colval, double **obj, double **lower,      */
/* double **upper, char ***names)                                    */
/* int mpq_QSget_columns_list (mpq_QSdata *p, int num, int *collist,             */
/* int **colcnt, int **colbeg, int **colind, double **colval,        */
/* double **obj, double **lower, double **upper, char ***names)      */
/* int mpq_QSget_rows (mpq_QSdata *p, int **rowcnt, int **rowbeg, int
   **rowind,  */
/* double **rowval, double **rhs, char **sense, char ***names)       */
/* int mpq_QSget_rows_list (mpq_QSdata *p, int num, int *rowlist, int
   **rowcnt,  */
/* int **rowbeg, int **rowind, double **rowval, double **rhs,        */
/* char **sense, char ***names)                                      */
/* int mpq_QSget_column_index (mpq_QSdata *p, const char *name, int
   *colindex)   */
/* int mpq_QSget_row_index (mpq_QSdata *p, const char *name, int *rowindex)      */
/* int mpq_QSget_rownames (mpq_QSdata *p, char **rownames)                       */
/* int mpq_QSget_colnames (mpq_QSdata *p, char **colnames)                       */
/* int mpq_QSget_bound (mpq_QSdata *p, int colindex, char lu, double *bound)     */
/* int mpq_QSget_bounds (mpq_QSdata *p, double *lower, double *upper)            */
/* int mpq_QSget_intcount (mpq_QSdata *p, int *count)                            */
/* int mpq_QSget_intflags (mpq_QSdata *p, int *intflags)                         */
/* int mpq_QScompute_row_norms (mpq_QSdata *p)                                   */
/* void mpq_QSfree_prob (mpq_QSdata *p)                                          */
/* void mpq_QSfree_basis (QSbasis *B)                                        */
/* int mpq_QSwrite_prob (mpq_QSdata *p, const char *filename,                    */
/* const char *filetype)                                             */
/* int mpq_QSwrite_prob_file (mpq_QSdata *p, FILE *file, const char
   *filetype)   */
/* int mpq_QSset_param (mpq_QSdata *p, int whichparam, int newvalue)             */
/* int QSset_param_double (mpq_QSdata *p, int whichparam, double newvalue)   */
/* int mpq_QSget_param (mpq_QSdata *p, int whichparam, int *value)               */
/* int QSget_param_double (mpq_QSdata *p, int whichparam, double *value)     */
/* int mpq_QStest_row_norms (mpq_QSdata *p)                                      */
/* int mpq_QSopt_strongbranch (mpq_QSdata *p, int ncand, int *candidatelist,     */
/* double *xlist, double *down_vals, double *up_vals,                */
/* int iterations, double objbound)                                  */
/* int mpq_QSopt_pivotin_row (mpq_QSdata *p, int rcnt, int *rlist)               */
/* int mpq_QSopt_pivotin_col (mpq_QSdata *p, int ccnt, int *clist)               */
/* void mpq_QSfree (void *ptr)                                               */
/* void mpq_QSstart (void)                                                   */
/* void mpq_QSend (void)                                                     */
/* char *mpq_QSversion (void))                                               */
/* */
/* mpq_NEW FUNCTIONS - Add to Docs                                           */
/* */
/* char *mpq_QSversion (void))                                               */
/* int mpq_QSget_objsense (mpq_QSdata *p, int *objsense)                         */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdata.h"
#include "mpq_lpdefs.h"
#include "mpq_simplex.h"
#include "mpq_price.h"
#include "mpq_qstruct.h"
#include "mpq_qsopt.h"
#include "mpq_lib.h"
#include "mpq_mps.h"
#include "mpq_lp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
void mpq_QSset_precision (const unsigned prec)
{
    EGlpNumSetPrecision (prec);
    mpq_ILLchange_precision ();
    /* change the numbers */
}

static void mpq_init_basis (QSbasis * B),
  mpq_free_cache (mpq_QSdata * p);

int mpq_grab_cache (mpq_QSdata * p,
      int status);
static int mpq_opt_work (mpq_QSdata * p,
      int *status,
      int primal_or_dual),
  mpq_qsbasis_to_illbasis (QSbasis * qB,
      mpq_ILLlp_basis * B),
  mpq_illbasis_to_qsbasis (mpq_ILLlp_basis * B,
      QSbasis * qB),
  mpq_grab_basis (mpq_QSdata * p),
  mpq_check_qsdata_pointer (mpq_QSdata * p);


mpq_QSLIB_INTERFACE int mpq_QSopt_primal (mpq_QSdata * p,
      int *status)
{
    int rval = 0;

    if (status)
	*status = QS_LP_UNSOLVED;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    /* If both the basis and the cache exist, then skip the optimization */

    if (!p->basis || !p->cache) {
	rval = mpq_opt_work (p, status, 0);
	ILL_CLEANUP_IF (rval);
    } else {
	if (status)
	    *status = p->cache->status;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSopt_primal");
}

mpq_QSLIB_INTERFACE int mpq_QSopt_dual (mpq_QSdata * p,
      int *status)
{
    int rval = 0;

    if (status)
	*status = QS_LP_UNSOLVED;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (!p->basis || !p->cache || !p->factorok) {
	rval = mpq_opt_work (p, status, 1);
	ILL_CLEANUP_IF (rval);
    } else {
	if (status)
	    *status = p->cache->status;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSopt_dual");
}

static int mpq_opt_work (mpq_QSdata * p,
      int *status,
      int primal_or_dual)
{
    int rval = 0;
    int rstatus = QS_LP_UNSOLVED;
    mpq_QSdata *p2 = 0;

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

	mpq_ILLprice_free_pricing_info (p->pricing);	/* Just to be sure  */
	p->factorok = 0;	/* that p is clean. */

	p2 = mpq_QScopy_prob (p, "scaled_lp");
	if (p2 == 0)
	    goto CLEANUP;

	rval = mpq_ILLlp_scale (p2->qslp);
	ILL_CLEANUP_IF (rval);

	if (primal_or_dual == 0) {
	    rval = mpq_ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
		PRIMAL_SIMPLEX, 0, p2->simplex_display);
	} else {
	    rval = mpq_ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
		DUAL_SIMPLEX, 0, p2->simplex_display);
	}
	ILL_CLEANUP_IF (rval);

	rval = mpq_grab_basis (p2);
	ILL_CLEANUP_IF (rval);

	if (p->basis) {
	    mpq_ILLlp_basis_free (p->basis);
	    ILL_IFFREE (p->basis, mpq_ILLlp_basis);
	}
	p->basis = p2->basis;
	p2->basis = 0;
	mpq_QSfree_prob (p2);
	p2 = 0;
    }
    if (primal_or_dual == 0) {
	if (p->factorok == 0) {
	    if (p->basis == 0)
		p->lp->basisid = -1;
	    rval = mpq_ILLlib_optimize (p->lp, p->basis, p->pricing, PRIMAL_SIMPLEX,
		&rstatus, p->simplex_display);
	} else {
	    mpq_ILLprice_free_pricing_info (p->pricing);
	    if (p->lp->basisid != -1)
		p->lp->fbasisid = p->lp->basisid;
	    rval = mpq_ILLlib_optimize (p->lp, 0, p->pricing,
		PRIMAL_SIMPLEX, &rstatus, p->simplex_display);
	}
    } else {
	if (p->factorok == 0) {
	    if (p->basis == 0)
		p->lp->basisid = -1;
	    rval = mpq_ILLlib_optimize (p->lp, p->basis, p->pricing, DUAL_SIMPLEX,
		&rstatus, p->simplex_display);
	} else {
	    /* The factorization and rownorms should be up-to-date */
	    if (p->lp->basisid != -1) {
		p->lp->fbasisid = p->lp->basisid;
	    } else {
		mpq_ILLprice_free_pricing_info (p->pricing);
	    }
	    rval = mpq_ILLlib_optimize (p->lp, 0, p->pricing,
		DUAL_SIMPLEX, &rstatus, p->simplex_display);
	}
    }
    ILL_CLEANUP_IF (rval);

    rval = mpq_grab_basis (p);
    ILL_CLEANUP_IF (rval);

    if (rstatus == QS_LP_OPTIMAL) {
	rval = mpq_grab_cache (p, rstatus);
	ILL_CLEANUP_IF (rval);
    } else {
	mpq_free_cache (p);
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
	mpq_QSfree_prob (p2);
    ILL_RETURN (rval, "mpq_opt_work");
}

mpq_QSLIB_INTERFACE int mpq_QSopt_pivotin_row (mpq_QSdata * p,
      int rcnt,
      int *rlist)
{
    int basismod = 0;
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing == 0) {
	ILL_ERROR (rval, "pricing info not available in mpq_QSopt_pivotin_row\n");
    }
    rval = mpq_ILLsimplex_pivotin (p->lp, p->pricing, rcnt, rlist,
	SIMPLEX_PIVOTINROW, &basismod);
    ILL_CLEANUP_IF (rval);

    rval = mpq_grab_basis (p);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSopt_pivotin_row");
}

mpq_QSLIB_INTERFACE int mpq_QSopt_pivotin_col (mpq_QSdata * p,
      int ccnt,
      int *clist)
{
    int basismod = 0;
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing == 0) {
	ILL_ERROR (rval, "pricing info not available in QSopt_pivotin\n");
    }
    rval = mpq_ILLsimplex_pivotin (p->lp, p->pricing, ccnt, clist,
	SIMPLEX_PIVOTINCOL, &basismod);
    ILL_CLEANUP_IF (rval);

    rval = mpq_grab_basis (p);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSopt_pivotin_col");
}


mpq_QSLIB_INTERFACE int mpq_QSopt_strongbranch (mpq_QSdata * p,
      int ncand,
      int *candidatelist,
      mpq_t * xlist,
      mpq_t * down_vals,
      mpq_t * up_vals,
      int iterations,
      mpq_t objbound)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing == 0) {
	rval = 1;
	ILL_CLEANUP_IF (rval);
    }
    rval = mpq_ILLlib_strongbranch (p->lp, p->pricing, candidatelist, ncand,
	xlist, down_vals, up_vals, iterations, objbound);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;
    mpq_free_cache (p);
    p->qstatus = QS_LP_UNSOLVED;/* Was set to MODIFIED in mpq_free_cache () */

CLEANUP:

    ILL_RETURN (rval, "mpq_QSopt_strongbranch");
}

mpq_QSLIB_INTERFACE mpq_QSdata *mpq_QScreate_prob (const char *name,
      int objsense)
{
    int rval = 0;
    mpq_QSdata *p = 0;
    int len;

    ILL_SAFE_MALLOC (p, 1, mpq_QSdata);
    if (!p) {
	fprintf (stderr, "out of memory in mpq_QScreate_prob\n");
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

    ILL_SAFE_MALLOC (p->qslp, 1, mpq_ILLlpdata);
    if (!p->qslp) {
	fprintf (stderr, "out of memory in mpq_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    mpq_ILLlpdata_init (p->qslp);

    ILL_SAFE_MALLOC (p->lp, 1, mpq_lpinfo);
    if (!p->lp) {
	fprintf (stderr, "out of memory in mpq_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    mpq_EGlpNumInitVar (p->lp->objval);
    mpq_EGlpNumInitVar (p->lp->pobjval);
    mpq_EGlpNumInitVar (p->lp->dobjval);
    mpq_EGlpNumInitVar (p->lp->pinfeas);
    mpq_EGlpNumInitVar (p->lp->dinfeas);
    mpq_EGlpNumInitVar (p->lp->objbound);
    mpq_EGlpNumInitVar (p->lp->upd.piv);
    mpq_EGlpNumInitVar (p->lp->upd.dty);
    mpq_EGlpNumInitVar (p->lp->upd.c_obj);
    mpq_EGlpNumInitVar (p->lp->upd.tz);
    mpq_ILLsimplex_init_lpinfo (p->lp);
    mpq_ILLsimplex_load_lpinfo (p->qslp, p->lp);

    ILL_SAFE_MALLOC (p->pricing, 1, mpq_price_info);
    if (!p->pricing) {
	fprintf (stderr, "out of memory in mpq_QScreate_prob\n");
	rval = 1;
	goto CLEANUP;
    }
    mpq_EGlpNumInitVar (p->pricing->htrigger);
    mpq_ILLprice_init_pricing_info (p->pricing);
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
	mpq_QSfree_prob (p);
	p = 0;
    }
    return p;
}

mpq_QSLIB_INTERFACE mpq_QSdata *mpq_QSread_prob (const char *filename,
      const char *filetype)
{
    mpq_QSdata *p = 0;
    FILE *file = 0;
    mpq_QSline_reader reader;

    if ((file = fopen (filename, "r")) == 0) {
	perror (filename);
	fprintf (stderr, "Unable to open \"%s\" for input.\n", filename);
    }
    if (file == NULL)
	goto CLEANUP;

    reader = mpq_ILLline_reader_new ((mpq_qsread_line_fct) fgets, file);
    p = mpq_QSget_prob (reader, filename, filetype);
    mpq_QSline_reader_free (reader);	/* Bico - 040723 */

CLEANUP:
    if (file != NULL) {
	fclose (file);
    }
    return p;
}

mpq_QSLIB_INTERFACE mpq_QSdata *mpq_QSload_prob (const char *probname,
      int ncols,
      int nrows,
      int *cmatcnt,
      int *cmatbeg,
      int *cmatind,
      mpq_t * cmatval,
      int objsense,
      mpq_t * obj,
      mpq_t * rhs,
      char *sense,
      mpq_t * lower,
      mpq_t * upper,
      const char **colnames,
      const char **rownames)
{
    int rval = 0;
    mpq_QSdata *p = 0;

    p = mpq_QScreate_prob (probname, objsense);
    if (p == 0)
	goto CLEANUP;

    rval = mpq_ILLlib_newrows (p->lp, 0, nrows, rhs, sense, 0, rownames);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_addcols (p->lp, 0, ncols, cmatcnt, cmatbeg, cmatind,
	cmatval, obj, lower, upper, colnames, 0);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;

CLEANUP:

    if (rval) {
	mpq_QSfree_prob (p);
	p = 0;
    }
    return p;
}

mpq_QSLIB_INTERFACE mpq_QSdata *mpq_QScopy_prob (mpq_QSdata * p,
      const char *newname)
{
    int rval = 0;
    int j, col, beg, pindex, hit;
    mpq_QSdata *p2 = 0;
    char *coln;
    char buf[ILL_namebufsize];

    /* printf ("mpq_QScopy_prob ...\n"); fflush (stdout); */

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    p2 = mpq_QScreate_prob (newname, p->qslp->objsense);
    if (p2 == 0)
	goto CLEANUP;

    rval = mpq_ILLlib_newrows (p2->lp, 0, p->qslp->nrows,
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
	/* code in mpq_lib.c.                                  */

	rval = mpq_ILLlib_addcol (p2->lp, 0,
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
    mpq_EGlpNumClearVar (p2->pricing->htrigger);
    *(p2->pricing) = *(p->pricing);
    /* I added this line because copying the mpq_heap (as a pointer) doesn't
       make any sense ! */
    mpq_ILLheap_init (&(p2->pricing->h));
    mpq_EGlpNumInitVar (p2->pricing->htrigger);
    mpq_EGlpNumCopy (p2->pricing->htrigger, p->pricing->htrigger);

    if (p->qslp->intmarker != 0) {
	ILL_SAFE_MALLOC (p2->qslp->intmarker, p->qslp->nstruct, char);
	for (j = 0; j < p->qslp->nstruct; j++) {
	    p2->qslp->intmarker[j] = p->qslp->intmarker[j];
	}
    }
    if (p->qslp->objname != 0) {
	mpq_ILL_UTIL_STR (p2->qslp->objname, p->qslp->objname);
    } else {
	strcpy (buf, "obj");
	rval = ILLsymboltab_uname (&p2->qslp->rowtab, buf, "", NULL);
	ILL_CLEANUP_IF (rval);
	mpq_ILL_UTIL_STR (p2->qslp->objname, buf);
    }
    rval = ILLsymboltab_register (&p2->qslp->rowtab, p2->qslp->objname,
	-1, &pindex, &hit);
    rval = rval || hit;
    ILL_CLEANUP_IF (rval);

    ILLstring_reporter_copy (&p2->qslp->reporter, &p->qslp->reporter);

CLEANUP:

    if (rval) {
	mpq_QSfree_prob (p2);
	p2 = 0;
    }
    return p2;
}

mpq_QSLIB_INTERFACE int mpq_QSchange_objsense (mpq_QSdata * p,
      int newsense)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (newsense != QS_MIN && newsense != QS_MAX) {
	fprintf (stderr, "Illegal objective sense %d\n", newsense);
	rval = 1;
	goto CLEANUP;
    }
    if (p->qslp->objsense != newsense) {
	p->qslp->objsense = newsense;
	mpq_free_cache (p);
    }
CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_objsense");
}

mpq_QSLIB_INTERFACE int mpq_QSget_objsense (mpq_QSdata * p,
      int *objsense)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (objsense)
	*objsense = p->qslp->objsense;

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_objsense");
}


mpq_QSLIB_INTERFACE int mpq_QSnew_col (mpq_QSdata * p,
      mpq_t obj,
      mpq_t lower,
      mpq_t upper,
      const char *name)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_newcol (p->lp, p->basis, obj, lower, upper, name, p->factorok);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSnew_col");
}

mpq_QSLIB_INTERFACE int mpq_QSadd_cols (mpq_QSdata * p,
      int num,
      int *cmatcnt,
      int *cmatbeg,
      int *cmatind,
      mpq_t * cmatval,
      mpq_t * obj,
      mpq_t * lower,
      mpq_t * upper,
      const char **names)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_addcols (p->lp, p->basis, num, cmatcnt, cmatbeg,
	cmatind, cmatval, obj, lower, upper, names,
	p->factorok);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSadd_cols");
}

mpq_QSLIB_INTERFACE int mpq_QSadd_col (mpq_QSdata * p,
      int cnt,
      int *cmatind,
      mpq_t * cmatval,
      mpq_t obj,
      mpq_t lower,
      mpq_t upper,
      const char *name)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_addcol (p->lp, p->basis, cnt, cmatind, cmatval,
	obj, lower, upper, name, p->factorok);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSadd_col");
}

mpq_QSLIB_INTERFACE int mpq_QSnew_row (mpq_QSdata * p,
      mpq_t rhs,
      int sense,
      const char *name)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_newrow (p->lp, p->basis, rhs, sense, mpq_zeroLpNum, name);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;
    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSnew_row");
}

mpq_QSLIB_INTERFACE int mpq_QSadd_rows (mpq_QSdata * p,
      int num,
      int *rmatcnt,
      int *rmatbeg,
      int *rmatind,
      mpq_t * rmatval,
      mpq_t * rhs,
      char *sense,
      const char **names)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_addrows (p->lp, p->basis, num, rmatcnt, rmatbeg,
	rmatind, rmatval, rhs, sense, 0,
	names, &(p->factorok));
    ILL_CLEANUP_IF (rval);

    if (p->factorok == 1 && p->basis->rownorms) {
	rval = mpq_ILLlib_loadrownorms (p->lp, p->pricing, p->basis->rownorms);
	ILL_CLEANUP_IF (rval);
	/* This really should go inside of mpq_ILLlib_addrows, once pinf is  */
	/* is moved into the lp struct.                                  */
    }
    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSadd_rows");
}

mpq_QSLIB_INTERFACE int mpq_QSadd_row (mpq_QSdata * p,
      int cnt,
      int *rmatind,
      mpq_t * rmatval,
      mpq_t * rhs,
      int sense,
      const char *name)
{
    int rval = 0;
    int vmatcnt[1];
    int vmatbeg[1];
    char vsense[1];
    const char *vnames[1];

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vmatcnt[0] = cnt;
    vmatbeg[0] = 0;
    vsense[0] = sense;
    vnames[0] = name;

    rval = mpq_QSadd_rows (p, 1, vmatcnt, vmatbeg, rmatind, rmatval, rhs, vsense,
	vnames);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSadd_row");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_rows (mpq_QSdata * p,
      int num,
      int *dellist)
{
    int rval = 0;
    int basis_ok = 0;
    int cache_ok = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_delrows (p->lp, p->basis, p->cache, num, dellist, &basis_ok,
	&cache_ok);
    ILL_CLEANUP_IF (rval);

    /* For now, just remove the basis - wait for pivotin */

    if (p->basis && !basis_ok) {
	mpq_ILLlp_basis_free (p->basis);
	ILL_IFFREE (p->basis, mpq_ILLlp_basis);
    }
    p->factorok = 0;

    if (!p->basis || !basis_ok || !cache_ok) {
	/* Note: If we only delete basic rows then cached soln is valid */
	mpq_free_cache (p);
    }
CLEANUP:

    ILL_RETURN (rval, "mpq_QSdelete_rows");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_row (mpq_QSdata * p,
      int rowindex)
{
    int rval = 0;
    int vdellist[1];

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = rowindex;

    rval = mpq_QSdelete_rows (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSdelete_row");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_setrows (mpq_QSdata * p,
      int *flags)
{
    int rval = 0;
    int j, num = 0;
    int *dellist = 0;
    int nrows;

    rval = mpq_check_qsdata_pointer (p);
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

	rval = mpq_QSdelete_rows (p, num, dellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (dellist, int);
    ILL_RETURN (rval, "mpq_QSdelete_setrows");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_named_row (mpq_QSdata * p,
      const char *rowname)
{
    int rval = 0;
    int i, vdellist[1];

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_QSget_row_index (p, rowname, &i);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = i;

    rval = mpq_QSdelete_rows (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSdelete_named_row");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_named_rows_list (mpq_QSdata * p,
      int num,
      const char **rownames)
{
    int rval = 0;
    int i, k;
    int *vdellist = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (num > 0) {
	ILL_SAFE_MALLOC (vdellist, num, int);
	for (k = 0; k < num; k++) {
	    rval = mpq_QSget_row_index (p, rownames[k], &i);
	    ILL_CLEANUP_IF (rval);
	    vdellist[k] = i;
	}

	rval = mpq_QSdelete_rows (p, num, vdellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (vdellist, int);
    ILL_RETURN (rval, "mpq_QSdelete_named_rows_list");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_cols (mpq_QSdata * p,
      int num,
      int *dellist)
{
    int rval = 0;
    int basis_ok;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_delcols (p->lp, p->basis, num, dellist, &basis_ok);
    ILL_CLEANUP_IF (rval);

    /* For now, just remove the basis - wait for pivotout */

    if (p->basis && !basis_ok) {
	mpq_ILLlp_basis_free (p->basis);
	ILL_IFFREE (p->basis, mpq_ILLlp_basis);
    }
    p->factorok = 0;
    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSdelete_cols");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_col (mpq_QSdata * p,
      int colindex)
{
    int rval = 0;
    int vdellist[1];

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = colindex;

    rval = mpq_QSdelete_cols (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSdelete_col");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_setcols (mpq_QSdata * p,
      int *flags)
{
    int rval = 0;
    int j, num = 0;
    int *dellist = 0;
    int ncols;

    rval = mpq_check_qsdata_pointer (p);
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

	rval = mpq_QSdelete_cols (p, num, dellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (dellist, int);
    ILL_RETURN (rval, "mpq_QSdelete_setcols");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_named_column (mpq_QSdata * p,
      const char *colname)
{
    int rval = 0;
    int j, vdellist[1];

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_QSget_column_index (p, colname, &j);
    ILL_CLEANUP_IF (rval);

    vdellist[0] = j;

    rval = mpq_QSdelete_cols (p, 1, vdellist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSdelete_named_column");
}

mpq_QSLIB_INTERFACE int mpq_QSdelete_named_columns_list (mpq_QSdata * p,
      int num,
      const char **colnames)
{
    int rval = 0;
    int i, j;
    int *vdellist = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (num > 0) {
	ILL_SAFE_MALLOC (vdellist, num, int);
	for (i = 0; i < num; i++) {
	    rval = mpq_QSget_column_index (p, colnames[i], &j);
	    ILL_CLEANUP_IF (rval);
	    vdellist[i] = j;
	}

	rval = mpq_QSdelete_cols (p, num, vdellist);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (vdellist, int);
    ILL_RETURN (rval, "mpq_QSdelete_named_columns_list");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_senses (mpq_QSdata * p,
      int num,
      int *rowlist,
      char *sense)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_chgsense (p->lp, num, rowlist, sense);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_senses");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_sense (mpq_QSdata * p,
      int rowindex,
      int sense)
{
    int rval = 0;
    int vrowlist[1];
    char vsenselist[1];

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    vrowlist[0] = rowindex;
    vsenselist[0] = sense;

    rval = mpq_QSchange_senses (p, 1, vrowlist, vsenselist);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_sense");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_coef (mpq_QSdata * p,
      int rowindex,
      int colindex,
      mpq_t coef)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_chgcoef (p->lp, rowindex, colindex, coef);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_coef");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_objcoef (mpq_QSdata * p,
      int indx,
      mpq_t coef)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_chgobj (p->lp, indx, coef);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_objcoef");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_rhscoef (mpq_QSdata * p,
      int indx,
      mpq_t coef)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_chgrhs (p->lp, indx, coef);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_rhscoef");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_bounds (mpq_QSdata * p,
      int num,
      int *collist,
      char *lu,
      mpq_t * bounds)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_chgbnds (p->lp, num, collist, lu, bounds);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_bounds");
}

mpq_QSLIB_INTERFACE int mpq_QSchange_bound (mpq_QSdata * p,
      int indx,
      int lu,
      mpq_t bound)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_chgbnd (p->lp, indx, lu, bound);
    ILL_CLEANUP_IF (rval);

    mpq_free_cache (p);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSchange_bound");
}

#if 0
/*
 * Bico - I removed this on 02.04.22.  I don't think we need to support
 * this type of interface (the loading via arrays can do the job)
 */
mpq_QSLIB_INTERFACE QSbasis *QScreate_basis (int nstruct,
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
	mpq_QSfree_basis (B);
	B = 0;
    }
    return B;
}
#endif

mpq_QSLIB_INTERFACE QSbasis *mpq_QSread_basis (mpq_QSdata * p,
      const char *filename)
{
    int rval = 0;
    QSbasis *qB = 0;
    mpq_ILLlp_basis B;

    mpq_ILLlp_basis_init (&B);

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ILL_NEW (qB, QSbasis);
    mpq_init_basis (qB);

    rval = mpq_ILLlib_readbasis (p->lp, &B, filename);
    ILL_CLEANUP_IF (rval);

    rval = mpq_illbasis_to_qsbasis (&B, qB);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    if (rval && qB) {
	mpq_QSfree_basis (qB);
	qB = 0;
    }
    mpq_ILLlp_basis_free (&B);

    return qB;
}

mpq_QSLIB_INTERFACE int mpq_QSload_basis (mpq_QSdata * p,
      QSbasis * B)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (B->nstruct != p->qslp->nstruct || B->nrows != p->qslp->nrows) {
	fprintf (stderr, "size of basis does not match lp\n");
	rval = 1;
	goto CLEANUP;
    }
    if (p->basis == 0) {
	ILL_SAFE_MALLOC (p->basis, 1, mpq_ILLlp_basis);
	mpq_ILLlp_basis_init (p->basis);
    } else {
	mpq_ILLlp_basis_free (p->basis);
    }

    rval = mpq_qsbasis_to_illbasis (B, p->basis);
    ILL_CLEANUP_IF (rval);

    p->factorok = 0;

CLEANUP:

    ILL_RETURN (rval, "mpq_QSload_basis");
}

mpq_QSLIB_INTERFACE int mpq_QSread_and_load_basis (mpq_QSdata * p,
      const char *filename)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->basis == 0) {
	ILL_SAFE_MALLOC (p->basis, 1, mpq_ILLlp_basis);
	mpq_ILLlp_basis_init (p->basis);
    } else {
	mpq_ILLlp_basis_free (p->basis);
    }

    rval = mpq_ILLlib_readbasis (p->lp, p->basis, filename);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    return rval;
}

mpq_QSLIB_INTERFACE int mpq_QSload_basis_array (mpq_QSdata * p,
      char *cstat,
      char *rstat)
{
    int rval = 0;
    int i;
    mpq_ILLlp_basis *B;
    mpq_ILLlpdata *qslp;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    qslp = p->qslp;

    if (qslp->nstruct > 0 && cstat == 0) {
	fprintf (stderr, "mpq_QSload_basis_array called without cstat\n");
	rval = 1;
	goto CLEANUP;
    }
    if (qslp->nrows > 0 && rstat == 0) {
	fprintf (stderr, "mpq_QSload_basis_array called without rstat\n");
	rval = 1;
	goto CLEANUP;
    }
    if (p->basis == 0) {
	ILL_SAFE_MALLOC (p->basis, 1, mpq_ILLlp_basis);
	mpq_ILLlp_basis_init (p->basis);
    } else {
	mpq_ILLlp_basis_free (p->basis);
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

    ILL_RETURN (rval, "mpq_QSload_basis_array");
}

mpq_QSLIB_INTERFACE int mpq_QSload_basis_and_row_norms_array (mpq_QSdata * p,
      char *cstat,
      char *rstat,
      mpq_t * rownorms)
{
    int rval = 0;
    int i, nrows;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = p->qslp->nrows;

    rval = mpq_QSload_basis_array (p, cstat, rstat);
    ILL_CLEANUP_IF (rval);
    p->basis->rownorms = mpq_EGlpNumAllocArray (nrows);

    for (i = 0; i < nrows; i++) {
	mpq_EGlpNumCopy (p->basis->rownorms[i], rownorms[i]);
    }

    p->factorok = 0;

CLEANUP:

    ILL_RETURN (rval, "mpq_QSload_basis_and_row_norms_array");
}

mpq_QSLIB_INTERFACE int mpq_QSwrite_basis (mpq_QSdata * p,
      QSbasis * B,
      const char *filename)
{
    int rval = 0;
    mpq_ILLlp_basis iB, *basis = 0;

    mpq_ILLlp_basis_init (&iB);

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (B) {
	rval = mpq_qsbasis_to_illbasis (B, &iB);
	ILL_CLEANUP_IF (rval);
	basis = &iB;
    } else {
	if (p->basis == 0) {
	    fprintf (stderr, "no basis available in mpq_QSwrite_basis\n");
	    rval = 1;
	    goto CLEANUP;
	}
	basis = p->basis;
    }

    rval = mpq_ILLlib_writebasis (p->lp, basis, filename);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    mpq_ILLlp_basis_free (basis);
    ILL_RETURN (rval, "mpq_QSwrite_basis");
}

mpq_QSLIB_INTERFACE QSbasis *mpq_QSget_basis (mpq_QSdata * p)
{
    int rval = 0;
    QSbasis *B = 0;

    if (p->basis == 0) {
	fprintf (stderr, "no basis available in mpq_QSget_basis\n");
	rval = 1;
	goto CLEANUP;
    }
    ILL_SAFE_MALLOC (B, 1, QSbasis);
    mpq_init_basis (B);
    rval = mpq_illbasis_to_qsbasis (p->basis, B);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    if (rval) {
	mpq_QSfree_basis (B);
	B = 0;
    }
    return B;
}

mpq_QSLIB_INTERFACE int mpq_QSget_basis_array (mpq_QSdata * p,
      char *cstat,
      char *rstat)
{
    int rval = 0;
    int i;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->basis == 0) {
	fprintf (stderr, "no basis available in mpq_QSget_basis_array\n");
	rval = 1;
	goto CLEANUP;
    }
    for (i = 0; i < p->basis->nstruct; i++)
	cstat[i] = p->basis->cstat[i];
    for (i = 0; i < p->basis->nrows; i++)
	rstat[i] = p->basis->rstat[i];

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_basis_array");
}

mpq_QSLIB_INTERFACE int mpq_QSget_basis_and_row_norms_array (mpq_QSdata * p,
      char *cstat,
      char *rstat,
      mpq_t * rownorms)
{
    int rval = 0;
    int i;

    rval = mpq_check_qsdata_pointer (p);
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
	    mpq_EGlpNumCopy (rownorms[i], p->basis->rownorms[i]);
	}
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_basis_and_row_norms_array");
}

static int mpq_illbasis_to_qsbasis (mpq_ILLlp_basis * B,
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

    ILL_RETURN (rval, "mpq_illbasis_to_qsbasis");
}

static int mpq_qsbasis_to_illbasis (QSbasis * qB,
      mpq_ILLlp_basis * B)
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

    ILL_RETURN (rval, "mpq_qsbasis_to_illbasis");
}

static int mpq_grab_basis (mpq_QSdata * p)
{
    int rval = 0;
    mpq_ILLlp_basis *B = p->basis;
    int nstruct = p->qslp->nstruct;
    int nrows = p->qslp->nrows;

    if (!B) {
	ILL_SAFE_MALLOC (p->basis, 1, mpq_ILLlp_basis);
	mpq_ILLlp_basis_init (p->basis);
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
    rval = mpq_ILLlib_getbasis (p->lp, B->cstat, B->rstat);
    ILL_CLEANUP_IF (rval);

    mpq_EGlpNumFreeArray (B->rownorms);
    mpq_EGlpNumFreeArray (B->colnorms);

    if (p->pricing->dII_price == QS_PRICE_DSTEEP) {
	B->rownorms = mpq_EGlpNumAllocArray (nrows);
	rval = mpq_ILLlib_getrownorms (p->lp, p->pricing, B->rownorms);
	if (rval) {
	    /*
	                fprintf (stderr, "no mpq_edge norms, continue anyway\n");
	    */
	    mpq_EGlpNumFreeArray (B->rownorms);
	    rval = 0;
	}
    }
CLEANUP:

    if (rval) {
	if (B) {
	    mpq_ILLlp_basis_free (B);
	    ILL_IFFREE (p->basis, mpq_ILLlp_basis);
	}
    }
    ILL_RETURN (rval, "mpq_grab_basis");
}

int mpq_grab_cache (mpq_QSdata * p,
      int status)
{
    int rval = 0;
    mpq_ILLlp_cache *C = p->cache;
    int nstruct = p->qslp->nstruct;
    int nrows = p->qslp->nrows;

    if (C == 0) {
	ILL_SAFE_MALLOC (p->cache, 1, mpq_ILLlp_cache);
	mpq_EGlpNumInitVar (p->cache->val);
	mpq_ILLlp_cache_init (p->cache);
	C = p->cache;
    }
    if (nstruct != C->nstruct || nrows != C->nrows) {
	mpq_ILLlp_cache_free (C);
	rval = mpq_ILLlp_cache_alloc (C, nstruct, nrows);
	ILL_CLEANUP_IF (rval);
    }
    rval = mpq_ILLlib_cache_solution (p->lp, C);
    ILL_CLEANUP_IF (rval);

    C->status = status;

CLEANUP:

    if (rval) {
	if (C) {
	    mpq_ILLlp_cache_free (C);
	    mpq_EGlpNumClearVar (p->cache->val);
	    ILL_IFFREE (p->cache, mpq_ILLlp_cache);
	}
    }
    ILL_RETURN (rval, "mpq_grab_cache");
}

void mpq_free_cache (mpq_QSdata * p)
{
    if (p->cache) {
	mpq_ILLlp_cache_free (p->cache);
	mpq_EGlpNumClearVar (p->cache->val);
	ILL_IFFREE (p->cache, mpq_ILLlp_cache);
    }
    p->qstatus = QS_LP_MODIFIED;
}

mpq_QSLIB_INTERFACE int mpq_QSget_binv_row (mpq_QSdata * p,
      int indx,
      mpq_t * binvrow)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "LP has not been optimized in mpq_QSget_binv_row\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_tableau (p->lp, indx, binvrow, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_binv_row");
}

mpq_QSLIB_INTERFACE int mpq_QSget_tableau_row (mpq_QSdata * p,
      int indx,
      mpq_t * tableaurow)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "LP has not been optimized in mpq_QSget_tableau_row\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_tableau (p->lp, indx, 0, tableaurow);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_tableau_row");
}

mpq_QSLIB_INTERFACE int mpq_QSget_basis_order (mpq_QSdata * p,
      int *basorder)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "LP has not been optimized in mpq_QSget_basis_order\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_basis_order (p->lp, basorder);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_basis_order");
}

mpq_QSLIB_INTERFACE int mpq_QScompute_row_norms (mpq_QSdata * p)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->pricing->dII_price != QS_PRICE_DSTEEP) {
	fprintf (stderr, "not using dual steepest mpq_edge\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_recompute_rownorms (p->lp, p->pricing);
    ILL_CLEANUP_IF (rval);

    rval = mpq_grab_basis (p);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QScompute_row_norms");
}

mpq_QSLIB_INTERFACE void mpq_QSfree_prob (mpq_QSdata * p)
{
    if (p) {
	if (p->qslp) {
	    mpq_ILLlpdata_free (p->qslp);
	    ILL_IFFREE (p->qslp, mpq_ILLlpdata);
	}
	if (p->lp) {
	    mpq_ILLsimplex_free_lpinfo (p->lp);
	    mpq_EGlpNumClearVar (p->lp->objval);
	    mpq_EGlpNumClearVar (p->lp->pobjval);
	    mpq_EGlpNumClearVar (p->lp->dobjval);
	    mpq_EGlpNumClearVar (p->lp->pinfeas);
	    mpq_EGlpNumClearVar (p->lp->dinfeas);
	    mpq_EGlpNumClearVar (p->lp->objbound);
	    mpq_EGlpNumClearVar (p->lp->upd.piv);
	    mpq_EGlpNumClearVar (p->lp->upd.dty);
	    mpq_EGlpNumClearVar (p->lp->upd.c_obj);
	    mpq_EGlpNumClearVar (p->lp->upd.tz);
	    ILL_IFFREE (p->lp, mpq_lpinfo);
	}
	if (p->basis) {
	    mpq_ILLlp_basis_free (p->basis);
	    ILL_IFFREE (p->basis, mpq_ILLlp_basis);
	}
	if (p->cache) {
	    mpq_ILLlp_cache_free (p->cache);
	    mpq_EGlpNumClearVar (p->cache->val);
	    ILL_IFFREE (p->cache, mpq_ILLlp_cache);
	}
	if (p->pricing) {
	    mpq_EGlpNumClearVar (p->pricing->htrigger);
	    mpq_ILLprice_free_pricing_info (p->pricing);
	    ILL_IFFREE (p->pricing, mpq_price_info);
	}
	ILL_IFFREE (p->name, char);
	ILL_IFFREE (p, mpq_QSdata);
    }
}

mpq_QSLIB_INTERFACE void mpq_QSfree_basis (QSbasis * B)
{
    if (B) {
	ILL_IFFREE (B->rstat, char);
	ILL_IFFREE (B->cstat, char);
	ILL_IFFREE (B, QSbasis);
    }
}

static void mpq_init_basis (QSbasis * B)
{
    if (B) {
	B->nstruct = 0;
	B->nrows = 0;
	B->cstat = 0;
	B->rstat = 0;
    }
}

mpq_QSLIB_INTERFACE int mpq_QSget_status (mpq_QSdata * p,
      int *status)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (status)
	*status = p->qstatus;

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_status");
}

mpq_QSLIB_INTERFACE int mpq_QSget_solution (mpq_QSdata * p,
      mpq_t * value,
      mpq_t * x,
      mpq_t * pi,
      mpq_t * slack,
      mpq_t * rc)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_solution\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_solution (p->lp, p->cache, value, x, pi, slack, rc);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_solution");
}

mpq_QSLIB_INTERFACE int mpq_QSget_objval (mpq_QSdata * p,
      mpq_t * value)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    /* Want to get objval after limited number of pivots. */

    if (p->qstatus == QS_LP_MODIFIED) {
	fprintf (stderr, "QSmsg: LP has been modified since last solve.\n");
	rval = 1;
	ILL_CLEANUP;
    }
    rval = mpq_ILLlib_objval (p->lp, p->cache, value);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_objval");
}

mpq_QSLIB_INTERFACE int mpq_QSget_x_array (mpq_QSdata * p,
      mpq_t * x)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_x_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_get_x (p->lp, p->cache, x);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_x_array");
}

mpq_QSLIB_INTERFACE int mpq_QSget_slack_array (mpq_QSdata * p,
      mpq_t * slack)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_slack_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_get_slack (p->lp, p->cache, slack);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_slack_array");
}

mpq_QSLIB_INTERFACE int mpq_QSget_rc_array (mpq_QSdata * p,
      mpq_t * rc)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_rc_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_solution (p->lp, p->cache, 0, 0, 0, 0, rc);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_rc_array");
}

mpq_QSLIB_INTERFACE int mpq_QSget_pi_array (mpq_QSdata * p,
      mpq_t * pi)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_pi_array\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_ILLlib_solution (p->lp, p->cache, 0, 0, pi, 0, 0);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_pi_array");
}

mpq_QSLIB_INTERFACE int mpq_QSget_infeas_array (mpq_QSdata * p,
      mpq_t * pi)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (pi == 0) {
	ILL_ERROR (rval, "QS_get_infeas_array called with NULL pi vector\n");
    }
    rval = mpq_ILLsimplex_infcertificate (p->lp, pi);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_infeas_array");
}

mpq_QSLIB_INTERFACE int mpq_QSget_named_x (mpq_QSdata * p,
      const char *colname,
      mpq_t * val)
{
    int rval = 0;
    int j;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_named_x\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_QSget_column_index (p, colname, &j);
    ILL_CLEANUP_IF (rval);

    if (j != -1) {
	mpq_EGlpNumCopy (*val, p->cache->x[j]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_named_x");
}

mpq_QSLIB_INTERFACE int mpq_QSget_named_rc (mpq_QSdata * p,
      const char *colname,
      mpq_t * val)
{
    int rval = 0;
    int j;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_named_rc\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_QSget_column_index (p, colname, &j);
    ILL_CLEANUP_IF (rval);

    if (j != -1) {
	mpq_EGlpNumCopy (*val, p->cache->rc[j]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_named_rc");
}

mpq_QSLIB_INTERFACE int mpq_QSget_named_pi (mpq_QSdata * p,
      const char *rowname,
      mpq_t * val)
{
    int rval = 0;
    int i;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_named_pi\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_QSget_row_index (p, rowname, &i);
    ILL_CLEANUP_IF (rval);

    if (i != -1) {
	mpq_EGlpNumCopy (*val, p->cache->pi[i]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_named_pi");
}

mpq_QSLIB_INTERFACE int mpq_QSget_named_slack (mpq_QSdata * p,
      const char *rowname,
      mpq_t * val)
{
    int rval = 0;
    int i;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->cache == 0) {
	fprintf (stderr, "no solution available in mpq_QSget_named_slack\n");
	rval = 1;
	goto CLEANUP;
    }
    rval = mpq_QSget_row_index (p, rowname, &i);
    ILL_CLEANUP_IF (rval);

    if (i != -1) {
	mpq_EGlpNumCopy (*val, p->cache->slack[i]);
    } else {
	rval = 1;
    }

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_named_slack");
}

mpq_QSLIB_INTERFACE int mpq_QSget_colcount (mpq_QSdata * p)
{
    int cnt;

    if (mpq_check_qsdata_pointer (p))
	cnt = 0;
    else
	cnt = p->qslp->nstruct;

    return cnt;
}

mpq_QSLIB_INTERFACE int mpq_QSget_rowcount (mpq_QSdata * p)
{
    int cnt;

    if (mpq_check_qsdata_pointer (p))
	cnt = 0;
    else
	cnt = p->qslp->nrows;

    return cnt;
}

mpq_QSLIB_INTERFACE int mpq_QSget_nzcount (mpq_QSdata * p)
{
    int cnt;

    if (mpq_check_qsdata_pointer (p))
	cnt = 0;
    else
	cnt = p->qslp->nzcount - p->qslp->nrows;

    return cnt;
}

mpq_QSLIB_INTERFACE int mpq_QStest_row_norms (mpq_QSdata * p)
{
    int yesno;

    if (mpq_check_qsdata_pointer (p)) {
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

mpq_QSLIB_INTERFACE int mpq_QSget_obj (mpq_QSdata * p,
      mpq_t * obj)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_getobj (p->lp, obj);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_obj");
}

mpq_QSLIB_INTERFACE int mpq_QSget_rhs (mpq_QSdata * p,
      mpq_t * rhs)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_getrhs (p->lp, rhs);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_rhs");
}

mpq_QSLIB_INTERFACE int mpq_QSget_rows_list (mpq_QSdata * p,
      int num,
      int *rowlist,
      int **rowcnt,
      int **rowbeg,
      int **rowind,
      mpq_t ** rowval,
      mpq_t ** rhs,
      char **sense,
      char ***names)
{
    int rval = 0;
    int i, nrows;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = mpq_QSget_rowcount (p);
    for (i = 0; i < num; i++) {
	if (rowlist[i] < 0 || rowlist[i] >= nrows) {
	    fprintf (stderr, "entry %d in rowlist out of range\n", i);
	    rval = 1;
	    goto CLEANUP;
	}
    }

    rval = mpq_ILLlib_getrows (p->lp, num, rowlist, rowcnt, rowbeg, rowind,
	rowval, rhs, sense, names);
    ILL_CLEANUP_IF (rval);


CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_rows_list");
}

mpq_QSLIB_INTERFACE int mpq_QSget_rows (mpq_QSdata * p,
      int **rowcnt,
      int **rowbeg,
      int **rowind,
      mpq_t ** rowval,
      mpq_t ** rhs,
      char **sense,
      char ***names)
{
    int rval = 0;
    int i, nrows;
    int *rowlist = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    nrows = mpq_QSget_rowcount (p);
    if (nrows > 0) {
	ILL_SAFE_MALLOC (rowlist, nrows, int);
	for (i = 0; i < nrows; i++) {
	    rowlist[i] = i;
	}
	rval = mpq_ILLlib_getrows (p->lp, nrows, rowlist, rowcnt, rowbeg, rowind,
	    rowval, rhs, sense, names);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (rowlist, int);
    ILL_RETURN (rval, "mpq_QSget_rows");
}

mpq_QSLIB_INTERFACE int mpq_QSget_columns_list (mpq_QSdata * p,
      int num,
      int *collist,
      int **colcnt,
      int **colbeg,
      int **colind,
      mpq_t ** colval,
      mpq_t ** obj,
      mpq_t ** lower,
      mpq_t ** upper,
      char ***names)
{
    int rval = 0;
    int j, ncols;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = mpq_QSget_colcount (p);
    for (j = 0; j < num; j++) {
	if (collist[j] < 0 || collist[j] >= ncols) {
	    fprintf (stderr, "entry %d in collist out of range\n", j);
	    rval = 1;
	    goto CLEANUP;
	}
    }

    rval = mpq_ILLlib_getcols (p->lp, num, collist, colcnt, colbeg, colind,
	colval, obj, lower, upper, names);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_columns_list");
}

mpq_QSLIB_INTERFACE int mpq_QSget_columns (mpq_QSdata * p,
      int **colcnt,
      int **colbeg,
      int **colind,
      mpq_t ** colval,
      mpq_t ** obj,
      mpq_t ** lower,
      mpq_t ** upper,
      char ***names)
{
    int rval = 0;
    int j, ncols;
    int *collist = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = mpq_QSget_colcount (p);
    if (ncols > 0) {
	ILL_SAFE_MALLOC (collist, ncols, int);
	for (j = 0; j < ncols; j++) {
	    collist[j] = j;
	}
	rval = mpq_ILLlib_getcols (p->lp, ncols, collist, colcnt, colbeg, colind,
	    colval, obj, lower, upper, names);
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:

    ILL_IFFREE (collist, int);
    ILL_RETURN (rval, "mpq_QSget_columns");
}

mpq_QSLIB_INTERFACE char *mpq_QSget_probname (mpq_QSdata * p)
{
    int rval = 0;
    char *name = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    mpq_ILL_UTIL_STR (name, p->name);

CLEANUP:
    ILL_RETURN_PTR (name, "mpq_QSget_probname");
}

mpq_QSLIB_INTERFACE char *mpq_QSget_objname (mpq_QSdata * p)
{
    int rval = 0;
    char *name = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (p->qslp->objname != 0) {
	mpq_ILL_UTIL_STR (name, p->qslp->objname);
    }
CLEANUP:
    ILL_RETURN_PTR (name, "mpq_QSget_objname");
}

mpq_QSLIB_INTERFACE int mpq_QSget_rownames (mpq_QSdata * p,
      char **rownames)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_rownames (p->lp, rownames);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_rownames");
}

mpq_QSLIB_INTERFACE int mpq_QSget_colnames (mpq_QSdata * p,
      char **colnames)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_colnames (p->lp, colnames);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_colnames");
}

mpq_QSLIB_INTERFACE int mpq_QSget_bound (mpq_QSdata * p,
      int colindex,
      int lu,
      mpq_t * bound)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_getbnd (p->lp, colindex, lu, bound);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_bound");
}

mpq_QSLIB_INTERFACE int mpq_QSget_bounds (mpq_QSdata * p,
      mpq_t * lower,
      mpq_t * upper)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_getbnds (p->lp, lower, upper);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_bounds");
}

mpq_QSLIB_INTERFACE int mpq_QSget_intflags (mpq_QSdata * p,
      int *intflags)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (intflags == 0) {
	rval = 1;
	ILL_CLEANUP;
    }
    rval = mpq_ILLlib_getintflags (p->lp, intflags);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_intflags");
}

mpq_QSLIB_INTERFACE int mpq_QSget_intcount (mpq_QSdata * p,
      int *count)
{
    int j, ncols, cnt = 0, rval = 0;
    int *intflags = 0;

    *count = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    ncols = mpq_QSget_colcount (p);

    if (ncols > 0) {
	ILL_SAFE_MALLOC (intflags, ncols, int);

	rval = mpq_ILLlib_getintflags (p->lp, intflags);
	ILL_CLEANUP_IF (rval);

	for (j = 0; j < ncols; j++) {
	    if (intflags[j] > 0)
		cnt++;
	}
    }
CLEANUP:

    *count = cnt;
    ILL_IFFREE (intflags, int);
    ILL_RETURN (rval, "mpq_QSget_intcount");
}

mpq_QSLIB_INTERFACE int mpq_QSget_column_index (mpq_QSdata * p,
      const char *name,
      int *colindex)
{
    int rval = 0;

    *colindex = -1;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_colindex (p->lp, name, colindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_column_index");
}

mpq_QSLIB_INTERFACE int mpq_QSget_row_index (mpq_QSdata * p,
      const char *name,
      int *rowindex)
{
    int rval = 0;

    *rowindex = -1;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLlib_rowindex (p->lp, name, rowindex);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_RETURN (rval, "mpq_QSget_row_index");
}

mpq_QSLIB_INTERFACE int mpq_QSwrite_prob (mpq_QSdata * p,
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
    rval = mpq_QSwrite_prob_file (p, file, filetype);
    if (file && (file != stdout) && (file != stderr)) {
	fclose (file);
    }
CLEANUP:
    ILL_RETURN (rval, "mpq_QSwrite_prob_file");
}

mpq_QSLIB_INTERFACE int mpq_QSwrite_prob_file (mpq_QSdata * p,
      FILE * out,
      const char *filetype)
{
    int rval = 0;
    qsstring_reporter rep;

    ILLstring_reporter_copy (&rep, &p->qslp->reporter);
    ILLstring_reporter_init (&p->qslp->reporter,
	(qsreport_string_fct) fprintf, out);
    rval = mpq_QSreport_prob (p, filetype, NULL);
    ILLstring_reporter_copy (&p->qslp->reporter, &rep);
    ILL_RESULT (rval, "mpq_QSwrite_prob_file");
}

mpq_QSLIB_INTERFACE void mpq_QSfree (void *ptr)
{
    ILL_IFFREE (ptr, void);
}

mpq_QSLIB_INTERFACE int mpq_QSset_param (mpq_QSdata * p,
      int whichparam,
      int newvalue)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
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

    ILL_RETURN (rval, "mpq_QSset_param");
}

mpq_QSLIB_INTERFACE int mpq_QSset_param_EGlpNum (mpq_QSdata * p,
      int whichparam,
      mpq_t newvalue)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    switch (whichparam) {
    case QS_PARAM_SIMPLEX_MAX_TIME:
	if (mpq_EGlpNumIsLess (mpq_zeroLpNum, newvalue)) {
	    p->lp->maxtime = mpq_EGlpNumToLf (newvalue);
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

mpq_QSLIB_INTERFACE int mpq_QSget_param (mpq_QSdata * p,
      int whichparam,
      int *value)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (!value) {
	fprintf (stderr, "mpq_QSget_param call without a value pointer\n");
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

    ILL_RETURN (rval, "mpq_QSget_param");
}

mpq_QSLIB_INTERFACE int mpq_QSget_param_EGlpNum (mpq_QSdata * p,
      int whichparam,
      mpq_t * value)
{
    int rval = 0;

    rval = mpq_check_qsdata_pointer (p);
    ILL_CLEANUP_IF (rval);

    if (!value) {
	fprintf (stderr, "QSget_param_double call without a value pointer\n");
	rval = 1;
	goto CLEANUP;
    }
    switch (whichparam) {
    case QS_PARAM_SIMPLEX_MAX_TIME:
	mpq_EGlpNumSet (*value, p->lp->maxtime);
	break;
    default:
	fprintf (stderr, "unknown parameter: %d\n", whichparam);
	rval = 1;
	goto CLEANUP;
    }

CLEANUP:

    ILL_RETURN (rval, "QSget_param_double");
}

static int mpq_check_qsdata_pointer (mpq_QSdata * p)
{
    if (p == NULL) {
	fprintf (stderr, "NULL mpq_QSprob pointer\n");
	return 1;
    } else {
	return 0;
    }
}

static int mpq_formatIsMps (const char *filetype,
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

static void mpq_check_pointer (void *p,
      const char *fct,
      const char *param)
{
    if (p == NULL)
	fprintf (stderr, "NULL %s argument to %s\n", param, fct);
}

/* mpq_QSline_reader: used by mps/lp reader to get input lines by default
   input is read froma FILE* via fgets */
mpq_QSLIB_INTERFACE mpq_QSline_reader mpq_QSline_reader_new (void *fct,
      void *data_src)
{
    mpq_check_pointer (fct, "mpq_QSline_reader_new", "fct");
    mpq_check_pointer (data_src, "mpq_QSline_reader_new", "data_src");
    return mpq_ILLline_reader_new ((mpq_qsread_line_fct) fct, data_src);
}

mpq_QSLIB_INTERFACE void mpq_QSline_reader_set_error_collector (mpq_QSline_reader reader,
      mpq_QSerror_collector
      collector)
{
    mpq_check_pointer (reader, "mpq_QSline_reader_set_error_collector", "reader");
    mpq_check_pointer (collector, "mpq_QSline_reader_set_error_collector", "collector");
    reader->error_collector = collector;
}

mpq_QSLIB_INTERFACE void mpq_QSline_reader_free (mpq_QSline_reader reader)
{
    mpq_ILLline_reader_free (reader);
}

mpq_QSLIB_INTERFACE char *mpq_QSline_reader_get (mpq_QSline_reader reader,
      char *s,
      int size)
{
    mpq_check_pointer (reader, "mpq_QSline_reader_get", "reader");
    mpq_check_pointer (s, "mpq_QSline_reader_get", "s");
    return mpq_ILLline_reader_get (s, size, reader);
}


mpq_QSLIB_INTERFACE mpq_QSerror_collector mpq_QSerror_collector_new (void *fct,
      void *dest)
{
    mpq_check_pointer (fct, "mpq_QSerror_collector_new", "fct");
    return mpq_ILLerror_collector_new ((mpq_qsadd_error_fct) fct, dest);
}

mpq_QSLIB_INTERFACE
mpq_QSerror_collector mpq_QSerror_memory_collector_new (mpq_QSerror_memory mem)
{
    mpq_check_pointer (mem, "mpq_QSerror_memory_collector_new", "mem");
    return mpq_ILLerror_memory_collector_new (mem);
}

mpq_QSLIB_INTERFACE void mpq_QSerror_collector_free (mpq_QSerror_collector c)
{
    mpq_ILLerror_collector_free (c);
}

mpq_QSLIB_INTERFACE mpq_QSdata *mpq_QSget_prob (mpq_QSline_reader reader,
      const char *probname,
      const char *filetype)
{
    int isMps, rval = 0;
    mpq_QSdata *p = 0;

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

    p = mpq_ILLread (reader, probname, isMps);
    ILL_CHECKnull (p, NULL);

    ILL_FAILfalse (p->qslp != NULL, "If there's a p there must be a p-qslp");
    ILL_IFFREE (p->name, char);
    mpq_ILL_UTIL_STR (p->name, p->qslp->probname);
    mpq_ILLsimplex_load_lpinfo (p->qslp, p->lp);

CLEANUP:

    if (rval != 0) {
	mpq_QSfree_prob (p);
	p = 0;
    }
    return p;
}

mpq_QSLIB_INTERFACE int mpq_QSreport_prob (mpq_QSdata * p,
      const char *filetype,
      mpq_qserror_collector * c)
{
    int isMps, rval = 0;

    rval = mpq_formatIsMps (filetype, &isMps);
    ILL_CLEANUP_IF (rval);
    if (isMps) {
	rval = mpq_ILLwrite_mps (p->qslp, c);
    } else {
	rval = mpq_ILLwrite_lp (p->qslp, c);
    }
CLEANUP:
    ILL_RESULT (rval, "mpq_QSreport_prob");
}

mpq_QSLIB_INTERFACE char *mpq_QSversion (void)
{
    char *name = 0;
    name = EGsMalloc (char, 256);
    snprintf (name, (size_t) 255, "QSopt_ex 2.5.0 (build %s-%s)", __DATE__, __TIME__);
    return name;
}

/* QSstring_reporter: used by solver code to report mpq_feedback by default
   mpq_feedback is sent to stdout via fprintf */
mpq_QSLIB_INTERFACE void mpq_QSset_reporter (mpq_QSprob prob,
      int skip,
      void *fct,
      void *dest)
{
    int rval = 0;
    rval = mpq_check_qsdata_pointer (prob);
    if (rval != 0)
	return;

    mpq_check_pointer (fct, "mpq_QSset_reporter", "fct");

    ILL_FAILtrue (prob->lp == NULL, "mpq_QSprob internal error: prob->lp == NULL");
    ILLstring_reporter_init (&prob->qslp->reporter,
	(qsreport_string_fct) fct, dest);

    prob->lp->iterskip = skip;
CLEANUP:
    return;
}

mpq_QSLIB_INTERFACE const char *mpq_QSformat_error_type_string (int tp)
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

mpq_QSLIB_INTERFACE int mpq_QSerror_get_type (mpq_QSformat_error error)
{
    mpq_check_pointer (error, "mpq_QSerror_get_type", "error");
    return error->type;
}

mpq_QSLIB_INTERFACE const char *mpq_QSerror_get_desc (mpq_QSformat_error error)
{
    mpq_check_pointer (error, "mpq_QSerror_get_desc", "error");
    return error->desc;
}

mpq_QSLIB_INTERFACE int mpq_QSerror_get_line_number (mpq_QSformat_error error)
{
    mpq_check_pointer (error, "mpq_QSerror_get_line_number", "error");
    return error->lineNumber;
}

mpq_QSLIB_INTERFACE int mpq_QSerror_get_pos (mpq_QSformat_error error)
{
    mpq_check_pointer (error, "mpq_QSerror_get_pos", "error");
    return error->at;
}

mpq_QSLIB_INTERFACE const char *mpq_QSerror_get_line (mpq_QSformat_error error)
{
    mpq_check_pointer (error, "mpq_QSerror_get_line", "error");
    return error->theLine;
}

mpq_QSLIB_INTERFACE void mpq_QSerror_print (FILE * f,
      mpq_QSformat_error error)
{
    mpq_check_pointer (f, "mpq_QSerror_print", "f");
    if (error == NULL) {
	fprintf (stderr, "0\n");
    } else {
	mpq_ILLformat_error_print (f, error);
    }
}

mpq_QSLIB_INTERFACE mpq_QSerror_memory mpq_QSerror_memory_create (int takeErrorLines)
{
    return mpq_ILLerror_memory_create (takeErrorLines);
}

mpq_QSLIB_INTERFACE void mpq_QSerror_memory_free (mpq_QSerror_memory mem)
{
    if (mem != NULL) {
	mpq_ILLerror_memory_free (mem);
    }
}

mpq_QSLIB_INTERFACE int mpq_QSerror_memory_get_nerrors (mpq_QSerror_memory mem)
{
    mpq_check_pointer (mem, "mpq_QSerror_memory_get_nerrors", "mem");
    return mem->nerror;
}

mpq_QSLIB_INTERFACE int mpq_QSerror_memory_get_nof (mpq_QSerror_memory mem,
      int type)
{
    mpq_check_pointer (mem, "mpq_QSerror_memory_get_nerrors", "mem");
    if (0 <= type && type < QS_INPUT_NERROR) {
	return mem->has_error[type];
    } else {
	ILL_REPRT ("bad error type");
	return 0;
    }
}

mpq_QSLIB_INTERFACE mpq_QSformat_error
  mpq_QSerror_memory_get_last_error (mpq_QSerror_memory mem)
{
    mpq_check_pointer (mem, "QSerror_memory_get_last_errors", "mem");
    return mem->error_list;
}

mpq_QSLIB_INTERFACE mpq_QSformat_error mpq_QSerror_memory_get_prev_error (mpq_QSformat_error e)
{
    mpq_check_pointer (e, "QSerror_memory_get_prev_errors", "e");
    if (e != NULL)
	e = e->next;
    return e;
}
