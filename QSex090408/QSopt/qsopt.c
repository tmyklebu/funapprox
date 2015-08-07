/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza.  Sanjeeb Dash's ownership of copyright in          */
/*  QSopt_ex is derived from his copyright in QSopt.                        */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */ 
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

/* RCS_INFO = "$RCSfile: qsopt.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/*                                                                          */
/*                     User-level Functions                                 */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int QSopt_primal (QSdata *p, int *status)                             */
/*    int QSopt_dual (QSdata *p, int *status)                               */
/*    QSdata *QScreate_prob (const char *name, int objsense)                */
/*    QSdata *QSread_prob (const char *filename, const char *filetype)      */
/*    QSdata *QSload_prob (const char *probname, int ncols, int nrows,      */
/*        int *cmatcnt, int *cmatbeg, int *cmatind, double *cmatval,        */
/*        int objsense, double *obj, double *rhs, char *sense,              */
/*        double *lower, double *upper, const char **colnames,              */
/*        const char **rownames)                                            */
/*    QSdata *QScopy_prob (QSdata *p, const char *newname)                  */
/*    int QSchange_objsense (QSdata *p, int newsense)                       */
/*    int QSget_objsense (QSdata *p, int *objsense)                         */
/*    int QSnew_col (QSdata *p, double obj, double lower, double upper,     */
/*        const char *name)                                                 */
/*    int QSadd_cols (QSdata *p, int num, int *cmatcnt, int *cmatbeg,       */
/*        int *cmatind, double *cmatval, double *obj, double *lower,        */
/*        double *upper, const char **names)                                */
/*    int QSadd_col (QSdata *p, int cnt, int *cmatind, double *cmatval,     */
/*        double obj, double lower, double upper, const char *name)         */
/*    int QSnew_row (QSdata *p, double rhs, const char sense, char *name)   */
/*    int QSadd_rows (QSdata *p, int num, int *rmatcnt, int *rmatbeg,       */
/*        int *rmatind, double *rmatval, double *rhs, char *sense,          */
/*        char **names)                                                     */
/*    int QSadd_row (QSdata *p, int cnt, int *rmatind, double *rmatval,     */
/*        double rhs, char sense, const char *name)                         */
/*    int QSdelete_rows (QSdata *p, int num, int *dellist)                  */
/*    int QSdelete_row (QSdata *p, int rowindex)                            */
/*    int QSdelete_setrows (QSdata *p, int *flags)                          */
/*    int QSdelete_cols (QSdata *p, int num, int *dellist)                  */
/*    int QSdelete_col (QSdata *p, int colindex)                            */
/*    int QSdelete_setcols (QSdata *p, int *flags)                          */
/*    int QSdelete_named_column (QSdata *p, const char *colname)            */
/*    int QSdelete_named_columns_list (QSdata *p, int num,                  */
/*        const char **colnames)                                            */
/*    int QSdelete_named_row (QSdata *p, const char *rowname)               */
/*    int QSdelete_named_rows_list (QSdata *p, int num,                     */
/*        const char **rownames)                                            */
/*    int QSchange_senses (QSdata *p, int num, int *rowlist, char *sense)   */
/*    int QSchange_sense (QSdata *p, int rowindex, char sense)              */
/*    int QSchange_coef (QSdata *p, int rowindex, int colindex,             */
/*        double coef)                                                      */
/*    int QSchange_objcoef (QSdata *p, int indx, double coef)               */
/*    int QSchange_rhscoef (QSdata *p, int indx, double coef)               */
/*    int QSchange_bounds (QSdata *p, int num, int *collist, char *lu,      */
/*        double *bounds)                                                   */
/*    int QSchange_bound (QSdata *p, int indx, char lu, double bound)       */
/*    int QSwrite_basis (QSdata *p, QSbasis *B, const char *filename)       */
/*    QSbasis *QSget_basis (QSdata *p)                                      */
/*    QSbasis *QSread_basis (QSdata *p, const char *filename)               */
/*    int QSload_basis (QSdata *p, QSbasis *B)                              */
/*    int QSread_and_load_basis (QSdata *p, const char *filename)           */
/*    int QSload_basis_array (QSdata *p, char *cstat, char *rstat)          */
/*    int QSload_basis_and_row_norms_array (QSdata *p, char *cstat,         */
/*      char *rstat, double *rownorms)                                      */
/*    int QSget_basis_array (QSdata *p, char *cstat, char *rstat)           */
/*    int QSget_basis_and_row_norms_array (QSdata *p, char *cstat,          */
/*      char *rstat, double *rownorms)                                      */
/*    int QSget_binv_row (QSdata *p, int indx, double *binvrow)             */
/*    int QSget_tableau_row (QSdata *p, int indx, double *tableaurow)       */
/*    int QSget_basis_order (QSdata *p, int *basorder)                      */
/*    int QSget_status (QSdata *p, int *status)                             */
/*    int QSget_solution (QSdata *p, double *value, double *x,              */
/*        double *pi, double *slack, double *rc),                           */
/*    int QSget_objval (QSdata *p, double *value)                           */
/*    int QSget_x_array (QSdata *p, double *x)                              */
/*    int QSget_rc_array (QSdata *p, double *rc)                            */
/*    int QSget_pi_array (QSdata *p, double *pi)                            */
/*    int QSget_slack_array (QSdata *p, double *slack)                      */
/*    int QSget_infeas_array (QSdata *p, double *pi)                        */
/*    int QSget_named_x (QSdata *p, const char *colname, double *val)       */
/*    int QSget_named_rc (QSdata *p, const char *colname, double *val)      */
/*    int QSget_named_pi (QSdata *p, const char *rowname, double *val)      */
/*    int QSget_named_slack (QSdata *p, const char *rowname, double *val)   */
/*    int QSget_colcount (QSdata *p)                                        */
/*    int QSget_rowcount (QSdata *p)                                        */
/*    int QSget_nzcount (QSdata *p)                                         */
/*    int QSget_obj (QSdata *p, double *obj),                               */
/*    int QSget_rhs (QSdata *p, double *rhs)                                */
/*    char* QSget_probname (QSdata *p)                                      */
/*    char* QSget_objname (QSdata *p)                                       */
/*    int QSget_columns (QSdata *p, int **colcnt, int **colbeg,             */
/*        int **colind, double **colval, double **obj, double **lower,      */
/*        double **upper, char ***names)                                    */
/*    int QSget_columns_list (QSdata *p, int num, int *collist,             */
/*        int **colcnt, int **colbeg, int **colind, double **colval,        */
/*        double **obj, double **lower, double **upper, char ***names)      */
/*    int QSget_rows (QSdata *p, int **rowcnt, int **rowbeg, int **rowind,  */
/*        double **rowval, double **rhs, char **sense, char ***names)       */
/*    int QSget_rows_list (QSdata *p, int num, int *rowlist, int **rowcnt,  */
/*        int **rowbeg, int **rowind, double **rowval, double **rhs,        */
/*        char **sense, char ***names)                                      */
/*    int QSget_column_index (QSdata *p, const char *name, int *colindex)   */
/*    int QSget_row_index (QSdata *p, const char *name, int *rowindex)      */
/*    int QSget_rownames (QSdata *p, char **rownames)                       */
/*    int QSget_colnames (QSdata *p, char **colnames)                       */
/*    int QSget_bound (QSdata *p, int colindex, char lu, double *bound)     */
/*    int QSget_bounds (QSdata *p, double *lower, double *upper)            */
/*    int QSget_intcount (QSdata *p, int *count)                            */
/*    int QSget_intflags (QSdata *p, int *intflags)                         */
/*    int QScompute_row_norms (QSdata *p)                                   */
/*    void QSfree_prob (QSdata *p)                                          */
/*    void QSfree_basis (QSbasis *B)                                        */
/*    int QSwrite_prob (QSdata *p, const char *filename,                    */
/*        const char *filetype)                                             */
/*    int QSwrite_prob_file (QSdata *p, FILE *file, const char *filetype)   */
/*    int QSset_param (QSdata *p, int whichparam, int newvalue)             */
/*    int QSset_param_double (QSdata *p, int whichparam, double newvalue)   */
/*    int QSget_param (QSdata *p, int whichparam, int *value)               */
/*    int QSget_param_double (QSdata *p, int whichparam, double *value)     */
/*    int QStest_row_norms (QSdata *p)                                      */
/*    int QSopt_strongbranch (QSdata *p, int ncand, int *candidatelist,     */
/*        double *xlist, double *down_vals, double *up_vals,                */
/*        int iterations, double objbound)                                  */
/*    int QSopt_pivotin_row (QSdata *p, int rcnt, int *rlist)               */
/*    int QSopt_pivotin_col (QSdata *p, int ccnt, int *clist)               */
/*    void QSfree (void *ptr)                                               */
/*    void QSstart (void)                                                   */
/*    void QSend (void)                                                     */
/*    char *QSversion (void))                                               */
/*                                                                          */
/*    NEW FUNCTIONS - Add to Docs                                           */
/*                                                                          */
/*    char *QSversion (void))                                               */
/*    int QSget_objsense (QSdata *p, int *objsense)                         */
/*                                                                          */
/****************************************************************************/

#include "econfig.h"
#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"
#include "simplex.h"
#include "price.h"
#include "qstruct.h"
#include "qsopt.h"
#include "lib.h"
#include "mps.h"
#include "lp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
void QSset_precision (const unsigned prec)
{
	EGlpNumSetPrecision (prec);
	ILLchange_precision ();
	/* change the numbers */
}

static void init_basis (QSbasis * B),
  free_cache (QSdata * p);

int grab_cache (QSdata * p,
								int status);
static int opt_work (QSdata * p,
										 int *status,
										 int primal_or_dual),
  qsbasis_to_illbasis (QSbasis * qB,
											 ILLlp_basis * B),
  illbasis_to_qsbasis (ILLlp_basis * B,
											 QSbasis * qB),
  grab_basis (QSdata * p),
  check_qsdata_pointer (QSdata * p);


QSLIB_INTERFACE int QSopt_primal (QSdata * p,
																	int *status)
{
	int rval = 0;

	if (status)
		*status = QS_LP_UNSOLVED;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	/* If both the basis and the cache exist, then skip the optimization */

	if (!p->basis || !p->cache)
	{
		rval = opt_work (p, status, 0);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		if (status)
			*status = p->cache->status;
	}

CLEANUP:

	ILL_RETURN (rval, "QSopt_primal");
}

QSLIB_INTERFACE int QSopt_dual (QSdata * p,
																int *status)
{
	int rval = 0;

	if (status)
		*status = QS_LP_UNSOLVED;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (!p->basis || !p->cache || !p->factorok)
	{
		rval = opt_work (p, status, 1);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		if (status)
			*status = p->cache->status;
	}

CLEANUP:

	ILL_RETURN (rval, "QSopt_dual");
}

static int opt_work (QSdata * p,
										 int *status,
										 int primal_or_dual)
{
	int rval = 0;
	int rstatus = QS_LP_UNSOLVED;
	QSdata *p2 = 0;

	if (p->basis)
	{
		if (p->basis->nstruct != p->qslp->nstruct ||
				p->basis->nrows != p->qslp->nrows)
		{
			fprintf (stderr, "Size of basis does not match LP\n");
			rval = 1;
			goto CLEANUP;
		}
	}

	if (!p->basis && p->lp->basisid == -1 && p->simplex_scaling == 1)
	{
		/* Try scaling by copying the LP and solving */

		ILLprice_free_pricing_info (p->pricing);	/* Just to be sure  */
		p->factorok = 0;						/* that p is clean. */

		p2 = QScopy_prob (p, "scaled_lp");
		if (p2 == 0)
			goto CLEANUP;

		rval = ILLlp_scale (p2->qslp);
		ILL_CLEANUP_IF (rval);

		if (primal_or_dual == 0)
		{
			rval = ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
															PRIMAL_SIMPLEX, 0, p2->simplex_display);
		}
		else
		{
			rval = ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
															DUAL_SIMPLEX, 0, p2->simplex_display);
		}
		ILL_CLEANUP_IF (rval);

		rval = grab_basis (p2);
		ILL_CLEANUP_IF (rval);

		if (p->basis)
		{
			ILLlp_basis_free (p->basis);
			ILL_IFFREE (p->basis, ILLlp_basis);
		}
		p->basis = p2->basis;
		p2->basis = 0;
		QSfree_prob (p2);
		p2 = 0;
	}

	if (primal_or_dual == 0)
	{
		if (p->factorok == 0)
		{
			if (p->basis == 0)
				p->lp->basisid = -1;
			rval = ILLlib_optimize (p->lp, p->basis, p->pricing, PRIMAL_SIMPLEX,
															&rstatus, p->simplex_display);
		}
		else
		{
			ILLprice_free_pricing_info (p->pricing);
			if (p->lp->basisid != -1)
				p->lp->fbasisid = p->lp->basisid;
			rval = ILLlib_optimize (p->lp, 0, p->pricing,
															PRIMAL_SIMPLEX, &rstatus, p->simplex_display);
		}
	}
	else
	{
		if (p->factorok == 0)
		{
			if (p->basis == 0)
				p->lp->basisid = -1;
			rval = ILLlib_optimize (p->lp, p->basis, p->pricing, DUAL_SIMPLEX,
															&rstatus, p->simplex_display);
		}
		else
		{
			/* The factorization and rownorms should be up-to-date */
			if (p->lp->basisid != -1)
			{
				p->lp->fbasisid = p->lp->basisid;
			}
			else
			{
				ILLprice_free_pricing_info (p->pricing);
			}
			rval = ILLlib_optimize (p->lp, 0, p->pricing,
															DUAL_SIMPLEX, &rstatus, p->simplex_display);
		}
	}
	ILL_CLEANUP_IF (rval);

	rval = grab_basis (p);
	ILL_CLEANUP_IF (rval);

	if (rstatus == QS_LP_OPTIMAL)
	{
		rval = grab_cache (p, rstatus);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		free_cache (p);
	}

	p->factorok = 1;

#if 0
	p->lp->basisid = -1;					/* This will cause the basis to be reloaded at the */
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
		QSfree_prob (p2);
	ILL_RETURN (rval, "opt_work");
}

QSLIB_INTERFACE int QSopt_pivotin_row (QSdata * p,
																			 int rcnt,
																			 int *rlist)
{
	int basismod = 0;
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->pricing == 0)
	{
		ILL_ERROR (rval, "pricing info not available in QSopt_pivotin_row\n");
	}

	rval = ILLsimplex_pivotin (p->lp, p->pricing, rcnt, rlist,
														 SIMPLEX_PIVOTINROW, &basismod);
	ILL_CLEANUP_IF (rval);

	rval = grab_basis (p);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSopt_pivotin_row");
}

QSLIB_INTERFACE int QSopt_pivotin_col (QSdata * p,
																			 int ccnt,
																			 int *clist)
{
	int basismod = 0;
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->pricing == 0)
	{
		ILL_ERROR (rval, "pricing info not available in QSopt_pivotin\n");
	}

	rval = ILLsimplex_pivotin (p->lp, p->pricing, ccnt, clist,
														 SIMPLEX_PIVOTINCOL, &basismod);
	ILL_CLEANUP_IF (rval);

	rval = grab_basis (p);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSopt_pivotin_col");
}


QSLIB_INTERFACE int QSopt_strongbranch (QSdata * p,
																				int ncand,
																				int *candidatelist,
																				EGlpNum_t * xlist,
																				EGlpNum_t * down_vals,
																				EGlpNum_t * up_vals,
																				int iterations,
																				EGlpNum_t objbound)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->pricing == 0)
	{
		rval = 1;
		ILL_CLEANUP_IF (rval);
	}

	rval = ILLlib_strongbranch (p->lp, p->pricing, candidatelist, ncand,
															xlist, down_vals, up_vals, iterations, objbound);
	ILL_CLEANUP_IF (rval);

	p->factorok = 0;
	free_cache (p);
	p->qstatus = QS_LP_UNSOLVED;	/* Was set to MODIFIED in free_cache () */

CLEANUP:

	ILL_RETURN (rval, "QSopt_strongbranch");
}

QSLIB_INTERFACE QSdata *QScreate_prob (const char *name,
																			 int objsense)
{
	int rval = 0;
	QSdata *p = 0;
	int len;

	ILL_SAFE_MALLOC (p, 1, QSdata);
	if (!p)
	{
		fprintf (stderr, "out of memory in QScreate_prob\n");
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

	ILL_SAFE_MALLOC (p->qslp, 1, ILLlpdata);
	if (!p->qslp)
	{
		fprintf (stderr, "out of memory in QScreate_prob\n");
		rval = 1;
		goto CLEANUP;
	}
	ILLlpdata_init (p->qslp);

	ILL_SAFE_MALLOC (p->lp, 1, lpinfo);
	if (!p->lp)
	{
		fprintf (stderr, "out of memory in QScreate_prob\n");
		rval = 1;
		goto CLEANUP;
	}
	EGlpNumInitVar (p->lp->objval);
	EGlpNumInitVar (p->lp->pobjval);
	EGlpNumInitVar (p->lp->dobjval);
	EGlpNumInitVar (p->lp->pinfeas);
	EGlpNumInitVar (p->lp->dinfeas);
	EGlpNumInitVar (p->lp->objbound);
	EGlpNumInitVar (p->lp->upd.piv);
	EGlpNumInitVar (p->lp->upd.dty);
	EGlpNumInitVar (p->lp->upd.c_obj);
	EGlpNumInitVar (p->lp->upd.tz);
	ILLsimplex_init_lpinfo (p->lp);
	ILLsimplex_load_lpinfo (p->qslp, p->lp);

	ILL_SAFE_MALLOC (p->pricing, 1, price_info);
	if (!p->pricing)
	{
		fprintf (stderr, "out of memory in QScreate_prob\n");
		rval = 1;
		goto CLEANUP;
	}
	EGlpNumInitVar (p->pricing->htrigger);
	ILLprice_init_pricing_info (p->pricing);
	p->pricing->pI_price = QS_DEFAULT_PRICE_PI;
	p->pricing->pII_price = QS_DEFAULT_PRICE_PII;
	p->pricing->dI_price = QS_DEFAULT_PRICE_DI;
	p->pricing->dII_price = QS_DEFAULT_PRICE_DII;

	if (name)
	{
		len = strlen (name) + 1;
		ILL_SAFE_MALLOC (p->name, len, char);
		strcpy (p->name, name);
	}
	else
	{
		ILL_SAFE_MALLOC (p->name, 7, char);
		sprintf (p->name, "noname");
	}

	len = strlen (p->name) + 1;
	ILL_SAFE_MALLOC (p->qslp->probname, len, char);
	strcpy (p->qslp->probname, p->name);

	if (objsense == QS_MAX)
	{
		p->qslp->objsense = QS_MAX;
	}

CLEANUP:

	if (rval)
	{
		QSfree_prob (p);
		p = 0;
	}

	return p;
}

QSLIB_INTERFACE QSdata *QSread_prob (const char *filename,
																		 const char *filetype)
{
	QSdata *p = 0;
	FILE *file = 0;
	QSline_reader reader;

	if ((file = fopen (filename, "r")) == 0)
	{
		perror (filename);
		fprintf (stderr, "Unable to open \"%s\" for input.\n", filename);
	}
	if (file == NULL)
		goto CLEANUP;

	reader = ILLline_reader_new ((qsread_line_fct) fgets, file);
	p = QSget_prob (reader, filename, filetype);
	QSline_reader_free (reader);	/* Bico - 040723 */

CLEANUP:
	if (file != NULL)
	{
		fclose (file);
	}
	return p;
}

QSLIB_INTERFACE QSdata *QSload_prob (const char *probname,
																		 int ncols,
																		 int nrows,
																		 int *cmatcnt,
																		 int *cmatbeg,
																		 int *cmatind,
																		 EGlpNum_t * cmatval,
																		 int objsense,
																		 EGlpNum_t * obj,
																		 EGlpNum_t * rhs,
																		 char *sense,
																		 EGlpNum_t * lower,
																		 EGlpNum_t * upper,
																		 const char **colnames,
																		 const char **rownames)
{
	int rval = 0;
	QSdata *p = 0;

	p = QScreate_prob (probname, objsense);
	if (p == 0)
		goto CLEANUP;

	rval = ILLlib_newrows (p->lp, 0, nrows, rhs, sense, 0, rownames);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_addcols (p->lp, 0, ncols, cmatcnt, cmatbeg, cmatind,
												 cmatval, obj, lower, upper, colnames, 0);
	ILL_CLEANUP_IF (rval);

	p->factorok = 0;

CLEANUP:

	if (rval)
	{
		QSfree_prob (p);
		p = 0;
	}

	return p;
}

QSLIB_INTERFACE QSdata *QScopy_prob (QSdata * p,
																		 const char *newname)
{
	int rval = 0;
	int j,
	  col,
	  beg,
	  pindex,
	  hit;
	QSdata *p2 = 0;
	char *coln;
	char buf[ILL_namebufsize];

	/* printf ("QScopy_prob ...\n"); fflush (stdout); */

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	p2 = QScreate_prob (newname, p->qslp->objsense);
	if (p2 == 0)
		goto CLEANUP;

	rval = ILLlib_newrows (p2->lp, 0, p->qslp->nrows,
												 p->qslp->rhs, p->qslp->sense, p->qslp->rangeval,
												 (const char **) p->qslp->rownames);
	ILL_CLEANUP_IF (rval);

	for (j = 0; j < p->qslp->nstruct; j++)
	{
		col = p->qslp->structmap[j];
		if (p->qslp->colnames)
			coln = p->qslp->colnames[j];
		else
			coln = 0;
		beg = p->qslp->A.matbeg[col];

		/* Monika: Note that Java will need to handle these arrays */
		/*         without using the beg offset.  The easiest way  */
		/*         may be to copy the arrays, as in the getcols()  */
		/*         code in lib.c.                                  */

		rval = ILLlib_addcol (p2->lp, 0,
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
	EGlpNumClearVar (p2->pricing->htrigger);
	*(p2->pricing) = *(p->pricing);
	/* I added this line because copying the heap (as a pointer) doesn't make any
	 * sense ! */
	ILLheap_init (&(p2->pricing->h));
	EGlpNumInitVar (p2->pricing->htrigger);
	EGlpNumCopy (p2->pricing->htrigger, p->pricing->htrigger);

	if (p->qslp->intmarker != 0)
	{
		ILL_SAFE_MALLOC (p2->qslp->intmarker, p->qslp->nstruct, char);
		for (j = 0; j < p->qslp->nstruct; j++)
		{
			p2->qslp->intmarker[j] = p->qslp->intmarker[j];
		}
	}

	if (p->qslp->objname != 0)
	{
		ILL_UTIL_STR (p2->qslp->objname, p->qslp->objname);
	}
	else
	{
		strcpy (buf, "obj");
		rval = ILLsymboltab_uname (&p2->qslp->rowtab, buf, "", NULL);
		ILL_CLEANUP_IF (rval);
		ILL_UTIL_STR (p2->qslp->objname, buf);
	}
	rval = ILLsymboltab_register (&p2->qslp->rowtab, p2->qslp->objname,
																-1, &pindex, &hit);
	rval = rval || hit;
	ILL_CLEANUP_IF (rval);

	ILLstring_reporter_copy (&p2->qslp->reporter, &p->qslp->reporter);

CLEANUP:

	if (rval)
	{
		QSfree_prob (p2);
		p2 = 0;
	}

	return p2;
}

QSLIB_INTERFACE int QSchange_objsense (QSdata * p,
																			 int newsense)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (newsense != QS_MIN && newsense != QS_MAX)
	{
		fprintf (stderr, "Illegal objective sense %d\n", newsense);
		rval = 1;
		goto CLEANUP;
	}

	if (p->qslp->objsense != newsense)
	{
		p->qslp->objsense = newsense;
		free_cache (p);
	}

CLEANUP:

	ILL_RETURN (rval, "QSchange_objsense");
}

QSLIB_INTERFACE int QSget_objsense (QSdata * p,
																		int *objsense)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (objsense)
		*objsense = p->qslp->objsense;

CLEANUP:

	ILL_RETURN (rval, "QSget_objsense");
}


QSLIB_INTERFACE int QSnew_col (QSdata * p,
															 EGlpNum_t obj,
															 EGlpNum_t lower,
															 EGlpNum_t upper,
															 const char *name)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_newcol (p->lp, p->basis, obj, lower, upper, name, p->factorok);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSnew_col");
}

QSLIB_INTERFACE int QSadd_cols (QSdata * p,
																int num,
																int *cmatcnt,
																int *cmatbeg,
																int *cmatind,
																EGlpNum_t * cmatval,
																EGlpNum_t * obj,
																EGlpNum_t * lower,
																EGlpNum_t * upper,
																const char **names)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_addcols (p->lp, p->basis, num, cmatcnt, cmatbeg,
												 cmatind, cmatval, obj, lower, upper, names,
												 p->factorok);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSadd_cols");
}

QSLIB_INTERFACE int QSadd_col (QSdata * p,
															 int cnt,
															 int *cmatind,
															 EGlpNum_t * cmatval,
															 EGlpNum_t obj,
															 EGlpNum_t lower,
															 EGlpNum_t upper,
															 const char *name)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_addcol (p->lp, p->basis, cnt, cmatind, cmatval,
												obj, lower, upper, name, p->factorok);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSadd_col");
}

QSLIB_INTERFACE int QSnew_row (QSdata * p,
															 EGlpNum_t rhs,
															 int sense,
															 const char *name)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_newrow (p->lp, p->basis, rhs, sense, zeroLpNum, name);
	ILL_CLEANUP_IF (rval);

	p->factorok = 0;
	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSnew_row");
}

QSLIB_INTERFACE int QSadd_rows (QSdata * p,
																int num,
																int *rmatcnt,
																int *rmatbeg,
																int *rmatind,
																EGlpNum_t * rmatval,
																EGlpNum_t * rhs,
																char *sense,
																const char **names)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_addrows (p->lp, p->basis, num, rmatcnt, rmatbeg,
												 rmatind, rmatval, rhs, sense, 0,
												 names, &(p->factorok));
	ILL_CLEANUP_IF (rval);

	if (p->factorok == 1 && p->basis->rownorms)
	{
		rval = ILLlib_loadrownorms (p->lp, p->pricing, p->basis->rownorms);
		ILL_CLEANUP_IF (rval);
		/* This really should go inside of ILLlib_addrows, once pinf is  */
		/* is moved into the lp struct.                                  */
	}

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSadd_rows");
}

QSLIB_INTERFACE int QSadd_row (QSdata * p,
															 int cnt,
															 int *rmatind,
															 EGlpNum_t * rmatval,
															 EGlpNum_t * rhs,
															 int sense,
															 const char *name)
{
	int rval = 0;
	int vmatcnt[1];
	int vmatbeg[1];
	char vsense[1];
	const char *vnames[1];

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	vmatcnt[0] = cnt;
	vmatbeg[0] = 0;
	vsense[0] = sense;
	vnames[0] = name;

	rval = QSadd_rows (p, 1, vmatcnt, vmatbeg, rmatind, rmatval, rhs, vsense,
										 vnames);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSadd_row");
}

QSLIB_INTERFACE int QSdelete_rows (QSdata * p,
																	 int num,
																	 int *dellist)
{
	int rval = 0;
	int basis_ok = 0;
	int cache_ok = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_delrows (p->lp, p->basis, p->cache, num, dellist, &basis_ok,
												 &cache_ok);
	ILL_CLEANUP_IF (rval);

	/* For now, just remove the basis - wait for pivotin */

	if (p->basis && !basis_ok)
	{
		ILLlp_basis_free (p->basis);
		ILL_IFFREE (p->basis, ILLlp_basis);
	}

	p->factorok = 0;

	if (!p->basis || !basis_ok || !cache_ok)
	{
		/* Note: If we only delete basic rows then cached soln is valid */
		free_cache (p);
	}

CLEANUP:

	ILL_RETURN (rval, "QSdelete_rows");
}

QSLIB_INTERFACE int QSdelete_row (QSdata * p,
																	int rowindex)
{
	int rval = 0;
	int vdellist[1];

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	vdellist[0] = rowindex;

	rval = QSdelete_rows (p, 1, vdellist);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSdelete_row");
}

QSLIB_INTERFACE int QSdelete_setrows (QSdata * p,
																			int *flags)
{
	int rval = 0;
	int j,
	  num = 0;
	int *dellist = 0;
	int nrows;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	nrows = p->qslp->nrows;

	for (j = 0; j < nrows; j++)
	{
		if (flags[j] == 1)
			num++;
	}

	if (num > 0)
	{
		ILL_SAFE_MALLOC (dellist, num, int);

		for (j = 0, num = 0; j < nrows; j++)
		{
			if (flags[j] == 1)
			{
				dellist[num++] = j;
			}
		}

		rval = QSdelete_rows (p, num, dellist);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	ILL_IFFREE (dellist, int);
	ILL_RETURN (rval, "QSdelete_setrows");
}

QSLIB_INTERFACE int QSdelete_named_row (QSdata * p,
																				const char *rowname)
{
	int rval = 0;
	int i,
	  vdellist[1];

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = QSget_row_index (p, rowname, &i);
	ILL_CLEANUP_IF (rval);

	vdellist[0] = i;

	rval = QSdelete_rows (p, 1, vdellist);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSdelete_named_row");
}

QSLIB_INTERFACE int QSdelete_named_rows_list (QSdata * p,
																							int num,
																							const char **rownames)
{
	int rval = 0;
	int i,
	  k;
	int *vdellist = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (num > 0)
	{
		ILL_SAFE_MALLOC (vdellist, num, int);
		for (k = 0; k < num; k++)
		{
			rval = QSget_row_index (p, rownames[k], &i);
			ILL_CLEANUP_IF (rval);
			vdellist[k] = i;
		}

		rval = QSdelete_rows (p, num, vdellist);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	ILL_IFFREE (vdellist, int);
	ILL_RETURN (rval, "QSdelete_named_rows_list");
}

QSLIB_INTERFACE int QSdelete_cols (QSdata * p,
																	 int num,
																	 int *dellist)
{
	int rval = 0;
	int basis_ok;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_delcols (p->lp, p->basis, num, dellist, &basis_ok);
	ILL_CLEANUP_IF (rval);

	/* For now, just remove the basis - wait for pivotout */

	if (p->basis && !basis_ok)
	{
		ILLlp_basis_free (p->basis);
		ILL_IFFREE (p->basis, ILLlp_basis);
	}

	p->factorok = 0;
	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSdelete_cols");
}

QSLIB_INTERFACE int QSdelete_col (QSdata * p,
																	int colindex)
{
	int rval = 0;
	int vdellist[1];

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	vdellist[0] = colindex;

	rval = QSdelete_cols (p, 1, vdellist);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSdelete_col");
}

QSLIB_INTERFACE int QSdelete_setcols (QSdata * p,
																			int *flags)
{
	int rval = 0;
	int j,
	  num = 0;
	int *dellist = 0;
	int ncols;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	ncols = p->qslp->nstruct;

	for (j = 0; j < ncols; j++)
	{
		if (flags[j] == 1)
			num++;
	}

	if (num > 0)
	{
		ILL_SAFE_MALLOC (dellist, num, int);

		for (j = 0, num = 0; j < ncols; j++)
		{
			if (flags[j] == 1)
			{
				dellist[num++] = j;
			}
		}

		rval = QSdelete_cols (p, num, dellist);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	ILL_IFFREE (dellist, int);
	ILL_RETURN (rval, "QSdelete_setcols");
}

QSLIB_INTERFACE int QSdelete_named_column (QSdata * p,
																					 const char *colname)
{
	int rval = 0;
	int j,
	  vdellist[1];

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = QSget_column_index (p, colname, &j);
	ILL_CLEANUP_IF (rval);

	vdellist[0] = j;

	rval = QSdelete_cols (p, 1, vdellist);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSdelete_named_column");
}

QSLIB_INTERFACE int QSdelete_named_columns_list (QSdata * p,
																								 int num,
																								 const char **colnames)
{
	int rval = 0;
	int i,
	  j;
	int *vdellist = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (num > 0)
	{
		ILL_SAFE_MALLOC (vdellist, num, int);
		for (i = 0; i < num; i++)
		{
			rval = QSget_column_index (p, colnames[i], &j);
			ILL_CLEANUP_IF (rval);
			vdellist[i] = j;
		}

		rval = QSdelete_cols (p, num, vdellist);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	ILL_IFFREE (vdellist, int);
	ILL_RETURN (rval, "QSdelete_named_columns_list");
}

QSLIB_INTERFACE int QSchange_senses (QSdata * p,
																		 int num,
																		 int *rowlist,
																		 char *sense)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_chgsense (p->lp, num, rowlist, sense);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSchange_senses");
}

QSLIB_INTERFACE int QSchange_sense (QSdata * p,
																		int rowindex,
																		int sense)
{
	int rval = 0;
	int vrowlist[1];
	char vsenselist[1];

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	vrowlist[0] = rowindex;
	vsenselist[0] = sense;

	rval = QSchange_senses (p, 1, vrowlist, vsenselist);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSchange_sense");
}

QSLIB_INTERFACE int QSchange_coef (QSdata * p,
																	 int rowindex,
																	 int colindex,
																	 EGlpNum_t coef)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_chgcoef (p->lp, rowindex, colindex, coef);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSchange_coef");
}

QSLIB_INTERFACE int QSchange_objcoef (QSdata * p,
																			int indx,
																			EGlpNum_t coef)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_chgobj (p->lp, indx, coef);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSchange_objcoef");
}

QSLIB_INTERFACE int QSchange_rhscoef (QSdata * p,
																			int indx,
																			EGlpNum_t coef)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_chgrhs (p->lp, indx, coef);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSchange_rhscoef");
}

QSLIB_INTERFACE int QSchange_bounds (QSdata * p,
																		 int num,
																		 int *collist,
																		 char *lu,
																		 EGlpNum_t * bounds)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_chgbnds (p->lp, num, collist, lu, bounds);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSchange_bounds");
}

QSLIB_INTERFACE int QSchange_bound (QSdata * p,
																		int indx,
																		int lu,
																		EGlpNum_t bound)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_chgbnd (p->lp, indx, lu, bound);
	ILL_CLEANUP_IF (rval);

	free_cache (p);

CLEANUP:

	ILL_RETURN (rval, "QSchange_bound");
}

#if 0
		/*
		 * Bico - I removed this on 02.04.22.  I don't think we need to support
		 * this type of interface (the loading via arrays can do the job)
		 */
QSLIB_INTERFACE QSbasis *QScreate_basis (int nstruct,
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

	if (nstruct)
	{
		ILL_SAFE_MALLOC (B->cstat, nstruct, char);
	}

	if (nrows)
	{
		ILL_SAFE_MALLOC (B->rstat, nrows, char);
	}

	for (i = 0; i < nstruct; i++)
		B->cstat[i] = 0;
	for (i = 0; i < nrows; i++)
		B->rstat[i] = 0;

CLEANUP:

	if (rval)
	{
		QSfree_basis (B);
		B = 0;
	}

	return B;
}
#endif

QSLIB_INTERFACE QSbasis *QSread_basis (QSdata * p,
																			 const char *filename)
{
	int rval = 0;
	QSbasis *qB = 0;
	ILLlp_basis B;

	ILLlp_basis_init (&B);

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	ILL_NEW (qB, QSbasis);
	init_basis (qB);

	rval = ILLlib_readbasis (p->lp, &B, filename);
	ILL_CLEANUP_IF (rval);

	rval = illbasis_to_qsbasis (&B, qB);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	if (rval && qB)
	{
		QSfree_basis (qB);
		qB = 0;
	}
	ILLlp_basis_free (&B);

	return qB;
}

QSLIB_INTERFACE int QSload_basis (QSdata * p,
																	QSbasis * B)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (B->nstruct != p->qslp->nstruct || B->nrows != p->qslp->nrows)
	{
		fprintf (stderr, "size of basis does not match lp\n");
		rval = 1;
		goto CLEANUP;
	}

	if (p->basis == 0)
	{
		ILL_SAFE_MALLOC (p->basis, 1, ILLlp_basis);
		ILLlp_basis_init (p->basis);
	}
	else
	{
		ILLlp_basis_free (p->basis);
	}

	rval = qsbasis_to_illbasis (B, p->basis);
	ILL_CLEANUP_IF (rval);

	p->factorok = 0;

CLEANUP:

	ILL_RETURN (rval, "QSload_basis");
}

QSLIB_INTERFACE int QSread_and_load_basis (QSdata * p,
																					 const char *filename)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->basis == 0)
	{
		ILL_SAFE_MALLOC (p->basis, 1, ILLlp_basis);
		ILLlp_basis_init (p->basis);
	}
	else
	{
		ILLlp_basis_free (p->basis);
	}

	rval = ILLlib_readbasis (p->lp, p->basis, filename);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	return rval;
}

QSLIB_INTERFACE int QSload_basis_array (QSdata * p,
																				char *cstat,
																				char *rstat)
{
	int rval = 0;
	int i;
	ILLlp_basis *B;
	ILLlpdata *qslp;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	qslp = p->qslp;

	if (qslp->nstruct > 0 && cstat == 0)
	{
		fprintf (stderr, "QSload_basis_array called without cstat\n");
		rval = 1;
		goto CLEANUP;
	}

	if (qslp->nrows > 0 && rstat == 0)
	{
		fprintf (stderr, "QSload_basis_array called without rstat\n");
		rval = 1;
		goto CLEANUP;
	}

	if (p->basis == 0)
	{
		ILL_SAFE_MALLOC (p->basis, 1, ILLlp_basis);
		ILLlp_basis_init (p->basis);
	}
	else
	{
		ILLlp_basis_free (p->basis);
	}

	B = p->basis;

	B->nstruct = qslp->nstruct;
	B->nrows = qslp->nrows;
	ILL_SAFE_MALLOC (B->cstat, qslp->nstruct, char);
	ILL_SAFE_MALLOC (B->rstat, qslp->nrows, char);

	for (i = 0; i < qslp->nstruct; i++)
	{
		B->cstat[i] = cstat[i];
	}

	for (i = 0; i < qslp->nrows; i++)
	{
		B->rstat[i] = rstat[i];
	}

	p->factorok = 0;

CLEANUP:

	ILL_RETURN (rval, "QSload_basis_array");
}

QSLIB_INTERFACE int QSload_basis_and_row_norms_array (QSdata * p,
																											char *cstat,
																											char *rstat,
																											EGlpNum_t * rownorms)
{
	int rval = 0;
	int i,
	  nrows;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	nrows = p->qslp->nrows;

	rval = QSload_basis_array (p, cstat, rstat);
	ILL_CLEANUP_IF (rval);
	p->basis->rownorms = EGlpNumAllocArray (nrows);

	for (i = 0; i < nrows; i++)
	{
		EGlpNumCopy (p->basis->rownorms[i], rownorms[i]);
	}

	p->factorok = 0;

CLEANUP:

	ILL_RETURN (rval, "QSload_basis_and_row_norms_array");
}

QSLIB_INTERFACE int QSwrite_basis (QSdata * p,
																	 QSbasis * B,
																	 const char *filename)
{
	int rval = 0;
	ILLlp_basis iB,
	 *basis = 0;

	ILLlp_basis_init (&iB);

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (B)
	{
		rval = qsbasis_to_illbasis (B, &iB);
		ILL_CLEANUP_IF (rval);
		basis = &iB;
	}
	else
	{
		if (p->basis == 0)
		{
			fprintf (stderr, "no basis available in QSwrite_basis\n");
			rval = 1;
			goto CLEANUP;
		}
		basis = p->basis;
	}

	rval = ILLlib_writebasis (p->lp, basis, filename);
	ILL_CLEANUP_IF (rval);


CLEANUP:

	ILLlp_basis_free (basis);
	ILL_RETURN (rval, "QSwrite_basis");
}

QSLIB_INTERFACE QSbasis *QSget_basis (QSdata * p)
{
	int rval = 0;
	QSbasis *B = 0;

	if (p->basis == 0)
	{
		fprintf (stderr, "no basis available in QSget_basis\n");
		rval = 1;
		goto CLEANUP;
	}

	ILL_SAFE_MALLOC (B, 1, QSbasis);
	init_basis (B);
	rval = illbasis_to_qsbasis (p->basis, B);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	if (rval)
	{
		QSfree_basis (B);
		B = 0;
	}

	return B;
}

QSLIB_INTERFACE int QSget_basis_array (QSdata * p,
																			 char *cstat,
																			 char *rstat)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->basis == 0)
	{
		fprintf (stderr, "no basis available in QSget_basis_array\n");
		rval = 1;
		goto CLEANUP;
	}

	for (i = 0; i < p->basis->nstruct; i++)
		cstat[i] = p->basis->cstat[i];
	for (i = 0; i < p->basis->nrows; i++)
		rstat[i] = p->basis->rstat[i];

CLEANUP:

	ILL_RETURN (rval, "QSget_basis_array");
}

QSLIB_INTERFACE int QSget_basis_and_row_norms_array (QSdata * p,
																										 char *cstat,
																										 char *rstat,
																										 EGlpNum_t * rownorms)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->basis == 0)
	{
		fprintf (stderr, "no basis available\n");
		rval = 1;
		goto CLEANUP;
	}

	for (i = 0; i < p->basis->nstruct; i++)
		cstat[i] = p->basis->cstat[i];
	for (i = 0; i < p->basis->nrows; i++)
		rstat[i] = p->basis->rstat[i];

	if (p->basis->rownorms == 0)
	{
		fprintf (stderr, "no row norms available\n");
		rval = 1;
		goto CLEANUP;
	}
	else
	{
		for (i = 0; i < p->basis->nrows; i++)
		{
			EGlpNumCopy (rownorms[i], p->basis->rownorms[i]);
		}
	}

CLEANUP:

	ILL_RETURN (rval, "QSget_basis_and_row_norms_array");
}

static int illbasis_to_qsbasis (ILLlp_basis * B,
																QSbasis * qB)
{
	int rval = 0;
	int i;

	qB->nstruct = B->nstruct;
	qB->nrows = B->nrows;
	ILL_SAFE_MALLOC (qB->cstat, B->nstruct, char);
	ILL_SAFE_MALLOC (qB->rstat, B->nrows, char);

	for (i = 0; i < B->nstruct; i++)
	{
		qB->cstat[i] = B->cstat[i];
	}

	for (i = 0; i < B->nrows; i++)
	{
		qB->rstat[i] = B->rstat[i];
	}

CLEANUP:

	ILL_RETURN (rval, "illbasis_to_qsbasis");
}

static int qsbasis_to_illbasis (QSbasis * qB,
																ILLlp_basis * B)
{
	int rval = 0;
	int i;

	B->nstruct = qB->nstruct;
	B->nrows = qB->nrows;
	ILL_SAFE_MALLOC (B->cstat, qB->nstruct, char);
	ILL_SAFE_MALLOC (B->rstat, qB->nrows, char);

	for (i = 0; i < qB->nstruct; i++)
	{
		B->cstat[i] = qB->cstat[i];
	}

	for (i = 0; i < qB->nrows; i++)
	{
		B->rstat[i] = qB->rstat[i];
	}

CLEANUP:

	ILL_RETURN (rval, "qsbasis_to_illbasis");
}

static int grab_basis (QSdata * p)
{
	int rval = 0;
	ILLlp_basis *B = p->basis;
	int nstruct = p->qslp->nstruct;
	int nrows = p->qslp->nrows;

	if (!B)
	{
		ILL_SAFE_MALLOC (p->basis, 1, ILLlp_basis);
		ILLlp_basis_init (p->basis);
		B = p->basis;
	}

	if (nstruct != B->nstruct)
	{
		ILL_IFFREE (B->cstat, char);
		ILL_SAFE_MALLOC (B->cstat, nstruct, char);
		B->nstruct = nstruct;
	}

	if (nrows != B->nrows)
	{
		ILL_IFFREE (B->rstat, char);
		ILL_SAFE_MALLOC (B->rstat, nrows, char);
		B->nrows = nrows;
	}

	rval = ILLlib_getbasis (p->lp, B->cstat, B->rstat);
	ILL_CLEANUP_IF (rval);

	EGlpNumFreeArray (B->rownorms);
	EGlpNumFreeArray (B->colnorms);

	if (p->pricing->dII_price == QS_PRICE_DSTEEP)
	{
		B->rownorms = EGlpNumAllocArray (nrows);
		rval = ILLlib_getrownorms (p->lp, p->pricing, B->rownorms);
		if (rval)
		{
/*
            fprintf (stderr, "no edge norms, continue anyway\n");
*/
			EGlpNumFreeArray (B->rownorms);
			rval = 0;
		}
	}

CLEANUP:

	if (rval)
	{
		if (B)
		{
			ILLlp_basis_free (B);
			ILL_IFFREE (p->basis, ILLlp_basis);
		}
	}

	ILL_RETURN (rval, "grab_basis");
}

int grab_cache (QSdata * p,
								int status)
{
	int rval = 0;
	ILLlp_cache *C = p->cache;
	int nstruct = p->qslp->nstruct;
	int nrows = p->qslp->nrows;

	if (C == 0)
	{
		ILL_SAFE_MALLOC (p->cache, 1, ILLlp_cache);
		EGlpNumInitVar (p->cache->val);
		ILLlp_cache_init (p->cache);
		C = p->cache;
	}

	if (nstruct != C->nstruct || nrows != C->nrows)
	{
		ILLlp_cache_free (C);
		rval = ILLlp_cache_alloc (C, nstruct, nrows);
		ILL_CLEANUP_IF (rval);
	}

	rval = ILLlib_cache_solution (p->lp, C);
	ILL_CLEANUP_IF (rval);

	C->status = status;

CLEANUP:

	if (rval)
	{
		if (C)
		{
			ILLlp_cache_free (C);
			EGlpNumClearVar (p->cache->val);
			ILL_IFFREE (p->cache, ILLlp_cache);
		}
	}

	ILL_RETURN (rval, "grab_cache");
}

void free_cache (QSdata * p)
{
	if (p->cache)
	{
		ILLlp_cache_free (p->cache);
		EGlpNumClearVar (p->cache->val);
		ILL_IFFREE (p->cache, ILLlp_cache);
	}
	p->qstatus = QS_LP_MODIFIED;
}

QSLIB_INTERFACE int QSget_binv_row (QSdata * p,
																		int indx,
																		EGlpNum_t * binvrow)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "LP has not been optimized in QSget_binv_row\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_tableau (p->lp, indx, binvrow, 0);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_binv_row");
}

QSLIB_INTERFACE int QSget_tableau_row (QSdata * p,
																			 int indx,
																			 EGlpNum_t * tableaurow)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "LP has not been optimized in QSget_tableau_row\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_tableau (p->lp, indx, 0, tableaurow);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_tableau_row");
}

QSLIB_INTERFACE int QSget_basis_order (QSdata * p,
																			 int *basorder)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "LP has not been optimized in QSget_basis_order\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_basis_order (p->lp, basorder);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_basis_order");
}

QSLIB_INTERFACE int QScompute_row_norms (QSdata * p)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->pricing->dII_price != QS_PRICE_DSTEEP)
	{
		fprintf (stderr, "not using dual steepest edge\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_recompute_rownorms (p->lp, p->pricing);
	ILL_CLEANUP_IF (rval);

	rval = grab_basis (p);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QScompute_row_norms");
}

QSLIB_INTERFACE void QSfree_prob (QSdata * p)
{
	if (p)
	{
		if (p->qslp)
		{
			ILLlpdata_free (p->qslp);
			ILL_IFFREE (p->qslp, ILLlpdata);
		}
		if (p->lp)
		{
			ILLsimplex_free_lpinfo (p->lp);
			EGlpNumClearVar (p->lp->objval);
			EGlpNumClearVar (p->lp->pobjval);
			EGlpNumClearVar (p->lp->dobjval);
			EGlpNumClearVar (p->lp->pinfeas);
			EGlpNumClearVar (p->lp->dinfeas);
			EGlpNumClearVar (p->lp->objbound);
			EGlpNumClearVar (p->lp->upd.piv);
			EGlpNumClearVar (p->lp->upd.dty);
			EGlpNumClearVar (p->lp->upd.c_obj);
			EGlpNumClearVar (p->lp->upd.tz);
			ILL_IFFREE (p->lp, lpinfo);
		}
		if (p->basis)
		{
			ILLlp_basis_free (p->basis);
			ILL_IFFREE (p->basis, ILLlp_basis);
		}
		if (p->cache)
		{
			ILLlp_cache_free (p->cache);
			EGlpNumClearVar (p->cache->val);
			ILL_IFFREE (p->cache, ILLlp_cache);
		}
		if (p->pricing)
		{
			EGlpNumClearVar (p->pricing->htrigger);
			ILLprice_free_pricing_info (p->pricing);
			ILL_IFFREE (p->pricing, price_info);
		}
		ILL_IFFREE (p->name, char);
		ILL_IFFREE (p, QSdata);
	}
}

QSLIB_INTERFACE void QSfree_basis (QSbasis * B)
{
	if (B)
	{
		ILL_IFFREE (B->rstat, char);
		ILL_IFFREE (B->cstat, char);
		ILL_IFFREE (B, QSbasis);
	}
}

static void init_basis (QSbasis * B)
{
	if (B)
	{
		B->nstruct = 0;
		B->nrows = 0;
		B->cstat = 0;
		B->rstat = 0;
	}
}

QSLIB_INTERFACE int QSget_status (QSdata * p,
																	int *status)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (status)
		*status = p->qstatus;

CLEANUP:

	ILL_RETURN (rval, "QSget_status");
}

QSLIB_INTERFACE int QSget_solution (QSdata * p,
																		EGlpNum_t * value,
																		EGlpNum_t * x,
																		EGlpNum_t * pi,
																		EGlpNum_t * slack,
																		EGlpNum_t * rc)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_solution\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_solution (p->lp, p->cache, value, x, pi, slack, rc);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_solution");
}

QSLIB_INTERFACE int QSget_objval (QSdata * p,
																	EGlpNum_t * value)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

/* Want to get objval after limited number of pivots. */

	if (p->qstatus == QS_LP_MODIFIED)
	{
		fprintf (stderr, "QSmsg: LP has been modified since last solve.\n");
		rval = 1;
		ILL_CLEANUP;
	}

	rval = ILLlib_objval (p->lp, p->cache, value);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_objval");
}

QSLIB_INTERFACE int QSget_x_array (QSdata * p,
																	 EGlpNum_t * x)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_x_array\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_get_x (p->lp, p->cache, x);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_x_array");
}

QSLIB_INTERFACE int QSget_slack_array (QSdata * p,
																			 EGlpNum_t * slack)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_slack_array\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_get_slack (p->lp, p->cache, slack);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_slack_array");
}

QSLIB_INTERFACE int QSget_rc_array (QSdata * p,
																		EGlpNum_t * rc)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_rc_array\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_solution (p->lp, p->cache, 0, 0, 0, 0, rc);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_rc_array");
}

QSLIB_INTERFACE int QSget_pi_array (QSdata * p,
																		EGlpNum_t * pi)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_pi_array\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = ILLlib_solution (p->lp, p->cache, 0, 0, pi, 0, 0);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_pi_array");
}

QSLIB_INTERFACE int QSget_infeas_array (QSdata * p,
																				EGlpNum_t * pi)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (pi == 0)
	{
		ILL_ERROR (rval, "QS_get_infeas_array called with NULL pi vector\n");
	}

	rval = ILLsimplex_infcertificate (p->lp, pi);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_infeas_array");
}

QSLIB_INTERFACE int QSget_named_x (QSdata * p,
																	 const char *colname,
																	 EGlpNum_t * val)
{
	int rval = 0;
	int j;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_named_x\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = QSget_column_index (p, colname, &j);
	ILL_CLEANUP_IF (rval);

	if (j != -1)
	{
		EGlpNumCopy (*val, p->cache->x[j]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	ILL_RETURN (rval, "QSget_named_x");
}

QSLIB_INTERFACE int QSget_named_rc (QSdata * p,
																		const char *colname,
																		EGlpNum_t * val)
{
	int rval = 0;
	int j;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_named_rc\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = QSget_column_index (p, colname, &j);
	ILL_CLEANUP_IF (rval);

	if (j != -1)
	{
		EGlpNumCopy (*val, p->cache->rc[j]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	ILL_RETURN (rval, "QSget_named_rc");
}

QSLIB_INTERFACE int QSget_named_pi (QSdata * p,
																		const char *rowname,
																		EGlpNum_t * val)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_named_pi\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = QSget_row_index (p, rowname, &i);
	ILL_CLEANUP_IF (rval);

	if (i != -1)
	{
		EGlpNumCopy (*val, p->cache->pi[i]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	ILL_RETURN (rval, "QSget_named_pi");
}

QSLIB_INTERFACE int QSget_named_slack (QSdata * p,
																			 const char *rowname,
																			 EGlpNum_t * val)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->cache == 0)
	{
		fprintf (stderr, "no solution available in QSget_named_slack\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = QSget_row_index (p, rowname, &i);
	ILL_CLEANUP_IF (rval);

	if (i != -1)
	{
		EGlpNumCopy (*val, p->cache->slack[i]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	ILL_RETURN (rval, "QSget_named_slack");
}

QSLIB_INTERFACE int QSget_colcount (QSdata * p)
{
	int cnt;

	if (check_qsdata_pointer (p))
		cnt = 0;
	else
		cnt = p->qslp->nstruct;

	return cnt;
}

QSLIB_INTERFACE int QSget_rowcount (QSdata * p)
{
	int cnt;

	if (check_qsdata_pointer (p))
		cnt = 0;
	else
		cnt = p->qslp->nrows;

	return cnt;
}

QSLIB_INTERFACE int QSget_nzcount (QSdata * p)
{
	int cnt;

	if (check_qsdata_pointer (p))
		cnt = 0;
	else
		cnt = p->qslp->nzcount - p->qslp->nrows;

	return cnt;
}

QSLIB_INTERFACE int QStest_row_norms (QSdata * p)
{
	int yesno;

	if (check_qsdata_pointer (p))
	{
		yesno = 0;
	}
	else
	{
		if (p->basis && p->basis->rownorms)
		{
			yesno = 1;
		}
		else
		{
			yesno = 0;
		}
	}

	return yesno;
}

QSLIB_INTERFACE int QSget_obj (QSdata * p,
															 EGlpNum_t * obj)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_getobj (p->lp, obj);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_obj");
}

QSLIB_INTERFACE int QSget_rhs (QSdata * p,
															 EGlpNum_t * rhs)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_getrhs (p->lp, rhs);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_rhs");
}

QSLIB_INTERFACE int QSget_rows_list (QSdata * p,
																		 int num,
																		 int *rowlist,
																		 int **rowcnt,
																		 int **rowbeg,
																		 int **rowind,
																		 EGlpNum_t ** rowval,
																		 EGlpNum_t ** rhs,
																		 char **sense,
																		 char ***names)
{
	int rval = 0;
	int i,
	  nrows;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	nrows = QSget_rowcount (p);
	for (i = 0; i < num; i++)
	{
		if (rowlist[i] < 0 || rowlist[i] >= nrows)
		{
			fprintf (stderr, "entry %d in rowlist out of range\n", i);
			rval = 1;
			goto CLEANUP;
		}
	}

	rval = ILLlib_getrows (p->lp, num, rowlist, rowcnt, rowbeg, rowind,
												 rowval, rhs, sense, names);
	ILL_CLEANUP_IF (rval);


CLEANUP:

	ILL_RETURN (rval, "QSget_rows_list");
}

QSLIB_INTERFACE int QSget_rows (QSdata * p,
																int **rowcnt,
																int **rowbeg,
																int **rowind,
																EGlpNum_t ** rowval,
																EGlpNum_t ** rhs,
																char **sense,
																char ***names)
{
	int rval = 0;
	int i,
	  nrows;
	int *rowlist = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	nrows = QSget_rowcount (p);
	if (nrows > 0)
	{
		ILL_SAFE_MALLOC (rowlist, nrows, int);
		for (i = 0; i < nrows; i++)
		{
			rowlist[i] = i;
		}
		rval = ILLlib_getrows (p->lp, nrows, rowlist, rowcnt, rowbeg, rowind,
													 rowval, rhs, sense, names);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	ILL_IFFREE (rowlist, int);
	ILL_RETURN (rval, "QSget_rows");
}

QSLIB_INTERFACE int QSget_columns_list (QSdata * p,
																				int num,
																				int *collist,
																				int **colcnt,
																				int **colbeg,
																				int **colind,
																				EGlpNum_t ** colval,
																				EGlpNum_t ** obj,
																				EGlpNum_t ** lower,
																				EGlpNum_t ** upper,
																				char ***names)
{
	int rval = 0;
	int j,
	  ncols;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	ncols = QSget_colcount (p);
	for (j = 0; j < num; j++)
	{
		if (collist[j] < 0 || collist[j] >= ncols)
		{
			fprintf (stderr, "entry %d in collist out of range\n", j);
			rval = 1;
			goto CLEANUP;
		}
	}

	rval = ILLlib_getcols (p->lp, num, collist, colcnt, colbeg, colind,
												 colval, obj, lower, upper, names);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_columns_list");
}

QSLIB_INTERFACE int QSget_columns (QSdata * p,
																	 int **colcnt,
																	 int **colbeg,
																	 int **colind,
																	 EGlpNum_t ** colval,
																	 EGlpNum_t ** obj,
																	 EGlpNum_t ** lower,
																	 EGlpNum_t ** upper,
																	 char ***names)
{
	int rval = 0;
	int j,
	  ncols;
	int *collist = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	ncols = QSget_colcount (p);
	if (ncols > 0)
	{
		ILL_SAFE_MALLOC (collist, ncols, int);
		for (j = 0; j < ncols; j++)
		{
			collist[j] = j;
		}
		rval = ILLlib_getcols (p->lp, ncols, collist, colcnt, colbeg, colind,
													 colval, obj, lower, upper, names);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	ILL_IFFREE (collist, int);
	ILL_RETURN (rval, "QSget_columns");
}

QSLIB_INTERFACE char *QSget_probname (QSdata * p)
{
	int rval = 0;
	char *name = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	ILL_UTIL_STR (name, p->name);

CLEANUP:
	ILL_RETURN_PTR (name, "QSget_probname");
}

QSLIB_INTERFACE char *QSget_objname (QSdata * p)
{
	int rval = 0;
	char *name = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (p->qslp->objname != 0)
	{
		ILL_UTIL_STR (name, p->qslp->objname);
	}

CLEANUP:
	ILL_RETURN_PTR (name, "QSget_objname");
}

QSLIB_INTERFACE int QSget_rownames (QSdata * p,
																		char **rownames)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_rownames (p->lp, rownames);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_rownames");
}

QSLIB_INTERFACE int QSget_colnames (QSdata * p,
																		char **colnames)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_colnames (p->lp, colnames);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_colnames");
}

QSLIB_INTERFACE int QSget_bound (QSdata * p,
																 int colindex,
																 int lu,
																 EGlpNum_t * bound)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_getbnd (p->lp, colindex, lu, bound);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_bound");
}

QSLIB_INTERFACE int QSget_bounds (QSdata * p,
																	EGlpNum_t * lower,
																	EGlpNum_t * upper)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_getbnds (p->lp, lower, upper);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_bounds");
}

QSLIB_INTERFACE int QSget_intflags (QSdata * p,
																		int *intflags)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (intflags == 0)
	{
		rval = 1;
		ILL_CLEANUP;
	}

	rval = ILLlib_getintflags (p->lp, intflags);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_intflags");
}

QSLIB_INTERFACE int QSget_intcount (QSdata * p,
																		int *count)
{
	int j,
	  ncols,
	  cnt = 0,
	  rval = 0;
	int *intflags = 0;

	*count = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	ncols = QSget_colcount (p);

	if (ncols > 0)
	{
		ILL_SAFE_MALLOC (intflags, ncols, int);

		rval = ILLlib_getintflags (p->lp, intflags);
		ILL_CLEANUP_IF (rval);

		for (j = 0; j < ncols; j++)
		{
			if (intflags[j] > 0)
				cnt++;
		}
	}

CLEANUP:

	*count = cnt;
	ILL_IFFREE (intflags, int);
	ILL_RETURN (rval, "QSget_intcount");
}

QSLIB_INTERFACE int QSget_column_index (QSdata * p,
																				const char *name,
																				int *colindex)
{
	int rval = 0;

	*colindex = -1;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_colindex (p->lp, name, colindex);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_column_index");
}

QSLIB_INTERFACE int QSget_row_index (QSdata * p,
																		 const char *name,
																		 int *rowindex)
{
	int rval = 0;

	*rowindex = -1;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_rowindex (p->lp, name, rowindex);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "QSget_row_index");
}

QSLIB_INTERFACE int QSwrite_prob (QSdata * p,
																	const char *filename,
																	const char *filetype)
{
	FILE *file = NULL;
	int rval = 0;

	if (filename == NULL)
	{
		file = stdout;
	}
	else
	{
		if ((file = fopen (filename, "w")) == 0)
		{
			perror (filename);
			fprintf (stderr, "Unable to open \"%s\" for output.\n", filename);
		}
	}
	ILL_CHECKnull (file, NULL);
	rval = QSwrite_prob_file (p, file, filetype);
	if (file && (file != stdout) && (file != stderr))
	{
		fclose (file);
	}
CLEANUP:
	ILL_RETURN (rval, "QSwrite_prob_file");
}

QSLIB_INTERFACE int QSwrite_prob_file (QSdata * p,
																			 FILE * out,
																			 const char *filetype)
{
	int rval = 0;
	qsstring_reporter rep;

	ILLstring_reporter_copy (&rep, &p->qslp->reporter);
	ILLstring_reporter_init (&p->qslp->reporter,
													 (qsreport_string_fct) fprintf, out);
	rval = QSreport_prob (p, filetype, NULL);
	ILLstring_reporter_copy (&p->qslp->reporter, &rep);
	ILL_RESULT (rval, "QSwrite_prob_file");
}

QSLIB_INTERFACE void QSfree (void *ptr)
{
	ILL_IFFREE (ptr, void);
}

QSLIB_INTERFACE int QSset_param (QSdata * p,
																 int whichparam,
																 int newvalue)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	switch (whichparam)
	{
	case QS_PARAM_PRIMAL_PRICING:
		if (newvalue == QS_PRICE_PDANTZIG ||
				newvalue == QS_PRICE_PDEVEX ||
				newvalue == QS_PRICE_PSTEEP || newvalue == QS_PRICE_PMULTPARTIAL)
		{
			p->pricing->pI_price = newvalue;
			p->pricing->pII_price = newvalue;
		}
		else
		{
			fprintf (stderr, "illegal value for QS_PARAM_PRIMAL_PRICING\n");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_DUAL_PRICING:
		if (newvalue == QS_PRICE_DDANTZIG ||
				newvalue == QS_PRICE_DSTEEP ||
				newvalue == QS_PRICE_DMULTPARTIAL || newvalue == QS_PRICE_DDEVEX)
		{
			p->pricing->dI_price = newvalue;
			p->pricing->dII_price = newvalue;
		}
		else
		{
			fprintf (stderr, "illegal value for QS_PARAM_DUAL_PRICING\n");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_SIMPLEX_DISPLAY:
		if (newvalue == 0 || (newvalue > 0 && newvalue < 4))
		{
			p->simplex_display = newvalue;
		}
		else
		{
			fprintf (stderr, "illegal value for QS_PARAM_SIMPLEX_DISPLAY\n");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_SIMPLEX_MAX_ITERATIONS:
		if (newvalue > 0)
		{
			p->lp->maxiter = newvalue;
		}
		else
		{
			fprintf (stderr, "illegal value for QS_PARAM_SIMPLEX_MAX_ITERATIONS\n");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_SIMPLEX_SCALING:
		if (newvalue == 0 || newvalue == 1)
		{
			p->simplex_scaling = newvalue;
		}
		else
		{
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

	ILL_RETURN (rval, "QSset_param");
}

QSLIB_INTERFACE int QSset_param_EGlpNum (QSdata * p,
																				 int whichparam,
																				 EGlpNum_t newvalue)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	switch (whichparam)
	{
	case QS_PARAM_SIMPLEX_MAX_TIME:
		if (EGlpNumIsLess (zeroLpNum, newvalue))
		{
			p->lp->maxtime = EGlpNumToLf (newvalue);
		}
		else
		{
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

QSLIB_INTERFACE int QSget_param (QSdata * p,
																 int whichparam,
																 int *value)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (!value)
	{
		fprintf (stderr, "QSget_param call without a value pointer\n");
		rval = 1;
		goto CLEANUP;
	}

	switch (whichparam)
	{
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

	ILL_RETURN (rval, "QSget_param");
}

QSLIB_INTERFACE int QSget_param_EGlpNum (QSdata * p,
																				 int whichparam,
																				 EGlpNum_t * value)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	ILL_CLEANUP_IF (rval);

	if (!value)
	{
		fprintf (stderr, "QSget_param_double call without a value pointer\n");
		rval = 1;
		goto CLEANUP;
	}

	switch (whichparam)
	{
	case QS_PARAM_SIMPLEX_MAX_TIME:
		EGlpNumSet (*value, p->lp->maxtime);
		break;
	default:
		fprintf (stderr, "unknown parameter: %d\n", whichparam);
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	ILL_RETURN (rval, "QSget_param_double");
}

static int check_qsdata_pointer (QSdata * p)
{
	if (p == NULL)
	{
		fprintf (stderr, "NULL QSprob pointer\n");
		return 1;
	}
	else
	{
		return 0;
	}
}

static int formatIsMps (const char *filetype,
												int *isMps)
{
	int rval = 0;

	if (!strcasecmp (filetype, "MPS"))
	{
		*isMps = 1;
	}
	else if (!strcasecmp (filetype, "LP"))
	{
		*isMps = 0;
	}
	else
	{
		fprintf (stderr, "Unknown prob-file type: %s\n", filetype);
		rval = 1;
		ILL_CLEANUP;
	}
CLEANUP:
	return rval;
}


/****************************************************************************/
/* 
 * undocumentyed functions 
 */

static void check_pointer (void *p,
													 const char *fct,
													 const char *param)
{
	if (p == NULL)
		fprintf (stderr, "NULL %s argument to %s\n", param, fct);
}

/* QSline_reader: 
 *    used by mps/lp reader to get input lines 
 *    by default input is read froma FILE* via fgets 
 */
QSLIB_INTERFACE QSline_reader QSline_reader_new (void *fct,
																								 void *data_src)
{
	check_pointer (fct, "QSline_reader_new", "fct");
	check_pointer (data_src, "QSline_reader_new", "data_src");
	return ILLline_reader_new ((qsread_line_fct) fct, data_src);
}

QSLIB_INTERFACE void QSline_reader_set_error_collector (QSline_reader reader,
																												QSerror_collector
																												collector)
{
	check_pointer (reader, "QSline_reader_set_error_collector", "reader");
	check_pointer (collector, "QSline_reader_set_error_collector", "collector");
	reader->error_collector = collector;
}

QSLIB_INTERFACE void QSline_reader_free (QSline_reader reader)
{
	ILLline_reader_free (reader);
}

QSLIB_INTERFACE char *QSline_reader_get (QSline_reader reader,
																				 char *s,
																				 int size)
{
	check_pointer (reader, "QSline_reader_get", "reader");
	check_pointer (s, "QSline_reader_get", "s");
	return ILLline_reader_get (s, size, reader);
}


QSLIB_INTERFACE QSerror_collector QSerror_collector_new (void *fct,
																												 void *dest)
{
	check_pointer (fct, "QSerror_collector_new", "fct");
	return ILLerror_collector_new ((qsadd_error_fct) fct, dest);
}

QSLIB_INTERFACE
	QSerror_collector QSerror_memory_collector_new (QSerror_memory mem)
{
	check_pointer (mem, "QSerror_memory_collector_new", "mem");
	return ILLerror_memory_collector_new (mem);
}

QSLIB_INTERFACE void QSerror_collector_free (QSerror_collector c)
{
	ILLerror_collector_free (c);
}

QSLIB_INTERFACE QSdata *QSget_prob (QSline_reader reader,
																		const char *probname,
																		const char *filetype)
{
	int isMps,
	  rval = 0;
	QSdata *p = 0;

	if ((filetype != NULL) && !strcasecmp (filetype, "MPS"))
	{
		isMps = 1;
	}
	else if ((filetype != NULL) && !strcasecmp (filetype, "LP"))
	{
		isMps = 0;
	}
	else
	{
		fprintf (stderr, "Unknown prob-file type: %s\n",
						 (filetype != NULL) ? filetype : "NULL");
		rval = 1;
		ILL_CLEANUP;
	}

	p = ILLread (reader, probname, isMps);
	ILL_CHECKnull (p, NULL);

	ILL_FAILfalse (p->qslp != NULL, "If there's a p there must be a p-qslp");
	ILL_IFFREE (p->name, char);
	ILL_UTIL_STR (p->name, p->qslp->probname);
	ILLsimplex_load_lpinfo (p->qslp, p->lp);

CLEANUP:

	if (rval != 0)
	{
		QSfree_prob (p);
		p = 0;
	}
	return p;
}

QSLIB_INTERFACE int QSreport_prob (QSdata * p,
																	 const char *filetype,
																	 qserror_collector * c)
{
	int isMps,
	  rval = 0;

	rval = formatIsMps (filetype, &isMps);
	ILL_CLEANUP_IF (rval);
	if (isMps)
	{
		rval = ILLwrite_mps (p->qslp, c);
	}
	else
	{
		rval = ILLwrite_lp (p->qslp, c);
	}
CLEANUP:
	ILL_RESULT (rval, "QSreport_prob");
}

QSLIB_INTERFACE char *QSversion (void)
{
	char *name = 0;
	name = EGsMalloc(char,256);
	snprintf(name,(size_t)255,"QSopt_ex 2.5.0 (build %s-%s)",__DATE__,__TIME__);
	return name;
}

/* QSstring_reporter: 
 *    used by solver code to report feedback 
 *    by default feedback is sent to stdout via fprintf
 */
QSLIB_INTERFACE void QSset_reporter (QSprob prob,
																		 int skip,
																		 void *fct,
																		 void *dest)
{
	int rval = 0;
	rval = check_qsdata_pointer (prob);
	if (rval != 0)
		return;

	check_pointer (fct, "QSset_reporter", "fct");

	ILL_FAILtrue (prob->lp == NULL, "QSprob internal error: prob->lp == NULL");
	ILLstring_reporter_init (&prob->qslp->reporter,
													 (qsreport_string_fct) fct, dest);

	prob->lp->iterskip = skip;
CLEANUP:
	return;
}

QSLIB_INTERFACE const char *QSformat_error_type_string (int tp)
{
	const char *type = "Error";

	if (tp == QS_DATA_ERROR)
	{
		type = "Data Error";
	}
	if (tp == QS_DATA_WARN)
	{
		type = "Data Warning";
	}
	if (tp == QS_MPS_FORMAT_ERROR)
	{
		type = "MPS Error";
	}
	if (tp == QS_MPS_FORMAT_WARN)
	{
		type = "MPS Warning";
	}
	if (tp == QS_LP_FORMAT_ERROR)
	{
		type = "LP Error";
	}
	if (tp == QS_LP_FORMAT_WARN)
	{
		type = "LP Warning";
	}
	return type;
}

QSLIB_INTERFACE int QSerror_get_type (QSformat_error error)
{
	check_pointer (error, "QSerror_get_type", "error");
	return error->type;
}

QSLIB_INTERFACE const char *QSerror_get_desc (QSformat_error error)
{
	check_pointer (error, "QSerror_get_desc", "error");
	return error->desc;
}

QSLIB_INTERFACE int QSerror_get_line_number (QSformat_error error)
{
	check_pointer (error, "QSerror_get_line_number", "error");
	return error->lineNumber;
}

QSLIB_INTERFACE int QSerror_get_pos (QSformat_error error)
{
	check_pointer (error, "QSerror_get_pos", "error");
	return error->at;
}

QSLIB_INTERFACE const char *QSerror_get_line (QSformat_error error)
{
	check_pointer (error, "QSerror_get_line", "error");
	return error->theLine;
}

QSLIB_INTERFACE void QSerror_print (FILE * f,
																		QSformat_error error)
{
	check_pointer (f, "QSerror_print", "f");
	if (error == NULL)
	{
		fprintf (stderr, "0\n");
	}
	else
	{
		ILLformat_error_print (f, error);
	}
}

QSLIB_INTERFACE QSerror_memory QSerror_memory_create (int takeErrorLines)
{
	return ILLerror_memory_create (takeErrorLines);
}

QSLIB_INTERFACE void QSerror_memory_free (QSerror_memory mem)
{
	if (mem != NULL)
	{
		ILLerror_memory_free (mem);
	}
}

QSLIB_INTERFACE int QSerror_memory_get_nerrors (QSerror_memory mem)
{
	check_pointer (mem, "QSerror_memory_get_nerrors", "mem");
	return mem->nerror;
}

QSLIB_INTERFACE int QSerror_memory_get_nof (QSerror_memory mem,
																						int type)
{
	check_pointer (mem, "QSerror_memory_get_nerrors", "mem");
	if (0 <= type && type < QS_INPUT_NERROR)
	{
		return mem->has_error[type];
	}
	else
	{
		ILL_REPRT ("bad error type");
		return 0;
	}
}

QSLIB_INTERFACE QSformat_error
QSerror_memory_get_last_error (QSerror_memory mem)
{
	check_pointer (mem, "QSerror_memory_get_last_errors", "mem");
	return mem->error_list;
}

QSLIB_INTERFACE QSformat_error QSerror_memory_get_prev_error (QSformat_error e)
{
	check_pointer (e, "QSerror_memory_get_prev_errors", "e");
	if (e != NULL)
		e = e->next;
	return e;
}
