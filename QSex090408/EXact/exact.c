/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
   and algorithms commons in solving MIP's

Copyright (C) 2005 Daniel Espinoza.

This library is free software; you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

You should have received a copy of the GNU Lesser General Public License
   along with this library; if not, write to the Free Software Foundation,
   Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */
/* ========================================================================= */

#define EGLPNUM_TYPE 9
#include "exact.h"
#include "util.h"
extern mpq_t mpq_ILL_MAXDOUBLE;
extern mpq_t mpq_ILL_MINDOUBLE;
extern void mpq_ILLfct_set_variable_type (mpq_lpinfo * lp);
extern void mpq_ILLfct_compute_piz (mpq_lpinfo * lp);
extern void mpq_ILLfct_compute_dz (mpq_lpinfo * lp);
extern void mpq_ILLfct_compute_xbz (mpq_lpinfo * lp);
extern void mpq_ILLfct_check_pfeasible (mpq_lpinfo * lp,
      mpq_feas_info * fs,
      mpq_t ftol);
extern void mpq_ILLfct_check_dfeasible (mpq_lpinfo * lp,
      mpq_feas_info * fs,
      mpq_t ftol);
extern void mpq_ILLfct_set_status_values (mpq_lpinfo * lp,
      int pstatus,
      int dstatus,
      int ptype,
      int dtype);
extern int mpq_grab_cache (mpq_QSdata * p,
      int status);
extern void mpq_ILLfct_compute_phaseI_piz (mpq_lpinfo * lp);

/* ========================================================================= */
/* this copies (in a loose way) a mpq pricing structure to a double *        */
/* ========================================================================= */
#define QScopy_pricing_mpq_dbl(p2,p) ({\
    p2->pricing->p_strategy = p->pricing->p_strategy;\
    p2->pricing->d_strategy = p->pricing->d_strategy;\
    p2->pricing->pI_price = p->pricing->pI_price;\
    p2->pricing->pII_price = p->pricing->pII_price;\
    p2->pricing->dI_price = p->pricing->dI_price;\
    p2->pricing->dII_price = p->pricing->dII_price;\
    p2->pricing->cur_price = p->pricing->cur_price;\
    p2->pricing->p_scaleinf = QScopy_array_mpq_dbl(p->pricing->p_scaleinf);\
    p2->pricing->d_scaleinf = QScopy_array_mpq_dbl(p->pricing->d_scaleinf);\
    p2->pricing->pdinfo.ninit = p->pricing->pdinfo.ninit;\
    p2->pricing->pdinfo.norms = QScopy_array_mpq_dbl(p->pricing->pdinfo.norms);\
    p2->pricing->pdinfo.refframe = p->pricing->pdinfo.refframe;\
    p2->pricing->psinfo.norms = QScopy_array_mpq_dbl(p->pricing->psinfo.norms);\
    p2->pricing->pmpinfo.k = p->pricing->pmpinfo.k;\
    p2->pricing->pmpinfo.cgroup = p->pricing->pmpinfo.cgroup;\
    p2->pricing->pmpinfo.ngroups = p->pricing->pmpinfo.ngroups;\
    p2->pricing->pmpinfo.gstart = p->pricing->pmpinfo.gstart;\
    p2->pricing->pmpinfo.gshift = p->pricing->pmpinfo.gshift;\
    p2->pricing->pmpinfo.gsize = p->pricing->pmpinfo.gsize;\
    p2->pricing->pmpinfo.bsize = p->pricing->pmpinfo.bsize;\
    p2->pricing->pmpinfo.bucket = p->pricing->pmpinfo.bucket;\
    p2->pricing->pmpinfo.perm = p->pricing->pmpinfo.perm;\
    p2->pricing->pmpinfo.infeas = QScopy_array_mpq_dbl(p->pricing->pmpinfo.infeas);\
    p2->pricing->ddinfo.ninit = p->pricing->ddinfo.ninit;\
    p2->pricing->ddinfo.norms = QScopy_array_mpq_dbl(p->pricing->ddinfo.norms);\
    p2->pricing->ddinfo.refframe = p->pricing->ddinfo.refframe;\
    p2->pricing->dsinfo.norms = QScopy_array_mpq_dbl(p->pricing->dsinfo.norms);\
    p2->pricing->dmpinfo.k = p->pricing->dmpinfo.k;\
    p2->pricing->dmpinfo.cgroup = p->pricing->dmpinfo.cgroup;\
    p2->pricing->dmpinfo.ngroups = p->pricing->dmpinfo.ngroups;\
    p2->pricing->dmpinfo.gstart = p->pricing->dmpinfo.gstart;\
    p2->pricing->dmpinfo.gshift = p->pricing->dmpinfo.gshift;\
    p2->pricing->dmpinfo.gsize = p->pricing->dmpinfo.gsize;\
    p2->pricing->dmpinfo.bsize = p->pricing->dmpinfo.bsize;\
    p2->pricing->dmpinfo.bucket = p->pricing->dmpinfo.bucket;\
    p2->pricing->dmpinfo.perm = p->pricing->dmpinfo.perm;\
    p2->pricing->dmpinfo.infeas = QScopy_array_mpq_dbl(p->pricing->dmpinfo.infeas);\
    dbl_ILLheap_init(&(p2->pricing->h));\
    p2->pricing->htrigger = mpq_get_d(p->pricing->htrigger);\
    p2->pricing->hineff = p->pricing->hineff;\
    p2->pricing->init = p->pricing->init;})

/* ========================================================================= */
/* copy (in a loose way) a mpq pricing structure to a mpf                    */
/* ========================================================================= */
#define QScopy_pricing_mpq_mpf(p2,p) ({\
    p2->pricing->p_strategy = p->pricing->p_strategy;\
    p2->pricing->d_strategy = p->pricing->d_strategy;\
    p2->pricing->pI_price = p->pricing->pI_price;\
    p2->pricing->pII_price = p->pricing->pII_price;\
    p2->pricing->dI_price = p->pricing->dI_price;\
    p2->pricing->dII_price = p->pricing->dII_price;\
    p2->pricing->cur_price = p->pricing->cur_price;\
    p2->pricing->p_scaleinf = QScopy_array_mpq_mpf(p->pricing->p_scaleinf);\
    p2->pricing->d_scaleinf = QScopy_array_mpq_mpf(p->pricing->d_scaleinf);\
    p2->pricing->pdinfo.ninit = p->pricing->pdinfo.ninit;\
    p2->pricing->pdinfo.norms = QScopy_array_mpq_mpf(p->pricing->pdinfo.norms);\
    p2->pricing->pdinfo.refframe = p->pricing->pdinfo.refframe;\
    p2->pricing->psinfo.norms = QScopy_array_mpq_mpf(p->pricing->psinfo.norms);\
    p2->pricing->pmpinfo.k = p->pricing->pmpinfo.k;\
    p2->pricing->pmpinfo.cgroup = p->pricing->pmpinfo.cgroup;\
    p2->pricing->pmpinfo.ngroups = p->pricing->pmpinfo.ngroups;\
    p2->pricing->pmpinfo.gstart = p->pricing->pmpinfo.gstart;\
    p2->pricing->pmpinfo.gshift = p->pricing->pmpinfo.gshift;\
    p2->pricing->pmpinfo.gsize = p->pricing->pmpinfo.gsize;\
    p2->pricing->pmpinfo.bsize = p->pricing->pmpinfo.bsize;\
    p2->pricing->pmpinfo.bucket = p->pricing->pmpinfo.bucket;\
    p2->pricing->pmpinfo.perm = p->pricing->pmpinfo.perm;\
    p2->pricing->pmpinfo.infeas = QScopy_array_mpq_mpf(p->pricing->pmpinfo.infeas);\
    p2->pricing->ddinfo.ninit = p->pricing->ddinfo.ninit;\
    p2->pricing->ddinfo.norms = QScopy_array_mpq_mpf(p->pricing->ddinfo.norms);\
    p2->pricing->ddinfo.refframe = p->pricing->ddinfo.refframe;\
    p2->pricing->dsinfo.norms = QScopy_array_mpq_mpf(p->pricing->dsinfo.norms);\
    p2->pricing->dmpinfo.k = p->pricing->dmpinfo.k;\
    p2->pricing->dmpinfo.cgroup = p->pricing->dmpinfo.cgroup;\
    p2->pricing->dmpinfo.ngroups = p->pricing->dmpinfo.ngroups;\
    p2->pricing->dmpinfo.gstart = p->pricing->dmpinfo.gstart;\
    p2->pricing->dmpinfo.gshift = p->pricing->dmpinfo.gshift;\
    p2->pricing->dmpinfo.gsize = p->pricing->dmpinfo.gsize;\
    p2->pricing->dmpinfo.bsize = p->pricing->dmpinfo.bsize;\
    p2->pricing->dmpinfo.bucket = p->pricing->dmpinfo.bucket;\
    p2->pricing->dmpinfo.perm = p->pricing->dmpinfo.perm;\
    p2->pricing->dmpinfo.infeas = QScopy_array_mpq_mpf(p->pricing->dmpinfo.infeas);\
    mpf_ILLheap_init(&(p2->pricing->h));\
    mpq_EGlpNumSet_mpf(p->pricing->htrigger,p2->pricing->htrigger );\
    p2->pricing->hineff = p->pricing->hineff;\
    p2->pricing->init = p->pricing->init;})

/* ========================================================================= */
int QSexact_print_sol (mpq_QSdata * p, FILE * out_f)
{
    int rval = 0, status;
    const int ncols = mpq_QSget_colcount (p);
    const int nrows = mpq_QSget_rowcount (p);
    mpq_t *x = mpq_EGlpNumAllocArray (ncols);
    mpq_t *rc = mpq_EGlpNumAllocArray (ncols);
    mpq_t *slack = mpq_EGlpNumAllocArray (nrows);
    mpq_t *pi = mpq_EGlpNumAllocArray (nrows);
    mpq_t value;
    register int i;
    char *str1 = 0;
    mpq_init (value);
    rval = mpq_QSget_status (p, &status);
    TESTG (rval, CLEANUP, "mpq_QSget_status failed");
    rval = mpq_QSget_x_array (p, x);
    if (rval)
	mpq_EGlpNumFreeArray (x);
    rval = mpq_QSget_slack_array (p, slack);
    if (rval)
	mpq_EGlpNumFreeArray (slack);
    rval = mpq_QSget_pi_array (p, pi);
    if (rval)
	mpq_EGlpNumFreeArray (pi);
    rval = mpq_QSget_rc_array (p, rc);
    if (rval)
	mpq_EGlpNumFreeArray (rc);
    rval = 0;

    switch (status) {
    case QS_LP_OPTIMAL:
	rval = mpq_QSget_objval (p, &value);
	TESTG (rval, CLEANUP, "mpq_QSget_objvail failed");
	str1 = mpq_EGlpNumGetStr (value);
	fprintf (out_f, "status OPTIMAL\n\tValue = %s\n", str1);
	free (str1);
	str1 = 0;
	break;
    case QS_LP_INFEASIBLE:
	fprintf (out_f, "status INFEASIBLE\n");
	break;
    case QS_LP_UNBOUNDED:
	fprintf (out_f, "status UNBOUNDED\n");
	break;
    case QS_LP_ITER_LIMIT:
    case QS_LP_TIME_LIMIT:
    case QS_LP_UNSOLVED:
    case QS_LP_ABORTED:
    case QS_LP_MODIFIED:
	fprintf (out_f, "status NOT_SOLVED\n");
	break;
    }
    if (x) {
	fprintf (out_f, "VARS:\n");
	for (i = 0; i < ncols; i++)
	    if (!mpq_equal (x[i], __zeroLpNum_mpq__)) {
		str1 = mpq_EGlpNumGetStr (x[i]);
		fprintf (out_f, "%s = %s\n", p->qslp->colnames[i], str1);
		free (str1);
	    }
    }
    if (rc) {
	fprintf (out_f, "REDUCED COST:\n");
	for (i = 0; i < ncols; i++)
	    if (!mpq_equal (rc[i], __zeroLpNum_mpq__)) {
		str1 = mpq_EGlpNumGetStr (rc[i]);
		fprintf (out_f, "%s = %s\n", p->qslp->colnames[i], str1);
		free (str1);
	    }
    }
    if (pi) {
	fprintf (out_f, "PI:\n");
	for (i = 0; i < nrows; i++)
	    if (!mpq_equal (pi[i], __zeroLpNum_mpq__)) {
		str1 = mpq_EGlpNumGetStr (pi[i]);
		fprintf (out_f, "%s = %s\n", p->qslp->rownames[i], str1);
		free (str1);
	    }
    }
    if (slack) {
	fprintf (out_f, "SLACK:\n");
	for (i = 0; i < nrows; i++)
	    if (!mpq_equal (slack[i], __zeroLpNum_mpq__)) {
		str1 = mpq_EGlpNumGetStr (slack[i]);
		fprintf (out_f, "%s = %s\n", p->qslp->rownames[i], str1);
		free (str1);
	    }
    }
    /* ending */
CLEANUP:
    if (x)
	mpq_EGlpNumFreeArray (x);
    if (pi)
	mpq_EGlpNumFreeArray (pi);
    if (rc)
	mpq_EGlpNumFreeArray (rc);
    if (slack)
	mpq_EGlpNumFreeArray (slack);
    mpq_clear (value);
    return rval;
}

/* ========================================================================= */
dbl_QSdata *QScopy_prob_mpq_dbl (mpq_QSdata * p, const char *newname)
{
    int rval = 0;
    int j, col, beg, pindex, hit;
    dbl_QSdata *p2 = 0;
    char *coln;
    char buf[ILL_namebufsize];
    double *rhs = 0, *rangeval = 0;

    /* test that the input is not NULL */
    ADVCHECKRVAL (!p, 0);

    p2 = dbl_QScreate_prob (newname, p->qslp->objsense);
    if (!p2)
	goto CLEANUP;

    rhs = QScopy_array_mpq_dbl (p->qslp->rhs);
    rangeval = QScopy_array_mpq_dbl (p->qslp->rangeval);
    rval = dbl_ILLlib_newrows (p2->lp, 0, p->qslp->nrows, rhs, p->qslp->sense,
	rangeval, (const char **) p->qslp->rownames);
    dbl_EGlpNumFreeArray (rhs);
    dbl_EGlpNumFreeArray (rangeval);
    CHECKRVALG (rval, CLEANUP);

    rhs = QScopy_array_mpq_dbl (p->qslp->A.matval);
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
	/* code in lib.c.                                  */
	rval = dbl_ILLlib_addcol (p2->lp, 0, p->qslp->A.matcnt[col],
	    p->qslp->A.matind + beg, rhs + beg,
	    mpq_get_d (p->qslp->obj[col]),
	    mpq_equal (p->qslp->lower[col],
		mpq_ILL_MINDOUBLE) ? dbl_ILL_MINDOUBLE
	    : mpq_get_d (p->qslp->lower[col]),
	    mpq_equal (p->qslp->upper[col],
		mpq_ILL_MAXDOUBLE) ? dbl_ILL_MAXDOUBLE
	    : mpq_get_d (p->qslp->upper[col]), coln, 0);
	ILL_CLEANUP_IF (rval);
    }
    dbl_EGlpNumFreeArray (rhs);

    p2->qslp->objsense = p->qslp->objsense;

    p2->factorok = 0;
    p2->simplex_display = p->simplex_display;
    p2->simplex_scaling = p->simplex_scaling;
    QScopy_pricing_mpq_dbl (p2, p);

    if (p->qslp->objname != 0) {
	ILL_UTIL_STR (p2->qslp->objname, p->qslp->objname);
    } else {
	strcpy (buf, "obj");
	rval = ILLsymboltab_uname (&p2->qslp->rowtab, buf, "", NULL);
	CHECKRVALG (rval, CLEANUP);
	ILL_UTIL_STR (p2->qslp->objname, buf);
    }
    rval = ILLsymboltab_register (&p2->qslp->rowtab, p2->qslp->objname,
	-1, &pindex, &hit);
    rval = rval || hit;
    CHECKRVALG (rval, CLEANUP);

    ILLstring_reporter_copy (&p2->qslp->reporter, &p->qslp->reporter);

    if (p->qslp->intmarker != 0) {
	ILL_SAFE_MALLOC (p2->qslp->intmarker, p->qslp->nstruct, char);
	for (j = 0; j < p->qslp->nstruct; j++) {
	    p2->qslp->intmarker[j] = p->qslp->intmarker[j];
	}
    }
CLEANUP:

    if (rval) {
	dbl_QSfree_prob (p2);
	p2 = 0;
    }
#if QSEXACT_SAVE_INT
    else {
	dbl_QSwrite_prob (p2, "prob.dbl.lp", "LP");
    }
#endif

    return p2;
}

/* ========================================================================= */
mpf_QSdata *QScopy_prob_mpq_mpf (mpq_QSdata * p, const char *newname)
{
    int rval = 0;
    int j, col, beg, pindex, hit;
    mpf_QSdata *p2 = 0;
    char *coln;
    char buf[ILL_namebufsize];
    mpf_t *rhs = 0, *rangeval = 0, obj, lower, upper;
    mpf_init (obj);
    mpf_init (lower);
    mpf_init (upper);

    /* test that the input is not NULL */
    ADVCHECKRVAL (!p, 0);

    p2 = mpf_QScreate_prob (newname, p->qslp->objsense);
    if (!p2)
	goto CLEANUP;

    rhs = QScopy_array_mpq_mpf (p->qslp->rhs);
    rangeval = QScopy_array_mpq_mpf (p->qslp->rangeval);
    rval = mpf_ILLlib_newrows (p2->lp, 0, p->qslp->nrows,
	rhs, p->qslp->sense, rangeval,
	(const char **) p->qslp->rownames);
    mpf_EGlpNumFreeArray (rhs);
    mpf_EGlpNumFreeArray (rangeval);
    CHECKRVALG (rval, CLEANUP);

    rhs = QScopy_array_mpq_mpf (p->qslp->A.matval);
    for (j = 0; j < p->qslp->nstruct; j++) {
	col = p->qslp->structmap[j];
	if (p->qslp->colnames)
	    coln = p->qslp->colnames[j];
	else
	    coln = 0;
	beg = p->qslp->A.matbeg[col];
	mpf_set_q (obj, p->qslp->obj[col]);
	if (mpq_equal (p->qslp->lower[col], mpq_ILL_MINDOUBLE))
	    mpf_set (lower, mpf_ILL_MINDOUBLE);
	else
	    mpf_set_q (lower, p->qslp->lower[col]);
	if (mpq_equal (p->qslp->upper[col], mpq_ILL_MAXDOUBLE))
	    mpf_set (upper, mpf_ILL_MAXDOUBLE);
	else
	    mpf_set_q (upper, p->qslp->upper[col]);

	rval = mpf_ILLlib_addcol (p2->lp, 0,
	    p->qslp->A.matcnt[col],
	    p->qslp->A.matind + beg, rhs + beg,
	    obj, lower, upper, coln, 0);
	ILL_CLEANUP_IF (rval);
    }
    mpf_EGlpNumFreeArray (rhs);

    p2->qslp->objsense = p->qslp->objsense;

    p2->factorok = 0;
    p2->simplex_display = p->simplex_display;
    p2->simplex_scaling = p->simplex_scaling;
    QScopy_pricing_mpq_mpf (p2, p);

    if (p->qslp->objname != 0) {
	ILL_UTIL_STR (p2->qslp->objname, p->qslp->objname);
    } else {
	strcpy (buf, "obj");
	rval = ILLsymboltab_uname (&p2->qslp->rowtab, buf, "", NULL);
	CHECKRVALG (rval, CLEANUP);
	ILL_UTIL_STR (p2->qslp->objname, buf);
    }
    rval = ILLsymboltab_register (&p2->qslp->rowtab, p2->qslp->objname,
	-1, &pindex, &hit);
    rval = rval || hit;
    CHECKRVALG (rval, CLEANUP);

    ILLstring_reporter_copy (&p2->qslp->reporter, &p->qslp->reporter);

    if (p->qslp->intmarker != 0) {
	ILL_SAFE_MALLOC (p2->qslp->intmarker, p->qslp->nstruct, char);
	for (j = 0; j < p->qslp->nstruct; j++) {
	    p2->qslp->intmarker[j] = p->qslp->intmarker[j];
	}
    }
CLEANUP:

    if (rval) {
	mpf_QSfree_prob (p2);
	p2 = 0;
    }
#if QSEXACT_SAVE_INT
    else {
	mpf_QSwrite_prob (p2, "prob.mpf.lp", "LP");
    }
#endif
    mpf_clear (obj);
    mpf_clear (lower);
    mpf_clear (upper);

    return p2;
}

#if QSEXACT_SAVE_OPTIMAL
/* ========================================================================= */
/* used to enumerate the generated optimal tests */
static int QSEXACT_SAVE_OPTIMAL_IND = 0;
#endif

/* ========================================================================= */
int QSexact_optimal_test (mpq_QSdata * p, mpq_t * p_sol, mpq_t * d_sol,
      QSbasis * basis)
{
    /* local variables */
    register int i, j;
    mpq_ILLlpdata *qslp = p->lp->O;
    int *iarr1 = 0, *rowmap = qslp->rowmap, *structmap = qslp->structmap, col;
    mpq_t *arr1 = 0, *arr2 = 0, *arr3 = 0, *arr4 = 0, *rhs_copy = 0;
    mpq_t *dz = 0;
    char *str1 = 0, *str2 = 0, *str3 = 0, *str4 = 0;
    int objsense = (qslp->objsense == QS_MIN) ? 1 : -1;
    int const msg_lvl = 100000 * (1 - p->simplex_display);
    int rval = 1;		/* store whether or not the solution is
				   optimal, we start assuming it is. */
    mpq_t num1, num2, num3, p_obj, d_obj;
    mpq_init (num1);
    mpq_init (num2);
    mpq_init (num3);
    mpq_init (p_obj);
    mpq_init (d_obj);
    mpq_set_ui (p_obj, 0UL, 1UL);
    mpq_set_ui (d_obj, 0UL, 1UL);

    /* now check if the given basis is the optimal basis */
    arr3 = qslp->lower;
    arr4 = qslp->upper;
    if (mpq_QSload_basis (p, basis)) {
	rval = 0;
	MESSAGE (msg_lvl, "QSload_basis failed");
	goto CLEANUP;
    }
    for (i = basis->nstruct; i--;) {
	/* check that the upper and lower bound define a non-empty space */
	if (mpq_cmp (arr3[structmap[i]], arr4[structmap[i]]) > 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (arr3[structmap[i]]);
	    str2 = mpq_EGlpNumGetStr (arr4[structmap[i]]);
	    MESSAGE (msg_lvl, "variable %s has empty feasible range [%s,%s]",
		qslp->colnames[i], str1, str2);
	    EGfree (str1);
	    EGfree (str2);
	    goto CLEANUP;
	}
	/* set the variable to its apropiate values, depending its status */
	switch (basis->cstat[i]) {
	case QS_COL_BSTAT_FREE:
	case QS_COL_BSTAT_BASIC:
	    if (mpq_cmp (p_sol[i], arr4[structmap[i]]) > 0)
		mpq_set (p_sol[i], arr4[structmap[i]]);
	    else if (mpq_cmp (p_sol[i], arr3[structmap[i]]) < 0)
		mpq_set (p_sol[i], arr3[structmap[i]]);
	    break;
	case QS_COL_BSTAT_UPPER:
	    mpq_set (p_sol[i], arr4[structmap[i]]);
	    break;
	case QS_COL_BSTAT_LOWER:
	    mpq_set (p_sol[i], arr3[structmap[i]]);
	    break;
	default:
	    rval = 0;
	    MESSAGE (msg_lvl, "Unknown Variable basic status %d, for variable "
		"(%s,%d)", basis->cstat[i], qslp->colnames[i], i);
	    goto CLEANUP;
	    break;
	}
    }
    for (i = basis->nrows; i--;) {
	/* check that the upper and lower bound define a non-empty space */
	if (mpq_cmp (arr3[rowmap[i]], arr4[rowmap[i]]) > 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (arr3[rowmap[i]]);
	    str2 = mpq_EGlpNumGetStr (arr4[rowmap[i]]);
	    MESSAGE (msg_lvl, "constraint %s logical has empty feasible range "
		"[%s,%s]", qslp->rownames[i], str1, str2);
	    EGfree (str1);
	    EGfree (str2);
	    goto CLEANUP;
	}
	/* set the variable to its apropiate values, depending its status */
	switch (basis->rstat[i]) {
	case QS_ROW_BSTAT_BASIC:
	    if (mpq_cmp (p_sol[i + basis->nstruct], arr4[rowmap[i]]) > 0)
		mpq_set (p_sol[i + basis->nstruct], arr4[rowmap[i]]);
	    else if (mpq_cmp (p_sol[i + basis->nstruct], arr3[rowmap[i]]) < 0)
		mpq_set (p_sol[i + basis->nstruct], arr3[rowmap[i]]);
	    break;
	case QS_ROW_BSTAT_UPPER:
	    mpq_set (p_sol[i + basis->nstruct], arr4[rowmap[i]]);
	    break;
	case QS_ROW_BSTAT_LOWER:
	    mpq_set (p_sol[i + basis->nstruct], arr3[rowmap[i]]);
	    break;
	default:
	    rval = 0;
	    MESSAGE (msg_lvl, "Unknown Variable basic status %d, for constraint "
		"(%s,%d)", basis->cstat[i], qslp->rownames[i], i);
	    goto CLEANUP;
	    break;
	}
    }

    /* compute the actual RHS */
    rhs_copy = mpq_EGlpNumAllocArray (qslp->nrows);
    for (i = qslp->nstruct; i--;) {
	if (!mpq_equal (p_sol[i], mpq_zeroLpNum)) {
	    arr1 = qslp->A.matval + qslp->A.matbeg[structmap[i]];
	    iarr1 = qslp->A.matind + qslp->A.matbeg[structmap[i]];
	    for (j = qslp->A.matcnt[structmap[i]]; j--;) {
		mpq_mul (num1, arr1[j], p_sol[i]);
		mpq_add (rhs_copy[iarr1[j]], rhs_copy[iarr1[j]], num1);
	    }
	}
    }

    /* now check if both rhs and copy_rhs are equal */
    arr4 = qslp->upper;
    arr1 = qslp->rhs;
    arr2 = qslp->lower;
    for (i = qslp->nrows; i--;) {
	mpq_mul (num1, arr1[i], d_sol[i]);
	mpq_add (d_obj, d_obj, num1);
	mpq_sub (num2, arr1[i], rhs_copy[i]);
	EXIT (qslp->A.matcnt[rowmap[i]] != 1, "Imposible!");
	if (basis->rstat[i] == QS_ROW_BSTAT_BASIC)
	    mpq_div (p_sol[qslp->nstruct + i], num2,
		qslp->A.matval[qslp->A.matbeg[rowmap[i]]]);
	else {
	    mpq_mul (num1, p_sol[qslp->nstruct + i],
		qslp->A.matval[qslp->A.matbeg[rowmap[i]]]);
	    if (!mpq_equal (num1, num2)) {
		rval = 0;
		MESSAGE (msg_lvl, "solution is infeasible for constraint %s, violation"
		    " %lg", qslp->rownames[i],
		    mpq_get_d (num1) - mpq_get_d (num2));
		goto CLEANUP;
	    }
	}
	mpq_set (num2, p_sol[qslp->nstruct + i]);
	/* now we check the bounds on the logical variables */
	if (mpq_cmp (num2, arr2[rowmap[i]]) < 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (num2);
	    str2 = mpq_EGlpNumGetStr (arr2[rowmap[i]]);
	    str3 = mpq_EGlpNumGetStr (rhs_copy[i]);
	    str4 = mpq_EGlpNumGetStr (arr1[i]);
	    MESSAGE (msg_lvl, "constraint %s artificial (%lg,%s) bellow lower"
		" bound (%lg,%s), actual LHS (%lg,%s), actual RHS (%lg,%s)",
		qslp->rownames[i], mpq_get_d (num2), str1,
		mpq_get_d (arr2[rowmap[i]]), str2, mpq_get_d (rhs_copy[i]), str3,
		mpq_get_d (arr1[i]), str4);
	    free (str1);
	    free (str2);
	    free (str3);
	    free (str4);
	    goto CLEANUP;
	} else if (mpq_cmp (num2, arr4[rowmap[i]]) > 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (num2);
	    str2 = mpq_EGlpNumGetStr (arr4[rowmap[i]]);
	    MESSAGE (msg_lvl, "constraint %s artificial (%lg,%s) bellow lower bound"
		" (%lg,%s)", qslp->rownames[i], mpq_get_d (num2), str1,
		mpq_get_d (arr4[rowmap[i]]), str2);
	    free (str1);
	    free (str2);
	    goto CLEANUP;
	}
    }

    /* compute the upper and lower bound dual variables, note that dl is the
       dual of the lower bounds, and du the dual of the upper bound, dl >= 0
       and du <= 0 and A^t y + Idl + Idu = c, and the dual objective value is
       max y*b + l*dl + u*du, we colapse both vector dl and du into dz, note
       that if we are maximizing, then dl <= 0 and du >=0 */
    dz = mpq_EGlpNumAllocArray (qslp->ncols);
    arr2 = qslp->obj;
    arr3 = qslp->lower;
    arr4 = qslp->upper;
    for (i = qslp->nstruct; i--;) {
	col = structmap[i];
	mpq_mul (num1, arr2[col], p_sol[i]);
	mpq_add (p_obj, p_obj, num1);
	arr1 = qslp->A.matval + qslp->A.matbeg[col];
	iarr1 = qslp->A.matind + qslp->A.matbeg[col];
	mpq_set (num1, arr2[col]);
	for (j = qslp->A.matcnt[col]; j--;) {
	    mpq_mul (num2, arr1[j], d_sol[iarr1[j]]);
	    mpq_sub (num1, num1, num2);
	}
	mpq_set (dz[col], num1);
	/* objective update */
	if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0) {
	    mpq_mul (num3, dz[col], arr3[col]);
	    mpq_add (d_obj, d_obj, num3);
	} else {
	    mpq_mul (num3, dz[col], arr4[col]);
	    mpq_add (d_obj, d_obj, num3);
	}
	/* now we check that only when the logical is tight then the dual
	   variable may be non-zero, also check for primal feasibility with
	   respect to lower/upper bounds. */
	mpq_set_ui (num2, 0UL, 1UL);
	if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0) {
	    mpq_sub (num1, p_sol[i], arr3[col]);
	    mpq_mul (num2, num1, dz[col]);
	}
	if (mpq_cmp_ui (num2, 0UL, 1UL) != 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (num1);
	    str2 = mpq_EGlpNumGetStr (dz[col]);
	    MESSAGE (msg_lvl, "lower bound (%s,%d) slack (%s) and dual variable (%s)"
		" don't satisfy complementary slacknes %s", qslp->colnames[i],
		i, str1, str2, "(real)");
	    free (str1);
	    free (str2);
	    goto CLEANUP;
	}
	mpq_set_ui (num2, 0UL, 1UL);
	if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) < 0) {
	    mpq_sub (num1, p_sol[i], arr4[col]);
	    mpq_mul (num2, num1, dz[col]);
	}
	if (mpq_cmp_ui (num2, 0UL, 1UL) != 0) {
	    rval = 0;
	    str3 = mpq_EGlpNumGetStr (arr4[col]);
	    str1 = mpq_EGlpNumGetStr (p_sol[i]);
	    str2 = mpq_EGlpNumGetStr (dz[col]);
	    MESSAGE (msg_lvl, "upper bound (%s) variable (%s) and dual variable (%s) "
		"don't satisfy complementary slacknes for variable (%s,%d) %s",
		str3, str1, str2, qslp->colnames[i], i, "(real)");
	    free (str1);
	    free (str2);
	    free (str3);
	    goto CLEANUP;
	}
    }
    /* complenetary slackness checked, now update the same for the logical
       variables */
    for (i = qslp->nrows; i--;) {
	col = rowmap[i];
	mpq_mul (num1, arr2[col], p_sol[i + qslp->nstruct]);
	WARNING (mpq_cmp (arr2[col], mpq_zeroLpNum), "logical variable %s with "
	    "non-zero objective function %lf", qslp->rownames[i],
	    mpq_get_d (arr2[col]));
	mpq_add (p_obj, p_obj, num1);
	arr1 = qslp->A.matval + qslp->A.matbeg[col];
	iarr1 = qslp->A.matind + qslp->A.matbeg[col];
	mpq_set (num1, arr2[col]);
	for (j = qslp->A.matcnt[col]; j--;) {
	    mpq_mul (num2, arr1[j], d_sol[iarr1[j]]);
	    mpq_sub (num1, num1, num2);
	}
	mpq_set (dz[col], num1);
	/* objective update */
	if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0) {
	    mpq_mul (num3, dz[col], arr3[col]);
	    mpq_add (d_obj, d_obj, num3);
	} else {
	    mpq_mul (num3, dz[col], arr4[col]);
	    mpq_add (d_obj, d_obj, num3);
	}
	/* now we check that only when the primal variable is tight then the
	   dual variable may be non-zero, also check for primal feasibility
	   with respect to lower/upper bounds. */
	mpq_set_ui (num2, 0UL, 1UL);
	if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) > 0) {
	    mpq_sub (num1, p_sol[i + qslp->nstruct], arr3[col]);
	    mpq_mul (num2, num1, dz[col]);
	}
	if (mpq_cmp_ui (num2, 0UL, 1UL) != 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (num1);
	    str2 = mpq_EGlpNumGetStr (dz[col]);
	    MESSAGE (msg_lvl, "lower bound (%s,%d) slack (%s) and dual variable (%s)"
		" don't satisfy complementary slacknes %s", qslp->colnames[col],
		i, str1, str2, "(real)");
	    free (str1);
	    free (str2);
	    goto CLEANUP;
	}
	mpq_set_ui (num2, 0UL, 1UL);
	if (objsense * mpq_cmp_ui (dz[col], 0UL, 1UL) < 0) {
	    mpq_sub (num1, p_sol[i + qslp->nstruct], arr4[col]);
	    mpq_mul (num2, num1, dz[col]);
	}
	if (mpq_cmp_ui (num2, 0UL, 1UL) != 0) {
	    rval = 0;
	    str3 = mpq_EGlpNumGetStr (arr4[col]);
	    str1 = mpq_EGlpNumGetStr (p_sol[i + qslp->nstruct]);
	    str2 = mpq_EGlpNumGetStr (dz[col]);
	    MESSAGE (msg_lvl, "upper bound (%s) variable (%s) and dual variable (%s) "
		"don't satisfy complementary slacknes for variable (%s,%d) %s",
		str3, str1, str2, qslp->colnames[col], i, "(real)");
	    free (str1);
	    free (str2);
	    free (str3);
	    goto CLEANUP;
	}
    }

    /* now check the objective values */
    if (mpq_cmp (p_obj, d_obj) != 0) {
	rval = 0;
	str1 = mpq_EGlpNumGetStr (p_obj);
	str2 = mpq_EGlpNumGetStr (d_obj);
	MESSAGE (msg_lvl, "primal and dual objective value differ %s %s", str1,
	    str2);
	free (str1);
	free (str2);
	goto CLEANUP;
    }
    str1 = mpq_EGlpNumGetStr (p_obj);
    MESSAGE (msg_lvl, "Problem solved to optimality, LP value %s", str1);
    free (str1);

    /* now we load into cache the solution */
    if (!p->cache) {
	p->cache = EGsMalloc (mpq_ILLlp_cache, 1);
	EGlpNumInitVar (p->cache->val);
	mpq_ILLlp_cache_init (p->cache);
    }
    if (qslp->nrows != p->cache->nrows || qslp->nstruct != p->cache->nstruct) {
	mpq_ILLlp_cache_free (p->cache);
	rval = mpq_ILLlp_cache_alloc (p->cache, qslp->nstruct, qslp->nrows);
	if (rval)
	    EXIT (1, "mpq_ILLlp_cache_alloc failed");
    }
    p->cache->status = QS_LP_OPTIMAL;
    p->qstatus = QS_LP_OPTIMAL;
    p->lp->basisstat.optimal = 1;
    mpq_set (p->cache->val, p_obj);
    for (i = qslp->nstruct; i--;) {
	mpq_set (p->cache->x[i], p_sol[i]);
	mpq_set (p->cache->rc[i], dz[structmap[i]]);
    }
    for (i = qslp->nrows; i--;) {
	mpq_set (p->cache->slack[i], p_sol[i + qslp->nstruct]);
	mpq_set (p->cache->pi[i], d_sol[i]);
    }

    /* save the problem and solution if enablred */
#if QSEXACT_SAVE_OPTIMAL
    {
	char stmp[1024];
	FILE *out_f = 0;
	snprintf (stmp, 1023, "%s-opt%03d.lp", p->name ? p->name : "UNNAMED",
	    QSEXACT_SAVE_OPTIMAL_IND);
	if (mpq_QSwrite_prob (p, stmp, "LP")) {
	    rval = 0;
	    MESSAGE (0, "Couldn't write output problem %s", stmp);
	    goto CLEANUP;
	}
	snprintf (stmp, 1023, "%s-opt%03d.sol", p->name ? p->name : "UNNAMED",
	    QSEXACT_SAVE_OPTIMAL_IND);
	if (!(out_f = fopen (stmp, "w+"))) {
	    rval = 0;
	    MESSAGE (0, "Couldn't open solution file %s", stmp);
	    goto CLEANUP;
	}
	if (QSexact_print_sol (p, out_f)) {
	    rval = 0;
	    MESSAGE (0, "Couldn't write output solution %s", stmp);
	    goto CLEANUP;
	}
	fclose (out_f);
	QSEXACT_SAVE_OPTIMAL_IND++;
    }
#endif
    rval = 1;

    /* ending */
CLEANUP:
    mpq_EGlpNumFreeArray (dz);
    mpq_EGlpNumFreeArray (rhs_copy);
    mpq_clear (num1);
    mpq_clear (num2);
    mpq_clear (num3);
    mpq_clear (p_obj);
    mpq_clear (d_obj);
    return rval;
}

/* ========================================================================= */
int QSexact_infeasible_test (mpq_QSdata * p, mpq_t * d_sol)
{
    /* local variables */
    register int i, j;
    int *iarr1;
    mpq_ILLlpdata *qslp = p->lp->O;
    mpq_t *arr1, *arr2, *arr3, *arr4;
    mpq_t *dl = 0, *du = 0;
    char *str1;
    int const msg_lvl = (1 - p->simplex_display) * 10000;
    int rval = 1;		/* store whether or not the solution is
				   optimal, we start assuming it is. */
    mpq_t num1, num2, num3, d_obj;
    mpq_init (num1);
    mpq_init (num2);
    mpq_init (num3);
    mpq_init (d_obj);
    mpq_set_ui (d_obj, 0UL, 1UL);

    /* compute the dual objective value */
    arr2 = qslp->rhs;
    for (i = qslp->nrows; i--;) {
	mpq_mul (num1, arr2[i], d_sol[i]);
	mpq_add (d_obj, d_obj, num1);
    }

    /* compute the upper and lower bound dual variables, note that dl is the
       dual of the lower bounds, and du the dual of the upper bound, dl <= 0
       and du >= 0 and A^t y + Idl + Idu = c, and the dual objective value is
       max y*b + l*dl + u*du */
    du = mpq_EGlpNumAllocArray (qslp->ncols);
    dl = mpq_EGlpNumAllocArray (qslp->ncols);
    arr3 = qslp->lower;
    arr4 = qslp->upper;
    for (i = qslp->ncols; i--;) {
	arr1 = qslp->A.matval + qslp->A.matbeg[i];
	iarr1 = qslp->A.matind + qslp->A.matbeg[i];
	mpq_set_ui (num1, 0UL, 1UL);
	mpq_set_ui (du[i], 0UL, 1UL);
	mpq_set_ui (dl[i], 0UL, 1UL);
	for (j = qslp->A.matcnt[i]; j--;) {
	    mpq_mul (num2, arr1[j], d_sol[iarr1[j]]);
	    mpq_sub (num1, num1, num2);
	}
	if (mpq_cmp_ui (num1, 0UL, 1UL) < 0)
	    mpq_set (du[i], num1);
	else
	    mpq_set (dl[i], num1);
	if (mpq_equal (arr4[i], mpq_ILL_MAXDOUBLE) &&
	    mpq_cmp_ui (du[i], 0UL, 1UL) != 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (du[i]);
	    MESSAGE (msg_lvl, "upper bound of variable is INFTY, and it's dual is "
		"non-zero %s", str1);
	    free (str1);
	    goto CLEANUP;
	}
	if (mpq_equal (arr3[i], mpq_ILL_MINDOUBLE) &&
	    mpq_cmp_ui (dl[i], 0UL, 1UL) != 0) {
	    rval = 0;
	    str1 = mpq_EGlpNumGetStr (dl[i]);
	    MESSAGE (msg_lvl, "lower bound of variable is -INFTY, and it's dual is "
		"non-zero %s", str1);
	    free (str1);
	    goto CLEANUP;
	}
	mpq_mul (num3, dl[i], arr3[i]);
	mpq_add (d_obj, d_obj, num3);
	mpq_mul (num3, du[i], arr4[i]);
	mpq_add (d_obj, d_obj, num3);
    }
    /* now check the objective values */
    if (mpq_cmp_ui (d_obj, 0UL, 1UL) <= 0) {
	rval = 0;
	str1 = mpq_EGlpNumGetStr (d_obj);
	MESSAGE (msg_lvl, "dual ray is feasible, but objective is non positive %s",
	    str1);
	free (str1);
	goto CLEANUP;
    }
    p->qstatus = QS_LP_INFEASIBLE;

    /* ending */
CLEANUP:
    mpq_EGlpNumFreeArray (dl);
    mpq_EGlpNumFreeArray (du);
    mpq_clear (num1);
    mpq_clear (num2);
    mpq_clear (num3);
    mpq_clear (d_obj);
    return rval;
}

/* ========================================================================= */
/* Used as separator while printing output to the screen (controled by
 * enabling simplex_display in the mpq_QSdata */
/* ========================================================================= */
static const char __sp[81] =
"================================================================================";

/* ========================================================================= */
/* print into screen (if enable) a message indicating that we have
 * successfully prove infeasibility, and save (if y is non
 * NULL ) the dual ray solution provided in y_mpq.
 * @param p_mpq the problem data.
 * @param y where to store the optimal dual solution (if not null).
 * @param y_mpq  the optimal dual solution.
 * */
/* ========================================================================= */
static void infeasible_output (mpq_QSdata * p_mpq, mpq_t * const y,
      mpq_t * y_mpq)
{
    if (p_mpq->simplex_display)
	fprintf (stderr, "%s\n\tProblem Is Infeasible\n%s\n", __sp, __sp);
    if (y) {
	unsigned sz = __EGlpNumArraySize (y_mpq);
	while (sz--)
	    mpq_set (y[sz], y_mpq[sz]);
    }
}

/* ========================================================================= */
/* print into screen (if enable) a message indicating that we have
 * successfully solved the problem at optimality, and save (if x and y are non
 * NULL respectivelly) the optimal primal/dual solution provided in x_mpq and
 * y_mpq.
 * param p_mpq the problem data.
 * param x where to store the optimal primal solution (if not null).
 * param y where to store the optimal dual solution (if not null).
 * param x_mpq  the optimal primal solution.
 * param y_mpq  the optimal dual solution.
 * */
/* ========================================================================= */
static void optimal_output (mpq_QSdata * p_mpq, mpq_t * const x,
      mpq_t * const y, mpq_t * x_mpq, mpq_t * y_mpq)
{
    if (p_mpq->simplex_display)
	fprintf (stderr, "%s\n\tProblem Solved Exactly\n%s\n", __sp, __sp);
    if (y) {
	unsigned sz = __EGlpNumArraySize (y_mpq);
	while (sz--)
	    mpq_set (y[sz], y_mpq[sz]);
    }
    if (x) {
	unsigned sz = __EGlpNumArraySize (x_mpq);
	while (sz--)
	    mpq_set (x[sz], x_mpq[sz]);
    }
}

/* ========================================================================= */
/* get the status for a given basis in rational arithmetic, it should
 * also leave everything set to get primal/dual solutions when needed.
 * */
static int QSexact_basis_status (mpq_QSdata * p_mpq, int *status,
      QSbasis * const basis, const int msg_lvl, int *const simplexalgo)
{
    int rval = 0, singular;
    mpq_feas_info fi;
    double szeit;

    mpq_EGlpNumInitVar (fi.totinfeas);
    szeit = ILLutil_zeit ();

    rval = mpq_QSload_basis (p_mpq, basis);
    CHECKRVALG (rval, CLEANUP);
    if (p_mpq->cache) {
	mpq_ILLlp_cache_free (p_mpq->cache);
	mpq_clear (p_mpq->cache->val);
	ILL_IFFREE (p_mpq->cache, mpq_ILLlp_cache);
    }
    p_mpq->qstatus = QS_LP_MODIFIED;
    if (p_mpq->qslp->sinfo) {
	mpq_ILLlp_sinfo_free (p_mpq->qslp->sinfo);
	ILL_IFFREE (p_mpq->qslp->sinfo, mpq_ILLlp_sinfo);
    }
    if (p_mpq->qslp->rA) {
	mpq_ILLlp_rows_clear (p_mpq->qslp->rA);
	ILL_IFFREE (p_mpq->qslp->rA, mpq_ILLlp_rows);
    }
    mpq_free_internal_lpinfo (p_mpq->lp);
    mpq_init_internal_lpinfo (p_mpq->lp);
    rval = mpq_build_internal_lpinfo (p_mpq->lp);
    CHECKRVALG (rval, CLEANUP);
    mpq_ILLfct_set_variable_type (p_mpq->lp);
    rval = mpq_ILLbasis_load (p_mpq->lp, p_mpq->basis);
    CHECKRVALG (rval, CLEANUP);
    rval = mpq_ILLbasis_factor (p_mpq->lp, &singular);
    CHECKRVALG (rval, CLEANUP);
    memset (&(p_mpq->lp->basisstat), 0, sizeof (mpq_lp_status_info));
    mpq_ILLfct_compute_piz (p_mpq->lp);
    mpq_ILLfct_compute_dz (p_mpq->lp);
    mpq_ILLfct_compute_xbz (p_mpq->lp);
    mpq_ILLfct_check_pfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
    mpq_ILLfct_check_dfeasible (p_mpq->lp, &fi, mpq_zeroLpNum);
    mpq_ILLfct_set_status_values (p_mpq->lp, fi.pstatus, fi.dstatus, PHASEII,
	PHASEII);
    if (p_mpq->lp->basisstat.optimal) {
	*status = QS_LP_OPTIMAL;
	rval = mpq_grab_cache (p_mpq, QS_LP_OPTIMAL);
	CHECKRVALG (rval, CLEANUP);
    } else if (p_mpq->lp->basisstat.primal_infeasible
	|| p_mpq->lp->basisstat.dual_unbounded) {
	if (*status == QS_LP_INFEASIBLE)
	    *simplexalgo = PRIMAL_SIMPLEX;
	*status = QS_LP_INFEASIBLE;
	p_mpq->lp->final_phase = PRIMAL_PHASEI;
	p_mpq->lp->pIpiz = mpq_EGlpNumAllocArray (p_mpq->lp->nrows);
	mpq_ILLfct_compute_phaseI_piz (p_mpq->lp);
    } else if (p_mpq->lp->basisstat.primal_unbounded)
	*status = QS_LP_UNBOUNDED;
    else
	*status = QS_LP_UNSOLVED;
    MESSAGE (msg_lvl, "Performing Rational Basic Solve on %s, %s, check"
	" done in %lg seconds, PS %s, DS %s", p_mpq->name, (*status ==
	    QS_LP_OPTIMAL) ? "RAT_optimal" : ((*status ==
		QS_LP_INFEASIBLE) ? "RAT_infeasible" : ((*status ==
		    QS_LP_UNBOUNDED) ? "RAT_unbounded" : "RAT_unsolved")),
	ILLutil_zeit () - szeit, p_mpq->lp->basisstat.primal_feasible ?
	"F" : (p_mpq->lp->basisstat.primal_infeasible ? "I" : "U"),
	p_mpq->lp->basisstat.dual_feasible ?
	"F" : (p_mpq->lp->basisstat.dual_infeasible ? "I" : "U"));
CLEANUP:
    mpq_EGlpNumClearVar (fi.totinfeas);
    return rval;
}

/* ========================================================================= */
int QSexact_solver (mpq_QSdata * p_mpq, mpq_t * const x, mpq_t * const y,
      QSbasis * const ebasis, int simplexalgo, int *status)
{
    /* local variables */
    QSbasis *basis = 0;
    unsigned precision = EGLPNUM_PRECISION;
    int rval = 0, it = QS_EXACT_MAX_ITER;
    dbl_QSdata *p_dbl = 0;
    mpf_QSdata *p_mpf = 0;
    double *x_dbl = 0, *y_dbl = 0;
    mpq_t *x_mpq = 0, *y_mpq = 0;
    mpf_t *x_mpf = 0, *y_mpf = 0;
    int const msg_lvl = (1 - p_mpq->simplex_display) * 10000;
    *status = 0;

    /* try first with doubles */
    if (p_mpq->simplex_display)
	fprintf (stderr, "%s\n\tTrying double precision\n%s\n", __sp, __sp);
    p_dbl = QScopy_prob_mpq_dbl (p_mpq, "dbl_problem");
    if (ebasis && ebasis->nstruct)
	dbl_QSload_basis (p_dbl, ebasis);
    rval = dbl_ILLeditor_solve (p_dbl, simplexalgo);
    if (rval) {
	if (p_mpq->simplex_display)
	    fprintf (stderr,
		"double approximation failed, code %d, continuing in extended"
		" precision\n", rval);
	rval = 0;
	goto MPF_PRECISION;
    }
    rval = dbl_QSget_status (p_dbl, status);
    CHECKRVALG (rval, CLEANUP);
    if ((*status == QS_LP_INFEASIBLE) &&
	(p_dbl->lp->final_phase != PRIMAL_PHASEI) &&
	(p_dbl->lp->final_phase != DUAL_PHASEII))
	dbl_QSopt_primal (p_dbl, status);
    rval = dbl_QSget_status (p_dbl, status);
    CHECKRVALG (rval, CLEANUP);
    /* deal with the problem depending on status we get from our optimizer */
    switch (*status) {
    case QS_LP_OPTIMAL:
	x_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->ncols);
	y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
	rval = dbl_QSget_x_array (p_dbl, x_dbl);
	CHECKRVALG (rval, CLEANUP);
	rval = dbl_QSget_pi_array (p_dbl, y_dbl);
	CHECKRVALG (rval, CLEANUP);
	x_mpq = QScopy_array_dbl_mpq (x_dbl);
	y_mpq = QScopy_array_dbl_mpq (y_dbl);
	dbl_EGlpNumFreeArray (x_dbl);
	dbl_EGlpNumFreeArray (y_dbl);
	basis = dbl_QSget_basis (p_dbl);
	rval = QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis);
	if (!rval) {
	    rval = QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo);
	    CHECKRVALG (rval, CLEANUP);
	    if (*status == QS_LP_OPTIMAL) {
		MESSAGE (msg_lvl, "Retesting solution");
		rval = mpq_QSget_x_array (p_mpq, x_mpq);
		CHECKRVALG (rval, CLEANUP);
		rval = mpq_QSget_pi_array (p_mpq, y_mpq);
		CHECKRVALG (rval, CLEANUP);
		rval = QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis);
		if (rval) {
		    rval = 0;
		    optimal_output (p_mpq, x, y, x_mpq, y_mpq);
		    goto CLEANUP;
		}
	    } else
		MESSAGE (msg_lvl, "Status is not optimal, but %d", *status);
	}
	if (rval) {
	    rval = 0;
	    optimal_output (p_mpq, x, y, x_mpq, y_mpq);
	    goto CLEANUP;
	}
	mpq_EGlpNumFreeArray (x_mpq);
	mpq_EGlpNumFreeArray (y_mpq);
	break;
    case QS_LP_INFEASIBLE:
	y_dbl = dbl_EGlpNumAllocArray (p_dbl->qslp->nrows);
	rval = dbl_QSget_infeas_array (p_dbl, y_dbl);
	if (rval) {
	    if (p_mpq->simplex_display)
		fprintf (stderr, "double approximation failed, code %d, continuing "
		    "in extended precision\n", rval);
	    rval = 0;
	    goto MPF_PRECISION;
	}
	y_mpq = QScopy_array_dbl_mpq (y_dbl);
	dbl_EGlpNumFreeArray (y_dbl);
	rval = QSexact_infeasible_test (p_mpq, y_mpq);
	if (rval) {
	    rval = 0;
	    infeasible_output (p_mpq, y, y_mpq);
	    goto CLEANUP;
	} else {
	    MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
	    basis = dbl_QSget_basis (p_dbl);
	    rval = QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo);
	    CHECKRVALG (rval, CLEANUP);
#if 0
	    mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_MAX_ITERATIONS, 1);
	    mpq_QSload_basis (p_mpq, basis);
	    mpq_QSfree_basis (basis);
	    rval = mpq_ILLeditor_solve (p_mpq, simplexalgo);
	    CHECKRVALG (rval, CLEANUP);
	    rval = mpq_QSget_status (p_mpq, status);
	    CHECKRVALG (rval, CLEANUP);
#endif
	    if (*status == QS_LP_INFEASIBLE) {
		mpq_EGlpNumFreeArray (y_mpq);
		y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
		rval = mpq_QSget_infeas_array (p_mpq, y_mpq);
		CHECKRVALG (rval, CLEANUP);
		rval = QSexact_infeasible_test (p_mpq, y_mpq);
		if (rval) {
		    infeasible_output (p_mpq, y, y_mpq);
		    rval = 0;
		    goto CLEANUP;
		}
	    }
	}
	mpq_EGlpNumFreeArray (y_mpq);
	break;
    case QS_LP_UNBOUNDED:
	if (p_mpq->simplex_display)
	    fprintf (stderr, "%s\n\tUnbounded Problem found, not implemented "
		"to deal with this\n%s\n", __sp, __sp);
	break;
    default:
	break;
    }
    /* if we reach this point, then we have to keep going, we use the
       previous basis ONLY if the previous precision think that it has the
       optimal solution, otherwise we start from scratch. */
MPF_PRECISION:
    dbl_QSfree_prob (p_dbl);
    p_dbl = 0;
    /* try with multiple precision floating points */
    for (; it--; precision = (unsigned) (precision * 1.5)) {
	QSexact_set_precision (precision);
	if (p_mpq->simplex_display)
	    fprintf (stderr, "%s\n\tTrying mpf with %u bits\n%s\n", __sp, precision,
		__sp);
	p_mpf = QScopy_prob_mpq_mpf (p_mpq, "mpf_problem");
	if (basis) {
	    rval = mpf_QSload_basis (p_mpf, basis);
	    mpf_QSfree_basis (basis);
	    basis = 0;
	} else if (ebasis && ebasis->nstruct)
	    mpf_QSload_basis (p_mpf, ebasis);
	rval = mpf_ILLeditor_solve (p_mpf, simplexalgo);
	if (rval) {
	    if (p_mpq->simplex_display)
		fprintf (stderr,
		    "mpf_%u precision falied, error code %d, continuing with "
		    "next precision", precision, rval);
	    rval = 0;
	    goto NEXT_PRECISION;
	}
	rval = mpf_QSget_status (p_mpf, status);
	CHECKRVALG (rval, CLEANUP);
	if ((*status == QS_LP_INFEASIBLE) &&
	    (p_mpf->lp->final_phase != PRIMAL_PHASEI) &&
	    (p_mpf->lp->final_phase != DUAL_PHASEII))
	    mpf_QSopt_primal (p_mpf, status);
	rval = mpf_QSget_status (p_mpf, status);
	CHECKRVALG (rval, CLEANUP);
	/* deal with the problem depending on status we got from our
	   optimizer */
	switch (*status) {
	case QS_LP_OPTIMAL:
	    basis = mpf_QSget_basis (p_mpf);
	    x_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->ncols);
	    y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
	    rval = mpf_QSget_x_array (p_mpf, x_mpf);
	    CHECKRVALG (rval, CLEANUP);
	    rval = mpf_QSget_pi_array (p_mpf, y_mpf);
	    CHECKRVALG (rval, CLEANUP);
	    x_mpq = QScopy_array_mpf_mpq (x_mpf);
	    y_mpq = QScopy_array_mpf_mpq (y_mpf);
	    mpf_EGlpNumFreeArray (x_mpf);
	    mpf_EGlpNumFreeArray (y_mpf);
	    rval = QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis);
	    if (!rval) {
		rval =
		    QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo);
		CHECKRVALG (rval, CLEANUP);
		if (*status == QS_LP_OPTIMAL) {
		    MESSAGE (msg_lvl, "Retesting solution");
		    rval = mpq_QSget_x_array (p_mpq, x_mpq);
		    CHECKRVALG (rval, CLEANUP);
		    rval = mpq_QSget_pi_array (p_mpq, y_mpq);
		    CHECKRVALG (rval, CLEANUP);
		    rval = QSexact_optimal_test (p_mpq, x_mpq, y_mpq, basis);
		    if (rval) {
			rval = 0;
			optimal_output (p_mpq, x, y, x_mpq, y_mpq);
			goto CLEANUP;
		    }
		} else
		    MESSAGE (msg_lvl, "Status is not optimal, but %d", *status);
	    }
	    if (rval) {
		rval = 0;
		optimal_output (p_mpq, x, y, x_mpq, y_mpq);
		goto CLEANUP;
	    }
	    mpq_EGlpNumFreeArray (x_mpq);
	    mpq_EGlpNumFreeArray (y_mpq);
	    break;
	case QS_LP_INFEASIBLE:
	    y_mpf = mpf_EGlpNumAllocArray (p_mpf->qslp->nrows);
	    rval = mpf_QSget_infeas_array (p_mpf, y_mpf);
	    CHECKRVALG (rval, CLEANUP);
	    y_mpq = QScopy_array_mpf_mpq (y_mpf);
	    mpf_EGlpNumFreeArray (y_mpf);
	    rval = QSexact_infeasible_test (p_mpq, y_mpq);
	    if (!rval) {
		MESSAGE (msg_lvl, "Retesting solution in exact arithmetic");
		basis = mpf_QSget_basis (p_mpf);
		rval =
		    QSexact_basis_status (p_mpq, status, basis, msg_lvl, &simplexalgo);
		CHECKRVALG (rval, CLEANUP);
#if 0
		mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_MAX_ITERATIONS, 1);
		mpq_QSload_basis (p_mpq, basis);
		mpq_QSfree_basis (basis);
		rval = mpq_ILLeditor_solve (p_mpq, simplexalgo);
		CHECKRVALG (rval, CLEANUP);
		rval = mpq_QSget_status (p_mpq, status);
		CHECKRVALG (rval, CLEANUP);
#endif
		if (*status == QS_LP_INFEASIBLE) {
		    mpq_EGlpNumFreeArray (y_mpq);
		    y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
		    rval = mpq_QSget_infeas_array (p_mpq, y_mpq);
		    CHECKRVALG (rval, CLEANUP);
		    rval = QSexact_infeasible_test (p_mpq, y_mpq);
		    if (rval) {
			infeasible_output (p_mpq, y, y_mpq);
			rval = 0;
			goto CLEANUP;
		    }
		}
	    }
	    if (rval) {
		infeasible_output (p_mpq, y, y_mpq);
		rval = 0;
		goto CLEANUP;
	    }
	    mpq_EGlpNumFreeArray (y_mpq);
	    break;
	case QS_LP_UNBOUNDED:
	    if (p_mpq->simplex_display)
		fprintf (stderr, "%s\n\tUnbounded Problem found, not implemented "
		    "to deal with this\n%s\n", __sp, __sp);
	    goto CLEANUP;
	    break;
	    break;
	default:
	    break;
	}
NEXT_PRECISION:
	mpf_QSfree_prob (p_mpf);
	p_mpf = 0;
    }
    /* ending */
CLEANUP:
    dbl_EGlpNumFreeArray (x_dbl);
    dbl_EGlpNumFreeArray (y_dbl);
    mpq_EGlpNumFreeArray (x_mpq);
    mpq_EGlpNumFreeArray (y_mpq);
    mpf_EGlpNumFreeArray (x_mpf);
    mpf_EGlpNumFreeArray (y_mpf);
    if (ebasis && basis) {
	ILL_IFFREE (ebasis->cstat, char);
	ILL_IFFREE (ebasis->rstat, char);
	ebasis->nstruct = basis->nstruct;
	ebasis->nrows = basis->nrows;
	ebasis->cstat = basis->cstat;
	ebasis->rstat = basis->rstat;
	basis->cstat = basis->rstat = 0;
    }
    mpq_QSfree_basis (basis);
    dbl_QSfree_prob (p_dbl);
    mpf_QSfree_prob (p_mpf);
    return rval;
}

