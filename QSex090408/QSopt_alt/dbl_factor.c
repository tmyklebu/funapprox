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

/* RCS_INFO = "$RCSfile: dbl_factor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/* implement a = max(a,abs(b)) and execute the extra code if the update is
   needed */
#define dbl_EGlpNumSetToMaxAbsAndDo(a,b,c) \
	if(dbl_EGlpNumIsLess(dbl_zeroLpNum,b))\
	{\
		if(dbl_EGlpNumIsLess(a,b)){\
			dbl_EGlpNumCopy(a,b);\
			c;\
			}\
	}\
	else\
	{\
		dbl_EGlpNumSign(a);\
		if(dbl_EGlpNumIsLess(b,a)){\
			dbl_EGlpNumCopy(a,b);\
			c;\
			}\
		dbl_EGlpNumSign(a);\
	}


#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdefs.h"
#include "dbl_factor.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


#undef  dbl_RECORD
#undef  dbl_DEBUG_FACTOR
#undef  dbl_SOLVE_DEBUG

#undef  dbl_FACTOR_DEBUG
#undef  dbl_UPDATE_DEBUG

#undef dbl_TRACK_FACTOR
#undef dbl_NOTICE_BLOWUP

#undef  dbl_FACTOR_STATS
#undef  dbl_UPDATE_STATS
#undef  dbl_GROWTH_STATS

#undef  dbl_UPDATE_STUDY

#undef  dbl_SORT_RESULTS

#ifdef dbl_UPDATE_STUDY
int dbl_nupdate = 0;
long int dbl_colspiketot = 0.0;
long int dbl_rowspiketot = 0.0;
long int dbl_permshifttot = 0.0;
long int dbl_leftetatot = 0.0;
#endif

void dbl_ILLfactor_init_factor_work (dbl_factor_work * f)
{
    f->max_k = 1000;		/* must be less than 46340 (2^15.5) */
    dbl_EGlpNumCopy (f->fzero_tol, dbl_SZERO_TOLER);	/* 2^-50 */
    dbl_EGlpNumCopy (f->szero_tol, dbl_SZERO_TOLER);	/* 2^-50 */
    dbl_EGlpNumCopy (f->partial_tol, dbl_OBJBND_TOLER);	/* 2^-7 */
    f->ur_space_mul = 2.0;
    f->uc_space_mul = 1.1;
    f->lc_space_mul = 1.1;
    f->er_space_mul = 1000.0;
    f->grow_mul = 1.5;
    f->p = 4;
    f->etamax = 100;
    f->minmult = 1e3;
    f->maxmult = 1e5;
    f->updmaxmult = 1e7;
    f->dense_fract = 0.25;
    f->dense_min = 25;
    dbl_EGlpNumCopy (f->partial_cur, f->partial_tol);
    f->work_coef = 0;
    f->work_indx = 0;
    f->uc_inf = 0;
    f->ur_inf = 0;
    f->lc_inf = 0;
    f->lr_inf = 0;
    f->er_inf = 0;
    f->ucindx = 0;
    f->ucrind = 0;
    f->uccoef = 0;
    f->urindx = 0;
    f->urcind = 0;
    f->urcoef = 0;
    f->lcindx = 0;
    f->lccoef = 0;
    f->lrindx = 0;
    f->lrcoef = 0;
    f->erindx = 0;
    f->ercoef = 0;
    f->rperm = 0;
    f->rrank = 0;
    f->cperm = 0;
    f->crank = 0;
    f->dmat = 0;
    dbl_ILLsvector_init (&f->xtmp);
}

void dbl_ILLfactor_free_factor_work (dbl_factor_work * f)
{
#ifdef dbl_UPDATE_STUDY
    if (dbl_nupdate) {
	printf
	("UPDATE STUDY: avg %d upd: %.2f col, %.2f row, %.2f lefteta, %.2f perm\n",
	    dbl_nupdate, ((double) dbl_colspiketot) / dbl_nupdate,
	    ((double) dbl_rowspiketot) / dbl_nupdate, ((double) dbl_leftetatot) / dbl_nupdate,
	    ((double) dbl_permshifttot) / dbl_nupdate);
    }
#endif
    dbl_EGlpNumFreeArray (f->work_coef);
    ILL_IFFREE (f->work_indx, int);
    ILL_IFFREE (f->uc_inf, dbl_uc_info);
    if (f->dim + f->max_k > 0 && f->ur_inf) {
	unsigned int i = f->dim + f->max_k + 1;
	while (i--)
	    dbl_EGlpNumClearVar (f->ur_inf[i].max);
    }
    ILL_IFFREE (f->ur_inf, dbl_ur_info);
    ILL_IFFREE (f->lc_inf, dbl_lc_info);
    ILL_IFFREE (f->lr_inf, dbl_lr_info);
    ILL_IFFREE (f->er_inf, dbl_er_info);
    ILL_IFFREE (f->ucindx, int);
    ILL_IFFREE (f->ucrind, int);
    dbl_EGlpNumFreeArray (f->uccoef);
    ILL_IFFREE (f->urindx, int);
    ILL_IFFREE (f->urcind, int);
    dbl_EGlpNumFreeArray (f->urcoef);
    ILL_IFFREE (f->lcindx, int);
    dbl_EGlpNumFreeArray (f->lccoef);
    ILL_IFFREE (f->lrindx, int);
    dbl_EGlpNumFreeArray (f->lrcoef);
    ILL_IFFREE (f->erindx, int);
    dbl_EGlpNumFreeArray (f->ercoef);
    ILL_IFFREE (f->rperm, int);
    ILL_IFFREE (f->rrank, int);
    ILL_IFFREE (f->cperm, int);
    ILL_IFFREE (f->crank, int);
    dbl_EGlpNumFreeArray (f->dmat);
    dbl_ILLsvector_free (&f->xtmp);
}

int dbl_ILLfactor_set_factor_iparam (dbl_factor_work * f,
      int param,
      int val)
{
    switch (param) {
    case QS_FACTOR_MAX_K:
	f->max_k = val;
	break;
    case QS_FACTOR_P:
	f->p = val;
	break;
    case QS_FACTOR_ETAMAX:
	f->etamax = val;
	break;
    case QS_FACTOR_DENSE_MIN:
	f->dense_min = val;
	break;
    default:
	fprintf (stderr, "Invalid param %d in dbl_ILLfactor_set_factor_iparam\n",
	    param);
	return 1;
    }
    return 0;
}

int dbl_ILLfactor_set_factor_dparam (dbl_factor_work * f,
      int param,
      double val)
{
    switch (param) {
    case QS_FACTOR_FZERO_TOL:
	dbl_EGlpNumCopy (f->fzero_tol, val);
	break;
    case QS_FACTOR_SZERO_TOL:
	dbl_EGlpNumCopy (f->szero_tol, val);
	break;
    case QS_FACTOR_UR_SPACE_MUL:
	f->ur_space_mul = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_UC_SPACE_MUL:
	f->uc_space_mul = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_LC_SPACE_MUL:
	f->lc_space_mul = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_LR_SPACE_MUL:
	f->lr_space_mul = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_ER_SPACE_MUL:
	f->er_space_mul = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_GROW_MUL:
	f->grow_mul = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_MAXMULT:
	f->maxmult = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_UPDMAXMULT:
	f->updmaxmult = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_DENSE_FRACT:
	f->dense_fract = dbl_EGlpNumToLf (val);
	break;
    case QS_FACTOR_PARTIAL_TOL:
	dbl_EGlpNumCopy (f->partial_tol, val);
	dbl_EGlpNumCopy (f->partial_cur, val);
	break;
    default:
	fprintf (stderr, "Invalid param %d in dbl_ILLfactor_set_factor_dparam\n",
	    param);
	return 1;
    }
    return 0;
}

int dbl_ILLfactor_create_factor_work (dbl_factor_work * f,
      int dim)
{
    int i;
    int rval;

    f->dim = dim;
    f->etacnt = 0;
    f->work_coef = dbl_EGlpNumAllocArray (dim);
    ILL_SAFE_MALLOC (f->work_indx, dim, int);
    ILL_SAFE_MALLOC (f->uc_inf, dim + (f->max_k + 1), dbl_uc_info);
    ILL_SAFE_MALLOC (f->ur_inf, dim + (f->max_k + 1), dbl_ur_info);
    ILL_SAFE_MALLOC (f->lc_inf, dim, dbl_lc_info);
    ILL_SAFE_MALLOC (f->lr_inf, dim, dbl_lr_info);
    ILL_SAFE_MALLOC (f->rperm, dim, int);
    ILL_SAFE_MALLOC (f->rrank, dim, int);
    ILL_SAFE_MALLOC (f->cperm, dim, int);
    ILL_SAFE_MALLOC (f->crank, dim, int);

    for (i = dim + f->max_k + 1; i--;)
	dbl_EGlpNumInitVar (f->ur_inf[i].max);

    for (i = 0; i < dim; i++) {
	dbl_EGlpNumZero (f->work_coef[i]);
	f->work_indx[i] = 0;
	f->uc_inf[i].nzcnt = 0;
	f->ur_inf[i].nzcnt = 0;
	f->lc_inf[i].nzcnt = 0;
	f->lr_inf[i].nzcnt = 0;
	f->rperm[i] = i;
	f->rrank[i] = i;
	f->cperm[i] = i;
	f->crank[i] = i;
    }
    for (i = 0; i <= f->max_k; i++) {
	f->uc_inf[dim + i].nzcnt = i;
	f->uc_inf[dim + i].next = dim + i;
	f->uc_inf[dim + i].prev = dim + i;
	f->ur_inf[dim + i].nzcnt = i;
	f->ur_inf[dim + i].next = dim + i;
	f->ur_inf[dim + i].prev = dim + i;
    }

    rval = dbl_ILLsvector_alloc (&f->xtmp, dim);
    ILL_CLEANUP_IF (rval);

    rval = 0;

CLEANUP:
    if (rval) {
	dbl_ILLfactor_free_factor_work (f);
    }
    ILL_RETURN (rval, "dbl_ILLfactor_create_factor_work");
}

#ifdef dbl_UPDATE_DEBUG
static void dbl_dump_matrix (dbl_factor_work * f,
      int remaining)
{
    int dim = f->dim;
    dbl_ur_info *ur_inf = f->ur_inf;
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_lc_info *lc_inf = f->lc_inf;
    dbl_lr_info *lr_inf = f->lr_inf;
    dbl_er_info *er_inf = f->er_inf;
    int nzcnt;
    int beg;

    int i;
    int j;

    for (i = 0; i < dim; i++) {
	if (!remaining || ur_inf[i].next >= 0) {
	    printf ("Row %d %d (max %.3f):", i, f->rrank[i],
		dbl_EGlpNumToLf (ur_inf[i].max));
	    nzcnt = ur_inf[i].nzcnt;
	    beg = ur_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		if (j == ur_inf[i].pivcnt) {
		    printf (" |");
		}
		printf (" %.3f*%d", dbl_EGlpNumToLf (f->urcoef[beg + j]),
		    f->urindx[beg + j]);
		if (f->urcind)
		    printf ("@%d", f->urcind[beg + j]);
	    }
	    printf ("\n");
	}
    }
    if (f->dmat) {
	int start = 0;
	if (remaining)
	    start = f->stage - f->dense_base;
	printf ("Dcols at %d %d - %d    :", f->stage - f->dense_base,
	    f->dense_base + start, f->nstages);
	for (j = start; j < f->dcols; j++) {
	    printf (" %5d", f->cperm[j + f->dense_base]);
	}
	printf ("\n");
	for (i = start; i < f->drows; i++) {
	    printf ("DRow %d %d (max %.3f):", i,
		f->rperm[i + f->dense_base],
		dbl_EGlpNumToLf (ur_inf[f->rperm[i + f->dense_base]].max));
	    for (j = start; j < f->dcols; j++) {
		if (j == f->drows) {
		    printf (" |");
		}
		printf (" %.3f", dbl_EGlpNumToLf (f->dmat[i * f->dcols + j]));
	    }
	    printf ("\n");
	}
    }
    if (!remaining) {
	for (i = 0; i < f->stage; i++) {
	    printf ("L col %d:", lc_inf[i].c);
	    nzcnt = lc_inf[i].nzcnt;
	    beg = lc_inf[i].cbeg;
	    for (j = 0; j < nzcnt; j++) {
		printf (" %.3f*%d", dbl_EGlpNumToLf (f->lccoef[beg + j]),
		    f->lcindx[beg + j]);
	    }
	    printf ("\n");
	}
	for (i = f->nstages; i < f->dim; i++) {
	    printf ("L col %d:", lc_inf[i].c);
	    nzcnt = lc_inf[i].nzcnt;
	    beg = lc_inf[i].cbeg;
	    for (j = 0; j < nzcnt; j++) {
		printf (" %.3f*%d", dbl_EGlpNumToLf (f->lccoef[beg + j]),
		    f->lcindx[beg + j]);
	    }
	    printf ("\n");
	}
	for (i = 0; i < f->dim; i++) {
	    if (!lr_inf[i].nzcnt)
		continue;
	    printf ("L row %d:", lr_inf[i].r);
	    nzcnt = lr_inf[i].nzcnt;
	    beg = lr_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		printf (" %.3f*%d", dbl_EGlpNumToLf (f->lrcoef[beg + j]),
		    f->lrindx[beg + j]);
	    }
	    printf ("\n");
	}
    }
    if (!remaining) {
	for (i = 0; i < f->etacnt; i++) {
	    printf ("Eta row %d:", f->er_inf[i].r);
	    nzcnt = er_inf[i].nzcnt;
	    beg = er_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		printf (" %.3f*%d", dbl_EGlpNumToLf (f->ercoef[beg + j]),
		    f->erindx[beg + j]);
	    }
	    printf ("\n");
	}
    }
    for (i = 0; i < dim; i++) {
	if (!remaining || uc_inf[i].next >= 0) {
	    printf ("Col %d %d:", i, f->crank[i]);
	    nzcnt = uc_inf[i].nzcnt;
	    beg = uc_inf[i].cbeg;
	    for (j = 0; j < nzcnt; j++) {
		if (f->uccoef != 0) {
		    printf (" %.3f*%d", dbl_EGlpNumToLf (f->uccoef[beg + j]),
			f->ucindx[beg + j]);
		    if (f->ucrind)
			printf ("@%d", f->ucrind[beg + j]);
		} else {
		    printf (" %d", f->ucindx[beg + j]);
		}
	    }
	    printf ("\n");
	}
    }

    if (!remaining) {
	printf ("rperm:");
	for (i = 0; i < dim; i++) {
	    if (i == f->nstages)
		printf ("|");
	    if (i == f->stage)
		printf ("|");
	    printf (" %d", f->rperm[i]);
	}
	printf ("\n");

	printf ("cperm:");
	for (i = 0; i < dim; i++) {
	    if (i == f->nstages)
		printf ("|");
	    if (i == f->stage)
		printf ("|");
	    printf (" %d", f->cperm[i]);
	}
	printf ("\n");
    }
    printf ("Rows by nzcnt:\n");
    for (i = 0; i <= f->max_k; i++) {
	if (ur_inf[dim + i].next != dim + i) {
	    printf ("%d:", i);
	    for (j = ur_inf[dim + i].next; j != dim + i; j = ur_inf[j].next) {
		printf (" %d", j);
	    }
	    printf ("\n");
	}
    }

    printf ("Cols by nzcnt:\n");
    for (i = 0; i <= f->max_k; i++) {
	if (uc_inf[dim + i].next != dim + i) {
	    printf ("%d:", i);
	    for (j = uc_inf[dim + i].next; j != dim + i; j = uc_inf[j].next) {
		printf (" %d", j);
	    }
	    printf ("\n");
	}
    }

    printf ("\n");
    fflush (stdout);
}
#endif

#ifdef dbl_SORT_RESULTS
static void dbl_sort_vector2 (int nzcnt,
      int *indx,
      double *coef)
{
    int i;
    int j;
    int itmp;
    double ctmp;
    dbl_EGlpNumInitVar (ctmp);

    for (i = 1; i < nzcnt; i++) {
	itmp = indx[i];
	dbl_EGlpNumCopy (ctmp, coef[i]);
	for (j = i; j >= 1 && indx[j - 1] > itmp; j--) {
	    indx[j] = indx[j - 1];
	    dbl_EGlpNumCopy (coef[j], coef[j - 1]);
	}
	indx[j] = itmp;
	dbl_EGlpNumCopy (coef[j], ctmp);
    }
    dbl_EGlpNumClearVar (ctmp);
}

static void dbl_sort_vector (dbl_svector * x)
{
    dbl_sort_vector2 (x->nzcnt, x->indx, x->coef);
}
#endif

#ifdef dbl_DEBUG_FACTOR
static int dbl_check_matrix (dbl_factor_work * f)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    dbl_uc_info *uc_inf = f->uc_inf;
    int rbeg;
    int nzcnt;
    int cbeg;
    int c;
    int r;
    int j;
    int nerr = 0;

    for (r = 0; r < f->dim; r++) {
	nzcnt = ur_inf[r].nzcnt;
	rbeg = ur_inf[r].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    c = f->urindx[rbeg + j];
	    cbeg = uc_inf[c].cbeg;
	    if (f->ucindx[cbeg + f->urcind[rbeg + j]] != r) {
		printf ("index mismatch, row %d column %d\n", r, c);
		nerr++;
	    }
	    if (f->uccoef[cbeg + f->urcind[rbeg + j]] != f->urcoef[rbeg + j]) {
		printf ("coef mismatch, row %d column %d\n", r, c);
		nerr++;
	    }
	}
    }
    if (f->urindx[f->ur_space] != 0) {
	printf ("last urindx entry %d != 0\n", f->urindx[f->ur_space]);
	nerr++;
    }
    for (c = 0; c < f->dim; c++) {
	nzcnt = uc_inf[c].nzcnt;
	cbeg = uc_inf[c].cbeg;
	for (j = 0; j < nzcnt; j++) {
	    r = f->ucindx[cbeg + j];
	    rbeg = ur_inf[r].rbeg;
	    if (f->urindx[rbeg + f->ucrind[cbeg + j]] != c) {
		printf ("index mismatch, column %d row %d\n", c, r);
		nerr++;
	    }
	    if (f->urcoef[rbeg + f->ucrind[cbeg + j]] != f->uccoef[cbeg + j]) {
		printf ("coef mismatch, column %d row %d\n", c, r);
		nerr++;
	    }
	}
    }
    if (f->ucindx[f->uc_space] != 0) {
	printf ("last ucindx entry %d != 0\n", f->ucindx[f->uc_space]);
    }
    if (nerr) {
	fflush (stdout);
	dbl_dump_matrix (f, 0);
	return E_CHECK_FAILED;
    }
    return 0;
}
#endif

#ifdef dbl_FACTOR_STATS
static void dbl_dump_factor_stats (dbl_factor_work * f)
{
    int dim = f->dim;
    int ecnt = f->etacnt;
    dbl_ur_info *ur_inf = f->ur_inf;
    dbl_lc_info *lc_inf = f->lc_inf;
    dbl_er_info *er_inf = f->er_inf;
    double *urcoef = f->urcoef;
    double *lccoef = f->lccoef;
    double *ercoef = f->ercoef;
    int lnzcnt = 0;
    int unzcnt = 0;
    int enzcnt = 0;
    int nzcnt;
    int beg;
    double umax;
    double lmax;
    double emax;
    int i;
    int j;
    dbl_EGlpNumInitVar (umax);
    dbl_EGlpNumInitVar (lmax);
    dbl_EGlpNumInitVar (emax);
    dbl_EGlpNumZero (umax);
    for (i = 0; i < dim; i++) {
	nzcnt = ur_inf[i].nzcnt;
	beg = ur_inf[i].rbeg;
	unzcnt += nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumSetToMaxAbs (umax, urcoef[beg + j]);
	}
    }
    dbl_EGlpNumZero (lmax);
    for (i = 0; i < dim; i++) {
	nzcnt = lc_inf[i].nzcnt;
	beg = lc_inf[i].cbeg;
	lnzcnt += nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumSetToMaxAbs (lmax, lccoef[beg + j]);
	}
    }
    dbl_EGlpNumZero (emax);
    for (i = 0; i < ecnt; i++) {
	nzcnt = er_inf[i].nzcnt;
	beg = er_inf[i].rbeg;
	enzcnt += nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumSetToMaxAbs (emax, ercoef[beg + j]);
	}
    }
    printf
	("factor U %d nzs %.3e max L %d nzs %.3e max E %d nzs %.3e max\n",
	unzcnt, dbl_EGlpNumToLf (umax), lnzcnt, dbl_EGlpNumToLf (lmax), enzcnt,
	dbl_EGlpNumToLf (emax));
    fflush (stdout);
    dbl_EGlpNumClearVar (umax);
    dbl_EGlpNumClearVar (lmax);
    dbl_EGlpNumClearVar (emax);
}
#endif

static void dbl_clear_work (dbl_factor_work * f)
{
    int i;
    int dim = f->dim;
    double *work_coef = f->work_coef;

    for (i = 0; i < dim; i++) {
	dbl_EGlpNumZero (work_coef[i]);
    }
}

static void dbl_load_row (dbl_factor_work * f,
      int r)
{
    double *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    double *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;

    for (i = 0; i < prow_nzcnt; i++) {
	j = prow_urindx[i];
	dbl_EGlpNumCopy (work_coef[j], prow_urcoef[i]);
	work_indx[j] = 1;
    }
}

static void dbl_clear_row (dbl_factor_work * f,
      int r)
{
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    double *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;

    for (i = 0; i < prow_nzcnt; i++) {
	j = prow_urindx[i];
	dbl_EGlpNumZero (work_coef[j]);
	work_indx[j] = 0;
    }
}

static int dbl_make_ur_space (dbl_factor_work * f,
      int space)
{
    double *new_urcoef = 0;
    int *new_urindx = 0;
    int *new_urcind = 0;
    double *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int minspace;
    dbl_ur_info *ur_inf = f->ur_inf;
    int dim = f->dim;
    int new_nzcnt = 0, old_nzcnt;
    int rbeg;
    int nzcnt;
    int i;
    int j;
    int rval;

    minspace = f->ur_space;
    nzcnt = space;
    for (i = 0; i < dim; i++)
	nzcnt += ur_inf[i].nzcnt;
    old_nzcnt = nzcnt;
    while (nzcnt * 2 >= minspace) {
	minspace = 1 + minspace * f->grow_mul;
    }

#ifdef dbl_GROWTH_STATS
    printf ("dbl_make_ur_space growing from %d to %d...", f->ur_space, minspace);
    fflush (stdout);
#endif
    new_urcoef = dbl_EGlpNumAllocArray (minspace);
    ILL_SAFE_MALLOC (new_urindx, minspace + 1, int);

    if (urcind) {
	ILL_SAFE_MALLOC (new_urcind, minspace, int);
    }
    if (urcind) {
	for (j = 0; j < dim; j++) {
	    rbeg = ur_inf[j].rbeg;
	    nzcnt = ur_inf[j].nzcnt;
	    ur_inf[j].rbeg = new_nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		new_urindx[new_nzcnt] = urindx[rbeg + i];
		dbl_EGlpNumCopy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
		new_urcind[new_nzcnt] = urcind[rbeg + i];
		new_nzcnt++;
	    }
	}
    } else {
	for (j = 0; j < dim; j++) {
	    rbeg = ur_inf[j].rbeg;
	    nzcnt = ur_inf[j].nzcnt;
	    ur_inf[j].rbeg = new_nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		new_urindx[new_nzcnt] = urindx[rbeg + i];
		dbl_EGlpNumCopy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
		new_nzcnt++;
	    }
	}
    }

    for (i = new_nzcnt; i < minspace; i++) {
	new_urindx[i] = -1;
    }
    new_urindx[minspace] = 0;
    dbl_EGlpNumFreeArray (f->urcoef);
    f->urcoef = new_urcoef;
    new_urcoef = 0;

    ILL_IFFREE (f->urindx, int);
    f->urindx = new_urindx;
    new_urindx = 0;

    ILL_IFFREE (f->urcind, int);
    f->urcind = new_urcind;
    new_urcind = 0;

    f->ur_freebeg = new_nzcnt;
    f->ur_space = minspace;

#ifdef dbl_GROWTH_STATS
    printf ("%d/%d nonzeros\n", new_nzcnt, old_nzcnt);
    fflush (stdout);
    dbl_dump_factor_stats (f);
#endif

    rval = 0;

CLEANUP:
    ILL_IFFREE (new_urcoef, double);
    ILL_IFFREE (new_urindx, int);
    ILL_IFFREE (new_urcind, int);
    ILL_RETURN (rval, "dbl_make_ur_space");
}

static int dbl_make_uc_space (dbl_factor_work * f,
      int space)
{
    double *new_uccoef = 0;
    int *new_ucindx = 0;
    int *new_ucrind = 0;
    int uc_freebeg = f->uc_freebeg;
    double *uccoef = f->uccoef;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int minspace = uc_freebeg + space;
    dbl_uc_info *uc_inf = f->uc_inf;
    int dim = f->dim;
    int new_nzcnt = 0;
    int cbeg;
    int nzcnt;
    int i;
    int j;
    int rval;
#if 0
    if (f->uc_space * f->grow_mul > minspace) {
        minspace = f->uc_space * f->grow_mul;
    }
#endif
    /* DAVE New -- 050420 */
    minspace = f->uc_space;
    nzcnt = space;
    for (i=0; i<dim; i++) nzcnt += uc_inf[i].nzcnt;
    while (nzcnt * 2 >= minspace) {
        minspace *= f->grow_mul;
    }
    /* End DAVE New -- 050420 */

#ifdef dbl_GROWTH_STATS
    printf ("dbl_make_uc_space growing from %d to %d...", f->uc_space, minspace);
    fflush (stdout);
#endif

    ILL_SAFE_MALLOC (new_ucindx, minspace + 1, int);

    if (ucrind) {
	new_uccoef = dbl_EGlpNumAllocArray (minspace);
	ILL_SAFE_MALLOC (new_ucrind, minspace, int);
    }
    if (ucrind) {
	for (j = 0; j < dim; j++) {
	    cbeg = uc_inf[j].cbeg;
	    nzcnt = uc_inf[j].nzcnt;
	    uc_inf[j].cbeg = new_nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		new_ucindx[new_nzcnt] = ucindx[cbeg + i];
		dbl_EGlpNumCopy (new_uccoef[new_nzcnt], uccoef[cbeg + i]);
		new_ucrind[new_nzcnt] = ucrind[cbeg + i];
		new_nzcnt++;
	    }
	}
    } else {
	for (j = 0; j < dim; j++) {
	    cbeg = uc_inf[j].cbeg;
	    nzcnt = uc_inf[j].nzcnt;
	    uc_inf[j].cbeg = new_nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		new_ucindx[new_nzcnt] = ucindx[cbeg + i];
		new_nzcnt++;
	    }
	}
    }

    for (i = new_nzcnt; i < minspace; i++) {
	new_ucindx[i] = -1;
    }
    new_ucindx[minspace] = 0;

    dbl_EGlpNumFreeArray (f->uccoef);
    f->uccoef = new_uccoef;
    new_uccoef = 0;

    ILL_IFFREE (f->ucindx, int);
    f->ucindx = new_ucindx;
    new_ucindx = 0;

    ILL_IFFREE (f->ucrind, int);
    f->ucrind = new_ucrind;
    new_ucrind = 0;

    f->uc_freebeg = new_nzcnt;
    f->uc_space = minspace;

#ifdef dbl_GROWTH_STATS
    printf ("%d nonzeros\n", new_nzcnt);
    fflush (stdout);
    dbl_dump_factor_stats (f);
#endif

    rval = 0;

CLEANUP:
    ILL_IFFREE (new_uccoef, double);
    ILL_IFFREE (new_ucindx, int);
    ILL_IFFREE (new_ucrind, int);
    ILL_RETURN (rval, "dbl_make_uc_space");
}

static int dbl_make_lc_space (dbl_factor_work * f,
      int space)
{
    double *new_lccoef = 0;
    int *new_lcindx = 0;
    int lc_freebeg = f->lc_freebeg;
    double *lccoef = f->lccoef;
    int *lcindx = f->lcindx;
    int minspace = lc_freebeg + space;
    int i;
    int rval;
    if (f->lc_space * f->grow_mul > minspace) {
	minspace = f->lc_space * f->grow_mul;
    }
#ifdef dbl_GROWTH_STATS
    printf ("dbl_make_lc_space growing from %d to %d...", f->lc_space, minspace);
    fflush (stdout);
#endif

    new_lccoef = dbl_EGlpNumAllocArray (minspace);
    ILL_SAFE_MALLOC (new_lcindx, minspace, int);

    for (i = 0; i < lc_freebeg; i++) {
	dbl_EGlpNumCopy (new_lccoef[i], lccoef[i]);
	new_lcindx[i] = lcindx[i];
    }

    dbl_EGlpNumFreeArray (lccoef);
    f->lccoef = new_lccoef;
    new_lccoef = 0;

    ILL_IFFREE (lcindx, int);
    f->lcindx = new_lcindx;
    new_lcindx = 0;

    f->lc_space = minspace;

#ifdef dbl_GROWTH_STATS
    printf ("done\n");
    fflush (stdout);
    dbl_dump_factor_stats (f);
#endif

    rval = 0;

CLEANUP:
    ILL_IFFREE (new_lccoef, double);
    ILL_IFFREE (new_lcindx, int);
    ILL_RETURN (rval, "dbl_make_lc_space");
}

static void dbl_set_col_nz (dbl_factor_work * f,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int nzcnt = uc_inf[c].nzcnt;
    int max_k = f->max_k;
    int dim = f->dim;

    if (uc_inf[c].next >= 0) {
	uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
	uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

	if (nzcnt >= max_k)
	    nzcnt = max_k;
	uc_inf[c].next = uc_inf[dim + nzcnt].next;
	uc_inf[c].prev = dim + nzcnt;
	uc_inf[dim + nzcnt].next = c;
	uc_inf[uc_inf[c].next].prev = c;
    }
}

static void dbl_set_row_nz (dbl_factor_work * f,
      int r)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int nzcnt = ur_inf[r].pivcnt;
    int max_k = f->max_k;
    int dim = f->dim;

    if (ur_inf[r].next >= 0) {
	ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
	ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

	if (nzcnt >= max_k)
	    nzcnt = max_k;
	ur_inf[r].next = ur_inf[dim + nzcnt].next;
	ur_inf[r].prev = dim + nzcnt;
	ur_inf[dim + nzcnt].next = r;
	ur_inf[ur_inf[r].next].prev = r;
    }
}

static void dbl_remove_col_nz (dbl_factor_work * f,
      int r,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int *ucindx = f->ucindx + uc_inf[c].cbeg;
    int nzcnt = uc_inf[c].nzcnt;
    int i;

    for (i = 0; i < nzcnt; i++) {
	if (ucindx[i] == r) {
	    --nzcnt;
	    ucindx[i] = ucindx[nzcnt];
	    ucindx[nzcnt] = -1;
	    break;
	}
    }
    uc_inf[c].nzcnt = nzcnt;

    dbl_set_col_nz (f, c);
}

static void dbl_remove_row_nz (dbl_factor_work * f,
      int r,
      int c)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int *urindx = f->urindx + ur_inf[r].rbeg;
    double *urcoef = f->urcoef + ur_inf[r].rbeg;
    int pivcnt = ur_inf[r].pivcnt;
    double max;
    int tind;
    double tcoef;
    int i;
    dbl_EGlpNumInitVar (tcoef);
    dbl_EGlpNumInitVar (max);
    dbl_EGlpNumZero (max);

    for (i = 0; i < pivcnt; i++) {
	if (urindx[i] == c) {
	    --pivcnt;
	    dbl_ILL_SWAP (urindx[i], urindx[pivcnt], tind);
	    dbl_EGLPNUM_SWAP (urcoef[i], urcoef[pivcnt], tcoef);
	    --i;
	} else {
	    dbl_EGlpNumSetToMaxAbs (max, urcoef[i]);
	}
    }
    ur_inf[r].pivcnt = pivcnt;
    dbl_EGlpNumCopy (ur_inf[r].max, max);
    dbl_set_row_nz (f, r);
    dbl_EGlpNumClearVar (max);
    dbl_EGlpNumClearVar (tcoef);
}

static int dbl_add_col_nz (dbl_factor_work * f,
      int r,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int cbeg = uc_inf[c].cbeg;
    int nzcnt = uc_inf[c].nzcnt;
    int uc_freebeg = f->uc_freebeg;
    int *ucindx = f->ucindx;
    int i;
    int rval = 0;

    if (uc_inf[c].next == -1) {
	return 0;
    }
    if (ucindx[cbeg + nzcnt] == -1) {
	ucindx[cbeg + nzcnt] = r;
	uc_inf[c].nzcnt++;
	if (nzcnt + cbeg == uc_freebeg) {
	    f->uc_freebeg = uc_freebeg + 1;
	}
    } else {
	if (uc_freebeg + nzcnt + 1 >= f->uc_space) {
	    rval = dbl_make_uc_space (f, nzcnt + 1);
	    ILL_CLEANUP_IF (rval);
	    uc_freebeg = f->uc_freebeg;
	    cbeg = uc_inf[c].cbeg;
	    ucindx = f->ucindx;
	}
	for (i = 0; i < nzcnt; i++) {
	    ucindx[uc_freebeg + i] = ucindx[cbeg + i];
	    ucindx[cbeg + i] = -1;
	}
	ucindx[uc_freebeg + nzcnt] = r;
	uc_inf[c].cbeg = uc_freebeg;
	uc_inf[c].nzcnt++;
	f->uc_freebeg = uc_freebeg + nzcnt + 1;
    }

    dbl_set_col_nz (f, c);
CLEANUP:
    ILL_RETURN (rval, "dbl_add_col_nz");
}

static void dbl_disable_col (dbl_factor_work * f,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;

    if (uc_inf[c].next >= 0) {
	uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
	uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

	uc_inf[c].next = -2;
	uc_inf[c].prev = -2;
    }
}

static void dbl_remove_col (dbl_factor_work * f,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int cbeg = uc_inf[c].cbeg;
    int nzcnt = uc_inf[c].nzcnt;
    int *ucindx = f->ucindx;
    int i;

    for (i = 0; i < nzcnt; i++) {
	ucindx[cbeg + i] = -1;
    }
    uc_inf[c].cbeg = 0;
    uc_inf[c].nzcnt = 0;

    if (uc_inf[c].next >= 0) {
	uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
	uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

	uc_inf[c].next = -1;
	uc_inf[c].prev = -1;
    }
}

static void dbl_remove_row (dbl_factor_work * f,
      int r)
{
    dbl_ur_info *ur_inf = f->ur_inf;

    if (ur_inf[r].next >= 0) {
	ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
	ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

	ur_inf[r].next = -1;
	ur_inf[r].prev = -1;
    }
}

static void dbl_find_coef (dbl_factor_work * f,
      int r,
      int c,
      double *coef)
{
    double *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int i;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    dbl_EGlpNumZero (*coef);
    for (i = 0; i < prow_nzcnt; i++) {
	if (prow_urindx[i] == c) {
	    dbl_EGlpNumCopy (*coef, prow_urcoef[i]);
	    return;
	}
    }
    fprintf (stderr, "Coefficient not found\n");
    return;
}

static int dbl_elim_row (dbl_factor_work * f,
      int elim_r,
      int r,
      int c,
      double *p_pivot_coef)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    double *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    double *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int prow_beg = ur_inf[r].rbeg;
    int prow_nzcnt = ur_inf[r].nzcnt;
    int prow_pivcnt = ur_inf[r].pivcnt;
    int fill = ur_inf[elim_r].nzcnt;
    int cancel = 0;
    double max;
    int erow_beg;
    int erow_nzcnt;
    int erow_pivcnt;
    double x;
    int i;
    int j;
    int rval = 0;
    double elim_coef;
    dbl_EGlpNumInitVar (max);
    dbl_EGlpNumInitVar (x);
    dbl_EGlpNumInitVar (elim_coef);
    dbl_EGlpNumZero (max);
    dbl_find_coef (f, r, c, &elim_coef);
    dbl_EGlpNumDivTo (elim_coef, work_coef[c]);
    dbl_EGlpNumCopy (*p_pivot_coef, elim_coef);

    for (i = 0; i < prow_nzcnt; i++) {
	j = urindx[prow_beg + i];
	if (work_indx[j] == 1) {
	    dbl_EGlpNumCopy (x, urcoef[prow_beg + i]);
	    dbl_EGlpNumSubInnProdTo (x, elim_coef, work_coef[j]);
	    if (dbl_EGlpNumIsEqual (x, dbl_zeroLpNum, f->fzero_tol) || j == c) {
		cancel++;
		if (j != c) {
		    dbl_remove_col_nz (f, r, j);
		}
		if (i < prow_pivcnt) {
		    prow_pivcnt--;
		    prow_nzcnt--;
		    urindx[prow_beg + i] = urindx[prow_beg + prow_pivcnt];
		    dbl_EGlpNumCopy (urcoef[prow_beg + i], urcoef[prow_beg + prow_pivcnt]);
		    if (prow_pivcnt != prow_nzcnt) {
			urindx[prow_beg + prow_pivcnt] = urindx[prow_beg + prow_nzcnt];
			dbl_EGlpNumCopy (urcoef[prow_beg + prow_pivcnt],
			    urcoef[prow_beg + prow_nzcnt]);
		    }
		} else {
		    prow_nzcnt--;
		    urindx[prow_beg + i] = urindx[prow_beg + prow_nzcnt];
		    dbl_EGlpNumCopy (urcoef[prow_beg + i], urcoef[prow_beg + prow_nzcnt]);
		}
		urindx[prow_beg + prow_nzcnt] = -1;
		i--;
	    } else {
		dbl_EGlpNumCopy (urcoef[prow_beg + i], x);
		if (i < prow_pivcnt) {
		    dbl_EGlpNumSetToMaxAbs (max, x);
		}
	    }
	    work_indx[j] = 0;
	    fill--;
	} else {
	    if (i < prow_pivcnt) {
		dbl_EGlpNumSetToMaxAbs (max, urcoef[prow_beg + i]);
	    }
	}
    }

    if (fill > 0) {
	ur_inf[r].nzcnt = prow_nzcnt;
	ur_inf[r].pivcnt = prow_pivcnt;
	if (fill > cancel) {
	    int ur_freebeg = f->ur_freebeg;

	    if (ur_freebeg + prow_nzcnt + fill >= f->ur_space) {
		rval = dbl_make_ur_space (f, prow_nzcnt + fill);
		ILL_CLEANUP_IF (rval);
		urcoef = f->urcoef;
		urindx = f->urindx;
		ur_freebeg = f->ur_freebeg;
		prow_beg = f->ur_inf[r].rbeg;
	    }
	    for (i = 0; i < prow_nzcnt; i++) {
		urindx[ur_freebeg + i] = urindx[prow_beg + i];
		dbl_EGlpNumCopy (urcoef[ur_freebeg + i], urcoef[prow_beg + i]);
		urindx[prow_beg + i] = -1;
	    }
	    ur_inf[r].rbeg = ur_freebeg;
	    f->ur_freebeg = ur_freebeg + prow_nzcnt + fill;
	    prow_beg = ur_freebeg;
	}
	erow_beg = ur_inf[elim_r].rbeg;
	erow_nzcnt = ur_inf[elim_r].nzcnt;
	erow_pivcnt = ur_inf[elim_r].pivcnt;

	for (i = 0; i < erow_pivcnt; i++) {
	    j = urindx[erow_beg + i];
	    if (work_indx[j] == 1) {
		dbl_EGlpNumCopyNeg (x, elim_coef);
		dbl_EGlpNumMultTo (x, urcoef[erow_beg + i]);
		if (dbl_EGlpNumIsNeqZero (x, f->fzero_tol)) {
		    rval = dbl_add_col_nz (f, r, j);
		    ILL_CLEANUP_IF (rval);
		    if (prow_pivcnt != prow_nzcnt) {
			urindx[prow_beg + prow_nzcnt] = urindx[prow_beg + prow_pivcnt];
			dbl_EGlpNumCopy (urcoef[prow_beg + prow_nzcnt],
			    urcoef[prow_beg + prow_pivcnt]);
		    }
		    urindx[prow_beg + prow_pivcnt] = j;
		    dbl_EGlpNumCopy (urcoef[prow_beg + prow_pivcnt], x);
		    dbl_EGlpNumSetToMaxAbs (max, x);
		    prow_pivcnt++;
		    prow_nzcnt++;
		}
	    } else {
		work_indx[j] = 1;
	    }
	}
	for (i = erow_pivcnt; i < erow_nzcnt; i++) {
	    j = urindx[erow_beg + i];
	    if (work_indx[j] == 1) {
		dbl_EGlpNumCopyNeg (x, elim_coef);
		dbl_EGlpNumMultTo (x, urcoef[erow_beg + i]);
		if (dbl_EGlpNumIsNeqZero (x, f->fzero_tol)) {
		    rval = dbl_add_col_nz (f, r, j);
		    ILL_CLEANUP_IF (rval);
		    urindx[prow_beg + prow_nzcnt] = j;
		    dbl_EGlpNumCopy (urcoef[prow_beg + prow_nzcnt], x);
		    prow_nzcnt++;
		}
	    } else {
		work_indx[j] = 1;
	    }
	}
    } else {
	erow_nzcnt = ur_inf[elim_r].nzcnt;
	erow_beg = ur_inf[elim_r].rbeg;
	for (i = 0; i < erow_nzcnt; i++) {
	    j = urindx[erow_beg + i];
	    work_indx[j] = 1;
	}
    }

    ur_inf[r].nzcnt = prow_nzcnt;
    ur_inf[r].pivcnt = prow_pivcnt;
    dbl_EGlpNumCopy (ur_inf[r].max, max);

    dbl_set_row_nz (f, r);
CLEANUP:
    dbl_EGlpNumClearVar (elim_coef);
    dbl_EGlpNumClearVar (x);
    dbl_EGlpNumClearVar (max);
    ILL_RETURN (rval, "dbl_elim_row");
}

#define dbl_SETPERM(f,s,r,c) {                    \
        f->rperm[f->rrank[r]] = f->rperm[s];  \
        f->rrank[f->rperm[s]] = f->rrank[r];  \
        f->rperm[s] = r;                      \
        f->rrank[r] = s;                      \
                                              \
        f->cperm[f->crank[c]] = f->cperm[s];  \
        f->crank[f->cperm[s]] = f->crank[c];  \
        f->cperm[s] = c;                      \
        f->crank[c] = s;                      \
}

static int dbl_elim (dbl_factor_work * f,
      int r,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_ur_info *ur_inf = f->ur_inf;
    dbl_lc_info *lc_inf = f->lc_inf;
    int *urindx;
    int *ucindx;
    int *lcindx;
    double *urcoef;
    double *lccoef;
    double pivot_coef;
    int nzcnt;
    int lc_freebeg;
    int s = f->stage;
    int i;
    int j;
    int rval = 0;
    dbl_EGlpNumInitVar (pivot_coef);

    if (uc_inf[c].nzcnt == 1) {
	/* col singleton */
	dbl_SETPERM (f, s, r, c);

	lc_inf[s].cbeg = -1;
	lc_inf[s].c = r;
	lc_inf[s].nzcnt = 0;
	f->stage++;

	urindx = f->urindx + ur_inf[r].rbeg;
	urcoef = f->urcoef + ur_inf[r].rbeg;
	nzcnt = ur_inf[r].nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    j = urindx[i];
	    dbl_remove_col_nz (f, r, j);
	    if (j == c) {
		urindx[i] = urindx[0];
		urindx[0] = c;
		dbl_EGLPNUM_SWAP (urcoef[0], urcoef[i], pivot_coef);
	    }
	}
	dbl_remove_row (f, r);
	dbl_remove_col (f, c);
    } else if (ur_inf[r].nzcnt == 1) {
	/* row singleton */
	--(f->nstages);
	dbl_SETPERM (f, f->nstages, r, c);

	lc_inf[f->nstages].cbeg = -1;
	lc_inf[f->nstages].c = r;
	lc_inf[f->nstages].nzcnt = 0;

	ucindx = f->ucindx + uc_inf[c].cbeg;
	nzcnt = uc_inf[c].nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    j = ucindx[i];
	    dbl_remove_row_nz (f, j, c);
	}
	dbl_remove_row (f, r);
	dbl_remove_col (f, c);
    } else {
	dbl_SETPERM (f, s, r, c);
	f->stage++;

	nzcnt = uc_inf[c].nzcnt;
	if (f->lc_freebeg + nzcnt >= f->lc_space) {
	    rval = dbl_make_lc_space (f, nzcnt);
	    ILL_CLEANUP_IF (rval);
	}
	lc_freebeg = f->lc_freebeg;
	lc_inf[s].cbeg = lc_freebeg;
	lc_inf[s].c = r;
	lcindx = f->lcindx;
	lccoef = f->lccoef;
	dbl_load_row (f, r);
	ucindx = f->ucindx + uc_inf[c].cbeg;
	for (i = 0; i < nzcnt; i++) {
	    j = f->ucindx[uc_inf[c].cbeg + i];
	    if (j != r) {
		rval = dbl_elim_row (f, r, j, c, &pivot_coef);
		ILL_CLEANUP_IF (rval);
		lcindx[lc_freebeg] = j;
		dbl_EGlpNumCopy (lccoef[lc_freebeg], pivot_coef);
		lc_freebeg++;
#ifdef dbl_TRACK_FACTOR
		dbl_EGlpNumSetToMaxAbs (f->maxelem_factor, pivot_coef);
		if (dbl_EGlpNumIsLess (f->maxelem_factor, ur_inf[r].max))
		    dbl_EGlpNumCopy (f->maxelem_factor, ur_inf[r].max);
#endif				/* dbl_TRACK_FACTOR */
	    }
	}
	lc_inf[s].nzcnt = lc_freebeg - lc_inf[s].cbeg;
	f->lc_freebeg = lc_freebeg;

	dbl_clear_row (f, r);

	urindx = f->urindx + ur_inf[r].rbeg;
	urcoef = f->urcoef + ur_inf[r].rbeg;
	nzcnt = ur_inf[r].nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    j = urindx[i];
	    dbl_remove_col_nz (f, r, j);
	    if (j == c) {
		urindx[i] = urindx[0];
		urindx[0] = c;
		dbl_EGLPNUM_SWAP (urcoef[0], urcoef[i], pivot_coef);
	    }
	}
	dbl_remove_row (f, r);
	dbl_remove_col (f, c);
    }
CLEANUP:
    dbl_EGlpNumClearVar (pivot_coef);
    ILL_RETURN (rval, "dbl_elim");
}

static void dbl_find_pivot_column (dbl_factor_work * f,
      int c,
      int *p_r)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *ucindx = f->ucindx;
    int nzcnt = uc_inf[c].nzcnt;
    int cbeg = uc_inf[c].cbeg;
    double num_tmp[2];
    int bestnz = -1;
    int i;
    int r;
    dbl_EGlpNumInitVar (num_tmp[0]);
    dbl_EGlpNumInitVar (num_tmp[1]);

    *p_r = -1;
    for (i = 0; i < nzcnt; i++) {
	r = ucindx[cbeg + i];
        if ((bestnz == -1 || ur_inf[r].pivcnt < bestnz)) {
            dbl_find_coef (f, r, c, &num_tmp[0]);
            if (dbl_EGlpNumIsLess (num_tmp[0], dbl_zeroLpNum))
                dbl_EGlpNumSign (num_tmp[0]);
            dbl_EGlpNumCopy (num_tmp[1], f->partial_cur);
            dbl_EGlpNumMultTo (num_tmp[1], ur_inf[r].max);
            if (dbl_EGlpNumIsLeq (num_tmp[1], num_tmp[0])) {
                bestnz = ur_inf[r].pivcnt;
                *p_r = r;
            }
	}
    }
    dbl_EGlpNumClearVar (num_tmp[0]);
    dbl_EGlpNumClearVar (num_tmp[1]);
}

static void dbl_find_pivot_row (dbl_factor_work * f,
      int r,
      int *p_c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *urindx = f->urindx;
    double *urcoef = f->urcoef;
    int pivcnt = ur_inf[r].pivcnt;
    int rbeg = ur_inf[r].rbeg;
    double thresh[2];
    int bestnz = -1;
    int i;
    int c;
    dbl_EGlpNumInitVar (thresh[0]);
    dbl_EGlpNumInitVar (thresh[1]);
    dbl_EGlpNumCopy (thresh[0], f->partial_cur);
    dbl_EGlpNumMultTo (thresh[0], ur_inf[r].max);
    *p_c = -1;
    for (i = 0; i < pivcnt; i++) {
        c = urindx[rbeg + i];
        if ((bestnz == -1 || uc_inf[c].nzcnt < bestnz)) {
            /* DAVID FIX */
            dbl_EGlpNumCopyAbs (thresh[1], urcoef[rbeg + i]);
            if (dbl_EGlpNumIsLeq (thresh[0], thresh[1])) {
                bestnz = uc_inf[c].nzcnt;
                *p_c = c;
            }
        }
    }
    dbl_EGlpNumClearVar (thresh[0]);
    dbl_EGlpNumClearVar (thresh[1]);
}

static int dbl_find_pivot (dbl_factor_work * f,
      int *p_r,
      int *p_c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_ur_info *ur_inf = f->ur_inf;
    int dim = f->dim;
    int max_k = f->max_k;
    int p = f->p;
    int c;
    int r;
    int mm = 0;
    int n = 0;
    int m;
    int k = 2;

    if (uc_inf[dim + 1].next != dim + 1) {
	c = uc_inf[dim + 1].next;
	r = f->ucindx[uc_inf[c].cbeg];
	*p_c = c;
	*p_r = r;
	return 0;
    } else if (ur_inf[dim + 1].next != dim + 1) {
	r = ur_inf[dim + 1].next;
	c = f->urindx[ur_inf[r].rbeg];
	*p_c = c;
	*p_r = r;
	return 0;
    }
    *p_r = -1;
    *p_c = -1;
    for (; k <= max_k && (mm == 0 || mm > (k - 1) * (k - 1)); k++) {
	if (uc_inf[dim + k].next != dim + k) {
	    for (c = uc_inf[dim + k].next; c != dim + k; c = uc_inf[c].next) {
		dbl_find_pivot_column (f, c, &r);
		if (r >= 0) {
		    m = (uc_inf[c].nzcnt - 1) * (ur_inf[r].pivcnt - 1);
		    if (mm == 0 || m < mm) {
			mm = m;
			*p_c = c;
			*p_r = r;
			if (mm <= (k - 1) * (k - 1)) {
			    return 0;
			}
		    }
		} else {
		    c = uc_inf[c].prev;
		    dbl_disable_col (f, uc_inf[c].next);
		}
		n++;
		if (n >= p && mm != 0) {
		    return 0;
		}
	    }
	}
	if (ur_inf[dim + k].next != dim + k) {
	    for (r = ur_inf[dim + k].next; r != dim + k; r = ur_inf[r].next) {
		dbl_find_pivot_row (f, r, &c);
		if (c >= 0) {
		    m = (uc_inf[c].nzcnt - 1) * (ur_inf[r].pivcnt - 1);
		    if (mm == 0 || m < mm) {
			mm = m;
			*p_c = c;
			*p_r = r;
			if (mm <= k * (k - 1)) {
			    return 0;
			}
		    }
		}
		n++;
		if (n >= p && mm != 0) {
		    return 0;
		}
	    }
	}
    }
    if (mm != 0) {
	return 0;
    } else {
	/* fprintf (stderr, "No acceptable pivot found\n"); */
	return E_NO_PIVOT;
    }
}

static int dbl_create_factor_space (dbl_factor_work * f)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_ur_info *ur_inf = f->ur_inf;
    int dim = f->dim;
    int nzcnt;
    int i;
    int rval;

    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	nzcnt += ur_inf[i].nzcnt;
    }

    if (f->ucindx == 0) {
	f->uc_space = nzcnt * f->uc_space_mul;
	ILL_SAFE_MALLOC (f->ucindx, f->uc_space + 1, int);
    }
    if (f->urindx == 0 || f->urcoef == 0) {
	ILL_IFFREE (f->urindx, int);
	dbl_EGlpNumFreeArray (f->urcoef);
	f->ur_space = nzcnt * f->ur_space_mul;
	ILL_SAFE_MALLOC (f->urindx, f->ur_space + 1, int);
	f->urcoef = dbl_EGlpNumAllocArray (f->ur_space);
    }
    if (f->lcindx == 0 || f->lccoef == 0) {
	ILL_IFFREE (f->lcindx, int);
	dbl_EGlpNumFreeArray (f->lccoef);
	f->lc_space = nzcnt * f->lc_space_mul;
	ILL_SAFE_MALLOC (f->lcindx, f->lc_space, int);
	f->lccoef = dbl_EGlpNumAllocArray (f->lc_space);
    }
    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	ur_inf[i].rbeg = nzcnt;
	nzcnt += ur_inf[i].nzcnt;
	ur_inf[i].nzcnt = ur_inf[i].rbeg;
    }
    f->ur_freebeg = nzcnt;

    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	uc_inf[i].cbeg = nzcnt;
	nzcnt += uc_inf[i].nzcnt;
	uc_inf[i].nzcnt = uc_inf[i].cbeg;
    }
    f->uc_freebeg = nzcnt;

    f->lc_freebeg = 0;

    rval = 0;
CLEANUP:
    ILL_RETURN (rval, "dbl_create_factor_space");
}

static int dbl_init_matrix (dbl_factor_work * f,
      int *basis,
      int *cbeg,
      int *clen,
      int *in_ucindx,
      double *in_uccoef)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    dbl_ur_info *ur_inf = f->ur_inf;
    int dim = f->dim;
    int max_k = f->max_k;
    int *ucindx;
    int *urindx;
    double *urcoef;
    int nzcnt;
    int beg;
    int i;
    int j;
    int r;
    int rval = 0;
    double v;
    double max;
    dbl_EGlpNumInitVar (v);
    dbl_EGlpNumInitVar (max);

    for (i = 0; i < dim; i++) {
	ur_inf[i].nzcnt = 0;
    }
    for (i = 0; i < dim; i++) {
	nzcnt = clen[basis[i]];
	beg = cbeg[basis[i]];
	uc_inf[i].nzcnt = nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    r = in_ucindx[beg + j];
	    ur_inf[r].nzcnt++;
	}
    }

    rval = dbl_create_factor_space (f);
    ILL_CLEANUP_IF (rval);

    urindx = f->urindx;
    ucindx = f->ucindx;
    urcoef = f->urcoef;

    for (i = 0; i < dim; i++) {
	nzcnt = clen[basis[i]];
	beg = cbeg[basis[i]];
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumCopy (v, in_uccoef[beg + j]);
	    if (dbl_EGlpNumIsEqual (v, dbl_zeroLpNum, f->fzero_tol))
		continue;
	    r = in_ucindx[beg + j];
	    ucindx[uc_inf[i].nzcnt++] = r;
	    urindx[ur_inf[r].nzcnt] = i;
	    dbl_EGlpNumCopy (urcoef[ur_inf[r].nzcnt], v);
	    ur_inf[r].nzcnt++;
	}
    }

    for (i = 0; i < dim; i++) {
	uc_inf[i].nzcnt -= uc_inf[i].cbeg;
	ur_inf[i].nzcnt -= ur_inf[i].rbeg;
    }

    j = f->uc_space;
    for (i = f->uc_freebeg; i < j; i++) {
	ucindx[i] = -1;
    }
    ucindx[j] = 0;

    j = f->ur_space;
    for (i = f->ur_freebeg; i < j; i++) {
	urindx[i] = -1;
    }
    urindx[j] = 0;

    for (i = 0; i < dim; i++) {
	nzcnt = ur_inf[i].nzcnt;
	ur_inf[i].pivcnt = nzcnt;
	beg = ur_inf[i].rbeg;
	dbl_EGlpNumZero (max);
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumSetToMaxAbs (max, urcoef[beg + j]);
	}
	dbl_EGlpNumCopy (ur_inf[i].max, max);
    }

    for (i = 0; i <= max_k; i++) {
	ur_inf[dim + i].next = dim + i;
	ur_inf[dim + i].prev = dim + i;
	uc_inf[dim + i].next = dim + i;
	uc_inf[dim + i].prev = dim + i;
    }

    for (i = 0; i < dim; i++) {
	nzcnt = uc_inf[i].nzcnt;
	if (nzcnt >= max_k)
	    nzcnt = max_k;
	uc_inf[i].next = uc_inf[dim + nzcnt].next;
	uc_inf[i].prev = dim + nzcnt;
	uc_inf[dim + nzcnt].next = i;
	uc_inf[uc_inf[i].next].prev = i;

	nzcnt = ur_inf[i].pivcnt;
	if (nzcnt >= max_k)
	    nzcnt = max_k;
	ur_inf[i].next = ur_inf[dim + nzcnt].next;
	ur_inf[i].prev = dim + nzcnt;
	ur_inf[dim + nzcnt].next = i;
	ur_inf[ur_inf[i].next].prev = i;
    }

#ifdef dbl_TRACK_FACTOR
    dbl_EGlpNumZero (max);
    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	if (dbl_EGlpNumIsLess (max, ur_inf[i].max))
	    dbl_EGlpNumCopy (max, ur_inf[i].max);
	nzcnt += ur_inf[i].nzcnt;
    }

    dbl_EGlpNumCopy (f->maxelem_orig, max);
    f->nzcnt_orig = nzcnt;
    dbl_EGlpNumCopy (f->maxelem_factor, f->maxelem_orig);
    f->nzcnt_factor = f->nzcnt_orig;
#endif				/* dbl_TRACK_FACTOR */

    /* sentinal for column space */
    ucindx[f->uc_space] = 0;

    dbl_clear_work (f);

CLEANUP:
    dbl_EGlpNumClearVar (max);
    dbl_EGlpNumClearVar (v);
    ILL_RETURN (rval, "dbl_init_matrix");
}

static int dbl_build_iteration_u_data (dbl_factor_work * f)
{
    int dim = f->dim;
    dbl_ur_info *ur_inf = f->ur_inf;
    dbl_uc_info *uc_inf = f->uc_inf;
    double *uccoef = 0;
    int *ucindx = 0;
    int *urindx = f->urindx;
    double *urcoef = f->urcoef;
    int *ucrind = 0;
    int *urcind = 0;
    int nzcnt;
    int beg;
    int cbeg;
    int cnzcnt;
    int uc_space = f->uc_space;
    int er_space;
    int i;
    int j;
    int k;
    int rval;
    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	nzcnt += ur_inf[i].nzcnt;
    }

#ifdef dbl_TRACK_FACTOR
    f->nzcnt_factor = nzcnt;
#endif				/* dbl_TRACK_FACTOR */

    dbl_EGlpNumFreeArray (f->uccoef);
    uccoef = dbl_EGlpNumAllocArray (nzcnt);
    f->uccoef = uccoef;

    ILL_IFFREE (f->ucrind, int);
    ILL_SAFE_MALLOC (ucrind, nzcnt, int);
    f->ucrind = ucrind;

    ILL_IFFREE (f->urcind, int);
    ILL_SAFE_MALLOC (urcind, f->ur_space, int);
    f->urcind = urcind;

    if (uc_space < nzcnt) {
	ILL_IFFREE (f->ucindx, int);
	ILL_SAFE_MALLOC (f->ucindx, nzcnt + 1, int);
    }
    f->uc_space = nzcnt;
    uc_space = nzcnt;
    ucindx = f->ucindx;

    for (i = 0; i < dim; i++) {
	uc_inf[i].nzcnt = 0;
    }

    for (i = 0; i < dim; i++) {
	nzcnt = ur_inf[i].nzcnt;
	beg = ur_inf[i].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    uc_inf[urindx[beg + j]].nzcnt++;
	}
	ur_inf[i].delay = 0;
    }

    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	uc_inf[i].cbeg = nzcnt;
	nzcnt += uc_inf[i].nzcnt;
	uc_inf[i].nzcnt = 0;
	uc_inf[i].delay = 0;
    }

    f->uc_freebeg = nzcnt;
    for (i = nzcnt; i < uc_space; i++) {
	ucindx[i] = -1;
    }
    ucindx[uc_space] = 0;

    for (i = 0; i < dim; i++) {
	nzcnt = ur_inf[i].nzcnt;
	beg = ur_inf[i].rbeg;
	k = urindx[beg];
	cbeg = uc_inf[k].cbeg;
	cnzcnt = uc_inf[k].nzcnt;
	if (cnzcnt != 0) {
	    ucindx[cbeg + cnzcnt] = ucindx[cbeg];
	    dbl_EGlpNumCopy (uccoef[cbeg + cnzcnt], uccoef[cbeg]);
	    ucrind[cbeg + cnzcnt] = ucrind[cbeg];
	    urcind[ur_inf[ucindx[cbeg]].rbeg + ucrind[cbeg]] = cnzcnt;
	}
	ucindx[cbeg] = i;
	dbl_EGlpNumCopy (uccoef[cbeg], urcoef[beg]);
	ucrind[cbeg] = 0;
	urcind[beg] = 0;
	uc_inf[k].nzcnt = cnzcnt + 1;
	for (j = 1; j < nzcnt; j++) {
	    k = urindx[beg + j];
	    cbeg = uc_inf[k].cbeg;
	    cnzcnt = uc_inf[k].nzcnt;
	    ucindx[cbeg + cnzcnt] = i;
	    dbl_EGlpNumCopy (uccoef[cbeg + cnzcnt], urcoef[beg + j]);
	    ucrind[cbeg + cnzcnt] = j;
	    urcind[beg + j] = cnzcnt;
	    uc_inf[k].nzcnt++;
	}
    }

    for (i = 0; i < dim; i++) {
	f->rrank[f->rperm[i]] = i;
    }

    nzcnt = f->ur_space;

    for (i = f->ur_freebeg; i < nzcnt; i++) {
	urindx[i] = -1;
    }
    urindx[nzcnt] = 0;

    dbl_clear_work (f);

    er_space = f->er_space_mul * f->etamax;
    ILL_SAFE_MALLOC (f->er_inf, f->etamax, dbl_er_info);
    ILL_SAFE_MALLOC (f->erindx, er_space, int);
    f->ercoef = dbl_EGlpNumAllocArray (er_space);
    f->etacnt = 0;
    f->er_freebeg = 0;
    f->er_space = er_space;

    rval = 0;

CLEANUP:
    ILL_RETURN (rval, "dbl_build_iteration_u_data");
}

static int dbl_build_iteration_l_data (dbl_factor_work * f)
{
    int dim = f->dim;
    dbl_lc_info *lc_inf = f->lc_inf;
    dbl_lr_info *lr_inf = f->lr_inf;
    double *lrcoef = 0;
    int *lrindx = 0;
    double *lccoef = f->lccoef;
    int *lcindx = f->lcindx;
    int nzcnt;
    int beg;
    int rnzcnt;
    int rbeg;
    int i;
    int j;
    int k;
    int c;
    int rval;

    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	nzcnt += lc_inf[i].nzcnt;
	lr_inf[i].nzcnt = 0;
	lr_inf[i].delay = 0;
	lc_inf[lc_inf[i].c].crank = i;
    }

    dbl_EGlpNumFreeArray (f->lrcoef);
    if (nzcnt) {
	lrcoef = dbl_EGlpNumAllocArray (nzcnt);
	f->lrcoef = lrcoef;
    }
    ILL_IFFREE (f->lrindx, int);
    ILL_SAFE_MALLOC (lrindx, nzcnt + 1, int);
    f->lrindx = lrindx;

    for (i = 0; i < dim; i++) {
	nzcnt = lc_inf[i].nzcnt;
	beg = lc_inf[i].cbeg;
	lc_inf[i].delay = 0;
	for (j = 0; j < nzcnt; j++) {
	    lr_inf[lc_inf[lcindx[beg + j]].crank].nzcnt++;
	}
    }

    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	lr_inf[i].rbeg = nzcnt;
	nzcnt += lr_inf[i].nzcnt;
	lr_inf[i].nzcnt = 0;
	lr_inf[i].r = lc_inf[i].c;
	lr_inf[lr_inf[i].r].rrank = i;
    }

    for (i = 0; i < dim; i++) {
	nzcnt = lc_inf[i].nzcnt;
	beg = lc_inf[i].cbeg;
	c = lc_inf[i].c;
	for (j = 0; j < nzcnt; j++) {
	    k = lc_inf[lcindx[beg + j]].crank;
	    rbeg = lr_inf[k].rbeg;
	    rnzcnt = lr_inf[k].nzcnt;
	    lrindx[rbeg + rnzcnt] = c;
	    dbl_EGlpNumCopy (lrcoef[rbeg + rnzcnt], lccoef[beg + j]);
	    lr_inf[k].nzcnt++;
	}
    }

#ifdef dbl_TRACK_FACTOR
    nzcnt = f->nzcnt_factor;
    for (i = 0; i < dim; i++) {
	nzcnt += lc_inf[i].nzcnt;
    }
    f->nzcnt_factor = nzcnt;

    dbl_EGlpNumCopy (f->maxelem_cur, f->maxelem_factor);
    f->nzcnt_cur = f->nzcnt_factor;

    /*
        dbl_dump_factor_stats (f);
        printf ("orig max  %e nzcnt %d\n", f->maxelem_orig, f->nzcnt_orig);
        printf ("f maxelem %e nzcnt %d\n", f->maxelem_cur, f->nzcnt_cur);
    */
#endif				/* dbl_TRACK_FACTOR */

    rval = 0;

CLEANUP:
    ILL_RETURN (rval, "dbl_build_iteration_l_data");
}

static int dbl_handle_singularity (dbl_factor_work * f)
{
    int rval = 0;
    int nsing;
    int *singr = 0;
    int *singc = 0;
    int i;

    if (f->p_nsing == 0 || f->p_singr == 0 || f->p_singc == 0) {
	fprintf (stderr, "singular basis, but no place for singularity data\n");
	return E_SING_NO_DATA;
    }
    nsing = f->nstages - f->stage;
    ILL_SAFE_MALLOC (singr, nsing, int);
    ILL_SAFE_MALLOC (singc, nsing, int);
    for (i = f->stage; i < f->nstages; i++) {
	singr[i - f->stage] = f->rperm[i];
	singc[i - f->stage] = f->cperm[i];
    }
    *f->p_nsing = nsing;
    *f->p_singr = singr;
    *f->p_singc = singc;
    singr = 0;
    singc = 0;

CLEANUP:
    ILL_IFFREE (singr, int);
    ILL_IFFREE (singc, int);
    ILL_RETURN (rval, "dbl_handle_singularity");
}

static int dbl_dense_build_matrix (dbl_factor_work * f)
{
    double *dmat = 0;
    int stage = f->stage;
    int drows = f->nstages - stage;
    int dcols = f->dim - stage;
    int dsize = drows * dcols;
    int *crank = f->crank;
    double *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int nzcnt;
    int beg;
    int i;
    int r;
    int j;
    int rval = 0;

    dmat = dbl_EGlpNumAllocArray (dsize);

    for (i = 0; i < dsize; i++)
	dbl_EGlpNumZero (dmat[i]);

    for (i = 0; i < drows; i++) {
	r = f->rperm[i + stage];
	nzcnt = f->ur_inf[r].nzcnt;
	beg = f->ur_inf[r].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumCopy (dmat[i * dcols - stage + crank[urindx[beg + j]]],
		urcoef[beg + j]);
	}
    }

    f->drows = drows;
    f->dcols = dcols;
    f->dense_base = f->stage;
    f->dmat = dmat;
    dmat = 0;

    /* CLEANUP: */
    dbl_EGlpNumFreeArray (dmat);
    ILL_RETURN (rval, "dbl_dense_build_matrix");
}

static int dbl_dense_find_pivot (dbl_factor_work * f,
      int *p_r,
      int *p_c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    double *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    double maxval;
    int max_r;
    int max_c;
    int i;
    dbl_EGlpNumInitVar (maxval);
    dbl_EGlpNumZero (maxval);
    max_r = -1;
    for (i = s; i < drows; i++) {
	if (dbl_EGlpNumIsLess (maxval, ur_inf[rperm[dense_base + i]].max)) {
	    dbl_EGlpNumCopy (maxval, ur_inf[rperm[dense_base + i]].max);
	    max_r = i;
	}
    }
    if (max_r == -1) {
	return E_NO_PIVOT;
    }
    dbl_EGlpNumZero (maxval);
    max_c = -1;
    for (i = s; i < drows; i++) {
	dbl_EGlpNumSetToMaxAbsAndDo (maxval, dmat[max_r * dcols + i], max_c = i);
    }
    if (max_c == -1) {
	return E_NO_PIVOT;
    }
    *p_r = max_r;
    *p_c = max_c;

    dbl_EGlpNumClearVar (maxval);
    return 0;
}

static void dbl_dense_swap (dbl_factor_work * f,
      int r,
      int c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    double *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    int i;
    double v;
    dbl_EGlpNumInitVar (v);

    if (r != s) {
	dbl_ILL_SWAP (f->rperm[dense_base + s], f->rperm[dense_base + r], i);
	f->rrank[f->rperm[dense_base + s]] = dense_base + s;
	f->rrank[f->rperm[dense_base + r]] = dense_base + r;
	for (i = 0; i < dcols; i++) {
	    dbl_EGLPNUM_SWAP (dmat[s * dcols + i], dmat[r * dcols + i], v);
	}
    }
    if (c != s) {
	dbl_ILL_SWAP (f->cperm[dense_base + s], f->cperm[dense_base + c], i);
	f->crank[f->cperm[dense_base + s]] = dense_base + s;
	f->crank[f->cperm[dense_base + c]] = dense_base + c;
	for (i = 0; i < drows; i++) {
	    dbl_EGLPNUM_SWAP (dmat[i * dcols + s], dmat[i * dcols + c], v);
	}
    }
    dbl_EGlpNumClearVar (v);
}

static void dbl_dense_elim (dbl_factor_work * f,
      int r,
      int c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    double *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int i;
    int j;
    double pivval;
    double max;
    double v;
    double w;
#ifdef dbl_TRACK_FACTOR
    double maxelem_factor;
    dbl_EGlpNumInitVar (maxelem_factor);
    dbl_EGlpNumCopy (maxelem_factor, f->maxelem_factor);
#endif
    dbl_EGlpNumInitVar (pivval);
    dbl_EGlpNumInitVar (max);
    dbl_EGlpNumInitVar (v);
    dbl_EGlpNumInitVar (w);

    dbl_dense_swap (f, r, c);
    f->stage++;
    dbl_EGlpNumCopyFrac (pivval, dbl_oneLpNum, dmat[s * dcols + s]);
    for (i = s + 1; i < drows; i++) {
	dbl_EGlpNumCopy (v, dmat[i * dcols + s]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    dbl_EGlpNumMultTo (v, pivval);
	    if (dbl_EGlpNumIsNeqZero (v, f->fzero_tol)) {
		dbl_EGlpNumCopy (dmat[i * dcols + s], v);
#ifdef dbl_TRACK_FACTOR
		dbl_EGlpNumSetToMaxAbs (maxelem_factor, v);
#endif
		dbl_EGlpNumZero (max);
		for (j = s + 1; j < drows; j++) {
		    dbl_EGlpNumCopy (w, dmat[i * dcols + j]);
		    dbl_EGlpNumSubInnProdTo (w, v, dmat[s * dcols + j]);
		    dbl_EGlpNumCopy (dmat[i * dcols + j], w);
		    dbl_EGlpNumSetToMaxAbs (max, w);
		}
		for (j = drows; j < dcols; j++) {
		    dbl_EGlpNumCopy (w, dmat[i * dcols + j]);
		    dbl_EGlpNumSubInnProdTo (w, v, dmat[s * dcols + j]);
		    dbl_EGlpNumCopy (dmat[i * dcols + j], w);
		}
		dbl_EGlpNumCopy (ur_inf[rperm[dense_base + i]].max, max);
#ifdef dbl_TRACK_FACTOR
		if (dbl_EGlpNumIsLess (maxelem_factor, max))
		    dbl_EGlpNumCopy (maxelem_factor, max);
#endif
	    } else {
		dbl_EGlpNumZero (dmat[i * dcols + s]);
	    }
	}
    }
#ifdef dbl_TRACK_FACTOR
    dbl_EGlpNumCopy (f->maxelem_factor, maxelem_factor);
    dbl_EGlpNumClearVar (maxelem_factor);
#endif
    dbl_EGlpNumClearVar (pivval);
    dbl_EGlpNumClearVar (max);
    dbl_EGlpNumClearVar (v);
    dbl_EGlpNumClearVar (w);
}

static int dbl_dense_replace_row (dbl_factor_work * f,
      int i)
{
    int dcols = f->dcols;
    int dense_base = f->dense_base;
    double *dmat = f->dmat + i * dcols;
    double *urcoef;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *cperm = f->cperm;
    int r = f->rperm[dense_base + i];
    int *urindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j = i; j < dcols; j++) {
	if (dbl_EGlpNumIsNeqZero (dmat[j], f->fzero_tol)) {
	    nzcnt++;
	}
    }
    if (nzcnt > ur_inf[r].nzcnt) {
	if (ur_inf[r].rbeg + ur_inf[r].nzcnt == f->ur_freebeg) {
	    f->ur_freebeg = ur_inf[r].rbeg;
	}
	ur_inf[r].nzcnt = 0;
	if (f->ur_freebeg + nzcnt > f->ur_space) {
	    rval = dbl_make_ur_space (f, nzcnt);
	    ILL_CLEANUP_IF (rval);
	}
	ur_inf[r].rbeg = f->ur_freebeg;
	f->ur_freebeg += nzcnt;
    }
    beg = ur_inf[r].rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    for (j = i; j < dcols; j++) {
	if (dbl_EGlpNumIsNeqZero (dmat[j], f->fzero_tol)) {
	    dbl_EGlpNumCopy (urcoef[beg], dmat[j]);
	    urindx[beg] = cperm[dense_base + j];
	    beg++;
	}
    }
    ur_inf[r].nzcnt = beg - ur_inf[r].rbeg;
CLEANUP:
    ILL_RETURN (rval, "dbl_dense_replace_row");
}

static int dbl_dense_create_col (dbl_factor_work * f,
      int i)
{
    int dcols = f->dcols;
    int drows = f->drows;
    int dense_base = f->dense_base;
    double *dmat = f->dmat;
    double *lccoef;
    dbl_lc_info *lc_inf = f->lc_inf;
    int *rperm = f->rperm;
    int *lcindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j = i + 1; j < drows; j++) {
	if (dbl_EGlpNumIsNeqZero (dmat[j * dcols + i], f->fzero_tol)) {
	    nzcnt++;
	}
    }

    if (f->lc_freebeg + nzcnt >= f->lc_space) {
	rval = dbl_make_lc_space (f, nzcnt);
	ILL_CLEANUP_IF (rval);
    }
    beg = f->lc_freebeg;
    lc_inf[dense_base + i].cbeg = beg;
    lc_inf[dense_base + i].c = rperm[dense_base + i];
    lcindx = f->lcindx;
    lccoef = f->lccoef;

    for (j = i + 1; j < drows; j++) {
	if (dbl_EGlpNumIsNeqZero (dmat[j * dcols + i], f->fzero_tol)) {
	    dbl_EGlpNumCopy (lccoef[beg], dmat[j * dcols + i]);
	    lcindx[beg] = rperm[dense_base + j];
	    beg++;
	}
    }
    lc_inf[dense_base + i].nzcnt = beg - lc_inf[dense_base + i].cbeg;
    f->lc_freebeg = beg;
CLEANUP:
    ILL_RETURN (rval, "create_col");
}

static int dbl_dense_replace (dbl_factor_work * f)
{
    int drows = f->drows;
    int rval = 0;
    int i;

    for (i = 0; i < drows; i++) {
	rval = dbl_dense_replace_row (f, i);
	ILL_CLEANUP_IF (rval);
	rval = dbl_dense_create_col (f, i);
	ILL_CLEANUP_IF (rval);
    }
    dbl_EGlpNumFreeArray (f->dmat);
    f->drows = 0;
    f->dcols = 0;
CLEANUP:
    ILL_RETURN (rval, "dbl_dense_replace");
}

static int dbl_dense_factor (dbl_factor_work * f)
{
    int r;
    int c;
    int rval = 0;
#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
    double tmpsize;
#endif
#endif

    /*
        printf ("dense kernel, %d rows, %d  cols...\n", f->nstages - f->stage,
                f->dim - f->stage);
        fflush (stdout);
    */

    rval = dbl_dense_build_matrix (f);
    ILL_CLEANUP_IF (rval);

#ifdef dbl_FACTOR_DEBUG
#if (dbl_FACTOR_DEBUG+0>1)
    printf ("before Dense dbl_ILLfactor\n");
    dbl_dump_matrix (f, 1);
    fflush (stdout);
#endif
#endif

    while (f->stage < f->nstages) {
	r = f->stage - f->dense_base;
	rval = dbl_dense_find_pivot (f, &r, &c);
	if (rval == E_NO_PIVOT) {
	    rval = dbl_handle_singularity (f);
	    ILL_CLEANUP_IF (rval);
	    return E_SINGULAR_INTERNAL;
	} else {
	    ILL_CLEANUP_IF (rval);
	}
#ifdef dbl_FACTOR_DEBUG
#if (dbl_FACTOR_DEBUG+0>2)
	printf ("dense pivot elem: %d %d\n", r, c);
	fflush (stdout);
#endif
#endif				/* dbl_FACTOR_DEBUG */
	dbl_dense_elim (f, r, c);

#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
	tmpsize = f->maxmult * dbl_EGlpNumToLf (f->maxelem_orig);
	if (tmpsize < dbl_EGlpNumToLf (f->maxelem_factor) &&
	    dbl_EGlpNumIsLess (f->partial_cur, dbl_oneLpNum)) {
	    return E_FACTOR_BLOWUP;
	}
#endif				/* dbl_NOTICE_BLOWUP */
#endif				/* dbl_TRACK_FACTOR */

#ifdef dbl_FACTOR_DEBUG
#if (dbl_FACTOR_DEBUG+0>1)
	printf ("After dense pivot stage %d (%d) of %d (%d)\n",
	    f->stage - f->dense_base, f->stage,
	    f->nstages - f->dense_base, f->nstages);
	fflush (stdout);
#endif
#if (dbl_FACTOR_DEBUG+0>2)
	dbl_dump_matrix (f, 1);
	fflush (stdout);
#endif
#endif				/* dbl_FACTOR_DEBUG */
    }

#ifdef dbl_FACTOR_DEBUG
    printf ("After dense dbl_ILLfactor:\n");
    dbl_dump_matrix (f, 0);
    fflush (stdout);
#endif				/* dbl_FACTOR_DEBUG */

    rval = dbl_dense_replace (f);
    ILL_CLEANUP_IF (rval);

#ifdef dbl_FACTOR_DEBUG
    printf ("After replacement:\n");
    dbl_dump_matrix (f, 0);
    fflush (stdout);
#endif				/* dbl_FACTOR_DEBUG */

CLEANUP:
    ILL_RETURN (rval, "dbl_dense_factor");
}

#ifdef dbl_RECORD
FILE *dbl_fsave = 0;
int dbl_fsavecnt = 0;
#endif				/* dbl_RECORD */

static int dbl_ILLfactor_try (dbl_factor_work * f,
      int *basis,
      int *cbeg,
      int *clen,
      int *cindx,
      double *ccoef)
{
    int rval = 0;
    int r;
    int c;
#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
    double tmpsize;
    dbl_EGlpNumInitVar (tmpsize);
#endif
#endif

#ifdef dbl_RECORD
    {
	int ncol = 0;
	int nzcnt = 0;
	int dim = f->dim;
	int i;
	int j;
	char fnambuf[40];

	for (i = 0; i < dim; i++) {
	    if (basis[i] > ncol)
		ncol = basis[i];
	}
	ncol++;
	for (i = 0; i < ncol; i++) {
	    nzcnt += clen[i];
	}
	if (dbl_fsave)
	    fclose (dbl_fsave);
	sprintf (fnambuf, "prob.mat.%d", dbl_fsavecnt);
	dbl_fsavecnt++;
	dbl_fsave = fopen (fnambuf, "w");
	fprintf (dbl_fsave, "%d %d %d\n", f->dim, ncol, nzcnt);
	for (i = 0; i < dim; i++) {
	    fprintf (dbl_fsave, "%d ", basis[i]);
	}
	fprintf (dbl_fsave, "\n");
	for (i = 0; i < ncol; i++) {
	    fprintf (dbl_fsave, "%d", clen[i]);
	    for (j = 0; j < clen[i]; j++) {
		fprintf (dbl_fsave, " %d %.16lg", cindx[cbeg[i] + j],
		    dbl_EGlpNumToLf (ccoef[cbeg[i] + j]));
	    }
	    fprintf (dbl_fsave, "\n");
	}
	fprintf (dbl_fsave, "\n");
	fflush (dbl_fsave);
    }
#endif				/* dbl_RECORD */

    rval = dbl_init_matrix (f, basis, cbeg, clen, cindx, ccoef);
    ILL_CLEANUP_IF (rval);

    f->stage = 0;
    f->nstages = f->dim;

#ifdef dbl_FACTOR_DEBUG
    printf ("Initial matrix:\n");
#if (dbl_FACTOR_DEBUG+0>1)
    dbl_dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif				/* dbl_FACTOR_DEBUG */
#ifdef dbl_FACTOR_STATS
    printf ("Initial matrix: ");
    dbl_dump_factor_stats (f);
#endif				/* dbl_FACTOR_STATS */

    while (f->stage < f->nstages) {
	rval = dbl_find_pivot (f, &r, &c);
	if (rval == E_NO_PIVOT) {
	    rval = dbl_handle_singularity (f);
	    ILL_CLEANUP_IF (rval);
	    return 0;
	} else {
	    ILL_CLEANUP_IF (rval);
	}
	if (f->ur_inf[r].pivcnt > f->dense_fract * (f->nstages - f->stage) &&
	    f->uc_inf[c].nzcnt > f->dense_fract * (f->nstages - f->stage) &&
	    f->nstages - f->stage > f->dense_min) {
	    rval = dbl_dense_factor (f);
	    if (rval == E_SINGULAR_INTERNAL)
		return 0;
	    if (rval)
		return rval;
	    break;
	}
#ifdef dbl_FACTOR_DEBUG
	printf ("pivot elem: %d %d\n", r, c);
	fflush (stdout);
#endif				/* dbl_FACTOR_DEBUG */
	rval = dbl_elim (f, r, c);
	ILL_CLEANUP_IF (rval);

#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
	dbl_EGlpNumSet (tmpsize, f->maxmult);
	dbl_EGlpNumMultTo (tmpsize, f->maxelem_orig);
	if (dbl_EGlpNumIsLess (tmpsize, f->maxelem_factor) &&
	    dbl_EGlpNumIsLess (f->partial_cur, dbl_oneLpNum)) {
	    return E_FACTOR_BLOWUP;
	}
#endif				/* dbl_NOTICE_BLOWUP */
#endif				/* dbl_TRACK_FACTOR */

#ifdef dbl_FACTOR_DEBUG
#if (dbl_FACTOR_DEBUG+0>3)
	printf ("After pivot stage %d of %d\n", f->stage, f->nstages);
	dbl_dump_matrix (f, 0);
	fflush (stdout);
#endif
#endif				/* dbl_FACTOR_DEBUG */
    }

    rval = dbl_build_iteration_u_data (f);
    ILL_CLEANUP_IF (rval);

    rval = dbl_build_iteration_l_data (f);
    ILL_CLEANUP_IF (rval);

#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
    dbl_EGlpNumSet (tmpsize, f->minmult);
    dbl_EGlpNumMultTo (tmpsize, f->maxelem_orig);
    if (dbl_EGlpNumIsLess (f->maxelem_factor, tmpsize) &&
	dbl_EGlpNumIsLess (f->partial_tol, f->partial_cur)) {
	if (dbl_EGlpNumIsGreaDbl (f->partial_cur, 0.5)) {
	    dbl_EGlpNumSet (f->partial_cur, 0.5);
	} else if (dbl_EGlpNumIsGreaDbl (f->partial_cur, 0.25)) {
	    dbl_EGlpNumSet (f->partial_cur, 0.25);
	} else if (dbl_EGlpNumIsGreaDbl (f->partial_cur, 0.1)) {
	    dbl_EGlpNumSet (f->partial_cur, 0.1);
	} else {
	    dbl_EGlpNumDivUiTo (f->partial_cur, 10);
	}
	if (dbl_EGlpNumIsLess (f->partial_cur, f->partial_tol)) {
	    dbl_EGlpNumCopy (f->partial_cur, f->partial_tol);
	}
	/* Bico - comment out for dist fprintf (stderr, "factor good,
	   lowering partial tolerance to %.2f\n", f->partial_cur); */
    }
#endif				/* dbl_NOTICE_BLOWUP */
#endif				/* dbl_TRACK_FACTOR */

#ifdef dbl_FACTOR_DEBUG
    printf ("Factored matrix:\n");
#if (dbl_FACTOR_DEBUG+0>1)
    dbl_dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif				/* dbl_FACTOR_DEBUG */

#ifdef dbl_FACTOR_STATS
    printf ("Factored matrix: ");
    dbl_dump_factor_stats (f);
#endif				/* dbl_FACTOR_STATS */
CLEANUP:
#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
    dbl_EGlpNumClearVar (tmpsize);
#endif
#endif
    ILL_RETURN (rval, "factor_try");
}

int dbl_ILLfactor (dbl_factor_work * f,
      int *basis,
      int *cbeg,
      int *clen,
      int *cindx,
      double *ccoef,
      int *p_nsing,
      int **p_singr,
      int **p_singc)
{
    int rval;
    f->p_nsing = p_nsing;
    f->p_singr = p_singr;
    f->p_singc = p_singc;
    *p_nsing = 0;

AGAIN:
    rval = dbl_ILLfactor_try (f, basis, cbeg, clen, cindx, ccoef);
    if (rval == E_FACTOR_BLOWUP) {
	if (dbl_EGlpNumIsLessDbl (f->partial_cur, 0.1)) {
	    dbl_EGlpNumMultUiTo (f->partial_cur, 10);
	} else if (dbl_EGlpNumIsLessDbl (f->partial_cur, 0.25)) {
	    dbl_EGlpNumSet (f->partial_cur, 0.25);
	} else if (dbl_EGlpNumIsLessDbl (f->partial_cur, 0.5)) {
	    dbl_EGlpNumSet (f->partial_cur, 0.5);
	} else if (dbl_EGlpNumIsLess (f->partial_cur, dbl_oneLpNum)) {
	    dbl_EGlpNumOne (f->partial_cur);
	} else {
	    ILL_RETURN (rval, "dbl_ILLfactor");
	}
	/* Bico - comment out for dist fprintf (stderr, "factor blowup,
	   changing partial tolerance to %.2f\n", f->partial_cur); */
	goto AGAIN;
    }
    ILL_RETURN (rval, "dbl_ILLfactor");
}

static void dbl_ILLfactor_ftranl (dbl_factor_work * f,
      double *a)
{
    int *lcindx = f->lcindx;
    dbl_lc_info *lc_inf = f->lc_inf;
    double *lccoef = f->lccoef;
    int dim = f->dim;
    int beg;
    int nzcnt;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = 0; i < dim; i++) {
	dbl_EGlpNumCopy (v, a[lc_inf[i].c]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    nzcnt = lc_inf[i].nzcnt;
	    beg = lc_inf[i].cbeg;
	    for (j = 0; j < nzcnt; j++) {
		dbl_EGlpNumSubInnProdTo (a[lcindx[beg + j]], v, lccoef[beg + j]);
	    }
	}
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 > 1)
	printf ("dbl_ILLfactor_ftran a after l %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", dbl_EGlpNumToLf (a[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    }
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 <= 1)
    printf ("dbl_ILLfactor_ftran a after l:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", dbl_EGlpNumToLf (a[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}

static void dbl_ftranl3_delay (dbl_factor_work * f,
      int c)
{
    dbl_lc_info *lc_inf = f->lc_inf;
    int nzcnt;
    int *indx;
    int i;

    c = lc_inf[c].crank;
    nzcnt = lc_inf[c].nzcnt;
    indx = f->lcindx + lc_inf[c].cbeg;
    for (i = 0; i < nzcnt; i++) {
	c = indx[i];
	if (lc_inf[c].delay++ == 0) {
	    dbl_ftranl3_delay (f, c);
	}
    }
}

static void dbl_ftranl3_delay2 (dbl_factor_work * f,
      int c)
{
    dbl_lc_info *lc_inf = f->lc_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
	c = lc_inf[c].crank;
	nzcnt = lc_inf[c].nzcnt;
	indx = f->lcindx + lc_inf[c].cbeg;
	last = -1;
	for (i = 0; i < nzcnt; i++) {
	    c = indx[i];
	    if (lc_inf[c].delay++ == 0) {
		if (last >= 0) {
		    dbl_ftranl3_delay2 (f, last);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
}

static void dbl_ftranl3_process (dbl_factor_work * f,
      int c,
      dbl_svector * x)
{
    dbl_lc_info *lc_inf = f->lc_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    int i;
    double *coef;
    double v;
    dbl_EGlpNumInitVar (v);

    dbl_EGlpNumCopy (v, work[c]);
    dbl_EGlpNumZero (work[c]);
    if (dbl_EGlpNumIsNeqqZero (v)) {
	x->indx[x->nzcnt] = c;
	dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    c = lc_inf[c].crank;
    nzcnt = lc_inf[c].nzcnt;
    indx = f->lcindx + lc_inf[c].cbeg;
    coef = f->lccoef + lc_inf[c].cbeg;
    for (i = 0; i < nzcnt; i++) {
	c = indx[i];
	dbl_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	if (--lc_inf[c].delay == 0) {
	    dbl_ftranl3_process (f, c, x);
	}
    }
    dbl_EGlpNumClearVar (v);
}

static void dbl_ftranl3_process2 (dbl_factor_work * f,
      int c,
      dbl_svector * x)
{
    dbl_lc_info *lc_inf = f->lc_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    int last;
    double v;
    dbl_EGlpNumInitVar (v);

    do {
	dbl_EGlpNumCopy (v, work[c]);
	dbl_EGlpNumZero (work[c]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    x->indx[x->nzcnt] = c;
	    dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	c = lc_inf[c].crank;
	nzcnt = lc_inf[c].nzcnt;
	indx = f->lcindx + lc_inf[c].cbeg;
	coef = f->lccoef + lc_inf[c].cbeg;
	last = -1;
	for (i = 0; i < nzcnt; i++) {
	    c = indx[i];
	    dbl_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	    if (--lc_inf[c].delay == 0) {
		if (last >= 0) {
		    dbl_ftranl3_process2 (f, last, x);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_ftranl3 (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * x)
{
    double *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    double *acoef = a->coef;
    dbl_lc_info *lc_inf = f->lc_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (lc_inf[aindx[i]].delay++ == 0) {
	    dbl_ftranl3_delay2 (f, aindx[i]);
	}
	dbl_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--lc_inf[aindx[i]].delay == 0) {
	    dbl_ftranl3_process2 (f, aindx[i], x);
	}
    }
#ifdef dbl_SOLVE_DEBUG
    printf ("dbl_ILLfactor_ftran x after l3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* dbl_SOLVE_DEBUG */
}

static void dbl_ILLfactor_ftrane (dbl_factor_work * f,
      double *a)
{
    int *erindx = f->erindx;
    double *ercoef = f->ercoef;
    dbl_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = 0; i < etacnt; i++) {
	dbl_EGlpNumCopy (v, a[er_inf[i].r]);
	nzcnt = er_inf[i].nzcnt;
	beg = er_inf[i].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumSubInnProdTo (v, ercoef[beg + j], a[erindx[beg + j]]);
	}
	dbl_EGlpNumCopy (a[er_inf[i].r], v);
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 > 1)
	printf ("dbl_ILLfactor_ftran a after eta %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", dbl_EGlpNumToLf (a[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    }
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 <= 1)
    printf ("dbl_ILLfactor_ftran a after eta:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", dbl_EGlpNumToLf (a[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_ftrane2 (dbl_factor_work * f,
      dbl_svector * a)
{
    int *erindx = f->erindx;
    double *ercoef = f->ercoef;
    dbl_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    double *acoef = a->coef;
    double *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;
    int r;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = 0; i < anzcnt; i++) {
	dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	work_indx[aindx[i]] = i + 1;
    }
    for (i = 0; i < etacnt; i++) {
	r = er_inf[i].r;
	dbl_EGlpNumCopy (v, work_coef[r]);
	nzcnt = er_inf[i].nzcnt;
	beg = er_inf[i].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    dbl_EGlpNumSubInnProdTo (v, ercoef[beg + j], work_coef[erindx[beg + j]]);
	}
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    dbl_EGlpNumCopy (work_coef[r], v);
	    if (work_indx[r] == 0) {
		dbl_EGlpNumCopy (acoef[anzcnt], v);
		aindx[anzcnt] = r;
		work_indx[r] = anzcnt + 1;
		anzcnt++;
	    } else {
		dbl_EGlpNumCopy (acoef[work_indx[r] - 1], v);
	    }
	} else {
	    dbl_EGlpNumZero (work_coef[r]);
	    if (work_indx[r]) {
		dbl_EGlpNumZero (acoef[work_indx[r] - 1]);
	    }
	}
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 > 1)
	printf ("dbl_ILLfactor_ftran a after eta2 %d:", i);
	for (j = 0; j < anzcnt; j++) {
	    printf (" %.3f*%d", dbl_EGlpNumToLf (acoef[j]), aindx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    }
    i = 0;
    while (i < anzcnt) {
	dbl_EGlpNumZero (work_coef[aindx[i]]);
	work_indx[aindx[i]] = 0;
	if (dbl_EGlpNumIsNeqZero (acoef[i], f->fzero_tol)) {
	    /* if (acoef[i] > fzero_tol || acoef[i] < -fzero_tol) */
	    i++;
	} else {
	    --anzcnt;
	    dbl_EGlpNumCopy (acoef[i], acoef[anzcnt]);
	    aindx[i] = aindx[anzcnt];
	}
    }
    a->nzcnt = anzcnt;

#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 <= 1)
    printf ("dbl_ILLfactor_ftran a after eta2:");
    for (j = 0; j < anzcnt; j++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (acoef[j]), aindx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_ftranu (dbl_factor_work * f,
      double *a,
      dbl_svector * x)
{
    int *ucindx = f->ucindx;
    double *uccoef = f->uccoef;
    dbl_uc_info *uc_inf = f->uc_inf;
    int *cperm = f->cperm;
    int *rperm = f->rperm;
    int dim = f->dim;
    int xnzcnt = 0;
    int *xindx = x->indx;
    double *xcoef = x->coef;
    int nzcnt;
    int beg;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = dim - 1; i >= 0; i--) {
	dbl_EGlpNumCopy (v, a[rperm[i]]);
	if (dbl_EGlpNumIsNeqqZero (v)) {	/* ((v = a[rperm[i]]) != 0.0) */
	    j = cperm[i];
	    beg = uc_inf[j].cbeg;
	    dbl_EGlpNumDivTo (v, uccoef[beg]);
	    if (dbl_EGlpNumIsNeqZero (v, f->szero_tol)) {
		/* if (v > szero_tol || v < -szero_tol) */
		xindx[xnzcnt] = j;
		dbl_EGlpNumCopy (xcoef[xnzcnt], v);
		xnzcnt++;
	    }
	    nzcnt = uc_inf[j].nzcnt;
	    for (j = 1; j < nzcnt; j++) {
		dbl_EGlpNumSubInnProdTo (a[ucindx[beg + j]], v, uccoef[beg + j]);
	    }
	    dbl_EGlpNumZero (a[rperm[i]]);
	}
    }
    x->nzcnt = xnzcnt;
#ifdef dbl_SOLVE_DEBUG
    printf ("dbl_ILLfactor_ftran x after u:");
    for (j = 0; j < x->nzcnt; j++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (x->coef[j]), x->indx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}


static void dbl_ftranu3_delay (dbl_factor_work * f,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int nzcnt;
    int *indx;
    int i;

    c = f->cperm[f->rrank[c]];
    nzcnt = uc_inf[c].nzcnt;
    indx = f->ucindx + uc_inf[c].cbeg;
    for (i = 1; i < nzcnt; i++) {
	c = indx[i];
	if (uc_inf[c].delay++ == 0) {
	    dbl_ftranu3_delay (f, c);
	}
    }
}

static void dbl_ftranu3_delay2 (dbl_factor_work * f,
      int c)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
	c = f->cperm[f->rrank[c]];
	nzcnt = uc_inf[c].nzcnt;
	indx = f->ucindx + uc_inf[c].cbeg;
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    c = indx[i];
	    if (uc_inf[c].delay++ == 0) {
		if (last >= 0) {
		    dbl_ftranu3_delay2 (f, last);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
}

static void dbl_ftranu3_process (dbl_factor_work * f,
      int c,
      dbl_svector * x)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    double v;
    dbl_EGlpNumInitVar (v);

    dbl_EGlpNumCopy (v, work[c]);
    dbl_EGlpNumZero (work[c]);
    c = f->cperm[f->rrank[c]];
    nzcnt = uc_inf[c].nzcnt;
    indx = f->ucindx + uc_inf[c].cbeg;
    coef = f->uccoef + uc_inf[c].cbeg;
    dbl_EGlpNumDivTo (v, coef[0]);
    if (dbl_EGlpNumIsNeqZero (v, f->szero_tol))
	/* if (v > szero_tol || v < -szero_tol) */
    {
	x->indx[x->nzcnt] = c;
	dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    for (i = 1; i < nzcnt; i++) {
	c = indx[i];
	dbl_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	if (--uc_inf[c].delay == 0) {
	    dbl_ftranu3_process (f, c, x);
	}
    }
    dbl_EGlpNumClearVar (v);
}

static void dbl_ftranu3_process2 (dbl_factor_work * f,
      int c,
      dbl_svector * x)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    int last;
    double v;
    dbl_EGlpNumInitVar (v);

    do {
	dbl_EGlpNumCopy (v, work[c]);
	dbl_EGlpNumZero (work[c]);
	c = f->cperm[f->rrank[c]];
	nzcnt = uc_inf[c].nzcnt;
	indx = f->ucindx + uc_inf[c].cbeg;
	coef = f->uccoef + uc_inf[c].cbeg;
	dbl_EGlpNumDivTo (v, coef[0]);
	if (dbl_EGlpNumIsNeqZero (v, f->szero_tol))
	    /* if (v > szero_tol || v < -szero_tol) */
	{
	    x->indx[x->nzcnt] = c;
	    dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    c = indx[i];
	    dbl_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	    if (--uc_inf[c].delay == 0) {
		if (last >= 0) {
		    dbl_ftranu3_process2 (f, last, x);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_ftranu3 (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * x)
{
    double *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    double *acoef = a->coef;
    dbl_uc_info *uc_inf = f->uc_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (uc_inf[aindx[i]].delay++ == 0) {
	    dbl_ftranu3_delay2 (f, aindx[i]);
	}
	dbl_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--uc_inf[aindx[i]].delay == 0) {
	    dbl_ftranu3_process2 (f, aindx[i], x);
	}
    }
#ifdef dbl_SOLVE_DEBUG
    printf ("dbl_ILLfactor_ftran x after u3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* dbl_SOLVE_DEBUG */
}

/* dbl_ILLfactor_ftran solves Bx=a for x */
void dbl_ILLfactor_ftran (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx;
    double *acoef;
    double *work_coef = f->work_coef;

#ifdef dbl_RECORD
    {
	fprintf (dbl_fsave, "f %d", a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (dbl_fsave, " %d %.16e", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	fprintf (dbl_fsave, "\n");
	fflush (dbl_fsave);
    }
#endif				/* dbl_RECORD */
#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_ftran a:");
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %d %la", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
	nzcnt = a->nzcnt;
	aindx = a->indx;
	acoef = a->coef;
	for (i = 0; i < nzcnt; i++) {
	    dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	dbl_ILLfactor_ftranl3 (f, a, &f->xtmp);
	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = f->xtmp.nzcnt;
	    aindx = f->xtmp.indx;
	    acoef = f->xtmp.coef;

	    for (i = 0; i < nzcnt; i++) {
		dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	dbl_ILLfactor_ftranl (f, work_coef);
    }

    if (sparse) {
	dbl_ILLfactor_ftrane2 (f, &f->xtmp);
	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = f->xtmp.nzcnt;
	    aindx = f->xtmp.indx;
	    acoef = f->xtmp.coef;

	    for (i = 0; i < nzcnt; i++) {
		dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	dbl_ILLfactor_ftrane (f, work_coef);
    }

    if (sparse) {
	dbl_ILLfactor_ftranu3 (f, &f->xtmp, x);
    } else {
	dbl_ILLfactor_ftranu (f, work_coef, x);
    }

#ifdef dbl_SORT_RESULTS
    dbl_sort_vector (x);
#endif

#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_ftran x:");
	for (i = 0; i < x->nzcnt; i++) {
	    printf (" %d %la", x->indx[i], dbl_EGlpNumToLf (x->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */
    return;
}

/* dbl_ILLfactor_ftran_update solves Bx=a for x, and also returns upd, where
   Ux=upd */
void dbl_ILLfactor_ftran_update (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * upd,
      dbl_svector * x)
{
    int i;
    int nzcnt;
    int dim;
    int sparse;
    int *aindx;
    double *acoef;
    double *work_coef = f->work_coef;

#ifdef dbl_RECORD
    {
	fprintf (dbl_fsave, "F %d", a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (dbl_fsave, " %d %.16e", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	fprintf (dbl_fsave, "\n");
	fflush (dbl_fsave);
    }
#endif				/* dbl_RECORD */
#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_ftran_update a:");
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %d %.3f", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
	aindx = a->indx;
	acoef = a->coef;
	nzcnt = a->nzcnt;

	for (i = 0; i < nzcnt; i++) {
	    dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	dbl_ILLfactor_ftranl3 (f, a, upd);
	if (upd->nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = upd->nzcnt;
	    aindx = upd->indx;
	    acoef = upd->coef;

	    for (i = 0; i < nzcnt; i++) {
		dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	dbl_ILLfactor_ftranl (f, work_coef);
    }

    if (sparse) {
	dbl_ILLfactor_ftrane2 (f, upd);
	if (upd->nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = upd->nzcnt;
	    aindx = upd->indx;
	    acoef = upd->coef;

	    for (i = 0; i < nzcnt; i++) {
		dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	dbl_ILLfactor_ftrane (f, work_coef);
	nzcnt = 0;
	dim = f->dim;
	aindx = upd->indx;
	acoef = upd->coef;
	for (i = 0; i < dim; i++) {
	    if (dbl_EGlpNumIsNeqqZero (work_coef[i])) {
		if (dbl_EGlpNumIsNeqZero (work_coef[i], f->szero_tol))
		    /* if(work_coef[i] > szero_tol || work_coef[i] <
		       -szero_tol) */
		{
		    aindx[nzcnt] = i;
		    dbl_EGlpNumCopy (acoef[nzcnt], work_coef[i]);
		    nzcnt++;
		}
	    }
	}
	upd->nzcnt = nzcnt;
    }

    if (sparse) {
	dbl_ILLfactor_ftranu3 (f, upd, x);
    } else {
	dbl_ILLfactor_ftranu (f, work_coef, x);
    }

#ifdef dbl_SORT_RESULTS
    dbl_sort_vector (upd);
    dbl_sort_vector (x);
#endif

#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_ftran update x:");
	for (i = 0; i < x->nzcnt; i++) {
	    printf (" %d %.3f", x->indx[i], dbl_EGlpNumToLf (x->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */
}


static void dbl_ILLfactor_btranl2 (dbl_factor_work * f,
      double *x)
{
    int *lrindx = f->lrindx;
    double *lrcoef = f->lrcoef;
    dbl_lr_info *lr_inf = f->lr_inf;
    int dim = f->dim;
    int nzcnt;
    int beg;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = dim - 1; i >= 0; i--) {
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 > 1)
	printf ("dbl_ILLfactor_btran x before l2 %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", dbl_EGlpNumToLf (x[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
	dbl_EGlpNumCopy (v, x[lr_inf[i].r]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    nzcnt = lr_inf[i].nzcnt;
	    beg = lr_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		dbl_EGlpNumSubInnProdTo (x[lrindx[beg + j]], v, lrcoef[beg + j]);
	    }
	}
    }
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 <= 1)
    printf ("dbl_ILLfactor_btran x after l2:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", dbl_EGlpNumToLf (x[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}

static void dbl_btranl3_delay (dbl_factor_work * f,
      int r)
{
    dbl_lr_info *lr_inf = f->lr_inf;
    int nzcnt;
    int *indx;
    int i;

    r = lr_inf[r].rrank;
    nzcnt = lr_inf[r].nzcnt;
    indx = f->lrindx + lr_inf[r].rbeg;
    for (i = 0; i < nzcnt; i++) {
	r = indx[i];
	if (lr_inf[r].delay++ == 0) {
	    dbl_btranl3_delay (f, r);
	}
    }
}

static void dbl_btranl3_delay2 (dbl_factor_work * f,
      int r)
{
    dbl_lr_info *lr_inf = f->lr_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
	r = lr_inf[r].rrank;
	nzcnt = lr_inf[r].nzcnt;
	indx = f->lrindx + lr_inf[r].rbeg;
	last = -1;
	for (i = 0; i < nzcnt; i++) {
	    r = indx[i];
	    if (lr_inf[r].delay++ == 0) {
		if (last >= 0) {
		    dbl_btranl3_delay2 (f, last);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
}

static void dbl_btranl3_process (dbl_factor_work * f,
      int r,
      dbl_svector * x)
{
    dbl_lr_info *lr_inf = f->lr_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    double v;
    dbl_EGlpNumInitVar (v);

    dbl_EGlpNumCopy (v, work[r]);
    dbl_EGlpNumZero (work[r]);
    if (dbl_EGlpNumIsNeqZero (v, f->szero_tol))
	/* if (v > szero_tol || v < -szero_tol) */
    {
	x->indx[x->nzcnt] = r;
	dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    r = lr_inf[r].rrank;
    nzcnt = lr_inf[r].nzcnt;
    indx = f->lrindx + lr_inf[r].rbeg;
    coef = f->lrcoef + lr_inf[r].rbeg;
    for (i = 0; i < nzcnt; i++) {
	r = indx[i];
	dbl_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	if (--lr_inf[r].delay == 0) {
	    dbl_btranl3_process (f, r, x);
	}
    }
    dbl_EGlpNumClearVar (v);
}

static void dbl_btranl3_process2 (dbl_factor_work * f,
      int r,
      dbl_svector * x)
{
    dbl_lr_info *lr_inf = f->lr_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    int last;
    double v;
    dbl_EGlpNumInitVar (v);

    do {
	dbl_EGlpNumCopy (v, work[r]);
	dbl_EGlpNumZero (work[r]);
	if (dbl_EGlpNumIsNeqZero (v, f->szero_tol))
	    /* if (v > szero_tol || v < -szero_tol) */
	{
	    x->indx[x->nzcnt] = r;
	    dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	r = lr_inf[r].rrank;
	nzcnt = lr_inf[r].nzcnt;
	indx = f->lrindx + lr_inf[r].rbeg;
	coef = f->lrcoef + lr_inf[r].rbeg;
	last = -1;
	for (i = 0; i < nzcnt; i++) {
	    r = indx[i];
	    dbl_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	    if (--lr_inf[r].delay == 0) {
		if (last >= 0) {
		    dbl_btranl3_process2 (f, last, x);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_btranl3 (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * x)
{
    double *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    double *acoef = a->coef;
    dbl_lr_info *lr_inf = f->lr_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (lr_inf[aindx[i]].delay++ == 0) {
	    dbl_btranl3_delay2 (f, aindx[i]);
	}
	dbl_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--lr_inf[aindx[i]].delay == 0) {
	    dbl_btranl3_process2 (f, aindx[i], x);
	}
    }
#ifdef dbl_SOLVE_DEBUG
    printf ("dbl_ILLfactor_btran x after l3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* dbl_SOLVE_DEBUG */
}

static void dbl_ILLfactor_btrane (dbl_factor_work * f,
      double *x)
{
    int *erindx = f->erindx;
    double *ercoef = f->ercoef;
    dbl_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = etacnt - 1; i >= 0; i--) {
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 > 1)
	printf ("dbl_ILLfactor_btran x before eta %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", dbl_EGlpNumToLf (x[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
	dbl_EGlpNumCopy (v, x[er_inf[i].r]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    nzcnt = er_inf[i].nzcnt;
	    beg = er_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		dbl_EGlpNumSubInnProdTo (x[erindx[beg + j]], v, ercoef[beg + j]);
	    }
	}
    }
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 <= 1)
    printf ("dbl_ILLfactor_btran x after eta:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", dbl_EGlpNumToLf (x[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_btrane2 (dbl_factor_work * f,
      dbl_svector * x)
{
    int *erindx = f->erindx;
    double *ercoef = f->ercoef;
    dbl_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int xnzcnt = x->nzcnt;
    int *xindx = x->indx;
    double *xcoef = x->coef;
    double *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = 0; i < xnzcnt; i++) {
	dbl_EGlpNumCopy (work_coef[xindx[i]], xcoef[i]);
	work_indx[xindx[i]] = i + 1;
    }
    for (i = etacnt - 1; i >= 0; i--) {
#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 > 1)
	printf ("dbl_ILLfactor_btran x before eta2 %d:", i);
	for (j = 0; j < xnzcnt; j++) {
	    printf (" %.3f*%d", dbl_EGlpNumToLf (work_coef[xindx[j]]), xindx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
	dbl_EGlpNumCopy (v, work_coef[er_inf[i].r]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    nzcnt = er_inf[i].nzcnt;
	    beg = er_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		if (work_indx[erindx[beg + j]] == 0) {
		    work_indx[erindx[beg + j]] = xnzcnt;
		    xindx[xnzcnt++] = erindx[beg + j];
		}
		dbl_EGlpNumSubInnProdTo (work_coef[erindx[beg + j]], v, ercoef[beg + j]);
	    }
	}
    }

    j = 0;
    while (j < xnzcnt) {
	dbl_EGlpNumCopy (xcoef[j], work_coef[xindx[j]]);
	dbl_EGlpNumZero (work_coef[xindx[j]]);
	work_indx[xindx[j]] = 0;
	if (dbl_EGlpNumIsEqqual (xcoef[j], dbl_zeroLpNum)) {
	    --xnzcnt;
	    xindx[j] = xindx[xnzcnt];
	} else {
	    j++;
	}
    }
    x->nzcnt = xnzcnt;

#ifdef dbl_SOLVE_DEBUG
#if (dbl_SOLVE_DEBUG+0 <= 1)
    printf ("dbl_ILLfactor_btran x after eta2:");
    for (j = 0; j < xnzcnt; j++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (xcoef[j]), xindx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_btranu (dbl_factor_work * f,
      double *a,
      dbl_svector * x)
{
    int *urindx = f->urindx;
    double *urcoef = f->urcoef;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int *cperm = f->cperm;
    int dim = f->dim;
    int xnzcnt = 0;
    int *xindx = x->indx;
    double *xcoef = x->coef;
    int nzcnt;
    int beg;
    int i;
    int j;
    double v;
    dbl_EGlpNumInitVar (v);

    for (i = 0; i < dim; i++) {
	dbl_EGlpNumCopy (v, a[cperm[i]]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    j = rperm[i];
	    beg = ur_inf[j].rbeg;
	    dbl_EGlpNumDivTo (v, urcoef[beg]);
	    if (dbl_EGlpNumIsNeqZero (v, f->szero_tol)) {	/*
																															 * if (v > szero_tol || v < -szero_tol) */
		xindx[xnzcnt] = j;
		dbl_EGlpNumCopy (xcoef[xnzcnt], v);
		xnzcnt++;
	    }
	    nzcnt = ur_inf[j].nzcnt;
	    for (j = 1; j < nzcnt; j++) {
		dbl_EGlpNumSubInnProdTo (a[urindx[beg + j]], v, urcoef[beg + j]);
	    }
	    dbl_EGlpNumZero (a[cperm[i]]);
	}
    }
    x->nzcnt = xnzcnt;
#ifdef dbl_SOLVE_DEBUG
    printf ("dbl_ILLfactor_btran x after u:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* dbl_SOLVE_DEBUG */
    dbl_EGlpNumClearVar (v);
}


static void dbl_btranu3_delay (dbl_factor_work * f,
      int r)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int nzcnt;
    int *indx;
    int i;

    r = f->rperm[f->crank[r]];
    nzcnt = ur_inf[r].nzcnt;
    indx = f->urindx + ur_inf[r].rbeg;
    for (i = 1; i < nzcnt; i++) {
	r = indx[i];
	if (ur_inf[r].delay++ == 0) {
	    dbl_btranu3_delay (f, r);
	}
    }
}

static void dbl_btranu3_delay2 (dbl_factor_work * f,
      int r)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    r = indx[i];
	    if (ur_inf[r].delay++ == 0) {
		if (last >= 0) {
		    dbl_btranu3_delay2 (f, last);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
}

static void dbl_btranu3_process (dbl_factor_work * f,
      int r,
      dbl_svector * x)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    double v;
    dbl_EGlpNumInitVar (v);

    dbl_EGlpNumCopy (v, work[r]);
    dbl_EGlpNumZero (work[r]);
    r = f->rperm[f->crank[r]];
    nzcnt = ur_inf[r].nzcnt;
    indx = f->urindx + ur_inf[r].rbeg;
    coef = f->urcoef + ur_inf[r].rbeg;
    dbl_EGlpNumDivTo (v, coef[0]);
    if (dbl_EGlpNumIsNeqqZero (v)) {
	x->indx[x->nzcnt] = r;
	dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    for (i = 1; i < nzcnt; i++) {
	r = indx[i];
	dbl_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	if (--ur_inf[r].delay == 0) {
	    dbl_btranu3_process (f, r, x);
	}
    }
    dbl_EGlpNumClearVar (v);
}

static void dbl_btranu3_process2 (dbl_factor_work * f,
      int r,
      dbl_svector * x)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    int last;
    double v;
    dbl_EGlpNumInitVar (v);

    do {
	dbl_EGlpNumCopy (v, work[r]);
	dbl_EGlpNumZero (work[r]);
	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	coef = f->urcoef + ur_inf[r].rbeg;
	dbl_EGlpNumDivTo (v, coef[0]);
	if (dbl_EGlpNumIsNeqqZero (v)) {
	    x->indx[x->nzcnt] = r;
	    dbl_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    r = indx[i];
	    dbl_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	    if (--ur_inf[r].delay == 0) {
		if (last >= 0) {
		    dbl_btranu3_process2 (f, last, x);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
    dbl_EGlpNumClearVar (v);
}

static void dbl_ILLfactor_btranu3 (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * x)
{
    double *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    double *acoef = a->coef;
    dbl_ur_info *ur_inf = f->ur_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (ur_inf[aindx[i]].delay++ == 0) {
	    dbl_btranu3_delay2 (f, aindx[i]);
	}
	dbl_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--ur_inf[aindx[i]].delay == 0) {
	    dbl_btranu3_process2 (f, aindx[i], x);
	}
    }
#ifdef dbl_SOLVE_DEBUG
    printf ("dbl_ILLfactor_btran x after u3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", dbl_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* dbl_SOLVE_DEBUG */
}

/* dbl_ILLfactor_btran solves x^tB=a^t (or, B^t x = a) for x */
void dbl_ILLfactor_btran (dbl_factor_work * f,
      dbl_svector * a,
      dbl_svector * x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx = a->indx;
    double *acoef = a->coef;
    double *work_coef = f->work_coef;
    int dim = f->dim;

#ifdef dbl_RECORD
    {
	fprintf (dbl_fsave, "b %d", a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (dbl_fsave, " %d %.16e", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	fprintf (dbl_fsave, "\n");
	fflush (dbl_fsave);
    }
#endif				/* dbl_RECORD */
#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_btran a:");
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %d %.3f", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
	aindx = a->indx;
	acoef = a->coef;
	work_coef = f->work_coef;
	nzcnt = a->nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	dbl_ILLfactor_btranu3 (f, a, &f->xtmp);
    } else {
	dbl_ILLfactor_btranu (f, work_coef, &f->xtmp);
    }

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	aindx = f->xtmp.indx;
	acoef = f->xtmp.coef;
	work_coef = f->work_coef;
	nzcnt = f->xtmp.nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	dbl_ILLfactor_btrane2 (f, &f->xtmp);
	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	    aindx = f->xtmp.indx;
	    acoef = f->xtmp.coef;
	    work_coef = f->work_coef;
	    nzcnt = f->xtmp.nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	dbl_ILLfactor_btrane (f, work_coef);
    }

    if (sparse) {
	dbl_ILLfactor_btranl3 (f, &f->xtmp, x);
    } else {
	dbl_ILLfactor_btranl2 (f, work_coef);
	dim = f->dim;
	nzcnt = 0;
	aindx = x->indx;
	acoef = x->coef;
	for (i = 0; i < dim; i++) {
	    if (dbl_EGlpNumIsNeqqZero (work_coef[i])) {
		if (dbl_EGlpNumIsNeqZero (work_coef[i], f->szero_tol))
		    /* if (work_coef[i] > szero_tol || work_coef[i] <
		       -szero_tol) */
		{
		    aindx[nzcnt] = i;
		    dbl_EGlpNumCopy (acoef[nzcnt], work_coef[i]);
		    nzcnt++;
		}
		dbl_EGlpNumZero (work_coef[i]);
	    }
	}
	x->nzcnt = nzcnt;
    }

#ifdef dbl_SORT_RESULTS
    dbl_sort_vector (x);
#endif

#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_btran x:");
	for (i = 0; i < x->nzcnt; i++) {
	    printf (" %d %.3f", x->indx[i], dbl_EGlpNumToLf (x->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */
    return;
}

static int dbl_expand_col (dbl_factor_work * f,
      int col)
{
    dbl_uc_info *uc_inf = f->uc_inf + col;
    int uc_freebeg = f->uc_freebeg;
    int nzcnt = uc_inf->nzcnt;
    int cbeg;
    double *uccoef;
    int *ucindx;
    int *ucrind;
    int i;
    int rval = 0;

    if (uc_freebeg + nzcnt + 1 >= f->uc_space) {
	rval = dbl_make_uc_space (f, nzcnt + 1);
	ILL_CLEANUP_IF (rval);
	uc_freebeg = f->uc_freebeg;
    }
    cbeg = uc_inf->cbeg;
    uccoef = f->uccoef;
    ucindx = f->ucindx;
    ucrind = f->ucrind;

    for (i = 0; i < nzcnt; i++) {
	dbl_EGlpNumCopy (uccoef[uc_freebeg + i], uccoef[cbeg + i]);
	ucindx[uc_freebeg + i] = ucindx[cbeg + i];
	ucrind[uc_freebeg + i] = ucrind[cbeg + i];
	ucindx[cbeg + i] = -1;
    }

    uc_inf->cbeg = uc_freebeg;
    f->uc_freebeg = uc_freebeg + nzcnt;
CLEANUP:
    ILL_RETURN (rval, "dbl_expand_col");
}

static int dbl_expand_row (dbl_factor_work * f,
      int row)
{
    dbl_ur_info *ur_inf = f->ur_inf + row;
    int ur_freebeg = f->ur_freebeg;
    int nzcnt = ur_inf->nzcnt;
    int rbeg;
    double *urcoef;
    int *urindx;
    int *urcind;
    int i;
    int rval = 0;

    if (ur_freebeg + nzcnt + 1 >= f->ur_space) {
	rval = dbl_make_ur_space (f, nzcnt + 1);
	ILL_CLEANUP_IF (rval);
	ur_freebeg = f->ur_freebeg;
    }
    rbeg = ur_inf->rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    urcind = f->urcind;

    for (i = 0; i < nzcnt; i++) {
	dbl_EGlpNumCopy (urcoef[ur_freebeg + i], urcoef[rbeg + i]);
	urindx[ur_freebeg + i] = urindx[rbeg + i];
	urcind[ur_freebeg + i] = urcind[rbeg + i];
	urindx[rbeg + i] = -1;
    }

    ur_inf->rbeg = ur_freebeg;
    f->ur_freebeg = ur_freebeg + nzcnt;
CLEANUP:
    ILL_RETURN (rval, "dbl_expand_row");
}

static int dbl_add_nonzero (dbl_factor_work * f,
      int row,
      int col,
      double val)
{
    dbl_ur_info *ur_inf = f->ur_inf + row;
    dbl_uc_info *uc_inf = f->uc_inf + col;
    int cnzcnt = uc_inf->nzcnt;
    int rnzcnt = ur_inf->nzcnt;
    int cloc = uc_inf->cbeg + cnzcnt;
    int rloc = ur_inf->rbeg + rnzcnt;
    int rval = 0;

    if (f->ucindx[cloc] != -1) {
	rval = dbl_expand_col (f, col);
	ILL_CLEANUP_IF (rval);
	cloc = uc_inf->cbeg + cnzcnt;
    }
    TESTG ((rval = (rloc < 0 || rloc > f->ur_space)), CLEANUP,
	"rloc %d outside boundaries [0:%d]", rloc, f->ur_space);
    if (f->urindx[rloc] != -1) {
	rval = dbl_expand_row (f, row);
	ILL_CLEANUP_IF (rval);
	rloc = ur_inf->rbeg + rnzcnt;
    }
    f->ucindx[cloc] = row;
    dbl_EGlpNumCopy (f->uccoef[cloc], val);
    f->ucrind[cloc] = rnzcnt;
    f->urindx[rloc] = col;
    dbl_EGlpNumCopy (f->urcoef[rloc], val);
    f->urcind[rloc] = cnzcnt;

    if (cloc == f->uc_freebeg)
	f->uc_freebeg++;
    if (rloc == f->ur_freebeg)
	f->ur_freebeg++;

    uc_inf->nzcnt = cnzcnt + 1;
    ur_inf->nzcnt = rnzcnt + 1;
CLEANUP:
    ILL_RETURN (rval, "dbl_add_nonzero");
}

static void dbl_delete_nonzero_row (dbl_factor_work * f,
      int row,
      int ind)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    double *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int *ucrind = f->ucrind;
    int rbeg = ur_inf[row].rbeg;
    int nzcnt = ur_inf[row].nzcnt - 1;
    int cbeg;

    if (ind != nzcnt) {
	dbl_EGlpNumCopy (urcoef[rbeg + ind], urcoef[rbeg + nzcnt]);
	urindx[rbeg + ind] = urindx[rbeg + nzcnt];
	urcind[rbeg + ind] = urcind[rbeg + nzcnt];
	cbeg = f->uc_inf[urindx[rbeg + nzcnt]].cbeg;
	ucrind[cbeg + urcind[rbeg + nzcnt]] = ind;
	urindx[rbeg + nzcnt] = -1;
    }
    ur_inf[row].nzcnt = nzcnt;
}

static void dbl_delete_nonzero_col (dbl_factor_work * f,
      int col,
      int ind)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    double *uccoef = f->uccoef;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int *urcind = f->urcind;
    int cbeg = uc_inf[col].cbeg;
    int nzcnt = uc_inf[col].nzcnt - 1;
    int rbeg;

    if (ind != nzcnt) {
	dbl_EGlpNumCopy (uccoef[cbeg + ind], uccoef[cbeg + nzcnt]);
	ucindx[cbeg + ind] = ucindx[cbeg + nzcnt];
	ucrind[cbeg + ind] = ucrind[cbeg + nzcnt];
	rbeg = f->ur_inf[ucindx[cbeg + nzcnt]].rbeg;
	urcind[rbeg + ucrind[cbeg + nzcnt]] = ind;
	ucindx[cbeg + nzcnt] = -1;
    }
    uc_inf[col].nzcnt = nzcnt;
}

static void dbl_delete_column (dbl_factor_work * f,
      int col)
{
    dbl_uc_info *uc_inf = f->uc_inf;
    int beg = uc_inf[col].cbeg;
    int nzcnt = uc_inf[col].nzcnt;
    int *ucindx = f->ucindx + beg;
    int *ucrind = f->ucrind + beg;
    int i;

    for (i = 0; i < nzcnt; i++) {
	dbl_delete_nonzero_row (f, ucindx[i], ucrind[i]);
	ucindx[i] = -1;
    }
    uc_inf[col].nzcnt = 0;

#ifdef dbl_TRACK_FACTOR
    f->nzcnt_cur -= nzcnt;
#endif

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_delete_column corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */
}

static void dbl_delete_row (dbl_factor_work * f,
      int row,
      dbl_svector * x)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int beg = ur_inf[row].rbeg;
    int nzcnt = ur_inf[row].nzcnt;
    int *urindx = f->urindx + beg;
    double *urcoef = f->urcoef + beg;
    int *urcind = f->urcind + beg;
    int i;

    for (i = 0; i < nzcnt; i++) {
	x->indx[i] = urindx[i];
	dbl_EGlpNumCopy (x->coef[i], urcoef[i]);
	dbl_delete_nonzero_col (f, urindx[i], urcind[i]);
	urindx[i] = -1;
    }
    x->nzcnt = nzcnt;
    ur_inf[row].nzcnt = 0;

#ifdef dbl_TRACK_FACTOR
    f->nzcnt_cur -= nzcnt;
#endif

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_delete_row corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */
}

static int dbl_create_column (dbl_factor_work * f,
      dbl_svector * a,
      int col,
      int *p_last_rank)
{
    int *rrank = f->rrank;
    int nzcnt = a->nzcnt;
    int *aindx = a->indx;
    double *acoef = a->coef;
    int i;
    int j;
    int rval = 0;
    int last_rank = -1;
#ifdef dbl_TRACK_FACTOR
    double max;
    dbl_EGlpNumInitVar (max);
    dbl_EGlpNumCopy (max, f->maxelem_cur);
#endif				/* dbl_TRACK_FACTOR */

    last_rank = 0;

    for (i = 0; i < nzcnt; i++) {
	rval = dbl_add_nonzero (f, aindx[i], col, acoef[i]);
	ILL_CLEANUP_IF (rval);
#ifdef dbl_TRACK_FACTOR
	dbl_EGlpNumSetToMaxAbs (max, acoef[i]);
#endif				/* dbl_TRACK_FACTOR */
	j = rrank[aindx[i]];
	if (j > last_rank)
	    last_rank = j;
    }
    *p_last_rank = last_rank;

#ifdef dbl_TRACK_FACTOR
    f->nzcnt_cur += nzcnt;
    dbl_EGlpNumCopy (f->maxelem_cur, max);
    dbl_EGlpNumClearVar (max);
#endif				/* dbl_TRACK_FACTOR */

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_create_column corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

CLEANUP:
    ILL_RETURN (rval, "dbl_create_column");
}

#ifdef dbl_UPDATE_STUDY
static int dbl_column_rank (dbl_factor_work * f,
      int col)
{
    int *cperm = f->cperm;
    int dim = f->dim;
    int i;

    for (i = 0; i < dim; i++) {
	if (cperm[i] == col) {
	    return i;
	}
    }
    return 0;
}
#endif

static void dbl_shift_permutations (dbl_factor_work * f,
      int rank_p,
      int rank_r)
{
    int *cperm = f->cperm;
    int *crank = f->crank;
    int *rperm = f->rperm;
    int *rrank = f->rrank;
    int col_p = cperm[rank_p];
    int row_p = rperm[rank_p];
    int i;

    for (i = rank_p; i < rank_r; i++) {
	cperm[i] = cperm[i + 1];
	crank[cperm[i]] = i;
	rperm[i] = rperm[i + 1];
	rrank[rperm[i]] = i;
    }
    cperm[rank_r] = col_p;
    crank[col_p] = rank_r;
    rperm[rank_r] = row_p;
    rrank[row_p] = rank_r;
}

static int dbl_eliminate_row (dbl_factor_work * f,
      int rank_p,
      int rank_r)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int *cperm = f->cperm;
    int *urindx = f->urindx;
    double *urcoef = f->urcoef;
    int *erindx = f->erindx;
    double *ercoef = f->ercoef;
    double *work_coef = f->work_coef;
    int er_freebeg = f->er_freebeg;
    int er_space = f->er_space;
    int beg;
    int nzcnt;
    int i;
    int j;
    int c;
    int r;
    double pivot_mul;
#ifdef dbl_TRACK_FACTOR
    double max;
    dbl_EGlpNumInitVar (max);
    dbl_EGlpNumCopy (max, f->maxelem_cur);
#endif				/* dbl_TRACK_FACTOR */
    dbl_EGlpNumInitVar (pivot_mul);

    for (i = rank_p; i < rank_r; i++) {
	c = cperm[i];
	if (dbl_EGlpNumIsNeqZero (work_coef[c], f->fzero_tol)) {	/*
																																					 * if (work_coef[c] > fzero_tol || work_coef[c] < -fzero_tol) */
	    r = rperm[i];
	    beg = ur_inf[r].rbeg;
	    nzcnt = ur_inf[r].nzcnt;
	    dbl_EGlpNumCopyFrac (pivot_mul, work_coef[c], urcoef[beg]);
	    dbl_EGlpNumZero (work_coef[c]);
	    for (j = 1; j < nzcnt; j++) {
		dbl_EGlpNumSubInnProdTo (work_coef[urindx[beg + j]], pivot_mul, urcoef[beg + j]);	/* 0.85 */
	    }
	    if (er_freebeg >= er_space) {
		/* fprintf (stderr, "no space in dbl_eliminate_row\n"); */
#ifdef dbl_TRACK_FACTOR
		dbl_EGlpNumClearVar (max);
#endif
		dbl_EGlpNumClearVar (pivot_mul);
		return E_UPDATE_NOSPACE;
	    }
	    erindx[er_freebeg] = r;
	    dbl_EGlpNumCopy (ercoef[er_freebeg], pivot_mul);
#ifdef dbl_TRACK_FACTOR
	    dbl_EGlpNumSetToMaxAbs (max, pivot_mul);
#endif				/* dbl_TRACK_FACTOR */
	    er_freebeg++;
	} else {
	    dbl_EGlpNumZero (work_coef[c]);
	}
    }
    f->er_freebeg = er_freebeg;
#ifdef dbl_TRACK_FACTOR
    dbl_EGlpNumCopy (f->maxelem_cur, max);
    dbl_EGlpNumClearVar (max);
#endif				/* dbl_TRACK_FACTOR */
    dbl_EGlpNumClearVar (pivot_mul);
    return 0;
}

static int dbl_create_row (dbl_factor_work * f,
      double *a,
      int row,
      int minrank)
{
    int *cperm = f->cperm;
    int dim = f->dim;
    int i;
    int j;
    int rval = 0;
#ifdef dbl_TRACK_FACTOR
    double max;
    dbl_EGlpNumInitVar (max);
    dbl_EGlpNumCopy (max, f->maxelem_cur);
#endif				/* dbl_TRACK_FACTOR */

    for (i = minrank; i < dim; i++) {
	if (dbl_EGlpNumIsNeqqZero (a[cperm[i]])) {
	    j = cperm[i];
	    if (dbl_EGlpNumIsNeqZero (a[j], f->fzero_tol)) {	/*
																																	 * if (a[j] > fzero_tol || a[j] < -fzero_tol) */
		rval = dbl_add_nonzero (f, row, j, a[j]);
		ILL_CLEANUP_IF (rval);
#ifdef dbl_TRACK_FACTOR
		dbl_EGlpNumSetToMaxAbs (max, a[j]);
#endif				/* dbl_TRACK_FACTOR */
	    }
	    dbl_EGlpNumZero (a[j]);
	}
    }

#ifdef dbl_TRACK_FACTOR
    f->nzcnt_cur += f->ur_inf[row].nzcnt;
    dbl_EGlpNumCopy (f->maxelem_cur, max);
    dbl_EGlpNumClearVar (max);
#endif				/* dbl_TRACK_FACTOR */

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_create_row corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */
CLEANUP:
    ILL_RETURN (rval, "dbl_create_row");
}

static void dbl_serow_delay (dbl_factor_work * f,
      int r,
      int rank_r)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    int *crank = f->crank;
    int nzcnt;
    int *indx;
    int i;
    int last;

    do {
	r = f->rperm[crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    r = indx[i];
	    if (ur_inf[r].delay++ == 0 && crank[r] < rank_r) {
		if (last >= 0) {
		    dbl_serow_delay (f, last, rank_r);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
}

static int dbl_serow_process (dbl_factor_work * f,
      int r,
      dbl_svector * newr,
      int rank_r)
{
    dbl_ur_info *ur_inf = f->ur_inf;
    double *work = f->work_coef;
    int nzcnt;
    int *indx;
    double *coef;
    int i;
    double v;
    int last;
    int rval;
    dbl_EGlpNumInitVar (v);

    do {
	dbl_EGlpNumCopy (v, work[r]);
	dbl_EGlpNumZero (work[r]);
	if (f->crank[r] >= rank_r) {
	    if (dbl_EGlpNumIsNeqZero (v, f->fzero_tol)) {	/*
																															 * if (v > fzero_tol || v < -fzero_tol) */
		/* stash this nonzero in the resulting row */
#ifdef dbl_TRACK_FACTOR
		dbl_EGlpNumSetToMaxAbs (f->maxelem_cur, v);
#endif				/* dbl_TRACK_FACTOR */
		newr->indx[newr->nzcnt] = r;
		dbl_EGlpNumCopy (newr->coef[newr->nzcnt], v);
		newr->nzcnt++;
		dbl_EGlpNumClearVar (v);
		return 0;
	    } else {
		dbl_EGlpNumClearVar (v);
		return 0;
	    }
	}
	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	coef = f->urcoef + ur_inf[r].rbeg;
	dbl_EGlpNumDivTo (v, coef[0]);
	if (dbl_EGlpNumIsNeqZero (v, f->fzero_tol)) {	/*
																													 * if (v > fzero_tol || v < -fzero_tol) */
	    /* stash v in eta */
	    if (f->er_freebeg >= f->er_space) {
		/* fprintf (stderr, "no space in dbl_eliminate_row\n"); */
		dbl_EGlpNumClearVar (v);
		return E_UPDATE_NOSPACE;
	    }
	    f->erindx[f->er_freebeg] = r;
	    dbl_EGlpNumCopy (f->ercoef[f->er_freebeg], v);
#ifdef dbl_TRACK_FACTOR
	    dbl_EGlpNumSetToMaxAbs (f->maxelem_cur, v);
#endif				/* dbl_TRACK_FACTOR */
	    f->er_freebeg++;
	}
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    r = indx[i];
	    dbl_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	    if (--ur_inf[r].delay == 0) {
		if (last >= 0) {
		    rval = dbl_serow_process (f, last, newr, rank_r);
		    if (rval) {
			dbl_EGlpNumClearVar (v);
			return rval;
		    }
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
    dbl_EGlpNumClearVar (v);
    return 0;
}

static int dbl_sparse_eliminate_row (dbl_factor_work * f,
      dbl_svector * x,
      int row_p,
      int rank_r)
{
    double *work = f->work_coef;
    int xnzcnt = x->nzcnt;
    int *xindx = x->indx;
    double *xcoef = x->coef;
    dbl_ur_info *ur_inf = f->ur_inf;
    int *crank = f->crank;
    int i;
    int j;
    int rval = 0;
    dbl_svector newr;

    newr.indx = 0;
    newr.coef = 0;

    for (i = 0; i < xnzcnt; i++) {
	j = xindx[i];
	if (ur_inf[j].delay++ == 0 && crank[j] < rank_r) {
	    dbl_serow_delay (f, j, rank_r);
	}
	dbl_EGlpNumCopy (work[j], xcoef[i]);
    }

    newr.nzcnt = 0;
    ILL_SAFE_MALLOC (newr.indx, f->dim, int);
    newr.coef = dbl_EGlpNumAllocArray (f->dim);

    for (i = 0; i < xnzcnt; i++) {
	j = xindx[i];
	if (--ur_inf[j].delay == 0) {
	    rval = dbl_serow_process (f, j, &newr, rank_r);
	    ILL_CLEANUP_IF (rval);
	}
    }

    for (i = 0; i < newr.nzcnt; i++) {
	rval = dbl_add_nonzero (f, row_p, newr.indx[i], newr.coef[i]);
	ILL_CLEANUP_IF (rval);
    }

#ifdef dbl_TRACK_FACTOR
    f->nzcnt_cur += newr.nzcnt;
#endif				/* dbl_TRACK_FACTOR */

CLEANUP:
    dbl_EGlpNumFreeArray (newr.coef);
    ILL_IFFREE (newr.indx, int);

    /* Bico 031210 - chg from ILL_RETURN */
    ILL_RESULT (rval, "dbl_sparse_eliminate_row");
}

static int dbl_move_pivot_row (dbl_factor_work * f,
      int r,
      int c)
{
    dbl_ur_info *ur_inf = f->ur_inf + r;
    dbl_uc_info *uc_inf = f->uc_inf;
    int beg = ur_inf->rbeg;
    int nzcnt = ur_inf->nzcnt;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int *ucrind = f->ucrind;
    double *urcoef = f->urcoef;
    double dt;
    int it;
    int i;

    if (urindx[beg] == c)
	return 0;
    dbl_EGlpNumInitVar (dt);

    for (i = 1; i < nzcnt; i++) {
	if (urindx[beg + i] == c) {
	    dbl_EGLPNUM_SWAP (urcoef[beg], urcoef[beg + i], dt);
	    dbl_ILL_SWAP (urcind[beg], urcind[beg + i], it);
	    urindx[beg + i] = urindx[beg];
	    urindx[beg] = c;
	    ucrind[uc_inf[c].cbeg + urcind[beg]] = 0;
	    ucrind[uc_inf[urindx[beg + i]].cbeg + urcind[beg + i]] = i;
	    dbl_EGlpNumClearVar (dt);
	    return 0;
	}
    }
    fprintf (stderr, "pivot row nonzero not found\n");
    dbl_EGlpNumClearVar (dt);
    return E_UPDATE_SINGULAR_ROW;
}

static int dbl_move_pivot_col (dbl_factor_work * f,
      int c,
      int r)
{
    dbl_uc_info *uc_inf = f->uc_inf + c;
    dbl_ur_info *ur_inf = f->ur_inf;
    int beg = uc_inf->cbeg;
    int nzcnt = uc_inf->nzcnt;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int *urcind = f->urcind;
    double *uccoef = f->uccoef;
    double dt;
    int i, it;

    if (ucindx[beg] == r)
	return 0;
    dbl_EGlpNumInitVar (dt);

    for (i = 1; i < nzcnt; i++) {
	if (ucindx[beg + i] == r) {
	    dbl_EGLPNUM_SWAP (uccoef[beg], uccoef[beg + i], dt);
	    dbl_ILL_SWAP (ucrind[beg], ucrind[beg + i], it);
	    ucindx[beg + i] = ucindx[beg];
	    ucindx[beg] = r;
	    urcind[ur_inf[r].rbeg + ucrind[beg]] = 0;
	    urcind[ur_inf[ucindx[beg + i]].rbeg + ucrind[beg + i]] = i;
	    dbl_EGlpNumClearVar (dt);
	    return 0;
	}
    }
    fprintf (stderr, "pivot col nonzero not found\n");
    dbl_EGlpNumClearVar (dt);
    return E_UPDATE_SINGULAR_COL;
}

static int dbl_move_pivot (dbl_factor_work * f,
      int rank_r)
{
    int r = f->rperm[rank_r];
    int c = f->cperm[rank_r];
    int rval = 0;

    rval = dbl_move_pivot_row (f, r, c);
    ILL_CLEANUP_IF (rval);

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_move_pivot_row corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

    rval = dbl_move_pivot_col (f, c, r);
    ILL_CLEANUP_IF (rval);

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_move_pivot_col corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

CLEANUP:
    ILL_RESULT (rval, "dbl_move_pivot");	/* Bico 031209 - chg from
						   RETURN */
}

int dbl_ILLfactor_update (dbl_factor_work * f,
      dbl_svector * a,
      int col_p,
      int *p_refact)
{
    int row_p;
    int rank_r = 0;
    int rank_p = 0;
    int rval = 0;
    int nzcnt;
    int *aindx;
    double *acoef;
    double *work_coef = f->work_coef;
#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
    double tmpsize;
#endif
#endif
    int i;

#ifdef dbl_RECORD
    {
	fprintf (dbl_fsave, "u %d %d", col_p, a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (dbl_fsave, " %d %.16e", a->indx[i], dbl_EGlpNumToLf (a->coef[i]));
	}
	fprintf (dbl_fsave, "\n");
	fflush (dbl_fsave);
    }
#endif				/* dbl_RECORD */

#ifdef dbl_DEBUG_FACTOR
    {
	printf ("dbl_ILLfactor_update col %d:", col_p);
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %.3f*%d", dbl_EGlpNumToLf (a->coef[i]), a->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

#ifdef dbl_DEBUG_FACTOR
    if (dbl_check_matrix (f)) {
	printf ("dbl_ILLfactor_update received corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* dbl_DEBUG_FACTOR */

    if (f->etacnt >= f->etamax) {
	*p_refact = 1;
	return 0;
    }
#ifdef dbl_UPDATE_STUDY
    dbl_nupdate++;
#endif

    row_p = f->ucindx[f->uc_inf[col_p].cbeg];

    dbl_delete_column (f, col_p);

    rval = dbl_create_column (f, a, col_p, &rank_r);
    /* if (rval) fprintf (stderr, "dbl_create_column failed\n"); */
    ILL_CLEANUP_IF (rval);

    rank_p = f->crank[col_p];
#ifdef dbl_UPDATE_STUDY
    if (rank_p != f->rrank[row_p] || rank_p != dbl_column_rank (f, col_p)) {
	printf ("rank_p %d rrank[row_p] %d dbl_column_rank(f,col_p) %d\n",
	    rank_p, f->rrank[row_p], dbl_column_rank (f, col_p));
    }
    if (rank_r > rank_p) {
	dbl_permshifttot += rank_r - rank_p;
    }
    for (i = 0; i < a->nzcnt; i++) {
	if (f->rrank[a->indx[i]] > rank_p)
	    dbl_colspiketot++;
    }
    for (i = 0; i < f->ur_inf[row_p].nzcnt; i++) {
	if (f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] <= rank_r &&
	    f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] != rank_p) {
	    dbl_rowspiketot++;
	}
    }
#endif

    dbl_shift_permutations (f, rank_p, rank_r);

    dbl_delete_row (f, row_p, &f->xtmp);

    f->er_inf[f->etacnt].rbeg = f->er_freebeg;
    f->er_inf[f->etacnt].r = row_p;

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	nzcnt = f->xtmp.nzcnt;
	aindx = f->xtmp.indx;
	acoef = f->xtmp.coef;

	for (i = 0; i < nzcnt; i++) {
	    dbl_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}

	rval = dbl_eliminate_row (f, rank_p, rank_r);
	/* if (rval) fprintf (stderr, "dbl_eliminate_row failed\n"); */
	ILL_CLEANUP_IF (rval);

	rval = dbl_create_row (f, f->work_coef, row_p, rank_r);
	/* if (rval) fprintf (stderr, "dbl_create_row failed\n"); */
	ILL_CLEANUP_IF (rval);
    } else {
	rval = dbl_sparse_eliminate_row (f, &f->xtmp, row_p, rank_r);
	/* if (rval) fprintf (stderr, "dbl_sparse_eliminate_row failed\n"); */
	ILL_CLEANUP_IF (rval);
    }

    if (f->er_freebeg - f->er_inf[f->etacnt].rbeg > 0) {
	f->er_inf[f->etacnt].nzcnt = f->er_freebeg - f->er_inf[f->etacnt].rbeg;
#ifdef dbl_TRACK_FACTOR
	f->nzcnt_cur += f->er_inf[f->etacnt].nzcnt;
#endif				/* dbl_TRACK_FACTOR */
#ifdef dbl_UPDATE_STUDY
	dbl_leftetatot += f->er_inf[f->etacnt].nzcnt;
#endif

#ifdef dbl_SORT_RESULTS
	dbl_sort_vector2 (f->er_inf[f->etacnt].nzcnt,
	    f->erindx + f->er_inf[f->etacnt].rbeg,
	    f->ercoef + f->er_inf[f->etacnt].rbeg);
#endif

	f->etacnt++;
    }
    rval = dbl_move_pivot (f, rank_r);
    /* if (rval) fprintf (stderr, "dbl_move_pivot failed\n"); */
    ILL_CLEANUP_IF (rval);

#ifdef dbl_UPDATE_DEBUG
    printf ("Updated factorization:\n");
#if (dbl_UPDATE_DEBUG+0>1)
    dbl_dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif				/* dbl_UPDATE_DEBUG */

#ifdef dbl_TRACK_FACTOR
#ifdef dbl_NOTICE_BLOWUP
    dbl_EGlpNumInitVar (tmpsize);
    dbl_EGlpNumSet (tmpsize, f->updmaxmult);
    dbl_EGlpNumMultTo (tmpsize, f->maxelem_orig);
    if (dbl_EGlpNumIsLess (tmpsize, f->maxelem_cur)) {
	/* Bico - comment out for dist fprintf (stderr, "factor_update blowup
	   max cur %e max orig %e\n", f->maxelem_cur, f->maxelem_orig); */
	dbl_EGlpNumClearVar (tmpsize);
	return E_FACTOR_BLOWUP;
    }
    dbl_EGlpNumClearVar (tmpsize);
#endif				/* dbl_NOTICE_BLOWUP */
#endif
#ifdef dbl_UPDATE_STATS
    dbl_dump_factor_stats (f);
#endif
CLEANUP:
    ILL_RESULT (rval, "dbl_ILLfactor_update");	/* Bico 031208 - chg from
						   RETURN */
}
