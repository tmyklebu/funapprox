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

/* RCS_INFO = "$RCSfile: mpq_factor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
static int TRACE = 0;

/* implement a = max(a,abs(b)) and execute the extra code if the update is
   needed */
#define mpq_EGlpNumSetToMaxAbsAndDo(a,b,c) \
	if(mpq_EGlpNumIsLess(mpq_zeroLpNum,b))\
	{\
		if(mpq_EGlpNumIsLess(a,b)){\
			mpq_EGlpNumCopy(a,b);\
			c;\
			}\
	}\
	else\
	{\
		mpq_EGlpNumSign(a);\
		if(mpq_EGlpNumIsLess(b,a)){\
			mpq_EGlpNumCopy(a,b);\
			c;\
			}\
		mpq_EGlpNumSign(a);\
	}


#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdefs.h"
#include "mpq_factor.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


#undef  mpq_RECORD
#undef  mpq_DEBUG_FACTOR
#undef  mpq_SOLVE_DEBUG

#undef  mpq_FACTOR_DEBUG
#undef  mpq_UPDATE_DEBUG

#undef mpq_TRACK_FACTOR
#undef mpq_NOTICE_BLOWUP

#undef  mpq_FACTOR_STATS
#undef  mpq_UPDATE_STATS
#undef  mpq_GROWTH_STATS

#undef  mpq_UPDATE_STUDY

#undef  mpq_SORT_RESULTS

#ifdef mpq_UPDATE_STUDY
int mpq_nupdate = 0;
long int mpq_colspiketot = 0.0;
long int mpq_rowspiketot = 0.0;
long int mpq_permshifttot = 0.0;
long int mpq_leftetatot = 0.0;
#endif

void mpq_ILLfactor_init_factor_work (mpq_factor_work * f)
{
    f->max_k = 1000;		/* must be less than 46340 (2^15.5) */
    mpq_EGlpNumCopy (f->fzero_tol, mpq_SZERO_TOLER);	/* 2^-50 */
    mpq_EGlpNumCopy (f->szero_tol, mpq_SZERO_TOLER);	/* 2^-50 */
    mpq_EGlpNumCopy (f->partial_tol, mpq_OBJBND_TOLER);	/* 2^-7 */
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
    mpq_EGlpNumCopy (f->partial_cur, f->partial_tol);
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
    mpq_ILLsvector_init (&f->xtmp);
}

void mpq_ILLfactor_free_factor_work (mpq_factor_work * f)
{
#ifdef mpq_UPDATE_STUDY
    if (mpq_nupdate) {
	printf
	("UPDATE STUDY: avg %d upd: %.2f col, %.2f row, %.2f lefteta, %.2f perm\n",
	    mpq_nupdate, ((double) mpq_colspiketot) / mpq_nupdate,
	    ((double) mpq_rowspiketot) / mpq_nupdate, ((double) mpq_leftetatot) / mpq_nupdate,
	    ((double) mpq_permshifttot) / mpq_nupdate);
    }
#endif
    mpq_EGlpNumFreeArray (f->work_coef);
    ILL_IFFREE (f->work_indx, int);
    ILL_IFFREE (f->uc_inf, mpq_uc_info);
    if (f->dim + f->max_k > 0 && f->ur_inf) {
	unsigned int i = f->dim + f->max_k + 1;
	while (i--)
	    mpq_EGlpNumClearVar (f->ur_inf[i].max);
    }
    ILL_IFFREE (f->ur_inf, mpq_ur_info);
    ILL_IFFREE (f->lc_inf, mpq_lc_info);
    ILL_IFFREE (f->lr_inf, mpq_lr_info);
    ILL_IFFREE (f->er_inf, mpq_er_info);
    ILL_IFFREE (f->ucindx, int);
    ILL_IFFREE (f->ucrind, int);
    mpq_EGlpNumFreeArray (f->uccoef);
    ILL_IFFREE (f->urindx, int);
    ILL_IFFREE (f->urcind, int);
    mpq_EGlpNumFreeArray (f->urcoef);
    ILL_IFFREE (f->lcindx, int);
    mpq_EGlpNumFreeArray (f->lccoef);
    ILL_IFFREE (f->lrindx, int);
    mpq_EGlpNumFreeArray (f->lrcoef);
    ILL_IFFREE (f->erindx, int);
    mpq_EGlpNumFreeArray (f->ercoef);
    ILL_IFFREE (f->rperm, int);
    ILL_IFFREE (f->rrank, int);
    ILL_IFFREE (f->cperm, int);
    ILL_IFFREE (f->crank, int);
    mpq_EGlpNumFreeArray (f->dmat);
    mpq_ILLsvector_free (&f->xtmp);
}

int mpq_ILLfactor_set_factor_iparam (mpq_factor_work * f,
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
	fprintf (stderr, "Invalid param %d in mpq_ILLfactor_set_factor_iparam\n",
	    param);
	return 1;
    }
    return 0;
}

int mpq_ILLfactor_set_factor_dparam (mpq_factor_work * f,
      int param,
      mpq_t val)
{
    switch (param) {
    case QS_FACTOR_FZERO_TOL:
	mpq_EGlpNumCopy (f->fzero_tol, val);
	break;
    case QS_FACTOR_SZERO_TOL:
	mpq_EGlpNumCopy (f->szero_tol, val);
	break;
    case QS_FACTOR_UR_SPACE_MUL:
	f->ur_space_mul = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_UC_SPACE_MUL:
	f->uc_space_mul = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_LC_SPACE_MUL:
	f->lc_space_mul = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_LR_SPACE_MUL:
	f->lr_space_mul = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_ER_SPACE_MUL:
	f->er_space_mul = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_GROW_MUL:
	f->grow_mul = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_MAXMULT:
	f->maxmult = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_UPDMAXMULT:
	f->updmaxmult = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_DENSE_FRACT:
	f->dense_fract = mpq_EGlpNumToLf (val);
	break;
    case QS_FACTOR_PARTIAL_TOL:
	mpq_EGlpNumCopy (f->partial_tol, val);
	mpq_EGlpNumCopy (f->partial_cur, val);
	break;
    default:
	fprintf (stderr, "Invalid param %d in mpq_ILLfactor_set_factor_dparam\n",
	    param);
	return 1;
    }
    return 0;
}

int mpq_ILLfactor_create_factor_work (mpq_factor_work * f,
      int dim)
{
    int i;
    int rval;

    f->dim = dim;
    f->etacnt = 0;
    f->work_coef = mpq_EGlpNumAllocArray (dim);
    ILL_SAFE_MALLOC (f->work_indx, dim, int);
    ILL_SAFE_MALLOC (f->uc_inf, dim + (f->max_k + 1), mpq_uc_info);
    ILL_SAFE_MALLOC (f->ur_inf, dim + (f->max_k + 1), mpq_ur_info);
    ILL_SAFE_MALLOC (f->lc_inf, dim, mpq_lc_info);
    ILL_SAFE_MALLOC (f->lr_inf, dim, mpq_lr_info);
    ILL_SAFE_MALLOC (f->rperm, dim, int);
    ILL_SAFE_MALLOC (f->rrank, dim, int);
    ILL_SAFE_MALLOC (f->cperm, dim, int);
    ILL_SAFE_MALLOC (f->crank, dim, int);

    for (i = dim + f->max_k + 1; i--;)
	mpq_EGlpNumInitVar (f->ur_inf[i].max);

    for (i = 0; i < dim; i++) {
	mpq_EGlpNumZero (f->work_coef[i]);
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

    rval = mpq_ILLsvector_alloc (&f->xtmp, dim);
    ILL_CLEANUP_IF (rval);

    rval = 0;

CLEANUP:
    if (rval) {
	mpq_ILLfactor_free_factor_work (f);
    }
    ILL_RETURN (rval, "mpq_ILLfactor_create_factor_work");
}

#ifdef mpq_UPDATE_DEBUG
static void mpq_dump_matrix (mpq_factor_work * f,
      int remaining)
{
    int dim = f->dim;
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_lc_info *lc_inf = f->lc_inf;
    mpq_lr_info *lr_inf = f->lr_inf;
    mpq_er_info *er_inf = f->er_inf;
    int nzcnt;
    int beg;

    int i;
    int j;

    for (i = 0; i < dim; i++) {
	if (!remaining || ur_inf[i].next >= 0) {
	    printf ("Row %d %d (max %.3f):", i, f->rrank[i],
		mpq_EGlpNumToLf (ur_inf[i].max));
	    nzcnt = ur_inf[i].nzcnt;
	    beg = ur_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		if (j == ur_inf[i].pivcnt) {
		    printf (" |");
		}
		printf (" %.3f*%d", mpq_EGlpNumToLf (f->urcoef[beg + j]),
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
		mpq_EGlpNumToLf (ur_inf[f->rperm[i + f->dense_base]].max));
	    for (j = start; j < f->dcols; j++) {
		if (j == f->drows) {
		    printf (" |");
		}
		printf (" %.3f", mpq_EGlpNumToLf (f->dmat[i * f->dcols + j]));
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
		printf (" %.3f*%d", mpq_EGlpNumToLf (f->lccoef[beg + j]),
		    f->lcindx[beg + j]);
	    }
	    printf ("\n");
	}
	for (i = f->nstages; i < f->dim; i++) {
	    printf ("L col %d:", lc_inf[i].c);
	    nzcnt = lc_inf[i].nzcnt;
	    beg = lc_inf[i].cbeg;
	    for (j = 0; j < nzcnt; j++) {
		printf (" %.3f*%d", mpq_EGlpNumToLf (f->lccoef[beg + j]),
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
		printf (" %.3f*%d", mpq_EGlpNumToLf (f->lrcoef[beg + j]),
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
		printf (" %.3f*%d", mpq_EGlpNumToLf (f->ercoef[beg + j]),
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
		    printf (" %.3f*%d", mpq_EGlpNumToLf (f->uccoef[beg + j]),
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

#ifdef mpq_SORT_RESULTS
static void mpq_sort_vector2 (int nzcnt,
      int *indx,
      mpq_t * coef)
{
    int i;
    int j;
    int itmp;
    mpq_t ctmp;
    mpq_EGlpNumInitVar (ctmp);

    for (i = 1; i < nzcnt; i++) {
	itmp = indx[i];
	mpq_EGlpNumCopy (ctmp, coef[i]);
	for (j = i; j >= 1 && indx[j - 1] > itmp; j--) {
	    indx[j] = indx[j - 1];
	    mpq_EGlpNumCopy (coef[j], coef[j - 1]);
	}
	indx[j] = itmp;
	mpq_EGlpNumCopy (coef[j], ctmp);
    }
    mpq_EGlpNumClearVar (ctmp);
}

static void mpq_sort_vector (mpq_svector * x)
{
    mpq_sort_vector2 (x->nzcnt, x->indx, x->coef);
}
#endif

#ifdef mpq_DEBUG_FACTOR
static int mpq_check_matrix (mpq_factor_work * f)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_uc_info *uc_inf = f->uc_inf;
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
	mpq_dump_matrix (f, 0);
	return E_CHECK_FAILED;
    }
    return 0;
}
#endif

#ifdef mpq_FACTOR_STATS
static void mpq_dump_factor_stats (mpq_factor_work * f)
{
    int dim = f->dim;
    int ecnt = f->etacnt;
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_lc_info *lc_inf = f->lc_inf;
    mpq_er_info *er_inf = f->er_inf;
    mpq_t *urcoef = f->urcoef;
    mpq_t *lccoef = f->lccoef;
    mpq_t *ercoef = f->ercoef;
    int lnzcnt = 0;
    int unzcnt = 0;
    int enzcnt = 0;
    int nzcnt;
    int beg;
    mpq_t umax;
    mpq_t lmax;
    mpq_t emax;
    int i;
    int j;
    mpq_EGlpNumInitVar (umax);
    mpq_EGlpNumInitVar (lmax);
    mpq_EGlpNumInitVar (emax);
    mpq_EGlpNumZero (umax);
    for (i = 0; i < dim; i++) {
	nzcnt = ur_inf[i].nzcnt;
	beg = ur_inf[i].rbeg;
	unzcnt += nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumSetToMaxAbs (umax, urcoef[beg + j]);
	}
    }
    mpq_EGlpNumZero (lmax);
    for (i = 0; i < dim; i++) {
	nzcnt = lc_inf[i].nzcnt;
	beg = lc_inf[i].cbeg;
	lnzcnt += nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumSetToMaxAbs (lmax, lccoef[beg + j]);
	}
    }
    mpq_EGlpNumZero (emax);
    for (i = 0; i < ecnt; i++) {
	nzcnt = er_inf[i].nzcnt;
	beg = er_inf[i].rbeg;
	enzcnt += nzcnt;
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumSetToMaxAbs (emax, ercoef[beg + j]);
	}
    }
    printf
	("factor U %d nzs %.3e max L %d nzs %.3e max E %d nzs %.3e max\n",
	unzcnt, mpq_EGlpNumToLf (umax), lnzcnt, mpq_EGlpNumToLf (lmax), enzcnt,
	mpq_EGlpNumToLf (emax));
    fflush (stdout);
    mpq_EGlpNumClearVar (umax);
    mpq_EGlpNumClearVar (lmax);
    mpq_EGlpNumClearVar (emax);
}
#endif

static void mpq_clear_work (mpq_factor_work * f)
{
    int i;
    int dim = f->dim;
    mpq_t *work_coef = f->work_coef;

    for (i = 0; i < dim; i++) {
	mpq_EGlpNumZero (work_coef[i]);
    }
}

static void mpq_load_row (mpq_factor_work * f,
      int r)
{
    mpq_t *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    mpq_t *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;

    for (i = 0; i < prow_nzcnt; i++) {
	j = prow_urindx[i];
	mpq_EGlpNumCopy (work_coef[j], prow_urcoef[i]);
	work_indx[j] = 1;
    }
}

static void mpq_clear_row (mpq_factor_work * f,
      int r)
{
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    mpq_t *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;

    for (i = 0; i < prow_nzcnt; i++) {
	j = prow_urindx[i];
	mpq_EGlpNumZero (work_coef[j]);
	work_indx[j] = 0;
    }
}

static int mpq_make_ur_space (mpq_factor_work * f,
      int space)
{
    mpq_t *new_urcoef = 0;
    int *new_urindx = 0;
    int *new_urcind = 0;
    mpq_t *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int minspace;
    mpq_ur_info *ur_inf = f->ur_inf;
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

#ifdef mpq_GROWTH_STATS
    printf ("mpq_make_ur_space growing from %d to %d...", f->ur_space, minspace);
    fflush (stdout);
#endif
    new_urcoef = mpq_EGlpNumAllocArray (minspace);
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
		mpq_EGlpNumCopy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
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
		mpq_EGlpNumCopy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
		new_nzcnt++;
	    }
	}
    }

    for (i = new_nzcnt; i < minspace; i++) {
	new_urindx[i] = -1;
    }
    new_urindx[minspace] = 0;
    mpq_EGlpNumFreeArray (f->urcoef);
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

#ifdef mpq_GROWTH_STATS
    printf ("%d/%d nonzeros\n", new_nzcnt, old_nzcnt);
    fflush (stdout);
    mpq_dump_factor_stats (f);
#endif

    rval = 0;

CLEANUP:
    ILL_IFFREE (new_urcoef, mpq_t);
    ILL_IFFREE (new_urindx, int);
    ILL_IFFREE (new_urcind, int);
    ILL_RETURN (rval, "mpq_make_ur_space");
}

static int mpq_make_uc_space (mpq_factor_work * f,
      int space)
{
    mpq_t *new_uccoef = 0;
    int *new_ucindx = 0;
    int *new_ucrind = 0;
    int uc_freebeg = f->uc_freebeg;
    mpq_t *uccoef = f->uccoef;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int minspace = uc_freebeg + space;
    mpq_uc_info *uc_inf = f->uc_inf;
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

#ifdef mpq_GROWTH_STATS
    printf ("mpq_make_uc_space growing from %d to %d...", f->uc_space, minspace);
    fflush (stdout);
#endif

    ILL_SAFE_MALLOC (new_ucindx, minspace + 1, int);

    if (ucrind) {
	new_uccoef = mpq_EGlpNumAllocArray (minspace);
	ILL_SAFE_MALLOC (new_ucrind, minspace, int);
    }
    if (ucrind) {
	for (j = 0; j < dim; j++) {
	    cbeg = uc_inf[j].cbeg;
	    nzcnt = uc_inf[j].nzcnt;
	    uc_inf[j].cbeg = new_nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		new_ucindx[new_nzcnt] = ucindx[cbeg + i];
		mpq_EGlpNumCopy (new_uccoef[new_nzcnt], uccoef[cbeg + i]);
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

    mpq_EGlpNumFreeArray (f->uccoef);
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

#ifdef mpq_GROWTH_STATS
    printf ("%d nonzeros\n", new_nzcnt);
    fflush (stdout);
    mpq_dump_factor_stats (f);
#endif

    rval = 0;

CLEANUP:
    ILL_IFFREE (new_uccoef, mpq_t);
    ILL_IFFREE (new_ucindx, int);
    ILL_IFFREE (new_ucrind, int);
    ILL_RETURN (rval, "mpq_make_uc_space");
}

static int mpq_make_lc_space (mpq_factor_work * f,
      int space)
{
    mpq_t *new_lccoef = 0;
    int *new_lcindx = 0;
    int lc_freebeg = f->lc_freebeg;
    mpq_t *lccoef = f->lccoef;
    int *lcindx = f->lcindx;
    int minspace = lc_freebeg + space;
    int i;
    int rval;
    if (f->lc_space * f->grow_mul > minspace) {
	minspace = f->lc_space * f->grow_mul;
    }
#ifdef mpq_GROWTH_STATS
    printf ("mpq_make_lc_space growing from %d to %d...", f->lc_space, minspace);
    fflush (stdout);
#endif

    new_lccoef = mpq_EGlpNumAllocArray (minspace);
    ILL_SAFE_MALLOC (new_lcindx, minspace, int);

    for (i = 0; i < lc_freebeg; i++) {
	mpq_EGlpNumCopy (new_lccoef[i], lccoef[i]);
	new_lcindx[i] = lcindx[i];
    }

    mpq_EGlpNumFreeArray (lccoef);
    f->lccoef = new_lccoef;
    new_lccoef = 0;

    ILL_IFFREE (lcindx, int);
    f->lcindx = new_lcindx;
    new_lcindx = 0;

    f->lc_space = minspace;

#ifdef mpq_GROWTH_STATS
    printf ("done\n");
    fflush (stdout);
    mpq_dump_factor_stats (f);
#endif

    rval = 0;

CLEANUP:
    ILL_IFFREE (new_lccoef, mpq_t);
    ILL_IFFREE (new_lcindx, int);
    ILL_RETURN (rval, "mpq_make_lc_space");
}

static void mpq_set_col_nz (mpq_factor_work * f,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
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

static void mpq_set_row_nz (mpq_factor_work * f,
      int r)
{
    mpq_ur_info *ur_inf = f->ur_inf;
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

static void mpq_remove_col_nz (mpq_factor_work * f,
      int r,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
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

    mpq_set_col_nz (f, c);
}

static void mpq_remove_row_nz (mpq_factor_work * f,
      int r,
      int c)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    int *urindx = f->urindx + ur_inf[r].rbeg;
    mpq_t *urcoef = f->urcoef + ur_inf[r].rbeg;
    int pivcnt = ur_inf[r].pivcnt;
    mpq_t max;
    int tind;
    mpq_t tcoef;
    int i;
    mpq_EGlpNumInitVar (tcoef);
    mpq_EGlpNumInitVar (max);
    mpq_EGlpNumZero (max);

    for (i = 0; i < pivcnt; i++) {
	if (urindx[i] == c) {
	    --pivcnt;
	    mpq_ILL_SWAP (urindx[i], urindx[pivcnt], tind);
	    mpq_EGLPNUM_SWAP (urcoef[i], urcoef[pivcnt], tcoef);
	    --i;
	} else {
	    mpq_EGlpNumSetToMaxAbs (max, urcoef[i]);
	}
    }
    ur_inf[r].pivcnt = pivcnt;
    mpq_EGlpNumCopy (ur_inf[r].max, max);
    mpq_set_row_nz (f, r);
    mpq_EGlpNumClearVar (max);
    mpq_EGlpNumClearVar (tcoef);
}

static int mpq_add_col_nz (mpq_factor_work * f,
      int r,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
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
	    rval = mpq_make_uc_space (f, nzcnt + 1);
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

    mpq_set_col_nz (f, c);
CLEANUP:
    ILL_RETURN (rval, "mpq_add_col_nz");
}

static void mpq_disable_col (mpq_factor_work * f,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;

    if (uc_inf[c].next >= 0) {
	uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
	uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

	uc_inf[c].next = -2;
	uc_inf[c].prev = -2;
    }
}

static void mpq_remove_col (mpq_factor_work * f,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
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

static void mpq_remove_row (mpq_factor_work * f,
      int r)
{
    mpq_ur_info *ur_inf = f->ur_inf;

    if (ur_inf[r].next >= 0) {
	ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
	ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

	ur_inf[r].next = -1;
	ur_inf[r].prev = -1;
    }
}

static void mpq_find_coef (mpq_factor_work * f,
      int r,
      int c,
      mpq_t * coef)
{
    mpq_t *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
    int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
    int i;
    int prow_nzcnt = f->ur_inf[r].nzcnt;
    mpq_EGlpNumZero (*coef);
    for (i = 0; i < prow_nzcnt; i++) {
	if (prow_urindx[i] == c) {
	    mpq_EGlpNumCopy (*coef, prow_urcoef[i]);
	    return;
	}
    }
    fprintf (stderr, "Coefficient not found\n");
    return;
}

static int mpq_elim_row (mpq_factor_work * f,
      int elim_r,
      int r,
      int c,
      mpq_t * p_pivot_coef)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_t *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    mpq_t *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int prow_beg = ur_inf[r].rbeg;
    int prow_nzcnt = ur_inf[r].nzcnt;
    int prow_pivcnt = ur_inf[r].pivcnt;
    int fill = ur_inf[elim_r].nzcnt;
    int cancel = 0;
    mpq_t max;
    int erow_beg;
    int erow_nzcnt;
    int erow_pivcnt;
    mpq_t x;
    int i;
    int j;
    int rval = 0;
    mpq_t elim_coef;
    mpq_EGlpNumInitVar (max);
    mpq_EGlpNumInitVar (x);
    mpq_EGlpNumInitVar (elim_coef);
    mpq_EGlpNumZero (max);
    mpq_find_coef (f, r, c, &elim_coef);
    mpq_EGlpNumDivTo (elim_coef, work_coef[c]);
    mpq_EGlpNumCopy (*p_pivot_coef, elim_coef);

    for (i = 0; i < prow_nzcnt; i++) {
	j = urindx[prow_beg + i];
	if (work_indx[j] == 1) {
	    mpq_EGlpNumCopy (x, urcoef[prow_beg + i]);
	    mpq_EGlpNumSubInnProdTo (x, elim_coef, work_coef[j]);
	    if (mpq_EGlpNumIsEqual (x, mpq_zeroLpNum, f->fzero_tol) || j == c) {
		cancel++;
		if (j != c) {
		    mpq_remove_col_nz (f, r, j);
		}
		if (i < prow_pivcnt) {
		    prow_pivcnt--;
		    prow_nzcnt--;
		    urindx[prow_beg + i] = urindx[prow_beg + prow_pivcnt];
		    mpq_EGlpNumCopy (urcoef[prow_beg + i], urcoef[prow_beg + prow_pivcnt]);
		    if (prow_pivcnt != prow_nzcnt) {
			urindx[prow_beg + prow_pivcnt] = urindx[prow_beg + prow_nzcnt];
			mpq_EGlpNumCopy (urcoef[prow_beg + prow_pivcnt],
			    urcoef[prow_beg + prow_nzcnt]);
		    }
		} else {
		    prow_nzcnt--;
		    urindx[prow_beg + i] = urindx[prow_beg + prow_nzcnt];
		    mpq_EGlpNumCopy (urcoef[prow_beg + i], urcoef[prow_beg + prow_nzcnt]);
		}
		urindx[prow_beg + prow_nzcnt] = -1;
		i--;
	    } else {
		mpq_EGlpNumCopy (urcoef[prow_beg + i], x);
		if (i < prow_pivcnt) {
		    mpq_EGlpNumSetToMaxAbs (max, x);
		}
	    }
	    work_indx[j] = 0;
	    fill--;
	} else {
	    if (i < prow_pivcnt) {
		mpq_EGlpNumSetToMaxAbs (max, urcoef[prow_beg + i]);
	    }
	}
    }

    if (fill > 0) {
	ur_inf[r].nzcnt = prow_nzcnt;
	ur_inf[r].pivcnt = prow_pivcnt;
	if (fill > cancel) {
	    int ur_freebeg = f->ur_freebeg;

	    if (ur_freebeg + prow_nzcnt + fill >= f->ur_space) {
		rval = mpq_make_ur_space (f, prow_nzcnt + fill);
		ILL_CLEANUP_IF (rval);
		urcoef = f->urcoef;
		urindx = f->urindx;
		ur_freebeg = f->ur_freebeg;
		prow_beg = f->ur_inf[r].rbeg;
	    }
	    for (i = 0; i < prow_nzcnt; i++) {
		urindx[ur_freebeg + i] = urindx[prow_beg + i];
		mpq_EGlpNumCopy (urcoef[ur_freebeg + i], urcoef[prow_beg + i]);
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
		mpq_EGlpNumCopyNeg (x, elim_coef);
		mpq_EGlpNumMultTo (x, urcoef[erow_beg + i]);
		if (mpq_EGlpNumIsNeqZero (x, f->fzero_tol)) {
		    rval = mpq_add_col_nz (f, r, j);
		    ILL_CLEANUP_IF (rval);
		    if (prow_pivcnt != prow_nzcnt) {
			urindx[prow_beg + prow_nzcnt] = urindx[prow_beg + prow_pivcnt];
			mpq_EGlpNumCopy (urcoef[prow_beg + prow_nzcnt],
			    urcoef[prow_beg + prow_pivcnt]);
		    }
		    urindx[prow_beg + prow_pivcnt] = j;
		    mpq_EGlpNumCopy (urcoef[prow_beg + prow_pivcnt], x);
		    mpq_EGlpNumSetToMaxAbs (max, x);
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
		mpq_EGlpNumCopyNeg (x, elim_coef);
		mpq_EGlpNumMultTo (x, urcoef[erow_beg + i]);
		if (mpq_EGlpNumIsNeqZero (x, f->fzero_tol)) {
		    rval = mpq_add_col_nz (f, r, j);
		    ILL_CLEANUP_IF (rval);
		    urindx[prow_beg + prow_nzcnt] = j;
		    mpq_EGlpNumCopy (urcoef[prow_beg + prow_nzcnt], x);
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
    mpq_EGlpNumCopy (ur_inf[r].max, max);

    mpq_set_row_nz (f, r);
CLEANUP:
    mpq_EGlpNumClearVar (elim_coef);
    mpq_EGlpNumClearVar (x);
    mpq_EGlpNumClearVar (max);
    ILL_RETURN (rval, "mpq_elim_row");
}

#define mpq_SETPERM(f,s,r,c) {                    \
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

static int mpq_elim (mpq_factor_work * f,
      int r,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_lc_info *lc_inf = f->lc_inf;
    int *urindx;
    int *ucindx;
    int *lcindx;
    mpq_t *urcoef;
    mpq_t *lccoef;
    mpq_t pivot_coef;
    int nzcnt;
    int lc_freebeg;
    int s = f->stage;
    int i;
    int j;
    int rval = 0;
    mpq_EGlpNumInitVar (pivot_coef);

    if (uc_inf[c].nzcnt == 1) {
	/* col singleton */
	mpq_SETPERM (f, s, r, c);

	lc_inf[s].cbeg = -1;
	lc_inf[s].c = r;
	lc_inf[s].nzcnt = 0;
	f->stage++;

	urindx = f->urindx + ur_inf[r].rbeg;
	urcoef = f->urcoef + ur_inf[r].rbeg;
	nzcnt = ur_inf[r].nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    j = urindx[i];
	    mpq_remove_col_nz (f, r, j);
	    if (j == c) {
		urindx[i] = urindx[0];
		urindx[0] = c;
		mpq_EGLPNUM_SWAP (urcoef[0], urcoef[i], pivot_coef);
	    }
	}
	mpq_remove_row (f, r);
	mpq_remove_col (f, c);
    } else if (ur_inf[r].nzcnt == 1) {
	/* row singleton */
	--(f->nstages);
	mpq_SETPERM (f, f->nstages, r, c);

	lc_inf[f->nstages].cbeg = -1;
	lc_inf[f->nstages].c = r;
	lc_inf[f->nstages].nzcnt = 0;

	ucindx = f->ucindx + uc_inf[c].cbeg;
	nzcnt = uc_inf[c].nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    j = ucindx[i];
	    mpq_remove_row_nz (f, j, c);
	}
	mpq_remove_row (f, r);
	mpq_remove_col (f, c);
    } else {
	mpq_SETPERM (f, s, r, c);
	f->stage++;

	nzcnt = uc_inf[c].nzcnt;
	if (f->lc_freebeg + nzcnt >= f->lc_space) {
	    rval = mpq_make_lc_space (f, nzcnt);
	    ILL_CLEANUP_IF (rval);
	}
	lc_freebeg = f->lc_freebeg;
	lc_inf[s].cbeg = lc_freebeg;
	lc_inf[s].c = r;
	lcindx = f->lcindx;
	lccoef = f->lccoef;
	mpq_load_row (f, r);
	ucindx = f->ucindx + uc_inf[c].cbeg;
	for (i = 0; i < nzcnt; i++) {
	    j = f->ucindx[uc_inf[c].cbeg + i];
	    if (j != r) {
		rval = mpq_elim_row (f, r, j, c, &pivot_coef);
		ILL_CLEANUP_IF (rval);
		lcindx[lc_freebeg] = j;
		mpq_EGlpNumCopy (lccoef[lc_freebeg], pivot_coef);
		lc_freebeg++;
#ifdef mpq_TRACK_FACTOR
		mpq_EGlpNumSetToMaxAbs (f->maxelem_factor, pivot_coef);
		if (mpq_EGlpNumIsLess (f->maxelem_factor, ur_inf[r].max))
		    mpq_EGlpNumCopy (f->maxelem_factor, ur_inf[r].max);
#endif				/* mpq_TRACK_FACTOR */
	    }
	}
	lc_inf[s].nzcnt = lc_freebeg - lc_inf[s].cbeg;
	f->lc_freebeg = lc_freebeg;

	mpq_clear_row (f, r);

	urindx = f->urindx + ur_inf[r].rbeg;
	urcoef = f->urcoef + ur_inf[r].rbeg;
	nzcnt = ur_inf[r].nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    j = urindx[i];
	    mpq_remove_col_nz (f, r, j);
	    if (j == c) {
		urindx[i] = urindx[0];
		urindx[0] = c;
		mpq_EGLPNUM_SWAP (urcoef[0], urcoef[i], pivot_coef);
	    }
	}
	mpq_remove_row (f, r);
	mpq_remove_col (f, c);
    }
CLEANUP:
    mpq_EGlpNumClearVar (pivot_coef);
    ILL_RETURN (rval, "mpq_elim");
}

static void mpq_find_pivot_column (mpq_factor_work * f,
      int c,
      int *p_r)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *ucindx = f->ucindx;
    int nzcnt = uc_inf[c].nzcnt;
    int cbeg = uc_inf[c].cbeg;
    mpq_t num_tmp[2];
    int bestnz = -1;
    int i;
    int r;
    mpq_EGlpNumInitVar (num_tmp[0]);
    mpq_EGlpNumInitVar (num_tmp[1]);

    *p_r = -1;
    for (i = 0; i < nzcnt; i++) {
	r = ucindx[cbeg + i];
        if ((bestnz == -1 || ur_inf[r].pivcnt < bestnz)) {
            mpq_find_coef (f, r, c, &num_tmp[0]);
            if (mpq_EGlpNumIsLess (num_tmp[0], mpq_zeroLpNum))
                mpq_EGlpNumSign (num_tmp[0]);
            mpq_EGlpNumCopy (num_tmp[1], f->partial_cur);
            mpq_EGlpNumMultTo (num_tmp[1], ur_inf[r].max);
            if (mpq_EGlpNumIsLeq (num_tmp[1], num_tmp[0])) {
                bestnz = ur_inf[r].pivcnt;
                *p_r = r;
            }
	}
    }
    mpq_EGlpNumClearVar (num_tmp[0]);
    mpq_EGlpNumClearVar (num_tmp[1]);
}

static void mpq_find_pivot_row (mpq_factor_work * f,
      int r,
      int *p_c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *urindx = f->urindx;
    mpq_t *urcoef = f->urcoef;
    int pivcnt = ur_inf[r].pivcnt;
    int rbeg = ur_inf[r].rbeg;
    mpq_t thresh[2];
    int bestnz = -1;
    int i;
    int c;
    mpq_EGlpNumInitVar (thresh[0]);
    mpq_EGlpNumInitVar (thresh[1]);
    mpq_EGlpNumCopy (thresh[0], f->partial_cur);
    mpq_EGlpNumMultTo (thresh[0], ur_inf[r].max);
    *p_c = -1;
    for (i = 0; i < pivcnt; i++) {
        c = urindx[rbeg + i];
        if ((bestnz == -1 || uc_inf[c].nzcnt < bestnz)) {
            /* DAVID FIX */
            mpq_EGlpNumCopyAbs (thresh[1], urcoef[rbeg + i]);
            if (mpq_EGlpNumIsLeq (thresh[0], thresh[1])) {
                bestnz = uc_inf[c].nzcnt;
                *p_c = c;
            }
        }
    }
    mpq_EGlpNumClearVar (thresh[0]);
    mpq_EGlpNumClearVar (thresh[1]);
}

static int mpq_find_pivot (mpq_factor_work * f,
      int *p_r,
      int *p_c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_ur_info *ur_inf = f->ur_inf;
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
		mpq_find_pivot_column (f, c, &r);
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
		    mpq_disable_col (f, uc_inf[c].next);
		}
		n++;
		if (n >= p && mm != 0) {
		    return 0;
		}
	    }
	}
	if (ur_inf[dim + k].next != dim + k) {
	    for (r = ur_inf[dim + k].next; r != dim + k; r = ur_inf[r].next) {
		mpq_find_pivot_row (f, r, &c);
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

static int mpq_create_factor_space (mpq_factor_work * f)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_ur_info *ur_inf = f->ur_inf;
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
	mpq_EGlpNumFreeArray (f->urcoef);
	f->ur_space = nzcnt * f->ur_space_mul;
	ILL_SAFE_MALLOC (f->urindx, f->ur_space + 1, int);
	f->urcoef = mpq_EGlpNumAllocArray (f->ur_space);
    }
    if (f->lcindx == 0 || f->lccoef == 0) {
	ILL_IFFREE (f->lcindx, int);
	mpq_EGlpNumFreeArray (f->lccoef);
	f->lc_space = nzcnt * f->lc_space_mul;
	ILL_SAFE_MALLOC (f->lcindx, f->lc_space, int);
	f->lccoef = mpq_EGlpNumAllocArray (f->lc_space);
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
    ILL_RETURN (rval, "mpq_create_factor_space");
}

static int mpq_init_matrix (mpq_factor_work * f,
      int *basis,
      int *cbeg,
      int *clen,
      int *in_ucindx,
      mpq_t * in_uccoef)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_ur_info *ur_inf = f->ur_inf;
    int dim = f->dim;
    int max_k = f->max_k;
    int *ucindx;
    int *urindx;
    mpq_t *urcoef;
    int nzcnt;
    int beg;
    int i;
    int j;
    int r;
    int rval = 0;
    mpq_t v;
    mpq_t max;
    mpq_EGlpNumInitVar (v);
    mpq_EGlpNumInitVar (max);

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

    rval = mpq_create_factor_space (f);
    ILL_CLEANUP_IF (rval);

    urindx = f->urindx;
    ucindx = f->ucindx;
    urcoef = f->urcoef;

    for (i = 0; i < dim; i++) {
	nzcnt = clen[basis[i]];
	beg = cbeg[basis[i]];
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumCopy (v, in_uccoef[beg + j]);
	    if (mpq_EGlpNumIsEqual (v, mpq_zeroLpNum, f->fzero_tol))
		continue;
	    r = in_ucindx[beg + j];
	    ucindx[uc_inf[i].nzcnt++] = r;
	    urindx[ur_inf[r].nzcnt] = i;
	    mpq_EGlpNumCopy (urcoef[ur_inf[r].nzcnt], v);
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
	mpq_EGlpNumZero (max);
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumSetToMaxAbs (max, urcoef[beg + j]);
	}
	mpq_EGlpNumCopy (ur_inf[i].max, max);
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

#ifdef mpq_TRACK_FACTOR
    mpq_EGlpNumZero (max);
    nzcnt = 0;
    for (i = 0; i < dim; i++) {
	if (mpq_EGlpNumIsLess (max, ur_inf[i].max))
	    mpq_EGlpNumCopy (max, ur_inf[i].max);
	nzcnt += ur_inf[i].nzcnt;
    }

    mpq_EGlpNumCopy (f->maxelem_orig, max);
    f->nzcnt_orig = nzcnt;
    mpq_EGlpNumCopy (f->maxelem_factor, f->maxelem_orig);
    f->nzcnt_factor = f->nzcnt_orig;
#endif				/* mpq_TRACK_FACTOR */

    /* sentinal for column space */
    ucindx[f->uc_space] = 0;

    mpq_clear_work (f);

CLEANUP:
    mpq_EGlpNumClearVar (max);
    mpq_EGlpNumClearVar (v);
    ILL_RETURN (rval, "mpq_init_matrix");
}

static int mpq_build_iteration_u_data (mpq_factor_work * f)
{
    int dim = f->dim;
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_t *uccoef = 0;
    int *ucindx = 0;
    int *urindx = f->urindx;
    mpq_t *urcoef = f->urcoef;
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

#ifdef mpq_TRACK_FACTOR
    f->nzcnt_factor = nzcnt;
#endif				/* mpq_TRACK_FACTOR */

    mpq_EGlpNumFreeArray (f->uccoef);
    uccoef = mpq_EGlpNumAllocArray (nzcnt);
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
	    mpq_EGlpNumCopy (uccoef[cbeg + cnzcnt], uccoef[cbeg]);
	    ucrind[cbeg + cnzcnt] = ucrind[cbeg];
	    urcind[ur_inf[ucindx[cbeg]].rbeg + ucrind[cbeg]] = cnzcnt;
	}
	ucindx[cbeg] = i;
	mpq_EGlpNumCopy (uccoef[cbeg], urcoef[beg]);
	ucrind[cbeg] = 0;
	urcind[beg] = 0;
	uc_inf[k].nzcnt = cnzcnt + 1;
	for (j = 1; j < nzcnt; j++) {
	    k = urindx[beg + j];
	    cbeg = uc_inf[k].cbeg;
	    cnzcnt = uc_inf[k].nzcnt;
	    ucindx[cbeg + cnzcnt] = i;
	    mpq_EGlpNumCopy (uccoef[cbeg + cnzcnt], urcoef[beg + j]);
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

    mpq_clear_work (f);

    er_space = f->er_space_mul * f->etamax;
    ILL_SAFE_MALLOC (f->er_inf, f->etamax, mpq_er_info);
    ILL_SAFE_MALLOC (f->erindx, er_space, int);
    f->ercoef = mpq_EGlpNumAllocArray (er_space);
    f->etacnt = 0;
    f->er_freebeg = 0;
    f->er_space = er_space;

    rval = 0;

CLEANUP:
    ILL_RETURN (rval, "mpq_build_iteration_u_data");
}

static int mpq_build_iteration_l_data (mpq_factor_work * f)
{
    int dim = f->dim;
    mpq_lc_info *lc_inf = f->lc_inf;
    mpq_lr_info *lr_inf = f->lr_inf;
    mpq_t *lrcoef = 0;
    int *lrindx = 0;
    mpq_t *lccoef = f->lccoef;
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

    mpq_EGlpNumFreeArray (f->lrcoef);
    if (nzcnt) {
	lrcoef = mpq_EGlpNumAllocArray (nzcnt);
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
	    mpq_EGlpNumCopy (lrcoef[rbeg + rnzcnt], lccoef[beg + j]);
	    lr_inf[k].nzcnt++;
	}
    }

#ifdef mpq_TRACK_FACTOR
    nzcnt = f->nzcnt_factor;
    for (i = 0; i < dim; i++) {
	nzcnt += lc_inf[i].nzcnt;
    }
    f->nzcnt_factor = nzcnt;

    mpq_EGlpNumCopy (f->maxelem_cur, f->maxelem_factor);
    f->nzcnt_cur = f->nzcnt_factor;

    /*
        mpq_dump_factor_stats (f);
        printf ("orig max  %e nzcnt %d\n", f->maxelem_orig, f->nzcnt_orig);
        printf ("f maxelem %e nzcnt %d\n", f->maxelem_cur, f->nzcnt_cur);
    */
#endif				/* mpq_TRACK_FACTOR */

    rval = 0;

CLEANUP:
    ILL_RETURN (rval, "mpq_build_iteration_l_data");
}

static int mpq_handle_singularity (mpq_factor_work * f)
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
    ILL_RETURN (rval, "mpq_handle_singularity");
}

static int mpq_dense_build_matrix (mpq_factor_work * f)
{
    mpq_t *dmat = 0;
    int stage = f->stage;
    int drows = f->nstages - stage;
    int dcols = f->dim - stage;
    int dsize = drows * dcols;
    int *crank = f->crank;
    mpq_t *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int nzcnt;
    int beg;
    int i;
    int r;
    int j;
    int rval = 0;

    dmat = mpq_EGlpNumAllocArray (dsize);

    for (i = 0; i < dsize; i++)
	mpq_EGlpNumZero (dmat[i]);

    for (i = 0; i < drows; i++) {
	r = f->rperm[i + stage];
	nzcnt = f->ur_inf[r].nzcnt;
	beg = f->ur_inf[r].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumCopy (dmat[i * dcols - stage + crank[urindx[beg + j]]],
		urcoef[beg + j]);
	}
    }

    f->drows = drows;
    f->dcols = dcols;
    f->dense_base = f->stage;
    f->dmat = dmat;
    dmat = 0;

    /* CLEANUP: */
    mpq_EGlpNumFreeArray (dmat);
    ILL_RETURN (rval, "mpq_dense_build_matrix");
}

static int mpq_dense_find_pivot (mpq_factor_work * f,
      int *p_r,
      int *p_c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    mpq_t *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    mpq_t maxval;
    int max_r;
    int max_c;
    int i;
    mpq_EGlpNumInitVar (maxval);
    mpq_EGlpNumZero (maxval);
    max_r = -1;
    for (i = s; i < drows; i++) {
	if (mpq_EGlpNumIsLess (maxval, ur_inf[rperm[dense_base + i]].max)) {
	    mpq_EGlpNumCopy (maxval, ur_inf[rperm[dense_base + i]].max);
	    max_r = i;
	}
    }
    if (max_r == -1) {
	return E_NO_PIVOT;
    }
    mpq_EGlpNumZero (maxval);
    max_c = -1;
    for (i = s; i < drows; i++) {
	mpq_EGlpNumSetToMaxAbsAndDo (maxval, dmat[max_r * dcols + i], max_c = i);
    }
    if (max_c == -1) {
	return E_NO_PIVOT;
    }
    *p_r = max_r;
    *p_c = max_c;

    mpq_EGlpNumClearVar (maxval);
    return 0;
}

static void mpq_dense_swap (mpq_factor_work * f,
      int r,
      int c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    mpq_t *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    int i;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    if (r != s) {
	mpq_ILL_SWAP (f->rperm[dense_base + s], f->rperm[dense_base + r], i);
	f->rrank[f->rperm[dense_base + s]] = dense_base + s;
	f->rrank[f->rperm[dense_base + r]] = dense_base + r;
	for (i = 0; i < dcols; i++) {
	    mpq_EGLPNUM_SWAP (dmat[s * dcols + i], dmat[r * dcols + i], v);
	}
    }
    if (c != s) {
	mpq_ILL_SWAP (f->cperm[dense_base + s], f->cperm[dense_base + c], i);
	f->crank[f->cperm[dense_base + s]] = dense_base + s;
	f->crank[f->cperm[dense_base + c]] = dense_base + c;
	for (i = 0; i < drows; i++) {
	    mpq_EGLPNUM_SWAP (dmat[i * dcols + s], dmat[i * dcols + c], v);
	}
    }
    mpq_EGlpNumClearVar (v);
}

static void mpq_dense_elim (mpq_factor_work * f,
      int r,
      int c)
{
    int dcols = f->dcols;
    int drows = f->drows;
    mpq_t *dmat = f->dmat;
    int dense_base = f->dense_base;
    int s = f->stage - dense_base;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int i;
    int j;
    mpq_t pivval;
    mpq_t max;
    mpq_t v;
    mpq_t w;
#ifdef mpq_TRACK_FACTOR
    mpq_t maxelem_factor;
    mpq_EGlpNumInitVar (maxelem_factor);
    mpq_EGlpNumCopy (maxelem_factor, f->maxelem_factor);
#endif
    mpq_EGlpNumInitVar (pivval);
    mpq_EGlpNumInitVar (max);
    mpq_EGlpNumInitVar (v);
    mpq_EGlpNumInitVar (w);

    mpq_dense_swap (f, r, c);
    f->stage++;
    mpq_EGlpNumCopyFrac (pivval, mpq_oneLpNum, dmat[s * dcols + s]);
    for (i = s + 1; i < drows; i++) {
	mpq_EGlpNumCopy (v, dmat[i * dcols + s]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    mpq_EGlpNumMultTo (v, pivval);
	    if (mpq_EGlpNumIsNeqZero (v, f->fzero_tol)) {
		mpq_EGlpNumCopy (dmat[i * dcols + s], v);
#ifdef mpq_TRACK_FACTOR
		mpq_EGlpNumSetToMaxAbs (maxelem_factor, v);
#endif
		mpq_EGlpNumZero (max);
		for (j = s + 1; j < drows; j++) {
		    mpq_EGlpNumCopy (w, dmat[i * dcols + j]);
		    mpq_EGlpNumSubInnProdTo (w, v, dmat[s * dcols + j]);
		    mpq_EGlpNumCopy (dmat[i * dcols + j], w);
		    mpq_EGlpNumSetToMaxAbs (max, w);
		}
		for (j = drows; j < dcols; j++) {
		    mpq_EGlpNumCopy (w, dmat[i * dcols + j]);
		    mpq_EGlpNumSubInnProdTo (w, v, dmat[s * dcols + j]);
		    mpq_EGlpNumCopy (dmat[i * dcols + j], w);
		}
		mpq_EGlpNumCopy (ur_inf[rperm[dense_base + i]].max, max);
#ifdef mpq_TRACK_FACTOR
		if (mpq_EGlpNumIsLess (maxelem_factor, max))
		    mpq_EGlpNumCopy (maxelem_factor, max);
#endif
	    } else {
		mpq_EGlpNumZero (dmat[i * dcols + s]);
	    }
	}
    }
#ifdef mpq_TRACK_FACTOR
    mpq_EGlpNumCopy (f->maxelem_factor, maxelem_factor);
    mpq_EGlpNumClearVar (maxelem_factor);
#endif
    mpq_EGlpNumClearVar (pivval);
    mpq_EGlpNumClearVar (max);
    mpq_EGlpNumClearVar (v);
    mpq_EGlpNumClearVar (w);
}

static int mpq_dense_replace_row (mpq_factor_work * f,
      int i)
{
    int dcols = f->dcols;
    int dense_base = f->dense_base;
    mpq_t *dmat = f->dmat + i * dcols;
    mpq_t *urcoef;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *cperm = f->cperm;
    int r = f->rperm[dense_base + i];
    int *urindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j = i; j < dcols; j++) {
	if (mpq_EGlpNumIsNeqZero (dmat[j], f->fzero_tol)) {
	    nzcnt++;
	}
    }
    if (nzcnt > ur_inf[r].nzcnt) {
	if (ur_inf[r].rbeg + ur_inf[r].nzcnt == f->ur_freebeg) {
	    f->ur_freebeg = ur_inf[r].rbeg;
	}
	ur_inf[r].nzcnt = 0;
	if (f->ur_freebeg + nzcnt > f->ur_space) {
	    rval = mpq_make_ur_space (f, nzcnt);
	    ILL_CLEANUP_IF (rval);
	}
	ur_inf[r].rbeg = f->ur_freebeg;
	f->ur_freebeg += nzcnt;
    }
    beg = ur_inf[r].rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    for (j = i; j < dcols; j++) {
	if (mpq_EGlpNumIsNeqZero (dmat[j], f->fzero_tol)) {
	    mpq_EGlpNumCopy (urcoef[beg], dmat[j]);
	    urindx[beg] = cperm[dense_base + j];
	    beg++;
	}
    }
    ur_inf[r].nzcnt = beg - ur_inf[r].rbeg;
CLEANUP:
    ILL_RETURN (rval, "mpq_dense_replace_row");
}

static int mpq_dense_create_col (mpq_factor_work * f,
      int i)
{
    int dcols = f->dcols;
    int drows = f->drows;
    int dense_base = f->dense_base;
    mpq_t *dmat = f->dmat;
    mpq_t *lccoef;
    mpq_lc_info *lc_inf = f->lc_inf;
    int *rperm = f->rperm;
    int *lcindx;
    int nzcnt;
    int beg;
    int j;
    int rval = 0;

    nzcnt = 0;
    for (j = i + 1; j < drows; j++) {
	if (mpq_EGlpNumIsNeqZero (dmat[j * dcols + i], f->fzero_tol)) {
	    nzcnt++;
	}
    }

    if (f->lc_freebeg + nzcnt >= f->lc_space) {
	rval = mpq_make_lc_space (f, nzcnt);
	ILL_CLEANUP_IF (rval);
    }
    beg = f->lc_freebeg;
    lc_inf[dense_base + i].cbeg = beg;
    lc_inf[dense_base + i].c = rperm[dense_base + i];
    lcindx = f->lcindx;
    lccoef = f->lccoef;

    for (j = i + 1; j < drows; j++) {
	if (mpq_EGlpNumIsNeqZero (dmat[j * dcols + i], f->fzero_tol)) {
	    mpq_EGlpNumCopy (lccoef[beg], dmat[j * dcols + i]);
	    lcindx[beg] = rperm[dense_base + j];
	    beg++;
	}
    }
    lc_inf[dense_base + i].nzcnt = beg - lc_inf[dense_base + i].cbeg;
    f->lc_freebeg = beg;
CLEANUP:
    ILL_RETURN (rval, "create_col");
}

static int mpq_dense_replace (mpq_factor_work * f)
{
    int drows = f->drows;
    int rval = 0;
    int i;

    for (i = 0; i < drows; i++) {
	rval = mpq_dense_replace_row (f, i);
	ILL_CLEANUP_IF (rval);
	rval = mpq_dense_create_col (f, i);
	ILL_CLEANUP_IF (rval);
    }
    mpq_EGlpNumFreeArray (f->dmat);
    f->drows = 0;
    f->dcols = 0;
CLEANUP:
    ILL_RETURN (rval, "mpq_dense_replace");
}

static int mpq_dense_factor (mpq_factor_work * f)
{
    int r;
    int c;
    int rval = 0;
#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
    double tmpsize;
#endif
#endif

    /*
        printf ("dense kernel, %d rows, %d  cols...\n", f->nstages - f->stage,
                f->dim - f->stage);
        fflush (stdout);
    */

    rval = mpq_dense_build_matrix (f);
    ILL_CLEANUP_IF (rval);

#ifdef mpq_FACTOR_DEBUG
#if (mpq_FACTOR_DEBUG+0>1)
    printf ("before Dense mpq_ILLfactor\n");
    mpq_dump_matrix (f, 1);
    fflush (stdout);
#endif
#endif

    while (f->stage < f->nstages) {
	r = f->stage - f->dense_base;
	rval = mpq_dense_find_pivot (f, &r, &c);
	if (rval == E_NO_PIVOT) {
	    rval = mpq_handle_singularity (f);
	    ILL_CLEANUP_IF (rval);
	    return E_SINGULAR_INTERNAL;
	} else {
	    ILL_CLEANUP_IF (rval);
	}
#ifdef mpq_FACTOR_DEBUG
#if (mpq_FACTOR_DEBUG+0>2)
	printf ("dense pivot elem: %d %d\n", r, c);
	fflush (stdout);
#endif
#endif				/* mpq_FACTOR_DEBUG */
	mpq_dense_elim (f, r, c);

#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
	tmpsize = f->maxmult * mpq_EGlpNumToLf (f->maxelem_orig);
	if (tmpsize < mpq_EGlpNumToLf (f->maxelem_factor) &&
	    mpq_EGlpNumIsLess (f->partial_cur, mpq_oneLpNum)) {
	    return E_FACTOR_BLOWUP;
	}
#endif				/* mpq_NOTICE_BLOWUP */
#endif				/* mpq_TRACK_FACTOR */

#ifdef mpq_FACTOR_DEBUG
#if (mpq_FACTOR_DEBUG+0>1)
	printf ("After dense pivot stage %d (%d) of %d (%d)\n",
	    f->stage - f->dense_base, f->stage,
	    f->nstages - f->dense_base, f->nstages);
	fflush (stdout);
#endif
#if (mpq_FACTOR_DEBUG+0>2)
	mpq_dump_matrix (f, 1);
	fflush (stdout);
#endif
#endif				/* mpq_FACTOR_DEBUG */
    }

#ifdef mpq_FACTOR_DEBUG
    printf ("After dense mpq_ILLfactor:\n");
    mpq_dump_matrix (f, 0);
    fflush (stdout);
#endif				/* mpq_FACTOR_DEBUG */

    rval = mpq_dense_replace (f);
    ILL_CLEANUP_IF (rval);

#ifdef mpq_FACTOR_DEBUG
    printf ("After replacement:\n");
    mpq_dump_matrix (f, 0);
    fflush (stdout);
#endif				/* mpq_FACTOR_DEBUG */

CLEANUP:
    ILL_RETURN (rval, "mpq_dense_factor");
}

#ifdef mpq_RECORD
FILE *mpq_fsave = 0;
int mpq_fsavecnt = 0;
#endif				/* mpq_RECORD */

static int mpq_ILLfactor_try (mpq_factor_work * f,
      int *basis,
      int *cbeg,
      int *clen,
      int *cindx,
      mpq_t * ccoef)
{
    int rval = 0;
    int r;
    int c;
#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
    mpq_t tmpsize;
    mpq_EGlpNumInitVar (tmpsize);
#endif
#endif

#ifdef mpq_RECORD
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
	if (mpq_fsave)
	    fclose (mpq_fsave);
	sprintf (fnambuf, "prob.mat.%d", mpq_fsavecnt);
	mpq_fsavecnt++;
	mpq_fsave = fopen (fnambuf, "w");
	fprintf (mpq_fsave, "%d %d %d\n", f->dim, ncol, nzcnt);
	for (i = 0; i < dim; i++) {
	    fprintf (mpq_fsave, "%d ", basis[i]);
	}
	fprintf (mpq_fsave, "\n");
	for (i = 0; i < ncol; i++) {
	    fprintf (mpq_fsave, "%d", clen[i]);
	    for (j = 0; j < clen[i]; j++) {
		fprintf (mpq_fsave, " %d %.16lg", cindx[cbeg[i] + j],
		    mpq_EGlpNumToLf (ccoef[cbeg[i] + j]));
	    }
	    fprintf (mpq_fsave, "\n");
	}
	fprintf (mpq_fsave, "\n");
	fflush (mpq_fsave);
    }
#endif				/* mpq_RECORD */

    rval = mpq_init_matrix (f, basis, cbeg, clen, cindx, ccoef);
    ILL_CLEANUP_IF (rval);

    f->stage = 0;
    f->nstages = f->dim;

#ifdef mpq_FACTOR_DEBUG
    printf ("Initial matrix:\n");
#if (mpq_FACTOR_DEBUG+0>1)
    mpq_dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif				/* mpq_FACTOR_DEBUG */
#ifdef mpq_FACTOR_STATS
    printf ("Initial matrix: ");
    mpq_dump_factor_stats (f);
#endif				/* mpq_FACTOR_STATS */

    while (f->stage < f->nstages) {
	rval = mpq_find_pivot (f, &r, &c);
	if (rval == E_NO_PIVOT) {
	    rval = mpq_handle_singularity (f);
	    ILL_CLEANUP_IF (rval);
	    return 0;
	} else {
	    ILL_CLEANUP_IF (rval);
	}
	if (f->ur_inf[r].pivcnt > f->dense_fract * (f->nstages - f->stage) &&
	    f->uc_inf[c].nzcnt > f->dense_fract * (f->nstages - f->stage) &&
	    f->nstages - f->stage > f->dense_min) {
	    rval = mpq_dense_factor (f);
	    if (rval == E_SINGULAR_INTERNAL)
		return 0;
	    if (rval)
		return rval;
	    break;
	}
#ifdef mpq_FACTOR_DEBUG
	printf ("pivot elem: %d %d\n", r, c);
	fflush (stdout);
#endif				/* mpq_FACTOR_DEBUG */
	rval = mpq_elim (f, r, c);
	ILL_CLEANUP_IF (rval);

#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
	mpq_EGlpNumSet (tmpsize, f->maxmult);
	mpq_EGlpNumMultTo (tmpsize, f->maxelem_orig);
	if (mpq_EGlpNumIsLess (tmpsize, f->maxelem_factor) &&
	    mpq_EGlpNumIsLess (f->partial_cur, mpq_oneLpNum)) {
	    return E_FACTOR_BLOWUP;
	}
#endif				/* mpq_NOTICE_BLOWUP */
#endif				/* mpq_TRACK_FACTOR */

#ifdef mpq_FACTOR_DEBUG
#if (mpq_FACTOR_DEBUG+0>3)
	printf ("After pivot stage %d of %d\n", f->stage, f->nstages);
	mpq_dump_matrix (f, 0);
	fflush (stdout);
#endif
#endif				/* mpq_FACTOR_DEBUG */
    }

    rval = mpq_build_iteration_u_data (f);
    ILL_CLEANUP_IF (rval);

    rval = mpq_build_iteration_l_data (f);
    ILL_CLEANUP_IF (rval);

#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
    mpq_EGlpNumSet (tmpsize, f->minmult);
    mpq_EGlpNumMultTo (tmpsize, f->maxelem_orig);
    if (mpq_EGlpNumIsLess (f->maxelem_factor, tmpsize) &&
	mpq_EGlpNumIsLess (f->partial_tol, f->partial_cur)) {
	if (mpq_EGlpNumIsGreaDbl (f->partial_cur, 0.5)) {
	    mpq_EGlpNumSet (f->partial_cur, 0.5);
	} else if (mpq_EGlpNumIsGreaDbl (f->partial_cur, 0.25)) {
	    mpq_EGlpNumSet (f->partial_cur, 0.25);
	} else if (mpq_EGlpNumIsGreaDbl (f->partial_cur, 0.1)) {
	    mpq_EGlpNumSet (f->partial_cur, 0.1);
	} else {
	    mpq_EGlpNumDivUiTo (f->partial_cur, 10);
	}
	if (mpq_EGlpNumIsLess (f->partial_cur, f->partial_tol)) {
	    mpq_EGlpNumCopy (f->partial_cur, f->partial_tol);
	}
	/* Bico - comment out for dist fprintf (stderr, "factor good,
	   lowering partial tolerance to %.2f\n", f->partial_cur); */
    }
#endif				/* mpq_NOTICE_BLOWUP */
#endif				/* mpq_TRACK_FACTOR */

#ifdef mpq_FACTOR_DEBUG
    printf ("Factored matrix:\n");
#if (mpq_FACTOR_DEBUG+0>1)
    mpq_dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif				/* mpq_FACTOR_DEBUG */

#ifdef mpq_FACTOR_STATS
    printf ("Factored matrix: ");
    mpq_dump_factor_stats (f);
#endif				/* mpq_FACTOR_STATS */
CLEANUP:
#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
    mpq_EGlpNumClearVar (tmpsize);
#endif
#endif
    ILL_RETURN (rval, "factor_try");
}

int mpq_ILLfactor (mpq_factor_work * f,
      int *basis,
      int *cbeg,
      int *clen,
      int *cindx,
      mpq_t * ccoef,
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
    rval = mpq_ILLfactor_try (f, basis, cbeg, clen, cindx, ccoef);
    if (rval == E_FACTOR_BLOWUP) {
	if (mpq_EGlpNumIsLessDbl (f->partial_cur, 0.1)) {
	    mpq_EGlpNumMultUiTo (f->partial_cur, 10);
	} else if (mpq_EGlpNumIsLessDbl (f->partial_cur, 0.25)) {
	    mpq_EGlpNumSet (f->partial_cur, 0.25);
	} else if (mpq_EGlpNumIsLessDbl (f->partial_cur, 0.5)) {
	    mpq_EGlpNumSet (f->partial_cur, 0.5);
	} else if (mpq_EGlpNumIsLess (f->partial_cur, mpq_oneLpNum)) {
	    mpq_EGlpNumOne (f->partial_cur);
	} else {
	    ILL_RETURN (rval, "mpq_ILLfactor");
	}
	/* Bico - comment out for dist fprintf (stderr, "factor blowup,
	   changing partial tolerance to %.2f\n", f->partial_cur); */
	goto AGAIN;
    }
    ILL_RETURN (rval, "mpq_ILLfactor");
}

static void mpq_ILLfactor_ftranl (mpq_factor_work * f,
      mpq_t * a)
{
    int *lcindx = f->lcindx;
    mpq_lc_info *lc_inf = f->lc_inf;
    mpq_t *lccoef = f->lccoef;
    int dim = f->dim;
    int beg;
    int nzcnt;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = 0; i < dim; i++) {
	mpq_EGlpNumCopy (v, a[lc_inf[i].c]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    nzcnt = lc_inf[i].nzcnt;
	    beg = lc_inf[i].cbeg;
	    for (j = 0; j < nzcnt; j++) {
		mpq_EGlpNumSubInnProdTo (a[lcindx[beg + j]], v, lccoef[beg + j]);
	    }
	}
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 > 1)
	printf ("mpq_ILLfactor_ftran a after l %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", mpq_EGlpNumToLf (a[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    }
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 <= 1)
    printf ("mpq_ILLfactor_ftran a after l:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", mpq_EGlpNumToLf (a[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}

static void mpq_ftranl3_delay (mpq_factor_work * f,
      int c)
{
    mpq_lc_info *lc_inf = f->lc_inf;
    int nzcnt;
    int *indx;
    int i;

    c = lc_inf[c].crank;
    nzcnt = lc_inf[c].nzcnt;
    indx = f->lcindx + lc_inf[c].cbeg;
    for (i = 0; i < nzcnt; i++) {
	c = indx[i];
	if (lc_inf[c].delay++ == 0) {
	    mpq_ftranl3_delay (f, c);
	}
    }
}

static void mpq_ftranl3_delay2 (mpq_factor_work * f,
      int c)
{
    mpq_lc_info *lc_inf = f->lc_inf;
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
		    mpq_ftranl3_delay2 (f, last);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
}

static void mpq_ftranl3_process (mpq_factor_work * f,
      int c,
      mpq_svector * x)
{
    mpq_lc_info *lc_inf = f->lc_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    int i;
    mpq_t *coef;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    mpq_EGlpNumCopy (v, work[c]);
    mpq_EGlpNumZero (work[c]);
    if (mpq_EGlpNumIsNeqqZero (v)) {
	x->indx[x->nzcnt] = c;
	mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    c = lc_inf[c].crank;
    nzcnt = lc_inf[c].nzcnt;
    indx = f->lcindx + lc_inf[c].cbeg;
    coef = f->lccoef + lc_inf[c].cbeg;
    for (i = 0; i < nzcnt; i++) {
	c = indx[i];
	mpq_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	if (--lc_inf[c].delay == 0) {
	    mpq_ftranl3_process (f, c, x);
	}
    }
    mpq_EGlpNumClearVar (v);
}

static void mpq_ftranl3_process2 (mpq_factor_work * f,
      int c,
      mpq_svector * x)
{
    mpq_lc_info *lc_inf = f->lc_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    int last;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    do {
	mpq_EGlpNumCopy (v, work[c]);
	mpq_EGlpNumZero (work[c]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    x->indx[x->nzcnt] = c;
	    mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	c = lc_inf[c].crank;
	nzcnt = lc_inf[c].nzcnt;
	indx = f->lcindx + lc_inf[c].cbeg;
	coef = f->lccoef + lc_inf[c].cbeg;
	last = -1;
	for (i = 0; i < nzcnt; i++) {
	    c = indx[i];
	    mpq_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	    if (--lc_inf[c].delay == 0) {
		if (last >= 0) {
		    mpq_ftranl3_process2 (f, last, x);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_ftranl3 (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * x)
{
    mpq_t *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    mpq_lc_info *lc_inf = f->lc_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (lc_inf[aindx[i]].delay++ == 0) {
	    mpq_ftranl3_delay2 (f, aindx[i]);
	}
	mpq_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--lc_inf[aindx[i]].delay == 0) {
	    mpq_ftranl3_process2 (f, aindx[i], x);
	}
    }
#ifdef mpq_SOLVE_DEBUG
    printf ("mpq_ILLfactor_ftran x after l3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* mpq_SOLVE_DEBUG */
}

static void mpq_ILLfactor_ftrane (mpq_factor_work * f,
      mpq_t * a)
{
    int *erindx = f->erindx;
    mpq_t *ercoef = f->ercoef;
    mpq_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = 0; i < etacnt; i++) {
	mpq_EGlpNumCopy (v, a[er_inf[i].r]);
	nzcnt = er_inf[i].nzcnt;
	beg = er_inf[i].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumSubInnProdTo (v, ercoef[beg + j], a[erindx[beg + j]]);
	}
	mpq_EGlpNumCopy (a[er_inf[i].r], v);
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 > 1)
	printf ("mpq_ILLfactor_ftran a after eta %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", mpq_EGlpNumToLf (a[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    }
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 <= 1)
    printf ("mpq_ILLfactor_ftran a after eta:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", mpq_EGlpNumToLf (a[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_ftrane2 (mpq_factor_work * f,
      mpq_svector * a)
{
    int *erindx = f->erindx;
    mpq_t *ercoef = f->ercoef;
    mpq_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    mpq_t *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;
    int r;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = 0; i < anzcnt; i++) {
	mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	work_indx[aindx[i]] = i + 1;
    }
    for (i = 0; i < etacnt; i++) {
	r = er_inf[i].r;
	mpq_EGlpNumCopy (v, work_coef[r]);
	nzcnt = er_inf[i].nzcnt;
	beg = er_inf[i].rbeg;
	for (j = 0; j < nzcnt; j++) {
	    mpq_EGlpNumSubInnProdTo (v, ercoef[beg + j], work_coef[erindx[beg + j]]);
	}
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    mpq_EGlpNumCopy (work_coef[r], v);
	    if (work_indx[r] == 0) {
		mpq_EGlpNumCopy (acoef[anzcnt], v);
		aindx[anzcnt] = r;
		work_indx[r] = anzcnt + 1;
		anzcnt++;
	    } else {
		mpq_EGlpNumCopy (acoef[work_indx[r] - 1], v);
	    }
	} else {
	    mpq_EGlpNumZero (work_coef[r]);
	    if (work_indx[r]) {
		mpq_EGlpNumZero (acoef[work_indx[r] - 1]);
	    }
	}
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 > 1)
	printf ("mpq_ILLfactor_ftran a after eta2 %d:", i);
	for (j = 0; j < anzcnt; j++) {
	    printf (" %.3f*%d", mpq_EGlpNumToLf (acoef[j]), aindx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    }
    i = 0;
    while (i < anzcnt) {
	mpq_EGlpNumZero (work_coef[aindx[i]]);
	work_indx[aindx[i]] = 0;
	if (mpq_EGlpNumIsNeqZero (acoef[i], f->fzero_tol)) {
	    /* if (acoef[i] > fzero_tol || acoef[i] < -fzero_tol) */
	    i++;
	} else {
	    --anzcnt;
	    mpq_EGlpNumCopy (acoef[i], acoef[anzcnt]);
	    aindx[i] = aindx[anzcnt];
	}
    }
    a->nzcnt = anzcnt;

#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 <= 1)
    printf ("mpq_ILLfactor_ftran a after eta2:");
    for (j = 0; j < anzcnt; j++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (acoef[j]), aindx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_ftranu (mpq_factor_work * f,
      mpq_t * a,
      mpq_svector * x)
{
    int *ucindx = f->ucindx;
    mpq_t *uccoef = f->uccoef;
    mpq_uc_info *uc_inf = f->uc_inf;
    int *cperm = f->cperm;
    int *rperm = f->rperm;
    int dim = f->dim;
    int xnzcnt = 0;
    int *xindx = x->indx;
    mpq_t *xcoef = x->coef;
    int nzcnt;
    int beg;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = dim - 1; i >= 0; i--) {
	mpq_EGlpNumCopy (v, a[rperm[i]]);
	if (mpq_EGlpNumIsNeqqZero (v)) {	/* ((v = a[rperm[i]]) != 0.0) */
	    j = cperm[i];
	    beg = uc_inf[j].cbeg;
	    mpq_EGlpNumDivTo (v, uccoef[beg]);
	    if (mpq_EGlpNumIsNeqZero (v, f->szero_tol)) {
		/* if (v > szero_tol || v < -szero_tol) */
		xindx[xnzcnt] = j;
		mpq_EGlpNumCopy (xcoef[xnzcnt], v);
		xnzcnt++;
	    }
	    nzcnt = uc_inf[j].nzcnt;
	    for (j = 1; j < nzcnt; j++) {
		mpq_EGlpNumSubInnProdTo (a[ucindx[beg + j]], v, uccoef[beg + j]);
	    }
	    mpq_EGlpNumZero (a[rperm[i]]);
	}
    }
    x->nzcnt = xnzcnt;
#ifdef mpq_SOLVE_DEBUG
    printf ("mpq_ILLfactor_ftran x after u:");
    for (j = 0; j < x->nzcnt; j++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (x->coef[j]), x->indx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}


static void mpq_ftranu3_delay (mpq_factor_work * f,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    int nzcnt;
    int *indx;
    int i;

    c = f->cperm[f->rrank[c]];
    nzcnt = uc_inf[c].nzcnt;
    indx = f->ucindx + uc_inf[c].cbeg;
    for (i = 1; i < nzcnt; i++) {
	c = indx[i];
	if (uc_inf[c].delay++ == 0) {
	    mpq_ftranu3_delay (f, c);
	}
    }
}

static void mpq_ftranu3_delay2 (mpq_factor_work * f,
      int c)
{
    mpq_uc_info *uc_inf = f->uc_inf;
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
		    mpq_ftranu3_delay2 (f, last);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
}

static void mpq_ftranu3_process (mpq_factor_work * f,
      int c,
      mpq_svector * x)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    mpq_EGlpNumCopy (v, work[c]);
    mpq_EGlpNumZero (work[c]);
    c = f->cperm[f->rrank[c]];
    nzcnt = uc_inf[c].nzcnt;
    indx = f->ucindx + uc_inf[c].cbeg;
    coef = f->uccoef + uc_inf[c].cbeg;
    mpq_EGlpNumDivTo (v, coef[0]);
    if (mpq_EGlpNumIsNeqZero (v, f->szero_tol))
	/* if (v > szero_tol || v < -szero_tol) */
    {
	x->indx[x->nzcnt] = c;
	mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    for (i = 1; i < nzcnt; i++) {
	c = indx[i];
	mpq_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	if (--uc_inf[c].delay == 0) {
	    mpq_ftranu3_process (f, c, x);
	}
    }
    mpq_EGlpNumClearVar (v);
}

static void mpq_ftranu3_process2 (mpq_factor_work * f,
      int c,
      mpq_svector * x)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    int last;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    do {
	mpq_EGlpNumCopy (v, work[c]);
	mpq_EGlpNumZero (work[c]);
	c = f->cperm[f->rrank[c]];
	nzcnt = uc_inf[c].nzcnt;
	indx = f->ucindx + uc_inf[c].cbeg;
	coef = f->uccoef + uc_inf[c].cbeg;
	mpq_EGlpNumDivTo (v, coef[0]);
	if (mpq_EGlpNumIsNeqZero (v, f->szero_tol))
	    /* if (v > szero_tol || v < -szero_tol) */
	{
	    x->indx[x->nzcnt] = c;
	    mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    c = indx[i];
	    mpq_EGlpNumSubInnProdTo (work[c], v, coef[i]);
	    if (--uc_inf[c].delay == 0) {
		if (last >= 0) {
		    mpq_ftranu3_process2 (f, last, x);
		}
		last = c;
	    }
	}
	c = last;
    } while (c >= 0);
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_ftranu3 (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * x)
{
    mpq_t *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    mpq_uc_info *uc_inf = f->uc_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (uc_inf[aindx[i]].delay++ == 0) {
	    mpq_ftranu3_delay2 (f, aindx[i]);
	}
	mpq_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--uc_inf[aindx[i]].delay == 0) {
	    mpq_ftranu3_process2 (f, aindx[i], x);
	}
    }
#ifdef mpq_SOLVE_DEBUG
    printf ("mpq_ILLfactor_ftran x after u3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* mpq_SOLVE_DEBUG */
}

/* mpq_ILLfactor_ftran solves Bx=a for x */
void mpq_ILLfactor_ftran (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx;
    mpq_t *acoef;
    mpq_t *work_coef = f->work_coef;

#ifdef mpq_RECORD
    {
	fprintf (mpq_fsave, "f %d", a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (mpq_fsave, " %d %.16e", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	fprintf (mpq_fsave, "\n");
	fflush (mpq_fsave);
    }
#endif				/* mpq_RECORD */
#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_ftran a:");
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %d %la", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
	nzcnt = a->nzcnt;
	aindx = a->indx;
	acoef = a->coef;
	for (i = 0; i < nzcnt; i++) {
	    mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	mpq_ILLfactor_ftranl3 (f, a, &f->xtmp);
	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = f->xtmp.nzcnt;
	    aindx = f->xtmp.indx;
	    acoef = f->xtmp.coef;

	    for (i = 0; i < nzcnt; i++) {
		mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	mpq_ILLfactor_ftranl (f, work_coef);
    }

    if (sparse) {
	mpq_ILLfactor_ftrane2 (f, &f->xtmp);
	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = f->xtmp.nzcnt;
	    aindx = f->xtmp.indx;
	    acoef = f->xtmp.coef;

	    for (i = 0; i < nzcnt; i++) {
		mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	mpq_ILLfactor_ftrane (f, work_coef);
    }

    if (sparse) {
	mpq_ILLfactor_ftranu3 (f, &f->xtmp, x);
    } else {
	mpq_ILLfactor_ftranu (f, work_coef, x);
    }

#ifdef mpq_SORT_RESULTS
    mpq_sort_vector (x);
#endif

#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_ftran x:");
	for (i = 0; i < x->nzcnt; i++) {
	    printf (" %d %la", x->indx[i], mpq_EGlpNumToLf (x->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */
    return;
}

/* mpq_ILLfactor_ftran_update solves Bx=a for x, and also returns upd, where
   Ux=upd */
void mpq_ILLfactor_ftran_update (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * upd,
      mpq_svector * x)
{
    int i;
    int nzcnt;
    int dim;
    int sparse;
    int *aindx;
    mpq_t *acoef;
    mpq_t *work_coef = f->work_coef;

#ifdef mpq_RECORD
    {
	fprintf (mpq_fsave, "F %d", a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (mpq_fsave, " %d %.16e", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	fprintf (mpq_fsave, "\n");
	fflush (mpq_fsave);
    }
#endif				/* mpq_RECORD */
#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_ftran_update a:");
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %d %.3f", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
	aindx = a->indx;
	acoef = a->coef;
	nzcnt = a->nzcnt;

	for (i = 0; i < nzcnt; i++) {
	    mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	mpq_ILLfactor_ftranl3 (f, a, upd);
	if (upd->nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = upd->nzcnt;
	    aindx = upd->indx;
	    acoef = upd->coef;

	    for (i = 0; i < nzcnt; i++) {
		mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	mpq_ILLfactor_ftranl (f, work_coef);
    }

    if (sparse) {
	mpq_ILLfactor_ftrane2 (f, upd);
	if (upd->nzcnt >= SPARSE_FACTOR * f->dim) {
	    nzcnt = upd->nzcnt;
	    aindx = upd->indx;
	    acoef = upd->coef;

	    for (i = 0; i < nzcnt; i++) {
		mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	mpq_ILLfactor_ftrane (f, work_coef);
	nzcnt = 0;
	dim = f->dim;
	aindx = upd->indx;
	acoef = upd->coef;
	for (i = 0; i < dim; i++) {
	    if (mpq_EGlpNumIsNeqqZero (work_coef[i])) {
		if (mpq_EGlpNumIsNeqZero (work_coef[i], f->szero_tol))
		    /* if(work_coef[i] > szero_tol || work_coef[i] <
		       -szero_tol) */
		{
		    aindx[nzcnt] = i;
		    mpq_EGlpNumCopy (acoef[nzcnt], work_coef[i]);
		    nzcnt++;
		}
	    }
	}
	upd->nzcnt = nzcnt;
    }

    if (sparse) {
	mpq_ILLfactor_ftranu3 (f, upd, x);
    } else {
	mpq_ILLfactor_ftranu (f, work_coef, x);
    }

#ifdef mpq_SORT_RESULTS
    mpq_sort_vector (upd);
    mpq_sort_vector (x);
#endif

#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_ftran update x:");
	for (i = 0; i < x->nzcnt; i++) {
	    printf (" %d %.3f", x->indx[i], mpq_EGlpNumToLf (x->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */
}


static void mpq_ILLfactor_btranl2 (mpq_factor_work * f,
      mpq_t * x)
{
    int *lrindx = f->lrindx;
    mpq_t *lrcoef = f->lrcoef;
    mpq_lr_info *lr_inf = f->lr_inf;
    int dim = f->dim;
    int nzcnt;
    int beg;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = dim - 1; i >= 0; i--) {
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 > 1)
	printf ("mpq_ILLfactor_btran x before l2 %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", mpq_EGlpNumToLf (x[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
	mpq_EGlpNumCopy (v, x[lr_inf[i].r]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    nzcnt = lr_inf[i].nzcnt;
	    beg = lr_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		mpq_EGlpNumSubInnProdTo (x[lrindx[beg + j]], v, lrcoef[beg + j]);
	    }
	}
    }
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 <= 1)
    printf ("mpq_ILLfactor_btran x after l2:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", mpq_EGlpNumToLf (x[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}

static void mpq_btranl3_delay (mpq_factor_work * f,
      int r)
{
    mpq_lr_info *lr_inf = f->lr_inf;
    int nzcnt;
    int *indx;
    int i;

    r = lr_inf[r].rrank;
    nzcnt = lr_inf[r].nzcnt;
    indx = f->lrindx + lr_inf[r].rbeg;
    for (i = 0; i < nzcnt; i++) {
	r = indx[i];
	if (lr_inf[r].delay++ == 0) {
	    mpq_btranl3_delay (f, r);
	}
    }
}

static void mpq_btranl3_delay2 (mpq_factor_work * f,
      int r)
{
    mpq_lr_info *lr_inf = f->lr_inf;
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
		    mpq_btranl3_delay2 (f, last);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
}

static void mpq_btranl3_process (mpq_factor_work * f,
      int r,
      mpq_svector * x)
{
    mpq_lr_info *lr_inf = f->lr_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    mpq_EGlpNumCopy (v, work[r]);
    mpq_EGlpNumZero (work[r]);
    if (mpq_EGlpNumIsNeqZero (v, f->szero_tol))
	/* if (v > szero_tol || v < -szero_tol) */
    {
	x->indx[x->nzcnt] = r;
	mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    r = lr_inf[r].rrank;
    nzcnt = lr_inf[r].nzcnt;
    indx = f->lrindx + lr_inf[r].rbeg;
    coef = f->lrcoef + lr_inf[r].rbeg;
    for (i = 0; i < nzcnt; i++) {
	r = indx[i];
	mpq_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	if (--lr_inf[r].delay == 0) {
	    mpq_btranl3_process (f, r, x);
	}
    }
    mpq_EGlpNumClearVar (v);
}

static void mpq_btranl3_process2 (mpq_factor_work * f,
      int r,
      mpq_svector * x)
{
    mpq_lr_info *lr_inf = f->lr_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    int last;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    do {
	mpq_EGlpNumCopy (v, work[r]);
	mpq_EGlpNumZero (work[r]);
	if (mpq_EGlpNumIsNeqZero (v, f->szero_tol))
	    /* if (v > szero_tol || v < -szero_tol) */
	{
	    x->indx[x->nzcnt] = r;
	    mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	r = lr_inf[r].rrank;
	nzcnt = lr_inf[r].nzcnt;
	indx = f->lrindx + lr_inf[r].rbeg;
	coef = f->lrcoef + lr_inf[r].rbeg;
	last = -1;
	for (i = 0; i < nzcnt; i++) {
	    r = indx[i];
	    mpq_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	    if (--lr_inf[r].delay == 0) {
		if (last >= 0) {
		    mpq_btranl3_process2 (f, last, x);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_btranl3 (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * x)
{
    mpq_t *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    mpq_lr_info *lr_inf = f->lr_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (lr_inf[aindx[i]].delay++ == 0) {
	    mpq_btranl3_delay2 (f, aindx[i]);
	}
	mpq_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--lr_inf[aindx[i]].delay == 0) {
	    mpq_btranl3_process2 (f, aindx[i], x);
	}
    }
#ifdef mpq_SOLVE_DEBUG
    printf ("mpq_ILLfactor_btran x after l3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* mpq_SOLVE_DEBUG */
}

static void mpq_ILLfactor_btrane (mpq_factor_work * f,
      mpq_t * x)
{
    int *erindx = f->erindx;
    mpq_t *ercoef = f->ercoef;
    mpq_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = etacnt - 1; i >= 0; i--) {
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 > 1)
	printf ("mpq_ILLfactor_btran x before eta %d:", i);
	for (j = 0; j < f->dim; j++) {
	    printf (" %.3f", mpq_EGlpNumToLf (x[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
	mpq_EGlpNumCopy (v, x[er_inf[i].r]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    nzcnt = er_inf[i].nzcnt;
	    beg = er_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		mpq_EGlpNumSubInnProdTo (x[erindx[beg + j]], v, ercoef[beg + j]);
	    }
	}
    }
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 <= 1)
    printf ("mpq_ILLfactor_btran x after eta:");
    for (j = 0; j < f->dim; j++) {
	printf (" %.3f", mpq_EGlpNumToLf (x[j]));
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_btrane2 (mpq_factor_work * f,
      mpq_svector * x)
{
    int *erindx = f->erindx;
    mpq_t *ercoef = f->ercoef;
    mpq_er_info *er_inf = f->er_inf;
    int etacnt = f->etacnt;
    int beg;
    int nzcnt;
    int xnzcnt = x->nzcnt;
    int *xindx = x->indx;
    mpq_t *xcoef = x->coef;
    mpq_t *work_coef = f->work_coef;
    int *work_indx = f->work_indx;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = 0; i < xnzcnt; i++) {
	mpq_EGlpNumCopy (work_coef[xindx[i]], xcoef[i]);
	work_indx[xindx[i]] = i + 1;
    }
    for (i = etacnt - 1; i >= 0; i--) {
#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 > 1)
	printf ("mpq_ILLfactor_btran x before eta2 %d:", i);
	for (j = 0; j < xnzcnt; j++) {
	    printf (" %.3f*%d", mpq_EGlpNumToLf (work_coef[xindx[j]]), xindx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
	mpq_EGlpNumCopy (v, work_coef[er_inf[i].r]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    nzcnt = er_inf[i].nzcnt;
	    beg = er_inf[i].rbeg;
	    for (j = 0; j < nzcnt; j++) {
		if (work_indx[erindx[beg + j]] == 0) {
		    work_indx[erindx[beg + j]] = xnzcnt;
		    xindx[xnzcnt++] = erindx[beg + j];
		}
		mpq_EGlpNumSubInnProdTo (work_coef[erindx[beg + j]], v, ercoef[beg + j]);
	    }
	}
    }

    j = 0;
    while (j < xnzcnt) {
	mpq_EGlpNumCopy (xcoef[j], work_coef[xindx[j]]);
	mpq_EGlpNumZero (work_coef[xindx[j]]);
	work_indx[xindx[j]] = 0;
	if (mpq_EGlpNumIsEqqual (xcoef[j], mpq_zeroLpNum)) {
	    --xnzcnt;
	    xindx[j] = xindx[xnzcnt];
	} else {
	    j++;
	}
    }
    x->nzcnt = xnzcnt;

#ifdef mpq_SOLVE_DEBUG
#if (mpq_SOLVE_DEBUG+0 <= 1)
    printf ("mpq_ILLfactor_btran x after eta2:");
    for (j = 0; j < xnzcnt; j++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (xcoef[j]), xindx[j]);
    }
    printf ("\n");
    fflush (stdout);
#endif
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_btranu (mpq_factor_work * f,
      mpq_t * a,
      mpq_svector * x)
{
    int *urindx = f->urindx;
    mpq_t *urcoef = f->urcoef;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int *cperm = f->cperm;
    int dim = f->dim;
    int xnzcnt = 0;
    int *xindx = x->indx;
    mpq_t *xcoef = x->coef;
    int nzcnt;
    int beg;
    int i;
    int j;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    for (i = 0; i < dim; i++) {
	mpq_EGlpNumCopy (v, a[cperm[i]]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    j = rperm[i];
	    beg = ur_inf[j].rbeg;
	    mpq_EGlpNumDivTo (v, urcoef[beg]);
	    if (mpq_EGlpNumIsNeqZero (v, f->szero_tol)) {	/*
																															 * if (v > szero_tol || v < -szero_tol) */
		xindx[xnzcnt] = j;
		mpq_EGlpNumCopy (xcoef[xnzcnt], v);
		xnzcnt++;
	    }
	    nzcnt = ur_inf[j].nzcnt;
	    for (j = 1; j < nzcnt; j++) {
		mpq_EGlpNumSubInnProdTo (a[urindx[beg + j]], v, urcoef[beg + j]);
	    }
	    mpq_EGlpNumZero (a[cperm[i]]);
	}
    }
    x->nzcnt = xnzcnt;
#ifdef mpq_SOLVE_DEBUG
    printf ("mpq_ILLfactor_btran x after u:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* mpq_SOLVE_DEBUG */
    mpq_EGlpNumClearVar (v);
}


static void mpq_btranu3_delay (mpq_factor_work * f,
      int r)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    int nzcnt;
    int *indx;
    int i;

    r = f->rperm[f->crank[r]];
    nzcnt = ur_inf[r].nzcnt;
    indx = f->urindx + ur_inf[r].rbeg;
    for (i = 1; i < nzcnt; i++) {
	r = indx[i];
	if (ur_inf[r].delay++ == 0) {
	    mpq_btranu3_delay (f, r);
	}
    }
}

static void mpq_btranu3_delay2 (mpq_factor_work * f,
      int r)
{
    mpq_ur_info *ur_inf = f->ur_inf;
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
		    mpq_btranu3_delay2 (f, last);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
}

static void mpq_btranu3_process (mpq_factor_work * f,
      int r,
      mpq_svector * x)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    mpq_EGlpNumCopy (v, work[r]);
    mpq_EGlpNumZero (work[r]);
    r = f->rperm[f->crank[r]];
    nzcnt = ur_inf[r].nzcnt;
    indx = f->urindx + ur_inf[r].rbeg;
    coef = f->urcoef + ur_inf[r].rbeg;
    mpq_EGlpNumDivTo (v, coef[0]);
    if (mpq_EGlpNumIsNeqqZero (v)) {
	x->indx[x->nzcnt] = r;
	mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	x->nzcnt++;
    }
    for (i = 1; i < nzcnt; i++) {
	r = indx[i];
	mpq_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	if (--ur_inf[r].delay == 0) {
	    mpq_btranu3_process (f, r, x);
	}
    }
    mpq_EGlpNumClearVar (v);
}

static void mpq_btranu3_process2 (mpq_factor_work * f,
      int r,
      mpq_svector * x)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    int last;
    mpq_t v;
    mpq_EGlpNumInitVar (v);

    do {
	mpq_EGlpNumCopy (v, work[r]);
	mpq_EGlpNumZero (work[r]);
	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	coef = f->urcoef + ur_inf[r].rbeg;
	mpq_EGlpNumDivTo (v, coef[0]);
	if (mpq_EGlpNumIsNeqqZero (v)) {
	    x->indx[x->nzcnt] = r;
	    mpq_EGlpNumCopy (x->coef[x->nzcnt], v);
	    x->nzcnt++;
	}
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    r = indx[i];
	    mpq_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	    if (--ur_inf[r].delay == 0) {
		if (last >= 0) {
		    mpq_btranu3_process2 (f, last, x);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
    mpq_EGlpNumClearVar (v);
}

static void mpq_ILLfactor_btranu3 (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * x)
{
    mpq_t *work = f->work_coef;
    int anzcnt = a->nzcnt;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    mpq_ur_info *ur_inf = f->ur_inf;
    int i;

    for (i = 0; i < anzcnt; i++) {
	if (ur_inf[aindx[i]].delay++ == 0) {
	    mpq_btranu3_delay2 (f, aindx[i]);
	}
	mpq_EGlpNumCopy (work[aindx[i]], acoef[i]);
    }
    x->nzcnt = 0;
    for (i = 0; i < anzcnt; i++) {
	if (--ur_inf[aindx[i]].delay == 0) {
	    mpq_btranu3_process2 (f, aindx[i], x);
	}
    }
#ifdef mpq_SOLVE_DEBUG
    printf ("mpq_ILLfactor_btran x after u3:");
    for (i = 0; i < x->nzcnt; i++) {
	printf (" %.3f*%d", mpq_EGlpNumToLf (x->coef[i]), x->indx[i]);
    }
    printf ("\n");
    fflush (stdout);
#endif				/* mpq_SOLVE_DEBUG */
}

/* mpq_ILLfactor_btran solves x^tB=a^t (or, B^t x = a) for x */
void mpq_ILLfactor_btran (mpq_factor_work * f,
      mpq_svector * a,
      mpq_svector * x)
{
    int i;
    int nzcnt;
    int sparse;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    mpq_t *work_coef = f->work_coef;
    int dim = f->dim;

#ifdef mpq_RECORD
    {
	fprintf (mpq_fsave, "b %d", a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (mpq_fsave, " %d %.16e", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	fprintf (mpq_fsave, "\n");
	fflush (mpq_fsave);
    }
#endif				/* mpq_RECORD */
#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_btran a:");
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %d %.3f", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

    if (a->nzcnt >= SPARSE_FACTOR * f->dim) {
	aindx = a->indx;
	acoef = a->coef;
	work_coef = f->work_coef;
	nzcnt = a->nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	mpq_ILLfactor_btranu3 (f, a, &f->xtmp);
    } else {
	mpq_ILLfactor_btranu (f, work_coef, &f->xtmp);
    }

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	aindx = f->xtmp.indx;
	acoef = f->xtmp.coef;
	work_coef = f->work_coef;
	nzcnt = f->xtmp.nzcnt;
	for (i = 0; i < nzcnt; i++) {
	    mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}
	sparse = 0;
    } else {
	sparse = 1;
    }

    if (sparse) {
	mpq_ILLfactor_btrane2 (f, &f->xtmp);
	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	    aindx = f->xtmp.indx;
	    acoef = f->xtmp.coef;
	    work_coef = f->work_coef;
	    nzcnt = f->xtmp.nzcnt;
	    for (i = 0; i < nzcnt; i++) {
		mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	    }
	    sparse = 0;
	}
    } else {
	mpq_ILLfactor_btrane (f, work_coef);
    }

    if (sparse) {
	mpq_ILLfactor_btranl3 (f, &f->xtmp, x);
    } else {
	mpq_ILLfactor_btranl2 (f, work_coef);
	dim = f->dim;
	nzcnt = 0;
	aindx = x->indx;
	acoef = x->coef;
	for (i = 0; i < dim; i++) {
	    if (mpq_EGlpNumIsNeqqZero (work_coef[i])) {
		if (mpq_EGlpNumIsNeqZero (work_coef[i], f->szero_tol))
		    /* if (work_coef[i] > szero_tol || work_coef[i] <
		       -szero_tol) */
		{
		    aindx[nzcnt] = i;
		    mpq_EGlpNumCopy (acoef[nzcnt], work_coef[i]);
		    nzcnt++;
		}
		mpq_EGlpNumZero (work_coef[i]);
	    }
	}
	x->nzcnt = nzcnt;
    }

#ifdef mpq_SORT_RESULTS
    mpq_sort_vector (x);
#endif

#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_btran x:");
	for (i = 0; i < x->nzcnt; i++) {
	    printf (" %d %.3f", x->indx[i], mpq_EGlpNumToLf (x->coef[i]));
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */
    return;
}

static int mpq_expand_col (mpq_factor_work * f,
      int col)
{
    mpq_uc_info *uc_inf = f->uc_inf + col;
    int uc_freebeg = f->uc_freebeg;
    int nzcnt = uc_inf->nzcnt;
    int cbeg;
    mpq_t *uccoef;
    int *ucindx;
    int *ucrind;
    int i;
    int rval = 0;

    if (uc_freebeg + nzcnt + 1 >= f->uc_space) {
	rval = mpq_make_uc_space (f, nzcnt + 1);
	ILL_CLEANUP_IF (rval);
	uc_freebeg = f->uc_freebeg;
    }
    cbeg = uc_inf->cbeg;
    uccoef = f->uccoef;
    ucindx = f->ucindx;
    ucrind = f->ucrind;

    for (i = 0; i < nzcnt; i++) {
	mpq_EGlpNumCopy (uccoef[uc_freebeg + i], uccoef[cbeg + i]);
	ucindx[uc_freebeg + i] = ucindx[cbeg + i];
	ucrind[uc_freebeg + i] = ucrind[cbeg + i];
	ucindx[cbeg + i] = -1;
    }

    uc_inf->cbeg = uc_freebeg;
    f->uc_freebeg = uc_freebeg + nzcnt;
CLEANUP:
    ILL_RETURN (rval, "mpq_expand_col");
}

static int mpq_expand_row (mpq_factor_work * f,
      int row)
{
    mpq_ur_info *ur_inf = f->ur_inf + row;
    int ur_freebeg = f->ur_freebeg;
    int nzcnt = ur_inf->nzcnt;
    int rbeg;
    mpq_t *urcoef;
    int *urindx;
    int *urcind;
    int i;
    int rval = 0;

    if (ur_freebeg + nzcnt + 1 >= f->ur_space) {
	rval = mpq_make_ur_space (f, nzcnt + 1);
	ILL_CLEANUP_IF (rval);
	ur_freebeg = f->ur_freebeg;
    }
    rbeg = ur_inf->rbeg;
    urcoef = f->urcoef;
    urindx = f->urindx;
    urcind = f->urcind;

    for (i = 0; i < nzcnt; i++) {
	mpq_EGlpNumCopy (urcoef[ur_freebeg + i], urcoef[rbeg + i]);
	urindx[ur_freebeg + i] = urindx[rbeg + i];
	urcind[ur_freebeg + i] = urcind[rbeg + i];
	urindx[rbeg + i] = -1;
    }

    ur_inf->rbeg = ur_freebeg;
    f->ur_freebeg = ur_freebeg + nzcnt;
CLEANUP:
    ILL_RETURN (rval, "mpq_expand_row");
}

static int mpq_add_nonzero (mpq_factor_work * f,
      int row,
      int col,
      mpq_t val)
{
    mpq_ur_info *ur_inf = f->ur_inf + row;
    mpq_uc_info *uc_inf = f->uc_inf + col;
    int cnzcnt = uc_inf->nzcnt;
    int rnzcnt = ur_inf->nzcnt;
    int cloc = uc_inf->cbeg + cnzcnt;
    int rloc = ur_inf->rbeg + rnzcnt;
    int rval = 0;

    if (f->ucindx[cloc] != -1) {
	rval = mpq_expand_col (f, col);
	ILL_CLEANUP_IF (rval);
	cloc = uc_inf->cbeg + cnzcnt;
    }
    TESTG ((rval = (rloc < 0 || rloc > f->ur_space)), CLEANUP,
	"rloc %d outside boundaries [0:%d]", rloc, f->ur_space);
    if (f->urindx[rloc] != -1) {
	rval = mpq_expand_row (f, row);
	ILL_CLEANUP_IF (rval);
	rloc = ur_inf->rbeg + rnzcnt;
    }
    f->ucindx[cloc] = row;
    mpq_EGlpNumCopy (f->uccoef[cloc], val);
    f->ucrind[cloc] = rnzcnt;
    f->urindx[rloc] = col;
    mpq_EGlpNumCopy (f->urcoef[rloc], val);
    f->urcind[rloc] = cnzcnt;

    if (cloc == f->uc_freebeg)
	f->uc_freebeg++;
    if (rloc == f->ur_freebeg)
	f->ur_freebeg++;

    uc_inf->nzcnt = cnzcnt + 1;
    ur_inf->nzcnt = rnzcnt + 1;
CLEANUP:
    ILL_RETURN (rval, "mpq_add_nonzero");
}

static void mpq_delete_nonzero_row (mpq_factor_work * f,
      int row,
      int ind)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_t *urcoef = f->urcoef;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int *ucrind = f->ucrind;
    int rbeg = ur_inf[row].rbeg;
    int nzcnt = ur_inf[row].nzcnt - 1;
    int cbeg;

    if (ind != nzcnt) {
	mpq_EGlpNumCopy (urcoef[rbeg + ind], urcoef[rbeg + nzcnt]);
	urindx[rbeg + ind] = urindx[rbeg + nzcnt];
	urcind[rbeg + ind] = urcind[rbeg + nzcnt];
	cbeg = f->uc_inf[urindx[rbeg + nzcnt]].cbeg;
	ucrind[cbeg + urcind[rbeg + nzcnt]] = ind;
	urindx[rbeg + nzcnt] = -1;
    }
    ur_inf[row].nzcnt = nzcnt;
}

static void mpq_delete_nonzero_col (mpq_factor_work * f,
      int col,
      int ind)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    mpq_t *uccoef = f->uccoef;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int *urcind = f->urcind;
    int cbeg = uc_inf[col].cbeg;
    int nzcnt = uc_inf[col].nzcnt - 1;
    int rbeg;

    if (ind != nzcnt) {
	mpq_EGlpNumCopy (uccoef[cbeg + ind], uccoef[cbeg + nzcnt]);
	ucindx[cbeg + ind] = ucindx[cbeg + nzcnt];
	ucrind[cbeg + ind] = ucrind[cbeg + nzcnt];
	rbeg = f->ur_inf[ucindx[cbeg + nzcnt]].rbeg;
	urcind[rbeg + ucrind[cbeg + nzcnt]] = ind;
	ucindx[cbeg + nzcnt] = -1;
    }
    uc_inf[col].nzcnt = nzcnt;
}

static void mpq_delete_column (mpq_factor_work * f,
      int col)
{
    mpq_uc_info *uc_inf = f->uc_inf;
    int beg = uc_inf[col].cbeg;
    int nzcnt = uc_inf[col].nzcnt;
    int *ucindx = f->ucindx + beg;
    int *ucrind = f->ucrind + beg;
    int i;

    for (i = 0; i < nzcnt; i++) {
	mpq_delete_nonzero_row (f, ucindx[i], ucrind[i]);
	ucindx[i] = -1;
    }
    uc_inf[col].nzcnt = 0;

#ifdef mpq_TRACK_FACTOR
    f->nzcnt_cur -= nzcnt;
#endif

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_delete_column corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */
}

static void mpq_delete_row (mpq_factor_work * f,
      int row,
      mpq_svector * x)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    int beg = ur_inf[row].rbeg;
    int nzcnt = ur_inf[row].nzcnt;
    int *urindx = f->urindx + beg;
    mpq_t *urcoef = f->urcoef + beg;
    int *urcind = f->urcind + beg;
    int i;

    for (i = 0; i < nzcnt; i++) {
	x->indx[i] = urindx[i];
	mpq_EGlpNumCopy (x->coef[i], urcoef[i]);
	mpq_delete_nonzero_col (f, urindx[i], urcind[i]);
	urindx[i] = -1;
    }
    x->nzcnt = nzcnt;
    ur_inf[row].nzcnt = 0;

#ifdef mpq_TRACK_FACTOR
    f->nzcnt_cur -= nzcnt;
#endif

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_delete_row corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */
}

static int mpq_create_column (mpq_factor_work * f,
      mpq_svector * a,
      int col,
      int *p_last_rank)
{
    int *rrank = f->rrank;
    int nzcnt = a->nzcnt;
    int *aindx = a->indx;
    mpq_t *acoef = a->coef;
    int i;
    int j;
    int rval = 0;
    int last_rank = -1;
#ifdef mpq_TRACK_FACTOR
    mpq_t max;
    mpq_EGlpNumInitVar (max);
    mpq_EGlpNumCopy (max, f->maxelem_cur);
#endif				/* mpq_TRACK_FACTOR */

    last_rank = 0;

    for (i = 0; i < nzcnt; i++) {
	rval = mpq_add_nonzero (f, aindx[i], col, acoef[i]);
	ILL_CLEANUP_IF (rval);
#ifdef mpq_TRACK_FACTOR
	mpq_EGlpNumSetToMaxAbs (max, acoef[i]);
#endif				/* mpq_TRACK_FACTOR */
	j = rrank[aindx[i]];
	if (j > last_rank)
	    last_rank = j;
    }
    *p_last_rank = last_rank;

#ifdef mpq_TRACK_FACTOR
    f->nzcnt_cur += nzcnt;
    mpq_EGlpNumCopy (f->maxelem_cur, max);
    mpq_EGlpNumClearVar (max);
#endif				/* mpq_TRACK_FACTOR */

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_create_column corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

CLEANUP:
    ILL_RETURN (rval, "mpq_create_column");
}

#ifdef mpq_UPDATE_STUDY
static int mpq_column_rank (mpq_factor_work * f,
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

static void mpq_shift_permutations (mpq_factor_work * f,
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

static int mpq_eliminate_row (mpq_factor_work * f,
      int rank_p,
      int rank_r)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    int *rperm = f->rperm;
    int *cperm = f->cperm;
    int *urindx = f->urindx;
    mpq_t *urcoef = f->urcoef;
    int *erindx = f->erindx;
    mpq_t *ercoef = f->ercoef;
    mpq_t *work_coef = f->work_coef;
    int er_freebeg = f->er_freebeg;
    int er_space = f->er_space;
    int beg;
    int nzcnt;
    int i;
    int j;
    int c;
    int r;
    mpq_t pivot_mul;
#ifdef mpq_TRACK_FACTOR
    mpq_t max;
    mpq_EGlpNumInitVar (max);
    mpq_EGlpNumCopy (max, f->maxelem_cur);
#endif				/* mpq_TRACK_FACTOR */
    mpq_EGlpNumInitVar (pivot_mul);

    for (i = rank_p; i < rank_r; i++) {
	c = cperm[i];
	if (mpq_EGlpNumIsNeqZero (work_coef[c], f->fzero_tol)) {	/*
																																					 * if (work_coef[c] > fzero_tol || work_coef[c] < -fzero_tol) */
	    r = rperm[i];
	    beg = ur_inf[r].rbeg;
	    nzcnt = ur_inf[r].nzcnt;
	    mpq_EGlpNumCopyFrac (pivot_mul, work_coef[c], urcoef[beg]);
	    mpq_EGlpNumZero (work_coef[c]);
	    for (j = 1; j < nzcnt; j++) {
		mpq_EGlpNumSubInnProdTo (work_coef[urindx[beg + j]], pivot_mul, urcoef[beg + j]);	/* 0.85 */
	    }
	    if (er_freebeg >= er_space) {
		/* fprintf (stderr, "no space in mpq_eliminate_row\n"); */
#ifdef mpq_TRACK_FACTOR
		mpq_EGlpNumClearVar (max);
#endif
		mpq_EGlpNumClearVar (pivot_mul);
		return E_UPDATE_NOSPACE;
	    }
	    erindx[er_freebeg] = r;
	    mpq_EGlpNumCopy (ercoef[er_freebeg], pivot_mul);
#ifdef mpq_TRACK_FACTOR
	    mpq_EGlpNumSetToMaxAbs (max, pivot_mul);
#endif				/* mpq_TRACK_FACTOR */
	    er_freebeg++;
	} else {
	    mpq_EGlpNumZero (work_coef[c]);
	}
    }
    f->er_freebeg = er_freebeg;
#ifdef mpq_TRACK_FACTOR
    mpq_EGlpNumCopy (f->maxelem_cur, max);
    mpq_EGlpNumClearVar (max);
#endif				/* mpq_TRACK_FACTOR */
    mpq_EGlpNumClearVar (pivot_mul);
    return 0;
}

static int mpq_create_row (mpq_factor_work * f,
      mpq_t * a,
      int row,
      int minrank)
{
    int *cperm = f->cperm;
    int dim = f->dim;
    int i;
    int j;
    int rval = 0;
#ifdef mpq_TRACK_FACTOR
    mpq_t max;
    mpq_EGlpNumInitVar (max);
    mpq_EGlpNumCopy (max, f->maxelem_cur);
#endif				/* mpq_TRACK_FACTOR */

    for (i = minrank; i < dim; i++) {
	if (mpq_EGlpNumIsNeqqZero (a[cperm[i]])) {
	    j = cperm[i];
	    if (mpq_EGlpNumIsNeqZero (a[j], f->fzero_tol)) {	/*
																																	 * if (a[j] > fzero_tol || a[j] < -fzero_tol) */
		rval = mpq_add_nonzero (f, row, j, a[j]);
		ILL_CLEANUP_IF (rval);
#ifdef mpq_TRACK_FACTOR
		mpq_EGlpNumSetToMaxAbs (max, a[j]);
#endif				/* mpq_TRACK_FACTOR */
	    }
	    mpq_EGlpNumZero (a[j]);
	}
    }

#ifdef mpq_TRACK_FACTOR
    f->nzcnt_cur += f->ur_inf[row].nzcnt;
    mpq_EGlpNumCopy (f->maxelem_cur, max);
    mpq_EGlpNumClearVar (max);
#endif				/* mpq_TRACK_FACTOR */

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_create_row corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */
CLEANUP:
    ILL_RETURN (rval, "mpq_create_row");
}

static void mpq_serow_delay (mpq_factor_work * f,
      int r,
      int rank_r)
{
    mpq_ur_info *ur_inf = f->ur_inf;
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
		    mpq_serow_delay (f, last, rank_r);
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
}

static int mpq_serow_process (mpq_factor_work * f,
      int r,
      mpq_svector * newr,
      int rank_r)
{
    mpq_ur_info *ur_inf = f->ur_inf;
    mpq_t *work = f->work_coef;
    int nzcnt;
    int *indx;
    mpq_t *coef;
    int i;
    mpq_t v;
    int last;
    int rval;
    mpq_EGlpNumInitVar (v);

    do {
	mpq_EGlpNumCopy (v, work[r]);
	mpq_EGlpNumZero (work[r]);
	if (f->crank[r] >= rank_r) {
	    if (mpq_EGlpNumIsNeqZero (v, f->fzero_tol)) {	/*
																															 * if (v > fzero_tol || v < -fzero_tol) */
		/* stash this nonzero in the resulting row */
#ifdef mpq_TRACK_FACTOR
		mpq_EGlpNumSetToMaxAbs (f->maxelem_cur, v);
#endif				/* mpq_TRACK_FACTOR */
		newr->indx[newr->nzcnt] = r;
		mpq_EGlpNumCopy (newr->coef[newr->nzcnt], v);
		newr->nzcnt++;
		mpq_EGlpNumClearVar (v);
		return 0;
	    } else {
		mpq_EGlpNumClearVar (v);
		return 0;
	    }
	}
	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	coef = f->urcoef + ur_inf[r].rbeg;
	mpq_EGlpNumDivTo (v, coef[0]);
	if (mpq_EGlpNumIsNeqZero (v, f->fzero_tol)) {	/*
																													 * if (v > fzero_tol || v < -fzero_tol) */
	    /* stash v in eta */
	    if (f->er_freebeg >= f->er_space) {
		/* fprintf (stderr, "no space in mpq_eliminate_row\n"); */
		mpq_EGlpNumClearVar (v);
		return E_UPDATE_NOSPACE;
	    }
	    f->erindx[f->er_freebeg] = r;
	    mpq_EGlpNumCopy (f->ercoef[f->er_freebeg], v);
#ifdef mpq_TRACK_FACTOR
	    mpq_EGlpNumSetToMaxAbs (f->maxelem_cur, v);
#endif				/* mpq_TRACK_FACTOR */
	    f->er_freebeg++;
	}
	last = -1;
	for (i = 1; i < nzcnt; i++) {
	    r = indx[i];
	    mpq_EGlpNumSubInnProdTo (work[r], v, coef[i]);
	    if (--ur_inf[r].delay == 0) {
		if (last >= 0) {
		    rval = mpq_serow_process (f, last, newr, rank_r);
		    if (rval) {
			mpq_EGlpNumClearVar (v);
			return rval;
		    }
		}
		last = r;
	    }
	}
	r = last;
    } while (r >= 0);
    mpq_EGlpNumClearVar (v);
    return 0;
}

static int mpq_sparse_eliminate_row (mpq_factor_work * f,
      mpq_svector * x,
      int row_p,
      int rank_r)
{
    mpq_t *work = f->work_coef;
    int xnzcnt = x->nzcnt;
    int *xindx = x->indx;
    mpq_t *xcoef = x->coef;
    mpq_ur_info *ur_inf = f->ur_inf;
    int *crank = f->crank;
    int i;
    int j;
    int rval = 0;
    mpq_svector newr;

    newr.indx = 0;
    newr.coef = 0;

    for (i = 0; i < xnzcnt; i++) {
	j = xindx[i];
	if (ur_inf[j].delay++ == 0 && crank[j] < rank_r) {
	    mpq_serow_delay (f, j, rank_r);
	}
	mpq_EGlpNumCopy (work[j], xcoef[i]);
    }

    newr.nzcnt = 0;
    ILL_SAFE_MALLOC (newr.indx, f->dim, int);
    newr.coef = mpq_EGlpNumAllocArray (f->dim);

    for (i = 0; i < xnzcnt; i++) {
	j = xindx[i];
	if (--ur_inf[j].delay == 0) {
	    rval = mpq_serow_process (f, j, &newr, rank_r);
	    ILL_CLEANUP_IF (rval);
	}
    }

    for (i = 0; i < newr.nzcnt; i++) {
	rval = mpq_add_nonzero (f, row_p, newr.indx[i], newr.coef[i]);
	ILL_CLEANUP_IF (rval);
    }

#ifdef mpq_TRACK_FACTOR
    f->nzcnt_cur += newr.nzcnt;
#endif				/* mpq_TRACK_FACTOR */

CLEANUP:
    mpq_EGlpNumFreeArray (newr.coef);
    ILL_IFFREE (newr.indx, int);

    /* Bico 031210 - chg from ILL_RETURN */
    ILL_RESULT (rval, "mpq_sparse_eliminate_row");
}

static int mpq_move_pivot_row (mpq_factor_work * f,
      int r,
      int c)
{
    mpq_ur_info *ur_inf = f->ur_inf + r;
    mpq_uc_info *uc_inf = f->uc_inf;
    int beg = ur_inf->rbeg;
    int nzcnt = ur_inf->nzcnt;
    int *urindx = f->urindx;
    int *urcind = f->urcind;
    int *ucrind = f->ucrind;
    mpq_t *urcoef = f->urcoef;
    mpq_t dt;
    int it;
    int i;

    if (urindx[beg] == c)
	return 0;
    mpq_EGlpNumInitVar (dt);

    for (i = 1; i < nzcnt; i++) {
	if (urindx[beg + i] == c) {
	    mpq_EGLPNUM_SWAP (urcoef[beg], urcoef[beg + i], dt);
	    mpq_ILL_SWAP (urcind[beg], urcind[beg + i], it);
	    urindx[beg + i] = urindx[beg];
	    urindx[beg] = c;
	    ucrind[uc_inf[c].cbeg + urcind[beg]] = 0;
	    ucrind[uc_inf[urindx[beg + i]].cbeg + urcind[beg + i]] = i;
	    mpq_EGlpNumClearVar (dt);
	    return 0;
	}
    }
    fprintf (stderr, "pivot row nonzero not found\n");
    mpq_EGlpNumClearVar (dt);
    return E_UPDATE_SINGULAR_ROW;
}

static int mpq_move_pivot_col (mpq_factor_work * f,
      int c,
      int r)
{
    mpq_uc_info *uc_inf = f->uc_inf + c;
    mpq_ur_info *ur_inf = f->ur_inf;
    int beg = uc_inf->cbeg;
    int nzcnt = uc_inf->nzcnt;
    int *ucindx = f->ucindx;
    int *ucrind = f->ucrind;
    int *urcind = f->urcind;
    mpq_t *uccoef = f->uccoef;
    mpq_t dt;
    int i, it;

    if (ucindx[beg] == r)
	return 0;
    mpq_EGlpNumInitVar (dt);

    for (i = 1; i < nzcnt; i++) {
	if (ucindx[beg + i] == r) {
	    mpq_EGLPNUM_SWAP (uccoef[beg], uccoef[beg + i], dt);
	    mpq_ILL_SWAP (ucrind[beg], ucrind[beg + i], it);
	    ucindx[beg + i] = ucindx[beg];
	    ucindx[beg] = r;
	    urcind[ur_inf[r].rbeg + ucrind[beg]] = 0;
	    urcind[ur_inf[ucindx[beg + i]].rbeg + ucrind[beg + i]] = i;
	    mpq_EGlpNumClearVar (dt);
	    return 0;
	}
    }
    fprintf (stderr, "pivot col nonzero not found\n");
    mpq_EGlpNumClearVar (dt);
    return E_UPDATE_SINGULAR_COL;
}

static int mpq_move_pivot (mpq_factor_work * f,
      int rank_r)
{
    int r = f->rperm[rank_r];
    int c = f->cperm[rank_r];
    int rval = 0;

    rval = mpq_move_pivot_row (f, r, c);
    ILL_CLEANUP_IF (rval);

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_move_pivot_row corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

    rval = mpq_move_pivot_col (f, c, r);
    ILL_CLEANUP_IF (rval);

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_move_pivot_col corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

CLEANUP:
    ILL_RESULT (rval, "mpq_move_pivot");	/* Bico 031209 - chg from
						   RETURN */
}

int mpq_ILLfactor_update (mpq_factor_work * f,
      mpq_svector * a,
      int col_p,
      int *p_refact)
{
    int row_p;
    int rank_r = 0;
    int rank_p = 0;
    int rval = 0;
    int nzcnt;
    int *aindx;
    mpq_t *acoef;
    mpq_t *work_coef = f->work_coef;
#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
    mpq_t tmpsize;
#endif
#endif
    int i;

#ifdef mpq_RECORD
    {
	fprintf (mpq_fsave, "u %d %d", col_p, a->nzcnt);
	for (i = 0; i < a->nzcnt; i++) {
	    fprintf (mpq_fsave, " %d %.16e", a->indx[i], mpq_EGlpNumToLf (a->coef[i]));
	}
	fprintf (mpq_fsave, "\n");
	fflush (mpq_fsave);
    }
#endif				/* mpq_RECORD */

#ifdef mpq_DEBUG_FACTOR
    {
	printf ("mpq_ILLfactor_update col %d:", col_p);
	for (i = 0; i < a->nzcnt; i++) {
	    printf (" %.3f*%d", mpq_EGlpNumToLf (a->coef[i]), a->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

#ifdef mpq_DEBUG_FACTOR
    if (mpq_check_matrix (f)) {
	printf ("mpq_ILLfactor_update received corrupted matrix\n");
	fflush (stdout);
    }
#endif				/* mpq_DEBUG_FACTOR */

    if (f->etacnt >= f->etamax) {
	*p_refact = 1;
	return 0;
    }
#ifdef mpq_UPDATE_STUDY
    mpq_nupdate++;
#endif

    row_p = f->ucindx[f->uc_inf[col_p].cbeg];

    mpq_delete_column (f, col_p);

    rval = mpq_create_column (f, a, col_p, &rank_r);
    /* if (rval) fprintf (stderr, "mpq_create_column failed\n"); */
    ILL_CLEANUP_IF (rval);

    rank_p = f->crank[col_p];
#ifdef mpq_UPDATE_STUDY
    if (rank_p != f->rrank[row_p] || rank_p != mpq_column_rank (f, col_p)) {
	printf ("rank_p %d rrank[row_p] %d mpq_column_rank(f,col_p) %d\n",
	    rank_p, f->rrank[row_p], mpq_column_rank (f, col_p));
    }
    if (rank_r > rank_p) {
	mpq_permshifttot += rank_r - rank_p;
    }
    for (i = 0; i < a->nzcnt; i++) {
	if (f->rrank[a->indx[i]] > rank_p)
	    mpq_colspiketot++;
    }
    for (i = 0; i < f->ur_inf[row_p].nzcnt; i++) {
	if (f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] <= rank_r &&
	    f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] != rank_p) {
	    mpq_rowspiketot++;
	}
    }
#endif

    mpq_shift_permutations (f, rank_p, rank_r);

    mpq_delete_row (f, row_p, &f->xtmp);

    f->er_inf[f->etacnt].rbeg = f->er_freebeg;
    f->er_inf[f->etacnt].r = row_p;

    if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim) {
	nzcnt = f->xtmp.nzcnt;
	aindx = f->xtmp.indx;
	acoef = f->xtmp.coef;

	for (i = 0; i < nzcnt; i++) {
	    mpq_EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
	}

	rval = mpq_eliminate_row (f, rank_p, rank_r);
	/* if (rval) fprintf (stderr, "mpq_eliminate_row failed\n"); */
	ILL_CLEANUP_IF (rval);

	rval = mpq_create_row (f, f->work_coef, row_p, rank_r);
	/* if (rval) fprintf (stderr, "mpq_create_row failed\n"); */
	ILL_CLEANUP_IF (rval);
    } else {
	rval = mpq_sparse_eliminate_row (f, &f->xtmp, row_p, rank_r);
	/* if (rval) fprintf (stderr, "mpq_sparse_eliminate_row failed\n"); */
	ILL_CLEANUP_IF (rval);
    }

    if (f->er_freebeg - f->er_inf[f->etacnt].rbeg > 0) {
	f->er_inf[f->etacnt].nzcnt = f->er_freebeg - f->er_inf[f->etacnt].rbeg;
#ifdef mpq_TRACK_FACTOR
	f->nzcnt_cur += f->er_inf[f->etacnt].nzcnt;
#endif				/* mpq_TRACK_FACTOR */
#ifdef mpq_UPDATE_STUDY
	mpq_leftetatot += f->er_inf[f->etacnt].nzcnt;
#endif

#ifdef mpq_SORT_RESULTS
	mpq_sort_vector2 (f->er_inf[f->etacnt].nzcnt,
	    f->erindx + f->er_inf[f->etacnt].rbeg,
	    f->ercoef + f->er_inf[f->etacnt].rbeg);
#endif

	f->etacnt++;
    }
    rval = mpq_move_pivot (f, rank_r);
    /* if (rval) fprintf (stderr, "mpq_move_pivot failed\n"); */
    ILL_CLEANUP_IF (rval);

#ifdef mpq_UPDATE_DEBUG
    printf ("Updated factorization:\n");
#if (mpq_UPDATE_DEBUG+0>1)
    mpq_dump_matrix (f, 0);
#endif
    fflush (stdout);
#endif				/* mpq_UPDATE_DEBUG */

#ifdef mpq_TRACK_FACTOR
#ifdef mpq_NOTICE_BLOWUP
    mpq_EGlpNumInitVar (tmpsize);
    mpq_EGlpNumSet (tmpsize, f->updmaxmult);
    mpq_EGlpNumMultTo (tmpsize, f->maxelem_orig);
    if (mpq_EGlpNumIsLess (tmpsize, f->maxelem_cur)) {
	/* Bico - comment out for dist fprintf (stderr, "factor_update blowup
	   max cur %e max orig %e\n", f->maxelem_cur, f->maxelem_orig); */
	mpq_EGlpNumClearVar (tmpsize);
	return E_FACTOR_BLOWUP;
    }
    mpq_EGlpNumClearVar (tmpsize);
#endif				/* mpq_NOTICE_BLOWUP */
#endif
#ifdef mpq_UPDATE_STATS
    mpq_dump_factor_stats (f);
#endif
CLEANUP:
    ILL_RESULT (rval, "mpq_ILLfactor_update");	/* Bico 031208 - chg from
						   RETURN */
}
