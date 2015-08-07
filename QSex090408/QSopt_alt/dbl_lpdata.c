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

/* RCS_INFO = "$RCSfile: dbl_lpdata.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
/****************************************************************************/
/* */
/* Routines for Manipulating and Writing LPdata               */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int ILLlpdata_buildrows (dbl_ILLlpdata *lp, int **rowbeg, int **rowcnt,     */
/* int **rowind, double **rowval, int include_logicals)          */
/* - include_logicals:  if nonzero, then logical variables will be     */
/* included in the row data                                        */
/* */
/* */
/* All _init routines initialize fields of allocated structure to        */
/* appropiate default values                                             */
/* The _free routines free structures contained in pareameter structure  */
/* but not the parameter itself.                                         */
/* The _alloc routines check whether given parameter is NULL; they either */
/* print an error message or fill structure with default values or the   */
/* given paremeter values.                                               */
/* */
/* void dbl_ILLlpdata_init (dbl_ILLlpdata *lp)                                   */
/* void dbl_ILLlpdata_free (dbl_ILLlpdata *lp)                                   */
/* */
/* void dbl_ILLlp_basis_init (dbl_ILLlp_basis *B)                                */
/* void dbl_ILLlp_basis_free (dbl_ILLlp_basis *B)                                */
/* int dbl_ILLlp_basis_alloc (dbl_ILLlp_basis *B, int nstruct, int nrows)        */
/* */
/* void dbl_ILLlp_cache_init (dbl_ILLlp_cache *C)                                */
/* void dbl_ILLlp_cache_free (dbl_ILLlp_cache *C)                                */
/* int dbl_ILLlp_cache_alloc (dbl_ILLlp_cache *C, int nstruct, int nrows)        */
/* */
/* void dbl_ILLlp_sinfo_init (dbl_ILLlp_sinfo *sinfo)                            */
/* void dbl_ILLlp_sinfo_free (dbl_ILLlp_sinfo *sinfo)                            */
/* */
/* int dbl_ILLlp_rows_init(dbl_ILLlp_rows *lprows, dbl_ILLlpdata *lp,                */
/* int include_logicals)          */
/* */
/****************************************************************************/

#include "econfig.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdata.h"
#include "dbl_qstruct.h"
#include "dbl_qsopt.h"
#include "dbl_lp.h"
#include "dbl_mps.h"
#include "dbl_rawlp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

double dbl_PARAM_IBASIS_RPIVOT;
double dbl_PARAM_IBASIS_RTRIANG;
double dbl_PARAM_MIN_DNORM;
double dbl_PFEAS_TOLER;
double dbl_BD_TOLER;
double dbl_DFEAS_TOLER;
double dbl_PIVOT_TOLER;
double dbl_SZERO_TOLER;
double dbl_PIVZ_TOLER;
double dbl_OBJBND_TOLER;
double dbl_DBNDPIV_TOLER;
double dbl_DBNDPIV_RATIO;
double dbl_PNSTEP_TOLER;
double dbl_DNSTEP_TOLER;
double dbl_ALTPIV_TOLER;
double dbl_DJZERO_TOLER;
double dbl_PROGRESS_ZERO;	/* 1e-7 */
double dbl_PROGRESS_THRESH;	/* 1e-5 */
double dbl_CB_EPS;
double dbl_CB_INF_RATIO;
double dbl_CB_PRI_RLIMIT;
double dbl_ILL_MAXDOUBLE;
double dbl_ILL_MINDOUBLE;

void dbl_ILLstart (void) __attribute__ ((constructor));
void dbl_ILLstart (void)
{
#if VERBOSE_LEVEL <= DEBUG
    char *strtmp = 0;
#endif
    if (!__EGlpNum_setup)
	EGlpNumStart ();
    dbl_EGlpNumInitVar (dbl_PARAM_IBASIS_RPIVOT);
    dbl_EGlpNumInitVar (dbl_PARAM_IBASIS_RTRIANG);
    dbl_EGlpNumInitVar (dbl_PARAM_MIN_DNORM);
    dbl_EGlpNumInitVar (dbl_PFEAS_TOLER);
    dbl_EGlpNumInitVar (dbl_BD_TOLER);
    dbl_EGlpNumInitVar (dbl_DFEAS_TOLER);
    dbl_EGlpNumInitVar (dbl_PIVOT_TOLER);
    dbl_EGlpNumInitVar (dbl_SZERO_TOLER);
    dbl_EGlpNumInitVar (dbl_PIVZ_TOLER);
    dbl_EGlpNumInitVar (dbl_OBJBND_TOLER);
    dbl_EGlpNumInitVar (dbl_DBNDPIV_TOLER);
    dbl_EGlpNumInitVar (dbl_DBNDPIV_RATIO);
    dbl_EGlpNumInitVar (dbl_PNSTEP_TOLER);
    dbl_EGlpNumInitVar (dbl_DNSTEP_TOLER);
    dbl_EGlpNumInitVar (dbl_ALTPIV_TOLER);
    dbl_EGlpNumInitVar (dbl_DJZERO_TOLER);
    dbl_EGlpNumInitVar (dbl_PROGRESS_ZERO);	/* 1e-7 */
    dbl_EGlpNumInitVar (dbl_PROGRESS_THRESH);	/* 1e-5 */
    dbl_EGlpNumInitVar (dbl_CB_PRI_RLIMIT);
    dbl_EGlpNumInitVar (dbl_CB_INF_RATIO);
    dbl_EGlpNumInitVar (dbl_CB_EPS);
    dbl_EGlpNumInitVar (dbl_ILL_MAXDOUBLE);
    dbl_EGlpNumInitVar (dbl_ILL_MINDOUBLE);
    /* parameters that do depend on the tolerance to zero */
    dbl_EGlpNumSet (dbl_PARAM_MIN_DNORM, 5e-9);
    dbl_EGlpNumMultTo (dbl_PARAM_MIN_DNORM, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PFEAS_TOLER, 5e9);
    dbl_EGlpNumMultTo (dbl_PFEAS_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_BD_TOLER, 5e10);
    dbl_EGlpNumMultTo (dbl_BD_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_DFEAS_TOLER, 5e9);
    dbl_EGlpNumMultTo (dbl_DFEAS_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PIVOT_TOLER, 5e5);
    dbl_EGlpNumMultTo (dbl_PIVOT_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_SZERO_TOLER, 5.0);
    dbl_EGlpNumMultTo (dbl_SZERO_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PIVZ_TOLER, 5e3);
    dbl_EGlpNumMultTo (dbl_PIVZ_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_OBJBND_TOLER, 5e13);
    dbl_EGlpNumMultTo (dbl_OBJBND_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_ALTPIV_TOLER, 5e7);
    dbl_EGlpNumMultTo (dbl_ALTPIV_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PROGRESS_ZERO, 5e8);
    dbl_EGlpNumMultTo (dbl_PROGRESS_ZERO, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PROGRESS_THRESH, 5e10);
    dbl_EGlpNumMultTo (dbl_PROGRESS_THRESH, dbl_epsLpNum);
#if VERBOSE_LEVEL <= DEBUG
    strtmp = dbl_EGlpNumGetStr (dbl_PARAM_MIN_DNORM);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_PARAM_MIN_DNORM to %lg %s", dbl_EGlpNumToLf (dbl_PARAM_MIN_DNORM), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_PFEAS_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_PFEAS_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_PFEAS_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_BD_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_BD_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_BD_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_DFEAS_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_DFEAS_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_DFEAS_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_PIVOT_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_PIVOT_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_PIVOT_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_SZERO_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_SZERO_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_SZERO_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_PIVZ_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_PIVZ_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_PIVZ_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_OBJBND_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_OBJBND_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_OBJBND_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_ALTPIV_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_ALTPIV_TOLER to %lg %s", dbl_EGlpNumToLf (dbl_ALTPIV_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_PROGRESS_ZERO);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_PROGRESS_ZERO to %lg %s", dbl_EGlpNumToLf (dbl_PROGRESS_ZERO), strtmp);
    EGfree (strtmp);
    strtmp = dbl_EGlpNumGetStr (dbl_PROGRESS_THRESH);
    MESSAGE (VERBOSE_LEVEL, "Setting dbl_PROGRESS_THRESH to %lg %s", dbl_EGlpNumToLf (dbl_PROGRESS_THRESH), strtmp);
    EGfree (strtmp);
#endif
    /* parameters that do not depend on the tolerance to zero */
    dbl_EGlpNumSet (dbl_ILL_MAXDOUBLE, 1e150);
    dbl_EGlpNumSet (dbl_ILL_MINDOUBLE, -1e150);
    dbl_EGlpNumSet (dbl_PARAM_IBASIS_RPIVOT, 0.98);
    dbl_EGlpNumSet (dbl_PARAM_IBASIS_RTRIANG, 0.01);
    dbl_EGlpNumSet (dbl_DBNDPIV_TOLER, 1e-3);
    dbl_EGlpNumSet (dbl_DBNDPIV_RATIO, 1e-2);
    dbl_EGlpNumSet (dbl_PNSTEP_TOLER, 1e-9);
    dbl_EGlpNumSet (dbl_DNSTEP_TOLER, 1e-9);
    dbl_EGlpNumSet (dbl_DJZERO_TOLER, 1e-8);
    dbl_EGlpNumSet (dbl_CB_EPS, 0.001);
    dbl_EGlpNumSet (dbl_CB_INF_RATIO, 10.0);
    dbl_EGlpNumSet (dbl_CB_PRI_RLIMIT, 0.25);
}

void dbl_ILLchange_precision (void)
{
    dbl_EGlpNumClearVar (dbl_PFEAS_TOLER);
    dbl_EGlpNumClearVar (dbl_BD_TOLER);
    dbl_EGlpNumClearVar (dbl_DFEAS_TOLER);
    dbl_EGlpNumClearVar (dbl_PIVOT_TOLER);
    dbl_EGlpNumClearVar (dbl_SZERO_TOLER);
    dbl_EGlpNumClearVar (dbl_PIVZ_TOLER);
    dbl_EGlpNumClearVar (dbl_OBJBND_TOLER);
    dbl_EGlpNumClearVar (dbl_ALTPIV_TOLER);
    dbl_EGlpNumClearVar (dbl_PARAM_MIN_DNORM);
    dbl_EGlpNumClearVar (dbl_PROGRESS_ZERO);
    dbl_EGlpNumClearVar (dbl_PROGRESS_THRESH);
    dbl_EGlpNumInitVar (dbl_PROGRESS_ZERO);
    dbl_EGlpNumInitVar (dbl_PROGRESS_THRESH);
    dbl_EGlpNumInitVar (dbl_PFEAS_TOLER);
    dbl_EGlpNumInitVar (dbl_BD_TOLER);
    dbl_EGlpNumInitVar (dbl_DFEAS_TOLER);
    dbl_EGlpNumInitVar (dbl_PIVOT_TOLER);
    dbl_EGlpNumInitVar (dbl_SZERO_TOLER);
    dbl_EGlpNumInitVar (dbl_PIVZ_TOLER);
    dbl_EGlpNumInitVar (dbl_OBJBND_TOLER);
    dbl_EGlpNumInitVar (dbl_ALTPIV_TOLER);
    dbl_EGlpNumInitVar (dbl_PARAM_MIN_DNORM);
    /* parameters that do depend on the tolerance to zero */
    dbl_EGlpNumSet (dbl_PARAM_MIN_DNORM, 5e-9);
    dbl_EGlpNumMultTo (dbl_PARAM_MIN_DNORM, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PFEAS_TOLER, 5e9);
    dbl_EGlpNumMultTo (dbl_PFEAS_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_BD_TOLER, 5e10);
    dbl_EGlpNumMultTo (dbl_BD_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_DFEAS_TOLER, 5e9);
    dbl_EGlpNumMultTo (dbl_DFEAS_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PIVOT_TOLER, 5e5);
    dbl_EGlpNumMultTo (dbl_PIVOT_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_SZERO_TOLER, 5.0);
    dbl_EGlpNumMultTo (dbl_SZERO_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PIVZ_TOLER, 5e3);
    dbl_EGlpNumMultTo (dbl_PIVZ_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_OBJBND_TOLER, 5e13);
    dbl_EGlpNumMultTo (dbl_OBJBND_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_ALTPIV_TOLER, 5e7);
    dbl_EGlpNumMultTo (dbl_ALTPIV_TOLER, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PROGRESS_ZERO, 5e8);
    dbl_EGlpNumMultTo (dbl_PROGRESS_ZERO, dbl_epsLpNum);
    dbl_EGlpNumSet (dbl_PROGRESS_THRESH, 5e10);
    dbl_EGlpNumMultTo (dbl_PROGRESS_THRESH, dbl_epsLpNum);
}

void dbl_ILLend (void) __attribute__ ((destructor));
void dbl_ILLend (void)
{
    dbl_EGlpNumClearVar (dbl_PARAM_IBASIS_RPIVOT);
    dbl_EGlpNumClearVar (dbl_PARAM_IBASIS_RTRIANG);
    dbl_EGlpNumClearVar (dbl_PARAM_MIN_DNORM);
    dbl_EGlpNumClearVar (dbl_PFEAS_TOLER);
    dbl_EGlpNumClearVar (dbl_BD_TOLER);
    dbl_EGlpNumClearVar (dbl_DFEAS_TOLER);
    dbl_EGlpNumClearVar (dbl_PIVOT_TOLER);
    dbl_EGlpNumClearVar (dbl_SZERO_TOLER);
    dbl_EGlpNumClearVar (dbl_PIVZ_TOLER);
    dbl_EGlpNumClearVar (dbl_OBJBND_TOLER);
    dbl_EGlpNumClearVar (dbl_DBNDPIV_TOLER);
    dbl_EGlpNumClearVar (dbl_DBNDPIV_RATIO);
    dbl_EGlpNumClearVar (dbl_PNSTEP_TOLER);
    dbl_EGlpNumClearVar (dbl_DNSTEP_TOLER);
    dbl_EGlpNumClearVar (dbl_ALTPIV_TOLER);
    dbl_EGlpNumClearVar (dbl_DJZERO_TOLER);
    dbl_EGlpNumClearVar (dbl_PROGRESS_ZERO);	/* 1e-7 */
    dbl_EGlpNumClearVar (dbl_PROGRESS_THRESH);	/* 1e-5 */
    dbl_EGlpNumClearVar (dbl_CB_EPS);
    dbl_EGlpNumClearVar (dbl_CB_INF_RATIO);
    dbl_EGlpNumClearVar (dbl_CB_PRI_RLIMIT);
    dbl_EGlpNumClearVar (dbl_ILL_MAXDOUBLE);
    dbl_EGlpNumClearVar (dbl_ILL_MINDOUBLE);
}

dbl_QSdata *dbl_ILLread (dbl_qsline_reader * file,
      const char *dbl_fname,
      int isMps)
{
    int rval = 0;
    dbl_QSdata *p = 0;
    dbl_ILLlpdata *lp;
    dbl_rawlpdata rawlp;

    ILL_FAILfalse (file != NULL, NULL);
    ILL_FAILfalse (dbl_fname != NULL, NULL);

    p = dbl_QScreate_prob (dbl_fname, QS_MIN);
    ILL_CHECKnull (p, NULL);
    ILL_IFFREE (p->qslp->probname, char);
    lp = p->qslp;

    dbl_ILLinit_rawlpdata (&rawlp, file->error_collector);
    dbl_ILLlpdata_init (lp);

    if (isMps != 0) {
	rval = dbl_ILLread_mps (file, dbl_fname, &rawlp);
    } else {
	rval = dbl_ILLread_lp (file, dbl_fname, &rawlp);
    }
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLrawlpdata_to_lpdata (&rawlp, lp);
    ILL_CLEANUP_IF (rval);

CLEANUP:
    dbl_ILLfree_rawlpdata (&rawlp);
    if (rval != 0) {
	dbl_QSfree_prob (p);
	p = 0;
    }
    return p;
}

void dbl_ILLlpdata_init (dbl_ILLlpdata * lp)
{
    if (lp) {
	lp->nrows = 0;
	lp->ncols = 0;
	lp->nstruct = 0;
	lp->nzcount = 0;
	lp->rowsize = 0;
	lp->colsize = 0;
	lp->structsize = 0;
	lp->objsense = dbl_ILL_MIN;
	lp->sense = 0;
	lp->obj = 0;
	lp->rhs = 0;
	lp->rangeval = 0;
	lp->lower = 0;
	lp->upper = 0;

	dbl_ILLmatrix_init (&lp->A);
	dbl_ILLmatrix_init (&lp->sos);
	lp->rA = 0;
	lp->is_sos_mem = NULL;
	lp->refrowname = NULL;
	lp->refind = -1;

	lp->colnames = 0;
	ILLsymboltab_init (&lp->coltab);
	lp->rownames = 0;
	ILLsymboltab_init (&lp->rowtab);
	lp->objname = 0;

	lp->probname = 0;
	lp->intmarker = 0;
	lp->structmap = 0;
	lp->rowmap = 0;
	lp->basis = 0;
	/* lp->presolve   = 0; */
	lp->sinfo = 0;

	ILLstring_reporter_init (&lp->reporter, ILL_fprintf, stdout);
    }
}

void dbl_ILLlpdata_free (dbl_ILLlpdata * lp)
{
    int i;

    if (lp) {
	ILL_IFFREE (lp->sense, char);
	dbl_EGlpNumFreeArray (lp->obj);
	dbl_EGlpNumFreeArray (lp->rhs);
	dbl_EGlpNumFreeArray (lp->rangeval);
	dbl_EGlpNumFreeArray (lp->lower);
	dbl_EGlpNumFreeArray (lp->upper);
	dbl_ILLmatrix_free (&lp->A);
	if (lp->rA) {
	    dbl_ILLlp_rows_clear (lp->rA);
	    ILL_IFFREE (lp->rA, dbl_ILLlp_rows);
	}
	ILL_IFFREE (lp->is_sos_mem, int);
	ILL_IFFREE (lp->refrowname, char);
	dbl_ILLmatrix_free (&lp->sos);
	if (lp->colnames) {
	    for (i = 0; i < lp->nstruct; i++) {
		ILL_IFFREE (lp->colnames[i], char);
	    }
	    ILL_IFFREE (lp->colnames, char *);
	}
	ILLsymboltab_free (&lp->coltab);
	if (lp->rownames) {
	    for (i = 0; i < lp->nrows; i++) {
		ILL_IFFREE (lp->rownames[i], char);
	    }
	    ILL_IFFREE (lp->rownames, char *);
	}
	ILLsymboltab_free (&lp->rowtab);
	ILL_IFFREE (lp->objname, char);
	ILL_IFFREE (lp->probname, char);
	ILL_IFFREE (lp->intmarker, char);
	ILL_IFFREE (lp->structmap, int);
	ILL_IFFREE (lp->rowmap, int);
	if (lp->sinfo) {
	    dbl_ILLlp_sinfo_free (lp->sinfo);
	    ILL_IFFREE (lp->sinfo, dbl_ILLlp_sinfo);
	}
	dbl_ILLlpdata_init (lp);
    }
}

void dbl_ILLlp_basis_init (dbl_ILLlp_basis * B)
{
    if (B) {
	B->cstat = 0;
	B->rstat = 0;
	B->rownorms = 0;
	B->colnorms = 0;
	B->nstruct = 0;
	B->nrows = 0;
    }
}

void dbl_ILLlp_basis_free (dbl_ILLlp_basis * B)
{
    if (B) {
	ILL_IFFREE (B->cstat, char);
	ILL_IFFREE (B->rstat, char);
	dbl_EGlpNumFreeArray (B->rownorms);
	dbl_EGlpNumFreeArray (B->colnorms);
	B->nstruct = 0;
	B->nrows = 0;
    }
}

int dbl_ILLlp_basis_alloc (dbl_ILLlp_basis * B,
      int nstruct,
      int nrows)
{
    int rval = 0;

    ILL_FAILtrue (B == NULL, "dbl_ILLlp_basis_alloc called without a basis");

    B->nstruct = nstruct;
    B->nrows = nrows;

    if (nstruct > 0) {
	ILL_SAFE_MALLOC (B->cstat, nstruct, char);
    }
    if (nrows > 0) {
	ILL_SAFE_MALLOC (B->rstat, nrows, char);
    }
CLEANUP:

    if (rval) {
	dbl_ILLlp_basis_free (B);
    }
    ILL_RETURN (rval, "dbl_ILLlp_basis_alloc");
}

void dbl_ILLlp_cache_init (dbl_ILLlp_cache * C)
{
    if (C) {
	C->x = 0;
	C->rc = 0;
	C->pi = 0;
	C->slack = 0;
	C->nstruct = 0;
	C->nrows = 0;
	C->status = 0;
	dbl_EGlpNumZero (C->val);
    }
}

void dbl_ILLlp_cache_free (dbl_ILLlp_cache * C)
{
    if (C) {
	dbl_EGlpNumFreeArray (C->x);
	dbl_EGlpNumFreeArray (C->rc);
	dbl_EGlpNumFreeArray (C->pi);
	dbl_EGlpNumFreeArray (C->slack);
	C->nstruct = 0;
	C->nrows = 0;
	C->status = 0;
    }
}

int dbl_ILLlp_cache_alloc (dbl_ILLlp_cache * C,
      int nstruct,
      int nrows)
{
    int rval = 0;

    ILL_FAILtrue (C == NULL, "dbl_ILLlp_cache_alloc called without a cache");

    C->nstruct = nstruct;
    C->nrows = nrows;

    if (nstruct > 0) {
	C->x = dbl_EGlpNumAllocArray (nstruct);
	C->rc = dbl_EGlpNumAllocArray (nstruct);
    }
    if (nrows > 0) {
	C->pi = dbl_EGlpNumAllocArray (nrows);
	C->slack = dbl_EGlpNumAllocArray (nrows);
    }
CLEANUP:

    if (rval) {
	dbl_ILLlp_cache_free (C);
    }
    ILL_RETURN (rval, "dbl_ILLlp_cache_alloc");
}


int dbl_ILLlp_rows_init (dbl_ILLlp_rows * lprows,
      dbl_ILLlpdata * lp,
      int include_logicals)
{
    int rval = 0;
    int i, k, st;
    int *beg, *cnt, *ind;
    double *val;
    dbl_ILLmatrix *A;
    char *hit = 0;
    int *inv_structmap = 0;

    /* If logicals are not included, then the columns are ordered as in */
    /* lp->structmap.  Otherwise, the columns are ordered as in the     */
    /* matrix structure.                                                */

    if (lprows != NULL) {
	lprows->rowbeg = 0;
	lprows->rowcnt = 0;
	lprows->rowind = 0;
	lprows->rowval = 0;
    }
    ILL_FAILfalse ((lp != NULL) && (lprows != NULL),
	"called with a NULL pointer");

    A = &lp->A;

    if (lp->nrows > 0) {
	if (include_logicals == 0) {
	    ILL_FAILtrue (lp->rowmap == NULL, "Programming error.");
	    ILL_SAFE_MALLOC (hit, lp->ncols, char);

	    for (i = 0; i < lp->ncols; i++) {
		hit[i] = 0;
	    }
	    for (i = 0; i < lp->nrows; i++) {
		hit[lp->rowmap[i]] = 1;
	    }

	    ILL_SAFE_MALLOC (inv_structmap, lp->ncols, int);

	    for (i = 0; i < lp->nstruct; i++) {
		inv_structmap[lp->structmap[i]] = i;
	    }
	}
	ILL_SAFE_MALLOC (lprows->rowbeg, lp->nrows, int);
	ILL_SAFE_MALLOC (lprows->rowcnt, lp->nrows, int);

	if (((include_logicals != 0) && lp->nzcount > 0) ||
	    ((include_logicals == 0) && lp->nzcount > lp->nrows)) {
	    if (include_logicals != 0) {
		ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount, int);
		lprows->rowval = dbl_EGlpNumAllocArray (lp->nzcount);
	    } else {
		ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount - lp->nrows, int);
		lprows->rowval = dbl_EGlpNumAllocArray (lp->nzcount - lp->nrows);
	    }
	}
	beg = lprows->rowbeg;
	cnt = lprows->rowcnt;
	ind = lprows->rowind;
	val = lprows->rowval;

	for (i = 0; i < lp->nrows; i++) {
	    cnt[i] = 0;
	}

	for (i = 0; i < lp->ncols; i++) {
	    if ((include_logicals != 0) || hit[i] == 0) {
		k = A->matbeg[i];
		st = k + A->matcnt[i];
		for (; k < st; k++) {
		    cnt[A->matind[k]]++;
		}
	    }
	}

	for (i = 0, k = 0; i < lp->nrows; i++) {
	    beg[i] = k;
	    k += cnt[i];
	}

	for (i = 0; i < lp->ncols; i++) {
	    if ((include_logicals != 0) || hit[i] == 0) {
		k = A->matbeg[i];
		st = k + A->matcnt[i];
		for (; k < st; k++) {
		    if (include_logicals != 0) {
			ind[beg[A->matind[k]]] = i;
		    } else {
			ind[beg[A->matind[k]]] = inv_structmap[i];
		    }
		    dbl_EGlpNumCopy (val[beg[A->matind[k]]], A->matval[k]);
		    beg[A->matind[k]]++;
		}
	    }
	}

	for (i = 0, k = 0; i < lp->nrows; i++) {
	    beg[i] = k;
	    k += cnt[i];
	}
    }
CLEANUP:

    if (rval) {
	dbl_ILLlp_rows_clear (lprows);
    }
    ILL_IFFREE (hit, char);
    ILL_IFFREE (inv_structmap, int);

    ILL_RETURN (rval, "dbl_ILLlp_rows_init");
}

void dbl_ILLlp_rows_clear (dbl_ILLlp_rows * lprows)
{
    if (lprows != NULL) {
	ILL_IFFREE (lprows->rowbeg, int);
	ILL_IFFREE (lprows->rowcnt, int);
	ILL_IFFREE (lprows->rowind, int);
	dbl_EGlpNumFreeArray (lprows->rowval);
    }
}

static int dbl_wr_line (dbl_ILLlpdata * lp,
      const char *format,
      va_list argptr)
{
    char buffer[ILL_namebufsize];
    int rval = 0;
    rval = vsprintf (buffer, format, argptr);
    if (rval > 0) {
	rval = fprintf (lp->reporter.dest, buffer);
	rval = (rval < 0) ? 1 : 0;
	/* rval = ILLstring_report (buffer, &lp->reporter); */
    }
    return rval;
}

int dbl_ILLprint_report (dbl_ILLlpdata * lp,
      const char *format,
    ...)
{
    va_list marker;
    int rval = 0;

    va_start (marker, format);	/* ANSI style */
    rval = dbl_wr_line (lp, format, marker);
    va_end (marker);		/* Reset variable arguments.      */
    return rval;
}
