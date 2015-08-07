/****************************************************************************/
/* */
/* This file is part of QSopt_ex.                                          */
/* */
/* (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/* and Daniel Espinoza                                                     */
/* */
/* This code may be used under the terms of the GNU General Public License */
/* (Version 2.1 or later) as published by the Free Software Foundation.    */
/* */
/* Alternatively, use is granted for research purposes only.               */
/* */
/* It is your choice of which of these two licenses you are operating      */
/* under.                                                                  */
/* */
/* We make no guarantees about the correctness or usefulness of this code. */
/* */
/****************************************************************************/

/* RCS_INFO = "$RCSfile: mpf_lpdata.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
/****************************************************************************/
/* */
/* Routines for Manipulating and Writing LPdata               */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int ILLlpdata_buildrows (mpf_ILLlpdata *lp, int **rowbeg, int **rowcnt,     */
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
/* void mpf_ILLlpdata_init (mpf_ILLlpdata *lp)                                   */
/* void mpf_ILLlpdata_free (mpf_ILLlpdata *lp)                                   */
/* */
/* void mpf_ILLlp_basis_init (mpf_ILLlp_basis *B)                                */
/* void mpf_ILLlp_basis_free (mpf_ILLlp_basis *B)                                */
/* int mpf_ILLlp_basis_alloc (mpf_ILLlp_basis *B, int nstruct, int nrows)        */
/* */
/* void mpf_ILLlp_cache_init (mpf_ILLlp_cache *C)                                */
/* void mpf_ILLlp_cache_free (mpf_ILLlp_cache *C)                                */
/* int mpf_ILLlp_cache_alloc (mpf_ILLlp_cache *C, int nstruct, int nrows)        */
/* */
/* void mpf_ILLlp_sinfo_init (mpf_ILLlp_sinfo *sinfo)                            */
/* void mpf_ILLlp_sinfo_free (mpf_ILLlp_sinfo *sinfo)                            */
/* */
/* int mpf_ILLlp_rows_init(mpf_ILLlp_rows *lprows, mpf_ILLlpdata *lp,                */
/* int include_logicals)          */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdata.h"
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#include "mpf_lp.h"
#include "mpf_mps.h"
#include "mpf_rawlp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

mpf_t mpf_PARAM_IBASIS_RPIVOT;
mpf_t mpf_PARAM_IBASIS_RTRIANG;
mpf_t mpf_PARAM_MIN_DNORM;
mpf_t mpf_PFEAS_TOLER;
mpf_t mpf_BD_TOLER;
mpf_t mpf_DFEAS_TOLER;
mpf_t mpf_PIVOT_TOLER;
mpf_t mpf_SZERO_TOLER;
mpf_t mpf_PIVZ_TOLER;
mpf_t mpf_OBJBND_TOLER;
mpf_t mpf_DBNDPIV_TOLER;
mpf_t mpf_DBNDPIV_RATIO;
mpf_t mpf_PNSTEP_TOLER;
mpf_t mpf_DNSTEP_TOLER;
mpf_t mpf_ALTPIV_TOLER;
mpf_t mpf_DJZERO_TOLER;
mpf_t mpf_PROGRESS_ZERO;	/* 1e-7 */
mpf_t mpf_PROGRESS_THRESH;	/* 1e-5 */
mpf_t mpf_CB_EPS;
mpf_t mpf_CB_INF_RATIO;
mpf_t mpf_CB_PRI_RLIMIT;
mpf_t mpf_ILL_MAXDOUBLE;
mpf_t mpf_ILL_MINDOUBLE;

void mpf_ILLstart (void) __attribute__ ((constructor));
void mpf_ILLstart (void)
{
#if VERBOSE_LEVEL <= DEBUG
    char *strtmp = 0;
#endif
    if (!__EGlpNum_setup)
	EGlpNumStart ();
    mpf_EGlpNumInitVar (mpf_PARAM_IBASIS_RPIVOT);
    mpf_EGlpNumInitVar (mpf_PARAM_IBASIS_RTRIANG);
    mpf_EGlpNumInitVar (mpf_PARAM_MIN_DNORM);
    mpf_EGlpNumInitVar (mpf_PFEAS_TOLER);
    mpf_EGlpNumInitVar (mpf_BD_TOLER);
    mpf_EGlpNumInitVar (mpf_DFEAS_TOLER);
    mpf_EGlpNumInitVar (mpf_PIVOT_TOLER);
    mpf_EGlpNumInitVar (mpf_SZERO_TOLER);
    mpf_EGlpNumInitVar (mpf_PIVZ_TOLER);
    mpf_EGlpNumInitVar (mpf_OBJBND_TOLER);
    mpf_EGlpNumInitVar (mpf_DBNDPIV_TOLER);
    mpf_EGlpNumInitVar (mpf_DBNDPIV_RATIO);
    mpf_EGlpNumInitVar (mpf_PNSTEP_TOLER);
    mpf_EGlpNumInitVar (mpf_DNSTEP_TOLER);
    mpf_EGlpNumInitVar (mpf_ALTPIV_TOLER);
    mpf_EGlpNumInitVar (mpf_DJZERO_TOLER);
    mpf_EGlpNumInitVar (mpf_PROGRESS_ZERO);	/* 1e-7 */
    mpf_EGlpNumInitVar (mpf_PROGRESS_THRESH);	/* 1e-5 */
    mpf_EGlpNumInitVar (mpf_CB_PRI_RLIMIT);
    mpf_EGlpNumInitVar (mpf_CB_INF_RATIO);
    mpf_EGlpNumInitVar (mpf_CB_EPS);
    mpf_EGlpNumInitVar (mpf_ILL_MAXDOUBLE);
    mpf_EGlpNumInitVar (mpf_ILL_MINDOUBLE);
    /* parameters that do depend on the tolerance to zero */
    mpf_EGlpNumSet (mpf_PARAM_MIN_DNORM, 5e-9);
    mpf_EGlpNumMultTo (mpf_PARAM_MIN_DNORM, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PFEAS_TOLER, 5e9);
    mpf_EGlpNumMultTo (mpf_PFEAS_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_BD_TOLER, 5e10);
    mpf_EGlpNumMultTo (mpf_BD_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_DFEAS_TOLER, 5e9);
    mpf_EGlpNumMultTo (mpf_DFEAS_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PIVOT_TOLER, 5e5);
    mpf_EGlpNumMultTo (mpf_PIVOT_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_SZERO_TOLER, 5.0);
    mpf_EGlpNumMultTo (mpf_SZERO_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PIVZ_TOLER, 5e3);
    mpf_EGlpNumMultTo (mpf_PIVZ_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_OBJBND_TOLER, 5e13);
    mpf_EGlpNumMultTo (mpf_OBJBND_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_ALTPIV_TOLER, 5e7);
    mpf_EGlpNumMultTo (mpf_ALTPIV_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PROGRESS_ZERO, 5e8);
    mpf_EGlpNumMultTo (mpf_PROGRESS_ZERO, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PROGRESS_THRESH, 5e10);
    mpf_EGlpNumMultTo (mpf_PROGRESS_THRESH, mpf_epsLpNum);
#if VERBOSE_LEVEL <= DEBUG
    strtmp = mpf_EGlpNumGetStr (mpf_PARAM_MIN_DNORM);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_PARAM_MIN_DNORM to %lg %s", mpf_EGlpNumToLf (mpf_PARAM_MIN_DNORM), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_PFEAS_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_PFEAS_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_PFEAS_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_BD_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_BD_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_BD_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_DFEAS_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_DFEAS_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_DFEAS_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_PIVOT_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_PIVOT_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_PIVOT_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_SZERO_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_SZERO_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_SZERO_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_PIVZ_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_PIVZ_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_PIVZ_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_OBJBND_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_OBJBND_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_OBJBND_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_ALTPIV_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_ALTPIV_TOLER to %lg %s", mpf_EGlpNumToLf (mpf_ALTPIV_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_PROGRESS_ZERO);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_PROGRESS_ZERO to %lg %s", mpf_EGlpNumToLf (mpf_PROGRESS_ZERO), strtmp);
    EGfree (strtmp);
    strtmp = mpf_EGlpNumGetStr (mpf_PROGRESS_THRESH);
    MESSAGE (VERBOSE_LEVEL, "Setting mpf_PROGRESS_THRESH to %lg %s", mpf_EGlpNumToLf (mpf_PROGRESS_THRESH), strtmp);
    EGfree (strtmp);
#endif
    /* parameters that do not depend on the tolerance to zero */
    mpf_EGlpNumSet (mpf_ILL_MAXDOUBLE, 1e150);
    mpf_EGlpNumSet (mpf_ILL_MINDOUBLE, -1e150);
    mpf_EGlpNumSet (mpf_PARAM_IBASIS_RPIVOT, 0.98);
    mpf_EGlpNumSet (mpf_PARAM_IBASIS_RTRIANG, 0.01);
    mpf_EGlpNumSet (mpf_DBNDPIV_TOLER, 1e-3);
    mpf_EGlpNumSet (mpf_DBNDPIV_RATIO, 1e-2);
    mpf_EGlpNumSet (mpf_PNSTEP_TOLER, 1e-9);
    mpf_EGlpNumSet (mpf_DNSTEP_TOLER, 1e-9);
    mpf_EGlpNumSet (mpf_DJZERO_TOLER, 1e-8);
    mpf_EGlpNumSet (mpf_CB_EPS, 0.001);
    mpf_EGlpNumSet (mpf_CB_INF_RATIO, 10.0);
    mpf_EGlpNumSet (mpf_CB_PRI_RLIMIT, 0.25);
}

void mpf_ILLchange_precision (void)
{
    mpf_EGlpNumClearVar (mpf_PFEAS_TOLER);
    mpf_EGlpNumClearVar (mpf_BD_TOLER);
    mpf_EGlpNumClearVar (mpf_DFEAS_TOLER);
    mpf_EGlpNumClearVar (mpf_PIVOT_TOLER);
    mpf_EGlpNumClearVar (mpf_SZERO_TOLER);
    mpf_EGlpNumClearVar (mpf_PIVZ_TOLER);
    mpf_EGlpNumClearVar (mpf_OBJBND_TOLER);
    mpf_EGlpNumClearVar (mpf_ALTPIV_TOLER);
    mpf_EGlpNumClearVar (mpf_PARAM_MIN_DNORM);
    mpf_EGlpNumClearVar (mpf_PROGRESS_ZERO);
    mpf_EGlpNumClearVar (mpf_PROGRESS_THRESH);
    mpf_EGlpNumInitVar (mpf_PROGRESS_ZERO);
    mpf_EGlpNumInitVar (mpf_PROGRESS_THRESH);
    mpf_EGlpNumInitVar (mpf_PFEAS_TOLER);
    mpf_EGlpNumInitVar (mpf_BD_TOLER);
    mpf_EGlpNumInitVar (mpf_DFEAS_TOLER);
    mpf_EGlpNumInitVar (mpf_PIVOT_TOLER);
    mpf_EGlpNumInitVar (mpf_SZERO_TOLER);
    mpf_EGlpNumInitVar (mpf_PIVZ_TOLER);
    mpf_EGlpNumInitVar (mpf_OBJBND_TOLER);
    mpf_EGlpNumInitVar (mpf_ALTPIV_TOLER);
    mpf_EGlpNumInitVar (mpf_PARAM_MIN_DNORM);
    /* parameters that do depend on the tolerance to zero */
    mpf_EGlpNumSet (mpf_PARAM_MIN_DNORM, 5e-9);
    mpf_EGlpNumMultTo (mpf_PARAM_MIN_DNORM, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PFEAS_TOLER, 5e9);
    mpf_EGlpNumMultTo (mpf_PFEAS_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_BD_TOLER, 5e10);
    mpf_EGlpNumMultTo (mpf_BD_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_DFEAS_TOLER, 5e9);
    mpf_EGlpNumMultTo (mpf_DFEAS_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PIVOT_TOLER, 5e5);
    mpf_EGlpNumMultTo (mpf_PIVOT_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_SZERO_TOLER, 5.0);
    mpf_EGlpNumMultTo (mpf_SZERO_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PIVZ_TOLER, 5e3);
    mpf_EGlpNumMultTo (mpf_PIVZ_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_OBJBND_TOLER, 5e13);
    mpf_EGlpNumMultTo (mpf_OBJBND_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_ALTPIV_TOLER, 5e7);
    mpf_EGlpNumMultTo (mpf_ALTPIV_TOLER, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PROGRESS_ZERO, 5e8);
    mpf_EGlpNumMultTo (mpf_PROGRESS_ZERO, mpf_epsLpNum);
    mpf_EGlpNumSet (mpf_PROGRESS_THRESH, 5e10);
    mpf_EGlpNumMultTo (mpf_PROGRESS_THRESH, mpf_epsLpNum);
}

void mpf_ILLend (void) __attribute__ ((destructor));
void mpf_ILLend (void)
{
    mpf_EGlpNumClearVar (mpf_PARAM_IBASIS_RPIVOT);
    mpf_EGlpNumClearVar (mpf_PARAM_IBASIS_RTRIANG);
    mpf_EGlpNumClearVar (mpf_PARAM_MIN_DNORM);
    mpf_EGlpNumClearVar (mpf_PFEAS_TOLER);
    mpf_EGlpNumClearVar (mpf_BD_TOLER);
    mpf_EGlpNumClearVar (mpf_DFEAS_TOLER);
    mpf_EGlpNumClearVar (mpf_PIVOT_TOLER);
    mpf_EGlpNumClearVar (mpf_SZERO_TOLER);
    mpf_EGlpNumClearVar (mpf_PIVZ_TOLER);
    mpf_EGlpNumClearVar (mpf_OBJBND_TOLER);
    mpf_EGlpNumClearVar (mpf_DBNDPIV_TOLER);
    mpf_EGlpNumClearVar (mpf_DBNDPIV_RATIO);
    mpf_EGlpNumClearVar (mpf_PNSTEP_TOLER);
    mpf_EGlpNumClearVar (mpf_DNSTEP_TOLER);
    mpf_EGlpNumClearVar (mpf_ALTPIV_TOLER);
    mpf_EGlpNumClearVar (mpf_DJZERO_TOLER);
    mpf_EGlpNumClearVar (mpf_PROGRESS_ZERO);	/* 1e-7 */
    mpf_EGlpNumClearVar (mpf_PROGRESS_THRESH);	/* 1e-5 */
    mpf_EGlpNumClearVar (mpf_CB_EPS);
    mpf_EGlpNumClearVar (mpf_CB_INF_RATIO);
    mpf_EGlpNumClearVar (mpf_CB_PRI_RLIMIT);
    mpf_EGlpNumClearVar (mpf_ILL_MAXDOUBLE);
    mpf_EGlpNumClearVar (mpf_ILL_MINDOUBLE);
}

mpf_QSdata *mpf_ILLread (mpf_qsline_reader * file,
      const char *mpf_fname,
      int isMps)
{
    int rval = 0;
    mpf_QSdata *p = 0;
    mpf_ILLlpdata *lp;
    mpf_rawlpdata rawlp;

    ILL_FAILfalse (file != NULL, NULL);
    ILL_FAILfalse (mpf_fname != NULL, NULL);

    p = mpf_QScreate_prob (mpf_fname, QS_MIN);
    ILL_CHECKnull (p, NULL);
    ILL_IFFREE (p->qslp->probname, char);
    lp = p->qslp;

    mpf_ILLinit_rawlpdata (&rawlp, file->error_collector);
    mpf_ILLlpdata_init (lp);

    if (isMps != 0) {
	rval = mpf_ILLread_mps (file, mpf_fname, &rawlp);
    } else {
	rval = mpf_ILLread_lp (file, mpf_fname, &rawlp);
    }
    ILL_CLEANUP_IF (rval);

    rval = mpf_ILLrawlpdata_to_lpdata (&rawlp, lp);
    ILL_CLEANUP_IF (rval);

CLEANUP:
    mpf_ILLfree_rawlpdata (&rawlp);
    if (rval != 0) {
	mpf_QSfree_prob (p);
	p = 0;
    }
    return p;
}

void mpf_ILLlpdata_init (mpf_ILLlpdata * lp)
{
    if (lp) {
	lp->nrows = 0;
	lp->ncols = 0;
	lp->nstruct = 0;
	lp->nzcount = 0;
	lp->rowsize = 0;
	lp->colsize = 0;
	lp->structsize = 0;
	lp->objsense = mpf_ILL_MIN;
	lp->sense = 0;
	lp->obj = 0;
	lp->rhs = 0;
	lp->rangeval = 0;
	lp->lower = 0;
	lp->upper = 0;

	mpf_ILLmatrix_init (&lp->A);
	mpf_ILLmatrix_init (&lp->sos);
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

void mpf_ILLlpdata_free (mpf_ILLlpdata * lp)
{
    int i;

    if (lp) {
	ILL_IFFREE (lp->sense, char);
	mpf_EGlpNumFreeArray (lp->obj);
	mpf_EGlpNumFreeArray (lp->rhs);
	mpf_EGlpNumFreeArray (lp->rangeval);
	mpf_EGlpNumFreeArray (lp->lower);
	mpf_EGlpNumFreeArray (lp->upper);
	mpf_ILLmatrix_free (&lp->A);
	if (lp->rA) {
	    mpf_ILLlp_rows_clear (lp->rA);
	    ILL_IFFREE (lp->rA, mpf_ILLlp_rows);
	}
	ILL_IFFREE (lp->is_sos_mem, int);
	ILL_IFFREE (lp->refrowname, char);
	mpf_ILLmatrix_free (&lp->sos);
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
	    mpf_ILLlp_sinfo_free (lp->sinfo);
	    ILL_IFFREE (lp->sinfo, mpf_ILLlp_sinfo);
	}
	mpf_ILLlpdata_init (lp);
    }
}

void mpf_ILLlp_basis_init (mpf_ILLlp_basis * B)
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

void mpf_ILLlp_basis_free (mpf_ILLlp_basis * B)
{
    if (B) {
	ILL_IFFREE (B->cstat, char);
	ILL_IFFREE (B->rstat, char);
	mpf_EGlpNumFreeArray (B->rownorms);
	mpf_EGlpNumFreeArray (B->colnorms);
	B->nstruct = 0;
	B->nrows = 0;
    }
}

int mpf_ILLlp_basis_alloc (mpf_ILLlp_basis * B,
      int nstruct,
      int nrows)
{
    int rval = 0;

    ILL_FAILtrue (B == NULL, "mpf_ILLlp_basis_alloc called without a basis");

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
	mpf_ILLlp_basis_free (B);
    }
    ILL_RETURN (rval, "mpf_ILLlp_basis_alloc");
}

void mpf_ILLlp_cache_init (mpf_ILLlp_cache * C)
{
    if (C) {
	C->x = 0;
	C->rc = 0;
	C->pi = 0;
	C->slack = 0;
	C->nstruct = 0;
	C->nrows = 0;
	C->status = 0;
	mpf_EGlpNumZero (C->val);
    }
}

void mpf_ILLlp_cache_free (mpf_ILLlp_cache * C)
{
    if (C) {
	mpf_EGlpNumFreeArray (C->x);
	mpf_EGlpNumFreeArray (C->rc);
	mpf_EGlpNumFreeArray (C->pi);
	mpf_EGlpNumFreeArray (C->slack);
	C->nstruct = 0;
	C->nrows = 0;
	C->status = 0;
    }
}

int mpf_ILLlp_cache_alloc (mpf_ILLlp_cache * C,
      int nstruct,
      int nrows)
{
    int rval = 0;

    ILL_FAILtrue (C == NULL, "mpf_ILLlp_cache_alloc called without a cache");

    C->nstruct = nstruct;
    C->nrows = nrows;

    if (nstruct > 0) {
	C->x = mpf_EGlpNumAllocArray (nstruct);
	C->rc = mpf_EGlpNumAllocArray (nstruct);
    }
    if (nrows > 0) {
	C->pi = mpf_EGlpNumAllocArray (nrows);
	C->slack = mpf_EGlpNumAllocArray (nrows);
    }
CLEANUP:

    if (rval) {
	mpf_ILLlp_cache_free (C);
    }
    ILL_RETURN (rval, "mpf_ILLlp_cache_alloc");
}


int mpf_ILLlp_rows_init (mpf_ILLlp_rows * lprows,
      mpf_ILLlpdata * lp,
      int include_logicals)
{
    int rval = 0;
    int i, k, st;
    int *beg, *cnt, *ind;
    mpf_t *val;
    mpf_ILLmatrix *A;
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
		lprows->rowval = mpf_EGlpNumAllocArray (lp->nzcount);
	    } else {
		ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount - lp->nrows, int);
		lprows->rowval = mpf_EGlpNumAllocArray (lp->nzcount - lp->nrows);
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
		    mpf_EGlpNumCopy (val[beg[A->matind[k]]], A->matval[k]);
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
	mpf_ILLlp_rows_clear (lprows);
    }
    ILL_IFFREE (hit, char);
    ILL_IFFREE (inv_structmap, int);

    ILL_RETURN (rval, "mpf_ILLlp_rows_init");
}

void mpf_ILLlp_rows_clear (mpf_ILLlp_rows * lprows)
{
    if (lprows != NULL) {
	ILL_IFFREE (lprows->rowbeg, int);
	ILL_IFFREE (lprows->rowcnt, int);
	ILL_IFFREE (lprows->rowind, int);
	mpf_EGlpNumFreeArray (lprows->rowval);
    }
}

static int mpf_wr_line (mpf_ILLlpdata * lp,
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

int mpf_ILLprint_report (mpf_ILLlpdata * lp,
      const char *format,
    ...)
{
    va_list marker;
    int rval = 0;

    va_start (marker, format);	/* ANSI style */
    rval = mpf_wr_line (lp, format, marker);
    va_end (marker);		/* Reset variable arguments.      */
    return rval;
}
