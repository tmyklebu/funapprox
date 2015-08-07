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

/* RCS_INFO = "$RCSfile: mpq_lpdata.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
/****************************************************************************/
/* */
/* Routines for Manipulating and Writing LPdata               */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int ILLlpdata_buildrows (mpq_ILLlpdata *lp, int **rowbeg, int **rowcnt,     */
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
/* void mpq_ILLlpdata_init (mpq_ILLlpdata *lp)                                   */
/* void mpq_ILLlpdata_free (mpq_ILLlpdata *lp)                                   */
/* */
/* void mpq_ILLlp_basis_init (mpq_ILLlp_basis *B)                                */
/* void mpq_ILLlp_basis_free (mpq_ILLlp_basis *B)                                */
/* int mpq_ILLlp_basis_alloc (mpq_ILLlp_basis *B, int nstruct, int nrows)        */
/* */
/* void mpq_ILLlp_cache_init (mpq_ILLlp_cache *C)                                */
/* void mpq_ILLlp_cache_free (mpq_ILLlp_cache *C)                                */
/* int mpq_ILLlp_cache_alloc (mpq_ILLlp_cache *C, int nstruct, int nrows)        */
/* */
/* void mpq_ILLlp_sinfo_init (mpq_ILLlp_sinfo *sinfo)                            */
/* void mpq_ILLlp_sinfo_free (mpq_ILLlp_sinfo *sinfo)                            */
/* */
/* int mpq_ILLlp_rows_init(mpq_ILLlp_rows *lprows, mpq_ILLlpdata *lp,                */
/* int include_logicals)          */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdata.h"
#include "mpq_qstruct.h"
#include "mpq_qsopt.h"
#include "mpq_lp.h"
#include "mpq_mps.h"
#include "mpq_rawlp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

mpq_t mpq_PARAM_IBASIS_RPIVOT;
mpq_t mpq_PARAM_IBASIS_RTRIANG;
mpq_t mpq_PARAM_MIN_DNORM;
mpq_t mpq_PFEAS_TOLER;
mpq_t mpq_BD_TOLER;
mpq_t mpq_DFEAS_TOLER;
mpq_t mpq_PIVOT_TOLER;
mpq_t mpq_SZERO_TOLER;
mpq_t mpq_PIVZ_TOLER;
mpq_t mpq_OBJBND_TOLER;
mpq_t mpq_DBNDPIV_TOLER;
mpq_t mpq_DBNDPIV_RATIO;
mpq_t mpq_PNSTEP_TOLER;
mpq_t mpq_DNSTEP_TOLER;
mpq_t mpq_ALTPIV_TOLER;
mpq_t mpq_DJZERO_TOLER;
mpq_t mpq_PROGRESS_ZERO;	/* 1e-7 */
mpq_t mpq_PROGRESS_THRESH;	/* 1e-5 */
mpq_t mpq_CB_EPS;
mpq_t mpq_CB_INF_RATIO;
mpq_t mpq_CB_PRI_RLIMIT;
mpq_t mpq_ILL_MAXDOUBLE;
mpq_t mpq_ILL_MINDOUBLE;

void mpq_ILLstart (void) __attribute__ ((constructor));
void mpq_ILLstart (void)
{
#if VERBOSE_LEVEL <= DEBUG
    char *strtmp = 0;
#endif
    if (!__EGlpNum_setup)
	EGlpNumStart ();
    mpq_EGlpNumInitVar (mpq_PARAM_IBASIS_RPIVOT);
    mpq_EGlpNumInitVar (mpq_PARAM_IBASIS_RTRIANG);
    mpq_EGlpNumInitVar (mpq_PARAM_MIN_DNORM);
    mpq_EGlpNumInitVar (mpq_PFEAS_TOLER);
    mpq_EGlpNumInitVar (mpq_BD_TOLER);
    mpq_EGlpNumInitVar (mpq_DFEAS_TOLER);
    mpq_EGlpNumInitVar (mpq_PIVOT_TOLER);
    mpq_EGlpNumInitVar (mpq_SZERO_TOLER);
    mpq_EGlpNumInitVar (mpq_PIVZ_TOLER);
    mpq_EGlpNumInitVar (mpq_OBJBND_TOLER);
    mpq_EGlpNumInitVar (mpq_DBNDPIV_TOLER);
    mpq_EGlpNumInitVar (mpq_DBNDPIV_RATIO);
    mpq_EGlpNumInitVar (mpq_PNSTEP_TOLER);
    mpq_EGlpNumInitVar (mpq_DNSTEP_TOLER);
    mpq_EGlpNumInitVar (mpq_ALTPIV_TOLER);
    mpq_EGlpNumInitVar (mpq_DJZERO_TOLER);
    mpq_EGlpNumInitVar (mpq_PROGRESS_ZERO);	/* 1e-7 */
    mpq_EGlpNumInitVar (mpq_PROGRESS_THRESH);	/* 1e-5 */
    mpq_EGlpNumInitVar (mpq_CB_PRI_RLIMIT);
    mpq_EGlpNumInitVar (mpq_CB_INF_RATIO);
    mpq_EGlpNumInitVar (mpq_CB_EPS);
    mpq_EGlpNumInitVar (mpq_ILL_MAXDOUBLE);
    mpq_EGlpNumInitVar (mpq_ILL_MINDOUBLE);
    /* parameters that do depend on the tolerance to zero */
    mpq_EGlpNumSet (mpq_PARAM_MIN_DNORM, 5e-9);
    mpq_EGlpNumMultTo (mpq_PARAM_MIN_DNORM, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PFEAS_TOLER, 5e9);
    mpq_EGlpNumMultTo (mpq_PFEAS_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_BD_TOLER, 5e10);
    mpq_EGlpNumMultTo (mpq_BD_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_DFEAS_TOLER, 5e9);
    mpq_EGlpNumMultTo (mpq_DFEAS_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PIVOT_TOLER, 5e5);
    mpq_EGlpNumMultTo (mpq_PIVOT_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_SZERO_TOLER, 5.0);
    mpq_EGlpNumMultTo (mpq_SZERO_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PIVZ_TOLER, 5e3);
    mpq_EGlpNumMultTo (mpq_PIVZ_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_OBJBND_TOLER, 5e13);
    mpq_EGlpNumMultTo (mpq_OBJBND_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_ALTPIV_TOLER, 5e7);
    mpq_EGlpNumMultTo (mpq_ALTPIV_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PROGRESS_ZERO, 5e8);
    mpq_EGlpNumMultTo (mpq_PROGRESS_ZERO, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PROGRESS_THRESH, 5e10);
    mpq_EGlpNumMultTo (mpq_PROGRESS_THRESH, mpq_epsLpNum);
#if VERBOSE_LEVEL <= DEBUG
    strtmp = mpq_EGlpNumGetStr (mpq_PARAM_MIN_DNORM);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_PARAM_MIN_DNORM to %lg %s", mpq_EGlpNumToLf (mpq_PARAM_MIN_DNORM), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_PFEAS_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_PFEAS_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_PFEAS_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_BD_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_BD_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_BD_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_DFEAS_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_DFEAS_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_DFEAS_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_PIVOT_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_PIVOT_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_PIVOT_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_SZERO_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_SZERO_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_SZERO_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_PIVZ_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_PIVZ_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_PIVZ_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_OBJBND_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_OBJBND_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_OBJBND_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_ALTPIV_TOLER);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_ALTPIV_TOLER to %lg %s", mpq_EGlpNumToLf (mpq_ALTPIV_TOLER), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_PROGRESS_ZERO);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_PROGRESS_ZERO to %lg %s", mpq_EGlpNumToLf (mpq_PROGRESS_ZERO), strtmp);
    EGfree (strtmp);
    strtmp = mpq_EGlpNumGetStr (mpq_PROGRESS_THRESH);
    MESSAGE (VERBOSE_LEVEL, "Setting mpq_PROGRESS_THRESH to %lg %s", mpq_EGlpNumToLf (mpq_PROGRESS_THRESH), strtmp);
    EGfree (strtmp);
#endif
    /* parameters that do not depend on the tolerance to zero */
    mpq_EGlpNumSet (mpq_ILL_MAXDOUBLE, 1e150);
    mpq_EGlpNumSet (mpq_ILL_MINDOUBLE, -1e150);
    mpq_EGlpNumSet (mpq_PARAM_IBASIS_RPIVOT, 0.98);
    mpq_EGlpNumSet (mpq_PARAM_IBASIS_RTRIANG, 0.01);
    mpq_EGlpNumSet (mpq_DBNDPIV_TOLER, 1e-3);
    mpq_EGlpNumSet (mpq_DBNDPIV_RATIO, 1e-2);
    mpq_EGlpNumSet (mpq_PNSTEP_TOLER, 1e-9);
    mpq_EGlpNumSet (mpq_DNSTEP_TOLER, 1e-9);
    mpq_EGlpNumSet (mpq_DJZERO_TOLER, 1e-8);
    mpq_EGlpNumSet (mpq_CB_EPS, 0.001);
    mpq_EGlpNumSet (mpq_CB_INF_RATIO, 10.0);
    mpq_EGlpNumSet (mpq_CB_PRI_RLIMIT, 0.25);
}

void mpq_ILLchange_precision (void)
{
    mpq_EGlpNumClearVar (mpq_PFEAS_TOLER);
    mpq_EGlpNumClearVar (mpq_BD_TOLER);
    mpq_EGlpNumClearVar (mpq_DFEAS_TOLER);
    mpq_EGlpNumClearVar (mpq_PIVOT_TOLER);
    mpq_EGlpNumClearVar (mpq_SZERO_TOLER);
    mpq_EGlpNumClearVar (mpq_PIVZ_TOLER);
    mpq_EGlpNumClearVar (mpq_OBJBND_TOLER);
    mpq_EGlpNumClearVar (mpq_ALTPIV_TOLER);
    mpq_EGlpNumClearVar (mpq_PARAM_MIN_DNORM);
    mpq_EGlpNumClearVar (mpq_PROGRESS_ZERO);
    mpq_EGlpNumClearVar (mpq_PROGRESS_THRESH);
    mpq_EGlpNumInitVar (mpq_PROGRESS_ZERO);
    mpq_EGlpNumInitVar (mpq_PROGRESS_THRESH);
    mpq_EGlpNumInitVar (mpq_PFEAS_TOLER);
    mpq_EGlpNumInitVar (mpq_BD_TOLER);
    mpq_EGlpNumInitVar (mpq_DFEAS_TOLER);
    mpq_EGlpNumInitVar (mpq_PIVOT_TOLER);
    mpq_EGlpNumInitVar (mpq_SZERO_TOLER);
    mpq_EGlpNumInitVar (mpq_PIVZ_TOLER);
    mpq_EGlpNumInitVar (mpq_OBJBND_TOLER);
    mpq_EGlpNumInitVar (mpq_ALTPIV_TOLER);
    mpq_EGlpNumInitVar (mpq_PARAM_MIN_DNORM);
    /* parameters that do depend on the tolerance to zero */
    mpq_EGlpNumSet (mpq_PARAM_MIN_DNORM, 5e-9);
    mpq_EGlpNumMultTo (mpq_PARAM_MIN_DNORM, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PFEAS_TOLER, 5e9);
    mpq_EGlpNumMultTo (mpq_PFEAS_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_BD_TOLER, 5e10);
    mpq_EGlpNumMultTo (mpq_BD_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_DFEAS_TOLER, 5e9);
    mpq_EGlpNumMultTo (mpq_DFEAS_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PIVOT_TOLER, 5e5);
    mpq_EGlpNumMultTo (mpq_PIVOT_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_SZERO_TOLER, 5.0);
    mpq_EGlpNumMultTo (mpq_SZERO_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PIVZ_TOLER, 5e3);
    mpq_EGlpNumMultTo (mpq_PIVZ_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_OBJBND_TOLER, 5e13);
    mpq_EGlpNumMultTo (mpq_OBJBND_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_ALTPIV_TOLER, 5e7);
    mpq_EGlpNumMultTo (mpq_ALTPIV_TOLER, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PROGRESS_ZERO, 5e8);
    mpq_EGlpNumMultTo (mpq_PROGRESS_ZERO, mpq_epsLpNum);
    mpq_EGlpNumSet (mpq_PROGRESS_THRESH, 5e10);
    mpq_EGlpNumMultTo (mpq_PROGRESS_THRESH, mpq_epsLpNum);
}

void mpq_ILLend (void) __attribute__ ((destructor));
void mpq_ILLend (void)
{
    mpq_EGlpNumClearVar (mpq_PARAM_IBASIS_RPIVOT);
    mpq_EGlpNumClearVar (mpq_PARAM_IBASIS_RTRIANG);
    mpq_EGlpNumClearVar (mpq_PARAM_MIN_DNORM);
    mpq_EGlpNumClearVar (mpq_PFEAS_TOLER);
    mpq_EGlpNumClearVar (mpq_BD_TOLER);
    mpq_EGlpNumClearVar (mpq_DFEAS_TOLER);
    mpq_EGlpNumClearVar (mpq_PIVOT_TOLER);
    mpq_EGlpNumClearVar (mpq_SZERO_TOLER);
    mpq_EGlpNumClearVar (mpq_PIVZ_TOLER);
    mpq_EGlpNumClearVar (mpq_OBJBND_TOLER);
    mpq_EGlpNumClearVar (mpq_DBNDPIV_TOLER);
    mpq_EGlpNumClearVar (mpq_DBNDPIV_RATIO);
    mpq_EGlpNumClearVar (mpq_PNSTEP_TOLER);
    mpq_EGlpNumClearVar (mpq_DNSTEP_TOLER);
    mpq_EGlpNumClearVar (mpq_ALTPIV_TOLER);
    mpq_EGlpNumClearVar (mpq_DJZERO_TOLER);
    mpq_EGlpNumClearVar (mpq_PROGRESS_ZERO);	/* 1e-7 */
    mpq_EGlpNumClearVar (mpq_PROGRESS_THRESH);	/* 1e-5 */
    mpq_EGlpNumClearVar (mpq_CB_EPS);
    mpq_EGlpNumClearVar (mpq_CB_INF_RATIO);
    mpq_EGlpNumClearVar (mpq_CB_PRI_RLIMIT);
    mpq_EGlpNumClearVar (mpq_ILL_MAXDOUBLE);
    mpq_EGlpNumClearVar (mpq_ILL_MINDOUBLE);
}

mpq_QSdata *mpq_ILLread (mpq_qsline_reader * file,
      const char *mpq_fname,
      int isMps)
{
    int rval = 0;
    mpq_QSdata *p = 0;
    mpq_ILLlpdata *lp;
    mpq_rawlpdata rawlp;

    ILL_FAILfalse (file != NULL, NULL);
    ILL_FAILfalse (mpq_fname != NULL, NULL);

    p = mpq_QScreate_prob (mpq_fname, QS_MIN);
    ILL_CHECKnull (p, NULL);
    ILL_IFFREE (p->qslp->probname, char);
    lp = p->qslp;

    mpq_ILLinit_rawlpdata (&rawlp, file->error_collector);
    mpq_ILLlpdata_init (lp);

    if (isMps != 0) {
	rval = mpq_ILLread_mps (file, mpq_fname, &rawlp);
    } else {
	rval = mpq_ILLread_lp (file, mpq_fname, &rawlp);
    }
    ILL_CLEANUP_IF (rval);

    rval = mpq_ILLrawlpdata_to_lpdata (&rawlp, lp);
    ILL_CLEANUP_IF (rval);

CLEANUP:
    mpq_ILLfree_rawlpdata (&rawlp);
    if (rval != 0) {
	mpq_QSfree_prob (p);
	p = 0;
    }
    return p;
}

void mpq_ILLlpdata_init (mpq_ILLlpdata * lp)
{
    if (lp) {
	lp->nrows = 0;
	lp->ncols = 0;
	lp->nstruct = 0;
	lp->nzcount = 0;
	lp->rowsize = 0;
	lp->colsize = 0;
	lp->structsize = 0;
	lp->objsense = mpq_ILL_MIN;
	lp->sense = 0;
	lp->obj = 0;
	lp->rhs = 0;
	lp->rangeval = 0;
	lp->lower = 0;
	lp->upper = 0;

	mpq_ILLmatrix_init (&lp->A);
	mpq_ILLmatrix_init (&lp->sos);
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

void mpq_ILLlpdata_free (mpq_ILLlpdata * lp)
{
    int i;

    if (lp) {
	ILL_IFFREE (lp->sense, char);
	mpq_EGlpNumFreeArray (lp->obj);
	mpq_EGlpNumFreeArray (lp->rhs);
	mpq_EGlpNumFreeArray (lp->rangeval);
	mpq_EGlpNumFreeArray (lp->lower);
	mpq_EGlpNumFreeArray (lp->upper);
	mpq_ILLmatrix_free (&lp->A);
	if (lp->rA) {
	    mpq_ILLlp_rows_clear (lp->rA);
	    ILL_IFFREE (lp->rA, mpq_ILLlp_rows);
	}
	ILL_IFFREE (lp->is_sos_mem, int);
	ILL_IFFREE (lp->refrowname, char);
	mpq_ILLmatrix_free (&lp->sos);
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
	    mpq_ILLlp_sinfo_free (lp->sinfo);
	    ILL_IFFREE (lp->sinfo, mpq_ILLlp_sinfo);
	}
	mpq_ILLlpdata_init (lp);
    }
}

void mpq_ILLlp_basis_init (mpq_ILLlp_basis * B)
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

void mpq_ILLlp_basis_free (mpq_ILLlp_basis * B)
{
    if (B) {
	ILL_IFFREE (B->cstat, char);
	ILL_IFFREE (B->rstat, char);
	mpq_EGlpNumFreeArray (B->rownorms);
	mpq_EGlpNumFreeArray (B->colnorms);
	B->nstruct = 0;
	B->nrows = 0;
    }
}

int mpq_ILLlp_basis_alloc (mpq_ILLlp_basis * B,
      int nstruct,
      int nrows)
{
    int rval = 0;

    ILL_FAILtrue (B == NULL, "mpq_ILLlp_basis_alloc called without a basis");

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
	mpq_ILLlp_basis_free (B);
    }
    ILL_RETURN (rval, "mpq_ILLlp_basis_alloc");
}

void mpq_ILLlp_cache_init (mpq_ILLlp_cache * C)
{
    if (C) {
	C->x = 0;
	C->rc = 0;
	C->pi = 0;
	C->slack = 0;
	C->nstruct = 0;
	C->nrows = 0;
	C->status = 0;
	mpq_EGlpNumZero (C->val);
    }
}

void mpq_ILLlp_cache_free (mpq_ILLlp_cache * C)
{
    if (C) {
	mpq_EGlpNumFreeArray (C->x);
	mpq_EGlpNumFreeArray (C->rc);
	mpq_EGlpNumFreeArray (C->pi);
	mpq_EGlpNumFreeArray (C->slack);
	C->nstruct = 0;
	C->nrows = 0;
	C->status = 0;
    }
}

int mpq_ILLlp_cache_alloc (mpq_ILLlp_cache * C,
      int nstruct,
      int nrows)
{
    int rval = 0;

    ILL_FAILtrue (C == NULL, "mpq_ILLlp_cache_alloc called without a cache");

    C->nstruct = nstruct;
    C->nrows = nrows;

    if (nstruct > 0) {
	C->x = mpq_EGlpNumAllocArray (nstruct);
	C->rc = mpq_EGlpNumAllocArray (nstruct);
    }
    if (nrows > 0) {
	C->pi = mpq_EGlpNumAllocArray (nrows);
	C->slack = mpq_EGlpNumAllocArray (nrows);
    }
CLEANUP:

    if (rval) {
	mpq_ILLlp_cache_free (C);
    }
    ILL_RETURN (rval, "mpq_ILLlp_cache_alloc");
}


int mpq_ILLlp_rows_init (mpq_ILLlp_rows * lprows,
      mpq_ILLlpdata * lp,
      int include_logicals)
{
    int rval = 0;
    int i, k, st;
    int *beg, *cnt, *ind;
    mpq_t *val;
    mpq_ILLmatrix *A;
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
		lprows->rowval = mpq_EGlpNumAllocArray (lp->nzcount);
	    } else {
		ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount - lp->nrows, int);
		lprows->rowval = mpq_EGlpNumAllocArray (lp->nzcount - lp->nrows);
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
		    mpq_EGlpNumCopy (val[beg[A->matind[k]]], A->matval[k]);
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
	mpq_ILLlp_rows_clear (lprows);
    }
    ILL_IFFREE (hit, char);
    ILL_IFFREE (inv_structmap, int);

    ILL_RETURN (rval, "mpq_ILLlp_rows_init");
}

void mpq_ILLlp_rows_clear (mpq_ILLlp_rows * lprows)
{
    if (lprows != NULL) {
	ILL_IFFREE (lprows->rowbeg, int);
	ILL_IFFREE (lprows->rowcnt, int);
	ILL_IFFREE (lprows->rowind, int);
	mpq_EGlpNumFreeArray (lprows->rowval);
    }
}

static int mpq_wr_line (mpq_ILLlpdata * lp,
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

int mpq_ILLprint_report (mpq_ILLlpdata * lp,
      const char *format,
    ...)
{
    va_list marker;
    int rval = 0;

    va_start (marker, format);	/* ANSI style */
    rval = mpq_wr_line (lp, format, marker);
    va_end (marker);		/* Reset variable arguments.      */
    return rval;
}
