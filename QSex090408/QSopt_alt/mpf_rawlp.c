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

/* RCS_INFO = "$RCSfile: mpf_rawlp.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */
/****************************************************************************/
/* DataStructure and routines to deal with raw lp information as read       */
/* from mps or lp files.                                                    */
/****************************************************************************/

#include "econfig.h"
#include "config.h"
#include "mpf_sortrus.h"
#include "mpf_iqsutil.h"
#include "mpf_rawlp.h"
#include "allocrus.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
static int TRACE = 0;

ILL_PTRWORLD_ROUTINES (mpf_colptr, colptralloc, colptr_bulkalloc, colptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (mpf_colptr, colptr_check_leaks, this, int)
const int mpf_ILL_SOS_TYPE1 = 1;
const int mpf_ILL_SOS_TYPE2 = 2;

static void mpf_ILLprt_EGlpNum (FILE * f,
      mpf_t * d)
{
    if (mpf_EGlpNumIsLeq (mpf_ILL_MAXDOUBLE, *d)) {
	fprintf (f, "MAX_DOUBLE");
    } else {
	if (mpf_EGlpNumIsLeq (*d, mpf_ILL_MINDOUBLE)) {
	    fprintf (f, "-MAX_DOUBLE");
	} else {
	    fprintf (f, "%f", mpf_EGlpNumToLf (*d));
	}
    }
}

static int mpf_ILLraw_check_bounds (mpf_rawlpdata * lp);

void mpf_ILLinit_rawlpdata (mpf_rawlpdata * lp,
      mpf_qserror_collector * collector)
{
    if (lp) {
	lp->name = 0;
	lp->ncols = 0;
	lp->nrows = 0;
	lp->cols = 0;
	lp->rowsense = 0;
	lp->rhsname = 0;
	lp->rhs = 0;
	lp->rhsind = 0;
	lp->rangesname = 0;
	lp->rangesind = 0;
	lp->ranges = 0;
	lp->boundsname = 0;
	lp->lbind = 0;
	lp->ubind = 0;
	lp->lower = 0;
	lp->upper = 0;
	lp->intmarker = 0;
	lp->colsize = 0;
	lp->sensesize = 0;
	lp->intsize = 0;
	lp->rhssize = 0;
	lp->refrow = NULL;
	lp->is_sos_size = 0;
	lp->is_sos_member = NULL;
	lp->nsos_member = 0;
	lp->sos_weight_size = 0;
	lp->sos_weight = NULL;
	lp->sos_col_size = 0;
	lp->sos_col = NULL;
	lp->nsos = 0;
	lp->sos_setsize = 0;
	lp->sos_set = NULL;
	ILLsymboltab_init (&lp->coltab);
	ILLsymboltab_init (&lp->rowtab);
	lp->objindex = -1;
	lp->objsense = mpf_ILL_MIN;
	lp->refrowind = -1;	/* undefined */
	ILLptrworld_init (&lp->ptrworld);
	lp->error_collector = collector;
    }
}

void mpf_ILLraw_clear_matrix (mpf_rawlpdata * lp)
{
    int i;
    mpf_colptr *next, *curr;
    if ((lp != NULL) && (lp->cols != NULL)) {
	for (i = 0; i < lp->ncols; i++) {
	    {
		curr = lp->cols[i];
		while (curr) {
		    next = curr->next;
		    mpf_EGlpNumClearVar ((curr->coef));
		    colptrfree (&(lp->ptrworld), curr);
		    curr = next;
		}
	    }
	    lp->cols[i] = NULL;
	}
    }
}

void mpf_ILLfree_rawlpdata (mpf_rawlpdata * lp)
{
    int total, onlist;
    mpf_colptr *next, *curr;

    if (lp) {
	ILL_IFFREE (lp->name, char);
	ILLsymboltab_free (&lp->rowtab);
	ILLsymboltab_free (&lp->coltab);
	ILL_IFFREE (lp->rowsense, char);
	mpf_ILLraw_clear_matrix (lp);
	ILL_IFFREE (lp->cols, mpf_colptr *);
	{
	    curr = lp->ranges;
	    while (curr) {
		next = curr->next;
		mpf_EGlpNumClearVar ((curr->coef));
		colptrfree (&(lp->ptrworld), curr);
		curr = next;
	    }
	}
	if (colptr_check_leaks (&lp->ptrworld, &total, &onlist)) {
	    fprintf (stderr, "WARNING: %d outstanding colptrs\n", total - onlist);
	}
	ILLptrworld_delete (&lp->ptrworld);
	ILL_IFFREE (lp->rhsname, char);
	mpf_EGlpNumFreeArray (lp->rhs);
	ILL_IFFREE (lp->rhsind, char);
	ILL_IFFREE (lp->rangesname, char);
	ILL_IFFREE (lp->rangesind, char);
	ILL_IFFREE (lp->boundsname, char);
	ILL_IFFREE (lp->lbind, char);
	ILL_IFFREE (lp->ubind, char);
	mpf_EGlpNumFreeArray (lp->lower);
	mpf_EGlpNumFreeArray (lp->upper);
	ILL_IFFREE (lp->intmarker, char);
	ILL_IFFREE (lp->refrow, char);
	ILL_IFFREE (lp->is_sos_member, int);
	mpf_EGlpNumFreeArray (lp->sos_weight);
	ILL_IFFREE (lp->sos_col, int);
	ILL_IFFREE (lp->sos_set, mpf_sosptr);
	mpf_ILLinit_rawlpdata (lp, NULL);
    }
}

const char *mpf_ILLraw_rowname (mpf_rawlpdata * lp,
      int i)
{
    const char *name = NULL;
    ILL_FAILfalse_no_rval ((i >= 0) && (i < lp->nrows), "index out of range");
    ILL_FAILfalse_no_rval (lp->nrows == lp->rowtab.tablesize,
	"tab and lp must be in synch");
    name = ILLsymboltab_get (&lp->rowtab, i);
CLEANUP:
    return name;
}
const char *mpf_ILLraw_colname (mpf_rawlpdata * lp,
      int i)
{
    const char *name = NULL;
    ILL_FAILfalse_no_rval ((i >= 0) && (i < lp->ncols), "index out of range");
    ILL_FAILfalse_no_rval (lp->ncols == lp->coltab.tablesize,
	"tab and lp must be in synch");
    name = ILLsymboltab_get (&lp->coltab, i);
CLEANUP:
    return name;
}

int mpf_ILLraw_add_col (mpf_rawlpdata * lp,
      const char *name,
      int intmarker)
{
    int rval = 0;
    int pindex, hit;

    rval = ILLsymboltab_register (&lp->coltab, name, -1, &pindex, &hit);
    rval = rval || hit;
    ILL_CLEANUP_IF (rval);
    if (lp->ncols >= lp->colsize) {
	lp->colsize *= 1.3;
	lp->colsize += 1000;
	if (lp->colsize < lp->ncols + 1)
	    lp->colsize = lp->ncols + 1;
	lp->cols = EGrealloc (lp->cols, lp->colsize * sizeof (mpf_colptr *));
	/*
	        rval = rval || ILLutil_reallocrus_scale (&lp->cols, &lp->colsize, lp->ncols + 1, 1.3, sizeof (mpf_colptr *));
	*/
    }
    if (lp->ncols >= lp->intsize) {
	lp->intsize *= 1.3;
	lp->intsize += 1000;
	if (lp->intsize < lp->ncols + 1)
	    lp->intsize = lp->ncols + 1;
	lp->intmarker = EGrealloc (lp->intmarker, lp->intsize * sizeof (char));
	/*
	        rval = rval || ILLutil_reallocrus_scale ((void **) &lp->intmarker, &lp->intsize, lp->ncols + 1, 1.3, sizeof (char));
	*/
    }
    if (lp->ncols >= lp->is_sos_size) {
	lp->is_sos_size *= 1.3;
	lp->is_sos_size += 1000;
	if (lp->is_sos_size < lp->ncols + 1)
	    lp->is_sos_size = lp->ncols + 1;
	lp->is_sos_member = EGrealloc (lp->is_sos_member,
	    sizeof (int) * lp->is_sos_size);
	/*
	        rval = rval || ILLutil_reallocrus_scale ((void **) &lp->is_sos_member, &lp->is_sos_size, lp->ncols + 1, 1.3, sizeof (int));
	*/
    }
    ILL_CLEANUP_IF (rval);
    lp->cols[lp->ncols] = 0;
    lp->is_sos_member[lp->ncols] = -1;
    lp->intmarker[lp->ncols] = intmarker;
    lp->ncols++;
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_add_col");
}

int mpf_ILLraw_init_rhs (mpf_rawlpdata * lp)
{
    int i, rval = 0;

    ILL_FAILfalse (lp->rhsind == NULL, "Should be called exactly once");
    if (lp->nrows > 0) {
	ILL_SAFE_MALLOC (lp->rhsind, lp->nrows, char);
	for (i = 0; i < lp->nrows; i++) {
	    lp->rhsind[i] = (char) 0;
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_init_rhs");
}

int mpf_ILLraw_init_ranges (mpf_rawlpdata * lp)
{
    int i, rval = 0;

    ILL_FAILfalse (lp->rangesind == NULL, "Should be called exactly once");
    if (lp->nrows > 0) {
	ILL_SAFE_MALLOC (lp->rangesind, lp->nrows, char);
	for (i = 0; i < lp->nrows; i++) {
	    lp->rangesind[i] = (char) 0;
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_init_ranges");
}

int mpf_ILLraw_add_col_coef (mpf_rawlpdata * lp,
      int colind,
      int rowind,
      mpf_t coef)
{
    mpf_colptr *cp = mpf_ILLcolptralloc (&lp->ptrworld);
    if (!cp) {
	return 1;
    }
    cp->this = rowind;
    mpf_EGlpNumCopy (cp->coef, coef);
    cp->next = lp->cols[colind];
    lp->cols[colind] = cp;
    return 0;
}


int mpf_ILLraw_add_ranges_coef (mpf_rawlpdata * lp,
      int rowind,
      mpf_t coef)
{
    mpf_colptr *cp = mpf_ILLcolptralloc (&lp->ptrworld);
    if (!cp) {
	return 1;
    }
    cp->this = rowind;
    mpf_EGlpNumCopy (cp->coef, coef);
    cp->next = lp->ranges;
    lp->ranges = cp;
    lp->rangesind[rowind] = (char) 1;
    return 0;
}

int mpf_ILLraw_add_sos (mpf_rawlpdata * lp,
      int tp)
{
    int rval = 0;
    mpf_sosptr *sos, *bef;

    if (lp->nsos >= lp->sos_setsize) {
	lp->sos_setsize *= 1.3;
	lp->sos_setsize += 1000;
	if (lp->sos_setsize < lp->nsos + 1)
	    lp->sos_setsize = lp->nsos + 1;
	lp->sos_set = EGrealloc (lp->sos_set, sizeof (mpf_sosptr *) * lp->sos_setsize);
	/*
	        if (ILLutil_reallocrus_scale ((void **) &lp->sos_set,
	                                                                    &lp->sos_setsize, lp->nsos + 1, 1.3,
	                                                                    sizeof (mpf_sosptr *)))
	        {
	            ILL_CLEANUP_IF (rval);
	        }
	*/
    }
    sos = lp->sos_set + lp->nsos;
    sos->nelem = 0;
    sos->type = tp;
    if (lp->nsos == 0) {
	sos->first = 0;
    } else {
	bef = &(lp->sos_set[lp->nsos - 1]);
	sos->first = bef->first + bef->nelem;
    }
    lp->nsos++;
    /* CLEANUP: */
    ILL_RETURN (rval, "mpf_ILLraw_add_sos");
}

int mpf_ILLraw_is_mem_other_sos (mpf_rawlpdata * lp,
      int colind)
{
    return (lp->is_sos_member[colind] >= 0) &&
    (lp->is_sos_member[colind] != (lp->nsos - 1));
}

int mpf_ILLraw_add_sos_member (mpf_rawlpdata * lp,
      int colind)
{
    int rval = 0;
    ILL_FAILfalse (lp->nsos > 0, "we should have called mpf_ILLraw_add_sos earlier");
    ILL_FAILtrue (mpf_ILLraw_is_mem_other_sos (lp, colind),
	"colind is member of another sos set");

    if (lp->is_sos_member[colind] == -1) {
	if (lp->nsos_member >= lp->sos_weight_size) {
	    lp->sos_weight_size *= 1.3;
	    lp->sos_weight_size += 1000;
	    if (lp->sos_weight_size < lp->nsos_member + 1)
		lp->sos_weight_size = lp->nsos_member + 1;
	    lp->sos_weight = EGrealloc (lp->sos_weight,
		lp->sos_weight_size * sizeof (double));
	    /*
	                if (ILLutil_reallocrus_scale ((void **) &lp->sos_weight,
	                                                                            &lp->sos_weight_size,
	                                                                            lp->nsos_member + 1, 1.3, sizeof (double)))
	                {
	                    ILL_CLEANUP_IF (rval);
	                }
	    */
	}
	if (lp->nsos_member >= lp->sos_col_size) {
	    lp->sos_col_size *= 1.3;
	    lp->sos_col_size += 1000;
	    if (lp->sos_col_size < lp->nsos_member + 1)
		lp->sos_col_size = lp->nsos_member + 1;
	    lp->sos_col = EGrealloc (lp->sos_col, sizeof (int) * lp->sos_col_size);
	    /*
	                if (ILLutil_reallocrus_scale ((void **) &lp->sos_col,
	                                                                            &lp->sos_col_size,
	                                                                            lp->nsos_member + 1, 1.3, sizeof (int)))
	                {
	                    ILL_CLEANUP_IF (rval);
	                }
	    */
	}
	lp->sos_col[lp->nsos_member] = colind;
	lp->sos_set[lp->nsos - 1].nelem++;
	lp->is_sos_member[colind] = lp->nsos - 1;
	lp->nsos_member++;
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_add_sos_member");
}


int mpf_ILLraw_add_row (mpf_rawlpdata * lp,
      const char *name,
      int sense,
      mpf_t rhs)
{
    int pindex, hit, rval = 0;

    rval = ILLsymboltab_register (&lp->rowtab, name, -1, &pindex, &hit);
    rval = rval || hit;
    ILL_CLEANUP_IF (rval);
    if (lp->nrows >= lp->sensesize) {
	lp->sensesize *= 1.3;
	lp->sensesize += 1000;
	if (lp->sensesize < lp->nrows + 1)
	    lp->sensesize = lp->nrows + 1;
	lp->rowsense = EGrealloc (lp->rowsense, sizeof (char) * lp->sensesize);
	/*
	        if (ILLutil_reallocrus_scale ((void **) &lp->rowsense,
	                                                                    &lp->sensesize, lp->nrows + 1,
	                                                                    1.3, sizeof (char)))
	        {
	            ILL_CLEANUP_IF (rval);
	        }
	*/
    }
    if (lp->nrows >= lp->rhssize) {
	if (lp->rhssize + 1000 < (lp->nrows + 1) * 1.3)
	    lp->rhssize = (lp->nrows + 1) * 1.3;
	else
	    lp->rhssize += 1000;
	mpf_EGlpNumReallocArray (&(lp->rhs), lp->rhssize);
    }
    lp->rowsense[lp->nrows] = sense;
    mpf_EGlpNumCopy (lp->rhs[lp->nrows], rhs);
    lp->nrows++;

CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_add_row");
}

static int mpf_ILLcheck_rawlpdata (mpf_rawlpdata * lp)
{
    int i, col, rval = 0;
    int si, *perm = NULL;
    const char *c1, *c2;
    mpf_sosptr *set;

    ILL_FAILfalse (lp, "lp must not be NULL");

    /* check    *) that there is at least one variable *) that the weights in
       all SOS sets are distinct *) all sos members are non integer variables
       ) sos set members have distint weights *) objindex is not -1 *)
       INVARIANT: rowname[objindex] != NULL *) INVARIANT: upper/lower arrays
       are filled in *) INVARIANT: if col or rownames != NULL then all their
       elements are not NULL */
    if (lp->ncols < 1) {
	return mpf_ILLdata_error (lp->error_collector, "There are no variables.");
    }
    if (lp->objindex == -1) {
	return mpf_ILLdata_error (lp->error_collector, "There is no objective fct.");
    }
    ILL_FAILfalse (mpf_ILLraw_rowname (lp, lp->objindex) != NULL,
	"must have objective name");
    if (lp->nsos_member > 1) {
	ILL_SAFE_MALLOC (perm, lp->nsos_member, int);
	for (si = 0; si < lp->nsos; si++) {
	    set = lp->sos_set + si;
	    for (i = 0; i < set->nelem; i++) {
		col = lp->sos_col[set->first + i];
		if (lp->intmarker[col]) {
		    rval = mpf_ILLdata_error (lp->error_collector,
			"SOS set member \"%s\" is an %s.\n",
			mpf_ILLraw_colname (lp, col),
			"integer/binary variable");
		}
	    }
	    if (set->nelem > 1) {
		for (i = 0; i < set->nelem; i++) {
		    perm[i] = set->first + i;
		}
		mpf_ILLutil_EGlpNum_perm_quicksort (perm, lp->sos_weight, set->nelem);
		for (i = 1; i < set->nelem; i++) {
		    if (mpf_EGlpNumIsEqqual
			(lp->sos_weight[perm[i - 1]], lp->sos_weight[perm[i]])) {
			c1 = mpf_ILLraw_colname (lp, lp->sos_col[perm[i]]);
			c2 = mpf_ILLraw_colname (lp, lp->sos_col[perm[i - 1]]);
			mpf_ILLdata_error (lp->error_collector,
			    "\"%s\" and \"%s\" both have %s %f.\n", c1, c2,
			    "SOS weight", lp->sos_weight[perm[i]]);
			rval = 1;
		    }
		}
	    }
	}
    }
    for (i = 0; i < lp->ncols; i++) {
	ILL_CHECKnull (mpf_ILLraw_colname (lp, i), "There is a NULL col name");
    }
    for (i = 0; i < lp->nrows; i++) {
	ILL_CHECKnull (mpf_ILLraw_rowname (lp, i), "There is a NULL row name");
    }
    ILL_FAILtrue ((lp->upper == NULL) | (lp->lower == NULL),
	"Upper/Lower arrays must be filled in.");

    rval += mpf_ILLraw_check_bounds (lp);
CLEANUP:
    ILL_IFFREE (perm, int);
    ILL_RESULT (rval, "mpf_ILLcheck_rawlpdata");
}

int mpf_ILLraw_init_bounds (mpf_rawlpdata * lp)
{
    int i, rval = 0;

    ILL_FAILfalse (lp->upper == NULL, "Should be called exactly once");
    ILL_FAILfalse (lp->lower == NULL, "Should be called exactly once");
    ILL_FAILfalse (lp->lbind == NULL, "Should be called exactly once");
    ILL_FAILfalse (lp->ubind == NULL, "Should be called exactly once");
    lp->upper = mpf_EGlpNumAllocArray (lp->ncols);
    lp->lower = mpf_EGlpNumAllocArray (lp->ncols);
    ILL_SAFE_MALLOC (lp->lbind, lp->ncols, char);
    ILL_SAFE_MALLOC (lp->ubind, lp->ncols, char);

    for (i = 0; i < lp->ncols; i++) {
	lp->lbind[i] = (char) 0;
	lp->ubind[i] = (char) 0;
	mpf_EGlpNumZero (lp->lower[i]);
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_init_bounds");
}

const char *mpf_ILLraw_set_lowerBound (mpf_rawlpdata * lp,
      int i,
      mpf_t bnd)
{
    ILL_FAILtrue_no_rval (i >= lp->ncols, "proper colind");
    if (lp->lbind[i]) {
	return "Using previous bound definition.";
    }
    mpf_EGlpNumCopy (lp->lower[i], bnd);
    lp->lbind[i] = (char) 1;
CLEANUP:
    return NULL;
}

const char *mpf_ILLraw_set_upperBound (mpf_rawlpdata * lp,
      int i,
      mpf_t bnd)
{
    ILL_FAILtrue_no_rval (i >= lp->ncols, "proper colind");
    if (lp->ubind[i]) {
	return "Using previous bound definition.";
    }
    mpf_EGlpNumCopy (lp->upper[i], bnd);
    lp->ubind[i] = (char) 1;
    if (mpf_EGlpNumIsEqqual (lp->lower[i], mpf_zeroLpNum) &&
	mpf_EGlpNumIsEqqual (bnd, mpf_zeroLpNum)) {
	return "0.0 upper bound fixes variable.";
    }
CLEANUP:
    return NULL;
}

const char *mpf_ILLraw_set_fixedBound (mpf_rawlpdata * lp,
      int i,
      mpf_t bnd)
{
    ILL_FAILtrue_no_rval (i >= lp->ncols, "proper colind");
    if (lp->ubind[i] || lp->lbind[i]) {
	return "Using previous bound definition.";
    }
    mpf_EGlpNumCopy (lp->lower[i], bnd);
    lp->lbind[i] = (char) 1;
    mpf_EGlpNumCopy (lp->upper[i], bnd);
    lp->ubind[i] = (char) 1;
CLEANUP:
    return NULL;
}

const char *mpf_ILLraw_set_unbound (mpf_rawlpdata * lp,
      int i)
{
    ILL_FAILtrue_no_rval (i >= lp->ncols, "proper colind");
    if (lp->lbind[i] || lp->ubind[i]) {
	return "Using previous bound definition.";
    }
    mpf_EGlpNumCopy (lp->lower[i], mpf_ILL_MINDOUBLE);
    mpf_EGlpNumCopy (lp->upper[i], mpf_ILL_MAXDOUBLE);
    lp->lbind[i] = 1;
    lp->ubind[i] = 1;
CLEANUP:
    return NULL;
}

const char *mpf_ILLraw_set_binaryBound (mpf_rawlpdata * lp,
      int i)
{
    ILL_FAILtrue_no_rval (i >= lp->ncols, "proper colind");
    if (lp->lbind[i] || lp->ubind[i]) {
	return "Using previous bound definition.";
    }
    mpf_EGlpNumZero (lp->lower[i]);
    mpf_EGlpNumOne (lp->upper[i]);
    lp->lbind[i] = 1;
    lp->ubind[i] = 1;
CLEANUP:
    return NULL;
}

int mpf_ILLraw_fill_in_bounds (mpf_rawlpdata * lp)
{
    int rval = 0, i;
    if (lp->lbind == NULL) {
	mpf_ILLraw_init_bounds (lp);
    }
    ILL_FAILtrue (lp->upper == NULL, "must all be there now");
    ILL_FAILtrue (lp->lower == NULL, "must all be there now");
    ILL_FAILtrue (lp->lbind == NULL, "must all be there now");
    ILL_FAILtrue (lp->ubind == NULL, "must all be there now");
    for (i = 0; i < lp->ncols; i++) {
	if (!lp->lbind[i]) {
	    if (lp->ubind[i] && mpf_EGlpNumIsLess (lp->upper[i], mpf_zeroLpNum)) {
		mpf_EGlpNumCopy (lp->lower[i], mpf_ILL_MINDOUBLE);
	    }
	}
	if (!lp->ubind[i]) {
	    /* int vars without bounds are binary                        */
	    /* all, also int vars                                        */
	    /* with explicit lower bound 0.0 are in [0.0,+inf]  */
	    if (((lp->intmarker != NULL) && lp->intmarker[i]) && !lp->lbind[i]) {
		mpf_EGlpNumOne (lp->upper[i]);
	    } else {
		mpf_EGlpNumCopy (lp->upper[i], mpf_ILL_MAXDOUBLE);
	    }
	}
    }

CLEANUP:
    if (rval) {
	mpf_EGlpNumFreeArray (lp->lower);
	mpf_EGlpNumFreeArray (lp->upper);
    }
    ILL_RETURN (rval, "mpf_ILLraw_fill_in_bounds");
}

static int mpf_ILLraw_check_bounds (mpf_rawlpdata * lp)
{
    int rval = 0, i;
    ILL_FAILtrue (lp->upper == NULL, "must all be there now");
    ILL_FAILtrue (lp->lower == NULL, "must all be there now");
    ILL_FAILtrue (lp->lbind == NULL, "must all be there now");
    ILL_FAILtrue (lp->ubind == NULL, "must all be there now");
    for (i = 0; i < lp->ncols; i++) {
	if (mpf_EGlpNumIsLess (lp->upper[i], lp->lower[i])) {
	    rval += mpf_ILLdata_error (lp->error_collector,
		"Lower bound is bigger than %s \"%s\".\n",
		"upper bound for", mpf_ILLraw_colname (lp, i));
	}
    }
    ILL_RESULT (rval, "mpf_ILLraw_check_bounds");
CLEANUP:
    ILL_RETURN (rval, "mpf_ILLraw_check_bounds");
}

int mpf_ILLraw_first_nondefault_bound (mpf_ILLlpdata * lp)
{
    int ri = lp->nstruct, i;
    ILL_FAILtrue_no_rval (lp->lower == NULL || lp->upper == NULL,
	"Should not call mpf_write_bounds when lower or upper are NULL");
    for (ri = 0; ri < lp->nstruct; ri++) {
	i = lp->structmap[ri];
	if (!mpf_ILLraw_default_lower (lp, i) || !mpf_ILLraw_default_upper (lp, i))
	    break;
    }
CLEANUP:
    return ri;
}

int mpf_ILLraw_default_lower (mpf_ILLlpdata * lp,
      int i)
{
    ILL_FAILtrue_no_rval (lp->lower == NULL || lp->upper == NULL,
	"Should not call mpf_write_bounds when lower or upper are NULL");
    ILL_FAILfalse_no_rval (lp->ncols > i, "i is not col index");
    if (mpf_EGlpNumIsEqqual (lp->lower[i], mpf_zeroLpNum) &&
	mpf_EGlpNumIsLeq (mpf_zeroLpNum, lp->upper[i])) {
	return 1;
    }
    if (mpf_EGlpNumIsEqqual (lp->lower[i], mpf_ILL_MINDOUBLE) &&
	mpf_EGlpNumIsLess (lp->upper[i], mpf_zeroLpNum)) {
	return 1;
    }
CLEANUP:
    return 0;
}

int mpf_ILLraw_default_upper (mpf_ILLlpdata * lp,
      int i)
{
    int isInt;

    ILL_FAILtrue_no_rval (lp->lower == NULL || lp->upper == NULL,
	"Should not call mpf_write_bounds when lower or upper are NULL");
    ILL_FAILfalse_no_rval (lp->ncols >= i, "i is not col index");
    isInt = (lp->intmarker != NULL) && lp->intmarker[i];
    if (isInt) {
	if (mpf_EGlpNumIsEqqual (lp->lower[i], mpf_zeroLpNum)) {
	    return (mpf_EGlpNumIsEqqual (lp->upper[i], mpf_oneLpNum));
	}
    }
    if (mpf_EGlpNumIsEqqual (lp->upper[i], mpf_ILL_MAXDOUBLE)) {
	return 1;
    }
CLEANUP:
    return 0;
}

int mpf_ILLraw_fill_in_rownames (mpf_rawlpdata * lp)
{
    int i, rval = 0;
    char uname[ILL_namebufsize];
    ILLsymboltab *rowtab;
    char first = 1;

    rowtab = &lp->rowtab;
    ILL_FAILtrue (lp->nrows != rowtab->tablesize, "must have same #entries");
    for (i = 0; (rval == 0) && i < lp->nrows; i++) {
	if (ILLsymboltab_get (rowtab, i) == NULL) {
	    if (first) {
		mpf_ILLdata_warn (lp->error_collector,
		    "Generating names for unnamed rows.");
		first = 0;
	    }
	    ILLsymboltab_unique_name (rowtab, i, "c", uname);
	    rval = ILLsymboltab_rename (rowtab, i, uname);
	    ILL_CLEANUP_IF (rval);
	}
    }
CLEANUP:
    ILL_RESULT (rval, "mpf_ILLraw_fill_in_rownames");
}

static int mpf_whichColsAreUsed (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *colindex)
{
    int rval = 0;
    int i, objind = raw->objindex;
    mpf_colptr *cp;
    char *colUsed = NULL;

    /* colUsed[i]  variable raw->colnames[i] is used in obj fct and/or
       equation(s) */
    ILL_SAFE_MALLOC (colUsed, raw->ncols, char);
    for (i = 0; i < raw->ncols; i++) {
	colUsed[i] = 0;
    }
    for (i = 0; i < raw->ncols; i++) {
	for (cp = raw->cols[i]; cp; cp = cp->next) {
	    if ((cp->this == objind) || (raw->rowsense[cp->this] != 'N')) {
		colUsed[i] = 1;
		break;
	    }
	}
    }

    /* colindex[i] = -1 for undefined, 0, 1, ... lp->ncol-1 lp->ncols <=
       raw->ncols */
    for (i = 0; i < raw->ncols; i++) {
	if (colUsed[i]) {
	    colindex[i] = lp->ncols++;
	} else {
	    colindex[i] = -1;
	    mpf_ILLdata_warn (raw->error_collector,
		"\"%s\" is used in non objective 'N' rows only.",
		mpf_ILLraw_colname (raw, i));
	}
    }
    if (lp->ncols < 1) {
	rval = mpf_ILLdata_error (raw->error_collector, "There are no variables.");
	ILL_CLEANUP_IF (rval);
    }
CLEANUP:
    ILL_IFFREE (colUsed, char);
    ILL_RESULT (rval, "mpf_whichColsAreUsed");
}

static int mpf_whichRowsAreUsed (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *rowindex)
{
    int i, rval = 0;

    /* only use non 'N' rows */
    for (i = 0; i < raw->nrows; i++) {
	if (raw->rowsense[i] != 'N') {
	    rowindex[i] = lp->nrows++;
	} else {
	    rowindex[i] = -1;
	}
    }
    if (lp->nrows == 0) {
	rval = mpf_ILLdata_error (raw->error_collector, "There are no constraints.");
    }
    ILL_RESULT (rval, "mpf_whichRowsAreUsed");
}


static int mpf_transferObjective (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *colindex)
{
    int rval = 0, i, ci, objind = raw->objindex;
    mpf_colptr *cp;
    int *coefWarn = NULL;

    /* transfer objective fct */
    lp->obj = mpf_EGlpNumAllocArray (lp->ncols);
    ILL_SAFE_MALLOC (coefWarn, lp->ncols, int);
    for (i = 0; i < lp->ncols; i++) {
	mpf_EGlpNumZero (lp->obj[i]);
	coefWarn[i] = 0;
    }
    for (i = 0; i < raw->ncols; i++) {
	for (cp = raw->cols[i]; cp; cp = cp->next) {
	    if (cp->this == objind) {
		ci = colindex[i];
		TESTG ((rval =
			(ci < 0
			    || ci >= lp->ncols)), CLEANUP, "ci %d is out of range [0,%d[",
		    ci, lp->ncols);
		ILL_FAILfalse (ci != -1,
		    "all vars in obj fct should be marked as useful");
		coefWarn[ci]++;
		if (mpf_EGlpNumIsNeqqZero (cp->coef))
		    mpf_EGlpNumAddTo (lp->obj[ci], cp->coef);
		if (coefWarn[ci] == 2) {
		    mpf_ILLdata_warn (raw->error_collector,
			"Multiple coefficients for \"%s\" in %s.",
			mpf_ILLraw_colname (raw, i), "objective function");
		}
	    }
	}
    }
CLEANUP:
    ILL_IFFREE (coefWarn, int);
    ILL_RETURN (rval, "mpf_transferObjective");
}

static int mpf_transferColNamesLowerUpperIntMarker (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *colindex)
{
    int i, ci, ind, pre, rval = 0;
    int hasIntVar;
    ILL_SAFE_MALLOC (lp->colnames, lp->ncols, char *);
    if (raw->upper)
	lp->upper = mpf_EGlpNumAllocArray (lp->ncols);
    if (raw->lower)
	lp->lower = mpf_EGlpNumAllocArray (lp->ncols);
    ILL_SAFE_MALLOC (lp->intmarker, lp->ncols, char);
    hasIntVar = 0;
    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i];
	if (ci != -1) {
	    ILL_FAILfalse ((ci >= 0) && (ci < lp->ncols), "colindex problem");
	    mpf_ILL_UTIL_STR (lp->colnames[ci], mpf_ILLraw_colname (raw, i));
	    rval = ILLsymboltab_register (&lp->coltab,
		lp->colnames[ci], -1, &ind, &pre);
	    ILL_FAILfalse ((rval == 0) && (ind == ci) && (pre == 0),
		"should have new entry");
	    if (raw->upper) {
		mpf_EGlpNumCopy (lp->upper[ci], raw->upper[i]);
	    }
	    if (raw->lower) {
		mpf_EGlpNumCopy (lp->lower[ci], raw->lower[i]);
	    }
	    lp->intmarker[ci] = raw->intmarker[i];
	    hasIntVar = hasIntVar || lp->intmarker[ci];
	    ILL_IFDOTRACE
	    {
		if (lp->lower) {
		    mpf_ILLprt_EGlpNum (stdout, &(lp->lower[ci]));
		    ILL_IFTRACE (" <= ");
		}
		ILL_IFTRACE ("%s", lp->colnames[ci]);
		if (lp->upper) {
		    ILL_IFTRACE (" <= ");
		    mpf_ILLprt_EGlpNum (stdout, &(lp->upper[ci]));
		}
		if (lp->intmarker[ci]) {
		    ILL_IFTRACE (" INTEGER ");
		}
		ILL_IFTRACE ("\n");
	    }
	}
    }
    if (!hasIntVar) {
	ILL_IFFREE (lp->intmarker, char);
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_transferColNamesLowerUpperIntMarker");
}

static void mpf_safeRegister (ILLsymboltab * tab,
      const char *name,
      int i)
{
    int ind, pre, rval;
    rval = ILLsymboltab_register (tab, name, -1, &ind, &pre);
    ILL_FAILfalse ((rval == 0) && (ind == i) && (pre == 0),
	"Pgming Error: should have new entry");
CLEANUP:
    return;
}

static int mpf_transferSenseRhsRowNames (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *rowindex)
{
    int i, ri, rval = 0;
    int objind = raw->objindex;

    /* transfer sense/rhs/rownames */
    if (lp->nrows > 0) {
	ILL_SAFE_MALLOC (lp->sense, lp->nrows, char);
	lp->rhs = mpf_EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (lp->rownames, lp->nrows, char *);

	ILL_FAILfalse (mpf_ILLraw_rowname (raw, raw->objindex), "NULL objname");
	mpf_safeRegister (&lp->rowtab, mpf_ILLraw_rowname (raw, raw->objindex), 0);

	ri = 0;
	for (i = 0; i < raw->nrows; i++) {
	    ri = rowindex[i];
	    if (i == raw->refrowind) {
		mpf_ILL_UTIL_STR (lp->refrowname, mpf_ILLraw_rowname (raw, i));
		lp->refind = ri;
	    }
	    if (raw->rowsense[i] != 'N') {
		ILL_FAILfalse (mpf_ILLraw_rowname (raw, i) != NULL,
		    "all rownames should be non NULL");
		mpf_ILL_UTIL_STR (lp->rownames[ri], mpf_ILLraw_rowname (raw, i));
		mpf_safeRegister (&lp->rowtab, lp->rownames[ri], ri + 1);
		lp->sense[ri] = raw->rowsense[i];
		mpf_EGlpNumCopy (lp->rhs[ri], raw->rhs[i]);
	    } else if (i == objind) {
		ILL_FAILfalse (lp->objname == NULL, "objname == NULL");
		mpf_ILL_UTIL_STR (lp->objname, mpf_ILLraw_rowname (raw, i));
	    } else {
		/* unused 'N' row */
	    }
	}
	ILL_FAILfalse ((lp->nrows + 1) == lp->rowtab.tablesize,
	    "problem with rowtab structure");
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_transferSenseRhsRowNames");
}

static int mpf_buildMatrix (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *rowindex,
      int *colindex)
{
    int i, ri, ci, k, nempty = 0, rval = 0;
    int *nRowsUsed = 0;
    int *coefSet = 0;
    int *coefWarn = 0;
    mpf_ILLmatrix *A = &lp->A;
    mpf_colptr *cp = NULL;

    /* put subjective fcts into matrix */
    ILL_SAFE_MALLOC (A->matcnt, lp->ncols, int);
    ILL_SAFE_MALLOC (A->matbeg, lp->ncols, int);
    ILL_SAFE_MALLOC (nRowsUsed, lp->nrows, int);

    ILL_SAFE_MALLOC (coefWarn, lp->ncols, int);
    for (i = 0; i < lp->ncols; i++) {
	coefWarn[i] = 0;
    }
    for (i = 0; i < lp->nrows; i++) {
	nRowsUsed[i] = -1;
    }
    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i];
	if (ci == -1)
	    continue;
	k = 0;
	for (cp = raw->cols[i]; cp; cp = cp->next) {
	    ri = rowindex[cp->this];
	    if (ri >= 0) {
		if (nRowsUsed[ri] != i) {
		    nRowsUsed[ri] = i;
		    k++;
		} else {
		    if (!coefWarn[ci]) {
			mpf_ILLdata_warn (raw->error_collector,
			    "Multiple coefficients for \"%s\" %s.",
			    lp->colnames[i], "in a row");
			coefWarn[ci] = 1;
		    }
		}
	    }
	}
	A->matcnt[ci] = k;
	A->matbeg[ci] = lp->nzcount + nempty;	/* mark empty cols */
	lp->nzcount += k;
	if (k == 0)
	    nempty++;
    }

    A->matrows = lp->nrows;
    A->matcols = lp->ncols;
    A->matcolsize = lp->ncols;
    A->matsize = lp->nzcount + nempty + 1;
    A->matfree = 1;
    ILL_SAFE_MALLOC (A->matind, A->matsize, int);
    A->matval = mpf_EGlpNumAllocArray (A->matsize);
    ILL_SAFE_MALLOC (coefSet, lp->nrows, int);

    for (k = 0; k < lp->nrows; k++) {
	coefSet[k] = -1;
    }

    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i];
	if (ci == -1)
	    continue;		/* unused variable */
	k = A->matbeg[ci];
	if (A->matcnt[ci] == 0) {
	    A->matind[k] = 1;	/* Used in addcols and addrows */
	} else {
	    for (cp = raw->cols[i]; cp; cp = cp->next) {
		ri = rowindex[cp->this];
		if (ri >= 0) {
		    if (coefSet[ri] == -1) {
			A->matind[k] = ri;
			mpf_EGlpNumCopy (A->matval[k], cp->coef);
			coefSet[ri] = k;
			k++;
		    } else {
			mpf_EGlpNumAddTo (A->matval[coefSet[ri]], cp->coef);
		    }
		}
	    }
	    if (k != A->matbeg[ci] + A->matcnt[ci]) {
		ILL_ERROR (rval, "problem with matrix");
	    }
	    for (k--; k >= A->matbeg[ci]; k--) {
		coefSet[A->matind[k]] = -1;
	    }
	}
    }
    A->matind[lp->nzcount + nempty] = -1;
CLEANUP:
    ILL_IFFREE (nRowsUsed, int);
    ILL_IFFREE (coefWarn, int);
    ILL_IFFREE (coefSet, int);
    ILL_RETURN (rval, "mpf_buildMatrix");
}

static int mpf_transferRanges (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *rowindex)
{
    int i, ri, rval = 0;
    mpf_colptr *cp;

    /*****************************************************/
    /* */
    /* Interpretation of RANGE values in MPS files      */
    /* */
    /* G    rhs           <= row <= rhs + |range|     */
    /* L    rhs - |range| <= row <= rhs               */
    /* E +  rhs           <= row <= rhs + range       */
    /* E -  rhs + range   <= row <= rhs               */
    /* */
    /* - where + and - refer to the sign of range    */
    /* and the letters refer to sense of the row.  */
    /* */
    /* We will store ranged rows as                   */
    /* */
    /* rhs  <= row  <= rhs + range                 */
    /* */
    /*****************************************************/


    lp->rangeval = mpf_EGlpNumAllocArray (lp->nrows);
    for (i = 0; i < lp->nrows; i++) {
	mpf_EGlpNumZero (lp->rangeval[i]);
    }
    for (cp = raw->ranges; cp; cp = cp->next) {
	i = cp->this;
	ri = rowindex[cp->this];
	switch (raw->rowsense[i]) {
	case 'N':
	    mpf_ILLdata_error (raw->error_collector, "No range for N-row.\n");
	    rval = 1;
	    goto CLEANUP;
	case 'G':
	    lp->sense[ri] = 'R';
	    mpf_EGlpNumCopyAbs (lp->rangeval[ri], cp->coef);
	    break;
	case 'L':
	    lp->sense[ri] = 'R';
	    mpf_EGlpNumCopyAbs (lp->rangeval[ri], cp->coef);
	    mpf_EGlpNumSubTo (lp->rhs[ri], lp->rangeval[ri]);
	    break;
	case 'E':
	    lp->sense[ri] = 'R';
	    if (mpf_EGlpNumIsLeq (mpf_zeroLpNum, cp->coef)) {
		mpf_EGlpNumCopy (lp->rangeval[ri], cp->coef);
	    } else {
		mpf_EGlpNumAddTo (lp->rhs[ri], cp->coef);
		mpf_EGlpNumCopyNeg (lp->rangeval[ri], cp->coef);
	    }
	    break;
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_transferRanges");
}

static int mpf_initStructmap (mpf_ILLlpdata * lp)
{
    int i, rval = 0;

    /* all vars are structural */
    ILL_SAFE_MALLOC (lp->structmap, lp->nstruct, int);
    for (i = 0; i < lp->nstruct; i++) {
	lp->structmap[i] = i;
    }

CLEANUP:
    ILL_RETURN (rval, "mpf_initStructmap");
}

static int mpf_buildSosInfo (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp,
      int *colindex)
{
    int i, ci, set, rval = 0;
    int nSosMem, nSetMem;

    /* build sos info */
    /* see comment in mpf_lpdata.h about mpf_ILLlpdata's sos and is_sos_mem
       fields and section of mpf_ILLprint_rawlpdata that prints SOS sets */

    ILL_SAFE_MALLOC (lp->is_sos_mem, lp->ncols, int);
    nSosMem = 0;
    for (i = 0; i < raw->ncols; i++) {
	ci = colindex[i];
	if (ci != -1) {
	    lp->is_sos_mem[ci] = raw->is_sos_member[i];
	    if (raw->is_sos_member[i] != -1)
		nSosMem++;
	}
    }
    if (nSosMem > 0) {
	lp->sos.matsize = nSosMem;
	lp->sos.matcols = raw->nsos;
	lp->sos.matcolsize = raw->nsos;
	lp->sos.matrows = lp->ncols;
	lp->sos.matfree = 0;
	lp->sos.matval = mpf_EGlpNumAllocArray (nSosMem);
	ILL_SAFE_MALLOC (lp->sos.matind, nSosMem, int);
	ILL_SAFE_MALLOC (lp->sos.matbeg, raw->nsos, int);
	ILL_SAFE_MALLOC (lp->sos.matcnt, raw->nsos, int);
	ILL_SAFE_MALLOC (lp->sos_type, raw->nsos, char);
	nSosMem = 0;
	for (set = 0; set < raw->nsos; set++) {
	    lp->sos_type[set] = raw->sos_set[set].type;
	    lp->sos.matbeg[set] = nSosMem;
	    nSetMem = 0;
	    for (i = raw->sos_set[set].first;
		i < raw->sos_set[set].first + raw->sos_set[set].nelem; i++) {
		ci = colindex[raw->sos_col[i]];
		if (ci != -1) {
		    lp->sos.matind[nSosMem + nSetMem] = ci;
		    mpf_EGlpNumCopy (lp->sos.matval[nSosMem + nSetMem], raw->sos_weight[i]);
		    nSetMem++;
		}
	    }
	    lp->sos.matcnt[set] = nSetMem;
	    nSosMem += nSetMem;
	}
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_buildSosInfo");
}

static int mpf_convert_rawlpdata_to_lpdata (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp)
/* only raw's non 'N' rows are converted to matrix entries in lp columns that
   are used in non objective 'N' rows only are not converted. That is they
   don't end up in lp's matrix, row/colnames, upper/lower bounds or SOS
   information. */
{
    int rval = 0;
    int *rowindex = 0;
    int *colindex = 0;

    ILL_FAILfalse ((raw && lp), "rawlpdata_to_lpdata called without input");
    if (raw->name == NULL) {
	mpf_ILLdata_warn (raw->error_collector, "Setting problem name to \"unnamed\".");
	mpf_ILL_UTIL_STR (raw->name, "unnamed");
    }
    rval = mpf_ILLcheck_rawlpdata (raw);
    ILL_CLEANUP_IF (rval);

    ILL_FAILfalse (raw->objindex != -1, "mpf_rawlpdata must have objective fct.");
    mpf_ILLlpdata_init (lp);

    ILL_IFFREE (lp->probname, char);
    lp->probname = raw->name;
    raw->name = 0;

    /* MINIMIZE or MAXIMIZE ? */
    lp->objsense = raw->objsense;
    if (lp->objsense != mpf_ILL_MIN && lp->objsense != mpf_ILL_MAX) {
	mpf_ILLdata_error (raw->error_collector, "Bad objsense.\n");
	rval = 1;
	goto CLEANUP;
    }
    ILL_SAFE_MALLOC (colindex, raw->ncols, int);
    ILL_SAFE_MALLOC (rowindex, raw->nrows, int);
    rval = mpf_whichColsAreUsed (raw, lp, colindex) ||
	mpf_whichRowsAreUsed (raw, lp, rowindex);
    ILL_CLEANUP_IF (rval);
    ILL_FAILtrue (lp->ncols == 0 || lp->nrows == 0, "we need rows and cols");

    /* array sizes */
    lp->rowsize = lp->nrows;
    lp->colsize = lp->ncols;
    lp->nstruct = lp->ncols;
    lp->structsize = lp->ncols;
    ILLsymboltab_create (&lp->rowtab, lp->nrows);
    ILLsymboltab_create (&lp->coltab, lp->ncols);

    rval = mpf_transferObjective (raw, lp, colindex);
    rval = rval || mpf_transferColNamesLowerUpperIntMarker (raw, lp, colindex);
    rval = rval || mpf_buildMatrix (raw, lp, rowindex, colindex);
    rval = rval || mpf_buildSosInfo (raw, lp, colindex);
    ILL_CLEANUP_IF (rval);
    ILL_IFDOTRACE
    {
	mpf_ILLmatrix_prt (stdout, &lp->A);
    }

    rval = mpf_transferSenseRhsRowNames (raw, lp, rowindex);
    if ((lp->nrows > 0) && raw->ranges) {
	rval = rval || mpf_transferRanges (raw, lp, rowindex);
    }
    ILL_CLEANUP_IF (rval);

    rval = mpf_initStructmap (lp);
    ILL_CLEANUP_IF (rval);

CLEANUP:

    ILL_IFFREE (rowindex, int);
    ILL_IFFREE (colindex, int);
    mpf_ILLfree_rawlpdata (raw);

    ILL_RESULT (rval, "mpf_convert_rawlpdata_to_lpdata");
}

int mpf_ILLrawlpdata_to_lpdata (mpf_rawlpdata * raw,
      mpf_ILLlpdata * lp)
{
    int rval = 0;

    ILL_IFDOTRACE
    {
	printf ("%s\n", __func__);
	mpf_ILLprint_rawlpdata (raw);
    }
    rval = mpf_convert_rawlpdata_to_lpdata (raw, lp);
    if (rval == 0) {
	rval = mpf_ILLlp_add_logicals (lp);
    }
    ILL_RESULT (rval, "mpf_ILLrawlpdata_to_lpdata");
}

static int mpf_set_field_name (char **field,
      const char *name,
      int *skip)
{
    int rval = 0;
    /* name is bounds/rhs/rangesname field from mpf_rawlpdata */
    *skip = 0;
    if (!*field) {
	mpf_ILL_UTIL_STR (*field, name);
    }
    if (strcmp (*field, name)) {
	/* not first specified RHS/BOUNDS - skip it */
	*skip = 1;
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_set_field_name");
}

int mpf_ILLraw_set_rhs_name (mpf_rawlpdata * lp,
      const char *name,
      int *skip)
{
    return mpf_set_field_name (&lp->rhsname, name, skip);
}

int mpf_ILLraw_set_bounds_name (mpf_rawlpdata * lp,
      const char *name,
      int *skip)
{
    return mpf_set_field_name (&lp->boundsname, name, skip);
}

int mpf_ILLraw_set_ranges_name (mpf_rawlpdata * lp,
      const char *name,
      int *skip)
{
    return mpf_set_field_name (&lp->rangesname, name, skip);
}

void mpf_ILLprint_rawlpdata (mpf_rawlpdata * lp)
{
    int i, cnt, si, m;
    char c;
    mpf_t d;
    mpf_colptr *cp;
    mpf_sosptr *set;
    mpf_EGlpNumInitVar (d);

    if (lp) {
	if (lp->name) {
	    printf ("PROBLEM  %s\n", lp->name);
	    fflush (stdout);
	}
	if (lp->rowsense && lp->rhs) {
	    printf ("Subject To\n");
	    for (i = 0; i < lp->nrows; i++) {
		switch (lp->rowsense[i]) {
		case 'E':
		    c = '=';
		    break;
		case 'L':
		    c = '<';
		    break;
		case 'G':
		    c = '>';
		    break;
		default:
		    c = '?';
		    break;
		}
		printf ("%s: %c %f\n", mpf_ILLraw_rowname (lp, i), c,
		    mpf_EGlpNumToLf (lp->rhs[i]));
	    }
	    printf ("\n");
	    fflush (stdout);
	}
	if (lp->ncols > 0) {
	    printf ("Columns\n");
	    for (i = 0; i < lp->ncols; i++) {
		for (cp = lp->cols[i]; cp; cp = cp->next) {
		    printf ("%s: ", mpf_ILLraw_rowname (lp, cp->this));
		    printf ("%c ", (mpf_EGlpNumIsLess (cp->coef, mpf_zeroLpNum)) ? '-' : '+');
		    mpf_EGlpNumCopyAbs (d, cp->coef);
		    if (mpf_EGlpNumIsNeqq (d, mpf_oneLpNum)) {
			printf (" %f ", mpf_EGlpNumToLf (d));
		    }
		    printf ("%s\n", mpf_ILLraw_colname (lp, i));
		}
		printf ("\n");
		fflush (stdout);
	    }
	}
	if (lp->rangesname) {
	    printf ("RANGES %s\n", lp->rangesname);
	    for (cp = lp->ranges; cp; cp = cp->next) {
		printf ("(%s, %f) ", mpf_ILLraw_rowname (lp, cp->this),
		    mpf_EGlpNumToLf (cp->coef));
	    }
	    printf ("\n");
	    fflush (stdout);
	}
	if (lp->boundsname) {
	    printf ("BOUNDS %s\n", lp->boundsname);
	    fflush (stdout);
	} else {
	    printf ("BOUNDS \n");
	    fflush (stdout);
	}
	if (lp->lower && lp->upper) {
	    for (i = 0; i < lp->ncols; i++) {
		mpf_ILLprt_EGlpNum (stdout, &(lp->lower[i]));
		printf (" <= %s <= ", mpf_ILLraw_colname (lp, i));
		mpf_ILLprt_EGlpNum (stdout, &(lp->upper[i]));
		printf ("\n");
	    }
	}
	if (lp->intmarker) {
	    printf ("Integer\n");
	    cnt = 0;
	    for (i = 0; i < lp->ncols; i++) {
		if (lp->intmarker[i]) {
		    printf ("%s", mpf_ILLraw_colname (lp, i));
		    cnt++;
		    if (cnt == 8) {
			printf ("\n    ");
			cnt = 0;
		    }
		}
	    }
	    printf ("\n");
	    fflush (stdout);
	}
	printf ("SOS-SETS\n");
	for (si = 0; si < lp->nsos; si++) {
	    set = lp->sos_set + si;
	    printf ("SOS-SET %d: %s; nelem=%d; first=%d;\n",
		si, ((set->type == mpf_ILL_SOS_TYPE1) ? "TYPE1" : "TYPE2"),
		set->nelem, set->first);
	    printf ("\t");
	    for (m = set->first; m < set->first + set->nelem; m++) {
		printf (" %s %f; ", mpf_ILLraw_colname (lp, lp->sos_col[m]),
		    mpf_EGlpNumToLf (lp->sos_weight[m]));
	    }
	    printf ("\n");
	}
	printf ("\n");
	fflush (stdout);
    }
    mpf_EGlpNumClearVar (d);
}

static int mpf_ILLmsg (mpf_qserror_collector * collector,
      int isError,
      const char *format,
      va_list args)
{
    const char *pre;
    int slen, errtype;
    mpf_qsformat_error error;
    char error_desc[256];

    vsprintf (error_desc, format, args);
    slen = strlen (error_desc);
    if ((slen > 0) && error_desc[slen - 1] != '\n') {
	error_desc[slen] = '\n';
	error_desc[slen + 1] = '\0';
    }
    if (collector != NULL) {
	errtype = (isError) ? QS_DATA_ERROR : QS_DATA_WARN;
	mpf_ILLformat_error_create (&error, errtype, error_desc, -1, NULL, -1);
	mpf_ILLformat_error (collector, &error);
	mpf_ILLformat_error_delete (&error);
    } else {
	pre = (isError) ? "Data Error" : "Data Warning";
	fprintf (stderr, "%s: %s", pre, error_desc);
    }
    return 1;
}

int mpf_ILLdata_error (mpf_qserror_collector * collector,
      const char *format,
    ...)
{
    va_list args;
    va_start (args, format);
    return mpf_ILLmsg (collector, mpf_TRUE, format, args);
}

void mpf_ILLdata_warn (mpf_qserror_collector * collector,
      const char *format,
    ...)
{
    va_list args;
    va_start (args, format);
    (void) mpf_ILLmsg (collector, mpf_FALSE, format, args);
}

mpf_colptr *mpf_ILLcolptralloc (ILLptrworld * p)
{
    mpf_colptr *sol = colptralloc (p);
    mpf_EGlpNumInitVar ((sol->coef));
    return sol;
}
