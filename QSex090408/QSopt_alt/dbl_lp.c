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

/* RCS_INFO = "$RCSfile: dbl_lp.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

/****************************************************************************/
/* */
/* Routines for Reading and Writing LP Files                  */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int dbl_ILLwrite_lp (FILE *out, dbl_ILLlpdata *lp)                        */
/* int dbl_ILLread_lp (dbl_qsline_reader f, const char *dbl_fname,
   dbl_rawlpdata *lp) */
/* int dbl_ILLis_lp_name_char (char c, int pos)                          */
/* int dbl_ILLread_constraint_expr(dbl_ILLread_lp_state *state,              */
/* dbl_rawlpdata *lp, int rowind, int allowNew)                        */
/* int dbl_ILLread_constraint_name (dbl_ILLread_lp_state *state,             */
/* char **rowname)                                                 */
/* int dbl_ILLread_one_constraint (dbl_ILLread_lp_state *state,              */
/* const char *rowname, dbl_rawlpdata *lp, int allowNewCols)           */
/* */
/****************************************************************************/

#include <stdio.h>
#include "econfig.h"
/* from the cplex manual:

not exceed 16 characters, all of which must be alphanumeric (a-z, A-Z, 0-9)
   or one of these symbols: ! " # $ % & ( ) / , . ; ? @ _ ` ' { } | ~. Longer
   names will be truncated to 16 characters. A variable name can not begin
   with a number or a period.

The letter E or e, alone or followed by other valid symbols, or followed by
   another E or e, should be avoided as this notation is reserved for
   exponential entries. Thus, variables can not be named e9, E-24, E8cats, or
   other names that could be interpreted as an exponent. Even variable names
   such as eels or example can cause a read error, depending on their
   placement in an input line.

Also, the following characters are not valid in variable names (in order to
   allow for quadratic objective information): ^, *, [ and ]. */
/* OUR DEFINTION: -) variables consist of a-z A-Z 0-9 !"#$%(),;.?@_`'{}|~
   don't start with a digit or '.' return  0 for variable return -1 for
   keyword return  1 otherwise */


int dbl_ILLis_lp_name_char (int c,
      int pos)
{
    return ((('a' <= c) && (c <= 'z')) ||
	(('A' <= c) && (c <= 'Z')) ||
	((pos > 0) && ('0' <= c) && (c <= '9')) ||
	((pos > 0) && (c == '.')) ||
	(strchr ("!\"#$%&()/,;?@_`'{}|~", c) != NULL));
}


/* -) anything after '\' is comment -) Problem is optional and comes first -)
   Minimize, Maximize, ... comes after Problem -) variables consist of a-z
   A-Z 0-9!"#$%(),;.?@_`'{}|~ */

/* NOTES                                                                   */
/* don't start with a digit or '.' */

static const int dbl_LINE_LEN = 64;

#include "dbl_iqsutil.h"
#include "dbl_lp.h"
#include "dbl_rawlp.h"
#include "dbl_read_lp.h"
#include "dbl_write_lp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
/* extern double dbl_SZERO_TOLER; */
static int TRACE = 0;

static int dbl_read_problem_name (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp);
static int dbl_read_minmax (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp);
static int dbl_read_objective (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp);
static int dbl_read_objective (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp);
static int dbl_read_constraints (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp,
      int allowNewCols);
static int dbl_read_colname (dbl_ILLread_lp_state * state,
      ILLsymboltab * coltab,
      int mustHave);
static int dbl_read_integer (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp);
static int dbl_read_bounds (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp);
static int dbl_add_var (dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state,
      double coef,
      int row,
      int allowNew);

/*------------------------------------------------------------------------
 * dbl_ILLwrite_lp and support routines
 */
static int dbl_fix_names ( /* dbl_qserror_collector * collector, */ char **names,
      int nnames, const char *extra, int prefix, char ***newnames);
static void dbl_write_objective (dbl_ILLlpdata * lp, const char *objname,
      char **colnames);
static int dbl_write_row (dbl_ILLlpdata * lp, dbl_ILLlp_rows * lprows,
      int i, char **rownames, char **colnames, int *colInRow, double *colCoef);
static int dbl_write_bounds (dbl_ILLlpdata * lp, char **colnames);
static void dbl_write_intvars (dbl_ILLlpdata * lp, char **colnames);

int dbl_ILLwrite_lp (dbl_ILLlpdata * lp, dbl_qserror_collector * collector)
{
    int rval = 0;
    int i;
    dbl_ILLlp_rows lp_rows, *lprows = NULL;
    char **colnames = (char **) NULL;
    char **rownames = (char **) NULL;
    double *colCoef = NULL;
    int *colInRow = NULL;
    const char *objname;

    ILL_FAILfalse (lp, "called without data\n");
    if (lp->nstruct == 0 || lp->nrows == 0) {
	ILL_RETURN (rval, "dbl_ILLwrite_lp");
    }
    ILL_FAILfalse (lp->colnames != NULL, "lp->colnames != NULL");
    ILL_FAILfalse (lp->rownames != NULL, "lp->rownames != NULL");
    ILL_FAILfalse (lp->nstruct == lp->coltab.tablesize,
	"lp coltab has nstruct entries");
    if (lp->objname == (char *) NULL) {
	ILL_FAILfalse (lp->nrows == lp->rowtab.tablesize,
	    "lp rowtab should have nrows entries");
    } else {
	ILL_FAILfalse (lp->nrows + 1 == lp->rowtab.tablesize,
	    "lp rowtab should have nrows+1 entries");
	ILL_FAILfalse (ILLsymboltab_contains (&lp->rowtab, lp->objname),
	    "rowtab must contain objname");
    }

    rval = dbl_fix_names (lp->colnames, lp->nstruct, NULL, 'x', &colnames);
    ILL_CLEANUP_IF (rval);

    rval = dbl_fix_names (lp->rownames, lp->nrows,
	(lp->objname) ? lp->objname : "obj", 'c', &rownames);
    ILL_CLEANUP_IF (rval);
    objname = rownames[lp->nrows];

    ILL_FAILtrue (objname == NULL, "OOps, that should never happen");
    ILL_CLEANUP_IF (rval);

    if (lp->sos.matcols > 0) {
	rval +=
	    dbl_ILLdata_error (collector, "Can't express SOS information in LP format.");
    }
    dbl_write_objective (lp, objname, colnames);

    /* Note, dbl_ILLlp_rows_init returns cols ordered by structmap, so we may
       use colnames[i] when pulling i from the matrix data.  */

    lprows = &lp_rows;
    if (dbl_ILLlp_rows_init (lprows, lp, 0) != 0) {
	rval += 1;
	ILL_FAILtrue (rval, "dbl_ILLlp_rows_init failed\n");
    }
    colCoef = dbl_EGlpNumAllocArray (lp->nstruct);
    ILL_SAFE_MALLOC (colInRow, lp->nstruct, int);
    for (i = 0; i < lp->nstruct; i++) {
	colInRow[i] = -1;
    }

    dbl_ILLprint_report (lp, "Subject To\n");
    for (i = 0; i < lp->nrows; i++) {
	if (lprows->rowcnt[i] == 0) {
	    /*
             * dbl_ILLdata_warn (collector, "Not printing  empty row \"%s\".", rownames[i]);
             */
	    continue;
	}
	rval += dbl_write_row (lp, lprows, i, rownames, colnames, colInRow, colCoef);
    }

    rval += dbl_write_bounds (lp, colnames);

    if (lp->intmarker != NULL) {
	dbl_write_intvars (lp, colnames);
    }
    dbl_ILLprint_report (lp, "End\n");
CLEANUP:
    if (lprows != NULL) {
	dbl_ILLlp_rows_clear (lprows);
    }
    ILLfree_names (colnames, lp->nstruct);
    ILLfree_names (rownames, lp->nrows + 1);
    dbl_EGlpNumFreeArray (colCoef);
    ILL_IFFREE (colInRow, int);
    ILL_RETURN (rval, "dbl_ILLwrite_lp");
}

static void dbl_write_objective (dbl_ILLlpdata * lp,
      const char *objname,
      char **colnames)
{
    int ri, i, k, var;
    dbl_ILLwrite_lp_state ln, *line = &ln;

    if (lp->probname != NULL) {
	dbl_ILLprint_report (lp, "Problem\n %s\n", lp->probname);
    }
    if (lp->objsense == dbl_ILL_MIN) {
	dbl_ILLprint_report (lp, "Minimize\n");
    } else {
	dbl_ILLprint_report (lp, "Maximize\n");
    }
    dbl_ILLwrite_lp_state_init (line, NULL);
    dbl_ILLwrite_lp_state_append (line, " ");
    dbl_ILLwrite_lp_state_append (line, objname);
    dbl_ILLwrite_lp_state_append (line, ": ");
    dbl_ILLwrite_lp_state_save_start (line);

    for (ri = 0, var = 0; ri < lp->nstruct; ri++) {
	i = lp->structmap[ri];
	if (dbl_EGlpNumIsNeqqZero (lp->obj[i])) {
	    dbl_ILLwrite_lp_state_append_coef (line, lp->obj[i], var);
	    dbl_ILLwrite_lp_state_append (line, " ");
	    dbl_ILLwrite_lp_state_append (line, colnames[ri]);
	    var++;

	    /* we put a least 4 terms on a line and then we stop after
	       dbl_LINE_LEN or more characters */
	    if ((line->total >= dbl_LINE_LEN) && (var >= 4)) {
		/* see whether there is another term if so append a '+' and
		   print line */
		k = ri + 1;
		while (k < lp->nstruct) {
		    if (dbl_EGlpNumIsLess (lp->obj[lp->structmap[k]], dbl_zeroLpNum)) {
			break;
		    } else {
			if (dbl_EGlpNumIsLess (dbl_zeroLpNum, lp->obj[lp->structmap[k]])) {
			    dbl_ILLwrite_lp_state_append (line, " +");
			    break;
			}
		    }
		    k++;
		}
		var = 0;	/* next line does not need to prefix coef
				   with '+' */
		dbl_ILLprint_report (lp, "%s\n", line->buf);
		dbl_ILLwrite_lp_state_start (line);
	    }
	}
    }
    if (var > 0) {
	dbl_ILLprint_report (lp, "%s\n", line->buf);
    }
}

static void dbl_write_the_expr (dbl_ILLlpdata * lp,
      dbl_ILLwrite_lp_state * line,
      char *rowname,
      dbl_ILLlp_rows * lprows,
      int row,
      char **colnames,
      int *colInRow,
      double *colCoef,
      int ncols)
{
    int var, firstVar, k, i;
    double *coef;

    dbl_ILLwrite_lp_state_init (line, NULL);
    if (rowname != NULL) {
	dbl_ILLwrite_lp_state_append (line, " ");
	dbl_ILLwrite_lp_state_append (line, rowname);
	dbl_ILLwrite_lp_state_append (line, ": ");
    } else {
	dbl_ILLwrite_lp_state_append (line, "   ");
    }
    dbl_ILLwrite_lp_state_save_start (line);

    for (k = lprows->rowbeg[row];
	k < lprows->rowbeg[row] + lprows->rowcnt[row]; k++) {
	i = lprows->rowind[k];
	colInRow[i] = row;
	dbl_EGlpNumCopy (colCoef[i], lprows->rowval[k]);
    }
    var = 0;
    firstVar = 1;
    for (i = 0; i < ncols; i++) {
	if (colInRow[i] == row) {
	    if (dbl_EGlpNumIsNeqqZero (colCoef[i])) {
		coef = &(colCoef[i]);
		if (line->total >= dbl_LINE_LEN) {
		    dbl_ILLprint_report (lp, "%s\n", line->buf);
		    dbl_ILLwrite_lp_state_start (line);
		    if ((!firstVar) && dbl_EGlpNumIsLeq (dbl_zeroLpNum, *coef)) {
			dbl_ILLwrite_lp_state_append (line, " +");
		    }
		    var = 0;
		}
		dbl_ILLwrite_lp_state_append_coef (line, *coef, var);
		dbl_ILLwrite_lp_state_append (line, " ");
		dbl_ILLwrite_lp_state_append (line, colnames[i]);
		var++;
		firstVar = 0;
	    }
	}
    }
}

static int dbl_write_row (dbl_ILLlpdata * lp,
      dbl_ILLlp_rows * lprows,
      int i,
      char **rownames,
      char **colnames,
      int *colInRow,
      double *colCoef)
{
    dbl_ILLwrite_lp_state ln, *line = &ln;
    int rval = 0;
    double ntmp;

    dbl_write_the_expr (lp, line, rownames[i], lprows, i, colnames,
	colInRow, colCoef, lp->nstruct);

    switch (lp->sense[i]) {
    case 'G':
	dbl_ILLwrite_lp_state_append (line, " >= ");
	dbl_ILLwrite_lp_state_append_number (line, lp->rhs[i]);
	break;
    case 'L':
	dbl_ILLwrite_lp_state_append (line, " <= ");
	dbl_ILLwrite_lp_state_append_number (line, lp->rhs[i]);
	break;
    case 'E':
	dbl_ILLwrite_lp_state_append (line, " = ");
	dbl_ILLwrite_lp_state_append_number (line, lp->rhs[i]);
	break;
    case 'R':
	ILL_FAILtrue (!lp->rangeval, "RANGE constraints without values\n");
	dbl_EGlpNumInitVar (ntmp);
	dbl_ILLwrite_lp_state_append (line, " >= ");
	dbl_ILLwrite_lp_state_append_number (line, lp->rhs[i]);

	dbl_ILLwrite_lp_state_append (line, " \t\\ RANGE (");
	dbl_ILLwrite_lp_state_append_number (line, lp->rhs[i]);
	dbl_ILLwrite_lp_state_append (line, ", ");
	dbl_EGlpNumCopySum (ntmp, lp->rhs[i], lp->rangeval[i]);
	dbl_ILLwrite_lp_state_append_number (line, ntmp);
	dbl_ILLwrite_lp_state_append (line, ")");
	dbl_ILLprint_report (lp, "%s\n", line->buf);

	dbl_write_the_expr (lp, line, NULL, lprows, i,
	    colnames, colInRow, colCoef, lp->nstruct);
	dbl_ILLwrite_lp_state_append (line, " <= ");
	dbl_ILLwrite_lp_state_append_number (line, ntmp);
	dbl_EGlpNumClearVar (ntmp);
	break;
    default:
	ILL_FAILtrue (1, "Unknown row sense\n");
    }

    dbl_ILLprint_report (lp, "%s\n", line->buf);
CLEANUP:
    ILL_RETURN (rval, "dbl_write_row");
}

static int dbl_write_bounds (dbl_ILLlpdata * lp,
      char **colnames)
{
    int ri, i, rval = 0;
    int prtLower, prtUpper;
    dbl_ILLwrite_lp_state l, *line = &l;

    ILL_FAILtrue (lp->lower == NULL || lp->upper == NULL,
	"Should not call dbl_write_bounds when lower or upper are NULL");
    ri = dbl_ILLraw_first_nondefault_bound (lp);
    if (ri != lp->nstruct) {
	dbl_ILLprint_report (lp, "Bounds\n");
	dbl_ILLwrite_lp_state_init (line, " ");
	dbl_ILLwrite_lp_state_save_start (line);

	for (ri = ri; ri < lp->nstruct; ri++) {
	    dbl_ILLwrite_lp_state_start (line);
	    i = lp->structmap[ri];
	    if (dbl_EGlpNumIsEqqual (lp->lower[i], lp->upper[i])) {
		dbl_ILLwrite_lp_state_append (line, " ");
		dbl_ILLwrite_lp_state_append (line, colnames[ri]);
		dbl_ILLwrite_lp_state_append (line, " = ");
		dbl_ILLwrite_lp_state_append_number (line, lp->upper[i]);
		dbl_ILLprint_report (lp, "%s\n", line->buf);
		continue;
	    }
	    if ((dbl_EGlpNumIsEqqual (lp->lower[i], dbl_ILL_MINDOUBLE)) &&
		(dbl_EGlpNumIsEqqual (lp->upper[i], dbl_ILL_MAXDOUBLE))) {
		dbl_ILLwrite_lp_state_append (line, colnames[ri]);
		dbl_ILLwrite_lp_state_append (line, " free");
		dbl_ILLprint_report (lp, "%s\n", line->buf);
		continue;
	    }
	    prtLower = !dbl_ILLraw_default_lower (lp, i);
	    prtUpper = !dbl_ILLraw_default_upper (lp, i);
	    if (prtLower || prtUpper) {
		if (prtLower) {
		    dbl_ILLwrite_lp_state_append_number (line, lp->lower[i]);
		    dbl_ILLwrite_lp_state_append (line, " <= ");
		}
		if (prtLower || prtUpper) {
		    dbl_ILLwrite_lp_state_append (line, colnames[ri]);
		}
		if (prtUpper) {
		    dbl_ILLwrite_lp_state_append (line, " <= ");
		    dbl_ILLwrite_lp_state_append_number (line, lp->upper[i]);
		}
		dbl_ILLprint_report (lp, "%s\n", line->buf);
	    }
	}
    }
CLEANUP:
    ILL_RETURN (rval, "dbl_write_bounds");
}

static void dbl_write_intvars (dbl_ILLlpdata * lp,
      char **colnames)
{
    dbl_ILLwrite_lp_state ln, *line = &ln;
    int var, j;

    dbl_ILLprint_report (lp, "Integer\n");
    dbl_ILLwrite_lp_state_init (line, " ");
    dbl_ILLwrite_lp_state_save_start (line);

    for (j = 0, var = 0; j < lp->nstruct; j++) {
	if (lp->intmarker[j]) {
	    if (var > 0) {
		dbl_ILLwrite_lp_state_append (line, " ");
	    }
	    dbl_ILLwrite_lp_state_append (line, colnames[j]);
	    var++;
	    if (line->total >= dbl_LINE_LEN) {
		dbl_ILLprint_report (lp, "%s\n", line->buf);
		dbl_ILLwrite_lp_state_init (line, " ");
		var = 0;
	    }
	}
    }
    if (var > 0) {
	dbl_ILLprint_report (lp, "%s\n", line->buf);
    }
}

/* ------------------------------------------------------------ */
/* fix up names that are numbers, i.e. prefix with x X x_ or X_ */

/* redefine names that start with [0-9]; i.e. give a prefix of "x" | "X" |
   "x_" | "X_"  or if all these are already taken prefix with "X"<number>
   make sure names contain solely the characters: [a-zA-Z0-9] and ! " # $ % &
   ( ) / , . ; ? @ _ ` ' { } | ~ rename names with 'bad' chars to <x X x_
   X_><number> */
static int dbl_fix_names ( /* dbl_qserror_collector * collector, */ char **names,
      int nnames, const char *extra, int pref, char ***newnames)
{
    ILLsymboltab symt, *symtab = NULL;
    int rval = 0, i, j, n, ind, hit;
    char **n_names = NULL;
    const char *old_name;
    char buf[ILL_namebufsize];
    char p1[2], p2[3];
    p1[0] = pref;
    p1[1] = '\0';
    p2[0] = pref;
    p2[1] = '_';
    p2[2] = '\0';

    ILL_SAFE_MALLOC (n_names, nnames + 1, char *);
    for (i = 0; i < nnames; i++) {
	n_names[i] = (char *) NULL;
    }

    for (i = 0; i <= nnames; i++) {
	if (i == nnames) {
	    if (extra == NULL)
		break;
	    old_name = extra;
	} else
	    old_name = names[i];

	n = strlen (old_name);
	strcpy (buf, old_name);
	if (!dbl_ILLis_lp_name_char (buf[0], 1)) {
	    sprintf (buf, "%d", i);
	} else {
	    for (j = 1; j < n; j++) {
		if (!dbl_ILLis_lp_name_char (buf[j], j)) {
		    sprintf (buf, "%d", i);
		    break;
		}
	    }
	}

	if (!dbl_ILLis_lp_name_char (buf[0], 0)) {
	    if (symtab == NULL) {
		symtab = &symt;
		ILLsymboltab_init (symtab);
		ILLsymboltab_create (symtab, nnames + 1);
		for (j = 0; j < nnames; j++) {
		    ILLsymboltab_register (symtab, names[j], -1, &ind, &hit);
		    ILL_FAILfalse (ind == j, "ind == j");
		}
		if (extra != NULL)
		    ILLsymboltab_register (symtab, extra, -1, &ind, &hit);
	    }
	    rval = ILLsymboltab_uname (symtab, buf, p1, p2);
	    ILL_CLEANUP_IF (rval);
	    rval = ILLsymboltab_rename (symtab, i, buf);
	    ILL_CLEANUP_IF (rval);

	    dbl_ILL_UTIL_STR (n_names[i], buf);
	    /*
             * dbl_ILLdata_warn (collector,
             * "\"%s\" is not a valid name in LP format; %s\"%s\".",
             * old_name, "renaiming to ", buf);
             */
	} else {
	    dbl_ILL_UTIL_STR (n_names[i], old_name);
	}
    }

CLEANUP:
    if (symtab != NULL) {
	ILLsymboltab_free (symtab);
    }
    *newnames = n_names;
    ILL_RETURN (rval, "dbl_fix_names");
}

/* end ILLlpdata_lpwrite
   ---------------------------------------------------------------------- */

int dbl_ILLread_lp (dbl_qsline_reader * file,
      const char *dbl_fname,
      dbl_rawlpdata * lp)
{
    /* file format: optional problem name (this is in addition to the bix
       format) min or max objective fct constraints (mandatory !) optional
       bounds optional integer (also new)
    
    as opposed to the official bix-format: no blanks in variable names there
       may be several bound defs on the same line; bound definitions may
       cross line boundaries constraints that have no name get generated
       names: if constraint i has not name we take the first name from the
       following list that is not in use: c<i>,  C<i>, or c<i>_0, c<i>_1,
       c<i>_2, ... */
    int rval = 0;
    dbl_ILLread_lp_state lpstate, *state = &lpstate;

    const char *bnds[3], *integer[3], *end[2];
    bnds[0] = "BOUNDS";
    bnds[1] = "BOUND";
    bnds[2] = NULL;
    integer[0] = "INTEGER";
    integer[1] = "INT";
    integer[2] = NULL;
    end[0] = "END";
    end[1] = NULL;

    rval = dbl_ILLread_lp_state_init (state, file, dbl_fname, 0);
    ILL_CLEANUP_IF (rval);

    dbl_ILLinit_rawlpdata (lp, file->error_collector);
    rval = ILLsymboltab_create (&lp->rowtab, 100) ||
	ILLsymboltab_create (&lp->coltab, 100);
    ILL_CLEANUP_IF (rval);

    if (dbl_ILLread_lp_state_next_field (state)) {
	rval = dbl_ILLlp_error (state, "Empty file.\n");
    }
    if (rval == 0)
	rval = dbl_read_problem_name (state, lp);
    if (rval == 0)
	rval = dbl_read_minmax (state, lp);
    if (rval == 0)
	rval = dbl_read_objective (state, lp);
    if (rval == 0)
	rval = dbl_read_constraints (state, lp, 1);
    if ((rval == 0) && (lp->ncols == 0 || lp->nrows == 0)) {
	rval = dbl_ILLlp_error (state,
	    "Problem must contain at least one %s.\n",
	    "non empty constraint");
    }
    ILL_CLEANUP_IF (rval);

    if (dbl_ILLread_lp_state_keyword (state, bnds) == 0) {
	rval = dbl_read_bounds (state, lp);
    }
    ILL_CLEANUP_IF (rval);

    if (dbl_ILLread_lp_state_keyword (state, integer) == 0) {
	rval = dbl_read_integer (state, lp);
    }
    ILL_CLEANUP_IF (rval);

    rval = dbl_ILLread_lp_state_keyword (state, end);
    if (rval != 0) {
	if (state->eof) {
	    rval = dbl_ILLlp_error (state, "Missing \"End\" at end of file.\n");
	} else {
	    rval = dbl_ILLlp_error (state, "\"%s\" unknown keyword\n", state->field);
	}
    }
    if (rval == 0) {
	rval = dbl_ILLraw_fill_in_rownames (lp) || dbl_ILLraw_fill_in_bounds (lp);
    }
CLEANUP:
    dbl_EGlpNumClearVar (lpstate.bound_val);
    ILL_RESULT (rval, "read_lp");
}

static int dbl_read_problem_name (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp)
{
    int rval = 0;

    if (!state->fieldOnFirstCol) {
	rval = dbl_ILLlp_error (state,
	    "Keyword \"%s\" not at beginning of line.\n",
	    state->field);
    }
    if (!dbl_ILLutil_strcasecmp (state->field, "PROBLEM") ||
	!dbl_ILLutil_strcasecmp (state->field, "PROB")) {
	if (dbl_ILLread_lp_state_next_field (state) != 0) {
	    rval = dbl_ILLlp_error (state, "No Problem name field.\n");
	} else {
	    ILL_IFFREE (lp->name, char);
	    dbl_ILL_UTIL_STR (lp->name, state->field);
	    ILL_IFTRACE ("ProblemName: %s\n", state->field);
	    (void) dbl_ILLread_lp_state_next_field (state);
	}
    }
CLEANUP:
    ILL_RESULT (rval, "dbl_read_problem_name");
}

static int dbl_read_minmax (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp)
{
    int rval = 0;
    if (!state->fieldOnFirstCol) {
	rval = dbl_ILLlp_error (state,
	    "Keyword \"%s\" not at beginning of line.\n",
	    state->field);
    }
    if (!dbl_ILLutil_strcasecmp (state->field, "MAX") ||
	!dbl_ILLutil_strcasecmp (state->field, "MAXIMUM") ||
	!dbl_ILLutil_strcasecmp (state->field, "MAXIMIZE")) {
	lp->objsense = dbl_ILL_MAX;
    } else {
	if (!dbl_ILLutil_strcasecmp (state->field, "MIN") ||
	    !dbl_ILLutil_strcasecmp (state->field, "MINIMUM") ||
	    !dbl_ILLutil_strcasecmp (state->field, "MINIMIZE")) {
	    lp->objsense = dbl_ILL_MIN;
	} else {
	    dbl_ILLread_lp_state_prev_field (state);
	    rval = dbl_ILLlp_error (state, "Expecting \"%s\" or \"%s\" keyword.\n",
		"Minimize", "Maximize");
	}
    }
    ILL_RESULT (rval, "dbl_read_minmax");
}

int dbl_ILLread_constraint_expr (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp,
      int rowind,
      int allowNew)
{
    int rval = 0;
    char firstTerm, haveCoef;
    const char *name;
    double sign, coef;
    double ntmp;
    dbl_EGlpNumInitVar (ntmp);
    dbl_EGlpNumInitVar (sign);
    dbl_EGlpNumInitVar (coef);

    firstTerm = 1;
    while (1) {
	if (dbl_ILLread_lp_state_sign (state, &sign) != 0) {
	    if (!firstTerm) {
		break;		/* we've ssen at least one term, this is the
				   constraint's end */
	    }
	}
	haveCoef = dbl_ILLread_lp_state_possible_coef (state, &coef, dbl_oneLpNum);
	if (dbl_ILLread_lp_state_next_var (state) == 0) {
	    dbl_EGlpNumCopy (ntmp, coef);
	    dbl_EGlpNumMultTo (ntmp, sign);
	    rval = dbl_add_var (lp, state, ntmp, rowind, allowNew);
	    ILL_CLEANUP_IF (rval);
	} else {
	    if (haveCoef == 0) {
		return dbl_ILLlp_error (state, "Coefficient without variable.\n");
	    } else {
		break;
	    }
	}
	firstTerm = 0;
    }
CLEANUP:
    if ((rval == 0) && firstTerm) {
	name = dbl_ILLraw_rowname (lp, rowind);
	if (name != NULL) {
	    dbl_ILLlp_warn (state,
		"No terms in constraint expression for \"%s\".\n", name);
	} else {
	    dbl_ILLlp_warn (state, "No terms in constraint expression.\n");
	}
    }
    dbl_EGlpNumClearVar (ntmp);
    dbl_EGlpNumClearVar (sign);
    dbl_EGlpNumClearVar (coef);
    ILL_RESULT (rval, "dbl_ILLread_constraint_expr");
}

static int dbl_read_objective (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp)
{
    int rval = 0;
    char objname[ILL_namebufsize];
    char *name;

    ILL_FAILfalse (lp->nrows == 0, "objective should be first row");
    dbl_ILLread_lp_state_skip_blanks (state, 1);
    if (dbl_ILLread_lp_state_has_colon (state)) {
	if (dbl_ILLread_lp_state_next_var (state) != 0) {
	    rval = dbl_ILLlp_error (state, "Bad objective function name.\n");
	}
	name = state->field;
	if (rval == 0) {
	    if (dbl_ILLread_lp_state_colon (state) != 0) {
		rval = dbl_ILLlp_error (state, "':' must follow constraint row name.\n");
	    }
	}
    } else {
	name = NULL;
    }

    if (rval == 0) {
	ILL_FAILfalse (lp->rowtab.tablesize == 0,
	    "objective row is first in symbol tab");
	if (name == NULL) {
	    strcpy (objname, "obj");
	    dbl_ILLlp_warn (state, "Empty obj name; using \"%s\".\n", objname);
	} else {
	    strcpy (objname, name);
	}
	rval = dbl_ILLraw_add_row (lp, objname, 'N', dbl_zeroLpNum);
	lp->objindex = lp->nrows - 1;
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLread_constraint_expr (state, lp, lp->objindex, 1);
    }
CLEANUP:
    ILL_RESULT (rval, "dbl_read_objective");
}

int dbl_ILLread_constraint_name (dbl_ILLread_lp_state * state,
      char **rowname)
{
    int rval = 0;
    *rowname = NULL;

    /* if there is a ':' on the line: look for constraint row name */
    if (dbl_ILLread_lp_state_has_colon (state)) {
	if (dbl_ILLread_lp_state_next_var (state) != 0) {
	    rval = dbl_ILLlp_error (state, "Bad constraint row name.\n");
	} else {
	    *rowname = state->field;
	    if (dbl_ILLread_lp_state_colon (state) != 0) {
		rval = dbl_ILLlp_error (state, "':' must follow constraint row name.\n");
	    }
	}
    }
    return rval;
}

int dbl_ILLread_one_constraint (dbl_ILLread_lp_state * state,
      const char *rowname,
      dbl_rawlpdata * lp,
      int allowNewCols)
{
    int rval = 0;
    int rowind;
    char sense;
    double d;
    dbl_EGlpNumInitVar (d);

    if ((rowname != NULL) &&
	(ILLsymboltab_lookup (&lp->rowtab, rowname, &rowind) == 0)) {
	rval = dbl_ILLlp_error (state, "Repeated row name \"%s\".\n", rowname);
	ILL_CLEANUP_IF (rval);
    }
    rowind = lp->nrows;
    rval = rval || dbl_ILLraw_add_row (lp, rowname, 'N', dbl_zeroLpNum);

    rval = rval || dbl_ILLread_constraint_expr (state, lp, rowind, allowNewCols);
    rval = rval || dbl_ILLread_lp_state_sense (state);
    sense = state->sense_val;
    if (rval == 0) {
	rval = dbl_ILLread_lp_state_value (state, &d);
	if (rval) {
	    (void) dbl_ILLlp_error (state, "No right hand side value in constraint.\n");
	}
    }
    if (rval == 0) {
	lp->rowsense[rowind] = sense;
	dbl_EGlpNumCopy (lp->rhs[rowind], d);
	ILL_IFTRACE ("SENSE \"%s\": %c %f\n",
	    dbl_ILLraw_rowname (lp, rowind), sense, dbl_EGlpNumToLf (d));
    }
CLEANUP:
    dbl_EGlpNumClearVar (d);
    ILL_RESULT (rval, "dbl_ILLread_one_constraint");
}

static int dbl_read_constraints (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp,
      int allowNewCols)
{
    int rval = 0;
    char *rowname = NULL;

    if (dbl_ILLcheck_subject_to (state) != 0) {
	return dbl_ILLlp_error (state, "Constraint section expected.\n");
    }
    while (rval == 0) {
	rval = dbl_ILLread_constraint_name (state, &rowname);
	if (rval == 0) {
	    rval = dbl_ILLread_one_constraint (state, rowname, lp, allowNewCols);
	}
	if (rval == 0) {
	    if (dbl_ILLread_lp_state_next_constraint (state) != 0) {
		break;
	    }
	}
    }
    dbl_ILLread_lp_state_next_field (state);
    ILL_RESULT (rval, "dbl_read_constraints");
}

/* return -2 iff next is not a variable and not a keyword return -1 iff next
   is a keyword and !mustHave return 1  iff unknown column name or mustHave
   and keyword return 0  for success */
static int dbl_read_colname (dbl_ILLread_lp_state * state,
      ILLsymboltab * coltab,
      int mustHave)
{
    int rval = 0;
    int colind = ILL_SYM_NOINDEX;
    rval = dbl_ILLread_lp_state_next_var (state);
    if (mustHave && (rval != 0)) {
	return dbl_ILLlp_error (state, "Expecting a column name.\n");
    }
    if (rval != 0) {
	return (rval == -1) ? rval : -2;
    }
    if (ILLsymboltab_lookup (coltab, state->field, &colind)) {
	dbl_ILLread_lp_state_prev_field (state);
	return dbl_ILLlp_error (state, "\"%s\" is not a column name.\n", state->field);
    }
    state->column_index = colind;
    return 0;
}

static int dbl_read_integer (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp)
{
    int rval = 0;
    ILLsymboltab *coltab = &lp->coltab;

    ILL_FAILfalse (lp->intmarker, "Programming error");

    while ((rval = dbl_read_colname (state, coltab, 0)) == 0) {
	ILL_FAILtrue (state->column_index == ILL_SYM_NOINDEX, "Programming error");
	lp->intmarker[state->column_index] = 1;
    }
CLEANUP:
    if (rval == -1) {		/* last try for a colname gave us a keyword */
	rval = 0;
    } else {
	rval = dbl_ILLlp_error (state, "Expecting a column name.");
    }
    dbl_ILLread_lp_state_next_field (state);
    ILL_RESULT (rval, "dbl_read_integer");
}

static int dbl_read_bounds (dbl_ILLread_lp_state * state,
      dbl_rawlpdata * lp)
{
    int rval = 0;
    int colind, haveBound;
    char sense;
    const char *msg;
    ILLsymboltab *coltab;

    dbl_ILLraw_init_bounds (lp);
    coltab = &lp->coltab;

    while (1) {
	colind = -1;
	haveBound = 0;
	if (dbl_ILLread_lp_state_possible_bound_value (state)) {
	    /* this must be for a lower bound */
	    dbl_ILLtest_lp_state_bound_sense (state);
	    if (state->sense_val != 'L') {
		rval = dbl_ILLlp_error (state, "Expecting \"<=\".\n");
		break;
	    }
	    rval = dbl_read_colname (state, coltab, 1);
	    if (rval != 0) {
		break;
	    }
	    colind = state->column_index;
	    /* add lower bound value */
	    msg = dbl_ILLraw_set_lowerBound (lp, colind, state->bound_val);
	    dbl_ILLlp_warn (state, msg);
	    haveBound = 1;
	}
	if (colind == -1) {
	    rval = dbl_read_colname (state, coltab, 0);
	    colind = state->column_index;
	    if (rval != 0) {
		if (rval == -1) {
		    rval = 0;	/* found a keyword and that's OK */
		} else if (rval == -2) {
		    rval = dbl_ILLlp_error (state, "Expecting a column name.\n");
		}
		break;
	    }
	}
	ILL_FAILtrue (colind == -1, "must have a valid colname");
	dbl_ILLtest_lp_state_bound_sense (state);
	if (state->sense_val != ' ') {
	    sense = state->sense_val;
	    if ((sense != 'L') && (sense != 'E')) {
		rval = dbl_ILLlp_error (state, "Expecting \"<=\" or \"=\".\n");
		break;
	    }
	    if (dbl_ILLread_lp_state_possible_bound_value (state)) {
		if (sense == 'E') {
		    msg = dbl_ILLraw_set_fixedBound (lp, colind, state->bound_val);
		} else {
		    msg = dbl_ILLraw_set_upperBound (lp, colind, state->bound_val);
		}
		dbl_ILLlp_warn (state, msg);
		haveBound = 1;
	    } else {
		rval = dbl_ILLlp_error (state, "Expecting bound value.\n");
		break;
	    }
	} else {
	    if (dbl_ILLtest_lp_state_next_is (state, "FREE")) {
		msg = dbl_ILLraw_set_unbound (lp, colind);
		dbl_ILLlp_warn (state, msg);
		haveBound = 1;
	    } else {
		if (!haveBound) {
		    rval = dbl_ILLlp_error (state, "Not a bound expression.\n");
		    break;
		}
	    }
	}
	ILL_IFTRACE ("BOUNDS: %f <= %s <= %f\n", dbl_EGlpNumToLf (lp->lower[colind]),
	    dbl_ILLraw_colname (lp, colind), dbl_EGlpNumToLf (lp->upper[colind]));
    }
    dbl_ILLread_lp_state_next_field (state);
CLEANUP:
    ILL_RESULT (rval, "dbl_read_bounds");
}

static int dbl_add_var (dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state,
      double coef,
      int row,
      int allowNew)
{
    char *var = state->field;
    int rval = 0;
    int colind;

    if (ILLsymboltab_lookup (&lp->coltab, var, &colind) != 0) {
	if (!allowNew) {
	    rval = dbl_ILLlp_error (state, "Unknown col name \"%s\".\n", var);
	}
	ILL_CLEANUP_IF (rval);
	rval = dbl_ILLraw_add_col (lp, var, 0 /* not an integer var */ );
	colind = lp->ncols - 1;
	ILL_CLEANUP_IF (rval);
    }
    ILL_IFTRACE ("dbl_add_var: \"%s\" coef=%f row=%s\n",
	var, dbl_EGlpNumToLf (coef), dbl_ILLraw_rowname (lp, row));
    rval = dbl_ILLraw_add_col_coef (lp, colind, row, coef);
CLEANUP:
    ILL_RESULT (rval, "dbl_add_var");
}
