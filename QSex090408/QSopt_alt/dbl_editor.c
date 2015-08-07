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

/* RCS_INFO = "$RCSfile: dbl_editor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

#include "econfig.h"
#include "dbl_qsopt.h"
#include "dbl_lpdata.h"
#include "dbl_qstruct.h"
#include "dbl_qsopt.h"
#include "dbl_editor.h"
#include "dbl_readline.h"
#include "dbl_rawlp.h"
#include "stddefs.h"		/* for MAX */
#include "dbl_read_lp.h"
#include "dbl_lp.h"
#include "dbl_lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

#define dbl_ILL_BREAK_BODY_IF(rval) if (rval != 0) goto CLEANUP
#define dbl_ILL_BREAK_BODY goto CLEANUP

static int dbl_transpose (dbl_rawlpdata * lp);
static int dbl_pull_info_from_p (dbl_QSdata * p,
      dbl_rawlpdata * lp);
static void dbl_add_row (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state);
/* static int new_row(dbl_QSdata *p, dbl_rawlpdata *lp, dbl_ILLread_lp_state
   state); */
static void dbl_del_row (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state);

static void dbl_add_col (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state);
static void dbl_del_col (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state);

#define dbl_NONE -1
#define dbl_QS_EXIT 0
#define dbl_ROW 1
#define dbl_COL 2
#define dbl_PLP 3
#define dbl_PRTX 4
#define dbl_SOLVE 5
#define dbl_PMPS 6
#define dbl_HELP 7
#define dbl_DEL 8
#define dbl_NEW 9
#define dbl_ADD 10
#define dbl_PRIMAL 11
#define dbl_DUAL 12
#define dbl_NCOMMAND 13
static const char *dbl_commands[dbl_NCOMMAND + 1];
static char dbl_hasSubCmd[dbl_NCOMMAND + 1];

void dbl_ILLeditor_init (void)
{
    dbl_commands[dbl_QS_EXIT] = "dbl_QS_EXIT";
    dbl_commands[dbl_ROW] = "dbl_ROW";
    dbl_commands[dbl_COL] = "dbl_COL";
    dbl_commands[dbl_PLP] = "LP";
    dbl_commands[dbl_PMPS] = "MPS";
    dbl_commands[dbl_SOLVE] = "dbl_SOLVE";
    dbl_commands[dbl_PRTX] = "PRT";
    dbl_commands[dbl_HELP] = "dbl_HELP";
    dbl_commands[dbl_ADD] = "dbl_ADD";
    dbl_commands[dbl_DEL] = "dbl_DEL";
    dbl_commands[dbl_NEW] = "dbl_NEW";
    dbl_commands[dbl_PRIMAL] = "dbl_PRIMAL";
    dbl_commands[dbl_DUAL] = "dbl_DUAL";
    dbl_commands[dbl_NCOMMAND] = NULL;

    dbl_hasSubCmd[dbl_QS_EXIT] = 0;
    dbl_hasSubCmd[dbl_ROW] = 1;
    dbl_hasSubCmd[dbl_COL] = 1;
    dbl_hasSubCmd[dbl_PLP] = 0;
    dbl_hasSubCmd[dbl_PMPS] = 0;
    dbl_hasSubCmd[dbl_SOLVE] = 1;
    dbl_hasSubCmd[dbl_PRTX] = 0;
    dbl_hasSubCmd[dbl_HELP] = 0;
    dbl_hasSubCmd[dbl_ADD] = 1;
    dbl_hasSubCmd[dbl_DEL] = 1;
    dbl_hasSubCmd[dbl_NEW] = 1;
    dbl_hasSubCmd[dbl_PRIMAL] = 1;
    dbl_hasSubCmd[dbl_DUAL] = 1;
    dbl_hasSubCmd[dbl_NCOMMAND] = 0;
}

static void dbl_ILLeditor_help_cmd (int cmd,
      int subcmd);

static void dbl_ILLeditor_help (void)
{
    dbl_ILLeditor_help_cmd (dbl_ROW, dbl_ADD);
    /* dbl_ILLeditor_help_cmd(dbl_ROW, dbl_NEW);  */
    dbl_ILLeditor_help_cmd (dbl_ROW, dbl_DEL);
    dbl_ILLeditor_help_cmd (dbl_COL, dbl_ADD);
    dbl_ILLeditor_help_cmd (dbl_COL, dbl_DEL);
    dbl_ILLeditor_help_cmd (dbl_SOLVE, dbl_NONE);
    dbl_ILLeditor_help_cmd (dbl_PRTX, dbl_NONE);
    dbl_ILLeditor_help_cmd (dbl_PLP, dbl_NONE);
    dbl_ILLeditor_help_cmd (dbl_PMPS, dbl_NONE);
    dbl_ILLeditor_help_cmd (dbl_QS_EXIT, dbl_NONE);
    dbl_ILLeditor_help_cmd (dbl_HELP, dbl_NONE);
}

static void dbl_ILLeditor_help_cmd (int cmd,
      int subcmd)
{
    if (cmd == dbl_ROW && subcmd == dbl_ADD)
	fprintf (stdout, "%s dbl_ADD:\t%s.\n",
	    dbl_commands[dbl_ROW], "add a row; enter in LP format");
    if (cmd == dbl_COL && subcmd == dbl_ADD)
	fprintf (stdout, "%s dbl_ADD:\t%s.\n",
	    dbl_commands[dbl_COL], "add a col; enter in LP format");
    /* if (cmd == dbl_ROW && subcmd == dbl_NEW) fprintf(stdout, "%s
       dbl_NEW:\t%s.\n", dbl_commands[dbl_ROW], "new row; enter rowname:
       sense rhs");  */
    if (cmd == dbl_ROW && subcmd == dbl_DEL)
	fprintf (stdout, "%s dbl_DEL:\t%s.\n",
	    dbl_commands[dbl_ROW], "delete a row; give rowname");
    if (cmd == dbl_COL && subcmd == dbl_DEL)
	fprintf (stdout, "%s dbl_DEL:\t%s.\n",
	    dbl_commands[dbl_COL], "delete a col; give colname");
    if (cmd == dbl_SOLVE)
	fprintf (stdout, "%s:\t%s.\n", dbl_commands[dbl_SOLVE], "solve problem");
    if (cmd == dbl_PRTX)
	fprintf (stdout, "%s:\t%s.\n",
	    dbl_commands[dbl_PRTX], "print variable values for optimal solution");
    if (cmd == dbl_PLP)
	fprintf (stdout, "%s [file]:\t%s.\n",
	    dbl_commands[dbl_PLP], "print problem in LP format to file or stdout");
    if (cmd == dbl_PMPS)
	fprintf (stdout, "%s [file]:\t%s.\n",
	    dbl_commands[dbl_PMPS], "print problem in MPS format to file or stdout");
    if (cmd == dbl_QS_EXIT)
	fprintf (stdout, "%s:\t%s.\n", dbl_commands[dbl_QS_EXIT], "dbl_QS_EXIT");
    if (cmd == dbl_HELP)
	fprintf (stdout, "%s:\t%s.\n", dbl_commands[dbl_HELP], "print this help");
}

static void dbl_getCmd (dbl_ILLread_lp_state * state,
      int *cmd,
      int *subcmd)
{
    const char *cmd_str, *subcmd_str;
    int tmp;

    *cmd = dbl_ILLutil_index (dbl_commands, state->field);
    *subcmd = -1;
    if (dbl_hasSubCmd[*cmd] && (dbl_ILLread_lp_state_next_field_on_line (state) == 0)) {
	*subcmd = dbl_ILLutil_index (dbl_commands, state->field);
	if ((*subcmd == dbl_ROW) || (*subcmd == dbl_COL) || (*subcmd == dbl_SOLVE)) {
	    dbl_ILL_SWAP (*subcmd, *cmd, tmp);
	}
    }
    cmd_str = (*cmd >= 0) ? dbl_commands[*cmd] : "???";
    subcmd_str = (*subcmd >= 0) ? dbl_commands[*subcmd] : "???";
    ILL_IFTRACE ("cmd = %s, subcmd = %s\n", cmd_str, subcmd_str);
}

void dbl_ILLeditor (dbl_QSdata * p)
{
    dbl_rawlpdata raw, *lp = &raw;
    int cmd, subcmd, tval, rval = 0;
    dbl_ILLread_lp_state lpstate, *state = &lpstate;
    dbl_qsline_reader *reader;

    ILL_IFTRACE ("dbl_ILLeditor\n");

    reader = dbl_ILLline_reader_new ((dbl_qsread_line_fct) fgets, stdin);
    rval = dbl_ILLread_lp_state_init (state, reader, "STDIN", 1);
    rval = rval || dbl_pull_info_from_p (p, lp);
    dbl_ILL_BREAK_BODY_IF (rval);

    while (dbl_ILLread_lp_state_next_field (state) == 0) {
	dbl_getCmd (state, &cmd, &subcmd);
	switch (cmd) {
	case dbl_QS_EXIT:
	    dbl_ILL_BREAK_BODY;

	case dbl_ROW:
	    {
		switch (subcmd) {
		case dbl_ADD:
		    dbl_add_row (p, lp, state);
		    break;
		    /* case dbl_NEW: rval = new_row(p, lp, state); break; */
		case dbl_DEL:
		    dbl_del_row (p, lp, state);
		    break;
		default:
		    dbl_ILLeditor_help ();
		    break;
		}
		break;
	    }
	case dbl_COL:
	    {
		switch (subcmd) {
		case dbl_ADD:
		    dbl_add_col (p, lp, state);
		    break;
		case dbl_DEL:
		    dbl_del_col (p, lp, state);
		    break;
		default:
		    dbl_ILLeditor_help ();
		    break;
		}
		break;
	    }

	case dbl_SOLVE:
	    {
		if (subcmd == dbl_PRIMAL) {
		    (void) dbl_ILLeditor_solve (p, PRIMAL_SIMPLEX);
		} else if (subcmd == dbl_DUAL) {
		    (void) dbl_ILLeditor_solve (p, DUAL_SIMPLEX);
		} else {
		    dbl_ILLeditor_help ();
		}
		break;
	    }

	case dbl_PRTX:
	    {
		if ((rval = dbl_ILLlib_print_x (stdout, p->lp, 0, 0, 1))) {
		    fprintf (stdout, "The problem may not be feasible.\n");
		}
		break;
	    }

	case dbl_PLP:
	case dbl_PMPS:
	    {
		if (dbl_ILLread_lp_state_next_field_on_line (state) == 0) {
		    if (cmd == dbl_PMPS) {
			tval = dbl_QSwrite_prob (p, state->field, "MPS");
		    } else {
			tval = dbl_QSwrite_prob (p, state->field, "LP");
		    }
		    if (tval) {
			fprintf (stdout, "Could not write problem to \"%s\".\n",
			    state->field);
		    } else {
			fprintf (stdout, "Saved to \"%s\".\n", state->field);
		    }
		} else {
		    if (cmd == dbl_PMPS) {
			(void) dbl_QSwrite_prob_file (p, stdout, "MPS");
		    } else {
			(void) dbl_QSwrite_prob_file (p, stdout, "LP");
		    }
		}
		break;
	    }

	case dbl_NONE:
	    fprintf (stdout, "Unknown command: %s\n", state->field);
	default:
	    dbl_ILLeditor_help ();
	    break;
	}
	fflush (stdout);
	dbl_ILLread_lp_state_next_line (state);
    }
CLEANUP:
    dbl_ILLline_reader_free (reader);
    dbl_ILLfree_rawlpdata (lp);
}

int dbl_ILLeditor_solve (dbl_QSdata * p,
      int salgo)
{
    int rval = 0;
    int status = 0;
    double val;
    dbl_EGlpNumInitVar (val);

    if (salgo == PRIMAL_SIMPLEX) {
	rval = dbl_QSopt_primal (p, &status);
    } else {
	rval = dbl_QSopt_dual (p, &status);
    }
    dbl_ILL_BREAK_BODY_IF (rval);
    rval = dbl_QSget_objval (p, &val);
    if (p->simplex_display)
	if (rval == 0) {
	    fprintf (stdout, "LP Value: %.6f, status %d\n", dbl_EGlpNumToLf (val),
		status);
	    fflush (stdout);
	}
CLEANUP:
    dbl_EGlpNumClearVar (val);
    ILL_RESULT (rval, "dbl_ILLeditor_solve");
}


static int dbl_pull_info_from_p (dbl_QSdata * p,
      dbl_rawlpdata * lp)
{
    int i, rval = 0;
    dbl_ILLlpdata *qslp = p->lp->O;
    int nrows, ncols;

    dbl_ILLinit_rawlpdata (lp, NULL);
    rval = ILLsymboltab_create (&lp->rowtab, 100) ||
	ILLsymboltab_create (&lp->coltab, 100);
    dbl_ILL_BREAK_BODY_IF (rval);

    nrows = qslp->nrows;
    ncols = qslp->nstruct;
    /* add rows to lp */
    dbl_ILLraw_add_row (lp, qslp->objname, 'N', dbl_zeroLpNum);
    for (i = 0; i < nrows; i++) {
	ILL_FAILfalse (qslp->rownames[i] != NULL, "should have no NULL names");
	dbl_ILLraw_add_row (lp, qslp->rownames[i], qslp->sense[i], qslp->rhs[i]);
    }

    /* add cols to coltab and lp */
    for (i = 0; i < ncols; i++) {
	ILL_FAILfalse (qslp->colnames[i] != NULL, "should have no NULL names");
	dbl_ILLraw_add_col (lp, qslp->colnames[i],
	    (qslp->intmarker) ? qslp->intmarker[i] : 0);
    }
CLEANUP:
    ILL_RETURN (rval, "dbl_pull_info_from_p");
}

static int dbl_transpose (dbl_rawlpdata * lp)
{
    int rval = 0;
    int tmp;
    ILLsymboltab tmptab;

    tmp = MAX (lp->nrows, lp->ncols);
    if (tmp >= lp->sensesize) {
	lp->sensesize *= 1.3;
	lp->sensesize += 1000;
	if (lp->sensesize < tmp + 1)
	    lp->sensesize = tmp + 1;
	lp->rowsense = EGrealloc (lp->rowsense, sizeof (char) * lp->sensesize);
	/*
			rval = ILLutil_reallocrus_scale ((void **) &lp->rowsense, &lp->sensesize, tmp + 1, 1.3, sizeof (char));
			ILL_CLEANUP_IF (rval);
	*/
    }
    if (tmp >= lp->rhssize) {
	lp->rhssize *= 1.3;
	lp->rhssize += 1000;
	if (lp->rhssize < tmp + 1)
	    lp->rhssize = tmp + 1;
	dbl_EGlpNumReallocArray (&(lp->rhs), lp->rhssize);
	/*
			lp->rhs = EGrealloc(lp->rhs, sizeof(double)*lp->rhssize);
			rval = ILLutil_reallocrus_scale ((void **) &lp->rhs, &lp->sensesize, tmp + 1, 1.3, sizeof (double));
			ILL_CLEANUP_IF (rval);
	*/
    }
    dbl_ILL_SWAP (lp->nrows, lp->ncols, tmp);
    dbl_ILL_SWAP (lp->rowtab, lp->coltab, tmptab);
    ILL_RETURN (rval, "dbl_transpose");
}

static char *dbl_get_row_col_name (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state,
      int doRow)
{
    int rval = 0;
    int ind;
    char *rname, *thename = NULL;
    char buf[ILL_namebufsize];
    ILLsymboltab *tab = (doRow) ? &lp->rowtab : &lp->coltab;
    int id = (doRow) ? lp->nrows : lp->ncols;
    id--;			/* in dbl_rawlpdata obj counts as a row */

    rval = dbl_ILLread_constraint_name (state, &rname);
    dbl_ILL_BREAK_BODY_IF (rval);

    if (rname == NULL) {
	dbl_ILLlib_findName (p->qslp, doRow /* forRow */ , rname, id, buf);
	dbl_ILL_UTIL_STR (thename, buf);
    } else {
	dbl_ILL_UTIL_STR (thename, rname);
    }
    if (ILLsymboltab_lookup (tab, thename, &ind) == 0) {
	rval = dbl_ILLlp_error (state, "\"%s\" already exists.", thename);
    }
CLEANUP:
    if (rval != 0) {
	ILL_IFFREE (thename, char);
    }
    return thename;
}

static int dbl_fill_matrix (dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state,
      dbl_ILLmatrix * m,
      double *obj,
      int n)
{
    int i, cnt, rval = 0;
    dbl_colptr *cp;
    double val;
    int newCol = (obj != NULL);
    dbl_EGlpNumInitVar (val);

    /* rely on fact that objective has rowindex 0 */

    m->matrows = lp->nrows;
    m->matcols = 1;
    m->matval = dbl_EGlpNumAllocArray (lp->ncols);
    ILL_SAFE_MALLOC (m->matind, lp->ncols, int);
    ILL_SAFE_MALLOC (m->matbeg, 1, int);
    ILL_SAFE_MALLOC (m->matcnt, 1, int);
    m->matsize = lp->ncols;
    m->matbeg[0] = 0;
    m->matcnt[0] = 0;
    for (i = 0; i < lp->ncols; i++) {
	cnt = 0;
	dbl_EGlpNumZero (val);
	for (cp = lp->cols[i]; cp != NULL; cp = cp->next) {
	    ILL_FAILfalse (cp->this == n, "n should be the only row around");
	    if (dbl_EGlpNumIsNeqq (cp->coef, dbl_zeroLpNum)) {
		dbl_EGlpNumAddTo (val, cp->coef);
		cnt++;
	    }
	}
	if (cnt > 1) {
	    dbl_ILLlp_warn (state, "Multiple coefficients for \"%s\".",
		dbl_ILLraw_colname (lp, i));
	}
	if (dbl_EGlpNumIsNeqqZero (val)) {
	    if ((i - newCol) >= 0) {
		dbl_EGlpNumCopy (m->matval[m->matcnt[0]], val);
		m->matind[m->matcnt[0]] = i - newCol;
		m->matcnt[0]++;
	    } else {
		dbl_EGlpNumCopy (obj[0], val);
	    }
	}
    }
    if (m->matcnt[0] == 0) {
	rval = dbl_ILLlp_error (state, "There are no non zero coefficients.");
    }
CLEANUP:
    dbl_EGlpNumClearVar (val);
    ILL_RESULT (rval, "dbl_fill_matrix");
}

static void dbl_add_row (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state)
{
    int rval = 0;
    int n;
    char *name;
    dbl_ILLmatrix m;
    char sense[1];

    dbl_ILLmatrix_init (&m);
    n = lp->nrows;
    name = dbl_get_row_col_name (p, lp, state, 1 /* doRow */ );

    if (name == NULL) {
	rval = 1;
    } else {
	rval = dbl_ILLread_one_constraint (state, name, lp, 0);

	/* adds row name to lp->rowtab the checks constraint expression  */
	if (rval != 0) {
	    /* failed because of error in expression => must remove name from
	       symbol table */
	    fprintf (stdout, "Incorrect expression.\n");
	} else {
	    ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");
	    ILL_IFTRACE ("ADDING row %s.\n", name);

	    sense[0] = lp->rowsense[n];

	    rval = dbl_fill_matrix (lp, state, &m, NULL, n);
	    dbl_ILL_BREAK_BODY_IF (rval);

	    dbl_QSadd_rows (p, 1, m.matcnt, m.matbeg, m.matind, m.matval,
		&(lp->rhs[n]), sense, (const char **) &name);
	}
    }
CLEANUP:
    dbl_ILLmatrix_free (&m);
    if (name != NULL) {
	if (rval != 0)
	    ILLsymboltab_delete (&lp->rowtab, name);
	ILL_IFFREE (name, char);
    }
    if (rval != 0) {
	lp->nrows = n;
    }
    dbl_ILLraw_clear_matrix (lp);
}

static void dbl_add_col (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state)
{
    int rval = 0;
    int n;
    char *name[1];
    int transposed = 1;
    dbl_ILLmatrix matrix, *m = &matrix;
    double obj[1], lower[1], upper[2];
    dbl_EGlpNumInitVar (*obj);
    dbl_EGlpNumInitVar (*lower);
    dbl_EGlpNumInitVar (upper[0]);
    dbl_EGlpNumInitVar (upper[1]);

    n = lp->ncols;
    dbl_ILLmatrix_init (m);
    name[0] = dbl_get_row_col_name (p, lp, state, 0 /* doRow */ );
    rval = (name[0] == NULL);
    dbl_ILL_BREAK_BODY_IF (rval);

    transposed = !dbl_transpose (lp);
    rval = dbl_ILLread_one_constraint (state, name[0], lp, 0);

    /* adds row name to lp->rowtab the checks constraint expression  */
    if (rval != 0) {
	/* failed because of error in expression => must remove name from
	   symbol table */
	fprintf (stdout, "Incorrect expression.\n");
    } else {
	ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");

	rval = dbl_fill_matrix (lp, state, m, obj, n);
	dbl_ILL_BREAK_BODY_IF (rval);

	fprintf (stdout, "lower ");
	rval = dbl_ILLread_lp_state_next_line (state) ||
	    dbl_ILLread_lp_state_value (state, &(lower[0]));
	dbl_ILL_BREAK_BODY_IF (rval);

	fprintf (stdout, "upper ");
	rval = dbl_ILLread_lp_state_next_line (state) ||
	    dbl_ILLread_lp_state_value (state, &(upper[0]));
	dbl_ILL_BREAK_BODY_IF (rval);

	ILL_IFTRACE ("ADDING col %s.\n", name[0]);

	dbl_QSadd_cols (p, 1, m->matcnt, m->matbeg, m->matind, m->matval,
	    obj, lower, upper, (const char **) name);

    }
CLEANUP:
    dbl_ILLmatrix_free (m);
    if (name[0] != NULL) {
	if (rval != 0)
	    ILLsymboltab_delete (&lp->rowtab, name[0]);
	ILL_IFFREE (name[0], char);
    }
    if (rval != 0) {
	lp->nrows = n;
    }
    dbl_ILLraw_clear_matrix (lp);
    if (transposed)
	dbl_transpose (lp);
    ILL_IFFREE (name[0], char);
    dbl_EGlpNumClearVar (*obj);
    dbl_EGlpNumClearVar (*lower);
    dbl_EGlpNumClearVar (upper[0]);
    dbl_EGlpNumClearVar (upper[1]);
}

#if 0
#ifndef JAVA_PORT
static void new_row (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state)
{
    int rval = 0;
    char *rowname = NULL, *rname = NULL;
    char sense;
    double d;
    int ind, hit;

    rval = dbl_ILLread_constraint_name (state, &rname);
    if (rname == NULL) {
	rval = 1;
	dbl_ILLeditor_help_cmd (dbl_ROW, dbl_NEW);
    }
    dbl_ILL_BREAK_BODY_IF (rval);

    ILLsymboltab_lookup (&lp->rowtab, rname, &ind);
    if (ind != ILL_SYM_NOINDEX) {
	rval = dbl_ILLlp_error (state, "\"%s\" is already defined.\n", rname);
	dbl_ILL_BREAK_BODY_IF (rval);
    }
    dbl_ILL_UTIL_STR (rowname, rname);

    rval = dbl_ILLread_lp_state_sense (state);
    sense = state->sense_val;
    dbl_ILL_BREAK_BODY_IF (rval);

    rval = dbl_ILLread_lp_state_value (state, &d);
    dbl_ILL_BREAK_BODY_IF (rval);

    rval = dbl_QSnew_row (p, d, sense, rowname);
    if (rval != 0) {
	fprintf (stderr, "could not add row\n");
    } else {
	ILLsymboltab_register (&lp->rowtab, rname, &ind, &hit);
    }
CLEANUP:
    ILL_IFFREE (rowname, char);
}
#endif
#endif

static int dbl_del_row_or_col (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state,
      int isRow)
{
    int i[1], rval = 0;
    char **names = (isRow) ? p->qslp->rownames : p->qslp->colnames;
    int nnames = (isRow) ? p->qslp->nrows : p->qslp->nstruct;
    ILLsymboltab *tab = (isRow) ? &lp->rowtab : &lp->coltab;

    rval = dbl_ILLread_lp_state_next_field_on_line (state);
    dbl_ILL_BREAK_BODY_IF (rval);

    i[0] = dbl_ILLutil_array_index (names, nnames, state->field);
    if (i[0] >= 0) {
	rval = (isRow) ? dbl_QSdelete_rows (p, 1, i) : dbl_QSdelete_cols (p, 1, i);
	if (rval == 0) {
	    ILLsymboltab_delete (tab, state->field);
	}
    } else {
	rval = dbl_ILLlp_error (state, "\"%s\" is not defined.\n", state->field);
    }

CLEANUP:
    ILL_RESULT (rval, "dbl_del_row_or_col");
}

static void dbl_del_row (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state)
{
    int rval = dbl_del_row_or_col (p, lp, state, 1);
    if (rval == 0) {
	lp->nrows--;
    }
}

static void dbl_del_col (dbl_QSdata * p,
      dbl_rawlpdata * lp,
      dbl_ILLread_lp_state * state)
{
    int rval = dbl_del_row_or_col (p, lp, state, 0);
    if (rval == 0) {
	lp->ncols--;
    }
}
