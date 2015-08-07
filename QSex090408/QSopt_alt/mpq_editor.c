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

/* RCS_INFO = "$RCSfile: mpq_editor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

#include "econfig.h"
#include "mpq_qsopt.h"
#include "mpq_lpdata.h"
#include "mpq_qstruct.h"
#include "mpq_qsopt.h"
#include "mpq_editor.h"
#include "mpq_readline.h"
#include "mpq_rawlp.h"
#include "stddefs.h"		/* for MAX */
#include "mpq_read_lp.h"
#include "mpq_lp.h"
#include "mpq_lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

#define mpq_ILL_BREAK_BODY_IF(rval) if (rval != 0) goto CLEANUP
#define mpq_ILL_BREAK_BODY goto CLEANUP

static int mpq_transpose (mpq_rawlpdata * lp);
static int mpq_pull_info_from_p (mpq_QSdata * p,
      mpq_rawlpdata * lp);
static void mpq_add_row (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state);
/* static int new_row(mpq_QSdata *p, mpq_rawlpdata *lp, mpq_ILLread_lp_state
   state); */
static void mpq_del_row (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state);

static void mpq_add_col (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state);
static void mpq_del_col (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state);

#define mpq_NONE -1
#define mpq_QS_EXIT 0
#define mpq_ROW 1
#define mpq_COL 2
#define mpq_PLP 3
#define mpq_PRTX 4
#define mpq_SOLVE 5
#define mpq_PMPS 6
#define mpq_HELP 7
#define mpq_DEL 8
#define mpq_NEW 9
#define mpq_ADD 10
#define mpq_PRIMAL 11
#define mpq_DUAL 12
#define mpq_NCOMMAND 13
static const char *mpq_commands[mpq_NCOMMAND + 1];
static char mpq_hasSubCmd[mpq_NCOMMAND + 1];

void mpq_ILLeditor_init (void)
{
    mpq_commands[mpq_QS_EXIT] = "mpq_QS_EXIT";
    mpq_commands[mpq_ROW] = "mpq_ROW";
    mpq_commands[mpq_COL] = "mpq_COL";
    mpq_commands[mpq_PLP] = "LP";
    mpq_commands[mpq_PMPS] = "MPS";
    mpq_commands[mpq_SOLVE] = "mpq_SOLVE";
    mpq_commands[mpq_PRTX] = "PRT";
    mpq_commands[mpq_HELP] = "mpq_HELP";
    mpq_commands[mpq_ADD] = "mpq_ADD";
    mpq_commands[mpq_DEL] = "mpq_DEL";
    mpq_commands[mpq_NEW] = "mpq_NEW";
    mpq_commands[mpq_PRIMAL] = "mpq_PRIMAL";
    mpq_commands[mpq_DUAL] = "mpq_DUAL";
    mpq_commands[mpq_NCOMMAND] = NULL;

    mpq_hasSubCmd[mpq_QS_EXIT] = 0;
    mpq_hasSubCmd[mpq_ROW] = 1;
    mpq_hasSubCmd[mpq_COL] = 1;
    mpq_hasSubCmd[mpq_PLP] = 0;
    mpq_hasSubCmd[mpq_PMPS] = 0;
    mpq_hasSubCmd[mpq_SOLVE] = 1;
    mpq_hasSubCmd[mpq_PRTX] = 0;
    mpq_hasSubCmd[mpq_HELP] = 0;
    mpq_hasSubCmd[mpq_ADD] = 1;
    mpq_hasSubCmd[mpq_DEL] = 1;
    mpq_hasSubCmd[mpq_NEW] = 1;
    mpq_hasSubCmd[mpq_PRIMAL] = 1;
    mpq_hasSubCmd[mpq_DUAL] = 1;
    mpq_hasSubCmd[mpq_NCOMMAND] = 0;
}

static void mpq_ILLeditor_help_cmd (int cmd,
      int subcmd);

static void mpq_ILLeditor_help (void)
{
    mpq_ILLeditor_help_cmd (mpq_ROW, mpq_ADD);
    /* mpq_ILLeditor_help_cmd(mpq_ROW, mpq_NEW);  */
    mpq_ILLeditor_help_cmd (mpq_ROW, mpq_DEL);
    mpq_ILLeditor_help_cmd (mpq_COL, mpq_ADD);
    mpq_ILLeditor_help_cmd (mpq_COL, mpq_DEL);
    mpq_ILLeditor_help_cmd (mpq_SOLVE, mpq_NONE);
    mpq_ILLeditor_help_cmd (mpq_PRTX, mpq_NONE);
    mpq_ILLeditor_help_cmd (mpq_PLP, mpq_NONE);
    mpq_ILLeditor_help_cmd (mpq_PMPS, mpq_NONE);
    mpq_ILLeditor_help_cmd (mpq_QS_EXIT, mpq_NONE);
    mpq_ILLeditor_help_cmd (mpq_HELP, mpq_NONE);
}

static void mpq_ILLeditor_help_cmd (int cmd,
      int subcmd)
{
    if (cmd == mpq_ROW && subcmd == mpq_ADD)
	fprintf (stdout, "%s mpq_ADD:\t%s.\n",
	    mpq_commands[mpq_ROW], "add a row; enter in LP format");
    if (cmd == mpq_COL && subcmd == mpq_ADD)
	fprintf (stdout, "%s mpq_ADD:\t%s.\n",
	    mpq_commands[mpq_COL], "add a col; enter in LP format");
    /* if (cmd == mpq_ROW && subcmd == mpq_NEW) fprintf(stdout, "%s
       mpq_NEW:\t%s.\n", mpq_commands[mpq_ROW], "new row; enter rowname:
       sense rhs");  */
    if (cmd == mpq_ROW && subcmd == mpq_DEL)
	fprintf (stdout, "%s mpq_DEL:\t%s.\n",
	    mpq_commands[mpq_ROW], "delete a row; give rowname");
    if (cmd == mpq_COL && subcmd == mpq_DEL)
	fprintf (stdout, "%s mpq_DEL:\t%s.\n",
	    mpq_commands[mpq_COL], "delete a col; give colname");
    if (cmd == mpq_SOLVE)
	fprintf (stdout, "%s:\t%s.\n", mpq_commands[mpq_SOLVE], "solve problem");
    if (cmd == mpq_PRTX)
	fprintf (stdout, "%s:\t%s.\n",
	    mpq_commands[mpq_PRTX], "print variable values for optimal solution");
    if (cmd == mpq_PLP)
	fprintf (stdout, "%s [file]:\t%s.\n",
	    mpq_commands[mpq_PLP], "print problem in LP format to file or stdout");
    if (cmd == mpq_PMPS)
	fprintf (stdout, "%s [file]:\t%s.\n",
	    mpq_commands[mpq_PMPS], "print problem in MPS format to file or stdout");
    if (cmd == mpq_QS_EXIT)
	fprintf (stdout, "%s:\t%s.\n", mpq_commands[mpq_QS_EXIT], "mpq_QS_EXIT");
    if (cmd == mpq_HELP)
	fprintf (stdout, "%s:\t%s.\n", mpq_commands[mpq_HELP], "print this help");
}

static void mpq_getCmd (mpq_ILLread_lp_state * state,
      int *cmd,
      int *subcmd)
{
    const char *cmd_str, *subcmd_str;
    int tmp;

    *cmd = mpq_ILLutil_index (mpq_commands, state->field);
    *subcmd = -1;
    if (mpq_hasSubCmd[*cmd] && (mpq_ILLread_lp_state_next_field_on_line (state) == 0)) {
	*subcmd = mpq_ILLutil_index (mpq_commands, state->field);
	if ((*subcmd == mpq_ROW) || (*subcmd == mpq_COL) || (*subcmd == mpq_SOLVE)) {
	    mpq_ILL_SWAP (*subcmd, *cmd, tmp);
	}
    }
    cmd_str = (*cmd >= 0) ? mpq_commands[*cmd] : "???";
    subcmd_str = (*subcmd >= 0) ? mpq_commands[*subcmd] : "???";
    ILL_IFTRACE ("cmd = %s, subcmd = %s\n", cmd_str, subcmd_str);
}

void mpq_ILLeditor (mpq_QSdata * p)
{
    mpq_rawlpdata raw, *lp = &raw;
    int cmd, subcmd, tval, rval = 0;
    mpq_ILLread_lp_state lpstate, *state = &lpstate;
    mpq_qsline_reader *reader;

    ILL_IFTRACE ("mpq_ILLeditor\n");

    reader = mpq_ILLline_reader_new ((mpq_qsread_line_fct) fgets, stdin);
    rval = mpq_ILLread_lp_state_init (state, reader, "STDIN", 1);
    rval = rval || mpq_pull_info_from_p (p, lp);
    mpq_ILL_BREAK_BODY_IF (rval);

    while (mpq_ILLread_lp_state_next_field (state) == 0) {
	mpq_getCmd (state, &cmd, &subcmd);
	switch (cmd) {
	case mpq_QS_EXIT:
	    mpq_ILL_BREAK_BODY;

	case mpq_ROW:
	    {
		switch (subcmd) {
		case mpq_ADD:
		    mpq_add_row (p, lp, state);
		    break;
		    /* case mpq_NEW: rval = new_row(p, lp, state); break; */
		case mpq_DEL:
		    mpq_del_row (p, lp, state);
		    break;
		default:
		    mpq_ILLeditor_help ();
		    break;
		}
		break;
	    }
	case mpq_COL:
	    {
		switch (subcmd) {
		case mpq_ADD:
		    mpq_add_col (p, lp, state);
		    break;
		case mpq_DEL:
		    mpq_del_col (p, lp, state);
		    break;
		default:
		    mpq_ILLeditor_help ();
		    break;
		}
		break;
	    }

	case mpq_SOLVE:
	    {
		if (subcmd == mpq_PRIMAL) {
		    (void) mpq_ILLeditor_solve (p, PRIMAL_SIMPLEX);
		} else if (subcmd == mpq_DUAL) {
		    (void) mpq_ILLeditor_solve (p, DUAL_SIMPLEX);
		} else {
		    mpq_ILLeditor_help ();
		}
		break;
	    }

	case mpq_PRTX:
	    {
		if ((rval = mpq_ILLlib_print_x (stdout, p->lp, 0, 0, 1))) {
		    fprintf (stdout, "The problem may not be feasible.\n");
		}
		break;
	    }

	case mpq_PLP:
	case mpq_PMPS:
	    {
		if (mpq_ILLread_lp_state_next_field_on_line (state) == 0) {
		    if (cmd == mpq_PMPS) {
			tval = mpq_QSwrite_prob (p, state->field, "MPS");
		    } else {
			tval = mpq_QSwrite_prob (p, state->field, "LP");
		    }
		    if (tval) {
			fprintf (stdout, "Could not write problem to \"%s\".\n",
			    state->field);
		    } else {
			fprintf (stdout, "Saved to \"%s\".\n", state->field);
		    }
		} else {
		    if (cmd == mpq_PMPS) {
			(void) mpq_QSwrite_prob_file (p, stdout, "MPS");
		    } else {
			(void) mpq_QSwrite_prob_file (p, stdout, "LP");
		    }
		}
		break;
	    }

	case mpq_NONE:
	    fprintf (stdout, "Unknown command: %s\n", state->field);
	default:
	    mpq_ILLeditor_help ();
	    break;
	}
	fflush (stdout);
	mpq_ILLread_lp_state_next_line (state);
    }
CLEANUP:
    mpq_ILLline_reader_free (reader);
    mpq_ILLfree_rawlpdata (lp);
}

int mpq_ILLeditor_solve (mpq_QSdata * p,
      int salgo)
{
    int rval = 0;
    int status = 0;
    mpq_t val;
    mpq_EGlpNumInitVar (val);

    if (salgo == PRIMAL_SIMPLEX) {
	rval = mpq_QSopt_primal (p, &status);
    } else {
	rval = mpq_QSopt_dual (p, &status);
    }
    mpq_ILL_BREAK_BODY_IF (rval);
    rval = mpq_QSget_objval (p, &val);
    if (p->simplex_display)
	if (rval == 0) {
	    fprintf (stdout, "LP Value: %.6f, status %d\n", mpq_EGlpNumToLf (val),
		status);
	    fflush (stdout);
	}
CLEANUP:
    mpq_EGlpNumClearVar (val);
    ILL_RESULT (rval, "mpq_ILLeditor_solve");
}


static int mpq_pull_info_from_p (mpq_QSdata * p,
      mpq_rawlpdata * lp)
{
    int i, rval = 0;
    mpq_ILLlpdata *qslp = p->lp->O;
    int nrows, ncols;

    mpq_ILLinit_rawlpdata (lp, NULL);
    rval = ILLsymboltab_create (&lp->rowtab, 100) ||
	ILLsymboltab_create (&lp->coltab, 100);
    mpq_ILL_BREAK_BODY_IF (rval);

    nrows = qslp->nrows;
    ncols = qslp->nstruct;
    /* add rows to lp */
    mpq_ILLraw_add_row (lp, qslp->objname, 'N', mpq_zeroLpNum);
    for (i = 0; i < nrows; i++) {
	ILL_FAILfalse (qslp->rownames[i] != NULL, "should have no NULL names");
	mpq_ILLraw_add_row (lp, qslp->rownames[i], qslp->sense[i], qslp->rhs[i]);
    }

    /* add cols to coltab and lp */
    for (i = 0; i < ncols; i++) {
	ILL_FAILfalse (qslp->colnames[i] != NULL, "should have no NULL names");
	mpq_ILLraw_add_col (lp, qslp->colnames[i],
	    (qslp->intmarker) ? qslp->intmarker[i] : 0);
    }
CLEANUP:
    ILL_RETURN (rval, "mpq_pull_info_from_p");
}

static int mpq_transpose (mpq_rawlpdata * lp)
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
	mpq_EGlpNumReallocArray (&(lp->rhs), lp->rhssize);
	/*
			lp->rhs = EGrealloc(lp->rhs, sizeof(double)*lp->rhssize);
			rval = ILLutil_reallocrus_scale ((void **) &lp->rhs, &lp->sensesize, tmp + 1, 1.3, sizeof (double));
			ILL_CLEANUP_IF (rval);
	*/
    }
    mpq_ILL_SWAP (lp->nrows, lp->ncols, tmp);
    mpq_ILL_SWAP (lp->rowtab, lp->coltab, tmptab);
    ILL_RETURN (rval, "mpq_transpose");
}

static char *mpq_get_row_col_name (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state,
      int doRow)
{
    int rval = 0;
    int ind;
    char *rname, *thename = NULL;
    char buf[ILL_namebufsize];
    ILLsymboltab *tab = (doRow) ? &lp->rowtab : &lp->coltab;
    int id = (doRow) ? lp->nrows : lp->ncols;
    id--;			/* in mpq_rawlpdata obj counts as a row */

    rval = mpq_ILLread_constraint_name (state, &rname);
    mpq_ILL_BREAK_BODY_IF (rval);

    if (rname == NULL) {
	mpq_ILLlib_findName (p->qslp, doRow /* forRow */ , rname, id, buf);
	mpq_ILL_UTIL_STR (thename, buf);
    } else {
	mpq_ILL_UTIL_STR (thename, rname);
    }
    if (ILLsymboltab_lookup (tab, thename, &ind) == 0) {
	rval = mpq_ILLlp_error (state, "\"%s\" already exists.", thename);
    }
CLEANUP:
    if (rval != 0) {
	ILL_IFFREE (thename, char);
    }
    return thename;
}

static int mpq_fill_matrix (mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state,
      mpq_ILLmatrix * m,
      mpq_t * obj,
      int n)
{
    int i, cnt, rval = 0;
    mpq_colptr *cp;
    mpq_t val;
    int newCol = (obj != NULL);
    mpq_EGlpNumInitVar (val);

    /* rely on fact that objective has rowindex 0 */

    m->matrows = lp->nrows;
    m->matcols = 1;
    m->matval = mpq_EGlpNumAllocArray (lp->ncols);
    ILL_SAFE_MALLOC (m->matind, lp->ncols, int);
    ILL_SAFE_MALLOC (m->matbeg, 1, int);
    ILL_SAFE_MALLOC (m->matcnt, 1, int);
    m->matsize = lp->ncols;
    m->matbeg[0] = 0;
    m->matcnt[0] = 0;
    for (i = 0; i < lp->ncols; i++) {
	cnt = 0;
	mpq_EGlpNumZero (val);
	for (cp = lp->cols[i]; cp != NULL; cp = cp->next) {
	    ILL_FAILfalse (cp->this == n, "n should be the only row around");
	    if (mpq_EGlpNumIsNeqq (cp->coef, mpq_zeroLpNum)) {
		mpq_EGlpNumAddTo (val, cp->coef);
		cnt++;
	    }
	}
	if (cnt > 1) {
	    mpq_ILLlp_warn (state, "Multiple coefficients for \"%s\".",
		mpq_ILLraw_colname (lp, i));
	}
	if (mpq_EGlpNumIsNeqqZero (val)) {
	    if ((i - newCol) >= 0) {
		mpq_EGlpNumCopy (m->matval[m->matcnt[0]], val);
		m->matind[m->matcnt[0]] = i - newCol;
		m->matcnt[0]++;
	    } else {
		mpq_EGlpNumCopy (obj[0], val);
	    }
	}
    }
    if (m->matcnt[0] == 0) {
	rval = mpq_ILLlp_error (state, "There are no non zero coefficients.");
    }
CLEANUP:
    mpq_EGlpNumClearVar (val);
    ILL_RESULT (rval, "mpq_fill_matrix");
}

static void mpq_add_row (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state)
{
    int rval = 0;
    int n;
    char *name;
    mpq_ILLmatrix m;
    char sense[1];

    mpq_ILLmatrix_init (&m);
    n = lp->nrows;
    name = mpq_get_row_col_name (p, lp, state, 1 /* doRow */ );

    if (name == NULL) {
	rval = 1;
    } else {
	rval = mpq_ILLread_one_constraint (state, name, lp, 0);

	/* adds row name to lp->rowtab the checks constraint expression  */
	if (rval != 0) {
	    /* failed because of error in expression => must remove name from
	       symbol table */
	    fprintf (stdout, "Incorrect expression.\n");
	} else {
	    ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");
	    ILL_IFTRACE ("ADDING row %s.\n", name);

	    sense[0] = lp->rowsense[n];

	    rval = mpq_fill_matrix (lp, state, &m, NULL, n);
	    mpq_ILL_BREAK_BODY_IF (rval);

	    mpq_QSadd_rows (p, 1, m.matcnt, m.matbeg, m.matind, m.matval,
		&(lp->rhs[n]), sense, (const char **) &name);
	}
    }
CLEANUP:
    mpq_ILLmatrix_free (&m);
    if (name != NULL) {
	if (rval != 0)
	    ILLsymboltab_delete (&lp->rowtab, name);
	ILL_IFFREE (name, char);
    }
    if (rval != 0) {
	lp->nrows = n;
    }
    mpq_ILLraw_clear_matrix (lp);
}

static void mpq_add_col (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state)
{
    int rval = 0;
    int n;
    char *name[1];
    int transposed = 1;
    mpq_ILLmatrix matrix, *m = &matrix;
    mpq_t obj[1], lower[1], upper[2];
    mpq_EGlpNumInitVar (*obj);
    mpq_EGlpNumInitVar (*lower);
    mpq_EGlpNumInitVar (upper[0]);
    mpq_EGlpNumInitVar (upper[1]);

    n = lp->ncols;
    mpq_ILLmatrix_init (m);
    name[0] = mpq_get_row_col_name (p, lp, state, 0 /* doRow */ );
    rval = (name[0] == NULL);
    mpq_ILL_BREAK_BODY_IF (rval);

    transposed = !mpq_transpose (lp);
    rval = mpq_ILLread_one_constraint (state, name[0], lp, 0);

    /* adds row name to lp->rowtab the checks constraint expression  */
    if (rval != 0) {
	/* failed because of error in expression => must remove name from
	   symbol table */
	fprintf (stdout, "Incorrect expression.\n");
    } else {
	ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");

	rval = mpq_fill_matrix (lp, state, m, obj, n);
	mpq_ILL_BREAK_BODY_IF (rval);

	fprintf (stdout, "lower ");
	rval = mpq_ILLread_lp_state_next_line (state) ||
	    mpq_ILLread_lp_state_value (state, &(lower[0]));
	mpq_ILL_BREAK_BODY_IF (rval);

	fprintf (stdout, "upper ");
	rval = mpq_ILLread_lp_state_next_line (state) ||
	    mpq_ILLread_lp_state_value (state, &(upper[0]));
	mpq_ILL_BREAK_BODY_IF (rval);

	ILL_IFTRACE ("ADDING col %s.\n", name[0]);

	mpq_QSadd_cols (p, 1, m->matcnt, m->matbeg, m->matind, m->matval,
	    obj, lower, upper, (const char **) name);

    }
CLEANUP:
    mpq_ILLmatrix_free (m);
    if (name[0] != NULL) {
	if (rval != 0)
	    ILLsymboltab_delete (&lp->rowtab, name[0]);
	ILL_IFFREE (name[0], char);
    }
    if (rval != 0) {
	lp->nrows = n;
    }
    mpq_ILLraw_clear_matrix (lp);
    if (transposed)
	mpq_transpose (lp);
    ILL_IFFREE (name[0], char);
    mpq_EGlpNumClearVar (*obj);
    mpq_EGlpNumClearVar (*lower);
    mpq_EGlpNumClearVar (upper[0]);
    mpq_EGlpNumClearVar (upper[1]);
}

#if 0
#ifndef JAVA_PORT
static void new_row (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state)
{
    int rval = 0;
    char *rowname = NULL, *rname = NULL;
    char sense;
    double d;
    int ind, hit;

    rval = mpq_ILLread_constraint_name (state, &rname);
    if (rname == NULL) {
	rval = 1;
	mpq_ILLeditor_help_cmd (mpq_ROW, mpq_NEW);
    }
    mpq_ILL_BREAK_BODY_IF (rval);

    ILLsymboltab_lookup (&lp->rowtab, rname, &ind);
    if (ind != ILL_SYM_NOINDEX) {
	rval = mpq_ILLlp_error (state, "\"%s\" is already defined.\n", rname);
	mpq_ILL_BREAK_BODY_IF (rval);
    }
    mpq_ILL_UTIL_STR (rowname, rname);

    rval = mpq_ILLread_lp_state_sense (state);
    sense = state->sense_val;
    mpq_ILL_BREAK_BODY_IF (rval);

    rval = mpq_ILLread_lp_state_value (state, &d);
    mpq_ILL_BREAK_BODY_IF (rval);

    rval = mpq_QSnew_row (p, d, sense, rowname);
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

static int mpq_del_row_or_col (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state,
      int isRow)
{
    int i[1], rval = 0;
    char **names = (isRow) ? p->qslp->rownames : p->qslp->colnames;
    int nnames = (isRow) ? p->qslp->nrows : p->qslp->nstruct;
    ILLsymboltab *tab = (isRow) ? &lp->rowtab : &lp->coltab;

    rval = mpq_ILLread_lp_state_next_field_on_line (state);
    mpq_ILL_BREAK_BODY_IF (rval);

    i[0] = mpq_ILLutil_array_index (names, nnames, state->field);
    if (i[0] >= 0) {
	rval = (isRow) ? mpq_QSdelete_rows (p, 1, i) : mpq_QSdelete_cols (p, 1, i);
	if (rval == 0) {
	    ILLsymboltab_delete (tab, state->field);
	}
    } else {
	rval = mpq_ILLlp_error (state, "\"%s\" is not defined.\n", state->field);
    }

CLEANUP:
    ILL_RESULT (rval, "mpq_del_row_or_col");
}

static void mpq_del_row (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state)
{
    int rval = mpq_del_row_or_col (p, lp, state, 1);
    if (rval == 0) {
	lp->nrows--;
    }
}

static void mpq_del_col (mpq_QSdata * p,
      mpq_rawlpdata * lp,
      mpq_ILLread_lp_state * state)
{
    int rval = mpq_del_row_or_col (p, lp, state, 0);
    if (rval == 0) {
	lp->ncols--;
    }
}
