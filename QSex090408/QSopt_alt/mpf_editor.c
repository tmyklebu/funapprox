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

/* RCS_INFO = "$RCSfile: mpf_editor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

#include "econfig.h"
#include "mpf_qsopt.h"
#include "mpf_lpdata.h"
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#include "mpf_editor.h"
#include "mpf_readline.h"
#include "mpf_rawlp.h"
#include "stddefs.h"		/* for MAX */
#include "mpf_read_lp.h"
#include "mpf_lp.h"
#include "mpf_lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

#define mpf_ILL_BREAK_BODY_IF(rval) if (rval != 0) goto CLEANUP
#define mpf_ILL_BREAK_BODY goto CLEANUP

static int mpf_transpose (mpf_rawlpdata * lp);
static int mpf_pull_info_from_p (mpf_QSdata * p,
      mpf_rawlpdata * lp);
static void mpf_add_row (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state);
/* static int new_row(mpf_QSdata *p, mpf_rawlpdata *lp, mpf_ILLread_lp_state
   state); */
static void mpf_del_row (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state);

static void mpf_add_col (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state);
static void mpf_del_col (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state);

#define mpf_NONE -1
#define mpf_QS_EXIT 0
#define mpf_ROW 1
#define mpf_COL 2
#define mpf_PLP 3
#define mpf_PRTX 4
#define mpf_SOLVE 5
#define mpf_PMPS 6
#define mpf_HELP 7
#define mpf_DEL 8
#define mpf_NEW 9
#define mpf_ADD 10
#define mpf_PRIMAL 11
#define mpf_DUAL 12
#define mpf_NCOMMAND 13
static const char *mpf_commands[mpf_NCOMMAND + 1];
static char mpf_hasSubCmd[mpf_NCOMMAND + 1];

void mpf_ILLeditor_init (void)
{
    mpf_commands[mpf_QS_EXIT] = "mpf_QS_EXIT";
    mpf_commands[mpf_ROW] = "mpf_ROW";
    mpf_commands[mpf_COL] = "mpf_COL";
    mpf_commands[mpf_PLP] = "LP";
    mpf_commands[mpf_PMPS] = "MPS";
    mpf_commands[mpf_SOLVE] = "mpf_SOLVE";
    mpf_commands[mpf_PRTX] = "PRT";
    mpf_commands[mpf_HELP] = "mpf_HELP";
    mpf_commands[mpf_ADD] = "mpf_ADD";
    mpf_commands[mpf_DEL] = "mpf_DEL";
    mpf_commands[mpf_NEW] = "mpf_NEW";
    mpf_commands[mpf_PRIMAL] = "mpf_PRIMAL";
    mpf_commands[mpf_DUAL] = "mpf_DUAL";
    mpf_commands[mpf_NCOMMAND] = NULL;

    mpf_hasSubCmd[mpf_QS_EXIT] = 0;
    mpf_hasSubCmd[mpf_ROW] = 1;
    mpf_hasSubCmd[mpf_COL] = 1;
    mpf_hasSubCmd[mpf_PLP] = 0;
    mpf_hasSubCmd[mpf_PMPS] = 0;
    mpf_hasSubCmd[mpf_SOLVE] = 1;
    mpf_hasSubCmd[mpf_PRTX] = 0;
    mpf_hasSubCmd[mpf_HELP] = 0;
    mpf_hasSubCmd[mpf_ADD] = 1;
    mpf_hasSubCmd[mpf_DEL] = 1;
    mpf_hasSubCmd[mpf_NEW] = 1;
    mpf_hasSubCmd[mpf_PRIMAL] = 1;
    mpf_hasSubCmd[mpf_DUAL] = 1;
    mpf_hasSubCmd[mpf_NCOMMAND] = 0;
}

static void mpf_ILLeditor_help_cmd (int cmd,
      int subcmd);

static void mpf_ILLeditor_help (void)
{
    mpf_ILLeditor_help_cmd (mpf_ROW, mpf_ADD);
    /* mpf_ILLeditor_help_cmd(mpf_ROW, mpf_NEW);  */
    mpf_ILLeditor_help_cmd (mpf_ROW, mpf_DEL);
    mpf_ILLeditor_help_cmd (mpf_COL, mpf_ADD);
    mpf_ILLeditor_help_cmd (mpf_COL, mpf_DEL);
    mpf_ILLeditor_help_cmd (mpf_SOLVE, mpf_NONE);
    mpf_ILLeditor_help_cmd (mpf_PRTX, mpf_NONE);
    mpf_ILLeditor_help_cmd (mpf_PLP, mpf_NONE);
    mpf_ILLeditor_help_cmd (mpf_PMPS, mpf_NONE);
    mpf_ILLeditor_help_cmd (mpf_QS_EXIT, mpf_NONE);
    mpf_ILLeditor_help_cmd (mpf_HELP, mpf_NONE);
}

static void mpf_ILLeditor_help_cmd (int cmd,
      int subcmd)
{
    if (cmd == mpf_ROW && subcmd == mpf_ADD)
	fprintf (stdout, "%s mpf_ADD:\t%s.\n",
	    mpf_commands[mpf_ROW], "add a row; enter in LP format");
    if (cmd == mpf_COL && subcmd == mpf_ADD)
	fprintf (stdout, "%s mpf_ADD:\t%s.\n",
	    mpf_commands[mpf_COL], "add a col; enter in LP format");
    /* if (cmd == mpf_ROW && subcmd == mpf_NEW) fprintf(stdout, "%s
       mpf_NEW:\t%s.\n", mpf_commands[mpf_ROW], "new row; enter rowname:
       sense rhs");  */
    if (cmd == mpf_ROW && subcmd == mpf_DEL)
	fprintf (stdout, "%s mpf_DEL:\t%s.\n",
	    mpf_commands[mpf_ROW], "delete a row; give rowname");
    if (cmd == mpf_COL && subcmd == mpf_DEL)
	fprintf (stdout, "%s mpf_DEL:\t%s.\n",
	    mpf_commands[mpf_COL], "delete a col; give colname");
    if (cmd == mpf_SOLVE)
	fprintf (stdout, "%s:\t%s.\n", mpf_commands[mpf_SOLVE], "solve problem");
    if (cmd == mpf_PRTX)
	fprintf (stdout, "%s:\t%s.\n",
	    mpf_commands[mpf_PRTX], "print variable values for optimal solution");
    if (cmd == mpf_PLP)
	fprintf (stdout, "%s [file]:\t%s.\n",
	    mpf_commands[mpf_PLP], "print problem in LP format to file or stdout");
    if (cmd == mpf_PMPS)
	fprintf (stdout, "%s [file]:\t%s.\n",
	    mpf_commands[mpf_PMPS], "print problem in MPS format to file or stdout");
    if (cmd == mpf_QS_EXIT)
	fprintf (stdout, "%s:\t%s.\n", mpf_commands[mpf_QS_EXIT], "mpf_QS_EXIT");
    if (cmd == mpf_HELP)
	fprintf (stdout, "%s:\t%s.\n", mpf_commands[mpf_HELP], "print this help");
}

static void mpf_getCmd (mpf_ILLread_lp_state * state,
      int *cmd,
      int *subcmd)
{
    const char *cmd_str, *subcmd_str;
    int tmp;

    *cmd = mpf_ILLutil_index (mpf_commands, state->field);
    *subcmd = -1;
    if (mpf_hasSubCmd[*cmd] && (mpf_ILLread_lp_state_next_field_on_line (state) == 0)) {
	*subcmd = mpf_ILLutil_index (mpf_commands, state->field);
	if ((*subcmd == mpf_ROW) || (*subcmd == mpf_COL) || (*subcmd == mpf_SOLVE)) {
	    mpf_ILL_SWAP (*subcmd, *cmd, tmp);
	}
    }
    cmd_str = (*cmd >= 0) ? mpf_commands[*cmd] : "???";
    subcmd_str = (*subcmd >= 0) ? mpf_commands[*subcmd] : "???";
    ILL_IFTRACE ("cmd = %s, subcmd = %s\n", cmd_str, subcmd_str);
}

void mpf_ILLeditor (mpf_QSdata * p)
{
    mpf_rawlpdata raw, *lp = &raw;
    int cmd, subcmd, tval, rval = 0;
    mpf_ILLread_lp_state lpstate, *state = &lpstate;
    mpf_qsline_reader *reader;

    ILL_IFTRACE ("mpf_ILLeditor\n");

    reader = mpf_ILLline_reader_new ((mpf_qsread_line_fct) fgets, stdin);
    rval = mpf_ILLread_lp_state_init (state, reader, "STDIN", 1);
    rval = rval || mpf_pull_info_from_p (p, lp);
    mpf_ILL_BREAK_BODY_IF (rval);

    while (mpf_ILLread_lp_state_next_field (state) == 0) {
	mpf_getCmd (state, &cmd, &subcmd);
	switch (cmd) {
	case mpf_QS_EXIT:
	    mpf_ILL_BREAK_BODY;

	case mpf_ROW:
	    {
		switch (subcmd) {
		case mpf_ADD:
		    mpf_add_row (p, lp, state);
		    break;
		    /* case mpf_NEW: rval = new_row(p, lp, state); break; */
		case mpf_DEL:
		    mpf_del_row (p, lp, state);
		    break;
		default:
		    mpf_ILLeditor_help ();
		    break;
		}
		break;
	    }
	case mpf_COL:
	    {
		switch (subcmd) {
		case mpf_ADD:
		    mpf_add_col (p, lp, state);
		    break;
		case mpf_DEL:
		    mpf_del_col (p, lp, state);
		    break;
		default:
		    mpf_ILLeditor_help ();
		    break;
		}
		break;
	    }

	case mpf_SOLVE:
	    {
		if (subcmd == mpf_PRIMAL) {
		    (void) mpf_ILLeditor_solve (p, PRIMAL_SIMPLEX);
		} else if (subcmd == mpf_DUAL) {
		    (void) mpf_ILLeditor_solve (p, DUAL_SIMPLEX);
		} else {
		    mpf_ILLeditor_help ();
		}
		break;
	    }

	case mpf_PRTX:
	    {
		if ((rval = mpf_ILLlib_print_x (stdout, p->lp, 0, 0, 1))) {
		    fprintf (stdout, "The problem may not be feasible.\n");
		}
		break;
	    }

	case mpf_PLP:
	case mpf_PMPS:
	    {
		if (mpf_ILLread_lp_state_next_field_on_line (state) == 0) {
		    if (cmd == mpf_PMPS) {
			tval = mpf_QSwrite_prob (p, state->field, "MPS");
		    } else {
			tval = mpf_QSwrite_prob (p, state->field, "LP");
		    }
		    if (tval) {
			fprintf (stdout, "Could not write problem to \"%s\".\n",
			    state->field);
		    } else {
			fprintf (stdout, "Saved to \"%s\".\n", state->field);
		    }
		} else {
		    if (cmd == mpf_PMPS) {
			(void) mpf_QSwrite_prob_file (p, stdout, "MPS");
		    } else {
			(void) mpf_QSwrite_prob_file (p, stdout, "LP");
		    }
		}
		break;
	    }

	case mpf_NONE:
	    fprintf (stdout, "Unknown command: %s\n", state->field);
	default:
	    mpf_ILLeditor_help ();
	    break;
	}
	fflush (stdout);
	mpf_ILLread_lp_state_next_line (state);
    }
CLEANUP:
    mpf_ILLline_reader_free (reader);
    mpf_ILLfree_rawlpdata (lp);
}

int mpf_ILLeditor_solve (mpf_QSdata * p,
      int salgo)
{
    int rval = 0;
    int status = 0;
    mpf_t val;
    mpf_EGlpNumInitVar (val);

    if (salgo == PRIMAL_SIMPLEX) {
	rval = mpf_QSopt_primal (p, &status);
    } else {
	rval = mpf_QSopt_dual (p, &status);
    }
    mpf_ILL_BREAK_BODY_IF (rval);
    rval = mpf_QSget_objval (p, &val);
    if (p->simplex_display)
	if (rval == 0) {
	    fprintf (stdout, "LP Value: %.6f, status %d\n", mpf_EGlpNumToLf (val),
		status);
	    fflush (stdout);
	}
CLEANUP:
    mpf_EGlpNumClearVar (val);
    ILL_RESULT (rval, "mpf_ILLeditor_solve");
}


static int mpf_pull_info_from_p (mpf_QSdata * p,
      mpf_rawlpdata * lp)
{
    int i, rval = 0;
    mpf_ILLlpdata *qslp = p->lp->O;
    int nrows, ncols;

    mpf_ILLinit_rawlpdata (lp, NULL);
    rval = ILLsymboltab_create (&lp->rowtab, 100) ||
	ILLsymboltab_create (&lp->coltab, 100);
    mpf_ILL_BREAK_BODY_IF (rval);

    nrows = qslp->nrows;
    ncols = qslp->nstruct;
    /* add rows to lp */
    mpf_ILLraw_add_row (lp, qslp->objname, 'N', mpf_zeroLpNum);
    for (i = 0; i < nrows; i++) {
	ILL_FAILfalse (qslp->rownames[i] != NULL, "should have no NULL names");
	mpf_ILLraw_add_row (lp, qslp->rownames[i], qslp->sense[i], qslp->rhs[i]);
    }

    /* add cols to coltab and lp */
    for (i = 0; i < ncols; i++) {
	ILL_FAILfalse (qslp->colnames[i] != NULL, "should have no NULL names");
	mpf_ILLraw_add_col (lp, qslp->colnames[i],
	    (qslp->intmarker) ? qslp->intmarker[i] : 0);
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_pull_info_from_p");
}

static int mpf_transpose (mpf_rawlpdata * lp)
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
	mpf_EGlpNumReallocArray (&(lp->rhs), lp->rhssize);
	/*
			lp->rhs = EGrealloc(lp->rhs, sizeof(double)*lp->rhssize);
			rval = ILLutil_reallocrus_scale ((void **) &lp->rhs, &lp->sensesize, tmp + 1, 1.3, sizeof (double));
			ILL_CLEANUP_IF (rval);
	*/
    }
    mpf_ILL_SWAP (lp->nrows, lp->ncols, tmp);
    mpf_ILL_SWAP (lp->rowtab, lp->coltab, tmptab);
    ILL_RETURN (rval, "mpf_transpose");
}

static char *mpf_get_row_col_name (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state,
      int doRow)
{
    int rval = 0;
    int ind;
    char *rname, *thename = NULL;
    char buf[ILL_namebufsize];
    ILLsymboltab *tab = (doRow) ? &lp->rowtab : &lp->coltab;
    int id = (doRow) ? lp->nrows : lp->ncols;
    id--;			/* in mpf_rawlpdata obj counts as a row */

    rval = mpf_ILLread_constraint_name (state, &rname);
    mpf_ILL_BREAK_BODY_IF (rval);

    if (rname == NULL) {
	mpf_ILLlib_findName (p->qslp, doRow /* forRow */ , rname, id, buf);
	mpf_ILL_UTIL_STR (thename, buf);
    } else {
	mpf_ILL_UTIL_STR (thename, rname);
    }
    if (ILLsymboltab_lookup (tab, thename, &ind) == 0) {
	rval = mpf_ILLlp_error (state, "\"%s\" already exists.", thename);
    }
CLEANUP:
    if (rval != 0) {
	ILL_IFFREE (thename, char);
    }
    return thename;
}

static int mpf_fill_matrix (mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state,
      mpf_ILLmatrix * m,
      mpf_t * obj,
      int n)
{
    int i, cnt, rval = 0;
    mpf_colptr *cp;
    mpf_t val;
    int newCol = (obj != NULL);
    mpf_EGlpNumInitVar (val);

    /* rely on fact that objective has rowindex 0 */

    m->matrows = lp->nrows;
    m->matcols = 1;
    m->matval = mpf_EGlpNumAllocArray (lp->ncols);
    ILL_SAFE_MALLOC (m->matind, lp->ncols, int);
    ILL_SAFE_MALLOC (m->matbeg, 1, int);
    ILL_SAFE_MALLOC (m->matcnt, 1, int);
    m->matsize = lp->ncols;
    m->matbeg[0] = 0;
    m->matcnt[0] = 0;
    for (i = 0; i < lp->ncols; i++) {
	cnt = 0;
	mpf_EGlpNumZero (val);
	for (cp = lp->cols[i]; cp != NULL; cp = cp->next) {
	    ILL_FAILfalse (cp->this == n, "n should be the only row around");
	    if (mpf_EGlpNumIsNeqq (cp->coef, mpf_zeroLpNum)) {
		mpf_EGlpNumAddTo (val, cp->coef);
		cnt++;
	    }
	}
	if (cnt > 1) {
	    mpf_ILLlp_warn (state, "Multiple coefficients for \"%s\".",
		mpf_ILLraw_colname (lp, i));
	}
	if (mpf_EGlpNumIsNeqqZero (val)) {
	    if ((i - newCol) >= 0) {
		mpf_EGlpNumCopy (m->matval[m->matcnt[0]], val);
		m->matind[m->matcnt[0]] = i - newCol;
		m->matcnt[0]++;
	    } else {
		mpf_EGlpNumCopy (obj[0], val);
	    }
	}
    }
    if (m->matcnt[0] == 0) {
	rval = mpf_ILLlp_error (state, "There are no non zero coefficients.");
    }
CLEANUP:
    mpf_EGlpNumClearVar (val);
    ILL_RESULT (rval, "mpf_fill_matrix");
}

static void mpf_add_row (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state)
{
    int rval = 0;
    int n;
    char *name;
    mpf_ILLmatrix m;
    char sense[1];

    mpf_ILLmatrix_init (&m);
    n = lp->nrows;
    name = mpf_get_row_col_name (p, lp, state, 1 /* doRow */ );

    if (name == NULL) {
	rval = 1;
    } else {
	rval = mpf_ILLread_one_constraint (state, name, lp, 0);

	/* adds row name to lp->rowtab the checks constraint expression  */
	if (rval != 0) {
	    /* failed because of error in expression => must remove name from
	       symbol table */
	    fprintf (stdout, "Incorrect expression.\n");
	} else {
	    ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");
	    ILL_IFTRACE ("ADDING row %s.\n", name);

	    sense[0] = lp->rowsense[n];

	    rval = mpf_fill_matrix (lp, state, &m, NULL, n);
	    mpf_ILL_BREAK_BODY_IF (rval);

	    mpf_QSadd_rows (p, 1, m.matcnt, m.matbeg, m.matind, m.matval,
		&(lp->rhs[n]), sense, (const char **) &name);
	}
    }
CLEANUP:
    mpf_ILLmatrix_free (&m);
    if (name != NULL) {
	if (rval != 0)
	    ILLsymboltab_delete (&lp->rowtab, name);
	ILL_IFFREE (name, char);
    }
    if (rval != 0) {
	lp->nrows = n;
    }
    mpf_ILLraw_clear_matrix (lp);
}

static void mpf_add_col (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state)
{
    int rval = 0;
    int n;
    char *name[1];
    int transposed = 1;
    mpf_ILLmatrix matrix, *m = &matrix;
    mpf_t obj[1], lower[1], upper[2];
    mpf_EGlpNumInitVar (*obj);
    mpf_EGlpNumInitVar (*lower);
    mpf_EGlpNumInitVar (upper[0]);
    mpf_EGlpNumInitVar (upper[1]);

    n = lp->ncols;
    mpf_ILLmatrix_init (m);
    name[0] = mpf_get_row_col_name (p, lp, state, 0 /* doRow */ );
    rval = (name[0] == NULL);
    mpf_ILL_BREAK_BODY_IF (rval);

    transposed = !mpf_transpose (lp);
    rval = mpf_ILLread_one_constraint (state, name[0], lp, 0);

    /* adds row name to lp->rowtab the checks constraint expression  */
    if (rval != 0) {
	/* failed because of error in expression => must remove name from
	   symbol table */
	fprintf (stdout, "Incorrect expression.\n");
    } else {
	ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");

	rval = mpf_fill_matrix (lp, state, m, obj, n);
	mpf_ILL_BREAK_BODY_IF (rval);

	fprintf (stdout, "lower ");
	rval = mpf_ILLread_lp_state_next_line (state) ||
	    mpf_ILLread_lp_state_value (state, &(lower[0]));
	mpf_ILL_BREAK_BODY_IF (rval);

	fprintf (stdout, "upper ");
	rval = mpf_ILLread_lp_state_next_line (state) ||
	    mpf_ILLread_lp_state_value (state, &(upper[0]));
	mpf_ILL_BREAK_BODY_IF (rval);

	ILL_IFTRACE ("ADDING col %s.\n", name[0]);

	mpf_QSadd_cols (p, 1, m->matcnt, m->matbeg, m->matind, m->matval,
	    obj, lower, upper, (const char **) name);

    }
CLEANUP:
    mpf_ILLmatrix_free (m);
    if (name[0] != NULL) {
	if (rval != 0)
	    ILLsymboltab_delete (&lp->rowtab, name[0]);
	ILL_IFFREE (name[0], char);
    }
    if (rval != 0) {
	lp->nrows = n;
    }
    mpf_ILLraw_clear_matrix (lp);
    if (transposed)
	mpf_transpose (lp);
    ILL_IFFREE (name[0], char);
    mpf_EGlpNumClearVar (*obj);
    mpf_EGlpNumClearVar (*lower);
    mpf_EGlpNumClearVar (upper[0]);
    mpf_EGlpNumClearVar (upper[1]);
}

#if 0
#ifndef JAVA_PORT
static void new_row (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state)
{
    int rval = 0;
    char *rowname = NULL, *rname = NULL;
    char sense;
    double d;
    int ind, hit;

    rval = mpf_ILLread_constraint_name (state, &rname);
    if (rname == NULL) {
	rval = 1;
	mpf_ILLeditor_help_cmd (mpf_ROW, mpf_NEW);
    }
    mpf_ILL_BREAK_BODY_IF (rval);

    ILLsymboltab_lookup (&lp->rowtab, rname, &ind);
    if (ind != ILL_SYM_NOINDEX) {
	rval = mpf_ILLlp_error (state, "\"%s\" is already defined.\n", rname);
	mpf_ILL_BREAK_BODY_IF (rval);
    }
    mpf_ILL_UTIL_STR (rowname, rname);

    rval = mpf_ILLread_lp_state_sense (state);
    sense = state->sense_val;
    mpf_ILL_BREAK_BODY_IF (rval);

    rval = mpf_ILLread_lp_state_value (state, &d);
    mpf_ILL_BREAK_BODY_IF (rval);

    rval = mpf_QSnew_row (p, d, sense, rowname);
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

static int mpf_del_row_or_col (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state,
      int isRow)
{
    int i[1], rval = 0;
    char **names = (isRow) ? p->qslp->rownames : p->qslp->colnames;
    int nnames = (isRow) ? p->qslp->nrows : p->qslp->nstruct;
    ILLsymboltab *tab = (isRow) ? &lp->rowtab : &lp->coltab;

    rval = mpf_ILLread_lp_state_next_field_on_line (state);
    mpf_ILL_BREAK_BODY_IF (rval);

    i[0] = mpf_ILLutil_array_index (names, nnames, state->field);
    if (i[0] >= 0) {
	rval = (isRow) ? mpf_QSdelete_rows (p, 1, i) : mpf_QSdelete_cols (p, 1, i);
	if (rval == 0) {
	    ILLsymboltab_delete (tab, state->field);
	}
    } else {
	rval = mpf_ILLlp_error (state, "\"%s\" is not defined.\n", state->field);
    }

CLEANUP:
    ILL_RESULT (rval, "mpf_del_row_or_col");
}

static void mpf_del_row (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state)
{
    int rval = mpf_del_row_or_col (p, lp, state, 1);
    if (rval == 0) {
	lp->nrows--;
    }
}

static void mpf_del_col (mpf_QSdata * p,
      mpf_rawlpdata * lp,
      mpf_ILLread_lp_state * state)
{
    int rval = mpf_del_row_or_col (p, lp, state, 0);
    if (rval == 0) {
	lp->ncols--;
    }
}
