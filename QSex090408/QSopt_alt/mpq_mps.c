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

/* RCS_INFO = "$RCSfile: mpq_mps.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

/****************************************************************************/
/* */
/* Routines for Reading and Writing MPS Files                 */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int  ILLlpdata_mpsread (mpq_ILLlpdata *lp, const char *filename);         */
/* int  ILLlpdata_mpswrite(mpq_ILLlpdata *lp, const char *filename);         */
/* */
/* NOTES                                                                   */
/* */
/* In the MPS reader, integer variables without an explicit bound are    */
/* set to binary; real variables without an explict lower bound and      */
/* either a nonnegative or free upperbound set to nonnegative.  (These   */
/* are standard settings.)                                               */
/* */
/* If a RHS is not specified for a row, it is set to 0.                  */
/* */
/* The MPS reader allows CPLEX's OBJSENSE extension to specify max or    */
/* min and the OBJNAME extension to specify an objective row.            */
/* */
/****************************************************************************/

#include "econfig.h"
#include "mpq_iqsutil.h"
#include "mpq_mps.h"
#include "mpq_rawlp.h"
#include "mpq_lpdata.h"
/* extern mpq_t mpq_SZERO_TOLER; */

static int TRACE = 0;

const char *mpq_ILLmps_section_name[ILL_MPS_N_SECTIONS + 2] = {
    "NAME", "OBJSENSE", "OBJNAME", "ROWS", "COLUMNS",
    "RHS", "RANGES", "BOUNDS", "REFROW", "ENDATA",
    NULL
};

static const char *mpq_mps_bound_name[] = {
    "LO", "UP", "FX", "FR", "MI", "PL", "BV", "UI", "LI", NULL
};

/****************************************************************************/
/* reading */
static int mpq_read_mps_section (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);

static int mpq_read_mps_name (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_read_mps_refrow (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_read_mps_objnamesense (ILLmps_section sec,
      mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_read_mps_objname (mpq_ILLread_mps_state * state);
static int mpq_read_mps_objsense (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);

static int mpq_read_mps_line_in_section (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);


static int mpq_add_row (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_add_col (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_add_rhs (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_add_ranges (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_add_bounds (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);

static int mpq_mps_read_marker_line (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);
static int mpq_is_marker_line (mpq_ILLread_mps_state * state);
static int mpq_mps_read_col_line (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp);

static int mpq_mps_fill_in (mpq_rawlpdata * lp,
      const char *obj);


static void mpq_mps_set_bound (mpq_rawlpdata * lp,
      mpq_ILLread_mps_state * state,
      int colind,
      const char *bndtype,
      mpq_t bnd);

int mpq_ILLread_mps (mpq_qsline_reader * file,
      const char *f,
      mpq_rawlpdata * lp)
{
    int rval = 0;
    int end = 0;
    mpq_ILLread_mps_state state;

    ILL_IFTRACE ("\tread_mps\n");
    if (ILLsymboltab_create (&lp->rowtab, 100) ||
	ILLsymboltab_create (&lp->coltab, 100)) {
	rval = 1;
    } else {
	rval = mpq_ILLmps_state_init (&state, file, f);
    }
    if (rval != 0) {
	goto CLEANUP;
    }
    while (mpq_ILLmps_next_line (&state) == 0) {
	if (mpq_ILLmps_empty_key (&state)) {
	    if (mpq_read_mps_line_in_section (&state, lp) != 0) {
		rval++;
	    }
	} else {
	    /* found a section indicator in col 1 */
	    if (!strcmp (state.key, mpq_ILLmps_section_name[ILL_MPS_ENDATA])) {
		end = 1;
		break;		/* done reading */
	    }
	    if (mpq_read_mps_section (&state, lp) != 0) {
		rval++;
	    }
	}
	if (rval == 50) {
	    (void) mpq_ILLmps_error (&state, "Too many errors.\n");
	}
    }

    if (!end) {
	mpq_ILLmps_warn (&state, "Missing ENDATA.");
    }
    if (!mpq_ILLmps_next_line (&state)) {
	mpq_ILLmps_warn (&state, "Ignoring text after ENDATA.");
    }
    if (rval == 0) {
	rval = mpq_mps_fill_in (lp, state.obj);
    }
CLEANUP:
    ILL_RESULT (rval, "read_mps");
}

static int mpq_check_section_order (mpq_ILLread_mps_state * state,
      int sec)
{
    switch (sec) {
    case ILL_MPS_REFROW:
	if (state->section[ILL_MPS_ROWS] == 1) {
	    return mpq_ILLmps_error (state, "%s section after ROWS section.\n",
		mpq_ILLmps_section_name[sec]);
	}
	break;

    case ILL_MPS_COLS:
    case ILL_MPS_RHS:
    case ILL_MPS_RANGES:
	if (state->section[ILL_MPS_ROWS] == 0) {
	    return mpq_ILLmps_error (state, "%s section before ROWS section.\n",
		mpq_ILLmps_section_name[sec]);
	};
	break;

    case ILL_MPS_BOUNDS:
	if (state->section[ILL_MPS_COLS] == 0) {
	    return mpq_ILLmps_error (state, "%s section before COLUMNS section.\n",
		mpq_ILLmps_section_name[sec]);
	}
	break;
    }
    return 0;
}

static int mpq_read_mps_section (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int sec;
    int rval = 0, r;

    ILL_FAILtrue (mpq_ILLmps_empty_key (state), "must have a key on this line");

    sec = mpq_ILLutil_index (mpq_ILLmps_section_name, state->key);
    if (sec < 0) {
	return mpq_ILLmps_error (state, "\"%s\" is not a key.\n", state->key);
    }
    rval = mpq_ILLmps_set_section (state, sec);

    state->active = ILL_MPS_NONE;
    rval = rval || mpq_check_section_order (state, sec);
    switch (sec) {
    case ILL_MPS_COLS:
    case ILL_MPS_ROWS:
	state->active = sec;
	break;

    case ILL_MPS_NAME:
	if (rval == 0) {
	    rval = mpq_read_mps_name (state, lp);
	}
	break;

    case ILL_MPS_RHS:
	if (rval == 0) {
	    rval = mpq_ILLraw_init_rhs (lp);
	}
	state->active = ILL_MPS_RHS;
	break;

    case ILL_MPS_RANGES:
	if (rval == 0) {
	    rval = mpq_ILLraw_init_ranges (lp);
	}
	state->active = ILL_MPS_RANGES;
	break;

    case ILL_MPS_BOUNDS:
	if (rval == 0) {
	    rval = mpq_ILLraw_init_bounds (lp);
	}
	state->active = ILL_MPS_BOUNDS;
	break;

    case ILL_MPS_OBJNAME:
    case ILL_MPS_OBJSENSE:
	r = mpq_read_mps_objnamesense (sec, state, lp);
	rval = r || rval;
	break;

    case ILL_MPS_REFROW:
	r = mpq_read_mps_refrow (state, lp);
	rval = r || rval;
	break;

    default:
	ILL_REPRT ("should never get here");
	goto CLEANUP;
    }
CLEANUP:
    ILL_RESULT (rval, "mpq_read_mps_section");
}

static int mpq_read_mps_name (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rval = 0;

    if (mpq_ILLmps_empty_field (state)) {
	mpq_ILLmps_warn (state, "Blank NAME.");
    } else {
	mpq_ILL_UTIL_STR (lp->name, state->field);
    }
CLEANUP:
    ILL_RESULT (rval, "mpq_read_mps_name");
}

static int mpq_read_mps_refrow (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rval = 0;

    rval = mpq_ILLmps_next_line (state);
    if (state->section[ILL_MPS_REFROW] > 1) {
	/* this is the second time we see this section; don't complain about
	   errors */
	return 0;
    }
    if (mpq_ILLmps_empty_key (state) && !mpq_ILLmps_empty_field (state)) {
	mpq_ILL_UTIL_STR (lp->refrow, state->field);
	return 0;
    } else {
	return mpq_ILLmps_error (state, "Bad row name in REFROW section.\n");
    }
CLEANUP:
    ILL_RETURN (rval, "mpq_read_mps_refrow");
}

static int mpq_read_mps_objnamesense (ILLmps_section sec,
      mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    if (state->section[sec] > 1) {
	/* this is the second time we see this section; just skip next line */
	(void) mpq_ILLmps_next_line (state);
	return 0;
    }
    if (mpq_ILLmps_next_line (state) != 0) {
	return mpq_ILLmps_error (state, "Missing %s line at end of file.\n",
	    mpq_ILLmps_section_name[sec]);
    }
    if ((!mpq_ILLmps_empty_key (state)) || mpq_ILLmps_empty_field (state)) {
	(void) mpq_ILLmps_error (state, "Bad %s in %s record.\n",
	    ((sec == ILL_MPS_OBJNAME) ? "row name"
		: "objective sense"), mpq_ILLmps_section_name[sec]);
	if (!mpq_ILLmps_empty_key (state)) {
	    (void) mpq_read_mps_section (state, lp);
	}
	return 1;
    }
    if (sec == ILL_MPS_OBJNAME) {
	if (mpq_read_mps_objname (state)) {
	    return 1;
	}
    } else {
	if (mpq_read_mps_objsense (state, lp) != 0) {
	    return 1;
	}
    }
    return 0;
}

static int mpq_read_mps_objsense (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rval = 0;
    char *objsense = state->field;
    ILL_FAILfalse (state->section[ILL_MPS_OBJSENSE] == 1, "should never happen");
    if (!strcmp (objsense, "MAX") ||
	!strcmp (objsense, "Max") ||
	!strcmp (objsense, "max") ||
	!strcmp (objsense, "MAXIMIZE") ||
	!strcmp (objsense, "Maximize") || !strcmp (objsense, "maximize")) {
	lp->objsense = mpq_ILL_MAX;
    } else if (!strcmp (objsense, "MIN") ||
	    !strcmp (objsense, "Min") ||
	    !strcmp (objsense, "min") ||
	    !strcmp (objsense, "MINIMIZE") ||
	!strcmp (objsense, "Minimize") || !strcmp (objsense, "minimize")) {
	lp->objsense = mpq_ILL_MIN;
    } else {
	return mpq_ILLmps_error (state, "\"%s\" is no OBJSENSE.\n", objsense);
    }
CLEANUP:
    ILL_RESULT (rval, "mpq_read_mps_objsense");
}

static int mpq_read_mps_objname (mpq_ILLread_mps_state * state)
{
    int rval = 0;
    ILL_FAILfalse (state->section[ILL_MPS_OBJNAME] == 1, "should never happen");
    mpq_ILL_UTIL_STR (state->obj, state->field);
CLEANUP:
    ILL_RETURN (rval, "mpq_read_mps_objname");
}

static int mpq_read_mps_line_in_section (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rval = 0;

    ILL_FAILtrue (!mpq_ILLmps_empty_key (state) || mpq_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    if (state->active == ILL_MPS_NONE) {
	return mpq_ILLmps_error (state, "Line is in no section.\n");
    } else {
	if (state->section[state->active] == 1) {
	    switch (state->active) {
	    case ILL_MPS_ROWS:
		rval = mpq_add_row (state, lp);
		break;
	    case ILL_MPS_COLS:
		rval = mpq_add_col (state, lp);
		break;
	    case ILL_MPS_RHS:
		rval = mpq_add_rhs (state, lp);
		break;
	    case ILL_MPS_RANGES:
		rval = mpq_add_ranges (state, lp);
		break;
	    case ILL_MPS_BOUNDS:
		rval = mpq_add_bounds (state, lp);
		break;
	    default:
		ILL_REPRT ("should never get here");
		ILL_CLEANUP;
	    }
	    if (rval == 0) {
		/* see whether there are extra fields on line */
		mpq_ILLmps_check_end_of_line (state);
	    }
	}
    }
CLEANUP:
    ILL_RESULT (rval, "mpq_read_mps_line_in_section");
}

static int mpq_add_row (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int ind, hit, rval = 0;
    char sense;

    ILL_FAILtrue (!mpq_ILLmps_empty_key (state) || mpq_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    /* field should contain exactly one character */
    if (state->field[1] == '\0') {
	sense = state->field[0];
	if (sense != 'L' && sense != 'G' && sense != 'E' && sense != 'N') {
	    return mpq_ILLmps_error (state,
		"Unknown rowsense '%c' in ROWS record.\n", sense);
	}
	if (mpq_ILLmps_next_field (state) == 0) {
	    hit = ILLsymboltab_lookup (&lp->rowtab, state->field, &ind);
	    if (!hit) {
		rval = mpq_ILLmps_error (state,
		    "Repeated row definition for \"%s\".\n",
		    state->field);
	    } else {
		rval = mpq_ILLraw_add_row (lp, state->field, sense, mpq_zeroLpNum);	/* default rhs is 0.0 */
	    }
	} else {
	    rval = mpq_ILLmps_error (state, "Missing rowname in ROWS record.\n");
	}
    } else {
	rval = mpq_ILLmps_error (state, "Unknown rowsense '%s' in ROWS record.\n",
	    state->field);
    }
CLEANUP:
    ILL_RESULT (rval, "mpq_add_row");
}

static int mpq_add_col (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rval = 0;

    ILL_FAILtrue (!mpq_ILLmps_empty_key (state) || mpq_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    if (mpq_is_marker_line (state)) {
	return mpq_mps_read_marker_line (state, lp);
    } else {
	return mpq_mps_read_col_line (state, lp);
    }
CLEANUP:
    ILL_RETURN (rval, "mpq_add_col");
}

static int mpq_mps_read_col_line (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int hit, colind, rowind, rval = 0, more, ind;
    mpq_t ncoef;
    mpq_EGlpNumInitVar (ncoef);

    ILL_FAILtrue (!mpq_ILLmps_empty_key (state) || mpq_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    hit = ILLsymboltab_lookup (&lp->coltab, state->field, &colind);
    if (hit) {
	rval = mpq_ILLraw_add_col (lp, state->field, state->intvar);
	ILL_CLEANUP_IF (rval);
	colind = lp->ncols - 1;
    } else {
	if (state->intvar) {
	    /* previously declared variable is in integer section */
	    lp->intmarker[colind] = 1;
	}
    }
#ifndef NDEBUG
    hit = ILLsymboltab_lookup (&lp->coltab, state->field, &ind);
    ILL_FAILfalse (colind == ind, "colind should be index of state->field");
#endif
    if (state->sosvar == 1) {
	if (mpq_ILLraw_is_mem_other_sos (lp, colind)) {
	    rval = mpq_ILLmps_error (state,
		"\"%s\" is a member of SOS set #%d.\n",
		mpq_ILLraw_colname (lp, colind),
		lp->is_sos_member[colind] + 1);
	} else {
	    rval = mpq_ILLraw_add_sos_member (lp, colind);
	}
	ILL_CLEANUP_IF (rval);
    }
    more = (mpq_ILLmps_next_field (state) == 0);
    if (!more) {
	return mpq_ILLmps_error (state, "Missing fields in COLUMNS record.\n");
    }
    for (more = 1; more; more = (mpq_ILLmps_next_field (state) == 0)) {
	hit = ILLsymboltab_lookup (&lp->rowtab, state->field, &rowind);
	if (hit) {
	    return mpq_ILLmps_error (state, "\"%s\" is not a row name.\n", state->field);
	}
	if (mpq_ILLmps_next_coef (state, &ncoef) != 0) {
	    return mpq_ILLmps_error (state,
		"Missing/Bad coefficient in COLUMNS record.\n");
	}
	rval = mpq_ILLraw_add_col_coef (lp, colind, rowind, ncoef);
    }
CLEANUP:
    mpq_EGlpNumClearVar (ncoef);
    ILL_RESULT (rval, "mpq_mps_read_col_line");
}

static int mpq_is_marker_line (mpq_ILLread_mps_state * state)
{
    const char *field = state->line;
    while ((field = mpq_ILLutil_strchr (field, '\''))) {
	if (strncmp (field, "'MARKER'", (size_t) 8) == 0) {
	    return 1;
	}
	while ((!mpq_ILL_ISBLANK (field)) && (*field != '\0')) {
	    field++;
	}
    }
    return 0;
}

static int mpq_mps_read_marker_line (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rval = 0;
    int sos_type = mpq_ILL_SOS_TYPE1;
    int sos_line = 0;
    int cur_sos_mode = state->sosvar;

    if (strcmp (state->field, "S2") == 0) {
	sos_type = mpq_ILL_SOS_TYPE2;
	sos_line = 1;
    } else if (strcmp (state->field, "S1") == 0) {
	sos_line = 1;
    }
    if (sos_line) {
	rval = mpq_ILLmps_next_field (state);
    }
    rval = rval || mpq_ILLmps_next_field (state);	/* swallow marker-id */
    if (strcmp (state->field, "'MARKER'")) {
	return mpq_ILLmps_error (state, "Bad 'MARKER' line.\n");
    }
    if (mpq_ILLmps_next_field (state)) {
	return mpq_ILLmps_error (state, "Missing field on 'MARKER' line.\n");
    }
    rval = mpq_ILLmps_int_sos_mode (state);
    if (rval == 0) {
	if (cur_sos_mode != state->sosvar) {
	    if (state->sosvar) {
		/* beginning of a new set */
		rval = mpq_ILLraw_add_sos (lp, sos_type);
	    }
	}
    }
    ILL_RESULT (rval, "mpq_mps_read_marker_line");
}

static int mpq_add_rhs (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    int rowind, more_fields, skip;
    const char *rhsname;
    mpq_t rhs;
    mpq_EGlpNumInitVar (rhs);

    rhsname = mpq_ILLmps_possibly_blank_name (state->field, state, &lp->rowtab);
    if (mpq_ILLraw_set_rhs_name (lp, rhsname, &skip)) {
	mpq_ILLmps_error (state, "Could not add right hand side.\n");
    }
    if (skip) {
	mpq_ILLmps_set_end_of_line (state);	/* to avoid warning about
						   extra fields */
    } else {
	if (strcmp (rhsname, " ")) {
	    /* field is non blank rhs name; advance to row name  */
	    if (mpq_ILLmps_next_field (state)) {
		return mpq_ILLmps_error (state, "Missing row name in RHS record.\n");
	    }
	}
	for (more_fields = 1; more_fields; more_fields = !mpq_ILLmps_next_field (state)) {
	    if (ILLsymboltab_lookup (&lp->rowtab, state->field, &rowind)) {
		return mpq_ILLmps_error (state, "\"%s\" is not a row name.\n",
		    state->field);
	    }
	    if (mpq_ILLmps_next_coef (state, &rhs)) {
		return mpq_ILLmps_error (state, "Missing/Bad coefficient in RHS record.\n");
	    }
	    if (lp->rhsind[rowind]) {
		return mpq_ILLmps_error (state, "Two rhs values for row \"%s\".\n",
		    state->field);
	    } else {
		if (lp->rowsense[rowind] == 'N') {
		    mpq_ILLmps_warn (state,
			"Ignoring right hand side for N-row \"%s\".",
			mpq_ILLraw_rowname (lp, rowind));
		} else {
		    mpq_EGlpNumCopy (lp->rhs[rowind], rhs);
		    lp->rhsind[rowind] = 1;
		}
	    }
	}
    }
    mpq_EGlpNumClearVar (rhs);
    return 0;
}

static int mpq_add_bounds (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    char bndtype[3];
    const char *bounds_name;
    int colind, skip, rval = 0;
    mpq_t bnd;
    mpq_EGlpNumInitVar (bnd);

    ILL_FAILtrue (!mpq_ILLmps_empty_key (state) || mpq_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    if (mpq_ILLutil_index (mpq_mps_bound_name, state->field) < 0) {
	return mpq_ILLmps_error (state, "\"%s\" is not a BOUNDS type.\n", state->field);
    }
    strcpy (bndtype, state->field);

    if (mpq_ILLmps_next_field (state) != 0) {
	return mpq_ILLmps_error (state,
	    "No bounds/column identifier in BOUNDS record.\n");
    }
    bounds_name = mpq_ILLmps_possibly_blank_name (state->field, state, &lp->coltab);
    if (bounds_name == NULL) {
	return 1;
    }
    if (mpq_ILLraw_set_bounds_name (lp, bounds_name, &skip)) {
	return 1;
    }
    if (skip) {
	mpq_ILLmps_set_end_of_line (state);	/* to avoid warning about
						   extra fields */
    } else {
	if (strcmp (bounds_name, " ")) {
	    /* non empty bounds_name ==> advance to col name field */
	    if (mpq_ILLmps_next_field (state)) {
		return mpq_ILLmps_error (state, "Missing column field in BOUNDS record.\n");
	    }
	}
	if (ILLsymboltab_lookup (&lp->coltab, state->field, &colind)) {
	    return mpq_ILLmps_error (state, "\"%s\" is not a column name.\n",
		state->field);
	}
	mpq_EGlpNumZero (bnd);
	if (strcmp (bndtype, "FR") && strcmp (bndtype, "BV") &&
	    strcmp (bndtype, "MI") && strcmp (bndtype, "PL")) {
	    /* neither "FR", "BV", "MI" nor "PL" ==> there should be a bound */
	    if (mpq_ILLmps_next_bound (state, &bnd)) {
		return mpq_ILLmps_error (state,
		    "Missing/Bad bound field in BOUNDS record.\n");
	    }
	}
	mpq_mps_set_bound (lp, state, colind, bndtype, bnd);
    }
CLEANUP:
    mpq_EGlpNumClearVar (bnd);
    ILL_RESULT (rval, "mpq_add_bounds");
}

static void mpq_mps_set_bound (mpq_rawlpdata * lp,
      mpq_ILLread_mps_state * state,
      int colind,
      const char *bndtype,
      mpq_t bnd)
{
    const char *msg = NULL;
    if (!strcmp (bndtype, "LO")) {
	msg = mpq_ILLraw_set_lowerBound (lp, colind, bnd);
    } else if (!strcmp (bndtype, "UP")) {
	msg = mpq_ILLraw_set_upperBound (lp, colind, bnd);
    } else if (!strcmp (bndtype, "FX")) {
	msg = mpq_ILLraw_set_fixedBound (lp, colind, bnd);
    } else if (!strcmp (bndtype, "FR")) {
	msg = mpq_ILLraw_set_unbound (lp, colind);
    } else if (!strcmp (bndtype, "BV")) {
	msg = mpq_ILLraw_set_binaryBound (lp, colind);
	if (msg == NULL) {
	    lp->intmarker[colind] = 1;
	}
    } else if (!strcmp (bndtype, "UI")) {
	msg = mpq_ILLraw_set_upperBound (lp, colind, bnd);
	if (msg == NULL) {
	    lp->intmarker[colind] = 1;
	}
    } else if (!strcmp (bndtype, "LI")) {
	msg = mpq_ILLraw_set_lowerBound (lp, colind, bnd);
	if (msg == NULL) {
	    lp->intmarker[colind] = 1;
	}
    } else if (!strcmp (bndtype, "MI")) {
	msg = mpq_ILLraw_set_lowerBound (lp, colind, mpq_ILL_MINDOUBLE);
    } else if (!strcmp (bndtype, "PL")) {
	msg = mpq_ILLraw_set_upperBound (lp, colind, mpq_ILL_MAXDOUBLE);
    } else {
	ILL_REPRT ("should never get here");
	ILL_CLEANUP;
    }
    mpq_ILLmps_warn (state, msg);
CLEANUP:
    return;
}

static int mpq_add_ranges (mpq_ILLread_mps_state * state,
      mpq_rawlpdata * lp)
{
    const char *rangesname;
    int skip, more_fields, rowind;
    mpq_t ntmp;
    mpq_EGlpNumInitVar (ntmp);

    rangesname = mpq_ILLmps_possibly_blank_name (state->field, state, &lp->rowtab);
    if (mpq_ILLraw_set_ranges_name (lp, rangesname, &skip)) {
	return mpq_ILLmps_error (state, "Could not add range.\n");
    }
    if (skip) {
	mpq_ILLmps_set_end_of_line (state);	/* to avoid warning about
						   extra fields */
    } else {
	if (strcmp (rangesname, " ")) {
	    /* field is non blank ranges name; advance to row name */
	    if (mpq_ILLmps_next_field (state)) {
		return mpq_ILLmps_error (state, "Missing row name in RANGES record.");
	    }
	}
	for (more_fields = 1; more_fields; more_fields = !mpq_ILLmps_next_field (state)) {
	    if (ILLsymboltab_lookup (&lp->rowtab, state->field, &rowind)) {
		return mpq_ILLmps_error (state, "\"%s\" is not a row name.\n",
		    state->field);
	    }
	    if (mpq_ILLmps_next_coef (state, &ntmp)) {
		return mpq_ILLmps_error (state,
		    "Missing/Bad coefficient in RANGES record.\n");
	    }
	    if (lp->rangesind[rowind]) {
		mpq_ILLmps_warn (state, "Ignoring second RANGE value %s \"%s\".",
		    "for row", mpq_ILLraw_rowname (lp, rowind));
	    } else {
		if (lp->rowsense[rowind] != 'N') {
		    if (mpq_ILLraw_add_ranges_coef (lp, rowind, ntmp))
			return 1;
		} else {
		    mpq_ILLmps_warn (state, "Ignoring RANGE value for N-row \"%s\".",
			mpq_ILLraw_rowname (lp, rowind));
		}
	    }
	}
    }
    mpq_EGlpNumClearVar (ntmp);
    return 0;
}


static int mpq_mps_fill_in (mpq_rawlpdata * lp,
      const char *obj)
{
    int i, hit, rval = 0;

    /* find the objective function -- the first N row if obj is not defined  */
    if (obj) {
	hit = ILLsymboltab_lookup (&lp->rowtab, obj, &lp->objindex);
	if (hit) {
	    return mpq_ILLdata_error (lp->error_collector,
		"Bad objective name \"%s\".\n", obj);
	} else if (lp->rowsense[lp->objindex] != 'N') {
	    mpq_ILLdata_warn (lp->error_collector,
		"Making objective row \"%s\" a N-row.", obj);
	}
	lp->rowsense[lp->objindex] = 'N';
    } else {
	for (i = 0; i < lp->nrows; i++) {
	    if (lp->rowsense[i] == 'N') {
		lp->objindex = i;
		break;
	    }
	}
	if (i == lp->nrows) {
	    return mpq_ILLdata_error (lp->error_collector,
		"No N-row in lp definition.\n");
	}
    }

    if (lp->ncols > 0) {
	rval = mpq_ILLraw_fill_in_bounds (lp);
    }
    /* set weights of sos set members */
    if (lp->refrow) {
	/* take the values from refrow */
	mpq_t weight;
	mpq_colptr *cp;
	mpq_EGlpNumInitVar (weight);

	hit = ILLsymboltab_lookup (&lp->rowtab, lp->refrow, &lp->refrowind);
	if (hit) {
	    return mpq_ILLdata_error (lp->error_collector,
		"REFROW \"%s\" is not a row name.\n", lp->refrow);
	}
	for (i = 0; i < lp->nsos_member; i++) {
	    for (cp = lp->cols[lp->sos_col[i]]; cp != NULL; cp = cp->next) {
		if (cp->this == lp->refrowind)
		    break;
	    }
	    if ((cp != NULL) && mpq_EGlpNumIsNeqqZero (cp->coef)) {
		mpq_EGlpNumCopy (weight, cp->coef);
	    } else {
		mpq_ILLdata_warn (lp->error_collector,
		    "\"%s\" has 0.0 coefficient in REFROW \"%s\".",
		    mpq_ILLraw_colname (lp, lp->sos_col[i]), lp->refrow);
		mpq_EGlpNumZero (weight);
	    }
	    mpq_EGlpNumCopy (lp->sos_weight[i], weight);
	}
	mpq_EGlpNumClearVar (weight);
    } else {			/* no refrow */
	/* set weights to 1, 2, 3, ... in order of definition */
	int si, w;
	mpq_sosptr *set;
	for (si = 0; si < lp->nsos; si++) {
	    set = lp->sos_set + si;
	    w = 1;
	    for (i = set->first; i < set->first + set->nelem; i++) {
		mpq_EGlpNumSet (lp->sos_weight[i], (double) w);
		w++;
	    }
	}
    }
    ILL_IFTRACE ("bound %lf <= x1 <= %lf\n", mpq_EGlpNumToLf (lp->lower[0]),
	mpq_EGlpNumToLf (lp->upper[0]));
    ILL_IFTRACE ("MAXDOUBLE %lf MINDOUBLE %lf\n", mpq_EGlpNumToLf (mpq_ILL_MAXDOUBLE),
	mpq_EGlpNumToLf (mpq_ILL_MINDOUBLE));
    ILL_RESULT (rval, "mpq_mps_fill_in");
}

/****************************************************************************/
/* writing */

static int mpq_mps_write_col (int i,
      int iorig,
      char *colname,
      mpq_ILLlpdata * lp,
      char **rownames,
      int intmode,
      char *objname);

int mpq_ILLwrite_mps (mpq_ILLlpdata * lp,
      mpq_qserror_collector * collector)
{
    int rval = 0;
    int marker = 0;
    int intmode = 0;
    int i, ri, set, el, empty, prtLower, prtUpper;
    char **colnames = (char **) NULL;
    char **rownames = (char **) NULL;
    mpq_ILLlp_rows lp_rows, *lprows = NULL;
    char buf[ILL_namebufsize];
    char *objname = NULL;

    ILL_CHECKnull (lp, "lp must not be null");
    ILL_FAILtrue (lp->probname == NULL, "oops should never happen");

    ILL_FAILfalse (lp->colnames != NULL, "colnames != NULL");
    ILL_FAILfalse (lp->rownames != NULL, "colnames != NULL");
    colnames = lp->colnames;
    rownames = lp->rownames;

    objname = lp->objname;
    if (objname == (char *) NULL) {
	strcpy (buf, "obj");
	rval = ILLsymboltab_uname (&lp->rowtab, buf, "", NULL);
	ILL_CLEANUP_IF (rval);
	mpq_ILL_UTIL_STR (objname, buf);
    }
    mpq_ILLprint_report (lp, "NAME    %s\n", lp->probname);
    mpq_ILLprint_report (lp, "OBJSENSE\n  %s\n",
	(lp->objsense == mpq_ILL_MIN) ? "MIN" : "MAX");
    mpq_ILLprint_report (lp, "OBJNAME\n  %s\n", objname);
    if (lp->refrowname) {
	mpq_ILLprint_report (lp, "REFROW\n");
	mpq_ILLprint_report (lp, " %s\n", lp->refrowname);
    }
    mpq_ILLprint_report (lp, "ROWS\n");
    mpq_ILLprint_report (lp, " N  %s\n", objname);
    if (lp->refrowname && (lp->refind == -1)) {
	mpq_ILLprint_report (lp, " N  %s\n", lp->refrowname);
    }
    lprows = &lp_rows;
    rval = mpq_ILLlp_rows_init (lprows, lp, 0);
    ILL_CLEANUP_IF (rval);
    for (i = 0; i < lp->nrows; i++) {
	if (lprows->rowcnt[i] == 0) {
	    mpq_ILLdata_warn (collector,
		"Deleting empty row \"%s\" from constraints.", rownames[i]);
	    continue;
	}
	switch (lp->sense[i]) {
	case 'G':
	    mpq_ILLprint_report (lp, " G  ");
	    break;
	case 'L':
	    mpq_ILLprint_report (lp, " L  ");
	    break;
	case 'E':
	    mpq_ILLprint_report (lp, " E  ");
	    break;
	case 'R':
	    mpq_ILLprint_report (lp, " G  ");
	    break;
	}
	mpq_ILLprint_report (lp, "%s\n", rownames[i]);
    }

    mpq_ILLprint_report (lp, "COLUMNS\n");
    for (set = 0; set < lp->sos.matcols; set++) {
	ILL_FAILtrue (lp->sos_type == NULL, "must have non NULL sos_type");
	ILL_FAILtrue (lp->is_sos_mem == NULL, "must have non NULL is_sos_mem");
	empty = 1;
	for (el = lp->sos.matbeg[set];
	    el < lp->sos.matbeg[set] + lp->sos.matcnt[set]; el++) {
	    if (empty) {
		mpq_ILLprint_report (lp, " %s SOS%dqs    'MARKER'    'SOSORG'\n",
		    ((lp->sos_type[set] == mpq_ILL_SOS_TYPE1)) ? "S1" : "S2",
		    marker++);
		empty = 0;
	    }
	    ri = lp->sos.matind[el];
	    i = lp->structmap[ri];
	    intmode = mpq_mps_write_col (i, ri, colnames[ri], lp, rownames,
		intmode, objname);
	    if (lp->refrowname && (lp->refind == -1)) {
		mpq_ILLprint_report (lp, "  %s    %s    %g\n",
		    colnames[ri], lp->refrowname, lp->sos.matval[el]);
	    }
	}
	if (!empty) {
	    mpq_ILLprint_report (lp, " SOS%dqs       'MARKER'    'SOSEND'\n", marker++);
	}
    }
    for (ri = 0; ri < lp->nstruct; ri++) {
	if (lp->is_sos_mem == (int *) NULL || lp->is_sos_mem[ri] == -1) {
	    i = lp->structmap[ri];
	    intmode = mpq_mps_write_col (i, ri, colnames[ri], lp, rownames,
		intmode, objname);
	}
    }
    if (intmode) {
	mpq_ILLprint_report (lp, " MARK%dqs      'MARKER'    'INTEND'\n", lp->nstruct);
    }
    mpq_ILLprint_report (lp, "RHS\n");
    for (i = 0; i < lp->nrows; i++) {
	if ((lprows->rowcnt[i] != 0) && mpq_EGlpNumIsNeqqZero (lp->rhs[i])) {
	    mpq_ILLprint_report (lp, " RHS    %s    %g\n", rownames[i],
		mpq_EGlpNumToLf (lp->rhs[i]));
	}
    }

    if (lp->rangeval) {
	mpq_ILLprint_report (lp, "RANGES\n");
	for (i = 0; i < lp->nrows; i++) {
	    if ((lprows->rowcnt[i] != 0) && mpq_EGlpNumIsNeqqZero (lp->rangeval[i])) {
		mpq_ILLprint_report (lp, " RANGE    %s    %g\n", rownames[i],
		    mpq_EGlpNumToLf (lp->rangeval[i]));
	    }
	}
    }
    ri = mpq_ILLraw_first_nondefault_bound (lp);
    if (ri != lp->nstruct) {
	mpq_ILLprint_report (lp, "BOUNDS\n");
	for (ri = ri; ri < lp->nstruct; ri++) {
	    i = lp->structmap[ri];
	    if (mpq_EGlpNumIsEqqual (lp->lower[i], lp->upper[i])) {
		mpq_ILLprint_report (lp, " FX BOUND    %s    %g\n", colnames[ri],
		    mpq_EGlpNumToLf (lp->lower[i]));
		continue;
	    }
	    if ((mpq_EGlpNumIsEqqual (lp->lower[i], mpq_ILL_MINDOUBLE)) &&
		(mpq_EGlpNumIsEqqual (lp->upper[i], mpq_ILL_MAXDOUBLE))) {
		mpq_ILLprint_report (lp, " FR BOUND    %s\n", colnames[ri]);
		continue;
	    }
	    prtLower = !mpq_ILLraw_default_lower (lp, i);
	    prtUpper = !mpq_ILLraw_default_upper (lp, i);
	    if (prtLower) {
		if (mpq_EGlpNumIsEqqual (lp->lower[i], mpq_ILL_MINDOUBLE)) {
		    mpq_ILLprint_report (lp, " MI BOUND    %s\n", colnames[ri]);
		} else {
		    mpq_ILLprint_report (lp, " LO BOUND    %s    %g\n", colnames[ri],
			mpq_EGlpNumToLf (lp->lower[i]));
		}
	    }
	    if (prtUpper) {
		if (mpq_EGlpNumIsEqqual (lp->upper[i], mpq_ILL_MAXDOUBLE)) {
		    mpq_ILLprint_report (lp, " PL BOUND    %s\n", colnames[ri]);
		} else {
		    mpq_ILLprint_report (lp, " UP BOUND    %s    %g\n", colnames[ri],
			mpq_EGlpNumToLf (lp->upper[i]));
		}
	    }
	}
    }
    mpq_ILLprint_report (lp, "ENDATA\n");

CLEANUP:
    mpq_ILLlp_rows_clear (lprows);
    if (!lp->colnames) {
	ILLfree_names (colnames, lp->ncols);
    }
    if (!lp->rownames) {
	ILLfree_names (rownames, lp->nrows);
    }
    if (objname != lp->objname) {
	ILL_IFFREE (objname, char);
    }
    ILL_RESULT (rval, "ILLlpdata_mpswrite");
}

static int mpq_mps_write_col (int i,
      int iorig,
      char *colname,
      mpq_ILLlpdata * lp,
      char **rownames,
      int intmode,
      char *objname)
{
    int row, k;
    mpq_ILLmatrix *A;

    A = &lp->A;
    if (lp->intmarker && (lp->intmarker[iorig] != intmode)) {
	mpq_ILLprint_report (lp, " MARK%dqs      'MARKER'    '%s'\n", iorig,
	    (lp->intmarker[iorig] ? "INTORG" : "INTEND"));
	intmode = lp->intmarker[iorig];
    }
    if (mpq_EGlpNumIsNeqqZero (lp->obj[i])) {
	mpq_ILLprint_report (lp, "  %s    %s    %g\n", colname, objname,
	    mpq_EGlpNumToLf (lp->obj[i]));
    }
    for (k = A->matbeg[i]; k < A->matbeg[i] + A->matcnt[i]; k++) {
	row = A->matind[k];
	mpq_ILLprint_report (lp, "  %s    %s    %g\n", colname, rownames[row],
	    mpq_EGlpNumToLf (A->matval[k]));
    }
    return intmode;
}
