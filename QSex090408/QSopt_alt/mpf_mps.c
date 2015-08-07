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

/* RCS_INFO = "$RCSfile: mpf_mps.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

/****************************************************************************/
/* */
/* Routines for Reading and Writing MPS Files                 */
/* */
/* EXPORTED FUNCTIONS                                                      */
/* */
/* int  ILLlpdata_mpsread (mpf_ILLlpdata *lp, const char *filename);         */
/* int  ILLlpdata_mpswrite(mpf_ILLlpdata *lp, const char *filename);         */
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
#include "mpf_iqsutil.h"
#include "mpf_mps.h"
#include "mpf_rawlp.h"
#include "mpf_lpdata.h"
/* extern mpf_t mpf_SZERO_TOLER; */

static int TRACE = 0;

const char *mpf_ILLmps_section_name[ILL_MPS_N_SECTIONS + 2] = {
    "NAME", "OBJSENSE", "OBJNAME", "ROWS", "COLUMNS",
    "RHS", "RANGES", "BOUNDS", "REFROW", "ENDATA",
    NULL
};

static const char *mpf_mps_bound_name[] = {
    "LO", "UP", "FX", "FR", "MI", "PL", "BV", "UI", "LI", NULL
};

/****************************************************************************/
/* reading */
static int mpf_read_mps_section (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);

static int mpf_read_mps_name (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_read_mps_refrow (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_read_mps_objnamesense (ILLmps_section sec,
      mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_read_mps_objname (mpf_ILLread_mps_state * state);
static int mpf_read_mps_objsense (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);

static int mpf_read_mps_line_in_section (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);


static int mpf_add_row (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_add_col (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_add_rhs (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_add_ranges (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_add_bounds (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);

static int mpf_mps_read_marker_line (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);
static int mpf_is_marker_line (mpf_ILLread_mps_state * state);
static int mpf_mps_read_col_line (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp);

static int mpf_mps_fill_in (mpf_rawlpdata * lp,
      const char *obj);


static void mpf_mps_set_bound (mpf_rawlpdata * lp,
      mpf_ILLread_mps_state * state,
      int colind,
      const char *bndtype,
      mpf_t bnd);

int mpf_ILLread_mps (mpf_qsline_reader * file,
      const char *f,
      mpf_rawlpdata * lp)
{
    int rval = 0;
    int end = 0;
    mpf_ILLread_mps_state state;

    ILL_IFTRACE ("\tread_mps\n");
    if (ILLsymboltab_create (&lp->rowtab, 100) ||
	ILLsymboltab_create (&lp->coltab, 100)) {
	rval = 1;
    } else {
	rval = mpf_ILLmps_state_init (&state, file, f);
    }
    if (rval != 0) {
	goto CLEANUP;
    }
    while (mpf_ILLmps_next_line (&state) == 0) {
	if (mpf_ILLmps_empty_key (&state)) {
	    if (mpf_read_mps_line_in_section (&state, lp) != 0) {
		rval++;
	    }
	} else {
	    /* found a section indicator in col 1 */
	    if (!strcmp (state.key, mpf_ILLmps_section_name[ILL_MPS_ENDATA])) {
		end = 1;
		break;		/* done reading */
	    }
	    if (mpf_read_mps_section (&state, lp) != 0) {
		rval++;
	    }
	}
	if (rval == 50) {
	    (void) mpf_ILLmps_error (&state, "Too many errors.\n");
	}
    }

    if (!end) {
	mpf_ILLmps_warn (&state, "Missing ENDATA.");
    }
    if (!mpf_ILLmps_next_line (&state)) {
	mpf_ILLmps_warn (&state, "Ignoring text after ENDATA.");
    }
    if (rval == 0) {
	rval = mpf_mps_fill_in (lp, state.obj);
    }
CLEANUP:
    ILL_RESULT (rval, "read_mps");
}

static int mpf_check_section_order (mpf_ILLread_mps_state * state,
      int sec)
{
    switch (sec) {
    case ILL_MPS_REFROW:
	if (state->section[ILL_MPS_ROWS] == 1) {
	    return mpf_ILLmps_error (state, "%s section after ROWS section.\n",
		mpf_ILLmps_section_name[sec]);
	}
	break;

    case ILL_MPS_COLS:
    case ILL_MPS_RHS:
    case ILL_MPS_RANGES:
	if (state->section[ILL_MPS_ROWS] == 0) {
	    return mpf_ILLmps_error (state, "%s section before ROWS section.\n",
		mpf_ILLmps_section_name[sec]);
	};
	break;

    case ILL_MPS_BOUNDS:
	if (state->section[ILL_MPS_COLS] == 0) {
	    return mpf_ILLmps_error (state, "%s section before COLUMNS section.\n",
		mpf_ILLmps_section_name[sec]);
	}
	break;
    }
    return 0;
}

static int mpf_read_mps_section (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int sec;
    int rval = 0, r;

    ILL_FAILtrue (mpf_ILLmps_empty_key (state), "must have a key on this line");

    sec = mpf_ILLutil_index (mpf_ILLmps_section_name, state->key);
    if (sec < 0) {
	return mpf_ILLmps_error (state, "\"%s\" is not a key.\n", state->key);
    }
    rval = mpf_ILLmps_set_section (state, sec);

    state->active = ILL_MPS_NONE;
    rval = rval || mpf_check_section_order (state, sec);
    switch (sec) {
    case ILL_MPS_COLS:
    case ILL_MPS_ROWS:
	state->active = sec;
	break;

    case ILL_MPS_NAME:
	if (rval == 0) {
	    rval = mpf_read_mps_name (state, lp);
	}
	break;

    case ILL_MPS_RHS:
	if (rval == 0) {
	    rval = mpf_ILLraw_init_rhs (lp);
	}
	state->active = ILL_MPS_RHS;
	break;

    case ILL_MPS_RANGES:
	if (rval == 0) {
	    rval = mpf_ILLraw_init_ranges (lp);
	}
	state->active = ILL_MPS_RANGES;
	break;

    case ILL_MPS_BOUNDS:
	if (rval == 0) {
	    rval = mpf_ILLraw_init_bounds (lp);
	}
	state->active = ILL_MPS_BOUNDS;
	break;

    case ILL_MPS_OBJNAME:
    case ILL_MPS_OBJSENSE:
	r = mpf_read_mps_objnamesense (sec, state, lp);
	rval = r || rval;
	break;

    case ILL_MPS_REFROW:
	r = mpf_read_mps_refrow (state, lp);
	rval = r || rval;
	break;

    default:
	ILL_REPRT ("should never get here");
	goto CLEANUP;
    }
CLEANUP:
    ILL_RESULT (rval, "mpf_read_mps_section");
}

static int mpf_read_mps_name (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rval = 0;

    if (mpf_ILLmps_empty_field (state)) {
	mpf_ILLmps_warn (state, "Blank NAME.");
    } else {
	mpf_ILL_UTIL_STR (lp->name, state->field);
    }
CLEANUP:
    ILL_RESULT (rval, "mpf_read_mps_name");
}

static int mpf_read_mps_refrow (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rval = 0;

    rval = mpf_ILLmps_next_line (state);
    if (state->section[ILL_MPS_REFROW] > 1) {
	/* this is the second time we see this section; don't complain about
	   errors */
	return 0;
    }
    if (mpf_ILLmps_empty_key (state) && !mpf_ILLmps_empty_field (state)) {
	mpf_ILL_UTIL_STR (lp->refrow, state->field);
	return 0;
    } else {
	return mpf_ILLmps_error (state, "Bad row name in REFROW section.\n");
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_read_mps_refrow");
}

static int mpf_read_mps_objnamesense (ILLmps_section sec,
      mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    if (state->section[sec] > 1) {
	/* this is the second time we see this section; just skip next line */
	(void) mpf_ILLmps_next_line (state);
	return 0;
    }
    if (mpf_ILLmps_next_line (state) != 0) {
	return mpf_ILLmps_error (state, "Missing %s line at end of file.\n",
	    mpf_ILLmps_section_name[sec]);
    }
    if ((!mpf_ILLmps_empty_key (state)) || mpf_ILLmps_empty_field (state)) {
	(void) mpf_ILLmps_error (state, "Bad %s in %s record.\n",
	    ((sec == ILL_MPS_OBJNAME) ? "row name"
		: "objective sense"), mpf_ILLmps_section_name[sec]);
	if (!mpf_ILLmps_empty_key (state)) {
	    (void) mpf_read_mps_section (state, lp);
	}
	return 1;
    }
    if (sec == ILL_MPS_OBJNAME) {
	if (mpf_read_mps_objname (state)) {
	    return 1;
	}
    } else {
	if (mpf_read_mps_objsense (state, lp) != 0) {
	    return 1;
	}
    }
    return 0;
}

static int mpf_read_mps_objsense (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rval = 0;
    char *objsense = state->field;
    ILL_FAILfalse (state->section[ILL_MPS_OBJSENSE] == 1, "should never happen");
    if (!strcmp (objsense, "MAX") ||
	!strcmp (objsense, "Max") ||
	!strcmp (objsense, "max") ||
	!strcmp (objsense, "MAXIMIZE") ||
	!strcmp (objsense, "Maximize") || !strcmp (objsense, "maximize")) {
	lp->objsense = mpf_ILL_MAX;
    } else if (!strcmp (objsense, "MIN") ||
	    !strcmp (objsense, "Min") ||
	    !strcmp (objsense, "min") ||
	    !strcmp (objsense, "MINIMIZE") ||
	!strcmp (objsense, "Minimize") || !strcmp (objsense, "minimize")) {
	lp->objsense = mpf_ILL_MIN;
    } else {
	return mpf_ILLmps_error (state, "\"%s\" is no OBJSENSE.\n", objsense);
    }
CLEANUP:
    ILL_RESULT (rval, "mpf_read_mps_objsense");
}

static int mpf_read_mps_objname (mpf_ILLread_mps_state * state)
{
    int rval = 0;
    ILL_FAILfalse (state->section[ILL_MPS_OBJNAME] == 1, "should never happen");
    mpf_ILL_UTIL_STR (state->obj, state->field);
CLEANUP:
    ILL_RETURN (rval, "mpf_read_mps_objname");
}

static int mpf_read_mps_line_in_section (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rval = 0;

    ILL_FAILtrue (!mpf_ILLmps_empty_key (state) || mpf_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    if (state->active == ILL_MPS_NONE) {
	return mpf_ILLmps_error (state, "Line is in no section.\n");
    } else {
	if (state->section[state->active] == 1) {
	    switch (state->active) {
	    case ILL_MPS_ROWS:
		rval = mpf_add_row (state, lp);
		break;
	    case ILL_MPS_COLS:
		rval = mpf_add_col (state, lp);
		break;
	    case ILL_MPS_RHS:
		rval = mpf_add_rhs (state, lp);
		break;
	    case ILL_MPS_RANGES:
		rval = mpf_add_ranges (state, lp);
		break;
	    case ILL_MPS_BOUNDS:
		rval = mpf_add_bounds (state, lp);
		break;
	    default:
		ILL_REPRT ("should never get here");
		ILL_CLEANUP;
	    }
	    if (rval == 0) {
		/* see whether there are extra fields on line */
		mpf_ILLmps_check_end_of_line (state);
	    }
	}
    }
CLEANUP:
    ILL_RESULT (rval, "mpf_read_mps_line_in_section");
}

static int mpf_add_row (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int ind, hit, rval = 0;
    char sense;

    ILL_FAILtrue (!mpf_ILLmps_empty_key (state) || mpf_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    /* field should contain exactly one character */
    if (state->field[1] == '\0') {
	sense = state->field[0];
	if (sense != 'L' && sense != 'G' && sense != 'E' && sense != 'N') {
	    return mpf_ILLmps_error (state,
		"Unknown rowsense '%c' in ROWS record.\n", sense);
	}
	if (mpf_ILLmps_next_field (state) == 0) {
	    hit = ILLsymboltab_lookup (&lp->rowtab, state->field, &ind);
	    if (!hit) {
		rval = mpf_ILLmps_error (state,
		    "Repeated row definition for \"%s\".\n",
		    state->field);
	    } else {
		rval = mpf_ILLraw_add_row (lp, state->field, sense, mpf_zeroLpNum);	/* default rhs is 0.0 */
	    }
	} else {
	    rval = mpf_ILLmps_error (state, "Missing rowname in ROWS record.\n");
	}
    } else {
	rval = mpf_ILLmps_error (state, "Unknown rowsense '%s' in ROWS record.\n",
	    state->field);
    }
CLEANUP:
    ILL_RESULT (rval, "mpf_add_row");
}

static int mpf_add_col (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rval = 0;

    ILL_FAILtrue (!mpf_ILLmps_empty_key (state) || mpf_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    if (mpf_is_marker_line (state)) {
	return mpf_mps_read_marker_line (state, lp);
    } else {
	return mpf_mps_read_col_line (state, lp);
    }
CLEANUP:
    ILL_RETURN (rval, "mpf_add_col");
}

static int mpf_mps_read_col_line (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int hit, colind, rowind, rval = 0, more, ind;
    mpf_t ncoef;
    mpf_EGlpNumInitVar (ncoef);

    ILL_FAILtrue (!mpf_ILLmps_empty_key (state) || mpf_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    hit = ILLsymboltab_lookup (&lp->coltab, state->field, &colind);
    if (hit) {
	rval = mpf_ILLraw_add_col (lp, state->field, state->intvar);
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
	if (mpf_ILLraw_is_mem_other_sos (lp, colind)) {
	    rval = mpf_ILLmps_error (state,
		"\"%s\" is a member of SOS set #%d.\n",
		mpf_ILLraw_colname (lp, colind),
		lp->is_sos_member[colind] + 1);
	} else {
	    rval = mpf_ILLraw_add_sos_member (lp, colind);
	}
	ILL_CLEANUP_IF (rval);
    }
    more = (mpf_ILLmps_next_field (state) == 0);
    if (!more) {
	return mpf_ILLmps_error (state, "Missing fields in COLUMNS record.\n");
    }
    for (more = 1; more; more = (mpf_ILLmps_next_field (state) == 0)) {
	hit = ILLsymboltab_lookup (&lp->rowtab, state->field, &rowind);
	if (hit) {
	    return mpf_ILLmps_error (state, "\"%s\" is not a row name.\n", state->field);
	}
	if (mpf_ILLmps_next_coef (state, &ncoef) != 0) {
	    return mpf_ILLmps_error (state,
		"Missing/Bad coefficient in COLUMNS record.\n");
	}
	rval = mpf_ILLraw_add_col_coef (lp, colind, rowind, ncoef);
    }
CLEANUP:
    mpf_EGlpNumClearVar (ncoef);
    ILL_RESULT (rval, "mpf_mps_read_col_line");
}

static int mpf_is_marker_line (mpf_ILLread_mps_state * state)
{
    const char *field = state->line;
    while ((field = mpf_ILLutil_strchr (field, '\''))) {
	if (strncmp (field, "'MARKER'", (size_t) 8) == 0) {
	    return 1;
	}
	while ((!mpf_ILL_ISBLANK (field)) && (*field != '\0')) {
	    field++;
	}
    }
    return 0;
}

static int mpf_mps_read_marker_line (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rval = 0;
    int sos_type = mpf_ILL_SOS_TYPE1;
    int sos_line = 0;
    int cur_sos_mode = state->sosvar;

    if (strcmp (state->field, "S2") == 0) {
	sos_type = mpf_ILL_SOS_TYPE2;
	sos_line = 1;
    } else if (strcmp (state->field, "S1") == 0) {
	sos_line = 1;
    }
    if (sos_line) {
	rval = mpf_ILLmps_next_field (state);
    }
    rval = rval || mpf_ILLmps_next_field (state);	/* swallow marker-id */
    if (strcmp (state->field, "'MARKER'")) {
	return mpf_ILLmps_error (state, "Bad 'MARKER' line.\n");
    }
    if (mpf_ILLmps_next_field (state)) {
	return mpf_ILLmps_error (state, "Missing field on 'MARKER' line.\n");
    }
    rval = mpf_ILLmps_int_sos_mode (state);
    if (rval == 0) {
	if (cur_sos_mode != state->sosvar) {
	    if (state->sosvar) {
		/* beginning of a new set */
		rval = mpf_ILLraw_add_sos (lp, sos_type);
	    }
	}
    }
    ILL_RESULT (rval, "mpf_mps_read_marker_line");
}

static int mpf_add_rhs (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    int rowind, more_fields, skip;
    const char *rhsname;
    mpf_t rhs;
    mpf_EGlpNumInitVar (rhs);

    rhsname = mpf_ILLmps_possibly_blank_name (state->field, state, &lp->rowtab);
    if (mpf_ILLraw_set_rhs_name (lp, rhsname, &skip)) {
	mpf_ILLmps_error (state, "Could not add right hand side.\n");
    }
    if (skip) {
	mpf_ILLmps_set_end_of_line (state);	/* to avoid warning about
						   extra fields */
    } else {
	if (strcmp (rhsname, " ")) {
	    /* field is non blank rhs name; advance to row name  */
	    if (mpf_ILLmps_next_field (state)) {
		return mpf_ILLmps_error (state, "Missing row name in RHS record.\n");
	    }
	}
	for (more_fields = 1; more_fields; more_fields = !mpf_ILLmps_next_field (state)) {
	    if (ILLsymboltab_lookup (&lp->rowtab, state->field, &rowind)) {
		return mpf_ILLmps_error (state, "\"%s\" is not a row name.\n",
		    state->field);
	    }
	    if (mpf_ILLmps_next_coef (state, &rhs)) {
		return mpf_ILLmps_error (state, "Missing/Bad coefficient in RHS record.\n");
	    }
	    if (lp->rhsind[rowind]) {
		return mpf_ILLmps_error (state, "Two rhs values for row \"%s\".\n",
		    state->field);
	    } else {
		if (lp->rowsense[rowind] == 'N') {
		    mpf_ILLmps_warn (state,
			"Ignoring right hand side for N-row \"%s\".",
			mpf_ILLraw_rowname (lp, rowind));
		} else {
		    mpf_EGlpNumCopy (lp->rhs[rowind], rhs);
		    lp->rhsind[rowind] = 1;
		}
	    }
	}
    }
    mpf_EGlpNumClearVar (rhs);
    return 0;
}

static int mpf_add_bounds (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    char bndtype[3];
    const char *bounds_name;
    int colind, skip, rval = 0;
    mpf_t bnd;
    mpf_EGlpNumInitVar (bnd);

    ILL_FAILtrue (!mpf_ILLmps_empty_key (state) || mpf_ILLmps_empty_field (state),
	"no key but at least one field on state->line");

    if (mpf_ILLutil_index (mpf_mps_bound_name, state->field) < 0) {
	return mpf_ILLmps_error (state, "\"%s\" is not a BOUNDS type.\n", state->field);
    }
    strcpy (bndtype, state->field);

    if (mpf_ILLmps_next_field (state) != 0) {
	return mpf_ILLmps_error (state,
	    "No bounds/column identifier in BOUNDS record.\n");
    }
    bounds_name = mpf_ILLmps_possibly_blank_name (state->field, state, &lp->coltab);
    if (bounds_name == NULL) {
	return 1;
    }
    if (mpf_ILLraw_set_bounds_name (lp, bounds_name, &skip)) {
	return 1;
    }
    if (skip) {
	mpf_ILLmps_set_end_of_line (state);	/* to avoid warning about
						   extra fields */
    } else {
	if (strcmp (bounds_name, " ")) {
	    /* non empty bounds_name ==> advance to col name field */
	    if (mpf_ILLmps_next_field (state)) {
		return mpf_ILLmps_error (state, "Missing column field in BOUNDS record.\n");
	    }
	}
	if (ILLsymboltab_lookup (&lp->coltab, state->field, &colind)) {
	    return mpf_ILLmps_error (state, "\"%s\" is not a column name.\n",
		state->field);
	}
	mpf_EGlpNumZero (bnd);
	if (strcmp (bndtype, "FR") && strcmp (bndtype, "BV") &&
	    strcmp (bndtype, "MI") && strcmp (bndtype, "PL")) {
	    /* neither "FR", "BV", "MI" nor "PL" ==> there should be a bound */
	    if (mpf_ILLmps_next_bound (state, &bnd)) {
		return mpf_ILLmps_error (state,
		    "Missing/Bad bound field in BOUNDS record.\n");
	    }
	}
	mpf_mps_set_bound (lp, state, colind, bndtype, bnd);
    }
CLEANUP:
    mpf_EGlpNumClearVar (bnd);
    ILL_RESULT (rval, "mpf_add_bounds");
}

static void mpf_mps_set_bound (mpf_rawlpdata * lp,
      mpf_ILLread_mps_state * state,
      int colind,
      const char *bndtype,
      mpf_t bnd)
{
    const char *msg = NULL;
    if (!strcmp (bndtype, "LO")) {
	msg = mpf_ILLraw_set_lowerBound (lp, colind, bnd);
    } else if (!strcmp (bndtype, "UP")) {
	msg = mpf_ILLraw_set_upperBound (lp, colind, bnd);
    } else if (!strcmp (bndtype, "FX")) {
	msg = mpf_ILLraw_set_fixedBound (lp, colind, bnd);
    } else if (!strcmp (bndtype, "FR")) {
	msg = mpf_ILLraw_set_unbound (lp, colind);
    } else if (!strcmp (bndtype, "BV")) {
	msg = mpf_ILLraw_set_binaryBound (lp, colind);
	if (msg == NULL) {
	    lp->intmarker[colind] = 1;
	}
    } else if (!strcmp (bndtype, "UI")) {
	msg = mpf_ILLraw_set_upperBound (lp, colind, bnd);
	if (msg == NULL) {
	    lp->intmarker[colind] = 1;
	}
    } else if (!strcmp (bndtype, "LI")) {
	msg = mpf_ILLraw_set_lowerBound (lp, colind, bnd);
	if (msg == NULL) {
	    lp->intmarker[colind] = 1;
	}
    } else if (!strcmp (bndtype, "MI")) {
	msg = mpf_ILLraw_set_lowerBound (lp, colind, mpf_ILL_MINDOUBLE);
    } else if (!strcmp (bndtype, "PL")) {
	msg = mpf_ILLraw_set_upperBound (lp, colind, mpf_ILL_MAXDOUBLE);
    } else {
	ILL_REPRT ("should never get here");
	ILL_CLEANUP;
    }
    mpf_ILLmps_warn (state, msg);
CLEANUP:
    return;
}

static int mpf_add_ranges (mpf_ILLread_mps_state * state,
      mpf_rawlpdata * lp)
{
    const char *rangesname;
    int skip, more_fields, rowind;
    mpf_t ntmp;
    mpf_EGlpNumInitVar (ntmp);

    rangesname = mpf_ILLmps_possibly_blank_name (state->field, state, &lp->rowtab);
    if (mpf_ILLraw_set_ranges_name (lp, rangesname, &skip)) {
	return mpf_ILLmps_error (state, "Could not add range.\n");
    }
    if (skip) {
	mpf_ILLmps_set_end_of_line (state);	/* to avoid warning about
						   extra fields */
    } else {
	if (strcmp (rangesname, " ")) {
	    /* field is non blank ranges name; advance to row name */
	    if (mpf_ILLmps_next_field (state)) {
		return mpf_ILLmps_error (state, "Missing row name in RANGES record.");
	    }
	}
	for (more_fields = 1; more_fields; more_fields = !mpf_ILLmps_next_field (state)) {
	    if (ILLsymboltab_lookup (&lp->rowtab, state->field, &rowind)) {
		return mpf_ILLmps_error (state, "\"%s\" is not a row name.\n",
		    state->field);
	    }
	    if (mpf_ILLmps_next_coef (state, &ntmp)) {
		return mpf_ILLmps_error (state,
		    "Missing/Bad coefficient in RANGES record.\n");
	    }
	    if (lp->rangesind[rowind]) {
		mpf_ILLmps_warn (state, "Ignoring second RANGE value %s \"%s\".",
		    "for row", mpf_ILLraw_rowname (lp, rowind));
	    } else {
		if (lp->rowsense[rowind] != 'N') {
		    if (mpf_ILLraw_add_ranges_coef (lp, rowind, ntmp))
			return 1;
		} else {
		    mpf_ILLmps_warn (state, "Ignoring RANGE value for N-row \"%s\".",
			mpf_ILLraw_rowname (lp, rowind));
		}
	    }
	}
    }
    mpf_EGlpNumClearVar (ntmp);
    return 0;
}


static int mpf_mps_fill_in (mpf_rawlpdata * lp,
      const char *obj)
{
    int i, hit, rval = 0;

    /* find the objective function -- the first N row if obj is not defined  */
    if (obj) {
	hit = ILLsymboltab_lookup (&lp->rowtab, obj, &lp->objindex);
	if (hit) {
	    return mpf_ILLdata_error (lp->error_collector,
		"Bad objective name \"%s\".\n", obj);
	} else if (lp->rowsense[lp->objindex] != 'N') {
	    mpf_ILLdata_warn (lp->error_collector,
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
	    return mpf_ILLdata_error (lp->error_collector,
		"No N-row in lp definition.\n");
	}
    }

    if (lp->ncols > 0) {
	rval = mpf_ILLraw_fill_in_bounds (lp);
    }
    /* set weights of sos set members */
    if (lp->refrow) {
	/* take the values from refrow */
	mpf_t weight;
	mpf_colptr *cp;
	mpf_EGlpNumInitVar (weight);

	hit = ILLsymboltab_lookup (&lp->rowtab, lp->refrow, &lp->refrowind);
	if (hit) {
	    return mpf_ILLdata_error (lp->error_collector,
		"REFROW \"%s\" is not a row name.\n", lp->refrow);
	}
	for (i = 0; i < lp->nsos_member; i++) {
	    for (cp = lp->cols[lp->sos_col[i]]; cp != NULL; cp = cp->next) {
		if (cp->this == lp->refrowind)
		    break;
	    }
	    if ((cp != NULL) && mpf_EGlpNumIsNeqqZero (cp->coef)) {
		mpf_EGlpNumCopy (weight, cp->coef);
	    } else {
		mpf_ILLdata_warn (lp->error_collector,
		    "\"%s\" has 0.0 coefficient in REFROW \"%s\".",
		    mpf_ILLraw_colname (lp, lp->sos_col[i]), lp->refrow);
		mpf_EGlpNumZero (weight);
	    }
	    mpf_EGlpNumCopy (lp->sos_weight[i], weight);
	}
	mpf_EGlpNumClearVar (weight);
    } else {			/* no refrow */
	/* set weights to 1, 2, 3, ... in order of definition */
	int si, w;
	mpf_sosptr *set;
	for (si = 0; si < lp->nsos; si++) {
	    set = lp->sos_set + si;
	    w = 1;
	    for (i = set->first; i < set->first + set->nelem; i++) {
		mpf_EGlpNumSet (lp->sos_weight[i], (double) w);
		w++;
	    }
	}
    }
    ILL_IFTRACE ("bound %lf <= x1 <= %lf\n", mpf_EGlpNumToLf (lp->lower[0]),
	mpf_EGlpNumToLf (lp->upper[0]));
    ILL_IFTRACE ("MAXDOUBLE %lf MINDOUBLE %lf\n", mpf_EGlpNumToLf (mpf_ILL_MAXDOUBLE),
	mpf_EGlpNumToLf (mpf_ILL_MINDOUBLE));
    ILL_RESULT (rval, "mpf_mps_fill_in");
}

/****************************************************************************/
/* writing */

static int mpf_mps_write_col (int i,
      int iorig,
      char *colname,
      mpf_ILLlpdata * lp,
      char **rownames,
      int intmode,
      char *objname);

int mpf_ILLwrite_mps (mpf_ILLlpdata * lp,
      mpf_qserror_collector * collector)
{
    int rval = 0;
    int marker = 0;
    int intmode = 0;
    int i, ri, set, el, empty, prtLower, prtUpper;
    char **colnames = (char **) NULL;
    char **rownames = (char **) NULL;
    mpf_ILLlp_rows lp_rows, *lprows = NULL;
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
	mpf_ILL_UTIL_STR (objname, buf);
    }
    mpf_ILLprint_report (lp, "NAME    %s\n", lp->probname);
    mpf_ILLprint_report (lp, "OBJSENSE\n  %s\n",
	(lp->objsense == mpf_ILL_MIN) ? "MIN" : "MAX");
    mpf_ILLprint_report (lp, "OBJNAME\n  %s\n", objname);
    if (lp->refrowname) {
	mpf_ILLprint_report (lp, "REFROW\n");
	mpf_ILLprint_report (lp, " %s\n", lp->refrowname);
    }
    mpf_ILLprint_report (lp, "ROWS\n");
    mpf_ILLprint_report (lp, " N  %s\n", objname);
    if (lp->refrowname && (lp->refind == -1)) {
	mpf_ILLprint_report (lp, " N  %s\n", lp->refrowname);
    }
    lprows = &lp_rows;
    rval = mpf_ILLlp_rows_init (lprows, lp, 0);
    ILL_CLEANUP_IF (rval);
    for (i = 0; i < lp->nrows; i++) {
	if (lprows->rowcnt[i] == 0) {
	    mpf_ILLdata_warn (collector,
		"Deleting empty row \"%s\" from constraints.", rownames[i]);
	    continue;
	}
	switch (lp->sense[i]) {
	case 'G':
	    mpf_ILLprint_report (lp, " G  ");
	    break;
	case 'L':
	    mpf_ILLprint_report (lp, " L  ");
	    break;
	case 'E':
	    mpf_ILLprint_report (lp, " E  ");
	    break;
	case 'R':
	    mpf_ILLprint_report (lp, " G  ");
	    break;
	}
	mpf_ILLprint_report (lp, "%s\n", rownames[i]);
    }

    mpf_ILLprint_report (lp, "COLUMNS\n");
    for (set = 0; set < lp->sos.matcols; set++) {
	ILL_FAILtrue (lp->sos_type == NULL, "must have non NULL sos_type");
	ILL_FAILtrue (lp->is_sos_mem == NULL, "must have non NULL is_sos_mem");
	empty = 1;
	for (el = lp->sos.matbeg[set];
	    el < lp->sos.matbeg[set] + lp->sos.matcnt[set]; el++) {
	    if (empty) {
		mpf_ILLprint_report (lp, " %s SOS%dqs    'MARKER'    'SOSORG'\n",
		    ((lp->sos_type[set] == mpf_ILL_SOS_TYPE1)) ? "S1" : "S2",
		    marker++);
		empty = 0;
	    }
	    ri = lp->sos.matind[el];
	    i = lp->structmap[ri];
	    intmode = mpf_mps_write_col (i, ri, colnames[ri], lp, rownames,
		intmode, objname);
	    if (lp->refrowname && (lp->refind == -1)) {
		mpf_ILLprint_report (lp, "  %s    %s    %g\n",
		    colnames[ri], lp->refrowname, lp->sos.matval[el]);
	    }
	}
	if (!empty) {
	    mpf_ILLprint_report (lp, " SOS%dqs       'MARKER'    'SOSEND'\n", marker++);
	}
    }
    for (ri = 0; ri < lp->nstruct; ri++) {
	if (lp->is_sos_mem == (int *) NULL || lp->is_sos_mem[ri] == -1) {
	    i = lp->structmap[ri];
	    intmode = mpf_mps_write_col (i, ri, colnames[ri], lp, rownames,
		intmode, objname);
	}
    }
    if (intmode) {
	mpf_ILLprint_report (lp, " MARK%dqs      'MARKER'    'INTEND'\n", lp->nstruct);
    }
    mpf_ILLprint_report (lp, "RHS\n");
    for (i = 0; i < lp->nrows; i++) {
	if ((lprows->rowcnt[i] != 0) && mpf_EGlpNumIsNeqqZero (lp->rhs[i])) {
	    mpf_ILLprint_report (lp, " RHS    %s    %g\n", rownames[i],
		mpf_EGlpNumToLf (lp->rhs[i]));
	}
    }

    if (lp->rangeval) {
	mpf_ILLprint_report (lp, "RANGES\n");
	for (i = 0; i < lp->nrows; i++) {
	    if ((lprows->rowcnt[i] != 0) && mpf_EGlpNumIsNeqqZero (lp->rangeval[i])) {
		mpf_ILLprint_report (lp, " RANGE    %s    %g\n", rownames[i],
		    mpf_EGlpNumToLf (lp->rangeval[i]));
	    }
	}
    }
    ri = mpf_ILLraw_first_nondefault_bound (lp);
    if (ri != lp->nstruct) {
	mpf_ILLprint_report (lp, "BOUNDS\n");
	for (ri = ri; ri < lp->nstruct; ri++) {
	    i = lp->structmap[ri];
	    if (mpf_EGlpNumIsEqqual (lp->lower[i], lp->upper[i])) {
		mpf_ILLprint_report (lp, " FX BOUND    %s    %g\n", colnames[ri],
		    mpf_EGlpNumToLf (lp->lower[i]));
		continue;
	    }
	    if ((mpf_EGlpNumIsEqqual (lp->lower[i], mpf_ILL_MINDOUBLE)) &&
		(mpf_EGlpNumIsEqqual (lp->upper[i], mpf_ILL_MAXDOUBLE))) {
		mpf_ILLprint_report (lp, " FR BOUND    %s\n", colnames[ri]);
		continue;
	    }
	    prtLower = !mpf_ILLraw_default_lower (lp, i);
	    prtUpper = !mpf_ILLraw_default_upper (lp, i);
	    if (prtLower) {
		if (mpf_EGlpNumIsEqqual (lp->lower[i], mpf_ILL_MINDOUBLE)) {
		    mpf_ILLprint_report (lp, " MI BOUND    %s\n", colnames[ri]);
		} else {
		    mpf_ILLprint_report (lp, " LO BOUND    %s    %g\n", colnames[ri],
			mpf_EGlpNumToLf (lp->lower[i]));
		}
	    }
	    if (prtUpper) {
		if (mpf_EGlpNumIsEqqual (lp->upper[i], mpf_ILL_MAXDOUBLE)) {
		    mpf_ILLprint_report (lp, " PL BOUND    %s\n", colnames[ri]);
		} else {
		    mpf_ILLprint_report (lp, " UP BOUND    %s    %g\n", colnames[ri],
			mpf_EGlpNumToLf (lp->upper[i]));
		}
	    }
	}
    }
    mpf_ILLprint_report (lp, "ENDATA\n");

CLEANUP:
    mpf_ILLlp_rows_clear (lprows);
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

static int mpf_mps_write_col (int i,
      int iorig,
      char *colname,
      mpf_ILLlpdata * lp,
      char **rownames,
      int intmode,
      char *objname)
{
    int row, k;
    mpf_ILLmatrix *A;

    A = &lp->A;
    if (lp->intmarker && (lp->intmarker[iorig] != intmode)) {
	mpf_ILLprint_report (lp, " MARK%dqs      'MARKER'    '%s'\n", iorig,
	    (lp->intmarker[iorig] ? "INTORG" : "INTEND"));
	intmode = lp->intmarker[iorig];
    }
    if (mpf_EGlpNumIsNeqqZero (lp->obj[i])) {
	mpf_ILLprint_report (lp, "  %s    %s    %g\n", colname, objname,
	    mpf_EGlpNumToLf (lp->obj[i]));
    }
    for (k = A->matbeg[i]; k < A->matbeg[i] + A->matcnt[i]; k++) {
	row = A->matind[k];
	mpf_ILLprint_report (lp, "  %s    %s    %g\n", colname, rownames[row],
	    mpf_EGlpNumToLf (A->matval[k]));
    }
    return intmode;
}
