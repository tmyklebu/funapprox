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

/*  RCS_INFO = "$RCSfile: read_lp_state.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef dbl_READ_LP_STATE_H
#define dbl_READ_LP_STATE_H

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading LP Files                       */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "dbl_iqsutil.h"
#include "dbl_readline.h"

typedef struct dbl_ILLread_lp_state
{
	dbl_qsline_reader *file;
	const char *file_name;
	char *p;
	double bound_val;
	int dbl_interactive;
	int line_num;
	int column_index;
	char realline[ILL_namebufsize];
	char line[ILL_namebufsize];
	char field[ILL_namebufsize + 1];
	char fieldOnFirstCol;
	char eof;
	char sense_val;
}
dbl_ILLread_lp_state;

extern int dbl_ILLread_lp_state_init (dbl_ILLread_lp_state * state,
																	dbl_qsline_reader * file,
																	const char *dbl_fname,
																	int interactve);
extern int dbl_ILLread_lp_state_next_line (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_next_var (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_keyword (dbl_ILLread_lp_state * state,
																		 const char **kwd);
extern int dbl_ILLread_lp_state_bad_keyword (dbl_ILLread_lp_state * state);
extern int dbl_ILLtest_lp_state_keyword (dbl_ILLread_lp_state * state,
																		 const char *kwd[]);
extern int dbl_ILLread_lp_state_next_field (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_next_field_on_line (dbl_ILLread_lp_state * state);
extern void dbl_ILLread_lp_state_prev_field (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_sign (dbl_ILLread_lp_state * state,
																	double * sign);
extern int dbl_ILLread_lp_state_possible_coef (dbl_ILLread_lp_state * state,
																					 double * coef,
																					 double defValue);
																				/* returns 1 iff found a number 
																				 * otherwise 0 */
extern int dbl_ILLread_lp_state_possible_bound_value (dbl_ILLread_lp_state * state);
																							 /* returns 1 iff found a number 
																							  * otherwise 0 */
extern int dbl_ILLread_lp_state_colon (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_has_colon (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_statxe_has_colon (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_next_constraint (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_sense (dbl_ILLread_lp_state * state);
extern int dbl_ILLtest_lp_state_sense (dbl_ILLread_lp_state * state,
																	 int all);
extern void dbl_ILLtest_lp_state_bound_sense (dbl_ILLread_lp_state * state);
extern int dbl_ILLread_lp_state_value (dbl_ILLread_lp_state * state,
																	 double * d);
extern int dbl_ILLtest_lp_state_next_is (dbl_ILLread_lp_state * state,
																		 const char *str);
extern int dbl_ILLread_lp_state_skip_blanks (dbl_ILLread_lp_state * state,
																				 int wrapLines);

extern int dbl_ILLcheck_subject_to (dbl_ILLread_lp_state * state);

/*---------------------------------------------------------------------------*/
/* errors and warnings 
 */
extern int dbl_ILLlp_error (dbl_ILLread_lp_state * state,
												const char *format,
												...);
extern void dbl_ILLlp_warn (dbl_ILLread_lp_state * state,
												const char *format,
												...);

/*---------------------------------------------------------------------------*/
/* shared with read_mps_state.c 
 */
extern int dbl_ILLget_value (char *line,
												 double * coef);

#endif
