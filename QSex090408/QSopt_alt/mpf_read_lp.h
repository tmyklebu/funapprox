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
#ifndef mpf_READ_LP_STATE_H
#define mpf_READ_LP_STATE_H

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

#include "mpf_iqsutil.h"
#include "mpf_readline.h"

typedef struct mpf_ILLread_lp_state
{
	mpf_qsline_reader *file;
	const char *file_name;
	char *p;
	mpf_t bound_val;
	int mpf_interactive;
	int line_num;
	int column_index;
	char realline[ILL_namebufsize];
	char line[ILL_namebufsize];
	char field[ILL_namebufsize + 1];
	char fieldOnFirstCol;
	char eof;
	char sense_val;
}
mpf_ILLread_lp_state;

extern int mpf_ILLread_lp_state_init (mpf_ILLread_lp_state * state,
																	mpf_qsline_reader * file,
																	const char *mpf_fname,
																	int interactve);
extern int mpf_ILLread_lp_state_next_line (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_next_var (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_keyword (mpf_ILLread_lp_state * state,
																		 const char **kwd);
extern int mpf_ILLread_lp_state_bad_keyword (mpf_ILLread_lp_state * state);
extern int mpf_ILLtest_lp_state_keyword (mpf_ILLread_lp_state * state,
																		 const char *kwd[]);
extern int mpf_ILLread_lp_state_next_field (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_next_field_on_line (mpf_ILLread_lp_state * state);
extern void mpf_ILLread_lp_state_prev_field (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_sign (mpf_ILLread_lp_state * state,
																	mpf_t * sign);
extern int mpf_ILLread_lp_state_possible_coef (mpf_ILLread_lp_state * state,
																					 mpf_t * coef,
																					 mpf_t defValue);
																				/* returns 1 iff found a number 
																				 * otherwise 0 */
extern int mpf_ILLread_lp_state_possible_bound_value (mpf_ILLread_lp_state * state);
																							 /* returns 1 iff found a number 
																							  * otherwise 0 */
extern int mpf_ILLread_lp_state_colon (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_has_colon (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_statxe_has_colon (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_next_constraint (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_sense (mpf_ILLread_lp_state * state);
extern int mpf_ILLtest_lp_state_sense (mpf_ILLread_lp_state * state,
																	 int all);
extern void mpf_ILLtest_lp_state_bound_sense (mpf_ILLread_lp_state * state);
extern int mpf_ILLread_lp_state_value (mpf_ILLread_lp_state * state,
																	 mpf_t * d);
extern int mpf_ILLtest_lp_state_next_is (mpf_ILLread_lp_state * state,
																		 const char *str);
extern int mpf_ILLread_lp_state_skip_blanks (mpf_ILLread_lp_state * state,
																				 int wrapLines);

extern int mpf_ILLcheck_subject_to (mpf_ILLread_lp_state * state);

/*---------------------------------------------------------------------------*/
/* errors and warnings 
 */
extern int mpf_ILLlp_error (mpf_ILLread_lp_state * state,
												const char *format,
												...);
extern void mpf_ILLlp_warn (mpf_ILLread_lp_state * state,
												const char *format,
												...);

/*---------------------------------------------------------------------------*/
/* shared with read_mps_state.c 
 */
extern int mpf_ILLget_value (char *line,
												 mpf_t * coef);

#endif
