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
#ifndef mpq_READ_LP_STATE_H
#define mpq_READ_LP_STATE_H

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

#include "mpq_iqsutil.h"
#include "mpq_readline.h"

typedef struct mpq_ILLread_lp_state
{
	mpq_qsline_reader *file;
	const char *file_name;
	char *p;
	mpq_t bound_val;
	int mpq_interactive;
	int line_num;
	int column_index;
	char realline[ILL_namebufsize];
	char line[ILL_namebufsize];
	char field[ILL_namebufsize + 1];
	char fieldOnFirstCol;
	char eof;
	char sense_val;
}
mpq_ILLread_lp_state;

extern int mpq_ILLread_lp_state_init (mpq_ILLread_lp_state * state,
																	mpq_qsline_reader * file,
																	const char *mpq_fname,
																	int interactve);
extern int mpq_ILLread_lp_state_next_line (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_next_var (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_keyword (mpq_ILLread_lp_state * state,
																		 const char **kwd);
extern int mpq_ILLread_lp_state_bad_keyword (mpq_ILLread_lp_state * state);
extern int mpq_ILLtest_lp_state_keyword (mpq_ILLread_lp_state * state,
																		 const char *kwd[]);
extern int mpq_ILLread_lp_state_next_field (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_next_field_on_line (mpq_ILLread_lp_state * state);
extern void mpq_ILLread_lp_state_prev_field (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_sign (mpq_ILLread_lp_state * state,
																	mpq_t * sign);
extern int mpq_ILLread_lp_state_possible_coef (mpq_ILLread_lp_state * state,
																					 mpq_t * coef,
																					 mpq_t defValue);
																				/* returns 1 iff found a number 
																				 * otherwise 0 */
extern int mpq_ILLread_lp_state_possible_bound_value (mpq_ILLread_lp_state * state);
																							 /* returns 1 iff found a number 
																							  * otherwise 0 */
extern int mpq_ILLread_lp_state_colon (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_has_colon (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_statxe_has_colon (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_next_constraint (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_sense (mpq_ILLread_lp_state * state);
extern int mpq_ILLtest_lp_state_sense (mpq_ILLread_lp_state * state,
																	 int all);
extern void mpq_ILLtest_lp_state_bound_sense (mpq_ILLread_lp_state * state);
extern int mpq_ILLread_lp_state_value (mpq_ILLread_lp_state * state,
																	 mpq_t * d);
extern int mpq_ILLtest_lp_state_next_is (mpq_ILLread_lp_state * state,
																		 const char *str);
extern int mpq_ILLread_lp_state_skip_blanks (mpq_ILLread_lp_state * state,
																				 int wrapLines);

extern int mpq_ILLcheck_subject_to (mpq_ILLread_lp_state * state);

/*---------------------------------------------------------------------------*/
/* errors and warnings 
 */
extern int mpq_ILLlp_error (mpq_ILLread_lp_state * state,
												const char *format,
												...);
extern void mpq_ILLlp_warn (mpq_ILLread_lp_state * state,
												const char *format,
												...);

/*---------------------------------------------------------------------------*/
/* shared with read_mps_state.c 
 */
extern int mpq_ILLget_value (char *line,
												 mpq_t * coef);

#endif
