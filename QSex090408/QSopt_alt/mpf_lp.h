/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */ 
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

/* RCSINFO $Id: mpf_lp.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpf_LP_H
#define mpf_LP_H

#include "mpf_readline.h"

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading and Writing LP Files           */
/*                                                                          */
/****************************************************************************/
/* 
 * -) anything after '\' is comment 
 * -) Problem is optional and comes first
 * -) Minimize, Maximize, ... comes after Problem
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "mpf_lpdata.h"
#include "mpf_rawlp.h"
#include "mpf_read_lp.h"
#include "mpf_write_lp.h"

extern int mpf_ILLread_lp (mpf_qsline_reader * file,
											 const char *mpf_fname,
											 mpf_rawlpdata * lp);
extern int mpf_ILLwrite_lp (mpf_ILLlpdata * l,
												mpf_qserror_collector * collector);
			/* write using current lp->reporter */
extern int mpf_ILLis_lp_name_char (int c,
															 int pos);

extern int mpf_ILLread_constraint_name (mpf_ILLread_lp_state * state,
																		char **rowname);
extern int mpf_ILLread_one_constraint (mpf_ILLread_lp_state * state,
																	 const char *rowname,
																	 mpf_rawlpdata * lp,
																	 int allowNewColsAddRow);
extern int mpf_ILLread_constraint_expr (mpf_ILLread_lp_state * state,
																		mpf_rawlpdata * lp,
																		int rowind,
																		int allowNew);

#endif
