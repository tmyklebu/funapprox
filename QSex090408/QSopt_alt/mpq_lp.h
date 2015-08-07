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

/* RCSINFO $Id: mpq_lp.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpq_LP_H
#define mpq_LP_H

#include "mpq_readline.h"

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

#include "mpq_lpdata.h"
#include "mpq_rawlp.h"
#include "mpq_read_lp.h"
#include "mpq_write_lp.h"

extern int mpq_ILLread_lp (mpq_qsline_reader * file,
											 const char *mpq_fname,
											 mpq_rawlpdata * lp);
extern int mpq_ILLwrite_lp (mpq_ILLlpdata * l,
												mpq_qserror_collector * collector);
			/* write using current lp->reporter */
extern int mpq_ILLis_lp_name_char (int c,
															 int pos);

extern int mpq_ILLread_constraint_name (mpq_ILLread_lp_state * state,
																		char **rowname);
extern int mpq_ILLread_one_constraint (mpq_ILLread_lp_state * state,
																	 const char *rowname,
																	 mpq_rawlpdata * lp,
																	 int allowNewColsAddRow);
extern int mpq_ILLread_constraint_expr (mpq_ILLread_lp_state * state,
																		mpq_rawlpdata * lp,
																		int rowind,
																		int allowNew);

#endif
