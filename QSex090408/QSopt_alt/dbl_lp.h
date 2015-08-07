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

/* RCSINFO $Id: dbl_lp.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef dbl_LP_H
#define dbl_LP_H

#include "dbl_readline.h"

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

#include "dbl_lpdata.h"
#include "dbl_rawlp.h"
#include "dbl_read_lp.h"
#include "dbl_write_lp.h"

extern int dbl_ILLread_lp (dbl_qsline_reader * file,
											 const char *dbl_fname,
											 dbl_rawlpdata * lp);
extern int dbl_ILLwrite_lp (dbl_ILLlpdata * l,
												dbl_qserror_collector * collector);
			/* write using current lp->reporter */
extern int dbl_ILLis_lp_name_char (int c,
															 int pos);

extern int dbl_ILLread_constraint_name (dbl_ILLread_lp_state * state,
																		char **rowname);
extern int dbl_ILLread_one_constraint (dbl_ILLread_lp_state * state,
																	 const char *rowname,
																	 dbl_rawlpdata * lp,
																	 int allowNewColsAddRow);
extern int dbl_ILLread_constraint_expr (dbl_ILLread_lp_state * state,
																		dbl_rawlpdata * lp,
																		int rowind,
																		int allowNew);

#endif
