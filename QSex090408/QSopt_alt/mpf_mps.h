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

/* RCSINFO $Id: mpf_mps.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpf_MPS_H
#define mpf_MPS_H

#include "mpf_readline.h"
#include "mpf_format.h"

/****************************************************************************/
/*                                                                          */
/*              Routines to support Reading and Writing MPS Files           */
/*                                                                          */
/****************************************************************************/
#include "basicdefs.h"
extern const char *mpf_ILLmps_section_name[ILL_MPS_N_SECTIONS + 2];

#include "mpf_lpdata.h"
#include "mpf_rawlp.h"
#include "mpf_read_mps.h"

extern int mpf_ILLread_mps (mpf_qsline_reader * file,
												const char *filename,
												mpf_rawlpdata * lp);

extern int mpf_ILLwrite_mps (mpf_ILLlpdata * lp,
												 mpf_qserror_collector * collector);
				/* use lp->reporter for output */

#endif
