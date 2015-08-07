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

/* RCSINFO $Id: mpq_mps.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpq_MPS_H
#define mpq_MPS_H

#include "mpq_readline.h"
#include "mpq_format.h"

/****************************************************************************/
/*                                                                          */
/*              Routines to support Reading and Writing MPS Files           */
/*                                                                          */
/****************************************************************************/
#include "basicdefs.h"
extern const char *mpq_ILLmps_section_name[ILL_MPS_N_SECTIONS + 2];

#include "mpq_lpdata.h"
#include "mpq_rawlp.h"
#include "mpq_read_mps.h"

extern int mpq_ILLread_mps (mpq_qsline_reader * file,
												const char *filename,
												mpq_rawlpdata * lp);

extern int mpq_ILLwrite_mps (mpq_ILLlpdata * lp,
												 mpq_qserror_collector * collector);
				/* use lp->reporter for output */

#endif
