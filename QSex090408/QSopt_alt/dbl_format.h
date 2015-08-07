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

/* RCSINFO $Id: dbl_format.h,v 1.3 2003/11/05 16:59:48 meven Exp $ */
#ifndef dbl_QS_FORMAT_ERROR_H
#define dbl_QS_FORMAT_ERROR_H

#include <stdio.h>
#include "dbl_qsopt.h"

/****************************************************************************/
/*
   The LP/MPS readers, writers, 
       dbl_ILLrawlpdata_to_lpdata, and 
   use dbl_ILLformat_error to report problems with their input iff
       the line reader used in reading the problem  or 
       the  dbl_qserror_collector pointer passed to dbl_ILLwrite_lp_file
   is not NULL.

   The QSgui code uses this feature to collect dbl_qsformat_error instances 
   which it uses after reading is done to insert error messages into the 
   input window. 
*/
/****************************************************************************/

/* 
for error type USE: 
          QS_DATA_ERROR			
          QS_DATA_WARN			
          QS_MPS_FORMAT_ERROR		
          QS_MPS_FORMAT_WARN		
          QS_LP_FORMAT_ERROR		
          QS_LP_FORMAT_WARN		
          QS_LP_OBJ_WARN			
          QS_GENERIC_ERROR		
*/

typedef struct dbl_qsformat_error
{
	char *desc;
	char *theLine;
	struct dbl_qsformat_error *next;
	int type;
	int lineNumber;								/* 1 based line counting */
	int at;
}
dbl_qsformat_error;

extern int dbl_ILLformat_error_create (dbl_qsformat_error * error,
																	 int mode,
																	 const char *desc,
																	 int lineNum,
																	 const char *theLine,
																	 int atPos);
extern void dbl_ILLformat_error_delete (dbl_qsformat_error * error);

extern void dbl_ILLformat_error_print (FILE * out,
																	 dbl_qsformat_error * e);



/*****************************************************************************
 * collecting error messages 
 * either with defining own qsad_error_fct and corresponding data structure 
 * or by using predefined dbl_ILLadd_error_to_memory fct with dbl_qserror_memory
 */

typedef int (*dbl_qsadd_error_fct) (void *dest,
																const dbl_qsformat_error * error);

typedef struct dbl_qserror_collector
{
	dbl_qsadd_error_fct add_error;
	void *dest;
}
dbl_qserror_collector;

typedef struct dbl_qserror_memory
{
	unsigned int nerror;
	dbl_qsformat_error *error_list;
	char has_error[QS_INPUT_NERROR];
	char hasErrorLines;
}
dbl_qserror_memory;


extern dbl_qserror_collector *dbl_ILLerror_collector_new (dbl_qsadd_error_fct fct,
																									void *dest);

dbl_qserror_collector *dbl_ILLerror_memory_collector_new (dbl_qserror_memory * dest);

extern void dbl_ILLerror_collector_free (dbl_qserror_collector * c);

#define dbl_ILLformat_error(collector, error)  \
	((collector)->add_error((collector)->dest, error))


extern int dbl_ILLadd_error_to_memory (void *dest,
																	 const dbl_qsformat_error * error);

extern dbl_qserror_memory *dbl_ILLerror_memory_create (int takeErrorLines);
extern void dbl_ILLerror_memory_free (dbl_qserror_memory * mem);

#endif
