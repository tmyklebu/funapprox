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

/*  $RCSfile: dbl_qstruct.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef dbl___QS_QSTRUCT_H
#define dbl___QS_QSTRUCT_H
#include "basicdefs.h"

typedef struct dbl_qsdata
{
	struct dbl_ILLlpdata *qslp;
	struct dbl_lpinfo *lp;
	struct dbl_price_info *pricing;
	struct dbl_ILLlp_basis *basis;
	struct dbl_ILLlp_cache *cache;
	char *name;
	int qstatus;									/* QS_LP_UNSOLVED or an opt status  */
	int factorok;									/* set to 0 if factorization is old */
	int simplex_display;					/* 0 off, 1 on */
	int simplex_scaling;					/* 0 off, 1 on */
}
dbl_QSdata;

#endif /* dbl___QS_QSTRUCT_H */
