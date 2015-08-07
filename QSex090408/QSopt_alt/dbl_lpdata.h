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

/* RCSINFO $Id: dbl_lpdata.h,v 1.4 2003/11/05 17:00:56 meven Exp $ */
#ifndef dbl_ILL_LPDATA_H
#define dbl_ILL_LPDATA_H

#include "config.h"
#include "dbl_qstruct.h"
#include "dbl_iqsutil.h"
#include "dbl_readline.h"
#include "reporter.h"
#include "dbl_format.h"
#include "dbl_dstruct.h"

extern double dbl_ILL_MAXDOUBLE;	/*  1267650600228229401496703205376.0 this is equal to  2^100 =  1.267e30 */
extern double dbl_ILL_MINDOUBLE;	/* -1267650600228229401496703205376.0 this is equal to -2^100 = -1.267e30 */
#define dbl_ILL_MAXINT    (2147483647)	/* this is equal to 2^31-1 */
#define dbl_ILL_MIN       (1)				/* Must be same as QS_MIN */
#define dbl_ILL_MAX       (-1)			/* Must be same as QS_MAX */

/*  Setting Alg in Presolve  */

#define dbl_ILL_PRE_SCALE           1
#define dbl_ILL_PRE_FIXED           2
#define dbl_ILL_PRE_SINGLE_ROW      4
#define dbl_ILL_PRE_FORCING         8
#define dbl_ILL_PRE_SINGLE_COL     16
#define dbl_ILL_PRE_DUPLICATE_ROW  32
#define dbl_ILL_PRE_DUPLICATE_COL  64
#define dbl_ILL_PRE_EMPTY_COL     128
#define dbl_ILL_PRE_ALL (dbl_ILL_PRE_SCALE | dbl_ILL_PRE_FIXED | dbl_ILL_PRE_SINGLE_ROW           \
                    dbl_ILL_PRE_FORCING | dbl_ILL_PRE_SINGLE_COL | dbl_ILL_PRE_DUPLICATE_ROW \
                   dbl_ILL_PRE_DUPLICATE_COL | dbl_ILL_PRE_EMPTY_COL)
#define dbl_ILL_PRE_SIMPLE (dbl_ILL_PRE_FIXED | dbl_ILL_PRE_EMPTY_COL)

typedef struct dbl_ILLlpdata
{																/* Complete LP data filled in by mpsread.  */
	int nrows;
	int ncols;
	int nstruct;									/* Not including logicals.                 */
	int nzcount;
	int rowsize;									/* Length of row arrays.                   */
	int colsize;									/* Length of col arrays.                   */
	int structsize;								/* Length of intmarker, structmap,         */
	/* colnames                                */
	int objsense;
	char *sense;									/* Original sense, not after logicals.     */
	double *obj;
	double *rhs;
	double *rangeval;
	double *lower;
	double *upper;
	dbl_ILLmatrix A;									/* The coef matrix.                        */
	struct dbl_ILLlp_rows *rA;				/* Coef matrix in row form.                */

	char **rownames;
	ILLsymboltab rowtab;					/* contains rownames in no particular order */
	char *objname;								/* if colname is not NULL it is entered into 
																 * the rowtab, see reader fcts in dbl_lp.c, dbl_mps.c*/

	char **colnames;							/* columns of struct variables */
	ILLsymboltab coltab;					/* contains colnames in no particular order */

	char *probname;
	char *intmarker;
	int *structmap;								/* Indices of structural variables         */
	int *rowmap;									/* Indices of logical and range variables  */
	struct dbl_ILLlp_basis *basis;
	struct dbl_ILLlp_predata *presolve;
	struct dbl_ILLlp_sinfo *sinfo;

	 /**************************************************************************/
	/* these fields are currently only set by dbl_mps.c reader fcts               */
	 /**************************************************************************/
	dbl_ILLmatrix sos;								/* columns are the sets, rows are the  
																 * problem's structural variables
																 * coefficients are the weights */

	char *sos_type;								/* type of each set */
	int *is_sos_mem;							/* for each structural variable contains 
																 *    -1 == not a set member
																 *     i == member of sos set i 
																 *          where 0 <= i < sos.matcols */
	char *refrowname;							/* name of reference row */
	int refind;										/* index of reference row 
																 *     -1 if refrow was a free row 
																 *          and weights are found only in the 
																 *          sos matrix 
																 *     index >=0 if refrow is also a lp-row */

	 /**************************************************************************
    * dbl_QSset_reporter initializes reporter 
    **************************************************************************/
	qsstring_reporter reporter;		/* used from within ILL fcts 
																 * to report dbl_feedback */
}
dbl_ILLlpdata;

typedef struct dbl_ILLlp_basis
{
	int nstruct;
	int nrows;
	int rownorms_size;
	int colnorms_size;
	char *cstat;
	char *rstat;
	double *rownorms;
	double *colnorms;
}
dbl_ILLlp_basis;

typedef struct dbl_ILLlp_cache
{
	int nstruct;
	int nrows;
	int status;
	double val;
	double *x;
	double *pi;
	double *rc;
	double *slack;
}
dbl_ILLlp_cache;

typedef struct dbl_ILLlp_sinfo
{																/* LP info returned by presolve            */
	int ncols;
	int nrows;
	int nzcount;
	int rowsize;
	int colsize;
	int objsense;

	double *obj;
	double *rhs;
	double *lower;
	double *upper;

	dbl_ILLmatrix A;

	char **colnames;							/* Just for debugging - not updated */
}
dbl_ILLlp_sinfo;

typedef struct dbl_ILLlp_preline
{
	double rhs;
	double obj;
	double lower;
	double upper;
	int count;
	int *ind;
	int row_or_col;								/* 0 is row, 1 is col */
	double *val;
}
dbl_ILLlp_preline;

typedef struct dbl_ILLlp_preop
{
	int ptype;
	int rowindex;
	int colindex;
	dbl_ILLlp_preline line;
}
dbl_ILLlp_preop;

typedef struct dbl_ILLlp_predata
{																/* Data needed in un-presolve.            */
	int opcount;
	int opsize;
	dbl_ILLlp_preop *oplist;
	int r_nrows;
	int r_ncols;
	int *colmap;
	int *rowmap;
	double *rowscale;
	double *colscale;
	double *colfixval;
	double *rowfixval;
}
dbl_ILLlp_predata;

typedef struct dbl_ILLlp_rows
{
	int *rowbeg;
	int *rowcnt;
	int *rowind;
	double *rowval;
}
dbl_ILLlp_rows;


/****************************************************************************/
/*                                                                          */
/*                             dbl_lpdata.c                                     */
/*                                                                          */
/****************************************************************************/

struct dbl_qsdata *dbl_ILLread (dbl_qsline_reader * file,
												const char *dbl_fname,
												int isMps);
void dbl_ILLstart (void);	/**< initialize dbl_ILL_MAXDOUBLE and other 
													 constants, this funtion should be callef AFTER 
													 EGlpNumStart() */
void dbl_ILLend (void);	/**< free any internal data asociated with variable 
												 dbl_precision numbers */
void dbl_ILLchange_precision (void);/**< This function re-compute the internal 
																		 variables dbl_precision to the (previously 
																		 set) EGLPNUM_PRECISION value (done with 
																		 EGlpNumSetPrecision) */
void dbl_ILLlpdata_init (dbl_ILLlpdata * lp);
void dbl_ILLlpdata_free (dbl_ILLlpdata * lp);
void dbl_ILLlp_basis_init (dbl_ILLlp_basis * B);
void dbl_ILLlp_basis_free (dbl_ILLlp_basis * B);
void dbl_ILLlp_cache_init (dbl_ILLlp_cache * C);
void dbl_ILLlp_cache_free (dbl_ILLlp_cache * C);
int dbl_ILLlp_basis_alloc (dbl_ILLlp_basis * B,
											 int ncols,
											 int nrows);
int dbl_ILLlp_cache_alloc (dbl_ILLlp_cache * C,
											 int ncols,
											 int nrows);

int dbl_ILLlp_rows_init (dbl_ILLlp_rows * lp_rows,
										 dbl_ILLlpdata * lp,
										 int include_logicals);
void dbl_ILLlp_rows_clear (dbl_ILLlp_rows * lp_rows);
int dbl_ILLprint_report (dbl_ILLlpdata * lp,
										 const char *format,
										 ...);
							/* print to lp->reporter */

/****************************************************************************/
/*                                                                          */
/*                             dbl_presolve.c                                   */
/*                                                                          */
/****************************************************************************/

void dbl_ILLlp_sinfo_init (dbl_ILLlp_sinfo * sinfo),
  dbl_ILLlp_sinfo_free (dbl_ILLlp_sinfo * sinfo),
  dbl_ILLlp_predata_init (dbl_ILLlp_predata * pre),
  dbl_ILLlp_predata_free (dbl_ILLlp_predata * pre);

int dbl_ILLlp_add_logicals (dbl_ILLlpdata * lp),
  dbl_ILLlp_scale (dbl_ILLlpdata * lp),
  dbl_ILLlp_presolve (dbl_ILLlpdata * lp,
									int pre_types);


#endif /* __ILL_LPDATA_H */
