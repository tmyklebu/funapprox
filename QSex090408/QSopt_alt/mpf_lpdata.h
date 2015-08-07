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

/* RCSINFO $Id: mpf_lpdata.h,v 1.4 2003/11/05 17:00:56 meven Exp $ */
#ifndef mpf_ILL_LPDATA_H
#define mpf_ILL_LPDATA_H

#include "config.h"
#include "mpf_qstruct.h"
#include "mpf_iqsutil.h"
#include "mpf_readline.h"
#include "reporter.h"
#include "mpf_format.h"
#include "mpf_dstruct.h"

extern mpf_t mpf_ILL_MAXDOUBLE;	/*  1267650600228229401496703205376.0 this is equal to  2^100 =  1.267e30 */
extern mpf_t mpf_ILL_MINDOUBLE;	/* -1267650600228229401496703205376.0 this is equal to -2^100 = -1.267e30 */
#define mpf_ILL_MAXINT    (2147483647)	/* this is equal to 2^31-1 */
#define mpf_ILL_MIN       (1)				/* Must be same as QS_MIN */
#define mpf_ILL_MAX       (-1)			/* Must be same as QS_MAX */

/*  Setting Alg in Presolve  */

#define mpf_ILL_PRE_SCALE           1
#define mpf_ILL_PRE_FIXED           2
#define mpf_ILL_PRE_SINGLE_ROW      4
#define mpf_ILL_PRE_FORCING         8
#define mpf_ILL_PRE_SINGLE_COL     16
#define mpf_ILL_PRE_DUPLICATE_ROW  32
#define mpf_ILL_PRE_DUPLICATE_COL  64
#define mpf_ILL_PRE_EMPTY_COL     128
#define mpf_ILL_PRE_ALL (mpf_ILL_PRE_SCALE | mpf_ILL_PRE_FIXED | mpf_ILL_PRE_SINGLE_ROW           \
                    mpf_ILL_PRE_FORCING | mpf_ILL_PRE_SINGLE_COL | mpf_ILL_PRE_DUPLICATE_ROW \
                   mpf_ILL_PRE_DUPLICATE_COL | mpf_ILL_PRE_EMPTY_COL)
#define mpf_ILL_PRE_SIMPLE (mpf_ILL_PRE_FIXED | mpf_ILL_PRE_EMPTY_COL)

typedef struct mpf_ILLlpdata
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
	mpf_t *obj;
	mpf_t *rhs;
	mpf_t *rangeval;
	mpf_t *lower;
	mpf_t *upper;
	mpf_ILLmatrix A;									/* The coef matrix.                        */
	struct mpf_ILLlp_rows *rA;				/* Coef matrix in row form.                */

	char **rownames;
	ILLsymboltab rowtab;					/* contains rownames in no particular order */
	char *objname;								/* if colname is not NULL it is entered into 
																 * the rowtab, see reader fcts in mpf_lp.c, mpf_mps.c*/

	char **colnames;							/* columns of struct variables */
	ILLsymboltab coltab;					/* contains colnames in no particular order */

	char *probname;
	char *intmarker;
	int *structmap;								/* Indices of structural variables         */
	int *rowmap;									/* Indices of logical and range variables  */
	struct mpf_ILLlp_basis *basis;
	struct mpf_ILLlp_predata *presolve;
	struct mpf_ILLlp_sinfo *sinfo;

	 /**************************************************************************/
	/* these fields are currently only set by mpf_mps.c reader fcts               */
	 /**************************************************************************/
	mpf_ILLmatrix sos;								/* columns are the sets, rows are the  
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
    * mpf_QSset_reporter initializes reporter 
    **************************************************************************/
	qsstring_reporter reporter;		/* used from within ILL fcts 
																 * to report mpf_feedback */
}
mpf_ILLlpdata;

typedef struct mpf_ILLlp_basis
{
	int nstruct;
	int nrows;
	int rownorms_size;
	int colnorms_size;
	char *cstat;
	char *rstat;
	mpf_t *rownorms;
	mpf_t *colnorms;
}
mpf_ILLlp_basis;

typedef struct mpf_ILLlp_cache
{
	int nstruct;
	int nrows;
	int status;
	mpf_t val;
	mpf_t *x;
	mpf_t *pi;
	mpf_t *rc;
	mpf_t *slack;
}
mpf_ILLlp_cache;

typedef struct mpf_ILLlp_sinfo
{																/* LP info returned by presolve            */
	int ncols;
	int nrows;
	int nzcount;
	int rowsize;
	int colsize;
	int objsense;

	mpf_t *obj;
	mpf_t *rhs;
	mpf_t *lower;
	mpf_t *upper;

	mpf_ILLmatrix A;

	char **colnames;							/* Just for debugging - not updated */
}
mpf_ILLlp_sinfo;

typedef struct mpf_ILLlp_preline
{
	mpf_t rhs;
	mpf_t obj;
	mpf_t lower;
	mpf_t upper;
	int count;
	int *ind;
	int row_or_col;								/* 0 is row, 1 is col */
	mpf_t *val;
}
mpf_ILLlp_preline;

typedef struct mpf_ILLlp_preop
{
	int ptype;
	int rowindex;
	int colindex;
	mpf_ILLlp_preline line;
}
mpf_ILLlp_preop;

typedef struct mpf_ILLlp_predata
{																/* Data needed in un-presolve.            */
	int opcount;
	int opsize;
	mpf_ILLlp_preop *oplist;
	int r_nrows;
	int r_ncols;
	int *colmap;
	int *rowmap;
	mpf_t *rowscale;
	mpf_t *colscale;
	mpf_t *colfixval;
	mpf_t *rowfixval;
}
mpf_ILLlp_predata;

typedef struct mpf_ILLlp_rows
{
	int *rowbeg;
	int *rowcnt;
	int *rowind;
	mpf_t *rowval;
}
mpf_ILLlp_rows;


/****************************************************************************/
/*                                                                          */
/*                             mpf_lpdata.c                                     */
/*                                                                          */
/****************************************************************************/

struct mpf_qsdata *mpf_ILLread (mpf_qsline_reader * file,
												const char *mpf_fname,
												int isMps);
void mpf_ILLstart (void);	/**< initialize mpf_ILL_MAXDOUBLE and other 
													 constants, this funtion should be callef AFTER 
													 EGlpNumStart() */
void mpf_ILLend (void);	/**< free any internal data asociated with variable 
												 mpf_precision numbers */
void mpf_ILLchange_precision (void);/**< This function re-compute the internal 
																		 variables mpf_precision to the (previously 
																		 set) EGLPNUM_PRECISION value (done with 
																		 EGlpNumSetPrecision) */
void mpf_ILLlpdata_init (mpf_ILLlpdata * lp);
void mpf_ILLlpdata_free (mpf_ILLlpdata * lp);
void mpf_ILLlp_basis_init (mpf_ILLlp_basis * B);
void mpf_ILLlp_basis_free (mpf_ILLlp_basis * B);
void mpf_ILLlp_cache_init (mpf_ILLlp_cache * C);
void mpf_ILLlp_cache_free (mpf_ILLlp_cache * C);
int mpf_ILLlp_basis_alloc (mpf_ILLlp_basis * B,
											 int ncols,
											 int nrows);
int mpf_ILLlp_cache_alloc (mpf_ILLlp_cache * C,
											 int ncols,
											 int nrows);

int mpf_ILLlp_rows_init (mpf_ILLlp_rows * lp_rows,
										 mpf_ILLlpdata * lp,
										 int include_logicals);
void mpf_ILLlp_rows_clear (mpf_ILLlp_rows * lp_rows);
int mpf_ILLprint_report (mpf_ILLlpdata * lp,
										 const char *format,
										 ...);
							/* print to lp->reporter */

/****************************************************************************/
/*                                                                          */
/*                             mpf_presolve.c                                   */
/*                                                                          */
/****************************************************************************/

void mpf_ILLlp_sinfo_init (mpf_ILLlp_sinfo * sinfo),
  mpf_ILLlp_sinfo_free (mpf_ILLlp_sinfo * sinfo),
  mpf_ILLlp_predata_init (mpf_ILLlp_predata * pre),
  mpf_ILLlp_predata_free (mpf_ILLlp_predata * pre);

int mpf_ILLlp_add_logicals (mpf_ILLlpdata * lp),
  mpf_ILLlp_scale (mpf_ILLlpdata * lp),
  mpf_ILLlp_presolve (mpf_ILLlpdata * lp,
									int pre_types);


#endif /* __ILL_LPDATA_H */
