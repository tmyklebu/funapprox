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

/* RCSINFO $Id: mpq_lpdata.h,v 1.4 2003/11/05 17:00:56 meven Exp $ */
#ifndef mpq_ILL_LPDATA_H
#define mpq_ILL_LPDATA_H

#include "config.h"
#include "mpq_qstruct.h"
#include "mpq_iqsutil.h"
#include "mpq_readline.h"
#include "reporter.h"
#include "mpq_format.h"
#include "mpq_dstruct.h"

extern mpq_t mpq_ILL_MAXDOUBLE;	/*  1267650600228229401496703205376.0 this is equal to  2^100 =  1.267e30 */
extern mpq_t mpq_ILL_MINDOUBLE;	/* -1267650600228229401496703205376.0 this is equal to -2^100 = -1.267e30 */
#define mpq_ILL_MAXINT    (2147483647)	/* this is equal to 2^31-1 */
#define mpq_ILL_MIN       (1)				/* Must be same as QS_MIN */
#define mpq_ILL_MAX       (-1)			/* Must be same as QS_MAX */

/*  Setting Alg in Presolve  */

#define mpq_ILL_PRE_SCALE           1
#define mpq_ILL_PRE_FIXED           2
#define mpq_ILL_PRE_SINGLE_ROW      4
#define mpq_ILL_PRE_FORCING         8
#define mpq_ILL_PRE_SINGLE_COL     16
#define mpq_ILL_PRE_DUPLICATE_ROW  32
#define mpq_ILL_PRE_DUPLICATE_COL  64
#define mpq_ILL_PRE_EMPTY_COL     128
#define mpq_ILL_PRE_ALL (mpq_ILL_PRE_SCALE | mpq_ILL_PRE_FIXED | mpq_ILL_PRE_SINGLE_ROW           \
                    mpq_ILL_PRE_FORCING | mpq_ILL_PRE_SINGLE_COL | mpq_ILL_PRE_DUPLICATE_ROW \
                   mpq_ILL_PRE_DUPLICATE_COL | mpq_ILL_PRE_EMPTY_COL)
#define mpq_ILL_PRE_SIMPLE (mpq_ILL_PRE_FIXED | mpq_ILL_PRE_EMPTY_COL)

typedef struct mpq_ILLlpdata
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
	mpq_t *obj;
	mpq_t *rhs;
	mpq_t *rangeval;
	mpq_t *lower;
	mpq_t *upper;
	mpq_ILLmatrix A;									/* The coef matrix.                        */
	struct mpq_ILLlp_rows *rA;				/* Coef matrix in row form.                */

	char **rownames;
	ILLsymboltab rowtab;					/* contains rownames in no particular order */
	char *objname;								/* if colname is not NULL it is entered into 
																 * the rowtab, see reader fcts in mpq_lp.c, mpq_mps.c*/

	char **colnames;							/* columns of struct variables */
	ILLsymboltab coltab;					/* contains colnames in no particular order */

	char *probname;
	char *intmarker;
	int *structmap;								/* Indices of structural variables         */
	int *rowmap;									/* Indices of logical and range variables  */
	struct mpq_ILLlp_basis *basis;
	struct mpq_ILLlp_predata *presolve;
	struct mpq_ILLlp_sinfo *sinfo;

	 /**************************************************************************/
	/* these fields are currently only set by mpq_mps.c reader fcts               */
	 /**************************************************************************/
	mpq_ILLmatrix sos;								/* columns are the sets, rows are the  
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
    * mpq_QSset_reporter initializes reporter 
    **************************************************************************/
	qsstring_reporter reporter;		/* used from within ILL fcts 
																 * to report mpq_feedback */
}
mpq_ILLlpdata;

typedef struct mpq_ILLlp_basis
{
	int nstruct;
	int nrows;
	int rownorms_size;
	int colnorms_size;
	char *cstat;
	char *rstat;
	mpq_t *rownorms;
	mpq_t *colnorms;
}
mpq_ILLlp_basis;

typedef struct mpq_ILLlp_cache
{
	int nstruct;
	int nrows;
	int status;
	mpq_t val;
	mpq_t *x;
	mpq_t *pi;
	mpq_t *rc;
	mpq_t *slack;
}
mpq_ILLlp_cache;

typedef struct mpq_ILLlp_sinfo
{																/* LP info returned by presolve            */
	int ncols;
	int nrows;
	int nzcount;
	int rowsize;
	int colsize;
	int objsense;

	mpq_t *obj;
	mpq_t *rhs;
	mpq_t *lower;
	mpq_t *upper;

	mpq_ILLmatrix A;

	char **colnames;							/* Just for debugging - not updated */
}
mpq_ILLlp_sinfo;

typedef struct mpq_ILLlp_preline
{
	mpq_t rhs;
	mpq_t obj;
	mpq_t lower;
	mpq_t upper;
	int count;
	int *ind;
	int row_or_col;								/* 0 is row, 1 is col */
	mpq_t *val;
}
mpq_ILLlp_preline;

typedef struct mpq_ILLlp_preop
{
	int ptype;
	int rowindex;
	int colindex;
	mpq_ILLlp_preline line;
}
mpq_ILLlp_preop;

typedef struct mpq_ILLlp_predata
{																/* Data needed in un-presolve.            */
	int opcount;
	int opsize;
	mpq_ILLlp_preop *oplist;
	int r_nrows;
	int r_ncols;
	int *colmap;
	int *rowmap;
	mpq_t *rowscale;
	mpq_t *colscale;
	mpq_t *colfixval;
	mpq_t *rowfixval;
}
mpq_ILLlp_predata;

typedef struct mpq_ILLlp_rows
{
	int *rowbeg;
	int *rowcnt;
	int *rowind;
	mpq_t *rowval;
}
mpq_ILLlp_rows;


/****************************************************************************/
/*                                                                          */
/*                             mpq_lpdata.c                                     */
/*                                                                          */
/****************************************************************************/

struct mpq_qsdata *mpq_ILLread (mpq_qsline_reader * file,
												const char *mpq_fname,
												int isMps);
void mpq_ILLstart (void);	/**< initialize mpq_ILL_MAXDOUBLE and other 
													 constants, this funtion should be callef AFTER 
													 EGlpNumStart() */
void mpq_ILLend (void);	/**< free any internal data asociated with variable 
												 mpq_precision numbers */
void mpq_ILLchange_precision (void);/**< This function re-compute the internal 
																		 variables mpq_precision to the (previously 
																		 set) EGLPNUM_PRECISION value (done with 
																		 EGlpNumSetPrecision) */
void mpq_ILLlpdata_init (mpq_ILLlpdata * lp);
void mpq_ILLlpdata_free (mpq_ILLlpdata * lp);
void mpq_ILLlp_basis_init (mpq_ILLlp_basis * B);
void mpq_ILLlp_basis_free (mpq_ILLlp_basis * B);
void mpq_ILLlp_cache_init (mpq_ILLlp_cache * C);
void mpq_ILLlp_cache_free (mpq_ILLlp_cache * C);
int mpq_ILLlp_basis_alloc (mpq_ILLlp_basis * B,
											 int ncols,
											 int nrows);
int mpq_ILLlp_cache_alloc (mpq_ILLlp_cache * C,
											 int ncols,
											 int nrows);

int mpq_ILLlp_rows_init (mpq_ILLlp_rows * lp_rows,
										 mpq_ILLlpdata * lp,
										 int include_logicals);
void mpq_ILLlp_rows_clear (mpq_ILLlp_rows * lp_rows);
int mpq_ILLprint_report (mpq_ILLlpdata * lp,
										 const char *format,
										 ...);
							/* print to lp->reporter */

/****************************************************************************/
/*                                                                          */
/*                             mpq_presolve.c                                   */
/*                                                                          */
/****************************************************************************/

void mpq_ILLlp_sinfo_init (mpq_ILLlp_sinfo * sinfo),
  mpq_ILLlp_sinfo_free (mpq_ILLlp_sinfo * sinfo),
  mpq_ILLlp_predata_init (mpq_ILLlp_predata * pre),
  mpq_ILLlp_predata_free (mpq_ILLlp_predata * pre);

int mpq_ILLlp_add_logicals (mpq_ILLlpdata * lp),
  mpq_ILLlp_scale (mpq_ILLlpdata * lp),
  mpq_ILLlp_presolve (mpq_ILLlpdata * lp,
									int pre_types);


#endif /* __ILL_LPDATA_H */
