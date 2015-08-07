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

/*  RCS_INFO = "$RCSfile: mpf_rawlp.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef mpf___ILL_RAWLP_H_
#define mpf___ILL_RAWLP_H_

/****************************************************************************/
/* DataStructure and Routines                                               */
/*          to deal with raw lp information as read from mps or lp files    */
/*          support scanning of input                                       */
/*          error reporting                                                 */
/****************************************************************************/

#include "trace.h"
#include "mpf_lpdata.h"
#include "mpf_iqsutil.h"
#include "mpf_format.h"
#include "mpf_lpdefs.h"

#define mpf_ILL_ISBLANK(p) \
             (((*(p))==' '||(*(p))=='\t'||(*(p))=='\r'||(*(p))=='\f') ? 1 : 0)

/* 
 * we rely on ILLsymboltab property:
 *   the ith name added can be retrieved by ILLsymboltab_get(table, i) 
 *   as long as we never delete names from the symbol table 
 */
typedef struct mpf_rawlpdata
{
	char *name;

	char *rhsname;
	char *rangesname;
	char *boundsname;

	int objsense;									/* maximize or minimize */
	int objindex;									/* index of objective row */

	int nrows;										/* number of rows in problem */
	ILLsymboltab rowtab;					/* ILLsymboltab_get(rowtab, i) name of ith row */

	int sensesize;								/* size of rowsense */
	char *rowsense;								/* rowsense[i] snese of row[i] */

	char *rhsind;									/* rhsind[i] == 1 we saw an rhs for row[i] */
	/* size is nrows */
	int rhssize;									/* size of rhs array */
	mpf_t *rhs;								/* rhs values for rows; size is nrows */
	char *rangesind;							/* ranges[i] == 1 we saw a range def for row[i] */
	struct mpf_colptr *ranges;				/* list of range values */

	int ncols;										/* number of cols in problem */
	ILLsymboltab coltab;					/* ILLsymboltab_get(coltab, i) name of ith col */
	int colsize;									/* size of cols array */
	struct mpf_colptr **cols;

	char *lbind;									/* lbind[i] == 1  we saw a lower bound for col[i] */
	char *ubind;									/* ubind[i] == 1  we saw a upper bound for col[i] */
	mpf_t *lower;							/* lower[i] = lower bound for col[i] */
	mpf_t *upper;							/* upper[i] = upper bound for col[i] */

	int intsize;									/* size of intmarker array */
	char *intmarker;							/* intmarker[i] == 1  col[i] is an int var */

	/* sos information is tranfered into mpf_ILLmatrix lpdata->sos */
	char *refrow;									/* name of reference row */
	int refrowind;								/* index of refrow or -1  */

	int is_sos_size;							/* size of is_sos_member array */
	int *is_sos_member;						/* for each col contains either               
																 *     -1 == no sos memeber 
																 *     i  == member of set #i */

	int nsos_member;							/* total number of sos set members */
	int sos_weight_size;					/* size of sos_weight array */
	mpf_t *sos_weight;				/* sos set elem i has weight of sos_weight[i] 
																 * value comes from refrow coeficients */
	int sos_col_size;							/* size of sos_col array */
	int *sos_col;									/* sos elem i is column sos_col[i] */

	int nsos;											/* number of sos sets */
	int sos_setsize;							/* size of sosset array */
	struct mpf_sosptr *sos_set;				/* type, size, first element of sos sets 
																 * first is index into sos_weight and sos_col 
																 * arrays */
	mpf_qserror_collector *error_collector;
	ILLptrworld ptrworld;
}
mpf_rawlpdata;

typedef struct mpf_colptr
{
	mpf_t coef;
	struct mpf_colptr *next;
	int this;											/* row index */
}
mpf_colptr;
extern mpf_colptr *mpf_ILLcolptralloc (ILLptrworld * p);

typedef struct mpf_sosptr
{
	int nelem;										/* number of set elements */
	int first;										/* index of first set element in sosmemeber */
	char type;										/* set type */
}
mpf_sosptr;
extern const int mpf_ILL_SOS_TYPE1;
extern const int mpf_ILL_SOS_TYPE2;

extern void mpf_ILLinit_rawlpdata (mpf_rawlpdata * lp,
															 mpf_qserror_collector * collector);
extern void mpf_ILLfree_rawlpdata (mpf_rawlpdata * lp);
extern void mpf_ILLraw_clear_matrix (mpf_rawlpdata * lp);

extern const char *mpf_ILLraw_rowname (mpf_rawlpdata * lp,
																	 int i);
extern const char *mpf_ILLraw_colname (mpf_rawlpdata * lp,
																	 int i);

extern int mpf_ILLraw_add_col (mpf_rawlpdata * lp,
													 const char *name,
													 int intmarker);
extern int mpf_ILLraw_add_row (mpf_rawlpdata * lp,
													 const char *name,
													 int sense,
													 mpf_t rhs);

extern int mpf_ILLraw_add_col_coef (mpf_rawlpdata * lp,
																int colind,
																int rowind,
																mpf_t coef);

extern int mpf_ILLraw_init_ranges (mpf_rawlpdata * lp);
extern int mpf_ILLraw_init_rhs (mpf_rawlpdata * lp);

extern int mpf_ILLraw_add_ranges_coef (mpf_rawlpdata * lp,
																	 int rowind,
																	 mpf_t coef);


extern int mpf_ILLraw_add_sos (mpf_rawlpdata * lp,
													 int sos_type);
																								/* add empty set with type */
extern int mpf_ILLraw_add_sos_member (mpf_rawlpdata * lp,
																	int colind);
																								/* add col to last set */
extern int mpf_ILLraw_is_mem_other_sos (mpf_rawlpdata * lp,
																		int colind);

extern int mpf_ILLraw_set_rhs_name (mpf_rawlpdata * lp,
																const char *name,
																int *skip);
extern int mpf_ILLraw_set_bounds_name (mpf_rawlpdata * lp,
																	 const char *name,
																	 int *skip);
extern int mpf_ILLraw_set_ranges_name (mpf_rawlpdata * lp,
																	 const char *name,
																	 int *skip);
extern void mpf_ILLprint_rawlpdata (mpf_rawlpdata * lp);

extern char *mpf_ILLraw_unique_name (ILLsymboltab * tab,
																 char *prefix,
																 int i);
extern int mpf_ILLraw_fill_in_rownames (mpf_rawlpdata * lp);

extern int mpf_ILLraw_init_bounds (mpf_rawlpdata * lp);

extern const char *mpf_ILLraw_set_lowerBound (mpf_rawlpdata * lp,
																					int i,
																					mpf_t bnd);
extern const char *mpf_ILLraw_set_upperBound (mpf_rawlpdata * lp,
																					int i,
																					mpf_t bnd);
extern const char *mpf_ILLraw_set_fixedBound (mpf_rawlpdata * lp,
																					int i,
																					mpf_t bnd);
extern const char *mpf_ILLraw_set_binaryBound (mpf_rawlpdata * lp,
																					 int i);
extern const char *mpf_ILLraw_set_unbound (mpf_rawlpdata * lp,
																			 int colind);
extern int mpf_ILLraw_fill_in_bounds (mpf_rawlpdata * lp);

extern int mpf_ILLraw_first_nondefault_bound (mpf_ILLlpdata * lp);
extern int mpf_ILLraw_default_lower (mpf_ILLlpdata * lp,
																 int i);
extern int mpf_ILLraw_default_upper (mpf_ILLlpdata * lp,
																 int i);

extern int mpf_ILLrawlpdata_to_lpdata (mpf_rawlpdata * raw,
																	 mpf_ILLlpdata * lp);

extern int mpf_ILLdata_error (mpf_qserror_collector * collector,
													const char *format,
													...);
extern void mpf_ILLdata_warn (mpf_qserror_collector * collector,
													const char *format,
													...);

#endif
