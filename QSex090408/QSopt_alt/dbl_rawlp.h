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

/*  RCS_INFO = "$RCSfile: dbl_rawlp.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef dbl___ILL_RAWLP_H_
#define dbl___ILL_RAWLP_H_

/****************************************************************************/
/* DataStructure and Routines                                               */
/*          to deal with raw lp information as read from mps or lp files    */
/*          support scanning of input                                       */
/*          error reporting                                                 */
/****************************************************************************/

#include "trace.h"
#include "dbl_lpdata.h"
#include "dbl_iqsutil.h"
#include "dbl_format.h"
#include "dbl_lpdefs.h"

#define dbl_ILL_ISBLANK(p) \
             (((*(p))==' '||(*(p))=='\t'||(*(p))=='\r'||(*(p))=='\f') ? 1 : 0)

/* 
 * we rely on ILLsymboltab property:
 *   the ith name added can be retrieved by ILLsymboltab_get(table, i) 
 *   as long as we never delete names from the symbol table 
 */
typedef struct dbl_rawlpdata
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
	double *rhs;								/* rhs values for rows; size is nrows */
	char *rangesind;							/* ranges[i] == 1 we saw a range def for row[i] */
	struct dbl_colptr *ranges;				/* list of range values */

	int ncols;										/* number of cols in problem */
	ILLsymboltab coltab;					/* ILLsymboltab_get(coltab, i) name of ith col */
	int colsize;									/* size of cols array */
	struct dbl_colptr **cols;

	char *lbind;									/* lbind[i] == 1  we saw a lower bound for col[i] */
	char *ubind;									/* ubind[i] == 1  we saw a upper bound for col[i] */
	double *lower;							/* lower[i] = lower bound for col[i] */
	double *upper;							/* upper[i] = upper bound for col[i] */

	int intsize;									/* size of intmarker array */
	char *intmarker;							/* intmarker[i] == 1  col[i] is an int var */

	/* sos information is tranfered into dbl_ILLmatrix lpdata->sos */
	char *refrow;									/* name of reference row */
	int refrowind;								/* index of refrow or -1  */

	int is_sos_size;							/* size of is_sos_member array */
	int *is_sos_member;						/* for each col contains either               
																 *     -1 == no sos memeber 
																 *     i  == member of set #i */

	int nsos_member;							/* total number of sos set members */
	int sos_weight_size;					/* size of sos_weight array */
	double *sos_weight;				/* sos set elem i has weight of sos_weight[i] 
																 * value comes from refrow coeficients */
	int sos_col_size;							/* size of sos_col array */
	int *sos_col;									/* sos elem i is column sos_col[i] */

	int nsos;											/* number of sos sets */
	int sos_setsize;							/* size of sosset array */
	struct dbl_sosptr *sos_set;				/* type, size, first element of sos sets 
																 * first is index into sos_weight and sos_col 
																 * arrays */
	dbl_qserror_collector *error_collector;
	ILLptrworld ptrworld;
}
dbl_rawlpdata;

typedef struct dbl_colptr
{
	double coef;
	struct dbl_colptr *next;
	int this;											/* row index */
}
dbl_colptr;
extern dbl_colptr *dbl_ILLcolptralloc (ILLptrworld * p);

typedef struct dbl_sosptr
{
	int nelem;										/* number of set elements */
	int first;										/* index of first set element in sosmemeber */
	char type;										/* set type */
}
dbl_sosptr;
extern const int dbl_ILL_SOS_TYPE1;
extern const int dbl_ILL_SOS_TYPE2;

extern void dbl_ILLinit_rawlpdata (dbl_rawlpdata * lp,
															 dbl_qserror_collector * collector);
extern void dbl_ILLfree_rawlpdata (dbl_rawlpdata * lp);
extern void dbl_ILLraw_clear_matrix (dbl_rawlpdata * lp);

extern const char *dbl_ILLraw_rowname (dbl_rawlpdata * lp,
																	 int i);
extern const char *dbl_ILLraw_colname (dbl_rawlpdata * lp,
																	 int i);

extern int dbl_ILLraw_add_col (dbl_rawlpdata * lp,
													 const char *name,
													 int intmarker);
extern int dbl_ILLraw_add_row (dbl_rawlpdata * lp,
													 const char *name,
													 int sense,
													 double rhs);

extern int dbl_ILLraw_add_col_coef (dbl_rawlpdata * lp,
																int colind,
																int rowind,
																double coef);

extern int dbl_ILLraw_init_ranges (dbl_rawlpdata * lp);
extern int dbl_ILLraw_init_rhs (dbl_rawlpdata * lp);

extern int dbl_ILLraw_add_ranges_coef (dbl_rawlpdata * lp,
																	 int rowind,
																	 double coef);


extern int dbl_ILLraw_add_sos (dbl_rawlpdata * lp,
													 int sos_type);
																								/* add empty set with type */
extern int dbl_ILLraw_add_sos_member (dbl_rawlpdata * lp,
																	int colind);
																								/* add col to last set */
extern int dbl_ILLraw_is_mem_other_sos (dbl_rawlpdata * lp,
																		int colind);

extern int dbl_ILLraw_set_rhs_name (dbl_rawlpdata * lp,
																const char *name,
																int *skip);
extern int dbl_ILLraw_set_bounds_name (dbl_rawlpdata * lp,
																	 const char *name,
																	 int *skip);
extern int dbl_ILLraw_set_ranges_name (dbl_rawlpdata * lp,
																	 const char *name,
																	 int *skip);
extern void dbl_ILLprint_rawlpdata (dbl_rawlpdata * lp);

extern char *dbl_ILLraw_unique_name (ILLsymboltab * tab,
																 char *prefix,
																 int i);
extern int dbl_ILLraw_fill_in_rownames (dbl_rawlpdata * lp);

extern int dbl_ILLraw_init_bounds (dbl_rawlpdata * lp);

extern const char *dbl_ILLraw_set_lowerBound (dbl_rawlpdata * lp,
																					int i,
																					double bnd);
extern const char *dbl_ILLraw_set_upperBound (dbl_rawlpdata * lp,
																					int i,
																					double bnd);
extern const char *dbl_ILLraw_set_fixedBound (dbl_rawlpdata * lp,
																					int i,
																					double bnd);
extern const char *dbl_ILLraw_set_binaryBound (dbl_rawlpdata * lp,
																					 int i);
extern const char *dbl_ILLraw_set_unbound (dbl_rawlpdata * lp,
																			 int colind);
extern int dbl_ILLraw_fill_in_bounds (dbl_rawlpdata * lp);

extern int dbl_ILLraw_first_nondefault_bound (dbl_ILLlpdata * lp);
extern int dbl_ILLraw_default_lower (dbl_ILLlpdata * lp,
																 int i);
extern int dbl_ILLraw_default_upper (dbl_ILLlpdata * lp,
																 int i);

extern int dbl_ILLrawlpdata_to_lpdata (dbl_rawlpdata * raw,
																	 dbl_ILLlpdata * lp);

extern int dbl_ILLdata_error (dbl_qserror_collector * collector,
													const char *format,
													...);
extern void dbl_ILLdata_warn (dbl_qserror_collector * collector,
													const char *format,
													...);

#endif
