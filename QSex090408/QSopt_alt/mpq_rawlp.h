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

/*  RCS_INFO = "$RCSfile: mpq_rawlp.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef mpq___ILL_RAWLP_H_
#define mpq___ILL_RAWLP_H_

/****************************************************************************/
/* DataStructure and Routines                                               */
/*          to deal with raw lp information as read from mps or lp files    */
/*          support scanning of input                                       */
/*          error reporting                                                 */
/****************************************************************************/

#include "trace.h"
#include "mpq_lpdata.h"
#include "mpq_iqsutil.h"
#include "mpq_format.h"
#include "mpq_lpdefs.h"

#define mpq_ILL_ISBLANK(p) \
             (((*(p))==' '||(*(p))=='\t'||(*(p))=='\r'||(*(p))=='\f') ? 1 : 0)

/* 
 * we rely on ILLsymboltab property:
 *   the ith name added can be retrieved by ILLsymboltab_get(table, i) 
 *   as long as we never delete names from the symbol table 
 */
typedef struct mpq_rawlpdata
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
	mpq_t *rhs;								/* rhs values for rows; size is nrows */
	char *rangesind;							/* ranges[i] == 1 we saw a range def for row[i] */
	struct mpq_colptr *ranges;				/* list of range values */

	int ncols;										/* number of cols in problem */
	ILLsymboltab coltab;					/* ILLsymboltab_get(coltab, i) name of ith col */
	int colsize;									/* size of cols array */
	struct mpq_colptr **cols;

	char *lbind;									/* lbind[i] == 1  we saw a lower bound for col[i] */
	char *ubind;									/* ubind[i] == 1  we saw a upper bound for col[i] */
	mpq_t *lower;							/* lower[i] = lower bound for col[i] */
	mpq_t *upper;							/* upper[i] = upper bound for col[i] */

	int intsize;									/* size of intmarker array */
	char *intmarker;							/* intmarker[i] == 1  col[i] is an int var */

	/* sos information is tranfered into mpq_ILLmatrix lpdata->sos */
	char *refrow;									/* name of reference row */
	int refrowind;								/* index of refrow or -1  */

	int is_sos_size;							/* size of is_sos_member array */
	int *is_sos_member;						/* for each col contains either               
																 *     -1 == no sos memeber 
																 *     i  == member of set #i */

	int nsos_member;							/* total number of sos set members */
	int sos_weight_size;					/* size of sos_weight array */
	mpq_t *sos_weight;				/* sos set elem i has weight of sos_weight[i] 
																 * value comes from refrow coeficients */
	int sos_col_size;							/* size of sos_col array */
	int *sos_col;									/* sos elem i is column sos_col[i] */

	int nsos;											/* number of sos sets */
	int sos_setsize;							/* size of sosset array */
	struct mpq_sosptr *sos_set;				/* type, size, first element of sos sets 
																 * first is index into sos_weight and sos_col 
																 * arrays */
	mpq_qserror_collector *error_collector;
	ILLptrworld ptrworld;
}
mpq_rawlpdata;

typedef struct mpq_colptr
{
	mpq_t coef;
	struct mpq_colptr *next;
	int this;											/* row index */
}
mpq_colptr;
extern mpq_colptr *mpq_ILLcolptralloc (ILLptrworld * p);

typedef struct mpq_sosptr
{
	int nelem;										/* number of set elements */
	int first;										/* index of first set element in sosmemeber */
	char type;										/* set type */
}
mpq_sosptr;
extern const int mpq_ILL_SOS_TYPE1;
extern const int mpq_ILL_SOS_TYPE2;

extern void mpq_ILLinit_rawlpdata (mpq_rawlpdata * lp,
															 mpq_qserror_collector * collector);
extern void mpq_ILLfree_rawlpdata (mpq_rawlpdata * lp);
extern void mpq_ILLraw_clear_matrix (mpq_rawlpdata * lp);

extern const char *mpq_ILLraw_rowname (mpq_rawlpdata * lp,
																	 int i);
extern const char *mpq_ILLraw_colname (mpq_rawlpdata * lp,
																	 int i);

extern int mpq_ILLraw_add_col (mpq_rawlpdata * lp,
													 const char *name,
													 int intmarker);
extern int mpq_ILLraw_add_row (mpq_rawlpdata * lp,
													 const char *name,
													 int sense,
													 mpq_t rhs);

extern int mpq_ILLraw_add_col_coef (mpq_rawlpdata * lp,
																int colind,
																int rowind,
																mpq_t coef);

extern int mpq_ILLraw_init_ranges (mpq_rawlpdata * lp);
extern int mpq_ILLraw_init_rhs (mpq_rawlpdata * lp);

extern int mpq_ILLraw_add_ranges_coef (mpq_rawlpdata * lp,
																	 int rowind,
																	 mpq_t coef);


extern int mpq_ILLraw_add_sos (mpq_rawlpdata * lp,
													 int sos_type);
																								/* add empty set with type */
extern int mpq_ILLraw_add_sos_member (mpq_rawlpdata * lp,
																	int colind);
																								/* add col to last set */
extern int mpq_ILLraw_is_mem_other_sos (mpq_rawlpdata * lp,
																		int colind);

extern int mpq_ILLraw_set_rhs_name (mpq_rawlpdata * lp,
																const char *name,
																int *skip);
extern int mpq_ILLraw_set_bounds_name (mpq_rawlpdata * lp,
																	 const char *name,
																	 int *skip);
extern int mpq_ILLraw_set_ranges_name (mpq_rawlpdata * lp,
																	 const char *name,
																	 int *skip);
extern void mpq_ILLprint_rawlpdata (mpq_rawlpdata * lp);

extern char *mpq_ILLraw_unique_name (ILLsymboltab * tab,
																 char *prefix,
																 int i);
extern int mpq_ILLraw_fill_in_rownames (mpq_rawlpdata * lp);

extern int mpq_ILLraw_init_bounds (mpq_rawlpdata * lp);

extern const char *mpq_ILLraw_set_lowerBound (mpq_rawlpdata * lp,
																					int i,
																					mpq_t bnd);
extern const char *mpq_ILLraw_set_upperBound (mpq_rawlpdata * lp,
																					int i,
																					mpq_t bnd);
extern const char *mpq_ILLraw_set_fixedBound (mpq_rawlpdata * lp,
																					int i,
																					mpq_t bnd);
extern const char *mpq_ILLraw_set_binaryBound (mpq_rawlpdata * lp,
																					 int i);
extern const char *mpq_ILLraw_set_unbound (mpq_rawlpdata * lp,
																			 int colind);
extern int mpq_ILLraw_fill_in_bounds (mpq_rawlpdata * lp);

extern int mpq_ILLraw_first_nondefault_bound (mpq_ILLlpdata * lp);
extern int mpq_ILLraw_default_lower (mpq_ILLlpdata * lp,
																 int i);
extern int mpq_ILLraw_default_upper (mpq_ILLlpdata * lp,
																 int i);

extern int mpq_ILLrawlpdata_to_lpdata (mpq_rawlpdata * raw,
																	 mpq_ILLlpdata * lp);

extern int mpq_ILLdata_error (mpq_qserror_collector * collector,
													const char *format,
													...);
extern void mpq_ILLdata_warn (mpq_qserror_collector * collector,
													const char *format,
													...);

#endif
