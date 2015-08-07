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

/*  $RCSfile: dbl_qsopt.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef dbl___QS_QSOPT_H
#define dbl___QS_QSOPT_H

#include <stdio.h>
#include "econfig.h"

#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define dbl_QSLIB_INTERFACE __declspec(dllexport)
#else
#define dbl_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define dbl_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct dbl_QSLIB_INTERFACE dbl_qsdata *dbl_QSprob;
typedef struct dbl_QSLIB_INTERFACE qsbasis *dbl_QSbas;
#else
typedef struct dbl_qsdata *dbl_QSprob;
typedef struct qsbasis *dbl_QSbas;
#endif

/****************************************************************************/
/*                                                                          */
/*                 PARAMETERS TO SPECIFY OBJECTIVE SENSE                    */
/*                                                                          */
/****************************************************************************/
#include "basicdefs.h"
/*
#define QS_LP_PRIMAL_FEASIBLE   11
#define QS_LP_PRIMAL_INFEASIBLE 12
#define QS_LP_PRIMAL_UNBOUNDED  13
#define QS_LP_DUAL_FEASIBLE     14
#define QS_LP_DUAL_INFEASIBLE   15
#define QS_LP_DUAL_UNBOUNDED    16
*/

/****************************************************************************/
/*                                                                          */
/*                      QSopt Library Functions                             */
/*                                                                          */
/****************************************************************************/
#ifdef  __cplusplus
extern "C"
{
#endif

#ifdef WIN32
/* 
 *  in WINDOWS we make 
 *     dbl_solver_main/dbl_reader_main part of DLL
 */
	dbl_QSLIB_INTERFACE int dbl_solver_main (int argc,
																	 char **argv);
	dbl_QSLIB_INTERFACE int dbl_reader_main (int argc,
																	 char **argv);
#endif

	dbl_QSLIB_INTERFACE void dbl_QSfree (void *ptr),
	  dbl_QSfree_prob (dbl_QSprob p),
	  dbl_QSfree_basis (dbl_QSbas B),
	  dbl_QSset_precision (const unsigned prec),/**< set the dbl_precision for floating 
																							 point numbers to the given 
																							 number of bits */

	  dbl_QSstart (void),	/**< whe we use non native numbers, we need to make 
												 some initializations before operating with the
												 library */

	  dbl_QSend (void);	/**< just to free any internal static data needed by
											 the variable dbl_precision numbers */

	dbl_QSLIB_INTERFACE int dbl_QSopt_primal (dbl_QSprob p,
																		int *status),
	  dbl_QSopt_dual (dbl_QSprob p,
								int *status),
	  dbl_QSopt_pivotin_col (dbl_QSprob p,
											 int ccnt,
											 int *clist),
	  dbl_QSopt_pivotin_row (dbl_QSprob p,
											 int rcnt,
											 int *rlist),
	  dbl_QSopt_strongbranch (dbl_QSprob p,
												int ncand,
												int *candidatelist,
												double * xlist,
												double * down_vals,
												double * up_vals,
												int iterations,
												double objbound),
	  dbl_QSchange_objsense (dbl_QSprob p,
											 int newsense),
	  dbl_QSget_objsense (dbl_QSprob p,
										int *newsense),
	  dbl_QSnew_col (dbl_QSprob p,
							 double obj,
							 double lower,
							 double upper,
							 const char *name),
	  dbl_QSadd_cols (dbl_QSprob p,
								int num,
								int *cmatcnt,
								int *cmatbeg,
								int *cmatind,
								double * cmatval,
								double * obj,
								double * lower,
								double * upper,
								const char **names),
	  dbl_QSadd_col (dbl_QSprob p,
							 int cnt,
							 int *cmatind,
							 double * cmatval,
							 double obj,
							 double lower,
							 double upper,
							 const char *name),
	  dbl_QSnew_row (dbl_QSprob p,
							 double rhs,
							 int sense,
							 const char *name),
	  dbl_QSadd_rows (dbl_QSprob p,
								int num,
								int *rmatcnt,
								int *rmatbeg,
								int *rmatind,
								double * rmatval,
								double * rhs,
								char *sense,
								const char **names),
	  dbl_QSadd_row (dbl_QSprob p,
							 int cnt,
							 int *rmatind,
							 double * rmatval,
							 double * rhs,
							 int sense,
							 const char *name),
	  dbl_QSdelete_rows (dbl_QSprob p,
									 int num,
									 int *dellist),
	  dbl_QSdelete_row (dbl_QSprob p,
									int rowindex),
	  dbl_QSdelete_setrows (dbl_QSprob p,
											int *flags),
	  dbl_QSdelete_named_row (dbl_QSprob p,
												const char *rowname),
	  dbl_QSdelete_named_rows_list (dbl_QSprob p,
															int num,
															const char **rownames),
	  dbl_QSdelete_cols (dbl_QSprob p,
									 int num,
									 int *dellist),
	  dbl_QSdelete_col (dbl_QSprob p,
									int colindex),
	  dbl_QSdelete_setcols (dbl_QSprob p,
											int *flags),
	  dbl_QSdelete_named_column (dbl_QSprob p,
													 const char *colname),
	  dbl_QSdelete_named_columns_list (dbl_QSprob p,
																 int num,
																 const char **colnames),
	  dbl_QSchange_senses (dbl_QSprob p,
										 int num,
										 int *rowlist,
										 char *sense),
	  dbl_QSchange_sense (dbl_QSprob p,
										int rowindex,
										int sense),
	  dbl_QSchange_coef (dbl_QSprob p,
									 int rowindex,
									 int colindex,
									 double coef),
	  dbl_QSchange_objcoef (dbl_QSprob p,
											int indx,
											double coef),
	  dbl_QSchange_rhscoef (dbl_QSprob p,
											int indx,
											double coef),
	  dbl_QSchange_bounds (dbl_QSprob p,
										 int num,
										 int *collist,
										 char *lu,
										 double * bounds),
	  dbl_QSchange_bound (dbl_QSprob p,
										int indx,
										int lu,
										double bound),
	  dbl_QSload_basis (dbl_QSprob p,
									dbl_QSbas B),
	  dbl_QSread_and_load_basis (dbl_QSprob p,
													 const char *filename),
	  dbl_QSload_basis_array (dbl_QSprob p,
												char *cstat,
												char *rstat),
	  dbl_QSload_basis_and_row_norms_array (dbl_QSprob p,
																			char *cstat,
																			char *rstat,
																			double * rownorms),
	  dbl_QSget_basis_array (dbl_QSprob p,
											 char *cstat,
											 char *rstat),
	  dbl_QSget_basis_and_row_norms_array (dbl_QSprob p,
																		 char *cstat,
																		 char *rstat,
																		 double * rownorms),
	  dbl_QSget_binv_row (dbl_QSprob p,
										int indx,
										double * binvrow),
	  dbl_QSget_tableau_row (dbl_QSprob p,
											 int indx,
											 double * tableaurow),
	  dbl_QSget_basis_order (dbl_QSprob p,
											 int *basorder),
	  dbl_QSget_status (dbl_QSprob p,
									int *status),
	  dbl_QSget_solution (dbl_QSprob p,
										double * value,
										double * x,
										double * pi,
										double * slack,
										double * rc),
	  dbl_QSget_objval (dbl_QSprob p,
									double * value),
	  dbl_QSget_pi_array (dbl_QSprob p,
										double * pi),
	  dbl_QSget_rc_array (dbl_QSprob p,
										double * rc),
	  dbl_QSget_x_array (dbl_QSprob p,
									 double * x),
	  dbl_QSget_slack_array (dbl_QSprob p,
											 double * slack),
	  dbl_QSget_infeas_array (dbl_QSprob p,
												double * pi),
	  dbl_QSget_colcount (dbl_QSprob p),
	  dbl_QSget_rowcount (dbl_QSprob p),
	  dbl_QSget_nzcount (dbl_QSprob p),
	  dbl_QSget_obj (dbl_QSprob p,
							 double * obj),
	  dbl_QSget_rhs (dbl_QSprob p,
							 double * rhs),
	  dbl_QSget_rows_list (dbl_QSprob p,
										 int num,
										 int *rowlist,
										 int **rowcnt,
										 int **rowbeg,
										 int **rowind,
										 double ** rowval,
										 double ** rhs,
										 char **sense,
										 char ***names),
	  dbl_QSget_rows (dbl_QSprob p,
								int **rowcnt,
								int **rowbeg,
								int **rowind,
								double ** rowval,
								double ** rhs,
								char **sense,
								char ***names),
	  dbl_QSget_columns_list (dbl_QSprob p,
												int num,
												int *collist,
												int **colcnt,
												int **colbeg,
												int **colind,
												double ** colval,
												double ** obj,
												double ** lower,
												double ** upper,
												char ***names),
	  dbl_QSget_columns (dbl_QSprob p,
									 int **colcnt,
									 int **colbeg,
									 int **colind,
									 double ** colval,
									 double ** obj,
									 double ** lower,
									 double ** upper,
									 char ***names),
	  dbl_QSget_rownames (dbl_QSprob p,
										char **rownames),
	  dbl_QSget_colnames (dbl_QSprob p,
										char **colnames),
	  dbl_QSget_bound (dbl_QSprob p,
								 int colindex,
								 int lu,
								 double * bound),
	  dbl_QSget_bounds (dbl_QSprob p,
									double * lower,
									double * upper),
	  dbl_QSget_intflags (dbl_QSprob p,
										int *intflags),
	  dbl_QSget_intcount (dbl_QSprob p,
										int *count),
	  dbl_QSget_column_index (dbl_QSprob p,
												const char *name,
												int *colindex),
	  dbl_QSget_row_index (dbl_QSprob p,
										 const char *name,
										 int *rowindex),
	  dbl_QSget_named_x (dbl_QSprob p,
									 const char *colname,
									 double * val),
	  dbl_QSget_named_rc (dbl_QSprob p,
										const char *colname,
										double * val),
	  dbl_QSget_named_pi (dbl_QSprob p,
										const char *rowname,
										double * val),
	  dbl_QSget_named_slack (dbl_QSprob p,
											 const char *rowname,
											 double * val),
	  dbl_QScompute_row_norms (dbl_QSprob p),
	  dbl_QSwrite_prob (dbl_QSprob p,
									const char *filename,
									const char *filetype),
	  dbl_QSwrite_prob_file (dbl_QSprob p,
											 FILE * file,
											 const char *filetype),
	  dbl_QSwrite_basis (dbl_QSprob p,
									 dbl_QSbas B,
									 const char *filename),
	  dbl_QStest_row_norms (dbl_QSprob p),
	  dbl_QSset_param (dbl_QSprob p,
								 int whichparam,
								 int newvalue),
	  dbl_QSset_param_EGlpNum (dbl_QSprob p,
												 int whichparam,
												 double newvalue),
	  dbl_QSget_param (dbl_QSprob p,
								 int whichparam,
								 int *value),
	  dbl_QSget_param_EGlpNum (dbl_QSprob p,
												 int whichparam,
												 double * value);

	dbl_QSLIB_INTERFACE char *dbl_QSget_probname (dbl_QSprob p);
	dbl_QSLIB_INTERFACE char *dbl_QSget_objname (dbl_QSprob p);
	dbl_QSLIB_INTERFACE char *dbl_QSversion (void);

	dbl_QSLIB_INTERFACE dbl_QSprob dbl_QScreate_prob (const char *name,
																				int objsense),
	  dbl_QSread_prob (const char *filename,
								 const char *filetype),
	  dbl_QSload_prob (const char *probname,
								 int ncols,
								 int nrows,
								 int *cmatcnt,
								 int *cmatbeg,
								 int *cmatind,
								 double * cmatval,
								 int objsense,
								 double * obj,
								 double * rhs,
								 char *sense,
								 double * lower,
								 double * upper,
								 const char **colnames,
								 const char **rownames),
	  dbl_QScopy_prob (dbl_QSprob p,
								 const char *newname);

	dbl_QSLIB_INTERFACE dbl_QSbas dbl_QSget_basis (dbl_QSprob p),
	  dbl_QSread_basis (dbl_QSprob p,
									const char *filename);

#ifdef  __cplusplus
}
#endif

/****************************************************************************
 *
 * This is the undocumented part of the QSlib interface 
 *
 ****************************************************************************/
/* 
 * functions to facilitate line by line reading from other sources than 
 * files from within MPS/LP parsers  
 * 
 * functions to facilitate the collection of error information instead of 
 * having the parsers print messages to stderr
 *                              by mps/lp format writers
 * 
 * a problem's reporter is used by the solver code to provide important 
 * dbl_feedback/progress information
 */

#ifdef WIN32
typedef struct dbl_QSLIB_INTERFACE dbl_qsline_reader *dbl_QSline_reader;
typedef struct dbl_QSLIB_INTERFACE dbl_qsformat_error *dbl_QSformat_error;
typedef struct dbl_QSLIB_INTERFACE dbl_qserror_collector *dbl_QSerror_collector;
typedef struct dbl_QSLIB_INTERFACE dbl_qserror_memory *dbl_QSerror_memory;
#else
typedef struct dbl_qsline_reader *dbl_QSline_reader;
typedef struct dbl_qsformat_error *dbl_QSformat_error;
typedef struct dbl_qserror_collector *dbl_QSerror_collector;
typedef struct dbl_qserror_memory *dbl_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
	dbl_QSLIB_INTERFACE const char *dbl_QSformat_error_type_string (int tp);

	dbl_QSLIB_INTERFACE int dbl_QSerror_get_type (dbl_QSformat_error error);
	dbl_QSLIB_INTERFACE const char *dbl_QSerror_get_desc (dbl_QSformat_error error);
	dbl_QSLIB_INTERFACE int dbl_QSerror_get_line_number (dbl_QSformat_error error);
	dbl_QSLIB_INTERFACE int dbl_QSerror_get_pos (dbl_QSformat_error error);
	dbl_QSLIB_INTERFACE const char *dbl_QSerror_get_line (dbl_QSformat_error error);
	dbl_QSLIB_INTERFACE void dbl_QSerror_print (FILE * f,
																			dbl_QSformat_error error);

	dbl_QSLIB_INTERFACE dbl_QSerror_collector dbl_QSerror_collector_new (void *fct,
																													 void *dest);
	  dbl_QSLIB_INTERFACE
		dbl_QSerror_collector dbl_QSerror_memory_collector_new (dbl_QSerror_memory mem);
	dbl_QSLIB_INTERFACE void dbl_QSerror_collector_free (dbl_QSerror_collector c);

/****************************************************************************
 * line reader 
 */
	dbl_QSLIB_INTERFACE dbl_QSline_reader dbl_QSline_reader_new (void *fct,
																									 void *data_src);
	/* reader->read_line_fct defaults to fgets */

	dbl_QSLIB_INTERFACE void dbl_QSline_reader_free (dbl_QSline_reader reader);

	dbl_QSLIB_INTERFACE void dbl_QSline_reader_set_error_collector (dbl_QSline_reader reader,
																													dbl_QSerror_collector
																													collector);

	dbl_QSLIB_INTERFACE char *dbl_QSline_reader_get (dbl_QSline_reader reader,
																					 char *s,
																					 int size);

	dbl_QSLIB_INTERFACE dbl_QSprob dbl_QSget_prob (dbl_QSline_reader reader,
																		 const char *probname,
																		 const char *filetype);
	/* the MPS and LP parsers uses the fct from reader 
	 * to get to next input line */


/****************************************************************************
 * error memory 
 */
	dbl_QSLIB_INTERFACE dbl_QSerror_memory dbl_QSerror_memory_create (int takeErrorLines);
	dbl_QSLIB_INTERFACE void dbl_QSerror_memory_free (dbl_QSerror_memory mem);

	dbl_QSLIB_INTERFACE int dbl_QSerror_memory_get_nof (dbl_QSerror_memory mem,
																							int error_type);
	dbl_QSLIB_INTERFACE int dbl_QSerror_memory_get_nerrors (dbl_QSerror_memory mem);

	dbl_QSLIB_INTERFACE dbl_QSformat_error
		dbl_QSerror_memory_get_last_error (dbl_QSerror_memory mem);
	dbl_QSLIB_INTERFACE dbl_QSformat_error
		dbl_QSerror_memory_get_prev_error (dbl_QSformat_error e);

/**************************************************************************** 
 * reporter for solver dbl_feedback 
 */
	dbl_QSLIB_INTERFACE void dbl_QSset_reporter (dbl_QSprob prob,
																			 int iterskip,
																			 void *fct,
																			 void *dest);

	dbl_QSLIB_INTERFACE int dbl_QSreport_prob (dbl_QSprob p,
																		 const char *filetype,
																		 dbl_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif													/* dbl___QS_QSOPT_H */
