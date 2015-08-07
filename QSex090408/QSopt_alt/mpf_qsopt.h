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

/*  $RCSfile: mpf_qsopt.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef mpf___QS_QSOPT_H
#define mpf___QS_QSOPT_H

#include <stdio.h>
#include "econfig.h"

#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define mpf_QSLIB_INTERFACE __declspec(dllexport)
#else
#define mpf_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define mpf_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct mpf_QSLIB_INTERFACE mpf_qsdata *mpf_QSprob;
typedef struct mpf_QSLIB_INTERFACE qsbasis *mpf_QSbas;
#else
typedef struct mpf_qsdata *mpf_QSprob;
typedef struct qsbasis *mpf_QSbas;
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
 *     mpf_solver_main/mpf_reader_main part of DLL
 */
	mpf_QSLIB_INTERFACE int mpf_solver_main (int argc,
																	 char **argv);
	mpf_QSLIB_INTERFACE int mpf_reader_main (int argc,
																	 char **argv);
#endif

	mpf_QSLIB_INTERFACE void mpf_QSfree (void *ptr),
	  mpf_QSfree_prob (mpf_QSprob p),
	  mpf_QSfree_basis (mpf_QSbas B),
	  mpf_QSset_precision (const unsigned prec),/**< set the mpf_precision for floating 
																							 point numbers to the given 
																							 number of bits */

	  mpf_QSstart (void),	/**< whe we use non native numbers, we need to make 
												 some initializations before operating with the
												 library */

	  mpf_QSend (void);	/**< just to free any internal static data needed by
											 the variable mpf_precision numbers */

	mpf_QSLIB_INTERFACE int mpf_QSopt_primal (mpf_QSprob p,
																		int *status),
	  mpf_QSopt_dual (mpf_QSprob p,
								int *status),
	  mpf_QSopt_pivotin_col (mpf_QSprob p,
											 int ccnt,
											 int *clist),
	  mpf_QSopt_pivotin_row (mpf_QSprob p,
											 int rcnt,
											 int *rlist),
	  mpf_QSopt_strongbranch (mpf_QSprob p,
												int ncand,
												int *candidatelist,
												mpf_t * xlist,
												mpf_t * down_vals,
												mpf_t * up_vals,
												int iterations,
												mpf_t objbound),
	  mpf_QSchange_objsense (mpf_QSprob p,
											 int newsense),
	  mpf_QSget_objsense (mpf_QSprob p,
										int *newsense),
	  mpf_QSnew_col (mpf_QSprob p,
							 mpf_t obj,
							 mpf_t lower,
							 mpf_t upper,
							 const char *name),
	  mpf_QSadd_cols (mpf_QSprob p,
								int num,
								int *cmatcnt,
								int *cmatbeg,
								int *cmatind,
								mpf_t * cmatval,
								mpf_t * obj,
								mpf_t * lower,
								mpf_t * upper,
								const char **names),
	  mpf_QSadd_col (mpf_QSprob p,
							 int cnt,
							 int *cmatind,
							 mpf_t * cmatval,
							 mpf_t obj,
							 mpf_t lower,
							 mpf_t upper,
							 const char *name),
	  mpf_QSnew_row (mpf_QSprob p,
							 mpf_t rhs,
							 int sense,
							 const char *name),
	  mpf_QSadd_rows (mpf_QSprob p,
								int num,
								int *rmatcnt,
								int *rmatbeg,
								int *rmatind,
								mpf_t * rmatval,
								mpf_t * rhs,
								char *sense,
								const char **names),
	  mpf_QSadd_row (mpf_QSprob p,
							 int cnt,
							 int *rmatind,
							 mpf_t * rmatval,
							 mpf_t * rhs,
							 int sense,
							 const char *name),
	  mpf_QSdelete_rows (mpf_QSprob p,
									 int num,
									 int *dellist),
	  mpf_QSdelete_row (mpf_QSprob p,
									int rowindex),
	  mpf_QSdelete_setrows (mpf_QSprob p,
											int *flags),
	  mpf_QSdelete_named_row (mpf_QSprob p,
												const char *rowname),
	  mpf_QSdelete_named_rows_list (mpf_QSprob p,
															int num,
															const char **rownames),
	  mpf_QSdelete_cols (mpf_QSprob p,
									 int num,
									 int *dellist),
	  mpf_QSdelete_col (mpf_QSprob p,
									int colindex),
	  mpf_QSdelete_setcols (mpf_QSprob p,
											int *flags),
	  mpf_QSdelete_named_column (mpf_QSprob p,
													 const char *colname),
	  mpf_QSdelete_named_columns_list (mpf_QSprob p,
																 int num,
																 const char **colnames),
	  mpf_QSchange_senses (mpf_QSprob p,
										 int num,
										 int *rowlist,
										 char *sense),
	  mpf_QSchange_sense (mpf_QSprob p,
										int rowindex,
										int sense),
	  mpf_QSchange_coef (mpf_QSprob p,
									 int rowindex,
									 int colindex,
									 mpf_t coef),
	  mpf_QSchange_objcoef (mpf_QSprob p,
											int indx,
											mpf_t coef),
	  mpf_QSchange_rhscoef (mpf_QSprob p,
											int indx,
											mpf_t coef),
	  mpf_QSchange_bounds (mpf_QSprob p,
										 int num,
										 int *collist,
										 char *lu,
										 mpf_t * bounds),
	  mpf_QSchange_bound (mpf_QSprob p,
										int indx,
										int lu,
										mpf_t bound),
	  mpf_QSload_basis (mpf_QSprob p,
									mpf_QSbas B),
	  mpf_QSread_and_load_basis (mpf_QSprob p,
													 const char *filename),
	  mpf_QSload_basis_array (mpf_QSprob p,
												char *cstat,
												char *rstat),
	  mpf_QSload_basis_and_row_norms_array (mpf_QSprob p,
																			char *cstat,
																			char *rstat,
																			mpf_t * rownorms),
	  mpf_QSget_basis_array (mpf_QSprob p,
											 char *cstat,
											 char *rstat),
	  mpf_QSget_basis_and_row_norms_array (mpf_QSprob p,
																		 char *cstat,
																		 char *rstat,
																		 mpf_t * rownorms),
	  mpf_QSget_binv_row (mpf_QSprob p,
										int indx,
										mpf_t * binvrow),
	  mpf_QSget_tableau_row (mpf_QSprob p,
											 int indx,
											 mpf_t * tableaurow),
	  mpf_QSget_basis_order (mpf_QSprob p,
											 int *basorder),
	  mpf_QSget_status (mpf_QSprob p,
									int *status),
	  mpf_QSget_solution (mpf_QSprob p,
										mpf_t * value,
										mpf_t * x,
										mpf_t * pi,
										mpf_t * slack,
										mpf_t * rc),
	  mpf_QSget_objval (mpf_QSprob p,
									mpf_t * value),
	  mpf_QSget_pi_array (mpf_QSprob p,
										mpf_t * pi),
	  mpf_QSget_rc_array (mpf_QSprob p,
										mpf_t * rc),
	  mpf_QSget_x_array (mpf_QSprob p,
									 mpf_t * x),
	  mpf_QSget_slack_array (mpf_QSprob p,
											 mpf_t * slack),
	  mpf_QSget_infeas_array (mpf_QSprob p,
												mpf_t * pi),
	  mpf_QSget_colcount (mpf_QSprob p),
	  mpf_QSget_rowcount (mpf_QSprob p),
	  mpf_QSget_nzcount (mpf_QSprob p),
	  mpf_QSget_obj (mpf_QSprob p,
							 mpf_t * obj),
	  mpf_QSget_rhs (mpf_QSprob p,
							 mpf_t * rhs),
	  mpf_QSget_rows_list (mpf_QSprob p,
										 int num,
										 int *rowlist,
										 int **rowcnt,
										 int **rowbeg,
										 int **rowind,
										 mpf_t ** rowval,
										 mpf_t ** rhs,
										 char **sense,
										 char ***names),
	  mpf_QSget_rows (mpf_QSprob p,
								int **rowcnt,
								int **rowbeg,
								int **rowind,
								mpf_t ** rowval,
								mpf_t ** rhs,
								char **sense,
								char ***names),
	  mpf_QSget_columns_list (mpf_QSprob p,
												int num,
												int *collist,
												int **colcnt,
												int **colbeg,
												int **colind,
												mpf_t ** colval,
												mpf_t ** obj,
												mpf_t ** lower,
												mpf_t ** upper,
												char ***names),
	  mpf_QSget_columns (mpf_QSprob p,
									 int **colcnt,
									 int **colbeg,
									 int **colind,
									 mpf_t ** colval,
									 mpf_t ** obj,
									 mpf_t ** lower,
									 mpf_t ** upper,
									 char ***names),
	  mpf_QSget_rownames (mpf_QSprob p,
										char **rownames),
	  mpf_QSget_colnames (mpf_QSprob p,
										char **colnames),
	  mpf_QSget_bound (mpf_QSprob p,
								 int colindex,
								 int lu,
								 mpf_t * bound),
	  mpf_QSget_bounds (mpf_QSprob p,
									mpf_t * lower,
									mpf_t * upper),
	  mpf_QSget_intflags (mpf_QSprob p,
										int *intflags),
	  mpf_QSget_intcount (mpf_QSprob p,
										int *count),
	  mpf_QSget_column_index (mpf_QSprob p,
												const char *name,
												int *colindex),
	  mpf_QSget_row_index (mpf_QSprob p,
										 const char *name,
										 int *rowindex),
	  mpf_QSget_named_x (mpf_QSprob p,
									 const char *colname,
									 mpf_t * val),
	  mpf_QSget_named_rc (mpf_QSprob p,
										const char *colname,
										mpf_t * val),
	  mpf_QSget_named_pi (mpf_QSprob p,
										const char *rowname,
										mpf_t * val),
	  mpf_QSget_named_slack (mpf_QSprob p,
											 const char *rowname,
											 mpf_t * val),
	  mpf_QScompute_row_norms (mpf_QSprob p),
	  mpf_QSwrite_prob (mpf_QSprob p,
									const char *filename,
									const char *filetype),
	  mpf_QSwrite_prob_file (mpf_QSprob p,
											 FILE * file,
											 const char *filetype),
	  mpf_QSwrite_basis (mpf_QSprob p,
									 mpf_QSbas B,
									 const char *filename),
	  mpf_QStest_row_norms (mpf_QSprob p),
	  mpf_QSset_param (mpf_QSprob p,
								 int whichparam,
								 int newvalue),
	  mpf_QSset_param_EGlpNum (mpf_QSprob p,
												 int whichparam,
												 mpf_t newvalue),
	  mpf_QSget_param (mpf_QSprob p,
								 int whichparam,
								 int *value),
	  mpf_QSget_param_EGlpNum (mpf_QSprob p,
												 int whichparam,
												 mpf_t * value);

	mpf_QSLIB_INTERFACE char *mpf_QSget_probname (mpf_QSprob p);
	mpf_QSLIB_INTERFACE char *mpf_QSget_objname (mpf_QSprob p);
	mpf_QSLIB_INTERFACE char *mpf_QSversion (void);

	mpf_QSLIB_INTERFACE mpf_QSprob mpf_QScreate_prob (const char *name,
																				int objsense),
	  mpf_QSread_prob (const char *filename,
								 const char *filetype),
	  mpf_QSload_prob (const char *probname,
								 int ncols,
								 int nrows,
								 int *cmatcnt,
								 int *cmatbeg,
								 int *cmatind,
								 mpf_t * cmatval,
								 int objsense,
								 mpf_t * obj,
								 mpf_t * rhs,
								 char *sense,
								 mpf_t * lower,
								 mpf_t * upper,
								 const char **colnames,
								 const char **rownames),
	  mpf_QScopy_prob (mpf_QSprob p,
								 const char *newname);

	mpf_QSLIB_INTERFACE mpf_QSbas mpf_QSget_basis (mpf_QSprob p),
	  mpf_QSread_basis (mpf_QSprob p,
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
 * mpf_feedback/progress information
 */

#ifdef WIN32
typedef struct mpf_QSLIB_INTERFACE mpf_qsline_reader *mpf_QSline_reader;
typedef struct mpf_QSLIB_INTERFACE mpf_qsformat_error *mpf_QSformat_error;
typedef struct mpf_QSLIB_INTERFACE mpf_qserror_collector *mpf_QSerror_collector;
typedef struct mpf_QSLIB_INTERFACE mpf_qserror_memory *mpf_QSerror_memory;
#else
typedef struct mpf_qsline_reader *mpf_QSline_reader;
typedef struct mpf_qsformat_error *mpf_QSformat_error;
typedef struct mpf_qserror_collector *mpf_QSerror_collector;
typedef struct mpf_qserror_memory *mpf_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
	mpf_QSLIB_INTERFACE const char *mpf_QSformat_error_type_string (int tp);

	mpf_QSLIB_INTERFACE int mpf_QSerror_get_type (mpf_QSformat_error error);
	mpf_QSLIB_INTERFACE const char *mpf_QSerror_get_desc (mpf_QSformat_error error);
	mpf_QSLIB_INTERFACE int mpf_QSerror_get_line_number (mpf_QSformat_error error);
	mpf_QSLIB_INTERFACE int mpf_QSerror_get_pos (mpf_QSformat_error error);
	mpf_QSLIB_INTERFACE const char *mpf_QSerror_get_line (mpf_QSformat_error error);
	mpf_QSLIB_INTERFACE void mpf_QSerror_print (FILE * f,
																			mpf_QSformat_error error);

	mpf_QSLIB_INTERFACE mpf_QSerror_collector mpf_QSerror_collector_new (void *fct,
																													 void *dest);
	  mpf_QSLIB_INTERFACE
		mpf_QSerror_collector mpf_QSerror_memory_collector_new (mpf_QSerror_memory mem);
	mpf_QSLIB_INTERFACE void mpf_QSerror_collector_free (mpf_QSerror_collector c);

/****************************************************************************
 * line reader 
 */
	mpf_QSLIB_INTERFACE mpf_QSline_reader mpf_QSline_reader_new (void *fct,
																									 void *data_src);
	/* reader->read_line_fct defaults to fgets */

	mpf_QSLIB_INTERFACE void mpf_QSline_reader_free (mpf_QSline_reader reader);

	mpf_QSLIB_INTERFACE void mpf_QSline_reader_set_error_collector (mpf_QSline_reader reader,
																													mpf_QSerror_collector
																													collector);

	mpf_QSLIB_INTERFACE char *mpf_QSline_reader_get (mpf_QSline_reader reader,
																					 char *s,
																					 int size);

	mpf_QSLIB_INTERFACE mpf_QSprob mpf_QSget_prob (mpf_QSline_reader reader,
																		 const char *probname,
																		 const char *filetype);
	/* the MPS and LP parsers uses the fct from reader 
	 * to get to next input line */


/****************************************************************************
 * error memory 
 */
	mpf_QSLIB_INTERFACE mpf_QSerror_memory mpf_QSerror_memory_create (int takeErrorLines);
	mpf_QSLIB_INTERFACE void mpf_QSerror_memory_free (mpf_QSerror_memory mem);

	mpf_QSLIB_INTERFACE int mpf_QSerror_memory_get_nof (mpf_QSerror_memory mem,
																							int error_type);
	mpf_QSLIB_INTERFACE int mpf_QSerror_memory_get_nerrors (mpf_QSerror_memory mem);

	mpf_QSLIB_INTERFACE mpf_QSformat_error
		mpf_QSerror_memory_get_last_error (mpf_QSerror_memory mem);
	mpf_QSLIB_INTERFACE mpf_QSformat_error
		mpf_QSerror_memory_get_prev_error (mpf_QSformat_error e);

/**************************************************************************** 
 * reporter for solver mpf_feedback 
 */
	mpf_QSLIB_INTERFACE void mpf_QSset_reporter (mpf_QSprob prob,
																			 int iterskip,
																			 void *fct,
																			 void *dest);

	mpf_QSLIB_INTERFACE int mpf_QSreport_prob (mpf_QSprob p,
																		 const char *filetype,
																		 mpf_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif													/* mpf___QS_QSOPT_H */
