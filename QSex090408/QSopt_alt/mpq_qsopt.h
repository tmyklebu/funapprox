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

/*  $RCSfile: mpq_qsopt.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef mpq___QS_QSOPT_H
#define mpq___QS_QSOPT_H

#include <stdio.h>
#include "econfig.h"

#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define mpq_QSLIB_INTERFACE __declspec(dllexport)
#else
#define mpq_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define mpq_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct mpq_QSLIB_INTERFACE mpq_qsdata *mpq_QSprob;
typedef struct mpq_QSLIB_INTERFACE qsbasis *mpq_QSbas;
#else
typedef struct mpq_qsdata *mpq_QSprob;
typedef struct qsbasis *mpq_QSbas;
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
 *     mpq_solver_main/mpq_reader_main part of DLL
 */
	mpq_QSLIB_INTERFACE int mpq_solver_main (int argc,
																	 char **argv);
	mpq_QSLIB_INTERFACE int mpq_reader_main (int argc,
																	 char **argv);
#endif

	mpq_QSLIB_INTERFACE void mpq_QSfree (void *ptr),
	  mpq_QSfree_prob (mpq_QSprob p),
	  mpq_QSfree_basis (mpq_QSbas B),
	  mpq_QSset_precision (const unsigned prec),/**< set the mpq_precision for floating 
																							 point numbers to the given 
																							 number of bits */

	  mpq_QSstart (void),	/**< whe we use non native numbers, we need to make 
												 some initializations before operating with the
												 library */

	  mpq_QSend (void);	/**< just to free any internal static data needed by
											 the variable mpq_precision numbers */

	mpq_QSLIB_INTERFACE int mpq_QSopt_primal (mpq_QSprob p,
																		int *status),
	  mpq_QSopt_dual (mpq_QSprob p,
								int *status),
	  mpq_QSopt_pivotin_col (mpq_QSprob p,
											 int ccnt,
											 int *clist),
	  mpq_QSopt_pivotin_row (mpq_QSprob p,
											 int rcnt,
											 int *rlist),
	  mpq_QSopt_strongbranch (mpq_QSprob p,
												int ncand,
												int *candidatelist,
												mpq_t * xlist,
												mpq_t * down_vals,
												mpq_t * up_vals,
												int iterations,
												mpq_t objbound),
	  mpq_QSchange_objsense (mpq_QSprob p,
											 int newsense),
	  mpq_QSget_objsense (mpq_QSprob p,
										int *newsense),
	  mpq_QSnew_col (mpq_QSprob p,
							 mpq_t obj,
							 mpq_t lower,
							 mpq_t upper,
							 const char *name),
	  mpq_QSadd_cols (mpq_QSprob p,
								int num,
								int *cmatcnt,
								int *cmatbeg,
								int *cmatind,
								mpq_t * cmatval,
								mpq_t * obj,
								mpq_t * lower,
								mpq_t * upper,
								const char **names),
	  mpq_QSadd_col (mpq_QSprob p,
							 int cnt,
							 int *cmatind,
							 mpq_t * cmatval,
							 mpq_t obj,
							 mpq_t lower,
							 mpq_t upper,
							 const char *name),
	  mpq_QSnew_row (mpq_QSprob p,
							 mpq_t rhs,
							 int sense,
							 const char *name),
	  mpq_QSadd_rows (mpq_QSprob p,
								int num,
								int *rmatcnt,
								int *rmatbeg,
								int *rmatind,
								mpq_t * rmatval,
								mpq_t * rhs,
								char *sense,
								const char **names),
	  mpq_QSadd_row (mpq_QSprob p,
							 int cnt,
							 int *rmatind,
							 mpq_t * rmatval,
							 mpq_t * rhs,
							 int sense,
							 const char *name),
	  mpq_QSdelete_rows (mpq_QSprob p,
									 int num,
									 int *dellist),
	  mpq_QSdelete_row (mpq_QSprob p,
									int rowindex),
	  mpq_QSdelete_setrows (mpq_QSprob p,
											int *flags),
	  mpq_QSdelete_named_row (mpq_QSprob p,
												const char *rowname),
	  mpq_QSdelete_named_rows_list (mpq_QSprob p,
															int num,
															const char **rownames),
	  mpq_QSdelete_cols (mpq_QSprob p,
									 int num,
									 int *dellist),
	  mpq_QSdelete_col (mpq_QSprob p,
									int colindex),
	  mpq_QSdelete_setcols (mpq_QSprob p,
											int *flags),
	  mpq_QSdelete_named_column (mpq_QSprob p,
													 const char *colname),
	  mpq_QSdelete_named_columns_list (mpq_QSprob p,
																 int num,
																 const char **colnames),
	  mpq_QSchange_senses (mpq_QSprob p,
										 int num,
										 int *rowlist,
										 char *sense),
	  mpq_QSchange_sense (mpq_QSprob p,
										int rowindex,
										int sense),
	  mpq_QSchange_coef (mpq_QSprob p,
									 int rowindex,
									 int colindex,
									 mpq_t coef),
	  mpq_QSchange_objcoef (mpq_QSprob p,
											int indx,
											mpq_t coef),
	  mpq_QSchange_rhscoef (mpq_QSprob p,
											int indx,
											mpq_t coef),
	  mpq_QSchange_bounds (mpq_QSprob p,
										 int num,
										 int *collist,
										 char *lu,
										 mpq_t * bounds),
	  mpq_QSchange_bound (mpq_QSprob p,
										int indx,
										int lu,
										mpq_t bound),
	  mpq_QSload_basis (mpq_QSprob p,
									mpq_QSbas B),
	  mpq_QSread_and_load_basis (mpq_QSprob p,
													 const char *filename),
	  mpq_QSload_basis_array (mpq_QSprob p,
												char *cstat,
												char *rstat),
	  mpq_QSload_basis_and_row_norms_array (mpq_QSprob p,
																			char *cstat,
																			char *rstat,
																			mpq_t * rownorms),
	  mpq_QSget_basis_array (mpq_QSprob p,
											 char *cstat,
											 char *rstat),
	  mpq_QSget_basis_and_row_norms_array (mpq_QSprob p,
																		 char *cstat,
																		 char *rstat,
																		 mpq_t * rownorms),
	  mpq_QSget_binv_row (mpq_QSprob p,
										int indx,
										mpq_t * binvrow),
	  mpq_QSget_tableau_row (mpq_QSprob p,
											 int indx,
											 mpq_t * tableaurow),
	  mpq_QSget_basis_order (mpq_QSprob p,
											 int *basorder),
	  mpq_QSget_status (mpq_QSprob p,
									int *status),
	  mpq_QSget_solution (mpq_QSprob p,
										mpq_t * value,
										mpq_t * x,
										mpq_t * pi,
										mpq_t * slack,
										mpq_t * rc),
	  mpq_QSget_objval (mpq_QSprob p,
									mpq_t * value),
	  mpq_QSget_pi_array (mpq_QSprob p,
										mpq_t * pi),
	  mpq_QSget_rc_array (mpq_QSprob p,
										mpq_t * rc),
	  mpq_QSget_x_array (mpq_QSprob p,
									 mpq_t * x),
	  mpq_QSget_slack_array (mpq_QSprob p,
											 mpq_t * slack),
	  mpq_QSget_infeas_array (mpq_QSprob p,
												mpq_t * pi),
	  mpq_QSget_colcount (mpq_QSprob p),
	  mpq_QSget_rowcount (mpq_QSprob p),
	  mpq_QSget_nzcount (mpq_QSprob p),
	  mpq_QSget_obj (mpq_QSprob p,
							 mpq_t * obj),
	  mpq_QSget_rhs (mpq_QSprob p,
							 mpq_t * rhs),
	  mpq_QSget_rows_list (mpq_QSprob p,
										 int num,
										 int *rowlist,
										 int **rowcnt,
										 int **rowbeg,
										 int **rowind,
										 mpq_t ** rowval,
										 mpq_t ** rhs,
										 char **sense,
										 char ***names),
	  mpq_QSget_rows (mpq_QSprob p,
								int **rowcnt,
								int **rowbeg,
								int **rowind,
								mpq_t ** rowval,
								mpq_t ** rhs,
								char **sense,
								char ***names),
	  mpq_QSget_columns_list (mpq_QSprob p,
												int num,
												int *collist,
												int **colcnt,
												int **colbeg,
												int **colind,
												mpq_t ** colval,
												mpq_t ** obj,
												mpq_t ** lower,
												mpq_t ** upper,
												char ***names),
	  mpq_QSget_columns (mpq_QSprob p,
									 int **colcnt,
									 int **colbeg,
									 int **colind,
									 mpq_t ** colval,
									 mpq_t ** obj,
									 mpq_t ** lower,
									 mpq_t ** upper,
									 char ***names),
	  mpq_QSget_rownames (mpq_QSprob p,
										char **rownames),
	  mpq_QSget_colnames (mpq_QSprob p,
										char **colnames),
	  mpq_QSget_bound (mpq_QSprob p,
								 int colindex,
								 int lu,
								 mpq_t * bound),
	  mpq_QSget_bounds (mpq_QSprob p,
									mpq_t * lower,
									mpq_t * upper),
	  mpq_QSget_intflags (mpq_QSprob p,
										int *intflags),
	  mpq_QSget_intcount (mpq_QSprob p,
										int *count),
	  mpq_QSget_column_index (mpq_QSprob p,
												const char *name,
												int *colindex),
	  mpq_QSget_row_index (mpq_QSprob p,
										 const char *name,
										 int *rowindex),
	  mpq_QSget_named_x (mpq_QSprob p,
									 const char *colname,
									 mpq_t * val),
	  mpq_QSget_named_rc (mpq_QSprob p,
										const char *colname,
										mpq_t * val),
	  mpq_QSget_named_pi (mpq_QSprob p,
										const char *rowname,
										mpq_t * val),
	  mpq_QSget_named_slack (mpq_QSprob p,
											 const char *rowname,
											 mpq_t * val),
	  mpq_QScompute_row_norms (mpq_QSprob p),
	  mpq_QSwrite_prob (mpq_QSprob p,
									const char *filename,
									const char *filetype),
	  mpq_QSwrite_prob_file (mpq_QSprob p,
											 FILE * file,
											 const char *filetype),
	  mpq_QSwrite_basis (mpq_QSprob p,
									 mpq_QSbas B,
									 const char *filename),
	  mpq_QStest_row_norms (mpq_QSprob p),
	  mpq_QSset_param (mpq_QSprob p,
								 int whichparam,
								 int newvalue),
	  mpq_QSset_param_EGlpNum (mpq_QSprob p,
												 int whichparam,
												 mpq_t newvalue),
	  mpq_QSget_param (mpq_QSprob p,
								 int whichparam,
								 int *value),
	  mpq_QSget_param_EGlpNum (mpq_QSprob p,
												 int whichparam,
												 mpq_t * value);

	mpq_QSLIB_INTERFACE char *mpq_QSget_probname (mpq_QSprob p);
	mpq_QSLIB_INTERFACE char *mpq_QSget_objname (mpq_QSprob p);
	mpq_QSLIB_INTERFACE char *mpq_QSversion (void);

	mpq_QSLIB_INTERFACE mpq_QSprob mpq_QScreate_prob (const char *name,
																				int objsense),
	  mpq_QSread_prob (const char *filename,
								 const char *filetype),
	  mpq_QSload_prob (const char *probname,
								 int ncols,
								 int nrows,
								 int *cmatcnt,
								 int *cmatbeg,
								 int *cmatind,
								 mpq_t * cmatval,
								 int objsense,
								 mpq_t * obj,
								 mpq_t * rhs,
								 char *sense,
								 mpq_t * lower,
								 mpq_t * upper,
								 const char **colnames,
								 const char **rownames),
	  mpq_QScopy_prob (mpq_QSprob p,
								 const char *newname);

	mpq_QSLIB_INTERFACE mpq_QSbas mpq_QSget_basis (mpq_QSprob p),
	  mpq_QSread_basis (mpq_QSprob p,
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
 * mpq_feedback/progress information
 */

#ifdef WIN32
typedef struct mpq_QSLIB_INTERFACE mpq_qsline_reader *mpq_QSline_reader;
typedef struct mpq_QSLIB_INTERFACE mpq_qsformat_error *mpq_QSformat_error;
typedef struct mpq_QSLIB_INTERFACE mpq_qserror_collector *mpq_QSerror_collector;
typedef struct mpq_QSLIB_INTERFACE mpq_qserror_memory *mpq_QSerror_memory;
#else
typedef struct mpq_qsline_reader *mpq_QSline_reader;
typedef struct mpq_qsformat_error *mpq_QSformat_error;
typedef struct mpq_qserror_collector *mpq_QSerror_collector;
typedef struct mpq_qserror_memory *mpq_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
	mpq_QSLIB_INTERFACE const char *mpq_QSformat_error_type_string (int tp);

	mpq_QSLIB_INTERFACE int mpq_QSerror_get_type (mpq_QSformat_error error);
	mpq_QSLIB_INTERFACE const char *mpq_QSerror_get_desc (mpq_QSformat_error error);
	mpq_QSLIB_INTERFACE int mpq_QSerror_get_line_number (mpq_QSformat_error error);
	mpq_QSLIB_INTERFACE int mpq_QSerror_get_pos (mpq_QSformat_error error);
	mpq_QSLIB_INTERFACE const char *mpq_QSerror_get_line (mpq_QSformat_error error);
	mpq_QSLIB_INTERFACE void mpq_QSerror_print (FILE * f,
																			mpq_QSformat_error error);

	mpq_QSLIB_INTERFACE mpq_QSerror_collector mpq_QSerror_collector_new (void *fct,
																													 void *dest);
	  mpq_QSLIB_INTERFACE
		mpq_QSerror_collector mpq_QSerror_memory_collector_new (mpq_QSerror_memory mem);
	mpq_QSLIB_INTERFACE void mpq_QSerror_collector_free (mpq_QSerror_collector c);

/****************************************************************************
 * line reader 
 */
	mpq_QSLIB_INTERFACE mpq_QSline_reader mpq_QSline_reader_new (void *fct,
																									 void *data_src);
	/* reader->read_line_fct defaults to fgets */

	mpq_QSLIB_INTERFACE void mpq_QSline_reader_free (mpq_QSline_reader reader);

	mpq_QSLIB_INTERFACE void mpq_QSline_reader_set_error_collector (mpq_QSline_reader reader,
																													mpq_QSerror_collector
																													collector);

	mpq_QSLIB_INTERFACE char *mpq_QSline_reader_get (mpq_QSline_reader reader,
																					 char *s,
																					 int size);

	mpq_QSLIB_INTERFACE mpq_QSprob mpq_QSget_prob (mpq_QSline_reader reader,
																		 const char *probname,
																		 const char *filetype);
	/* the MPS and LP parsers uses the fct from reader 
	 * to get to next input line */


/****************************************************************************
 * error memory 
 */
	mpq_QSLIB_INTERFACE mpq_QSerror_memory mpq_QSerror_memory_create (int takeErrorLines);
	mpq_QSLIB_INTERFACE void mpq_QSerror_memory_free (mpq_QSerror_memory mem);

	mpq_QSLIB_INTERFACE int mpq_QSerror_memory_get_nof (mpq_QSerror_memory mem,
																							int error_type);
	mpq_QSLIB_INTERFACE int mpq_QSerror_memory_get_nerrors (mpq_QSerror_memory mem);

	mpq_QSLIB_INTERFACE mpq_QSformat_error
		mpq_QSerror_memory_get_last_error (mpq_QSerror_memory mem);
	mpq_QSLIB_INTERFACE mpq_QSformat_error
		mpq_QSerror_memory_get_prev_error (mpq_QSformat_error e);

/**************************************************************************** 
 * reporter for solver mpq_feedback 
 */
	mpq_QSLIB_INTERFACE void mpq_QSset_reporter (mpq_QSprob prob,
																			 int iterskip,
																			 void *fct,
																			 void *dest);

	mpq_QSLIB_INTERFACE int mpq_QSreport_prob (mpq_QSprob p,
																		 const char *filetype,
																		 mpq_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif													/* mpq___QS_QSOPT_H */
