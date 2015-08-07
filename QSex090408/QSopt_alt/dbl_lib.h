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

/* RCSINFO $Id: dbl_lib.h,v 1.4 2003/11/05 17:00:26 meven Exp $ */
#ifndef dbl_ILL_LIB_H
#define dbl_ILL_LIB_H

#include "dbl_lpdefs.h"
#include "dbl_lpdata.h"
#include "dbl_price.h"
#include "basicdefs.h"

/****************************************************************************/
/*                                                                          */
/*                   Return Status for dbl_ILLlib_optimize                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                               dbl_lib.c                                      */
/*                                                                          */
/****************************************************************************/


int dbl_ILLlib_optimize (dbl_lpinfo * lp,
										 dbl_ILLlp_basis * B,
										 dbl_price_info * pinf,
										 int algo,
										 int *status,
										 int simplex_display),
  dbl_ILLlib_cache_solution (dbl_lpinfo * lp,
												 dbl_ILLlp_cache * C),
  dbl_ILLlib_solution (dbl_lpinfo * lp,
									 dbl_ILLlp_cache * C,
									 double * val,
									 double * x,
									 double * pi,
									 double * slack,
									 double * rc),
  dbl_ILLlib_get_x (dbl_lpinfo * lp,
								dbl_ILLlp_cache * C,
								double * x),
  dbl_ILLlib_get_slack (dbl_lpinfo * lp,
										dbl_ILLlp_cache * C,
										double * slack),
  dbl_ILLlib_objval (dbl_lpinfo * lp,
								 dbl_ILLlp_cache * C,
								 double * val),
  dbl_ILLlib_tableau (dbl_lpinfo * lp,
									int row,
									double * binv,
									double * tabrow),
  dbl_ILLlib_basis_order (dbl_lpinfo * lp,
											int *header),
  dbl_ILLlib_newrow (dbl_lpinfo * lp,
								 dbl_ILLlp_basis * B,
								 double rhs,
								 int sense,
								 double range,
								 const char *name),
  dbl_ILLlib_newrows (dbl_lpinfo * lp,
									dbl_ILLlp_basis * B,
									int num,
									double * rhs,
									char *sense,
									double * range,
									const char **names),
  dbl_ILLlib_addrow (dbl_lpinfo * lp,
								 dbl_ILLlp_basis * B,
								 int cnt,
								 int *ind,
								 double * val,
								 double rhs,
								 int sense,
								 double range,
								 const char *rowname),
  dbl_ILLlib_addrows (dbl_lpinfo * lp,
									dbl_ILLlp_basis * B,
									int num,
									int *rmatcnt,
									int *rmatbeg,
									int *rmatind,
									double * rmatval,
									double * rhs,
									char *sense,
									double * range,
									const char **names,
									int *nofactor),
  dbl_ILLlib_delrows (dbl_lpinfo * lp,
									dbl_ILLlp_basis * B,
									dbl_ILLlp_cache * C,
									int num,
									int *dellist,
									int *basis_ok,
									int *cache_ok),
  dbl_ILLlib_newcol (dbl_lpinfo * lp,
								 dbl_ILLlp_basis * B,
								 double obj,
								 double lower,
								 double upper,
								 const char *name,
								 int factorok),
  dbl_ILLlib_newcols (dbl_lpinfo * lp,
									dbl_ILLlp_basis * B,
									int num,
									double * obj,
									double * lower,
									double * upper,
									const char **names,
									int factorok),
  dbl_ILLlib_addcol (dbl_lpinfo * lp,
								 dbl_ILLlp_basis * B,
								 int cnt,
								 int *ind,
								 double * val,
								 double obj,
								 double lower,
								 double upper,
								 const char *name,
								 int factorok),
  dbl_ILLlib_addcols (dbl_lpinfo * lp,
									dbl_ILLlp_basis * B,
									int num,
									int *cmatcnt,
									int *cmatbeg,
									int *cmatind,
									double * cmatval,
									double * obj,
									double * lower,
									double * upper,
									const char **names,
									int factorok),
  dbl_ILLlib_delcols (dbl_lpinfo * lp,
									dbl_ILLlp_basis * B,
									int num,
									int *dellist,
									int *basis_ok),
  dbl_ILLlib_chgcoef (dbl_lpinfo * lp,
									int rowindex,
									int colindex,
									double coef),
  dbl_ILLlib_chgsense (dbl_lpinfo * lp,
									 int num,
									 int *rowlist,
									 char *sense),
  dbl_ILLlib_getrows (dbl_lpinfo * lp,
									int num,
									int *rowlist,
									int **rowcnt,
									int **rowbeg,
									int **rowind,
									double ** rowval,
									double ** rhs,
									char **sense,
									char ***names),
  dbl_ILLlib_getcols (dbl_lpinfo * lp,
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
  dbl_ILLlib_getobj (dbl_lpinfo * lp,
								 double * obj),
  dbl_ILLlib_chgobj (dbl_lpinfo * lp,
								 int indx,
								 double coef),
  dbl_ILLlib_getrhs (dbl_lpinfo * lp,
								 double * rhs),
  dbl_ILLlib_chgrhs (dbl_lpinfo * lp,
								 int indx,
								 double coef),
  dbl_ILLlib_getintflags (dbl_lpinfo * lp,
											int *intflags),
  dbl_ILLlib_rownames (dbl_lpinfo * lp,
									 char **rownames),
  dbl_ILLlib_colnames (dbl_lpinfo * lp,
									 char **colnames),
  dbl_ILLlib_colindex (dbl_lpinfo * lp,
									 const char *name,
									 int *colindex),
  dbl_ILLlib_rowindex (dbl_lpinfo * lp,
									 const char *name,
									 int *rowindex),
  dbl_ILLlib_chgbnd (dbl_lpinfo * lp,
								 int indx,
								 int lu,
								 double bnd),
  dbl_ILLlib_chgbnds (dbl_lpinfo * lp,
									int cnt,
									int *indx,
									char *lu,
									double * bnd),
  dbl_ILLlib_getbnd (dbl_lpinfo * lp,
								 int indx,
								 int lu,
								 double * bnd),
  dbl_ILLlib_getbnds (dbl_lpinfo * lp,
									double * lower,
									double * upper),
  dbl_ILLlib_strongbranch (dbl_lpinfo * lp,
											 dbl_price_info * pinf,
											 int *candidatelist,
											 int ncand,
											 double * xlist,
											 double * downpen,
											 double * uppen,
											 int iterations,
											 double objbound),
  dbl_ILLlib_getbasis (dbl_lpinfo * lp,
									 char *cstat,
									 char *rstat),
  dbl_ILLlib_loadbasis (dbl_ILLlp_basis * B,
										int nstruct,
										int nrows,
										char *cstat,
										char *rstat),
  dbl_ILLlib_readbasis (dbl_lpinfo * lp,
										dbl_ILLlp_basis * B,
										const char *dbl_fname),
  dbl_ILLlib_writebasis (dbl_lpinfo * lp,
										 dbl_ILLlp_basis * B,
										 const char *dbl_fname),
  dbl_ILLlib_getrownorms (dbl_lpinfo * lp,
											dbl_price_info * pinf,
											double * rownorms),
  dbl_ILLlib_loadrownorms (dbl_lpinfo * lp,
											 dbl_price_info * pinf,
											 double * rownorms),
  dbl_ILLlib_recompute_rownorms (dbl_lpinfo * lp,
														 dbl_price_info * pinf),
  dbl_ILLlib_iter (dbl_lpinfo * lp),
  dbl_ILLlib_print_x (FILE * fd,
									dbl_lpinfo * lp,
									dbl_ILLlp_cache * C,
									double * x,
									int nonZerosOnly),
  dbl_ILLwrite_lp_file (dbl_ILLlpdata * lp,
										FILE * eout,
										dbl_qserror_collector * c);


extern int dbl_ILLlib_findName (dbl_ILLlpdata * qslp,
														int forRow,
														const char *name,
														int id,
														char buf[ILL_namebufsize]);

/****************************************************************************/
/*                                                                          */
/*                           dbl_presolve.c                                     */
/*                                                                          */
/****************************************************************************/

int dbl_ILLpresolve_add_logicals (dbl_ILLlpdata * lp);


/****************************************************************************/
/*                                                                          */
/*                            dbl_binary.c                                      */
/*                                                                          */
/****************************************************************************/

int dbl_ILLmip_binary_dfs (dbl_lpinfo * lp);

#endif /* dbl_ILL_LIB_H */
