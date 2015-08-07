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

/* RCSINFO $Id: mpf_lib.h,v 1.4 2003/11/05 17:00:26 meven Exp $ */
#ifndef mpf_ILL_LIB_H
#define mpf_ILL_LIB_H

#include "mpf_lpdefs.h"
#include "mpf_lpdata.h"
#include "mpf_price.h"
#include "basicdefs.h"

/****************************************************************************/
/*                                                                          */
/*                   Return Status for mpf_ILLlib_optimize                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                               mpf_lib.c                                      */
/*                                                                          */
/****************************************************************************/


int mpf_ILLlib_optimize (mpf_lpinfo * lp,
										 mpf_ILLlp_basis * B,
										 mpf_price_info * pinf,
										 int algo,
										 int *status,
										 int simplex_display),
  mpf_ILLlib_cache_solution (mpf_lpinfo * lp,
												 mpf_ILLlp_cache * C),
  mpf_ILLlib_solution (mpf_lpinfo * lp,
									 mpf_ILLlp_cache * C,
									 mpf_t * val,
									 mpf_t * x,
									 mpf_t * pi,
									 mpf_t * slack,
									 mpf_t * rc),
  mpf_ILLlib_get_x (mpf_lpinfo * lp,
								mpf_ILLlp_cache * C,
								mpf_t * x),
  mpf_ILLlib_get_slack (mpf_lpinfo * lp,
										mpf_ILLlp_cache * C,
										mpf_t * slack),
  mpf_ILLlib_objval (mpf_lpinfo * lp,
								 mpf_ILLlp_cache * C,
								 mpf_t * val),
  mpf_ILLlib_tableau (mpf_lpinfo * lp,
									int row,
									mpf_t * binv,
									mpf_t * tabrow),
  mpf_ILLlib_basis_order (mpf_lpinfo * lp,
											int *header),
  mpf_ILLlib_newrow (mpf_lpinfo * lp,
								 mpf_ILLlp_basis * B,
								 mpf_t rhs,
								 int sense,
								 mpf_t range,
								 const char *name),
  mpf_ILLlib_newrows (mpf_lpinfo * lp,
									mpf_ILLlp_basis * B,
									int num,
									mpf_t * rhs,
									char *sense,
									mpf_t * range,
									const char **names),
  mpf_ILLlib_addrow (mpf_lpinfo * lp,
								 mpf_ILLlp_basis * B,
								 int cnt,
								 int *ind,
								 mpf_t * val,
								 mpf_t rhs,
								 int sense,
								 mpf_t range,
								 const char *rowname),
  mpf_ILLlib_addrows (mpf_lpinfo * lp,
									mpf_ILLlp_basis * B,
									int num,
									int *rmatcnt,
									int *rmatbeg,
									int *rmatind,
									mpf_t * rmatval,
									mpf_t * rhs,
									char *sense,
									mpf_t * range,
									const char **names,
									int *nofactor),
  mpf_ILLlib_delrows (mpf_lpinfo * lp,
									mpf_ILLlp_basis * B,
									mpf_ILLlp_cache * C,
									int num,
									int *dellist,
									int *basis_ok,
									int *cache_ok),
  mpf_ILLlib_newcol (mpf_lpinfo * lp,
								 mpf_ILLlp_basis * B,
								 mpf_t obj,
								 mpf_t lower,
								 mpf_t upper,
								 const char *name,
								 int factorok),
  mpf_ILLlib_newcols (mpf_lpinfo * lp,
									mpf_ILLlp_basis * B,
									int num,
									mpf_t * obj,
									mpf_t * lower,
									mpf_t * upper,
									const char **names,
									int factorok),
  mpf_ILLlib_addcol (mpf_lpinfo * lp,
								 mpf_ILLlp_basis * B,
								 int cnt,
								 int *ind,
								 mpf_t * val,
								 mpf_t obj,
								 mpf_t lower,
								 mpf_t upper,
								 const char *name,
								 int factorok),
  mpf_ILLlib_addcols (mpf_lpinfo * lp,
									mpf_ILLlp_basis * B,
									int num,
									int *cmatcnt,
									int *cmatbeg,
									int *cmatind,
									mpf_t * cmatval,
									mpf_t * obj,
									mpf_t * lower,
									mpf_t * upper,
									const char **names,
									int factorok),
  mpf_ILLlib_delcols (mpf_lpinfo * lp,
									mpf_ILLlp_basis * B,
									int num,
									int *dellist,
									int *basis_ok),
  mpf_ILLlib_chgcoef (mpf_lpinfo * lp,
									int rowindex,
									int colindex,
									mpf_t coef),
  mpf_ILLlib_chgsense (mpf_lpinfo * lp,
									 int num,
									 int *rowlist,
									 char *sense),
  mpf_ILLlib_getrows (mpf_lpinfo * lp,
									int num,
									int *rowlist,
									int **rowcnt,
									int **rowbeg,
									int **rowind,
									mpf_t ** rowval,
									mpf_t ** rhs,
									char **sense,
									char ***names),
  mpf_ILLlib_getcols (mpf_lpinfo * lp,
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
  mpf_ILLlib_getobj (mpf_lpinfo * lp,
								 mpf_t * obj),
  mpf_ILLlib_chgobj (mpf_lpinfo * lp,
								 int indx,
								 mpf_t coef),
  mpf_ILLlib_getrhs (mpf_lpinfo * lp,
								 mpf_t * rhs),
  mpf_ILLlib_chgrhs (mpf_lpinfo * lp,
								 int indx,
								 mpf_t coef),
  mpf_ILLlib_getintflags (mpf_lpinfo * lp,
											int *intflags),
  mpf_ILLlib_rownames (mpf_lpinfo * lp,
									 char **rownames),
  mpf_ILLlib_colnames (mpf_lpinfo * lp,
									 char **colnames),
  mpf_ILLlib_colindex (mpf_lpinfo * lp,
									 const char *name,
									 int *colindex),
  mpf_ILLlib_rowindex (mpf_lpinfo * lp,
									 const char *name,
									 int *rowindex),
  mpf_ILLlib_chgbnd (mpf_lpinfo * lp,
								 int indx,
								 int lu,
								 mpf_t bnd),
  mpf_ILLlib_chgbnds (mpf_lpinfo * lp,
									int cnt,
									int *indx,
									char *lu,
									mpf_t * bnd),
  mpf_ILLlib_getbnd (mpf_lpinfo * lp,
								 int indx,
								 int lu,
								 mpf_t * bnd),
  mpf_ILLlib_getbnds (mpf_lpinfo * lp,
									mpf_t * lower,
									mpf_t * upper),
  mpf_ILLlib_strongbranch (mpf_lpinfo * lp,
											 mpf_price_info * pinf,
											 int *candidatelist,
											 int ncand,
											 mpf_t * xlist,
											 mpf_t * downpen,
											 mpf_t * uppen,
											 int iterations,
											 mpf_t objbound),
  mpf_ILLlib_getbasis (mpf_lpinfo * lp,
									 char *cstat,
									 char *rstat),
  mpf_ILLlib_loadbasis (mpf_ILLlp_basis * B,
										int nstruct,
										int nrows,
										char *cstat,
										char *rstat),
  mpf_ILLlib_readbasis (mpf_lpinfo * lp,
										mpf_ILLlp_basis * B,
										const char *mpf_fname),
  mpf_ILLlib_writebasis (mpf_lpinfo * lp,
										 mpf_ILLlp_basis * B,
										 const char *mpf_fname),
  mpf_ILLlib_getrownorms (mpf_lpinfo * lp,
											mpf_price_info * pinf,
											mpf_t * rownorms),
  mpf_ILLlib_loadrownorms (mpf_lpinfo * lp,
											 mpf_price_info * pinf,
											 mpf_t * rownorms),
  mpf_ILLlib_recompute_rownorms (mpf_lpinfo * lp,
														 mpf_price_info * pinf),
  mpf_ILLlib_iter (mpf_lpinfo * lp),
  mpf_ILLlib_print_x (FILE * fd,
									mpf_lpinfo * lp,
									mpf_ILLlp_cache * C,
									mpf_t * x,
									int nonZerosOnly),
  mpf_ILLwrite_lp_file (mpf_ILLlpdata * lp,
										FILE * eout,
										mpf_qserror_collector * c);


extern int mpf_ILLlib_findName (mpf_ILLlpdata * qslp,
														int forRow,
														const char *name,
														int id,
														char buf[ILL_namebufsize]);

/****************************************************************************/
/*                                                                          */
/*                           mpf_presolve.c                                     */
/*                                                                          */
/****************************************************************************/

int mpf_ILLpresolve_add_logicals (mpf_ILLlpdata * lp);


/****************************************************************************/
/*                                                                          */
/*                            mpf_binary.c                                      */
/*                                                                          */
/****************************************************************************/

int mpf_ILLmip_binary_dfs (mpf_lpinfo * lp);

#endif /* mpf_ILL_LIB_H */
