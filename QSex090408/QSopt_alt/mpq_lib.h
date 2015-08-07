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

/* RCSINFO $Id: mpq_lib.h,v 1.4 2003/11/05 17:00:26 meven Exp $ */
#ifndef mpq_ILL_LIB_H
#define mpq_ILL_LIB_H

#include "mpq_lpdefs.h"
#include "mpq_lpdata.h"
#include "mpq_price.h"
#include "basicdefs.h"

/****************************************************************************/
/*                                                                          */
/*                   Return Status for mpq_ILLlib_optimize                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                               mpq_lib.c                                      */
/*                                                                          */
/****************************************************************************/


int mpq_ILLlib_optimize (mpq_lpinfo * lp,
										 mpq_ILLlp_basis * B,
										 mpq_price_info * pinf,
										 int algo,
										 int *status,
										 int simplex_display),
  mpq_ILLlib_cache_solution (mpq_lpinfo * lp,
												 mpq_ILLlp_cache * C),
  mpq_ILLlib_solution (mpq_lpinfo * lp,
									 mpq_ILLlp_cache * C,
									 mpq_t * val,
									 mpq_t * x,
									 mpq_t * pi,
									 mpq_t * slack,
									 mpq_t * rc),
  mpq_ILLlib_get_x (mpq_lpinfo * lp,
								mpq_ILLlp_cache * C,
								mpq_t * x),
  mpq_ILLlib_get_slack (mpq_lpinfo * lp,
										mpq_ILLlp_cache * C,
										mpq_t * slack),
  mpq_ILLlib_objval (mpq_lpinfo * lp,
								 mpq_ILLlp_cache * C,
								 mpq_t * val),
  mpq_ILLlib_tableau (mpq_lpinfo * lp,
									int row,
									mpq_t * binv,
									mpq_t * tabrow),
  mpq_ILLlib_basis_order (mpq_lpinfo * lp,
											int *header),
  mpq_ILLlib_newrow (mpq_lpinfo * lp,
								 mpq_ILLlp_basis * B,
								 mpq_t rhs,
								 int sense,
								 mpq_t range,
								 const char *name),
  mpq_ILLlib_newrows (mpq_lpinfo * lp,
									mpq_ILLlp_basis * B,
									int num,
									mpq_t * rhs,
									char *sense,
									mpq_t * range,
									const char **names),
  mpq_ILLlib_addrow (mpq_lpinfo * lp,
								 mpq_ILLlp_basis * B,
								 int cnt,
								 int *ind,
								 mpq_t * val,
								 mpq_t rhs,
								 int sense,
								 mpq_t range,
								 const char *rowname),
  mpq_ILLlib_addrows (mpq_lpinfo * lp,
									mpq_ILLlp_basis * B,
									int num,
									int *rmatcnt,
									int *rmatbeg,
									int *rmatind,
									mpq_t * rmatval,
									mpq_t * rhs,
									char *sense,
									mpq_t * range,
									const char **names,
									int *nofactor),
  mpq_ILLlib_delrows (mpq_lpinfo * lp,
									mpq_ILLlp_basis * B,
									mpq_ILLlp_cache * C,
									int num,
									int *dellist,
									int *basis_ok,
									int *cache_ok),
  mpq_ILLlib_newcol (mpq_lpinfo * lp,
								 mpq_ILLlp_basis * B,
								 mpq_t obj,
								 mpq_t lower,
								 mpq_t upper,
								 const char *name,
								 int factorok),
  mpq_ILLlib_newcols (mpq_lpinfo * lp,
									mpq_ILLlp_basis * B,
									int num,
									mpq_t * obj,
									mpq_t * lower,
									mpq_t * upper,
									const char **names,
									int factorok),
  mpq_ILLlib_addcol (mpq_lpinfo * lp,
								 mpq_ILLlp_basis * B,
								 int cnt,
								 int *ind,
								 mpq_t * val,
								 mpq_t obj,
								 mpq_t lower,
								 mpq_t upper,
								 const char *name,
								 int factorok),
  mpq_ILLlib_addcols (mpq_lpinfo * lp,
									mpq_ILLlp_basis * B,
									int num,
									int *cmatcnt,
									int *cmatbeg,
									int *cmatind,
									mpq_t * cmatval,
									mpq_t * obj,
									mpq_t * lower,
									mpq_t * upper,
									const char **names,
									int factorok),
  mpq_ILLlib_delcols (mpq_lpinfo * lp,
									mpq_ILLlp_basis * B,
									int num,
									int *dellist,
									int *basis_ok),
  mpq_ILLlib_chgcoef (mpq_lpinfo * lp,
									int rowindex,
									int colindex,
									mpq_t coef),
  mpq_ILLlib_chgsense (mpq_lpinfo * lp,
									 int num,
									 int *rowlist,
									 char *sense),
  mpq_ILLlib_getrows (mpq_lpinfo * lp,
									int num,
									int *rowlist,
									int **rowcnt,
									int **rowbeg,
									int **rowind,
									mpq_t ** rowval,
									mpq_t ** rhs,
									char **sense,
									char ***names),
  mpq_ILLlib_getcols (mpq_lpinfo * lp,
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
  mpq_ILLlib_getobj (mpq_lpinfo * lp,
								 mpq_t * obj),
  mpq_ILLlib_chgobj (mpq_lpinfo * lp,
								 int indx,
								 mpq_t coef),
  mpq_ILLlib_getrhs (mpq_lpinfo * lp,
								 mpq_t * rhs),
  mpq_ILLlib_chgrhs (mpq_lpinfo * lp,
								 int indx,
								 mpq_t coef),
  mpq_ILLlib_getintflags (mpq_lpinfo * lp,
											int *intflags),
  mpq_ILLlib_rownames (mpq_lpinfo * lp,
									 char **rownames),
  mpq_ILLlib_colnames (mpq_lpinfo * lp,
									 char **colnames),
  mpq_ILLlib_colindex (mpq_lpinfo * lp,
									 const char *name,
									 int *colindex),
  mpq_ILLlib_rowindex (mpq_lpinfo * lp,
									 const char *name,
									 int *rowindex),
  mpq_ILLlib_chgbnd (mpq_lpinfo * lp,
								 int indx,
								 int lu,
								 mpq_t bnd),
  mpq_ILLlib_chgbnds (mpq_lpinfo * lp,
									int cnt,
									int *indx,
									char *lu,
									mpq_t * bnd),
  mpq_ILLlib_getbnd (mpq_lpinfo * lp,
								 int indx,
								 int lu,
								 mpq_t * bnd),
  mpq_ILLlib_getbnds (mpq_lpinfo * lp,
									mpq_t * lower,
									mpq_t * upper),
  mpq_ILLlib_strongbranch (mpq_lpinfo * lp,
											 mpq_price_info * pinf,
											 int *candidatelist,
											 int ncand,
											 mpq_t * xlist,
											 mpq_t * downpen,
											 mpq_t * uppen,
											 int iterations,
											 mpq_t objbound),
  mpq_ILLlib_getbasis (mpq_lpinfo * lp,
									 char *cstat,
									 char *rstat),
  mpq_ILLlib_loadbasis (mpq_ILLlp_basis * B,
										int nstruct,
										int nrows,
										char *cstat,
										char *rstat),
  mpq_ILLlib_readbasis (mpq_lpinfo * lp,
										mpq_ILLlp_basis * B,
										const char *mpq_fname),
  mpq_ILLlib_writebasis (mpq_lpinfo * lp,
										 mpq_ILLlp_basis * B,
										 const char *mpq_fname),
  mpq_ILLlib_getrownorms (mpq_lpinfo * lp,
											mpq_price_info * pinf,
											mpq_t * rownorms),
  mpq_ILLlib_loadrownorms (mpq_lpinfo * lp,
											 mpq_price_info * pinf,
											 mpq_t * rownorms),
  mpq_ILLlib_recompute_rownorms (mpq_lpinfo * lp,
														 mpq_price_info * pinf),
  mpq_ILLlib_iter (mpq_lpinfo * lp),
  mpq_ILLlib_print_x (FILE * fd,
									mpq_lpinfo * lp,
									mpq_ILLlp_cache * C,
									mpq_t * x,
									int nonZerosOnly),
  mpq_ILLwrite_lp_file (mpq_ILLlpdata * lp,
										FILE * eout,
										mpq_qserror_collector * c);


extern int mpq_ILLlib_findName (mpq_ILLlpdata * qslp,
														int forRow,
														const char *name,
														int id,
														char buf[ILL_namebufsize]);

/****************************************************************************/
/*                                                                          */
/*                           mpq_presolve.c                                     */
/*                                                                          */
/****************************************************************************/

int mpq_ILLpresolve_add_logicals (mpq_ILLlpdata * lp);


/****************************************************************************/
/*                                                                          */
/*                            mpq_binary.c                                      */
/*                                                                          */
/****************************************************************************/

int mpq_ILLmip_binary_dfs (mpq_lpinfo * lp);

#endif /* mpq_ILL_LIB_H */
