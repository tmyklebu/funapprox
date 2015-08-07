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

/* RCSINFO $Id: dbl_lpdefs.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef dbl___QS_LPDEFS_H
#define dbl___QS_LPDEFS_H

#include "dbl_qsopt.h"
#include "dbl_lpdata.h"
#include "dbl_factor.h"

/* infinity and negative infinity */
#define dbl_INFTY  dbl_ILL_MAXDOUBLE
#define dbl_NINFTY dbl_ILL_MINDOUBLE

#include "basicdefs.h"
/* tolerances, these are initialized in dbl_ILLstart, file dbl_lpdata.c */
/* #if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE */
/* these three constants are defined in dbl_lpdata.c */
extern double dbl_PARAM_IBASIS_RPIVOT;	/*      0.98 */
extern double dbl_PARAM_IBASIS_RTRIANG;	/*     0.01 */
extern double dbl_PARAM_MIN_DNORM;	/*          1e-24 */
extern double dbl_PFEAS_TOLER;		/*            1e-6 */
extern double dbl_BD_TOLER;			/*            1e-7 */
extern double dbl_DFEAS_TOLER;		/*            1e-6 */
extern double dbl_PIVOT_TOLER;		/*            1e-10 */
extern double dbl_SZERO_TOLER;		/*            1e-15 */
extern double dbl_PIVZ_TOLER;		/*            1e-12 */
extern double dbl_OBJBND_TOLER;	/*            1e-2 */
extern double dbl_DBNDPIV_TOLER;	/*            1e-3 */
extern double dbl_DBNDPIV_RATIO;	/*            1e-2 */
extern double dbl_PNSTEP_TOLER;	/*            1e-9 */
extern double dbl_DNSTEP_TOLER;	/*            1e-9 */
extern double dbl_ALTPIV_TOLER;	/*            1e-8 */
extern double dbl_DJZERO_TOLER;	/*            1e-8 */
extern double dbl_PROGRESS_ZERO;	/*            1e-7 */
extern double dbl_PROGRESS_THRESH;	/*           1e-5 */
extern double dbl_CB_EPS;				/*           0.001 */
extern double dbl_CB_INF_RATIO;	/*           10.0 */
extern double dbl_CB_PRI_RLIMIT;	/*           0.25 */
/* structure for statistics */
typedef struct
{
	int ynz_cnt;									/* nz in entering columns */
	int num_y;
	double y_ravg;							/* weighted avg. of current & prior y */
	int znz_cnt;									/* nz in ith row of B^{-1}, ie z_i */
	int num_z;
	double z_ravg;							/* weighted avg. of current & prior z */
	int zanz_cnt;									/* nz in z^TA */
	int num_za;
	double za_ravg;						/* weighted avg. of current & prior za */
	int pnorm_cnt;								/* nz in columns for primal norms */
	int dnorm_cnt;								/* nz in rows for dual norms */
	int pinz_cnt;									/* nz in phase II pi (solve) */
	int num_pi;										/* # of pi solves */
	int pi1nz_cnt;								/* nz in phase I pi (solve) */
	int num_pi1;									/* # of phase I pi solves */
	int upnz_cnt;									/* nz in ftran update vector */
	int num_up;										/* # of ftran_updates */
	int pupv_cnt;									/* nz in primal steep updates */
	int dupv_cnt;									/* nz in dual steep updates */

	int start_slacks;							/* # slacks in beginning */
	int final_slacks;							/* # slacks in the end */
	int start_art;								/* # arts in beginning */
	int final_art;								/* # arts in the end */

	int pI_iter;									/* primal phase I iterations */
	int pII_iter;
	int dI_iter;									/* dual phase I iterations */
	int dII_iter;
	int tot_iter;

	int pivpI[10];								/* sizes of pivots */
	int pivpII[10];
	int pivdI[10];
	int pivdII[10];
}
dbl_count_struct;

/* structure for tolerances */
typedef struct
{
	double pfeas_tol;
	double dfeas_tol;
	double pivot_tol;
	double szero_tol;
	double ip_tol;							/* inner primal & dual feas toler */
	double id_tol;
}
dbl_tol_struct;

/* bound information */
typedef struct dbl_bndinfo
{
	double pbound;
	double cbound;
	int btype;
	int varnum;
	struct dbl_bndinfo *next;
}
dbl_bndinfo;

/* bound information */
typedef struct dbl_coefinfo
{
	double pcoef;
	double ccoef;
	int varnum;
	struct dbl_coefinfo *next;
}
dbl_coefinfo;

/* feasibility info */
typedef struct dbl_feas_info
{
	int pstatus;
	int dstatus;
	double totinfeas;
}
dbl_feas_info;

typedef struct dbl_lp_status_info
{
	char optimal;
	char primal_feasible;
	char primal_infeasible;
	char primal_unbounded;
	char dual_feasible;
	char dual_infeasible;
	char dual_unbounded;
	char padd;
}
dbl_lp_status_info;

typedef struct dbl_pI_uinfo
{
	int tctr;
	int i;
	int *perm;
	int *ix;
	int fs;
	double piv;
	double *t;
	double dty;
	double c_obj;
	double tz;
}
dbl_pI_uinfo;

extern void dbl_ILLlp_status_info_init (dbl_lp_status_info * ls);

/* structure for local lp information
 * contains lp obj values - status - dimensions - input data -
 * solution vecs - basis info - update vecs - work vecs - bound changes -
 * tolerances - time info - statistics 
 */
typedef struct dbl_lpinfo
{

	double objval;							/* obj info */
	double pobjval;						/* intermediate status info */
	double dobjval;
	double pinfeas;
	double dinfeas;
	double objbound;
	dbl_lp_status_info probstat;			/* final status */
	dbl_lp_status_info basisstat;			/* final status */
	int nrows;										/* input info follows; given in col format */
	int ncols;
	int *matcnt;
	int *matbeg;
	int *matind;
	double *matval;
	int matfree;
	int matsize;
	double *bz;
	double *lz;
	double *uz;
	double *cz;
	int localrows;								/* set to 1 if these are created locally */
	int *rowcnt;									/* row info follows, copy of col info */
	int *rowbeg;
	int *rowind;
	double *rowval;

	double *xbz;								/* output info x, pi, reduced cost */
	double *piz;
	double *dz;
	double *pIxbz;							/* output info (phase I) x, pi, reduced cost */
	double *pIpiz;
	double *pIdz;

	int final_phase;							/* final phase, inf & unboundedness info */
	int infub_ix;

	int basisid;									/* basis and variable info follows */
	int nnbasic;
	int *baz;
	int *nbaz;
	int *vstat;
	int *vindex;
	int fbasisid;
	dbl_factor_work *f;
	int *vtype;										/* internal var info */
	char *vclass;									/* structural or logical */

	dbl_svector zz;										/* local dbl_ILLfactor_update vectors z, yj, za */
	dbl_svector yjz;
	dbl_svector zA;
	dbl_svector work;									/* local work vector */
	dbl_svector srhs;									/* local vectors for lin. eq. solves */
	dbl_svector ssoln;
	int *iwork;										/* local work vector */
	dbl_pI_uinfo upd;									/* phase I update info */
	int *bfeas;										/* primal and dual infeasibility info */
	int *dfeas;

	dbl_tol_struct *tol;							/* tolerances */
	dbl_count_struct *cnts;						/* counts */
	int nbchange;									/* # bound shifts */
	int ncchange;									/* # obj coef shifts */
	dbl_bndinfo *bchanges;						/* list of bound shifts */
	dbl_coefinfo *cchanges;						/* list of coef shifts */
	int pIratio;									/* ratio tests */
	int pIIratio;
	int dIratio;
	int dIIratio;

	int maxiter;
	int iterskip;
	double maxtime;
	double starttime;
	struct dbl_ILLlpdata *O;
	ILLrandstate rstate;

}
dbl_lpinfo;

/* pricing structures */
typedef struct
{
	int ninit;
	double *norms;
	int *refframe;
}
dbl_p_devex_info;

typedef struct
{
	double *norms;
}
dbl_p_steep_info;

typedef struct
{
	int k;
	int cgroup;
	int ngroups;
	int *gstart;
	int *gshift;
	int *gsize;
	int bsize;
	int *bucket;
	int *perm;
	double *infeas;
}
dbl_mpart_info;

typedef struct
{
	int ninit;
	double *norms;
	int *refframe;
}
dbl_d_devex_info;

typedef struct
{
	double *norms;
}
dbl_d_steep_info;

/* pricing information */
typedef struct dbl_price_info
{
	int p_strategy;
	int d_strategy;
	int pI_price;
	int pII_price;
	int dI_price;
	int dII_price;
	int cur_price;
	double *p_scaleinf;
	double *d_scaleinf;
	dbl_p_devex_info pdinfo;
	dbl_p_steep_info psinfo;
	dbl_mpart_info pmpinfo;
	dbl_d_devex_info ddinfo;
	dbl_d_steep_info dsinfo;
	dbl_mpart_info dmpinfo;
	dbl_heap h;
	double htrigger;
	int hineff;
	char init;
}
dbl_price_info;

#endif /* dbl___QS_LPDEFS_H */
