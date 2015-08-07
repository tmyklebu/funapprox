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

/* RCSINFO $Id: mpf_lpdefs.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpf___QS_LPDEFS_H
#define mpf___QS_LPDEFS_H

#include "mpf_qsopt.h"
#include "mpf_lpdata.h"
#include "mpf_factor.h"

/* infinity and negative infinity */
#define mpf_INFTY  mpf_ILL_MAXDOUBLE
#define mpf_NINFTY mpf_ILL_MINDOUBLE

#include "basicdefs.h"
/* tolerances, these are initialized in mpf_ILLstart, file mpf_lpdata.c */
/* #if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE */
/* these three constants are defined in mpf_lpdata.c */
extern mpf_t mpf_PARAM_IBASIS_RPIVOT;	/*      0.98 */
extern mpf_t mpf_PARAM_IBASIS_RTRIANG;	/*     0.01 */
extern mpf_t mpf_PARAM_MIN_DNORM;	/*          1e-24 */
extern mpf_t mpf_PFEAS_TOLER;		/*            1e-6 */
extern mpf_t mpf_BD_TOLER;			/*            1e-7 */
extern mpf_t mpf_DFEAS_TOLER;		/*            1e-6 */
extern mpf_t mpf_PIVOT_TOLER;		/*            1e-10 */
extern mpf_t mpf_SZERO_TOLER;		/*            1e-15 */
extern mpf_t mpf_PIVZ_TOLER;		/*            1e-12 */
extern mpf_t mpf_OBJBND_TOLER;	/*            1e-2 */
extern mpf_t mpf_DBNDPIV_TOLER;	/*            1e-3 */
extern mpf_t mpf_DBNDPIV_RATIO;	/*            1e-2 */
extern mpf_t mpf_PNSTEP_TOLER;	/*            1e-9 */
extern mpf_t mpf_DNSTEP_TOLER;	/*            1e-9 */
extern mpf_t mpf_ALTPIV_TOLER;	/*            1e-8 */
extern mpf_t mpf_DJZERO_TOLER;	/*            1e-8 */
extern mpf_t mpf_PROGRESS_ZERO;	/*            1e-7 */
extern mpf_t mpf_PROGRESS_THRESH;	/*           1e-5 */
extern mpf_t mpf_CB_EPS;				/*           0.001 */
extern mpf_t mpf_CB_INF_RATIO;	/*           10.0 */
extern mpf_t mpf_CB_PRI_RLIMIT;	/*           0.25 */
/* structure for statistics */
typedef struct
{
	int ynz_cnt;									/* nz in entering columns */
	int num_y;
	mpf_t y_ravg;							/* weighted avg. of current & prior y */
	int znz_cnt;									/* nz in ith row of B^{-1}, ie z_i */
	int num_z;
	mpf_t z_ravg;							/* weighted avg. of current & prior z */
	int zanz_cnt;									/* nz in z^TA */
	int num_za;
	mpf_t za_ravg;						/* weighted avg. of current & prior za */
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
mpf_count_struct;

/* structure for tolerances */
typedef struct
{
	mpf_t pfeas_tol;
	mpf_t dfeas_tol;
	mpf_t pivot_tol;
	mpf_t szero_tol;
	mpf_t ip_tol;							/* inner primal & dual feas toler */
	mpf_t id_tol;
}
mpf_tol_struct;

/* bound information */
typedef struct mpf_bndinfo
{
	mpf_t pbound;
	mpf_t cbound;
	int btype;
	int varnum;
	struct mpf_bndinfo *next;
}
mpf_bndinfo;

/* bound information */
typedef struct mpf_coefinfo
{
	mpf_t pcoef;
	mpf_t ccoef;
	int varnum;
	struct mpf_coefinfo *next;
}
mpf_coefinfo;

/* feasibility info */
typedef struct mpf_feas_info
{
	int pstatus;
	int dstatus;
	mpf_t totinfeas;
}
mpf_feas_info;

typedef struct mpf_lp_status_info
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
mpf_lp_status_info;

typedef struct mpf_pI_uinfo
{
	int tctr;
	int i;
	int *perm;
	int *ix;
	int fs;
	mpf_t piv;
	mpf_t *t;
	mpf_t dty;
	mpf_t c_obj;
	mpf_t tz;
}
mpf_pI_uinfo;

extern void mpf_ILLlp_status_info_init (mpf_lp_status_info * ls);

/* structure for local lp information
 * contains lp obj values - status - dimensions - input data -
 * solution vecs - basis info - update vecs - work vecs - bound changes -
 * tolerances - time info - statistics 
 */
typedef struct mpf_lpinfo
{

	mpf_t objval;							/* obj info */
	mpf_t pobjval;						/* intermediate status info */
	mpf_t dobjval;
	mpf_t pinfeas;
	mpf_t dinfeas;
	mpf_t objbound;
	mpf_lp_status_info probstat;			/* final status */
	mpf_lp_status_info basisstat;			/* final status */
	int nrows;										/* input info follows; given in col format */
	int ncols;
	int *matcnt;
	int *matbeg;
	int *matind;
	mpf_t *matval;
	int matfree;
	int matsize;
	mpf_t *bz;
	mpf_t *lz;
	mpf_t *uz;
	mpf_t *cz;
	int localrows;								/* set to 1 if these are created locally */
	int *rowcnt;									/* row info follows, copy of col info */
	int *rowbeg;
	int *rowind;
	mpf_t *rowval;

	mpf_t *xbz;								/* output info x, pi, reduced cost */
	mpf_t *piz;
	mpf_t *dz;
	mpf_t *pIxbz;							/* output info (phase I) x, pi, reduced cost */
	mpf_t *pIpiz;
	mpf_t *pIdz;

	int final_phase;							/* final phase, inf & unboundedness info */
	int infub_ix;

	int basisid;									/* basis and variable info follows */
	int nnbasic;
	int *baz;
	int *nbaz;
	int *vstat;
	int *vindex;
	int fbasisid;
	mpf_factor_work *f;
	int *vtype;										/* internal var info */
	char *vclass;									/* structural or logical */

	mpf_svector zz;										/* local mpf_ILLfactor_update vectors z, yj, za */
	mpf_svector yjz;
	mpf_svector zA;
	mpf_svector work;									/* local work vector */
	mpf_svector srhs;									/* local vectors for lin. eq. solves */
	mpf_svector ssoln;
	int *iwork;										/* local work vector */
	mpf_pI_uinfo upd;									/* phase I update info */
	int *bfeas;										/* primal and dual infeasibility info */
	int *dfeas;

	mpf_tol_struct *tol;							/* tolerances */
	mpf_count_struct *cnts;						/* counts */
	int nbchange;									/* # bound shifts */
	int ncchange;									/* # obj coef shifts */
	mpf_bndinfo *bchanges;						/* list of bound shifts */
	mpf_coefinfo *cchanges;						/* list of coef shifts */
	int pIratio;									/* ratio tests */
	int pIIratio;
	int dIratio;
	int dIIratio;

	int maxiter;
	int iterskip;
	double maxtime;
	double starttime;
	struct mpf_ILLlpdata *O;
	ILLrandstate rstate;

}
mpf_lpinfo;

/* pricing structures */
typedef struct
{
	int ninit;
	mpf_t *norms;
	int *refframe;
}
mpf_p_devex_info;

typedef struct
{
	mpf_t *norms;
}
mpf_p_steep_info;

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
	mpf_t *infeas;
}
mpf_mpart_info;

typedef struct
{
	int ninit;
	mpf_t *norms;
	int *refframe;
}
mpf_d_devex_info;

typedef struct
{
	mpf_t *norms;
}
mpf_d_steep_info;

/* pricing information */
typedef struct mpf_price_info
{
	int p_strategy;
	int d_strategy;
	int pI_price;
	int pII_price;
	int dI_price;
	int dII_price;
	int cur_price;
	mpf_t *p_scaleinf;
	mpf_t *d_scaleinf;
	mpf_p_devex_info pdinfo;
	mpf_p_steep_info psinfo;
	mpf_mpart_info pmpinfo;
	mpf_d_devex_info ddinfo;
	mpf_d_steep_info dsinfo;
	mpf_mpart_info dmpinfo;
	mpf_heap h;
	mpf_t htrigger;
	int hineff;
	char init;
}
mpf_price_info;

#endif /* mpf___QS_LPDEFS_H */
