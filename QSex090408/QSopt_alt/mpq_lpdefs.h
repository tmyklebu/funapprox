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

/* RCSINFO $Id: mpq_lpdefs.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpq___QS_LPDEFS_H
#define mpq___QS_LPDEFS_H

#include "mpq_qsopt.h"
#include "mpq_lpdata.h"
#include "mpq_factor.h"

/* infinity and negative infinity */
#define mpq_INFTY  mpq_ILL_MAXDOUBLE
#define mpq_NINFTY mpq_ILL_MINDOUBLE

#include "basicdefs.h"
/* tolerances, these are initialized in mpq_ILLstart, file mpq_lpdata.c */
/* #if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE */
/* these three constants are defined in mpq_lpdata.c */
extern mpq_t mpq_PARAM_IBASIS_RPIVOT;	/*      0.98 */
extern mpq_t mpq_PARAM_IBASIS_RTRIANG;	/*     0.01 */
extern mpq_t mpq_PARAM_MIN_DNORM;	/*          1e-24 */
extern mpq_t mpq_PFEAS_TOLER;		/*            1e-6 */
extern mpq_t mpq_BD_TOLER;			/*            1e-7 */
extern mpq_t mpq_DFEAS_TOLER;		/*            1e-6 */
extern mpq_t mpq_PIVOT_TOLER;		/*            1e-10 */
extern mpq_t mpq_SZERO_TOLER;		/*            1e-15 */
extern mpq_t mpq_PIVZ_TOLER;		/*            1e-12 */
extern mpq_t mpq_OBJBND_TOLER;	/*            1e-2 */
extern mpq_t mpq_DBNDPIV_TOLER;	/*            1e-3 */
extern mpq_t mpq_DBNDPIV_RATIO;	/*            1e-2 */
extern mpq_t mpq_PNSTEP_TOLER;	/*            1e-9 */
extern mpq_t mpq_DNSTEP_TOLER;	/*            1e-9 */
extern mpq_t mpq_ALTPIV_TOLER;	/*            1e-8 */
extern mpq_t mpq_DJZERO_TOLER;	/*            1e-8 */
extern mpq_t mpq_PROGRESS_ZERO;	/*            1e-7 */
extern mpq_t mpq_PROGRESS_THRESH;	/*           1e-5 */
extern mpq_t mpq_CB_EPS;				/*           0.001 */
extern mpq_t mpq_CB_INF_RATIO;	/*           10.0 */
extern mpq_t mpq_CB_PRI_RLIMIT;	/*           0.25 */
/* structure for statistics */
typedef struct
{
	int ynz_cnt;									/* nz in entering columns */
	int num_y;
	mpq_t y_ravg;							/* weighted avg. of current & prior y */
	int znz_cnt;									/* nz in ith row of B^{-1}, ie z_i */
	int num_z;
	mpq_t z_ravg;							/* weighted avg. of current & prior z */
	int zanz_cnt;									/* nz in z^TA */
	int num_za;
	mpq_t za_ravg;						/* weighted avg. of current & prior za */
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
mpq_count_struct;

/* structure for tolerances */
typedef struct
{
	mpq_t pfeas_tol;
	mpq_t dfeas_tol;
	mpq_t pivot_tol;
	mpq_t szero_tol;
	mpq_t ip_tol;							/* inner primal & dual feas toler */
	mpq_t id_tol;
}
mpq_tol_struct;

/* bound information */
typedef struct mpq_bndinfo
{
	mpq_t pbound;
	mpq_t cbound;
	int btype;
	int varnum;
	struct mpq_bndinfo *next;
}
mpq_bndinfo;

/* bound information */
typedef struct mpq_coefinfo
{
	mpq_t pcoef;
	mpq_t ccoef;
	int varnum;
	struct mpq_coefinfo *next;
}
mpq_coefinfo;

/* feasibility info */
typedef struct mpq_feas_info
{
	int pstatus;
	int dstatus;
	mpq_t totinfeas;
}
mpq_feas_info;

typedef struct mpq_lp_status_info
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
mpq_lp_status_info;

typedef struct mpq_pI_uinfo
{
	int tctr;
	int i;
	int *perm;
	int *ix;
	int fs;
	mpq_t piv;
	mpq_t *t;
	mpq_t dty;
	mpq_t c_obj;
	mpq_t tz;
}
mpq_pI_uinfo;

extern void mpq_ILLlp_status_info_init (mpq_lp_status_info * ls);

/* structure for local lp information
 * contains lp obj values - status - dimensions - input data -
 * solution vecs - basis info - update vecs - work vecs - bound changes -
 * tolerances - time info - statistics 
 */
typedef struct mpq_lpinfo
{

	mpq_t objval;							/* obj info */
	mpq_t pobjval;						/* intermediate status info */
	mpq_t dobjval;
	mpq_t pinfeas;
	mpq_t dinfeas;
	mpq_t objbound;
	mpq_lp_status_info probstat;			/* final status */
	mpq_lp_status_info basisstat;			/* final status */
	int nrows;										/* input info follows; given in col format */
	int ncols;
	int *matcnt;
	int *matbeg;
	int *matind;
	mpq_t *matval;
	int matfree;
	int matsize;
	mpq_t *bz;
	mpq_t *lz;
	mpq_t *uz;
	mpq_t *cz;
	int localrows;								/* set to 1 if these are created locally */
	int *rowcnt;									/* row info follows, copy of col info */
	int *rowbeg;
	int *rowind;
	mpq_t *rowval;

	mpq_t *xbz;								/* output info x, pi, reduced cost */
	mpq_t *piz;
	mpq_t *dz;
	mpq_t *pIxbz;							/* output info (phase I) x, pi, reduced cost */
	mpq_t *pIpiz;
	mpq_t *pIdz;

	int final_phase;							/* final phase, inf & unboundedness info */
	int infub_ix;

	int basisid;									/* basis and variable info follows */
	int nnbasic;
	int *baz;
	int *nbaz;
	int *vstat;
	int *vindex;
	int fbasisid;
	mpq_factor_work *f;
	int *vtype;										/* internal var info */
	char *vclass;									/* structural or logical */

	mpq_svector zz;										/* local mpq_ILLfactor_update vectors z, yj, za */
	mpq_svector yjz;
	mpq_svector zA;
	mpq_svector work;									/* local work vector */
	mpq_svector srhs;									/* local vectors for lin. eq. solves */
	mpq_svector ssoln;
	int *iwork;										/* local work vector */
	mpq_pI_uinfo upd;									/* phase I update info */
	int *bfeas;										/* primal and dual infeasibility info */
	int *dfeas;

	mpq_tol_struct *tol;							/* tolerances */
	mpq_count_struct *cnts;						/* counts */
	int nbchange;									/* # bound shifts */
	int ncchange;									/* # obj coef shifts */
	mpq_bndinfo *bchanges;						/* list of bound shifts */
	mpq_coefinfo *cchanges;						/* list of coef shifts */
	int pIratio;									/* ratio tests */
	int pIIratio;
	int dIratio;
	int dIIratio;

	int maxiter;
	int iterskip;
	double maxtime;
	double starttime;
	struct mpq_ILLlpdata *O;
	ILLrandstate rstate;

}
mpq_lpinfo;

/* pricing structures */
typedef struct
{
	int ninit;
	mpq_t *norms;
	int *refframe;
}
mpq_p_devex_info;

typedef struct
{
	mpq_t *norms;
}
mpq_p_steep_info;

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
	mpq_t *infeas;
}
mpq_mpart_info;

typedef struct
{
	int ninit;
	mpq_t *norms;
	int *refframe;
}
mpq_d_devex_info;

typedef struct
{
	mpq_t *norms;
}
mpq_d_steep_info;

/* pricing information */
typedef struct mpq_price_info
{
	int p_strategy;
	int d_strategy;
	int pI_price;
	int pII_price;
	int dI_price;
	int dII_price;
	int cur_price;
	mpq_t *p_scaleinf;
	mpq_t *d_scaleinf;
	mpq_p_devex_info pdinfo;
	mpq_p_steep_info psinfo;
	mpq_mpart_info pmpinfo;
	mpq_d_devex_info ddinfo;
	mpq_d_steep_info dsinfo;
	mpq_mpart_info dmpinfo;
	mpq_heap h;
	mpq_t htrigger;
	int hineff;
	char init;
}
mpq_price_info;

#endif /* mpq___QS_LPDEFS_H */
