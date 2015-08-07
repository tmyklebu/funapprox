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

/* RCSINFO $Id: mpq_fct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef mpq___FUNCTIONS_H
#define mpq___FUNCTIONS_H
#include "basicdefs.h"
int mpq_ILLfct_compute_zA (mpq_lpinfo * lp,
                                             mpq_svector * z,
                                             mpq_svector * zA),
  mpq_ILLfct_compute_wz (mpq_lpinfo * lp,
                                         mpq_t * wz),
  mpq_ILLfct_adjust_viol_bounds (mpq_lpinfo * lp),
  mpq_ILLfct_perturb_bounds (mpq_lpinfo * lp),
  mpq_ILLfct_perturb_phaseI_bounds (mpq_lpinfo * lp),
  mpq_ILLfct_bound_shift (mpq_lpinfo * lp,
                                            int col,
                                            int bndtype,
                                            mpq_t newbnd),
  mpq_ILLfct_adjust_viol_coefs (mpq_lpinfo * lp),
  mpq_ILLfct_perturb_coefs (mpq_lpinfo * lp),
  mpq_ILLfct_coef_shift (mpq_lpinfo * lp,
                                         int col,
                                         mpq_t newcoef),
  mpq_ILLfct_test_pivot (mpq_lpinfo * lp,
                                         int indx,
                                         int indxtype,
                                         mpq_t piv_val);

void mpq_ILLfct_load_workvector (mpq_lpinfo * lp,
                                                         mpq_svector * s),
  mpq_ILLfct_zero_workvector (mpq_lpinfo * lp),
  mpq_ILLfct_set_variable_type (mpq_lpinfo * lp),
  mpq_ILLfct_compute_pobj (mpq_lpinfo * lp),
  mpq_ILLfct_compute_dobj (mpq_lpinfo * lp),
  mpq_ILLfct_compute_xbz (mpq_lpinfo * lp),
  mpq_ILLfct_compute_piz (mpq_lpinfo * lp),
  mpq_ILLfct_compute_phaseI_xbz (mpq_lpinfo * lp),
  mpq_ILLfct_compute_phaseI_piz (mpq_lpinfo * lp),
  mpq_ILLfct_compute_dz (mpq_lpinfo * lp),
  mpq_ILLfct_compute_phaseI_dz (mpq_lpinfo * lp),
  mpq_ILLfct_compute_yz (mpq_lpinfo * lp,
                                         mpq_svector * yz,
                                         mpq_svector * updz,
                                         int ecol),
  mpq_ILLfct_compute_zz (mpq_lpinfo * lp,
                                         mpq_svector * zz,
                                         int lindex),
  mpq_ILLfct_compute_binvrow (mpq_lpinfo * lp,
                                                    mpq_svector * zz,
                                                    int row,
                                                    mpq_t ztoler),
  mpq_ILLfct_compute_vA (mpq_lpinfo * lp,
                                         mpq_svector * v,
                                         mpq_t * vA),
  mpq_ILLfct_compute_psteep_upv (mpq_lpinfo * lp,
                                                         mpq_svector * swz),
  mpq_ILLfct_compute_dsteep_upv (mpq_lpinfo * lp,
                                                         mpq_svector * swz),
  mpq_ILLfct_update_basis_info (mpq_lpinfo * lp,
                                                        int eindex,
                                                        int lindex,
                                                        int lvstat),
  mpq_ILLfct_update_xz (mpq_lpinfo * lp,
                                        mpq_t tz,
                                        int eindex,
                                        int lindex),
  mpq_ILLfct_update_piz (mpq_lpinfo * lp,
                                         mpq_t alpha),
  mpq_ILLfct_update_pIpiz (mpq_lpinfo * lp,
                                             mpq_svector * z,
                                             mpq_t alpha),
  mpq_ILLfct_update_dz (mpq_lpinfo * lp,
                                        int eindex,
                                        mpq_t alpha),
  mpq_ILLfct_update_pIdz (mpq_lpinfo * lp,
                                            mpq_svector * zA,
                                            int eindex,
                                            mpq_t alpha),
    mpq_ILLfct_unroll_bound_change (mpq_lpinfo * lp),
    mpq_ILLfct_unroll_coef_change (mpq_lpinfo * lp),
    mpq_ILLfct_check_pfeasible (mpq_lpinfo * lp, mpq_feas_info * fs,
        mpq_t ftol),
    mpq_ILLfct_check_pIpfeasible (mpq_lpinfo * lp, mpq_feas_info * fs),
    mpq_ILLfct_check_dfeasible (mpq_lpinfo * lp, mpq_feas_info * fs),
    mpq_ILLfct_dual_adjust (mpq_lpinfo * lp),
    mpq_ILLfct_dphaseI_simple_update (mpq_lpinfo * lp),
    mpq_ILLfct_set_status_values (mpq_lpinfo * lp, int pstatus, int dstatus,
        int ptype, int dtype),
    mpq_ILLfct_init_counts (mpq_lpinfo * lp),
    mpq_ILLfct_update_counts (mpq_lpinfo * lp, int f, int upi, mpq_t upd),
    mpq_ILLfct_print_counts (mpq_lpinfo * lp),
    mpq_ILLfct_check_pIdfeasible (mpq_lpinfo * lp, mpq_feas_info * fs),
  mpq_ILLfct_update_pfeas (mpq_lpinfo * lp,
                                             int lindex,
                                             mpq_svector * srhs),
  mpq_ILLfct_compute_ppIzz (mpq_lpinfo * lp,
                                                mpq_svector * srhs,
                                                mpq_svector * ssoln),
  mpq_ILLfct_update_ppI_prices (mpq_lpinfo * lp,
                                                        mpq_price_info * pinf,
                                                        mpq_svector * srhs,
                                                        mpq_svector * ssoln,
                                                        int eindex,
                                                        int lindex,
                                                        mpq_t alpha),
  mpq_ILLfct_update_dfeas (mpq_lpinfo * lp,
                                             int eindex,
                                             mpq_svector * srhs),
  mpq_ILLfct_compute_dpIy (mpq_lpinfo * lp,
                                             mpq_svector * srhs,
                                             mpq_svector * ssoln),
  mpq_ILLfct_update_dpI_prices (mpq_lpinfo * lp,
                                                        mpq_price_info * pinf,
                                                        mpq_svector * srhs,
                                                        mpq_svector * ssoln,
                                                        int lindex,
                                                        mpq_t alpha),
  mpq_ILLfct_update_dIIfeas (mpq_lpinfo * lp,
                                                 int eindex,
                                                 mpq_svector * srhs),
  mpq_ILLfct_compute_dpIIy (mpq_lpinfo * lp,
                                                mpq_svector * srhs,
                                                mpq_svector * ssoln),
    mpq_ILLfct_update_dpII_prices (mpq_lpinfo * lp, mpq_price_info * pinf,
        mpq_svector * srhs, mpq_svector * ssoln, int lindex, mpq_t eval,
        mpq_t alpha);

void mpq_fct_test_workvector (mpq_lpinfo * lp),
  mpq_fct_test_pfeasible (mpq_lpinfo * lp),
  mpq_fct_test_dfeasible (mpq_lpinfo * lp),
  mpq_fct_test_pI_x (mpq_lpinfo * lp,
                                 mpq_price_info * p),
  mpq_fct_test_pII_x (mpq_lpinfo * lp,
                                    mpq_price_info * p),
  mpq_fct_test_pI_pi_dz (mpq_lpinfo * lp,
                                         mpq_price_info * p),
  mpq_fct_test_pII_pi_dz (mpq_lpinfo * lp,
                                            mpq_price_info * p);

mpq_bndinfo *mpq_ILLfct_new_bndinfo (void);
void mpq_ILLfct_free_bndinfo (mpq_bndinfo * binfo);

#endif /* mpq___FUNCTIONS_H */
