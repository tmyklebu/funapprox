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

/*  $RCSfile: mpf_price.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef mpf___PRICE_H
#define mpf___PRICE_H

#include "mpf_dstruct.h"
#include "basicdefs.h"

typedef struct mpf_price_res
{
    int eindex;
    int dir;
    int lindex;
    int lvstat;
    int price_stat;
    mpf_t dinfeas;
    mpf_t pinfeas;
}
mpf_price_res;

int mpf_ILLprice_test_for_heap (mpf_lpinfo * const lp,
                                                        mpf_price_info * const pinf,
                                                        int const nkeys,
                                                        mpf_t * keylist,
                                                        int const algo,
                                                        int const upd),
  mpf_ILLprice_build_heap (mpf_price_info * const pinf,
                                             int const nkeys,
                                             mpf_t * keylist),
  mpf_ILLprice_build_pricing_info (mpf_lpinfo * const lp,
                                                             mpf_price_info * const pinf,
                                                             int const phase),
  mpf_ILLprice_update_pricing_info (mpf_lpinfo * const lp,
                                                                mpf_price_info * const pinf,
                                                                int const phase,
                                                                mpf_svector * const wz,
                                                                int const eindex,
                                                                int const lindex,
                                                                mpf_t y),
  mpf_ILLprice_get_price (mpf_price_info * const p,
                                            int const phase),
  mpf_ILLprice_build_mpartial_info (mpf_lpinfo * const lp,
                                                                mpf_price_info * const pinf,
                                                                int const pricetype),
  mpf_ILLprice_build_pdevex_norms (mpf_lpinfo * const lp,
                                                             mpf_p_devex_info * const pdinfo,
                                                             int const reinit),
  mpf_ILLprice_update_pdevex_norms (mpf_lpinfo * const lp,
                                                                mpf_p_devex_info * const pdinfo,
                                                                int const eindex,
                                                                mpf_t yl),
  mpf_ILLprice_build_psteep_norms (mpf_lpinfo * const lp,
                                                             mpf_p_steep_info * const psinfo),
  mpf_ILLprice_build_ddevex_norms (mpf_lpinfo * const lp,
                                                             mpf_d_devex_info * const ddinfo,
                                                             int const reinit),
  mpf_ILLprice_update_ddevex_norms (mpf_lpinfo * const lp,
                                                                mpf_d_devex_info * const ddinfo,
                                                                int const eindex,
                                                                mpf_t yl),
  mpf_ILLprice_build_dsteep_norms (mpf_lpinfo * const lp,
                                                             mpf_d_steep_info * const dsinfo),
  mpf_ILLprice_get_dsteep_norms (mpf_lpinfo * const lp,
                                                         int const count,
                                                         int *constrowind,
                                                         mpf_t * const norms),
  mpf_ILLprice_get_rownorms (mpf_lpinfo * const lp,
                                                 mpf_price_info * const pinf,
                                                 mpf_t * const rnorms),
  mpf_ILLprice_get_colnorms (mpf_lpinfo * const lp,
                                                 mpf_price_info * const pinf,
                                                 mpf_t * const cnorms),
  mpf_ILLprice_get_newnorms (mpf_lpinfo * const lp,
                                                 int const nelems,
                                                 mpf_t * const norms,
                                                 int *const matcnt,
                                                 int *const matbeg,
                                                 int *const matind,
                                                 mpf_t * const matval,
                                                 int const option),
  mpf_ILLprice_get_new_rownorms (mpf_lpinfo * const lp,
                                                         int const newrows,
                                                         mpf_t * const rnorms,
                                                         int *const rmatcnt,
                                                         int *const rmatbeg,
                                                         int *const rmatind,
                                                         mpf_t * const rmatval),
  mpf_ILLprice_get_new_colnorms (mpf_lpinfo * const lp,
                                                         int const newrows,
                                                         mpf_t * const rnorms,
                                                         int *const matcnt,
                                                         int *const matbeg,
                                                         int *const matind,
                                                         mpf_t * const matval),
  mpf_ILLprice_load_rownorms (mpf_lpinfo * const lp,
                                                    mpf_t * const rnorms,
                                                    mpf_price_info * const pinf),
  mpf_ILLprice_load_colnorms (mpf_lpinfo * const lp,
                                                    mpf_t * const cnorms,
                                                    mpf_price_info * const pinf);


void mpf_ILLprice_free_heap (mpf_price_info * const pinf),
  mpf_ILLprice_init_pricing_info (mpf_price_info * const pinf),
  mpf_ILLprice_free_pricing_info (mpf_price_info * const pinf),
  mpf_ILLprice_free_mpartial_info (mpf_mpart_info * p),
  mpf_ILLprice_init_mpartial_price (mpf_lpinfo * const lp,
                                                                mpf_price_info * const pinf,
                                                                int const phase,
                                                                int const pricetype),
  mpf_ILLprice_update_mpartial_price (mpf_lpinfo * const lp,
                                                                    mpf_price_info * const pinf,
                                                                    int const phase,
                                                                    int const pricetype),
    mpf_ILLprice_delete_onempart_price ( mpf_price_info * const pinf,
        int const indx, int const pricetype),
  mpf_ILLprice_mpartial_group (mpf_lpinfo * const lp,
                                                     mpf_mpart_info * const p,
                                                     int const phase,
                                                     int const g,
                                                     int const pricetype),
  mpf_ILLprice_column (mpf_lpinfo * const lp,
                                     int const ix,
                                     int const phase,
                                     mpf_price_res * const pr),
  mpf_ILLprice_row (mpf_lpinfo * const lp,
                                int const ix,
                                int const phase,
                                mpf_price_res * const pr),
  mpf_ILLprice_update_psteep_norms (mpf_lpinfo * lp,
                                                                mpf_p_steep_info * psinfo,
                                                                mpf_svector * wz,
                                                                int eindex,
                                                                mpf_t yl),
  mpf_ILLprice_update_dsteep_norms (mpf_lpinfo * const lp,
                                                                mpf_d_steep_info * const dsinfo,
                                                                mpf_svector * const wz,
                                                                int const lindex,
                                                                mpf_t yl),
  mpf_ILLprice_compute_dual_inf (mpf_lpinfo * const lp,
                                                         mpf_price_info * const p,
                                                         int *const ix,
                                                         int const icnt,
                                                         int const phase),
  mpf_ILLprice_primal (mpf_lpinfo * const lp,
                                     mpf_price_info * const pinf,
                                     mpf_price_res * const pr,
                                     int const phase),
  mpf_ILLprice_compute_primal_inf (mpf_lpinfo * const lp,
                                                             mpf_price_info * const p,
                                                             int *const ix,
                                                             int const icnt,
                                                             int const phase),
  mpf_ILLprice_dual (mpf_lpinfo * const lp,
                                 mpf_price_info * const pinf,
                                 int const phase,
                                 mpf_price_res * const pr);

void mpf_test_dsteep_norms (mpf_lpinfo * const lp,
                                                mpf_price_info * const p);

#endif /* mpf___PRICE_H */
