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

/*  $RCSfile: mpq_price.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef mpq___PRICE_H
#define mpq___PRICE_H

#include "mpq_dstruct.h"
#include "basicdefs.h"

typedef struct mpq_price_res
{
    int eindex;
    int dir;
    int lindex;
    int lvstat;
    int price_stat;
    mpq_t dinfeas;
    mpq_t pinfeas;
}
mpq_price_res;

int mpq_ILLprice_test_for_heap (mpq_lpinfo * const lp,
                                                        mpq_price_info * const pinf,
                                                        int const nkeys,
                                                        mpq_t * keylist,
                                                        int const algo,
                                                        int const upd),
  mpq_ILLprice_build_heap (mpq_price_info * const pinf,
                                             int const nkeys,
                                             mpq_t * keylist),
  mpq_ILLprice_build_pricing_info (mpq_lpinfo * const lp,
                                                             mpq_price_info * const pinf,
                                                             int const phase),
  mpq_ILLprice_update_pricing_info (mpq_lpinfo * const lp,
                                                                mpq_price_info * const pinf,
                                                                int const phase,
                                                                mpq_svector * const wz,
                                                                int const eindex,
                                                                int const lindex,
                                                                mpq_t y),
  mpq_ILLprice_get_price (mpq_price_info * const p,
                                            int const phase),
  mpq_ILLprice_build_mpartial_info (mpq_lpinfo * const lp,
                                                                mpq_price_info * const pinf,
                                                                int const pricetype),
  mpq_ILLprice_build_pdevex_norms (mpq_lpinfo * const lp,
                                                             mpq_p_devex_info * const pdinfo,
                                                             int const reinit),
  mpq_ILLprice_update_pdevex_norms (mpq_lpinfo * const lp,
                                                                mpq_p_devex_info * const pdinfo,
                                                                int const eindex,
                                                                mpq_t yl),
  mpq_ILLprice_build_psteep_norms (mpq_lpinfo * const lp,
                                                             mpq_p_steep_info * const psinfo),
  mpq_ILLprice_build_ddevex_norms (mpq_lpinfo * const lp,
                                                             mpq_d_devex_info * const ddinfo,
                                                             int const reinit),
  mpq_ILLprice_update_ddevex_norms (mpq_lpinfo * const lp,
                                                                mpq_d_devex_info * const ddinfo,
                                                                int const eindex,
                                                                mpq_t yl),
  mpq_ILLprice_build_dsteep_norms (mpq_lpinfo * const lp,
                                                             mpq_d_steep_info * const dsinfo),
  mpq_ILLprice_get_dsteep_norms (mpq_lpinfo * const lp,
                                                         int const count,
                                                         int *constrowind,
                                                         mpq_t * const norms),
  mpq_ILLprice_get_rownorms (mpq_lpinfo * const lp,
                                                 mpq_price_info * const pinf,
                                                 mpq_t * const rnorms),
  mpq_ILLprice_get_colnorms (mpq_lpinfo * const lp,
                                                 mpq_price_info * const pinf,
                                                 mpq_t * const cnorms),
  mpq_ILLprice_get_newnorms (mpq_lpinfo * const lp,
                                                 int const nelems,
                                                 mpq_t * const norms,
                                                 int *const matcnt,
                                                 int *const matbeg,
                                                 int *const matind,
                                                 mpq_t * const matval,
                                                 int const option),
  mpq_ILLprice_get_new_rownorms (mpq_lpinfo * const lp,
                                                         int const newrows,
                                                         mpq_t * const rnorms,
                                                         int *const rmatcnt,
                                                         int *const rmatbeg,
                                                         int *const rmatind,
                                                         mpq_t * const rmatval),
  mpq_ILLprice_get_new_colnorms (mpq_lpinfo * const lp,
                                                         int const newrows,
                                                         mpq_t * const rnorms,
                                                         int *const matcnt,
                                                         int *const matbeg,
                                                         int *const matind,
                                                         mpq_t * const matval),
  mpq_ILLprice_load_rownorms (mpq_lpinfo * const lp,
                                                    mpq_t * const rnorms,
                                                    mpq_price_info * const pinf),
  mpq_ILLprice_load_colnorms (mpq_lpinfo * const lp,
                                                    mpq_t * const cnorms,
                                                    mpq_price_info * const pinf);


void mpq_ILLprice_free_heap (mpq_price_info * const pinf),
  mpq_ILLprice_init_pricing_info (mpq_price_info * const pinf),
  mpq_ILLprice_free_pricing_info (mpq_price_info * const pinf),
  mpq_ILLprice_free_mpartial_info (mpq_mpart_info * p),
  mpq_ILLprice_init_mpartial_price (mpq_lpinfo * const lp,
                                                                mpq_price_info * const pinf,
                                                                int const phase,
                                                                int const pricetype),
  mpq_ILLprice_update_mpartial_price (mpq_lpinfo * const lp,
                                                                    mpq_price_info * const pinf,
                                                                    int const phase,
                                                                    int const pricetype),
    mpq_ILLprice_delete_onempart_price (mpq_price_info * const pinf,
        int const indx, int const pricetype),
  mpq_ILLprice_mpartial_group (mpq_lpinfo * const lp,
                                                     mpq_mpart_info * const p,
                                                     int const phase,
                                                     int const g,
                                                     int const pricetype),
  mpq_ILLprice_column (mpq_lpinfo * const lp,
                                     int const ix,
                                     int const phase,
                                     mpq_price_res * const pr),
  mpq_ILLprice_row (mpq_lpinfo * const lp,
                                int const ix,
                                int const phase,
                                mpq_price_res * const pr),
  mpq_ILLprice_update_psteep_norms (mpq_lpinfo * lp,
                                                                mpq_p_steep_info * psinfo,
                                                                mpq_svector * wz,
                                                                int eindex,
                                                                mpq_t yl),
  mpq_ILLprice_update_dsteep_norms (mpq_lpinfo * const lp,
                                                                mpq_d_steep_info * const dsinfo,
                                                                mpq_svector * const wz,
                                                                int const lindex,
                                                                mpq_t yl),
  mpq_ILLprice_compute_dual_inf (mpq_lpinfo * const lp,
                                                         mpq_price_info * const p,
                                                         int *const ix,
                                                         int const icnt,
                                                         int const phase),
  mpq_ILLprice_primal (mpq_lpinfo * const lp,
                                     mpq_price_info * const pinf,
                                     mpq_price_res * const pr,
                                     int const phase),
  mpq_ILLprice_compute_primal_inf (mpq_lpinfo * const lp,
                                                             mpq_price_info * const p,
                                                             int *const ix,
                                                             int const icnt,
                                                             int const phase),
  mpq_ILLprice_dual (mpq_lpinfo * const lp,
                                 mpq_price_info * const pinf,
                                 int const phase,
                                 mpq_price_res * const pr);

void mpq_test_dsteep_norms (mpq_lpinfo * const lp,
                                                mpq_price_info * const p);

#endif /* mpq___PRICE_H */
