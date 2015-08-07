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

/*  $RCSfile: dbl_price.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef dbl___PRICE_H
#define dbl___PRICE_H

#include "dbl_dstruct.h"
#include "basicdefs.h"

typedef struct dbl_price_res
{
    int eindex;
    int dir;
    int lindex;
    int lvstat;
    int price_stat;
    double dinfeas;
    double pinfeas;
}
dbl_price_res;

int dbl_ILLprice_test_for_heap (dbl_lpinfo * const lp,
                                                        dbl_price_info * const pinf,
                                                        int const nkeys,
                                                        double * keylist,
                                                        int const algo,
                                                        int const upd),
  dbl_ILLprice_build_heap (dbl_price_info * const pinf,
                                             int const nkeys,
                                             double * keylist),
  dbl_ILLprice_build_pricing_info (dbl_lpinfo * const lp,
                                                             dbl_price_info * const pinf,
                                                             int const phase),
  dbl_ILLprice_update_pricing_info (dbl_lpinfo * const lp,
                                                                dbl_price_info * const pinf,
                                                                int const phase,
                                                                dbl_svector * const wz,
                                                                int const eindex,
                                                                int const lindex,
                                                                double y),
  dbl_ILLprice_get_price (dbl_price_info * const p,
                                            int const phase),
  dbl_ILLprice_build_mpartial_info (dbl_lpinfo * const lp,
                                                                dbl_price_info * const pinf,
                                                                int const pricetype),
  dbl_ILLprice_build_pdevex_norms (dbl_lpinfo * const lp,
                                                             dbl_p_devex_info * const pdinfo,
                                                             int const reinit),
  dbl_ILLprice_update_pdevex_norms (dbl_lpinfo * const lp,
                                                                dbl_p_devex_info * const pdinfo,
                                                                int const eindex,
                                                                double yl),
  dbl_ILLprice_build_psteep_norms (dbl_lpinfo * const lp,
                                                             dbl_p_steep_info * const psinfo),
  dbl_ILLprice_build_ddevex_norms (dbl_lpinfo * const lp,
                                                             dbl_d_devex_info * const ddinfo,
                                                             int const reinit),
  dbl_ILLprice_update_ddevex_norms (dbl_lpinfo * const lp,
                                                                dbl_d_devex_info * const ddinfo,
                                                                int const eindex,
                                                                double yl),
  dbl_ILLprice_build_dsteep_norms (dbl_lpinfo * const lp,
                                                             dbl_d_steep_info * const dsinfo),
  dbl_ILLprice_get_dsteep_norms (dbl_lpinfo * const lp,
                                                         int const count,
                                                         int *constrowind,
                                                         double * const norms),
  dbl_ILLprice_get_rownorms (dbl_lpinfo * const lp,
                                                 dbl_price_info * const pinf,
                                                 double * const rnorms),
  dbl_ILLprice_get_colnorms (dbl_lpinfo * const lp,
                                                 dbl_price_info * const pinf,
                                                 double * const cnorms),
  dbl_ILLprice_get_newnorms (dbl_lpinfo * const lp,
                                                 int const nelems,
                                                 double * const norms,
                                                 int *const matcnt,
                                                 int *const matbeg,
                                                 int *const matind,
                                                 double * const matval,
                                                 int const option),
  dbl_ILLprice_get_new_rownorms (dbl_lpinfo * const lp,
                                                         int const newrows,
                                                         double * const rnorms,
                                                         int *const rmatcnt,
                                                         int *const rmatbeg,
                                                         int *const rmatind,
                                                         double * const rmatval),
  dbl_ILLprice_get_new_colnorms (dbl_lpinfo * const lp,
                                                         int const newrows,
                                                         double * const rnorms,
                                                         int *const matcnt,
                                                         int *const matbeg,
                                                         int *const matind,
                                                         double * const matval),
  dbl_ILLprice_load_rownorms (dbl_lpinfo * const lp,
                                                    double * const rnorms,
                                                    dbl_price_info * const pinf),
  dbl_ILLprice_load_colnorms (dbl_lpinfo * const lp,
                                                    double * const cnorms,
                                                    dbl_price_info * const pinf);


void dbl_ILLprice_free_heap (dbl_price_info * const pinf),
  dbl_ILLprice_init_pricing_info (dbl_price_info * const pinf),
  dbl_ILLprice_free_pricing_info (dbl_price_info * const pinf),
  dbl_ILLprice_free_mpartial_info (dbl_mpart_info * p),
  dbl_ILLprice_init_mpartial_price (dbl_lpinfo * const lp,
                                                                dbl_price_info * const pinf,
                                                                int const phase,
                                                                int const pricetype),
  dbl_ILLprice_update_mpartial_price (dbl_lpinfo * const lp,
                                                                    dbl_price_info * const pinf,
                                                                    int const phase,
                                                                    int const pricetype),
    dbl_ILLprice_delete_onempart_price (dbl_price_info * const pinf,
        int const indx, int const pricetype),
  dbl_ILLprice_mpartial_group (dbl_lpinfo * const lp,
                                                     dbl_mpart_info * const p,
                                                     int const phase,
                                                     int const g,
                                                     int const pricetype),
  dbl_ILLprice_column (dbl_lpinfo * const lp,
                                     int const ix,
                                     int const phase,
                                     dbl_price_res * const pr),
  dbl_ILLprice_row (dbl_lpinfo * const lp,
                                int const ix,
                                int const phase,
                                dbl_price_res * const pr),
  dbl_ILLprice_update_psteep_norms (dbl_lpinfo * lp,
                                                                dbl_p_steep_info * psinfo,
                                                                dbl_svector * wz,
                                                                int eindex,
                                                                double yl),
  dbl_ILLprice_update_dsteep_norms (dbl_lpinfo * const lp,
                                                                dbl_d_steep_info * const dsinfo,
                                                                dbl_svector * const wz,
                                                                int const lindex,
                                                                double yl),
  dbl_ILLprice_compute_dual_inf (dbl_lpinfo * const lp,
                                                         dbl_price_info * const p,
                                                         int *const ix,
                                                         int const icnt,
                                                         int const phase),
  dbl_ILLprice_primal (dbl_lpinfo * const lp,
                                     dbl_price_info * const pinf,
                                     dbl_price_res * const pr,
                                     int const phase),
  dbl_ILLprice_compute_primal_inf (dbl_lpinfo * const lp,
                                                             dbl_price_info * const p,
                                                             int *const ix,
                                                             int const icnt,
                                                             int const phase),
dbl_ILLprice_dual (dbl_lpinfo * const lp, dbl_price_info * const pinf,
    int const phase, dbl_price_res * const pr);

void dbl_test_dsteep_norms (dbl_lpinfo * const lp,
                                                dbl_price_info * const p);

#endif /* dbl___PRICE_H */
