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

/*  $RCSfile: mpf_ratio.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef mpf___RATIO_H
#define mpf___RATIO_H
#include "basicdefs.h"
typedef struct mpf_ratio_res
{
    mpf_t tz;
    int eindex;
    int lindex;
    int lvstat;
    int ratio_stat;
    int boundch;
    int coeffch;
    mpf_t lbound;
    mpf_t ecoeff;
    mpf_t pivotval;
}
mpf_ratio_res;

void mpf_ILLratio_pI_test (mpf_lpinfo * const lp,
                                             int const eindex,
                                             int const dir,
                                             mpf_ratio_res * const rs),
  mpf_ILLratio_pII_test (mpf_lpinfo * const lp,
                                         int const eindex,
                                         int const dir,
                                         mpf_ratio_res * const rs),
  mpf_ILLratio_dI_test (mpf_lpinfo * const lp,
                                        int const lindex,
                                        int const lvstat,
                                        mpf_ratio_res * const rs),
    mpf_ILLratio_dII_test (mpf_lpinfo * const lp, int const lvstat,
        mpf_ratio_res * const rs),
  mpf_ILLratio_longdII_test (mpf_lpinfo * const lp,
                                                 int const lindex,
                                                 int const lvstat,
                                                 mpf_ratio_res * const rs),
  mpf_ILLratio_pivotin_test (mpf_lpinfo * const lp,
                                                 int *const rlist,
                                                 int const rcnt,
                                                 mpf_ratio_res * const rs);

#endif /* mpf___RATIO_H */
