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

/*  $RCSfile: mpq_ratio.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef mpq___RATIO_H
#define mpq___RATIO_H
#include "basicdefs.h"
typedef struct mpq_ratio_res
{
    mpq_t tz;
    int eindex;
    int lindex;
    int lvstat;
    int ratio_stat;
    int boundch;
    int coeffch;
    mpq_t lbound;
    mpq_t ecoeff;
    mpq_t pivotval;
}
mpq_ratio_res;

void mpq_ILLratio_pI_test (mpq_lpinfo * const lp,
                                             int const eindex,
                                             int const dir,
                                             mpq_ratio_res * const rs),
  mpq_ILLratio_pII_test (mpq_lpinfo * const lp,
                                         int const eindex,
                                         int const dir,
                                         mpq_ratio_res * const rs),
  mpq_ILLratio_dI_test (mpq_lpinfo * const lp,
                                        int const lindex,
                                        int const lvstat,
                                        mpq_ratio_res * const rs),
  mpq_ILLratio_dII_test (mpq_lpinfo * const lp,
                                         int const lvstat,
                                         mpq_ratio_res * const rs),
  mpq_ILLratio_longdII_test (mpq_lpinfo * const lp,
                                                 int const lindex,
                                                 int const lvstat,
                                                 mpq_ratio_res * const rs),
  mpq_ILLratio_pivotin_test (mpq_lpinfo * const lp,
                                                 int *const rlist,
                                                 int const rcnt,
                                                 mpq_ratio_res * const rs);

#endif /* mpq___RATIO_H */
