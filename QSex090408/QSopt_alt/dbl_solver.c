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

/* RCS_INFO = "$RCSfile: dbl_solver.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

#include "econfig.h"
#include "dbl_solver.h"
#include "dbl_iqsutil.h"
#include "dbl_lpdefs.h"		/* for PRIMAL_SIMPLEX */
#include "dbl_qstruct.h"
#include "dbl_qsopt.h"
#include "dbl_binary.h"
#include "dbl_editor.h"
#include "dbl_price.h"
#include "dbl_lib.h"		/* for dbl_ILLmip_binary_dfs */
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static char *dbl_fname = 0;
static int dbl_lpfile = 0;
static int dbl_solvemip = 0;
static int dbl_interactive = 0;
static int dbl_usescaling = 1;
static int dbl_showversion = 0;
static int dbl_simplexalgo = PRIMAL_SIMPLEX;
static int dbl_pstrategy = QS_PRICE_PSTEEP;
static int dbl_dstrategy = QS_PRICE_DSTEEP;
static unsigned dbl_precision = 128;
static int dbl_printsol = 0;
static char *dbl_readbasis = 0;
static char *dbl_writebasis = 0;

static void dbl_usage (char *s),
  dbl_get_ftype (char *name,
      int *ftype);
#ifdef TEST_FEEDBACK
static int dbl_feedback (FILE * dest,
      const char *str);
#endif
static int parseargs (int ac,
      char **av);

#ifndef WIN32
dbl_QSLIB_INTERFACE int main (int ac,
      char **av);
int main (int ac,
      char **av)
{
    int rval;
    rval = dbl_solver_main (ac, av);
    return rval;
}
#endif

dbl_QSLIB_INTERFACE int dbl_solver_main (int ac,
      char **av)
{
    int rval = 0;
    dbl_QSdata *p = 0;
    ILLutil_timer timer_solve;
    ILLutil_timer timer_read;
    int ftype = 0;		/* 0 mps, 1 lp */
    double *x = 0;
    double val;

    rval = parseargs (ac, av);
    if (rval)
	ILL_CLEANUP;

    dbl_QSset_precision (dbl_precision);
    dbl_EGlpNumInitVar (val);
    if (dbl_showversion) {
	char *buf = 0;
	buf = dbl_QSversion ();
	if (buf == 0) {
	    ILL_CLEANUP;
	} else {
	    printf ("%s\n", buf);
	    dbl_QSfree ((void *) buf);
	}
    }
    if (dbl_lpfile) {
	ftype = 1;
    } else {
	dbl_get_ftype (dbl_fname, &ftype);
    }

    ILLutil_init_timer (&timer_read, "SOLVER_READ");
    ILLutil_start_timer (&timer_read);
    if (ftype == 1) {
	p = dbl_QSread_prob ((const char *) dbl_fname, "LP");
	if (p == 0) {
	    fprintf (stderr, "Could not read lp file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	p = dbl_QSread_prob ((const char *) dbl_fname, "MPS");
	if (p == 0) {
	    fprintf (stderr, "Could not read mps file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    }

    if (dbl_readbasis) {
	rval = dbl_QSread_and_load_basis (p, (const char *) dbl_readbasis);
	ILL_CLEANUP_IF (rval);
    }
    ILLutil_stop_timer (&timer_read, 1);
    rval = dbl_QSset_param (p, QS_PARAM_SIMPLEX_DISPLAY, 1)
	|| dbl_QSset_param (p, QS_PARAM_PRIMAL_PRICING, dbl_pstrategy)
	|| dbl_QSset_param (p, QS_PARAM_DUAL_PRICING, dbl_dstrategy)
	|| dbl_QSset_param (p, QS_PARAM_SIMPLEX_SCALING, dbl_usescaling);
    ILL_CLEANUP_IF (rval);

    if (dbl_interactive) {
	dbl_ILLeditor_init ();
	dbl_ILLeditor (p);
	goto CLEANUP;
    }
    ILLutil_init_timer (&timer_solve, "SOLVER");
    ILLutil_start_timer (&timer_solve);

#ifdef TEST_FEEDBACK
    dbl_QSset_reporter (p, (void *) dbl_feedback, stdout);
#endif
    rval = dbl_ILLeditor_solve (p, dbl_simplexalgo);
    ILLutil_stop_timer (&timer_solve, 1);
    ILL_CLEANUP_IF (rval);

    if (dbl_printsol)
	x = dbl_EGlpNumAllocArray (p->lp->O->nstruct);

    if (dbl_solvemip) {
	rval = dbl_ILLmip_bfs (p->lp, &val, x);
	ILL_CLEANUP_IF (rval);
	printf ("MIP Objective Value: %.6f\n", dbl_EGlpNumToLf (val));
	fflush (stdout);
	if (dbl_printsol && dbl_EGlpNumIsNeqq (val, dbl_ILL_MAXDOUBLE) &&
	    dbl_EGlpNumIsNeqq (val, dbl_ILL_MINDOUBLE)) {
	    rval = dbl_ILLlib_print_x (stdout, p->lp, 0, x, 1);
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	if (dbl_writebasis) {
	    rval = dbl_QSwrite_basis (p, 0, dbl_writebasis);
	    ILL_CLEANUP_IF (rval);
	}
	if (dbl_printsol) {
	    rval = dbl_ILLlib_print_x (stdout, p->lp, 0, 0, 1);
	    ILL_CLEANUP_IF (rval);
	}
    }

CLEANUP:
    dbl_EGlpNumFreeArray (x);
    dbl_QSfree_prob (p);
    dbl_EGlpNumClearVar (val);
    return rval;		/* main return */
}

static void dbl_usage (char *s)
{
    char *buf = 0;

    buf = dbl_QSversion ();
    if (buf) {
	fprintf (stderr, "%s\n", buf);
	dbl_QSfree ((void *) buf);
    }
    fprintf (stderr, "Usage: %s [- below -] prob_file\n", s);
    fprintf (stderr, "   -b f  write basis to file f\n");
    fprintf (stderr, "   -B f  read initial basis from file f\n");
#if 0
    fprintf (stderr, "   -I    solve the MIP using BestBound\n");
    fprintf (stderr, "   -E    edit problem after solving initial version\n");
#endif
    fprintf (stderr, "   -L    input file is in lp format (default: mps)\n");
    fprintf (stderr, "   -O    print the final solution\n");
    fprintf (stderr, "   -p #  run primal simplex with pricing rule #\n");
    fprintf (stderr,
	"         (%d-Dantzig, %d-Devex, %d-Steep (default), %d-Partial\n",
	QS_PRICE_PDANTZIG, QS_PRICE_PDEVEX, QS_PRICE_PSTEEP,
	QS_PRICE_PMULTPARTIAL);
#if EGLPNUM_TYPE == GNU_MP_F
    fprintf (stderr,
	"   -P #  number of bits to use for the float representation (default: 128)\n");
#endif
    fprintf (stderr, "   -d #  run dual simplex with pricing rule #\n");
    fprintf (stderr, "         (%d-Dantzig, %d-Steep, %d-Partial, %d-Devex)\n",
	QS_PRICE_DDANTZIG, QS_PRICE_DSTEEP, QS_PRICE_DMULTPARTIAL,
	QS_PRICE_DDEVEX);
    fprintf (stderr, "   -S    do NOT scale the initial LP\n");
    fprintf (stderr, "   -v    print QSopt version number\n");
}

static int parseargs (int ac,
      char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = 0;

    while ((c =
	    ILLutil_bix_getopt (ac, av, "b:B:d:p:P:IELOSv", &boptind,
		&boptarg)) != EOF)
	switch (c) {
	case 'b':
	    dbl_writebasis = boptarg;
	    break;
	case 'B':
	    dbl_readbasis = boptarg;
	    break;
	case 'P':
	    dbl_precision = atoi (boptarg);
	    break;
	case 'd':
	    dbl_simplexalgo = DUAL_SIMPLEX;
	    dbl_dstrategy = atoi (boptarg);
	    break;
#if 0
	case 'E':
	    dbl_interactive = 1;
	    break;
	case 'I':
	    dbl_solvemip = 1;
	    break;
#endif
	case 'L':
	    dbl_lpfile = 1;
	    break;
	case 'O':
	    dbl_printsol = 1;
	    break;
	case 'p':
	    dbl_simplexalgo = PRIMAL_SIMPLEX;
	    dbl_pstrategy = atoi (boptarg);
	    break;
	case 'S':
	    dbl_usescaling = 0;
	    break;
	case 'v':
	    dbl_showversion = 1;
	    break;
	case '?':
	default:
	    dbl_usage (av[0]);
	    return 1;
	}
    if (boptind != (ac - 1)) {
	dbl_usage (av[0]);
	return 1;
    }
    dbl_fname = av[boptind++];
    fprintf (stderr, "Reading problem from %s\n", dbl_fname);
    return 0;
}

static void dbl_get_ftype (char *name,
      int *ftype)
{
    char *q;

    q = strrchr (name, '.');
    if (q) {
	q++;
	if (!strcmp (q, "lp") || !strcmp (q, "LP")) {
	    *ftype = 1;
	} else {
	    *ftype = 0;
	}
    }
}

#ifdef TEST_FEEDBACK
static int dbl_feedback (FILE * dest,
      const char *str)
{
    if (str != NULL) {
	int rc = fprintf ((FILE *) dest, "FEEDBACK: %s", str);
	fflush (dest);
	return rc;
    }
    return 0;
}
#endif
