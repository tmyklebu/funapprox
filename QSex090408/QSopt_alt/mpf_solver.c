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

/* RCS_INFO = "$RCSfile: mpf_solver.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

#include "econfig.h"
#include "mpf_solver.h"
#include "mpf_iqsutil.h"
#include "mpf_lpdefs.h"		/* for PRIMAL_SIMPLEX */
#include "mpf_qstruct.h"
#include "mpf_qsopt.h"
#include "mpf_binary.h"
#include "mpf_editor.h"
#include "mpf_price.h"
#include "mpf_lib.h"		/* for mpf_ILLmip_binary_dfs */
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static char *mpf_fname = 0;
static int mpf_lpfile = 0;
static int mpf_solvemip = 0;
static int mpf_interactive = 0;
static int mpf_usescaling = 1;
static int mpf_showversion = 0;
static int mpf_simplexalgo = PRIMAL_SIMPLEX;
static int mpf_pstrategy = QS_PRICE_PSTEEP;
static int mpf_dstrategy = QS_PRICE_DSTEEP;
static unsigned mpf_precision = 128;
static int mpf_printsol = 0;
static char *mpf_readbasis = 0;
static char *mpf_writebasis = 0;

static void mpf_usage (char *s),
  mpf_get_ftype (char *name,
      int *ftype);
#ifdef TEST_FEEDBACK
static int mpf_feedback (FILE * dest,
      const char *str);
#endif
static int parseargs (int ac,
      char **av);

#ifndef WIN32
mpf_QSLIB_INTERFACE int main (int ac,
      char **av);
int main (int ac,
      char **av)
{
    int rval;
    rval = mpf_solver_main (ac, av);
    return rval;
}
#endif

mpf_QSLIB_INTERFACE int mpf_solver_main (int ac,
      char **av)
{
    int rval = 0;
    mpf_QSdata *p = 0;
    ILLutil_timer timer_solve;
    ILLutil_timer timer_read;
    int ftype = 0;		/* 0 mps, 1 lp */
    mpf_t *x = 0;
    mpf_t val;

    rval = parseargs (ac, av);
    if (rval)
	ILL_CLEANUP;

    mpf_QSset_precision (mpf_precision);
    mpf_EGlpNumInitVar (val);
    if (mpf_showversion) {
	char *buf = 0;
	buf = mpf_QSversion ();
	if (buf == 0) {
	    ILL_CLEANUP;
	} else {
	    printf ("%s\n", buf);
	    mpf_QSfree ((void *) buf);
	}
    }
    if (mpf_lpfile) {
	ftype = 1;
    } else {
	mpf_get_ftype (mpf_fname, &ftype);
    }

    ILLutil_init_timer (&timer_read, "SOLVER_READ");
    ILLutil_start_timer (&timer_read);
    if (ftype == 1) {
	p = mpf_QSread_prob ((const char *) mpf_fname, "LP");
	if (p == 0) {
	    fprintf (stderr, "Could not read lp file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	p = mpf_QSread_prob ((const char *) mpf_fname, "MPS");
	if (p == 0) {
	    fprintf (stderr, "Could not read mps file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    }

    if (mpf_readbasis) {
	rval = mpf_QSread_and_load_basis (p, (const char *) mpf_readbasis);
	ILL_CLEANUP_IF (rval);
    }
    ILLutil_stop_timer (&timer_read, 1);
    rval = mpf_QSset_param (p, QS_PARAM_SIMPLEX_DISPLAY, 1)
	|| mpf_QSset_param (p, QS_PARAM_PRIMAL_PRICING, mpf_pstrategy)
	|| mpf_QSset_param (p, QS_PARAM_DUAL_PRICING, mpf_dstrategy)
	|| mpf_QSset_param (p, QS_PARAM_SIMPLEX_SCALING, mpf_usescaling);
    ILL_CLEANUP_IF (rval);

    if (mpf_interactive) {
	mpf_ILLeditor_init ();
	mpf_ILLeditor (p);
	goto CLEANUP;
    }
    ILLutil_init_timer (&timer_solve, "SOLVER");
    ILLutil_start_timer (&timer_solve);

#ifdef TEST_FEEDBACK
    mpf_QSset_reporter (p, (void *) mpf_feedback, stdout);
#endif
    rval = mpf_ILLeditor_solve (p, mpf_simplexalgo);
    ILLutil_stop_timer (&timer_solve, 1);
    ILL_CLEANUP_IF (rval);

    if (mpf_printsol)
	x = mpf_EGlpNumAllocArray (p->lp->O->nstruct);

    if (mpf_solvemip) {
	rval = mpf_ILLmip_bfs (p->lp, &val, x);
	ILL_CLEANUP_IF (rval);
	printf ("MIP Objective Value: %.6f\n", mpf_EGlpNumToLf (val));
	fflush (stdout);
	if (mpf_printsol && mpf_EGlpNumIsNeqq (val, mpf_ILL_MAXDOUBLE) &&
	    mpf_EGlpNumIsNeqq (val, mpf_ILL_MINDOUBLE)) {
	    rval = mpf_ILLlib_print_x (stdout, p->lp, 0, x, 1);
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	if (mpf_writebasis) {
	    rval = mpf_QSwrite_basis (p, 0, mpf_writebasis);
	    ILL_CLEANUP_IF (rval);
	}
	if (mpf_printsol) {
	    rval = mpf_ILLlib_print_x (stdout, p->lp, 0, 0, 1);
	    ILL_CLEANUP_IF (rval);
	}
    }

CLEANUP:
    mpf_EGlpNumFreeArray (x);
    mpf_QSfree_prob (p);
    mpf_EGlpNumClearVar (val);
    return rval;		/* main return */
}

static void mpf_usage (char *s)
{
    char *buf = 0;

    buf = mpf_QSversion ();
    if (buf) {
	fprintf (stderr, "%s\n", buf);
	mpf_QSfree ((void *) buf);
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
	    mpf_writebasis = boptarg;
	    break;
	case 'B':
	    mpf_readbasis = boptarg;
	    break;
	case 'P':
	    mpf_precision = atoi (boptarg);
	    break;
	case 'd':
	    mpf_simplexalgo = DUAL_SIMPLEX;
	    mpf_dstrategy = atoi (boptarg);
	    break;
#if 0
	case 'E':
	    mpf_interactive = 1;
	    break;
	case 'I':
	    mpf_solvemip = 1;
	    break;
#endif
	case 'L':
	    mpf_lpfile = 1;
	    break;
	case 'O':
	    mpf_printsol = 1;
	    break;
	case 'p':
	    mpf_simplexalgo = PRIMAL_SIMPLEX;
	    mpf_pstrategy = atoi (boptarg);
	    break;
	case 'S':
	    mpf_usescaling = 0;
	    break;
	case 'v':
	    mpf_showversion = 1;
	    break;
	case '?':
	default:
	    mpf_usage (av[0]);
	    return 1;
	}
    if (boptind != (ac - 1)) {
	mpf_usage (av[0]);
	return 1;
    }
    mpf_fname = av[boptind++];
    fprintf (stderr, "Reading problem from %s\n", mpf_fname);
    return 0;
}

static void mpf_get_ftype (char *name,
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
static int mpf_feedback (FILE * dest,
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
