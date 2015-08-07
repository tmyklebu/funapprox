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

/* RCS_INFO = "$RCSfile: mpq_solver.c,v $ $Revision: 1.2 $ $Date: 2003/11/05
   16:49:52 $"; */

#include "econfig.h"
#include "mpq_solver.h"
#include "mpq_iqsutil.h"
#include "mpq_lpdefs.h"		/* for PRIMAL_SIMPLEX */
#include "mpq_qstruct.h"
#include "mpq_qsopt.h"
#include "mpq_binary.h"
#include "mpq_editor.h"
#include "mpq_price.h"
#include "mpq_lib.h"		/* for mpq_ILLmip_binary_dfs */
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static char *mpq_fname = 0;
static int mpq_lpfile = 0;
static int mpq_solvemip = 0;
static int mpq_interactive = 0;
static int mpq_usescaling = 1;
static int mpq_showversion = 0;
static int mpq_simplexalgo = PRIMAL_SIMPLEX;
static int mpq_pstrategy = QS_PRICE_PSTEEP;
static int mpq_dstrategy = QS_PRICE_DSTEEP;
static unsigned mpq_precision = 128;
static int mpq_printsol = 0;
static char *mpq_readbasis = 0;
static char *mpq_writebasis = 0;

static void mpq_usage (char *s),
  mpq_get_ftype (char *name,
      int *ftype);
#ifdef TEST_FEEDBACK
static int mpq_feedback (FILE * dest,
      const char *str);
#endif
static int parseargs (int ac,
      char **av);

#ifndef WIN32
mpq_QSLIB_INTERFACE int main (int ac,
      char **av);
int main (int ac,
      char **av)
{
    int rval;
    rval = mpq_solver_main (ac, av);
    return rval;
}
#endif

mpq_QSLIB_INTERFACE int mpq_solver_main (int ac,
      char **av)
{
    int rval = 0;
    mpq_QSdata *p = 0;
    ILLutil_timer timer_solve;
    ILLutil_timer timer_read;
    int ftype = 0;		/* 0 mps, 1 lp */
    mpq_t *x = 0;
    mpq_t val;

    rval = parseargs (ac, av);
    if (rval)
	ILL_CLEANUP;

    mpq_QSset_precision (mpq_precision);
    mpq_EGlpNumInitVar (val);
    if (mpq_showversion) {
	char *buf = 0;
	buf = mpq_QSversion ();
	if (buf == 0) {
	    ILL_CLEANUP;
	} else {
	    printf ("%s\n", buf);
	    mpq_QSfree ((void *) buf);
	}
    }
    if (mpq_lpfile) {
	ftype = 1;
    } else {
	mpq_get_ftype (mpq_fname, &ftype);
    }

    ILLutil_init_timer (&timer_read, "SOLVER_READ");
    ILLutil_start_timer (&timer_read);
    if (ftype == 1) {
	p = mpq_QSread_prob ((const char *) mpq_fname, "LP");
	if (p == 0) {
	    fprintf (stderr, "Could not read lp file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	p = mpq_QSread_prob ((const char *) mpq_fname, "MPS");
	if (p == 0) {
	    fprintf (stderr, "Could not read mps file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    }

    if (mpq_readbasis) {
	rval = mpq_QSread_and_load_basis (p, (const char *) mpq_readbasis);
	ILL_CLEANUP_IF (rval);
    }
    ILLutil_stop_timer (&timer_read, 1);
    rval = mpq_QSset_param (p, QS_PARAM_SIMPLEX_DISPLAY, 1)
	|| mpq_QSset_param (p, QS_PARAM_PRIMAL_PRICING, mpq_pstrategy)
	|| mpq_QSset_param (p, QS_PARAM_DUAL_PRICING, mpq_dstrategy)
	|| mpq_QSset_param (p, QS_PARAM_SIMPLEX_SCALING, mpq_usescaling);
    ILL_CLEANUP_IF (rval);

    if (mpq_interactive) {
	mpq_ILLeditor_init ();
	mpq_ILLeditor (p);
	goto CLEANUP;
    }
    ILLutil_init_timer (&timer_solve, "SOLVER");
    ILLutil_start_timer (&timer_solve);

#ifdef TEST_FEEDBACK
    mpq_QSset_reporter (p, (void *) mpq_feedback, stdout);
#endif
    rval = mpq_ILLeditor_solve (p, mpq_simplexalgo);
    ILLutil_stop_timer (&timer_solve, 1);
    ILL_CLEANUP_IF (rval);

    if (mpq_printsol)
	x = mpq_EGlpNumAllocArray (p->lp->O->nstruct);

    if (mpq_solvemip) {
	rval = mpq_ILLmip_bfs (p->lp, &val, x);
	ILL_CLEANUP_IF (rval);
	printf ("MIP Objective Value: %.6f\n", mpq_EGlpNumToLf (val));
	fflush (stdout);
	if (mpq_printsol && mpq_EGlpNumIsNeqq (val, mpq_ILL_MAXDOUBLE) &&
	    mpq_EGlpNumIsNeqq (val, mpq_ILL_MINDOUBLE)) {
	    rval = mpq_ILLlib_print_x (stdout, p->lp, 0, x, 1);
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	if (mpq_writebasis) {
	    rval = mpq_QSwrite_basis (p, 0, mpq_writebasis);
	    ILL_CLEANUP_IF (rval);
	}
	if (mpq_printsol) {
	    rval = mpq_ILLlib_print_x (stdout, p->lp, 0, 0, 1);
	    ILL_CLEANUP_IF (rval);
	}
    }

CLEANUP:
    mpq_EGlpNumFreeArray (x);
    mpq_QSfree_prob (p);
    mpq_EGlpNumClearVar (val);
    return rval;		/* main return */
}

static void mpq_usage (char *s)
{
    char *buf = 0;

    buf = mpq_QSversion ();
    if (buf) {
	fprintf (stderr, "%s\n", buf);
	mpq_QSfree ((void *) buf);
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
	    mpq_writebasis = boptarg;
	    break;
	case 'B':
	    mpq_readbasis = boptarg;
	    break;
	case 'P':
	    mpq_precision = atoi (boptarg);
	    break;
	case 'd':
	    mpq_simplexalgo = DUAL_SIMPLEX;
	    mpq_dstrategy = atoi (boptarg);
	    break;
#if 0
	case 'E':
	    mpq_interactive = 1;
	    break;
	case 'I':
	    mpq_solvemip = 1;
	    break;
#endif
	case 'L':
	    mpq_lpfile = 1;
	    break;
	case 'O':
	    mpq_printsol = 1;
	    break;
	case 'p':
	    mpq_simplexalgo = PRIMAL_SIMPLEX;
	    mpq_pstrategy = atoi (boptarg);
	    break;
	case 'S':
	    mpq_usescaling = 0;
	    break;
	case 'v':
	    mpq_showversion = 1;
	    break;
	case '?':
	default:
	    mpq_usage (av[0]);
	    return 1;
	}
    if (boptind != (ac - 1)) {
	mpq_usage (av[0]);
	return 1;
    }
    mpq_fname = av[boptind++];
    fprintf (stderr, "Reading problem from %s\n", mpq_fname);
    return 0;
}

static void mpq_get_ftype (char *name,
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
static int mpq_feedback (FILE * dest,
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
