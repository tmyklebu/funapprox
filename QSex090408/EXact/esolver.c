/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures
   and algorithms commons in solving MIP's

Copyright (C) 2008 David Applegate, Bill Cook, Sanjeeb Dash, Daniel Espinoza.
Sanjeeb Dash's ownership of copyright in QSopt_ex is derived from his
copyright in QSopt.


This library is free software; you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

You should have received a copy of the GNU Lesser General Public License
   along with this library; if not, write to the Free Software Foundation,
   Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA */
/* ========================================================================= */
#include "econfig.h"
#include "qs_all.h"

static char *fname = 0;
static int lpfile = 0;
static int usescaling = 1;
static int showversion = 0;
static int simplexalgo = PRIMAL_SIMPLEX;
static int pstrategy = QS_PRICE_PSTEEP;
static int dstrategy = QS_PRICE_DSTEEP;
static unsigned precision = 128;
static int printsol = 0;
static char *readbasis = 0;
static char *writebasis = 0;

/* ========================================================================= */
/* brief Display options to the screen                                       */
/* ========================================================================= */
static void usage (char *s)
{
    fprintf (stderr, "Usage: %s [- below -] prob_file\n", s);
    fprintf (stderr, "   -b f  write basis to file f\n");
    fprintf (stderr, "   -B f  read initial basis from file f\n");
#if 0
    fprintf (stderr, "   -I    solve the MIP using BestBound\n");
    fprintf (stderr, "   -E    edit problem after solving initial version\n");
#endif
    fprintf (stderr, "   -L    input file is in lp format (default: mps)\n");
    fprintf (stderr, "   -O    write the final solution to a .sol file\n");
    fprintf (stderr, "   -p #  run primal simplex with pricing rule #\n");
    fprintf (stderr,
	"         (%d-Dantzig, %d-Devex, %d-Steep (default), %d-Partial\n",
	QS_PRICE_PDANTZIG, QS_PRICE_PDEVEX, QS_PRICE_PSTEEP,
	QS_PRICE_PMULTPARTIAL);
    fprintf (stderr,
	"   -P #  number of bits to use for the float representation (default: 128)\n");
    fprintf (stderr, "   -d #  run dual simplex with pricing rule #\n");
    fprintf (stderr, "         (%d-Dantzig, %d-Steep, %d-Partial, %d-Devex)\n",
	QS_PRICE_DDANTZIG, QS_PRICE_DSTEEP, QS_PRICE_DMULTPARTIAL,
	QS_PRICE_DDEVEX);
    fprintf (stderr, "   -S    do NOT scale the initial LP\n");
    fprintf (stderr, "   -v    print QSopt version number\n");
}

/* ========================================================================= */
/* brief decide if a given file is mps or lp (only by extension)             */
/* ========================================================================= */
static void get_ftype (char *name, int *ftype)
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

/* ========================================================================= */
/* brief parssing options for the program                                    */
/* ========================================================================= */
static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = 0;

    while ((c =
	    ILLutil_bix_getopt (ac, av, "b:B:d:p:P:IELOSv", &boptind,
		&boptarg)) != EOF)
	switch (c) {
	case 'b':
	    writebasis = boptarg;
	    break;
	case 'B':
	    readbasis = boptarg;
	    break;
	case 'P':
	    precision = atoi (boptarg);
	    break;
	case 'd':
	    simplexalgo = DUAL_SIMPLEX;
	    dstrategy = atoi (boptarg);
	    break;
	case 'L':
	    lpfile = 1;
	    break;
	case 'O':
	    printsol = 1;
	    break;
	case 'p':
	    simplexalgo = PRIMAL_SIMPLEX;
	    pstrategy = atoi (boptarg);
	    break;
	case 'S':
	    usescaling = 0;
	    break;
	case 'v':
	    showversion = 1;
	    break;
	case '?':
	default:
	    usage (av[0]);
	    return 1;
	}
    if ((boptind == ac) && (showversion)) {
	char *buf = 0;
	buf = mpq_QSversion ();
	printf ("%s\n", buf);
	mpq_QSfree ((void *) buf);
	exit (0);
    }
    if (boptind != (ac - 1)) {
	usage (av[0]);
	return 1;
    }
    fname = av[boptind++];
    fprintf (stderr, "Reading problem from %s\n", fname);
    return 0;
}

/* ========================================================================= */
/* the main thing!                                                           */
/* ========================================================================= */
int main (int ac, char **av)
{
    int rval = 0, status = 0;
    mpq_QSdata *p_mpq = 0;
    QSbasis *basis = 0;
    ILLutil_timer timer_solve;
    ILLutil_timer timer_read;
    int ftype = 0;		/* 0 mps, 1 lp */
    mpq_t *y_mpq = 0, *x_mpq = 0;

    /* parse arguments and initialize EGlpNum related things */
    rval = parseargs (ac, av);
    QSexact_set_precision (precision);
    if (rval)
	goto CLEANUP;
    if (writebasis) {
	basis = EGsMalloc (QSbasis, 1);
	memset (basis, 0, sizeof (QSbasis));
    }
    /* just for the bell's and wistle */
    if (showversion) {
	char *buf = 0;
	buf = mpq_QSversion ();
	if (buf == 0) {
	    ILL_CLEANUP;
	} else {
	    printf ("%s\n", buf);
	    mpq_QSfree ((void *) buf);
	}
    }
    /* get the file type */
    if (lpfile)
	ftype = 1;
    else
	get_ftype (fname, &ftype);

    /* read the mpq problem */
    ILLutil_init_timer (&timer_read, "SOLVER_READ_MPQ");
    ILLutil_start_timer (&timer_read);
    if (ftype == 1) {
	p_mpq = mpq_QSread_prob ((const char *) fname, "LP");
	if (p_mpq == 0) {
	    fprintf (stderr, "Could not read lp file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    } else {
	p_mpq = mpq_QSread_prob ((const char *) fname, "MPS");
	if (p_mpq == 0) {
	    fprintf (stderr, "Could not read mps file.\n");
	    rval = 1;
	    ILL_CLEANUP_IF (rval);
	}
    }

    /* and get the basis if needed */
    if (readbasis) {
	rval = mpq_QSread_and_load_basis (p_mpq, (const char *) readbasis);
	ILL_CLEANUP_IF (rval);
	if (basis)
	    mpq_QSfree_basis (basis);
	basis = mpq_QSget_basis (p_mpq);
    }
    ILLutil_stop_timer (&timer_read, 1);
    /* set the readed flags */
    rval = mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_DISPLAY, 1)
	|| mpq_QSset_param (p_mpq, QS_PARAM_PRIMAL_PRICING, pstrategy)
	|| mpq_QSset_param (p_mpq, QS_PARAM_DUAL_PRICING, dstrategy)
	|| mpq_QSset_param (p_mpq, QS_PARAM_SIMPLEX_SCALING, usescaling);
    ILL_CLEANUP_IF (rval);
    if (printsol) {
	x_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->ncols);
	y_mpq = mpq_EGlpNumAllocArray (p_mpq->qslp->nrows);
    }
    ILLutil_init_timer (&timer_solve, "SOLVER");
    ILLutil_start_timer (&timer_solve);
    rval = QSexact_solver (p_mpq, x_mpq, y_mpq, basis, simplexalgo, &status);
    ILL_CLEANUP_IF (rval);
    ILLutil_stop_timer (&timer_solve, 1);
    if (printsol) {
	char out_f_name[100];
	FILE *out_f;
	sprintf (out_f_name, "%s.sol", p_mpq->qslp->probname);
	out_f = fopen (out_f_name, "w+");
	switch (status) {
	case QS_LP_OPTIMAL:
	    fprintf (out_f, "status = OPTIMAL\n");
	    break;
	case QS_LP_INFEASIBLE:
	    fprintf (out_f, "status = INFEASIBLE\n");
	    break;
	case QS_LP_UNBOUNDED:
	    fprintf (out_f, "status = UNBOUNDED\n");
	    break;
	default:
	    fprintf (out_f, "status = UNDEFINED\n");
	    break;
	}
	QSexact_print_sol (p_mpq, out_f);
	fclose (out_f);
    }
    /* ending */
CLEANUP:
    mpq_EGlpNumFreeArray (x_mpq);
    mpq_EGlpNumFreeArray (y_mpq);
    /* free the last allocated basis, and if we wanted to save it, do so */
    if (basis) {
	if (writebasis)
	    rval = mpq_QSwrite_basis (p_mpq, 0, writebasis);
    }
    mpq_QSfree_basis (basis);
    mpq_QSfree_prob (p_mpq);
    return rval;		/* main return */
}
