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

/* RCS_INFO = "$RCSfile: format_error.c,v $ $Revision: 1.2 $ $Date:
   2003/11/05 16:49:52 $"; */

#include "econfig.h"
#include "mpq_format.h"
#include "mpq_qsopt.h"
#include "mpq_iqsutil.h"

int mpq_ILLformat_error_create (mpq_qsformat_error * error,
      int mode,
      const char *desc,
      int lineNum,
      const char *theLine,
      int atPos)
{
    int len;
    int rval = 0;

    error->theLine = NULL;
    error->desc = NULL;
    error->next = NULL;

    ILL_FAILtrue (desc == NULL, "non empty error desc please");
    ILL_FAILtrue (mode >= QS_INPUT_NERROR
	|| mode < 0, "0<= mode <=QS_INPUT_NERROR");
    error->type = mode;
    len = strlen (desc);
    ILL_SAFE_MALLOC (error->desc, len + 1, char);
    strcpy (error->desc, desc);
    error->lineNumber = lineNum;
    if (theLine != NULL) {
	len = strlen (theLine);
	ILL_SAFE_MALLOC (error->theLine, len + 2, char);
	strcpy (error->theLine, theLine);
	if (error->theLine[len - 1] != '\n') {
	    error->theLine[len] = '\n';
	    error->theLine[len + 1] = '\0';
	}
    }
    error->at = atPos;
CLEANUP:
    if (rval) {
	mpq_ILLformat_error_delete (error);
    }
    return rval;
}

void mpq_ILLformat_error_delete (mpq_qsformat_error * error)
{
    ILL_IFFREE (error->desc, char);
    ILL_IFFREE (error->theLine, char);
}

void mpq_ILLformat_error_print (FILE * out,
      mpq_qsformat_error * error)
{
    int at = error->at;
    int tp = error->type;
    const char *type = "Error";
    const char *line = NULL;
    int i;

    type = mpq_QSformat_error_type_string (tp);

    fprintf (out, "%s  line %d pos %d\n",
	type, mpq_QSerror_get_line_number (error), at);
    line = mpq_QSerror_get_line (error);
    if (line != NULL) {
	fprintf (out, "LINE %s", line);
	if (at >= 0) {
	    fprintf (out, ".....");
	    for (i = 0; i <= (at - 1); i++) {
		if (line[i] == '\t') {
		    fputc ('\t', out);
		} else {
		    fputc ('.', out);
		}
	    }
	    fprintf (out, "^\n");
	}
    } else {
	fprintf (out, "NO LINE\n");
    }
    fprintf (out, "MSG: %s\n", mpq_QSerror_get_desc (error));
}


mpq_qserror_collector *mpq_ILLerror_collector_new (mpq_qsadd_error_fct fct,
      void *dest)
{
    int rval = 0;
    mpq_qserror_collector *c = NULL;
    ILL_SAFE_MALLOC (c, 1, mpq_qserror_collector);
    c->add_error = fct;
    c->dest = dest;

CLEANUP:
    if (rval) {
	ILL_IFFREE (c, mpq_qserror_collector);
    }
    return c;
}

mpq_qserror_collector *mpq_ILLerror_memory_collector_new (mpq_qserror_memory * dest)
{
    return mpq_ILLerror_collector_new (mpq_ILLadd_error_to_memory, dest);
}

void mpq_ILLerror_collector_free (mpq_qserror_collector * c)
{
    ILL_IFFREE (c, mpq_qserror_collector);
}

mpq_qserror_memory *mpq_ILLerror_memory_create (int takeErrorLines)
{
    int rval = 0, i;
    mpq_qserror_memory *mem = NULL;
    ILL_SAFE_MALLOC (mem, 1, mpq_qserror_memory);
    for (i = 0; i < QS_INPUT_NERROR; i++) {
	mem->has_error[i] = 0;
    }
    mem->error_list = NULL;
    mem->nerror = 0;
    mem->hasErrorLines = takeErrorLines;
CLEANUP:
    return mem;
}

void mpq_ILLerror_memory_free (mpq_qserror_memory * mem)
{
    mpq_qsformat_error *ths, *nxt;

    if (mem != NULL) {
	ths = mem->error_list;
	while (ths != NULL) {
	    nxt = ths->next;
	    ILL_IFFREE (ths, mpq_qsformat_error);
	    ths = nxt;
	}
	ILL_IFFREE (mem, mpq_qserror_memory);
    }
}

int mpq_ILLadd_error_to_memory (void *dest,
      const mpq_qsformat_error * error)
{
    int rval = 0;
    mpq_qserror_memory *mem = (mpq_qserror_memory *) dest;
    mpq_qsformat_error *e = 0;
    ILL_CHECKnull (mem, "must give non NULL mpq_qserror_memory");

    ILL_SAFE_MALLOC (e, 1, mpq_qsformat_error);
    rval = mpq_ILLformat_error_create (e, error->type, error->desc,
	error->lineNumber,
	(mem->hasErrorLines) ? error->theLine : NULL,
	error->at);
    ILL_CLEANUP_IF (rval);
    e->next = mem->error_list;
    mem->error_list = e;
    mem->nerror++;
    mem->has_error[error->type]++;

CLEANUP:
    if (rval) {
	mpq_ILLformat_error_delete (e);
	ILL_IFFREE (e, mpq_qsformat_error);
    }
    return rval;
}
