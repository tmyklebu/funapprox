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

/* $RCSfile: rd_mps.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $ */
#ifndef dbl_READ_MPS_STATE_H
#define dbl_READ_MPS_STATE_H

#include "dbl_iqsutil.h"

#include "dbl_mps.h"

typedef struct dbl_ILLread_mps_state_struct
{
	int section[ILL_MPS_N_SECTIONS];
	ILLmps_section active;
	const char *file_name;
	dbl_qsline_reader *file;
	unsigned int line_num;
	unsigned int field_num;				/* number of successfully read fields on line */
	int intvar;
	int sosvar;
	char line[ILL_namebufsize];
	char key[ILL_namebufsize];
	char field[ILL_namebufsize];
	char *obj;
	char *p;											/* ptr to next 'unread' character */
}
dbl_ILLread_mps_state;

extern int dbl_ILLmps_state_init (dbl_ILLread_mps_state * state,
															dbl_qsline_reader * file,
															const char *dbl_fname);
extern void dbl_ILLmps_state_clear (dbl_ILLread_mps_state * state);
extern int dbl_ILLmps_set_section (dbl_ILLread_mps_state * state,
															 const ILLmps_section sec);

extern int dbl_ILLmps_next_line (dbl_ILLread_mps_state * state);
extern int dbl_ILLmps_next_field (dbl_ILLread_mps_state * state);
extern int dbl_ILLmps_next_coef (dbl_ILLread_mps_state * state,
														 double * coef);
extern int dbl_ILLmps_next_bound (dbl_ILLread_mps_state * state,
															double * coef);
extern void dbl_ILLmps_check_end_of_line (dbl_ILLread_mps_state * state);
extern void dbl_ILLmps_set_end_of_line (dbl_ILLread_mps_state * state);

extern int dbl_ILLmps_int_sos_mode (dbl_ILLread_mps_state * state);

extern const char *dbl_ILLmps_possibly_blank_name (const char *field,
																							 dbl_ILLread_mps_state * state,
																							 ILLsymboltab * tab);
extern int dbl_ILLmps_empty_key (dbl_ILLread_mps_state * state);
extern int dbl_ILLmps_empty_field (dbl_ILLread_mps_state * state);

extern int dbl_ILLmps_error (dbl_ILLread_mps_state * state,
												 const char *format,
												 ...);
extern void dbl_ILLmps_warn (dbl_ILLread_mps_state * state,
												 const char *format,
												 ...);

#endif
