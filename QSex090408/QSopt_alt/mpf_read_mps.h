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
#ifndef mpf_READ_MPS_STATE_H
#define mpf_READ_MPS_STATE_H

#include "mpf_iqsutil.h"

#include "mpf_mps.h"

typedef struct mpf_ILLread_mps_state_struct
{
	int section[ILL_MPS_N_SECTIONS];
	ILLmps_section active;
	const char *file_name;
	mpf_qsline_reader *file;
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
mpf_ILLread_mps_state;

extern int mpf_ILLmps_state_init (mpf_ILLread_mps_state * state,
															mpf_qsline_reader * file,
															const char *mpf_fname);
extern void mpf_ILLmps_state_clear (mpf_ILLread_mps_state * state);
extern int mpf_ILLmps_set_section (mpf_ILLread_mps_state * state,
															 const ILLmps_section sec);

extern int mpf_ILLmps_next_line (mpf_ILLread_mps_state * state);
extern int mpf_ILLmps_next_field (mpf_ILLread_mps_state * state);
extern int mpf_ILLmps_next_coef (mpf_ILLread_mps_state * state,
														 mpf_t * coef);
extern int mpf_ILLmps_next_bound (mpf_ILLread_mps_state * state,
															mpf_t * coef);
extern void mpf_ILLmps_check_end_of_line (mpf_ILLread_mps_state * state);
extern void mpf_ILLmps_set_end_of_line (mpf_ILLread_mps_state * state);

extern int mpf_ILLmps_int_sos_mode (mpf_ILLread_mps_state * state);

extern const char *mpf_ILLmps_possibly_blank_name (const char *field,
																							 mpf_ILLread_mps_state * state,
																							 ILLsymboltab * tab);
extern int mpf_ILLmps_empty_key (mpf_ILLread_mps_state * state);
extern int mpf_ILLmps_empty_field (mpf_ILLread_mps_state * state);

extern int mpf_ILLmps_error (mpf_ILLread_mps_state * state,
												 const char *format,
												 ...);
extern void mpf_ILLmps_warn (mpf_ILLread_mps_state * state,
												 const char *format,
												 ...);

#endif
