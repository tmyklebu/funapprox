/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */ 
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

/* RCSINFO $Id: dbl_editor.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef dbl_EDITOR_H
#define dbl_EDITOR_H

extern void dbl_ILLeditor_init (void);
extern void dbl_ILLeditor (dbl_QSdata * p);
extern int dbl_ILLeditor_solve (dbl_QSdata * p,
														int salgo);

#endif
