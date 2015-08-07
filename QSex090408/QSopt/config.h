/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza.  Sanjeeb Dash's ownership of copyright in          */
/*  QSopt_ex is derived from his copyright in QSopt.                        */
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

/* RCSINFO $Id: config.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef __CONFIG_H
#define __CONFIG_H

/* INCLUDE/config.h.  Generated automatically by configure.  */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/


/* Define if your compiler is missing the appropriate function prototype */

/* #undef ILL_PROTO_PRINTF */
/* #undef ILL_PROTO_GETHOSTNAME */
/* #undef ILL_PROTO_GETRUSAGE */

/* Define if you want to use posix threads */
/* #undef ILL_POSIXTHREADS */

/* Define to empty if the keyword `const' does not work.  */
/* #undef const */

/* Define to `int' if <sys/types.h> doesn't define.  */
/* #undef pid_t */

/* Define to `unsigned' if <sys/types.h> doesn't define.  */
/* #undef size_t */

/* Define to `unsigned char' if <sys/types.h> doesn't define.  */
/* #undef u_char */

/* Define to `int' if the builtin type `void' does not work.  */
/* #undef void */

/* Define if you have the gethostname function.  */

/* Define if you have the socket function.  */

/* Define if you have the strdup function.  */

/* Define if you have the getrusage function.  */
#ifndef WIN32
#ifndef HAVE_GETRUSAGE
#define HAVE_GETRUSAGE
#endif
#endif

/* Define if you have the times function.  */

/* Define if you have the clock function.  */
#ifndef HAVE_CLOCK
#define HAVE_CLOCK
#endif

/* Define if you have the sleep function.  */

/* Define if you have the <stdio.h> header file.  */
#define HAVE_STDIO_H 1

/* Define if you have the <stdarg.h> header file.  */
#define HAVE_STDARG_H 1

/* Define if you have the <stdlib.h> header file.  */
#define HAVE_STDLIB_H 1

/* Define if you have the <math.h> header file.  */
#define HAVE_MATH_H 1

/* Define if you have the <string.h> header file.  */
#define HAVE_STRING_H 1


/* Define if you have the <strings.h> header file.  */
#ifndef WIN32
#define HAVE_STRINGS_H 1
#endif

/* Define if you have the <errno.h> header file.  */
#define HAVE_ERRNO_H 1

/* Define if you have the <assert.h> header file.  */
#define HAVE_ASSERT_H 1

/* Define if you can safely include both <sys/time.h> and <time.h>.  */

/* Define if you have the <sys/time.h> header file.  */

/* Define if you have the <time.h> header file.  */
#ifndef HAVE_TIME_H
#define  HAVE_TIME_H
#endif

/* Define if you have the <stddef.h> header file.  */
#define HAVE_STDDEF_H 1

/* Define if you have the <unistd.h> header file.  */
#ifndef WIN32
#ifndef HAVE_UNISTD_H
#define HAVE_UNISTD_H
#endif
#endif

/* Define if you have the <malloc.h> header file.  */
#define HAVE_MALLOC_H 1

/* Define if you have the <sys/types.h> header file.  */

/* Define if you have the <sys/stat.h> header file.  */

/* Define if you have the <sys/resource.h> header file.  */
#ifndef WIN32
#ifndef HAVE_SYS_RESOURCE_H
#define HAVE_SYS_RESOURCE_H
#endif
#endif

/* Define if you have the <fcntl.h> header file.  */

/* Define if you have the <signal.h> header file.  */

/* Define if you have the <sys/socket.h> header file.  */

/* Define if you have the <netdb.h> header file.  */

/* Define if you have the <netinet/in.h> header file.  */

/* Define if you have the <sys/param.h> header file.  */

/* Define if you have the <sys/times.h> header file.  */

/* Define if your compiler supports attribute modifiers  */
/* such as __attribute__ ((unused)) (gcc 2.8.1 does)     */
#define ILL_ATTRIBUTE 1

/* Some machine (o/s) specific problems */

/* Define if unistd.h uses __vfork but does not prototype it */
/* This happens under Irix 6 */
/* #undef ILL_PROTO___VFORK */

#ifdef WIN32
#define HAVE_DOS_TIME
#endif

/* __CONFIG_H */
#endif
