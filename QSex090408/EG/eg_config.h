/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/* ========================================================================= */
/** Main Configuration for the library, as debug levels and so on
 * 
 * @par History:
 * - 2006-01-27
 *					- Handle some problems with stdint.h in SUN
 * - 2005-08-17
 * 					- Set memory aligment to 8 bits
 * - 2003-06-02
 * 					- First Implementation
 * @version 0.0.1 
 * */
/* ========================================================================= */
#ifndef __EG_CONGIH_H__
#define __EG_CONFIG_H__
#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

/* ========================================================================= */
/** @brief we try to detect what type of OS we are working in, this is 
 * important because different flavor of unixes have small tweaks that 
 * we must take care of. 
 * So far, we can choose LINUX or SUN */
#define LINUX 0
#define SUN 1
#define OSXMAC 2
#define CYGWIN 3
#define AIX 4
#ifndef OS
#define OS LINUX
#endif
#if ( (OS > AIX) || (OS < LINUX) )
#error OS configuration value unknown
#endif

/* ========================================================================= */
/** @brief Depending on the OS, we include some files and make some 
 * definitions */
#if OS == LINUX
#include <getopt.h>
#include <stdint.h>
#endif


/* ========================================================================= */
/** @brief assert Debug options definitions, by defoult set on */
#ifndef DEBUG
#warning you should define DEBUG, assuming it to be 1
#define DEBUG 1
#endif

/* ========================================================================= */
/** @brief assert Verbose options definition, by default set on */
#ifndef VERBOSE_LEVEL
#warning you should define VERBOSE_LEVEL, assuming it to be 1
#define VERBOSE_LEVEL 1
#endif

/* end eg_config.h */
#endif
