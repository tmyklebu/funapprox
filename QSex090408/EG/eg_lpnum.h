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
#ifndef __EG_LPNUM_H__
#define __EG_LPNUM_H__
/* ========================================================================= */
/** @defgroup EGlpNum EGlpNum
 *
 * Here we define a common interface to handle numbers in general, the idea is
 * to be able to work with infinite precicion numbers, plain doubles, floats,
 * integers, or fixed point numbers, without actually making different codes
 * for each of those types, rather we preffer to fix that at compyle time.
 *
 * @par History:
 * Revision 0.0.2
 *  We start doing the migration to gmp and 'natural' types, this means that we
 *  drop support for EGrat, this allow us to drop the requirement to use
 *  pointers, and instead we can just call the functions with the original
 *  parameters, still we have to be VERY carefull regarding changing
 *  local/external copies.
 *  - 2007-10-08
 *  					- Move EGswap, EGabs, EGmin and Egmax to eg_numutil.h
 *  - 2005-10-31
 *  					- Add EGswap to swap elements of any predefined type.
 *  - 2005-08-31
 *  					- Add EGmin and EGmax for built in types (i.e. for types where
 *  					the < comparison works as we want).
 *  - 2005-08-16
 *  					- Streamline mpq_EGlpNumGetStr
 *  					- Minor Fixes for zeroLpNum
 *  - 2005-07-29
 *  					- Add EGabs definition.
 *  - 2005-07-24
 *  					- Split eg_lpnum.h into different headers for each type of
 *  						suported numbers.
 *  					- Deprecate EGlpNumCOmpUFrac
 *  - 2005-05-26
 *  					- Add epsLpNum
 *  - 2005-05-17
 *  					- Add mpq_EGlpNumReadStrXc(mpq_t,str)
 *  - 2005-05-16
 *  					- Add mpq_EGlpNumSet_mpf(mpq,mpf)
 *  - 2005-05-12
 *  					- Add mpf_EGlpNumEpow(num,power)
 *  					- Add function to change precision of the numbers on the fly.
 *  					- Add EGlpNumReadStr to set a number from a given input string.
 *  					- Add EGlpNumGetStr to get (hopefully) the exact representation 
 *  						of the given input as string.
 *  - 2005-05-03
 *  					- Change the structure of the header so it provides an interface
 *  					in which a program can use all types of numbers at the same time,
 *  					this implies that we must define a start-up and clean-up function
 *  					that would initialize all constants for all numbers, and change
 *  					the naming scheme accordingly to support this.
 *  					- Deprecate EGlpNumInitArray
 *  					- Deprecate EGlpNumFreeIntArray
 *  					- Change all static-inline definitions to Statement Exprs style.
 *  - 2005-04-25
 *  					- Add EGlpNumIsSumLess(a,b,c)
 *  					- Add EGlpNumIsDiffLess(a,b,c)
 *  - 2005-04-13
 *  					- Add EGlpNumCopyDiffRatio(a,b,c,d)
 *  					- Add EGlpNumIsEqqual(a,b)
 *  					- Add EGlpNumIsEqual(a,b,error)
 *  					- Add EGlpNumCopy(a,b)
 *  					- Add EGlpNumIsLess(a,b)
 *  					- Add EGlpNumToLf(a)
 *  					- Add EGlpNumCopyDiff(a,b,c)
 *  					- Add EGlpNumCopyAbs(a,b)
 *  					- Add EGlpNumSubTo(a,b)
 *  					- Add EGlpNumAddTo(a,b)
 *  					- Add EGlpNumDivTo(a,b)
 *  					- Add EGlpNumMultTo(a,b)
 *  					- Add EGlpNumZero(a)
 *  					- Add EGlpNumOne(a)
 *  					- Add EGlpNumAddInnProdTo(a,b,c)
 *  					- Add EGlpNumSubInnProdTo(a,b,c)
 *  					- Add EGlpNumSign(a)
 *  					- Add EGlpNumCopyNeg(a,b)
 *  					- Add EGlpNumDivUiTo(a,b)
 *  					- Add EGlpNumMultUiTo(a,b)
 *  					- Add EGlpNumIsLeq(a,b)
 *  					- Add EGlpNumCopySqrOver(a,b,c)
 *  					- Add EGlpNumSet(a,b)
 *  					- Add EGlpNumCopyFrac(a,b,c)
 *  					- Add EGlpNumAddUiTo(a,b)
 *  					- Add EGlpNumSubUiTo(a,b)
 *  					- Add EGlpNumCopySum(a,b,c)
 *  					- Add EGlpNumInv(a)
 *  					- Add EGlpNumFloor(a)
 *  					- Add EGlpNumCeil(a)
 *  					- Add EGlpNumIsLessDbl(a,b)
 *  					- Add EGlpNumIsGreaDbl(a,b)
 *  					- Add EGlpNumSetToMaxAbs(a,b)
 *  					- Add EGlpNumAllocArray(size)
 *  					- Add EGlpNumFreeArray(a)
 *  					- Add EGlpNumReallocArray(a,size)
 *  					- Add EGlpNumInitVar(a)
 *  					- Add EGlpNumClearVar(a)
 *  					- Add EGlpNumIsNeqq(a,b)
 *  					- Add EGlpNumIsNeq(a,b,c)
 *  					- Add EGlpNumIsNeqqZero(a)
 *  					- Add EGlpNumIsNeqZero(a,b)
 *  					- Add EGlpNumSetToMinAbs(a,b)
 * Revision 0.0.1
 * - 2004-07-15
 * 						- Add support for GNU_MP_F types
 * - 2004-07-14
 * 						- Add support for GNU_MP_Q types
 * - 2004-07-12
 * 						- Add support for EG-rationals
 * - 2004-06-21
 * 						- First Implementation/Definition
 * */

/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <float.h>
/* ========================================================================= */
/** @name Number Types Definitions:
 * Define (as its name suggest) an internal identifier for the given 
 * type. this definitions are provided to select different types of data at 
 * compile  time, thus allowing us to provide limited template support. */
/*@{*/
/** C double type. */
#define DBL_TYPE 0
/** C float type. */
#define FLT_TYPE 1
/** C int type. */
#define INT_TYPE 2
/** EGlib #EGfp10_t type, this is an implementation of fixed precision
 * arithmetic with 10 bits for fractional representation. */
#define FP10_TYPE 3
/** EGlib #EGfp20_t type, this is an implementation of fixed precision
 * arithmetic with 20 bits for fractional representation. */
#define FP20_TYPE 4
/** EGlib #EGfp28_t type, this is an implementation of fixed precision
 * arithmetic with 28 bits for fractional representation. */
#define FP28_TYPE 5
/** EGlib #EGfp25_t type, this is an implementation of fixed precision
 * arithmetic with 25 bits for fractional representation. */
#define FP25_TYPE 6
/** GNU_MP library mpz_t type */
#define GNU_MP_Z 8
/** GNU_MP library mpq_t type */
#define GNU_MP_Q 9
/** GNU_MP library mpf_t type */
#define GNU_MP_F 10
/** C long double type */
#define LDBL_TYPE 11
/** C long long int type */
#define LLINT_TYPE 12
/** SoftFloat 128-bit floating point numbner */
#define FLOAT128_TYPE 13
/*@}*/


/* ========================================================================= */
#ifndef EGLPNUM_TYPE
/**
 * @brief default type for EGLPNUM_TYPE
 * @par Description:
 * If the type of number for EGlpNum is not defined beforehand (eg. eg_config.h
 * or through make.conf) we default to GNU_MP_Q.
 * available options:
 * 			- GNU_MP_Q
 * 			- GNU_MP_F
 * 			- DBL_TYPE
 * 			- LDBL_TYPE
 * */
#define EGLPNUM_TYPE DBL_TYPE
#endif
#include "eg_lpnum.dbl.h"
#include "eg_lpnum.int.h"

#include "gmp.h"
#include "eg_lpnum.mpz.h"
#include "eg_lpnum.mpq.h"
#include "eg_lpnum.mpf.h"

#include "eg_macros.h"
#include "eg_mem.h"
#include "eg_nummacros.h"
/* ========================================================================= */
/** @brief Debugging verbosity messages deped on the value of DEBUG (defined in
 * eg_configure.h) and on the value of EGLPNUM_DEBUGL macro defined here.
* */
#define EGLPNUM_DEBUGL 100

#ifndef EGLPNUM_MINEPS
/* ========================================================================= */
/** @brief This constant define the of the acuracy required while
 * converting doubles to rationals, a good number is 1e-5. More exactly, we 
 * stop the continued fraction method whenever the next e_i-[e_i] computed is
 * less than EGLPNUM_MINEPS. Note that this value can't be smaller than
 * 1/ULONG_MAX, otherwise we will have problems in the confertion step. */
#define EGLPNUM_MINEPS 0x1ep-20
#else
#if EGLPNUM_MINEPS < 3e-10
#undef EGLPNUM_MINEPS
#define EGLPNUM_MINEPS 3e-10
#endif
#endif

/* ========================================================================= */
/** @brief Set the default number of __BITS__ used in the precision of the
 * float point numbers (mpf_t), a normal double use up to 56-64 bits., the 
 * default precision is set to 128 */
extern unsigned long int EGLPNUM_PRECISION;

/* ========================================================================= */
/** @brief Change the default precision for mpf_t numbers. */
void EGlpNumSetPrecision (const unsigned prec);

/* ========================================================================= */
/** @brief Allocate an array of a given type and store (sizeof(size_t) bytes 
 * before the actual array) the size of the allocated array. 
 * @param type the type of the array to be returned.
 * @param size the length of the array to be returned, note that it can be
 * zero, in wich case no memory allocation is made and NULL is returned. */
#define __EGlpNumAllocArray(type,size) ({\
	size_t __sz = (size);\
	size_t *__utmp = __sz ? (size_t*) EGmalloc (sizeof(type) * __sz + sizeof(size_t)) : 0;\
	if(__sz) __utmp[0] = __sz;\
	(type*)(__sz ? (__utmp+1):0);})

/* ========================================================================= */
/** @brief Given an array allocated with __EGlpNumAllocArray, return the size of
 * the given array, if the array is null, return zero. 
 * @param array the array from where we need the size. */
/* ========================================================================= */
#define __EGlpNumArraySize(array) ({\
	size_t *__utmp = (size_t*)(array);\
	if(__utmp) __utmp--;\
	__utmp ? __utmp[0]:0;})

/* ========================================================================= */
/** @brief, given an array allocated by __EGlpNumAllocArray, free the allocated
 * memory.
 * @param array the array to be freed, it can be null. The given array will
 * always pooint to NULL when this function is done.
 * */
/* ========================================================================= */
#define __EGlpNumFreeArray(array) ({\
	size_t *__utmp = (size_t*)(array);\
	if(__utmp) free (__utmp-1);\
	(array) = 0;})

/* ========================================================================= */
/** @brief indicate if the global data has been initialized, if zero,
 * initialization routine should be called. This is provided to allow
 * syncronization of constructors between libraries */
extern int __EGlpNum_setup;
/* ========================================================================= */
/** @brief initialization routine for global data. This function is called as a
 * constructor, but calling it twice won't cause any problems, it is provided
 * to ensure that all EGlpnum globals are initialized at the beggining 
 * and in case they where not (__EGlpNMum_setup==0), then call the initializator */
extern void EGlpNumStart(void);
 /** @}*/
/* ========================================================================= */
#endif
