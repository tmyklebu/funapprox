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
#ifndef __EG_LPNUM_MPZ__
#define __EG_LPNUM_MPZ__
#include "eg_lpnum.h"
/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** extern definitions of constaants for different set-ups */
extern mpz_t __zeroLpNum_mpz__;
extern mpz_t __oneLpNum_mpz__;
extern mpz_t __lpnum_mpz__;
#define mpz_zeroLpNum __zeroLpNum_mpz__
#define mpz_oneLpNum  __oneLpNum_mpz__
#define mpz_epsLpNum  __zeroLpNum_mpz__

/* ========================================================================= */
/** @brief Read from a string a number and store it in the given mpz_t, return
 * the number of chars readed from the input string */
int mpz_EGlpNumReadStr (mpz_t a,
												const char *const str);

/* ========================================================================= */
/** @brief given a mpz_t, write it to a string (to be allocated internally), 
 * and return it. */
#define mpz_EGlpNumGetStr(a) ({\
	char *__str=0;\
	__str = mpz_get_str(__str,10,a);\
	__str;})

/* ========================================================================= */
/** @brief given an array of type mpz_t, free it, if the pointer is NULL
 * nothing happen. */
#define mpz_EGlpNumFreeArray(ea) ({\
	size_t __sz = __EGlpNumArraySize(ea);\
	mpz_t* __ptr__ = (ea);\
	while(__sz--) mpz_clear(__ptr__[__sz]);\
	__EGlpNumFreeArray(ea);})


/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * mpz_t and return it, if no more memory, exit(1) */
#define mpz_EGlpNumReallocArray(lptr, lsize) ({\
	mpz_t **__ptr__ = (lptr);\
	size_t *__ntmp__ = (size_t *) *__ptr__, __sz__ = (lsize);\
	size_t __psz__;\
	/* if no memory allocated before we just call the regular allocator */\
	if (!*__ptr__)\
		*__ptr__ = mpz_EGlpNumAllocArray (__sz__);\
	else\
	{\
		/* first check that the previous size is not larger than the current */\
		__ntmp__--;\
		__psz__ = __ntmp__[0];\
		if (__psz__ < __sz__)\
		{\
			/* now we have to do the reallocation */\
			*__ptr__ = (mpz_t *) __ntmp__;\
			*__ptr__ = EGrealloc(*__ptr__,sizeof(mpz_t) * __sz__ + sizeof(size_t));\
			__ntmp__ = (size_t *) *__ptr__;\
			__ntmp__[0] = __sz__;\
			__ntmp__++;\
			*__ptr__ = (mpz_t *) __ntmp__;\
			for (; __psz__ < __sz__; __psz__++) mpz_init ((*__ptr__)[__psz__]);\
		}\
	}\
})

/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type mpz_t
 * and return it, if no more memory, exit(1) */
#define mpz_EGlpNumAllocArray(size) ({\
	size_t __i__ = (size);\
	mpz_t *__res = __EGlpNumAllocArray(mpz_t,__i__);\
	while(__i__--) mpz_init(__res[__i__]);\
	__res;})

/* ========================================================================= */
/** @brief set the given number pointer, set its value to the given double.
 * @param var mpz_t where we will store the double value.
 * @param dbl double value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is a number and not a pointer to that value, be carefull with this
 * detail. Also, due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
#define mpz_EGlpNumSet (var, dbl) mpz_set_d(var,dbl)

/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(a,b) <==> a= ceil(b) */
#define mpz_EGlpNumCeil(a, b) ({\
	mpz_cdiv_qr (mpz_numref (a), mpz_denref (a), mpz_numref (b), mpz_denref (b));\
	mpz_set_ui (mpz_denref (a), 1);\
})

/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(a,b) <==> a= floor(b) */
#define mpz_EGlpNumFloor(a, b) ({\
	mpz_fdiv_q (mpz_numref (a), mpz_numref (b), mpz_denref (b));\
	mpz_set_ui (mpz_denref (a), 1);\
})

/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define mpz_EGlpNumInv(a) mpz_inv(a,a)

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a mpz_t first number to compare.
 * @param b mpz_t second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define mpz_EGlpNumIsEqqual(a,b) (mpz_equal(a,b))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a mpz_t first number to compare.
 * @param b mpz_t second number to compare.
 * @param error mpz_t maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance 'error',
 * return 1 if |a-b|<= error, otherwise it return 0.
 * */
#define mpz_EGlpNumIsEqual(a,b,error) (mpz_equal(a,b))
#define mpz_EGlpNumIsNeq(a,b,error)   (!(mpz_equal(a,b)))
#define mpz_EGlpNumIsNeqq(a,b)        (!(mpz_equal(a,b)))
#define mpz_EGlpNumIsNeqZero(a,error) (!(mpz_equal(a,mpz_zeroLpNum)))
#define mpz_EGlpNumIsNeqqZero(a)      (!(mpz_equal(a,mpz_zeroLpNum)))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpz_t the first number.
 * @param b mpz_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define mpz_EGlpNumIsLess(a,b) (mpz_cmp(a,b) < 0)

/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a mpz_t the first number.
 * @param b mpz_t the second number
 * @param c mpz_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define mpz_EGlpNumIsSumLess(a, b, c) ({\
	mpz_add (__lpnum_mpz__, a, b);\
	(mpz_cmp (__lpnum_mpz__, c) < 0);\
})

/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a mpz_t the first number.
 * @param b mpz_t the second number
 * @param c mpz_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define mpz_EGlpNumIsDiffLess(a, b, c) ({\
	mpz_sub (__lpnum_mpz__, a, b);\
	(mpz_cmp (__lpnum_mpz__, c) < 0);\
})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpz_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define mpz_EGlpNumIsLessDbl(a,b) (mpz_get_d(a) < b)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpz_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define mpz_EGlpNumIsGreaDbl(a,b) (mpz_get_d(a) > b)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpz_t the first number.
 * @param b mpz_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define mpz_EGlpNumIsLeq(a,b) (mpz_cmp(a,b) <= 0)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a mpz_t source number (it won't change value).
 * @param b mpz_t source number (it won't change value).
 * @param c mpz_t denominator of the difference (it won't change value).
 * @param d mpz_t where to store the value .
 * @par Description:
 * Set @f$a = \frac{b - c}{d} @f$ */
#define mpz_EGlpNumCopyDiffRatio(a, b, c, d) ({\
	mpz_sub (__lpnum_mpz__, b, c);\
	mpz_div (a, __lpnum_mpz__, d);\
})

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a mpz_t source number (it won't change value).
 * @param b mpz_t source number (it won't change value).
 * @param dest mpz_t where to store the value stored in 'orig'.
 * @par Description:
 * Set dest = a - b */
#define mpz_EGlpNumCopyDiff(dest,a,b) mpz_sub(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a mpz_t source number (it won't change value).
 * @param b mpz_t source number (it won't change value).
 * @param dest mpz_t where to store the sum.
 * @par Description:
 * Set dest = a + b */
#define mpz_EGlpNumCopySum(dest,a,b) mpz_add(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param orig mpz_t source number (it won't change value).
 * @param dest mpz_t where to store the value stored in 'orig'.
 * @par Description:
 * Given two numbers copy the values in 'orig', into 'dest'.
 * */
#define mpz_EGlpNumCopy(dest,orig) mpz_set(dest,orig)

/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param orig mpz_t source number (it won't change value).
 * @param dest mpz_t where to store the value stored in 'orig'.
 * @par Description:
 * implement dest = max(dest,abs(orig))
 * */
#define mpz_EGlpNumSetToMaxAbs(dest, orig) ({\
	mpz_abs (__lpnum_mpz__, orig);\
	if (mpz_cmp (dest, __lpnum_mpz__) < 0)\
		mpz_set (dest, __lpnum_mpz__);\
})

#define mpz_EGlpNumSetToMinAbs(dest, orig) ({\
	mpz_abs (__lpnum_mpz__, orig);\
	if (mpz_cmp (dest, __lpnum_mpz__) > 0)\
		mpz_set (dest, __lpnum_mpz__);\
})

/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param dest mpz_t where to store the result
 * @param orig mpz_t second parameter
 * @param den mpz_t third parameter
 * @par Description:
 * compute dest = (orig*orig)/den
 * */
#define mpz_EGlpNumCopySqrOver(dest, orig, den) ({\
	mpz_mul (__lpnum_mpz__, orig, orig);\
	mpz_div (dest, __lpnum_mpz__, den);\
})

/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param orig mpz_t source number (it won't change value).
 * @param dest mpz_t where to store the absolute value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy its absolute value to 'dest'. i.e.
 * dest = |orig|
 * */
#define mpz_EGlpNumCopyAbs(dest,orig) mpz_abs(dest,orig)

/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param orig mpz_t the source number (it won't change value).
 * @param dest mpz_t where to store minus the value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy minus the value to 'dest'. i.e.
 * dest = -orig
 * */
#define mpz_EGlpNumCopyNeg(dest,orig) mpz_neg(dest,orig)

/* ========================================================================= */
/** @brief Set des = op1/op2.
 * @param dest mpz_t where we will store the result.
 * @param op1 mpz_t numerator of the fraction (possibly non an integer)
 * @param op2 mpz_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des = op1/op2
 * */
#define mpz_EGlpNumCopyFrac(dest,op1,op2) mpz_div(dest,op1,op2)

/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param orig mpz_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param dest mpz_t* pointer to where to store the first 'size' values 
 * stored in 'orig'.
 * @param size unsigned int specifying how many values of 'orig' will be copied
 * onto 'dest'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define mpz_EGlpNumCopyArray(dest,orig,size) {\
	register unsigned int __i__ = size;\
	for(;__i__--;)\
	{\
		mpz_set(dest[__i__],orig[__i__]);\
	}\
}

/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a mpz_t the number that we are going to Sub to.
 * @param b mpz_t value to be multiplyed.
 * @param c mpz_t value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define mpz_EGlpNumSubInnProdTo(a, b, c) ({\
	mpz_mul (__lpnum_mpz__, b, c);\
	mpz_sub (a, a, __lpnum_mpz__);\
})

/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a mpz_t the number that we are going to add to.
 * @param b mpz_t value to be multiplyed.
 * @param c mpz_t value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define mpz_EGlpNumAddInnProdTo(a, b, c) ({\
	mpz_mul (__lpnum_mpz__, b, c);\
	mpz_add (a, a, __lpnum_mpz__);\
})

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a mpz_t the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumSubUiTo(a,b) (mpz_submul_ui(mpz_numref(a),mpz_denref(a),b),mpz_canonicalize(a))

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a mpz_t the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumAddUiTo(a,b) (mpz_addmul_ui(mpz_numref(a),mpz_denref(a),b),mpz_canonicalize(a))

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a mpz_t the number that we are going to add to.
 * @param b mpz_t value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumAddTo(a,b) mpz_add(a,a,b)

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a mpz_t the number that we are going to substract
 * from.
 * @param b mpz_t value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumSubTo(a,b) mpz_sub(a,a,b)

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a mpz_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b mpz_t value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumMultTo(a,b) mpz_mul(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a mpz_t the number that we are going to divide by
 * the second number and store the result.
 * @param b mpz_t value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumDivTo(a,b) mpz_div(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a mpz_t the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumDivUiTo(a,b) (mpz_mul_ui(mpz_denref(a),mpz_denref(a),b),mpz_canonicalize(a))

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a mpz_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpz_EGlpNumMultUiTo(a,b) (mpz_mul_ui(mpz_numref(a),mpz_numref(a),b),mpz_canonicalize(a))

/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a mpz_t the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define mpz_EGlpNumZero(a) mpz_set_ui(a,0U,1U)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a mpz_t value to be set to one.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 1;
 * */
#define mpz_EGlpNumOne(a) mpz_set_ui(a,1U,1U)

/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a mpz_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define mpz_EGlpNumSign(a) mpz_neg(a,a)

/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a mpz_t number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in a.
 * */
#define mpz_EGlpNumToLf(a) mpz_get_d(a)

/* ========================================================================= */
/** @brief initialize the internal memory of a given variable */
#define mpz_EGlpNumInitVar(a) mpz_init(a)

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define mpz_EGlpNumClearVar(a) mpz_clear(a)

/* ========================================================================= */
/** This are the definitions needed to use mpz_t numbers within our general
 * framework @{ */
#if EGLPNUM_TYPE==GNU_MP_Z
#ifndef EGLPNUM_MAXCONVORDER
#define EGLPNUM_MAXCONVORDER 2U
#endif
#define EGLPNUM_ISFP 0
#define epsLpNum 													mpz_epsLpNum
#define oneLpNum 													mpz_oneLpNum
#define zeroLpNum 												mpz_zeroLpNum
#define EGlpNum_t 												mpz_t
#define EGlpNumAddTo(a,b)									mpz_EGlpNumAddTo(a,b)
#define EGlpNumAddUiTo(a,b) 							mpz_EGlpNumAddUiTo(a,b)
#define EGlpNumAddInnProdTo(a, b, c) 			mpz_EGlpNumAddInnProdTo(a,b,c)
#define EGlpNumAllocArray(size) 					mpz_EGlpNumAllocArray(size)
#define EGlpNumCeil(a, b) 								mpz_EGlpNumCeil(a,b)
#define EGlpNumClearVar(a) 								mpz_EGlpNumClearVar(a)
#define EGlpNumCopy(dest,orig) 						mpz_EGlpNumCopy(dest,orig)
#define EGlpNumCopyArray(dest,orig,size) 	mpz_EGlpNumCopyArray(dest,orig,size)
#define EGlpNumCopyDiff(dest,a,b) 				mpz_EGlpNumCopyDiff(dest,a,b)
#define EGlpNumCopyDiffRatio(a, b, c, d) 	mpz_EGlpNumCopyDiffRatio(a,b,c,d)
#define EGlpNumCopySum(dest,a,b) 					mpz_EGlpNumCopySum(dest,a,b)
#define EGlpNumCopySqrOver(dest,orig,den) mpz_EGlpNumCopySqrOver(dest,orig,den)
#define EGlpNumCopyAbs(dest,orig) 				mpz_EGlpNumCopyAbs(dest,orig)
#define EGlpNumCopyNeg(dest,orig) 				mpz_EGlpNumCopyNeg(dest,orig)
#define EGlpNumCopyFrac(dest,op1,op2) 		mpz_EGlpNumCopyFrac(dest,op1,op2)
#define EGlpNumDivTo(a,b) 								mpz_EGlpNumDivTo(a,b)
#define EGlpNumDivUiTo(a,b) 							mpz_EGlpNumDivUiTo(a,b)
#define EGlpNumFloor(a, b) 								mpz_EGlpNumFloor(a,b)
#define EGlpNumFreeArray(ea) 							mpz_EGlpNumFreeArray(ea)
#define EGlpNumGetStr(a) 									mpz_EGlpNumGetStr(a)
#define EGlpNumInitVar(a) 								mpz_EGlpNumInitVar(a)
#define EGlpNumInv(a) 										mpz_EGlpNumInv(a,b)
#define EGlpNumIsDiffLess(a, b, c) 				mpz_EGlpNumIsDiffLess(a,b,c)
#define EGlpNumIsEqqual(a,b) 							mpz_EGlpNumIsEqqual(a,b)
#define EGlpNumIsEqual(a,b,error) 				mpz_EGlpNumIsEqual(a,b,error)
#define EGlpNumIsGreaDbl(a,b) 						mpz_EGlpNumIsGreaDbl(a,b)
#define EGlpNumIsLessDbl(a,b) 						mpz_EGlpNumIsLessDbl(a,b)
#define EGlpNumIsLess(a,b) 								mpz_EGlpNumIsLess(a,b)
#define EGlpNumIsLeq(a,b) 								mpz_EGlpNumIsLeq(a,b)
#define EGlpNumIsNeq(a,b,error)   				mpz_EGlpNumIsNeq(a,b)
#define EGlpNumIsNeqq(a,b)        				mpz_EGlpNumIsNeqq(a,b)
#define EGlpNumIsNeqZero(a,error) 				mpz_EGlpNumIsNeqZero(a,b)
#define EGlpNumIsNeqqZero(a)      				mpz_EGlpNumIsNeqqZero(a)
#define EGlpNumIsSumLess(a, b, c) 				mpz_EGlpNumIsSumLess(a,b,c)
#define EGlpNumMultTo(a,b) 								mpz_EGlpNumMultTo(a,b)
#define EGlpNumMultUiTo(a,b) 							mpz_EGlpNumMultUiTo(a,b)
#define EGlpNumOne(a) 										mpz_EGlpNumOne(a)
#define EGlpNumReadStr(a,str) 						mpz_EGlpNumReadStrXc(a,str)
#define EGlpNumReallocArray(lptr,size) 		mpz_EGlpNumReallocArray(lptr,size)
#define EGlpNumSet(var,dbl) 							mpz_EGlpNumSet(var,dbl)
#define EGlpNumSetToMaxAbs(dest, orig) 		mpz_EGlpNumSetToMaxAbs(dest,orig)
#define EGlpNumSetToMinAbs(dest, orig) 		mpz_EGlpNumSetToMinAbs(dest,orig)
#define EGlpNumSign(a) 										mpz_EGlpNumSign(a)
#define EGlpNumSubTo(a,b) 								mpz_EGlpNumSubTo(a,b)
#define EGlpNumSubUiTo(a,b) 							mpz_EGlpNumSubUiTo(a,b)
#define EGlpNumSubInnProdTo(a, b, c) 			mpz_EGlpNumSubInnProdTo(a,b,c)
#define EGlpNumToLf(a) 										mpz_EGlpNumToLf(a)
#define EGlpNumZero(a) 										mpz_EGlpNumZero(a)
#endif

/** @} */
/* ========================================================================= */
/** @} */
#endif
