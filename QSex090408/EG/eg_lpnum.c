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
#include "eg_lpnum.h"
/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** @brief type-dependant constants and helper numbers @{ */
mpz_t __zeroLpNum_mpz__;
mpz_t __oneLpNum_mpz__;
mpz_t __lpnum_mpz__;
mpq_t __zeroLpNum_mpq__;
mpq_t __oneLpNum_mpq__;
mpq_t __lpnum_mpq__;
mpf_t __zeroLpNum_mpf__;
mpf_t __oneLpNum_mpf__;
mpf_t __lpnum_mpf__;
mpf_t mpf_eps;
unsigned long int EGLPNUM_PRECISION = 128;
/** @} */

/* ========================================================================= */
/** @brief internal marker where we store wether or not we have already
 * setup-unsetup the static EGlpNumStart part */
int __EGlpNum_setup=0;

/* ========================================================================= */
/** @brief This function must be called once to initialize internal data needed
 * for all supported EGlpNum_t type, note that this is achieved through the
 * constructor attribute available in GCC.
 * @par Description:
 * This function initiaize all internal data needed to operate with the 
 * EGlpNum_t types, this function must be called before any operation on
 * EGlpNum_t types, and must be called only once.
 * */
void EGlpNumStart(void) __attribute__ ((constructor));
void EGlpNumStart(void)
{
	if(__EGlpNum_setup) return;
	mpf_set_default_prec (EGLPNUM_PRECISION);
	mpz_init (__zeroLpNum_mpz__);
	mpz_init (__oneLpNum_mpz__);
	mpz_set_ui (__zeroLpNum_mpz__, (unsigned long int)0);
	mpz_set_ui (__oneLpNum_mpz__, (unsigned long int)1);
	mpf_init (__zeroLpNum_mpf__);
	mpf_init (__oneLpNum_mpf__);
	mpf_set_ui (__oneLpNum_mpf__, (unsigned long int)1);
	mpf_set_ui (__zeroLpNum_mpf__, (unsigned long int)0);
	mpf_init (__lpnum_mpf__);
	mpf_init_set_ui (mpf_eps, (unsigned long int)1);
	mpf_div_2exp (mpf_eps, mpf_eps,  (unsigned long int)(EGLPNUM_PRECISION - 1));
	mpq_init (__zeroLpNum_mpq__);
	mpq_init (__oneLpNum_mpq__);
	mpq_init (__lpnum_mpq__);
	mpq_set_ui (__oneLpNum_mpq__, (unsigned long int)1, (unsigned long int)1);
	mpq_set_ui (__zeroLpNum_mpq__, (unsigned long int)0, (unsigned long int)1);
	__EGlpNum_setup=1;
}

/* ========================================================================= */
void EGlpNumSetPrecision (const unsigned prec)
{
	EGLPNUM_PRECISION = prec;
	mpf_set_default_prec (EGLPNUM_PRECISION);
	mpf_clear (__lpnum_mpf__);
	mpf_clear (mpf_eps);
	mpf_init (__lpnum_mpf__);
	mpf_init_set_ui (mpf_eps, (unsigned long int)1);
	mpf_div_2exp (mpf_eps, mpf_eps,  (unsigned long int)(EGLPNUM_PRECISION - 1));
}

/* ========================================================================= */
/** @brief This function must be called at the end of the program to free all
 * internal data used in the EGlpNum_t structures, once this function is called
 * any operation on EGlpNum_t types may fail. This behavior is now handled
 * through the destructor attribute available in GCC.
 * */
void EGlpNumExit(void) __attribute__ ((destructor));
void EGlpNumExit(void)
{
	if(!__EGlpNum_setup) return;
	mpf_clear (__zeroLpNum_mpf__);
	mpf_clear (__oneLpNum_mpf__);
	mpf_clear (__lpnum_mpf__);
	mpf_clear (mpf_eps);
	mpq_clear (__zeroLpNum_mpq__);
	mpq_clear (__oneLpNum_mpq__);
	mpq_clear (__lpnum_mpq__);
	mpz_clear (__zeroLpNum_mpz__);
	mpz_clear (__oneLpNum_mpz__);
	__EGlpNum_setup=0;
}

/* ========================================================================= */
void mpq_EGlpNumSet_mpf (mpq_t var,
												 mpf_t flt)
{
	/* local variables */
	unsigned long int __lsgn = mpf_cmp_ui (flt, (unsigned long int)0) < 0 ? (unsigned long int)1 : (unsigned long int)0;
	mpz_t __utmp,
	  __z[7],
	  max_den;
	long int __lexp = 0;
	int i;
	unsigned long int uexp,
	  cexp;
	mpf_t __cvl;
	/* check if the given number is zero, if so, set to zero var and return */
	if (mpf_cmp_ui (flt, (unsigned long int)0) == 0)
	{
		mpq_set_ui (var, (unsigned long int)0, (unsigned long int)1);
		return;
	}
	/* if not, then we have some work to do */
	/* now we initialize the internal numbers */
	mpf_init (__cvl);
	mpf_abs (__cvl, flt);
	mpz_init_set_ui (__utmp, (unsigned long int)0);
	for (i = 7; i--;)
		mpz_init_set_ui (__z[i], (unsigned long int)0);
	mpz_set_ui (__z[0], (unsigned long int)1);
	mpz_set_ui (__z[4], (unsigned long int)1);
	/* max_den is the maximum denominator that we want to see, this number should
	 * be sligtly larger than the square root of 2^EGLPNUM_PRECISION */
	mpz_init_set_ui (max_den, (unsigned long int)1);
	mpz_mul_2exp (max_den, max_den, EGLPNUM_PRECISION >> 1);
	/* first we compute the exponent stored in the limbs */
	__lexp = __cvl->_mp_exp * __GMP_BITS_PER_MP_LIMB;
	if (__lexp < 0)
	{
		uexp = -__lexp;
		mpf_mul_2exp (__cvl, __cvl, (unsigned long int) uexp);
	}
	else
	{
		uexp = __lexp;
		mpf_div_2exp (__cvl, __cvl, (unsigned long int) uexp);
	}
	/* now we compute the 2^n part needed to set this number between 0.5 and 1 */
	cexp = (unsigned long int)1 << 6;
	while (mpf_cmp_ui (__cvl, (unsigned long int)1) > 0 && cexp)
	{
		mpf_mul_2exp (__lpnum_mpf__, __oneLpNum_mpf__, (unsigned long int) cexp);
		if (mpf_cmp (__cvl, __lpnum_mpf__) > 0)
		{
			__lexp += cexp;
			mpf_div_2exp (__cvl, __cvl, (unsigned long int) cexp);
		}
		cexp = cexp >> 1;
	}
	cexp = (unsigned long int)1 << 6;
	while (mpf_cmp_d (__cvl, 0.5) < 0 && cexp)
	{
		mpf_div_2exp (__lpnum_mpf__, __oneLpNum_mpf__, (unsigned long int) cexp);
		if (mpf_cmp (__cvl, __lpnum_mpf__) < 0)
		{
			__lexp -= cexp;
			mpf_mul_2exp (__cvl, __cvl, (unsigned long int) cexp);
		}
		cexp = cexp >> 1;
	}
	/* now we loop until the next t's is more than mpf_eps */
	/* the formula is 
	 * p_i = t_i*p_{i-1} + p_{i-2}, and 
	 * q_i = t_i*q_{i-1} + q_{i-2} 
	 * note that |x-p_i/q_i|<1/q_i^2
	 * for us t_i = __utmp, and the current number is either [0,1,2] in the __z
	 * array, we use those popsitions ciclicly, and use the four position as a
	 * temporary number, __z+4 is used to store q's, at the beginning i = 1. */
	while (1)
	{
		if (mpf_cmp (__cvl, mpf_eps) < 0 || (mpz_cmp (__z[4], max_den) > 0))
		{
			mpz_set (mpq_denref (var), __z[4]);
			mpz_set (mpq_numref (var), __z[1]);
			break;
		}
		/* first run */
		mpf_ui_div (__cvl, (unsigned long int)1, __cvl);
		mpz_set_f (__utmp, __cvl);
		mpf_set_z (__lpnum_mpf__, __utmp);
		mpf_sub (__cvl, __cvl, __lpnum_mpf__);
		mpz_set (__z[6], __utmp);
		mpz_set (__z[2], __z[0]);
		mpz_addmul (__z[2], __z[1], __z[6]);
		mpz_set (__z[5], __z[3]);
		mpz_addmul (__z[5], __z[4], __z[6]);
		if (mpf_cmp (__cvl, mpf_eps) < 0 || (mpz_cmp (__z[5], max_den) > 0))
		{
			mpz_set (mpq_denref (var), __z[5]);
			mpz_set (mpq_numref (var), __z[2]);
			break;
		}
		/* second run */
		mpf_ui_div (__cvl, (unsigned long int)1, __cvl);
		mpz_set_f (__utmp, __cvl);
		mpf_set_z (__lpnum_mpf__, __utmp);
		mpf_sub (__cvl, __cvl, __lpnum_mpf__);
		mpz_set (__z[6], __utmp);
		mpz_set (__z[0], __z[1]);
		mpz_addmul (__z[0], __z[2], __z[6]);
		mpz_set (__z[3], __z[4]);
		mpz_addmul (__z[3], __z[5], __z[6]);
		if (mpf_cmp (__cvl, mpf_eps) < 0 || (mpz_cmp (__z[3], max_den) > 0))
		{
			mpz_set (mpq_denref (var), __z[3]);
			mpz_set (mpq_numref (var), __z[0]);
			break;
		}
		/* third run */
		mpf_ui_div (__cvl, (unsigned long int)1, __cvl);
		mpz_set_f (__utmp, __cvl);
		mpf_set_z (__lpnum_mpf__, __utmp);
		mpf_sub (__cvl, __cvl, __lpnum_mpf__);
		mpz_set (__z[6], __utmp);
		mpz_set (__z[1], __z[2]);
		mpz_addmul (__z[1], __z[0], __z[6]);
		mpz_set (__z[4], __z[5]);
		mpz_addmul (__z[4], __z[3], __z[6]);
	}

	/* ending */
	mpq_canonicalize (var);
	if (__lsgn)
		mpq_neg (var, var);
	if (__lexp > 0)
		mpq_mul_2exp (var, var, (unsigned long int) __lexp);
	if (__lexp < 0)
		mpq_div_2exp (var, var, (unsigned long int) (-__lexp));
	for (i = 7; i--;)
		mpz_clear (__z[i]);
	mpf_clear (__cvl);
	mpz_clear (max_den);
	mpz_clear (__utmp);
	return;
}

/* ========================================================================= */
void mpq_EGlpNumSet (mpq_t var,
										 double const dbl)
{
	/* local variables */
	double __dbl = dbl;
	unsigned __lsgn = __dbl > 0.0 ? (unsigned long int)0 : (unsigned long int)1;
	unsigned long __utmp = 0;
	int __lexp = 0;
	double __cvl = __dbl = fabs (__dbl);
	/* we use the first three numbers for p, and the last three numbers for q */
	/* first check that the dbl is not zero */
	if (__dbl < 1e-151)
	{
		mpq_set_ui (var, (unsigned long int)0, (unsigned long int)1);
		__lsgn = 0;
	}
	else if (__dbl > 1e151)
		mpq_set_d (var, __dbl);
	else
	{
		/* now we initialize the integer numbers */
		mpz_t __z[7];
		for (__utmp = 7; __utmp--;)
			mpz_init (__z[__utmp]);
		mpz_set_ui (__z[0], (unsigned long int)1);
		mpz_set_ui (__z[4], (unsigned long int)1);
		/* now we compute the 2^n part needed to set this number between 0 and 1 */
#define __HI_EXP(x,e,v,lv) {if( x >=v ){ e = e + lv; x /= v;}}
		if (__cvl > 1)
		{
			__HI_EXP (__cvl, __lexp,
								115792089237316195423570985008687907853269984665640564039457584007913129639936.0,
								256);
			__HI_EXP (__cvl, __lexp, 340282366920938463463374607431768211456.0, 128);
			__HI_EXP (__cvl, __lexp, 18446744073709551616.0, 64);
			__HI_EXP (__cvl, __lexp, 4294967296.0, 32);
			__HI_EXP (__cvl, __lexp, 65536.0, 16);
			__HI_EXP (__cvl, __lexp, 256.0, 8);
			__HI_EXP (__cvl, __lexp, 16.0, 4);
			__HI_EXP (__cvl, __lexp, 4.0, 2);
			__HI_EXP (__cvl, __lexp, 2.0, 1);
#undef __HI_EXP
		}
		else if (__cvl < 0.5)
		{
#define __LO_EXP(x,e,v,lv) {if( x < 1/v ) { e = e - lv; x *= v;}}
			__LO_EXP (__cvl, __lexp,
								115792089237316195423570985008687907853269984665640564039457584007913129639936.0,
								256);
			__LO_EXP (__cvl, __lexp, 340282366920938463463374607431768211456.0, 128);
			__LO_EXP (__cvl, __lexp, 18446744073709551616.0, 64);
			__LO_EXP (__cvl, __lexp, 4294967296.0, 32);
			__LO_EXP (__cvl, __lexp, 65536.0, 16);
			__LO_EXP (__cvl, __lexp, 256.0, 8);
			__LO_EXP (__cvl, __lexp, 16.0, 4);
			__LO_EXP (__cvl, __lexp, 4.0, 2);
			__LO_EXP (__cvl, __lexp, 2.0, 1);
#undef __LO_EXP
		}
		/* now we loop until the next t's is more than EGLPNUM_MINEPS */
		/* the formula is 
		 * p_i = t_i*p_{i-1} + p_{i-2}, and 
		 * q_i = t_i*q_{i-1} + q_{i-2} 
		 * note that |x-p_i/q_i|<1/q_i^2
		 * for us t_i = __utmp, and the current number is either [0,1,2] in the __z
		 * array, we use those popsitions ciclicly, and use the four position as a
		 * temporary number, __z+4 is used to store q's, at the beginning i = 1. */
		while (1)
		{
			if (__cvl < EGLPNUM_MINEPS || (mpz_cmp_ui (__z[4], (unsigned long int)0xfffffff) > 0))
			{
				mpz_set (mpq_denref (var), __z[4]);
				mpz_set (mpq_numref (var), __z[1]);
				break;
			}
			MESSAGE(MPQ_VERBOSE_CNT_FRAC,"cur approximate %ld/%ld, error %10.7lg", 
							mpz_get_si(__z[1]), mpz_get_si(__z[4]), __cvl);
			/* first run */
			__cvl = 1 / __cvl;
			__utmp = floor (__cvl);
			__cvl -= __utmp;
			mpz_set_ui (__z[6], __utmp);
			mpz_set (__z[2], __z[0]);
			mpz_addmul (__z[2], __z[1], __z[6]);
			mpz_set (__z[5], __z[3]);
			mpz_addmul (__z[5], __z[4], __z[6]);
			if (__cvl < EGLPNUM_MINEPS || (mpz_cmp_ui (__z[5], (unsigned long int)0xfffffff) > 0))
			{
				mpz_set (mpq_denref (var), __z[5]);
				mpz_set (mpq_numref (var), __z[2]);
				break;
			}
			MESSAGE(MPQ_VERBOSE_CNT_FRAC,"cur approximate %ld/%ld, error %10.7lg", 
							mpz_get_si(__z[2]), mpz_get_si(__z[5]), __cvl);
			/* second run */
			__cvl = 1 / __cvl;
			__utmp = floor (__cvl);
			__cvl -= __utmp;
			mpz_set_ui (__z[6], __utmp);
			mpz_set (__z[0], __z[1]);
			mpz_addmul (__z[0], __z[2], __z[6]);
			mpz_set (__z[3], __z[4]);
			mpz_addmul (__z[3], __z[5], __z[6]);
			if (__cvl < EGLPNUM_MINEPS || (mpz_cmp_ui (__z[3], (unsigned long int)0xfffffff) > 0))
			{
				mpz_set (mpq_denref (var), __z[3]);
				mpz_set (mpq_numref (var), __z[0]);
				break;
			}
			MESSAGE(MPQ_VERBOSE_CNT_FRAC,"cur approximate %ld/%ld, error %10.7lg", 
							mpz_get_si(__z[0]), mpz_get_si(__z[3]), __cvl);
			/* third run */
			__cvl = 1 / __cvl;
			__utmp = floor (__cvl);
			__cvl -= __utmp;
			mpz_set_ui (__z[6], __utmp);
			mpz_set (__z[1], __z[2]);
			mpz_addmul (__z[1], __z[0], __z[6]);
			mpz_set (__z[4], __z[5]);
			mpz_addmul (__z[4], __z[3], __z[6]);
		}
		for (__utmp = 7; __utmp--;)
			mpz_clear (__z[__utmp]);
	}
	/* ending */
	mpq_canonicalize (var);
	if (__lsgn)
		mpq_neg (var, var);
	if (__lexp > 0)
		mpq_mul_2exp (var, var, (unsigned long int) __lexp);
	if (__lexp < 0)
		mpq_div_2exp (var, var, (unsigned long int) (-__lexp));
	return;
}

/* ========================================================================= */
int mpz_EGlpNumReadStr (mpz_t var,
												char const *str)
{
	/* local variables */
	char unsigned a_sgn = 1;
	char unsigned sgn = 0;
	char c = 0;
	unsigned n_char = 0;
	/* now we read the string */
	c = str[n_char];
	mpz_set_ui (var, (unsigned long int)0);
	while ((('0' <= c) && (c <= '9')) ||	/* allow to read digits */
				 (a_sgn && (c == '+' || c == '-')) /* allow sign for exponent */ )
	{
		switch (c)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			mpz_mul_ui (var, var, (unsigned long int)10);
			mpz_add_ui (var, var,  (unsigned long int)(c - '0'));
			a_sgn = 0;
			break;
		case '-':
			sgn = 1;
		case '+':
			a_sgn = 0;
			break;
		}
		/* advance the reading character */
		c = str[++n_char];
	}
	if (sgn)
		mpz_neg (var, var);
	return n_char;
}

/* ========================================================================= */
int mpq_EGlpNumReadStrXc (mpq_t var,
													char const *str)
{
	/* local variables */
	char unsigned a_dot = 1,
	  a_exp = 0,
	  a_exp_sgn = 0,
	  a_sgn = 1,
	  a_div = 1;
	char c = 0;
	int l_exp = 0,
	  sgn = 0,
	  exp_sgn = 0;
	unsigned n_char = 0,
	  n_dig = 0,
	  cn = 0;
	mpq_t den[2];
	mpq_init (den[0]);
	mpq_init (den[1]);
	mpq_set_ui (den[1], (unsigned long int)1, (unsigned long int)1);
	mpq_set_ui (den[0], (unsigned long int)0, (unsigned long int)1);

	/* now we read the string */
	c = str[n_char];
	while ((('0' <= c) && (c <= '9')) ||	/* allow to read digits */
				 (a_dot && (c == '.')) ||	/* allow to read a dot point */
				 (a_exp && (c == 'e' || c == 'E')) ||	/* allow an exponent marker */
				 (a_sgn && (c == '+' || c == '-')) ||	/* allow a number sign */
				 (a_div && c == '/') ||	/* allow the division sign */
				 (a_exp_sgn && (c == '+' || c == '-')) /* allow sign for exponent */ )
	{
		switch (c)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			/* if we haven't read the exponent then the digits bellongs to the mantisa
			 * */
			if (a_exp || n_dig == 0)
			{
				if (!a_dot)
					mpz_mul_ui (mpq_denref (den[cn]), mpq_denref (den[cn]), (unsigned long int)10);
				mpz_mul_ui (mpq_numref (den[cn]), mpq_numref (den[cn]), (unsigned long int)10);
				mpz_add_ui (mpq_numref (den[cn]), mpq_numref (den[cn]),
										 (unsigned long int)(c - '0'));
				n_dig++;
				a_exp = 1;
			}
			/* otherwise, if we have read the exponent, the digits should go to the
			 * exponent */
			else
			{
				l_exp = 10 * l_exp + c - '0';
				a_exp_sgn = 0;
			}
			a_sgn = 0;
			break;
		case '.':
			a_sgn = 0;
			a_dot = 0;
			a_sgn = 0;
			break;
		case '-':
			if (a_sgn)
				sgn = 1;
			else
				exp_sgn = 1;
		case '+':
			if (a_sgn)
				a_sgn = 0;
			if (a_exp_sgn)
				a_exp_sgn = 0;
			break;
		case 'e':
		case 'E':
			a_sgn = 0;
			a_exp = 0;
			a_exp_sgn = 1;
			break;
		case '/':
			if (exp_sgn)
				l_exp = -l_exp;
			if (l_exp > 0)
				while (l_exp--)
					mpz_mul_ui (mpq_numref (den[0]), mpq_numref (den[0]), (unsigned long int)10);
			else if (l_exp < 0)
			{
				l_exp = -l_exp;
				while (l_exp--)
					mpz_mul_ui (mpq_denref (den[0]), mpq_denref (den[0]), (unsigned long int)10);
			}
			if (sgn)
				mpz_neg (mpq_numref (den[0]), mpq_numref (den[0]));
			mpq_canonicalize (den[0]);
			mpq_set_ui (den[1], (unsigned long int)0, (unsigned long int)1);
			sgn = 0;
			exp_sgn = 0;
			l_exp = 0;
			a_div = 0;
			n_dig = 0;
			a_dot = 1;
			a_exp = 0;
			a_exp_sgn = 0;
			a_sgn = 1;
			cn = 1;
			break;
		}
		/* advance the reading character */
		c = str[++n_char];
	}
	if (n_char)
	{
		/* now expand the exponent of the denominator */
		if (exp_sgn)
			l_exp = -l_exp;
		if (l_exp > 0)
			while (l_exp--)
				mpz_mul_ui (mpq_numref (den[cn]), mpq_numref (den[cn]), (unsigned long int)10);
		else if (l_exp < 0)
		{
			l_exp = -l_exp;
			while (l_exp--)
				mpz_mul_ui (mpq_denref (den[cn]), mpq_denref (den[cn]), (unsigned long int)10);
		}
		/* check the sign of the whole number */
		if (sgn)
			mpz_neg (mpq_numref (den[cn]), mpq_numref (den[cn]));
		/* ending */
		mpq_canonicalize (den[0]);
		mpq_canonicalize (den[1]);
		mpq_div (var, den[0], den[1]);
	}
	mpq_clear (den[0]);
	mpq_clear (den[1]);
	return n_char;
}

/* ========================================================================= */
/** @} */
