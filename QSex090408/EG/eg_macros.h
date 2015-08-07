/* ========================================================================= */
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
/** @defgroup EGmacros General Macros
 * global macros and types for EGlib
 *
 * @version 0.9.2 
 * @par History:
 * - 2007-12-07
 * 						- Add FTESTG and FTEST, that always test the parameters
 * 							regardless of the debug level
 * - 2007-01-19
 * 						- Delete EGosGetOffset
 * - 2006-09-28
 * 						- Add function that display basic process information, including
 * 						version and date of compilation of EGlib
 * - 2005-12-19
 * 						- Add float128 support ussing SoftFloat.
 * - 2005-10-28
 * 						- Add some status definitions for algorithms.
 * - 2005-06-14
 * 						- Add strdup definition, just for cleanliness when compiling
 * - 2005-05-23
 * 						- Add EGcontainerOf
 * - 2005-05-03
 * 						- Add typeof definition;
 * - 2004-07-14
 * 						- Add GNU_MP_Z GNU_MP_F and GNU_MP_Q to the type definitions.
 * - 2004-07-12
 * 						- Add EGRAT_TYPE to the type definitions.
 * - 2004-03-17
 * 						- Add TESTG that if recives something that is nonzero print an
 * 							user message and the go to the given location.
 * - 2004-02-05
 * 						- Add CHECKRVALG that checks a return value, display a mesage,
 * 							and then perform a goto.
 * - 2003-12-01
 * 						- Add definition of a 'copy' function and its MP version.
 * - 2003-11-20
 * 						- Add PTRTEST that check if a pointer points to the first 64Kb of
 * 							memory internal memory. Althought such situation may happend 
 * 							(if we work in kernel-related stuff), it is usually an error 
 * 							when we try to access such a memory.
 * - 2003-09-08
 * 						- Add ADVTESTL
 * - 2003-07-10
 * 						- Add MESSAGEF, ADVCHECKRVAL
 * - 2003-07-02
 * 						- Add EGosGetData, EGosSetData, EGosGetOffset
 * - 2003-06-16
 * 						- Add EXITL macro
 * - 2003-06-06 
 * 						- Add TESTL macro to test conditions but only when the debug
 * 							level is at least some value
 * - 2003-05-22 
 * 						- Add EXITRVAL
 * - 2003-05-15 
 * 						- Add CHECKRVAL MESSAGE and WARNING macros.
 * - 2003-05-08 
 * 						- Add support for variadic macros for EXIT and TEST
 * 						- Define EGRAND_MAX for SUN and LINUX acordingly, this is becouse 
 * 							for some reason the value of RAND_MAX in SUN is not as 
 * 							specified in stdlib.h but rather 1 << 31 - 1. Still I am not 
 * 							sure about the maximum value of rand() on sun... will fix that 
 * 							later on to.
 * 						- Add a mesage macro, it only print if the debug level is as high
 * 							as required by the first field. Again the definition is
 * 							variadric and if the debug level is 0 we reduce the macro to 
 * 							the empty instruction.
 *
 * */
/** @{*/
/** @file 
 * */
/* ========================================================================= */

#ifndef __EG_MACROS_H__
#define __EG_MACROS_H__
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "eg_config.h"

/* ========================================================================= */
/** @brief We define the GNU C extension typeof if necesary. */
#ifndef typeof
#define typeof __typeof__
#endif

/* ========================================================================= */
/** @brief return the offset of a member inside a structure.
 * @param type the type of the containing structure.
 * @param member the name of the member that we are interested in compute the
 * offset.
 * @return the number of bytes between the member and the beginning of the
 * structure. */
#define EGoffsetOf(type,member) ((size_t) &((type *)0)->member)

/* ========================================================================= */
/** @brief given a pointer to a member of a structure, return the pointer to
 * the head of the structure. (idea taken from Linux Kernel).
 * @param __ptr pointer to the member of the containing structure.
 * @param __type name type of the containing structure.
 * @param __member name of the given member in the containing structure.
 * @return pointer to the containing structure.
 * */
#define EGcontainerOf(__ptr,__type,__member) ({\
	typeof(((__type *)0)->__member) *const __EGcOf_ptr = (__ptr);\
	(__type *)( (char*)__EGcOf_ptr - ((size_t) &((__type *)0)->__member));})

/* ========================================================================= */
/** @name External C-library functions:
 * Here we define some C-library functions that for some weird reason can't 
 * be found at compile time but are present at linking time */
/*@{*/
#if OS != CYGWIN
extern void srandom (unsigned int);
extern long random (void);
extern double drand48 (void);
extern long lrand48 (void);
extern char *optarg;
extern int getopt (int,
								 char *const *,
									 const char *);
extern int finite (double);
#endif
/*@}*/

/* ========================================================================= */
/** @name Code Location Utility:
 * this are utility macros to print information about where we are.
 * @note For some reason __func__ don't work correctly for sun's cc in inline
 * functions, so we don't use in SUN architecture */
/* @{ */
#define __EG_PRINTLOCF__(F)  fprintf(((F==0)?stderr:F),", in %s (%s:%d)\n",__func__,__FILE__,__LINE__)
#define __EG_PRINTLOC__      __EG_PRINTLOCF__(stderr)
#define __EG_PRINTLOC2F__(F) fprintf(((F==0)?stderr:F),"in %s (%s:%d)\n",__func__,__FILE__,__LINE__)
#define __EG_PRINTLOC2__     __EG_PRINTLOC2F__(stderr)
/* @} */

#if DEBUG>=1
/* ========================================================================= */
/** @brief This macro is to print error messages and to return with value one
 * from the current function, it also print the file and line where this 
 * happend, but the condition is looked only if the debug level is at least L
 * */
#define EXITL(L,A,...) ({\
	if(L<=DEBUG){\
	if(A){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		exit(1);}}})

/* ========================================================================= */
/** @brief This macro is to print error messages and to return with value one 
 * from the current function, it also print the file and line where this 
 * happend, but the condition is looked only if the debug level is at least L 
 * */
#define TESTL(L,A,...) ({\
	if(L<=DEBUG){\
	if(A){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		return 1;}}})

/* ========================================================================= */
/** @brief this macro check if the value of a pointer is not bellow the first 
 * 64Kb, if so it return the given value */
#define PTRTEST(ptr,rval) {\
	if(ptr) ADVTESTL(0,((size_t)(ptr)) < (1U<<16),rval, \
								 "%s=%p is not a valid pointer",\
									#ptr, (void*)(ptr));}

/* ========================================================================= */
/** @brief This macro is to print error messages and to return with value one 
 * from the current function, it also print the file and line where this 
 * happend */
#define TESTG(A,B,...) ({\
	if(A){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		goto B;}})

/* ========================================================================= */
/** @brief This macro is to print error messages and to return with value one 
 * from the current function, it also print the file and line where this 
 * happend */
#define TEST(A,...) ({\
	if(A){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		return 1;}})

/* ========================================================================= */
/** @brief This macro print messages to the screen when the debug level is as 
 * big as the first parameter, if the debug level is zero we eliminate the 
 * code and reduce it to the empty instruction. */
#define MESSAGEF(A,F,...) ({\
	if(A <= DEBUG ){\
		fprintf(((F==0)?stderr:F),__VA_ARGS__);\
		__EG_PRINTLOCF__(F);}})

/* ========================================================================= */
/** @brief This macro print messages to the screen when the debug level is as 
 * big as the first parameter, if the debug level is zero we eliminate the 
 * code and reduce it to the empty instruction. */
#define MESSAGE(A,...) ({\
	if(A <= DEBUG ){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;}})

#if VERBOSE_LEVEL >= 1
/* ========================================================================= */
/** @brief This macro print messages to the screen when the verbose level is as 
 * big as the first parameter, if the verbose level is zero we eliminate the 
 * code and reduce it to the empty instruction. */
#define OUTPUT(A,...) ({\
	if(A <= VERBOSE_LEVEL ){\
		fprintf(stderr,__VA_ARGS__);}})
#else
#define OUTPUT(A,...) ;
#endif

/* ========================================================================= */
/** @brief This macro print messages to the screen when the condition A is 
 * true .if the debug level is one we don't print any warning message. if 
 * the debug level is zero we eliminate the code and reduce it to the empty 
 * instruction. */
#define WARNINGL(L,A,...) ({\
	if((A)&&(DEBUG>=L)){\
		fprintf(stderr,"WARNING: ");\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;}})

#else
#define TESTL(L,A,...) ;
#define EXITL(L,A,...) ;
#define TEST(A,...) ;
#define TESTG(A,B,...) ;
#define MESSAGE(A,...) ;
#define MESSAGEF(A,F,...) ;
#define WARNINGL(L,A,...) ;
#define PTRTEST(ptr,rval) ;
#endif

/* ========================================================================= */
/** @brief This macro is to print error messages and to return with value one 
 * from the current function, it also print the file and line where this 
 * happend */
#define FTEST(A,...) ({\
	if(A){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		return 1;}})

/* ========================================================================= */
/** @brief This macro is to print error messages and to return with value one 
 * from the current function, it also print the file and line where this 
 * happend */
#define FTESTG(A,B,...) ({\
	if(A){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		goto B;}})

/* ========================================================================= */
/** @brief This macro print messages to the screen when the condition A is 
 * true. */
#define WARNING(A,...) ({if(A){\
		fprintf(stderr,"WARNING: ");\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;}})

/* ========================================================================= */
/** @brief this macro test if a value is non zero, if it is it print where is 
 * it and exit 1. The idea is to use it to check return values of functions, 
 * and the calling function can't return a status, and then we are forced to 
 * exit. */
#define EXITRVAL(A) ({\
	if(A){\
		__EG_PRINTLOC2__;\
		exit(1);}})

/* ========================================================================= */
/** @brief This macro is to print error messages and exit the program with 
 * code one from the current function, it also print the file and line where 
 * this happend */
#define EXIT(A,...) ({if(A){\
		fprintf(stderr,"EXIT: ");\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		exit(1);}})

/* ========================================================================= */
/** @brief this macro test if a value is non zero, if it is it print where is 
 * it and return B. The idea is to use it to check return values of functions 
 * */
#define ADVCHECKRVAL(A,B) ({\
	if(A){\
		__EG_PRINTLOC2__;\
		return B;}})

/* ========================================================================= */
/** @brief This macro test a condition 'cond' when the debug level used at 
 * compile time is at least 'level'. If the condition is true, it print the 
 * message and return the 'rval' value. */
#define ADVTESTL(level,cond,rval,...) ({\
	if((DEBUG>=level)&&(cond)){\
		fprintf(stderr,__VA_ARGS__);\
		__EG_PRINTLOC__;\
		return rval;}})

/* ========================================================================= */
/** @brief this macro test if a value is non zero, if it is it print where is 
 * it and return 1. The idea is to use it to check return values of functions 
 * */
#define CHECKRVAL(A) ({\
	if(A){\
		__EG_PRINTLOC2__;\
		return A;}})

/* ========================================================================= */
/** @brief, if a non-zero value is given as an argument, check the errno stored
 * in the system, print the related message, and return the non-zero given
 * parameter, otherwise, do nothing.
 * @param __value if non-zero check systems errors, and return this value
 * */
#define TESTERRNOIF(__value) do{\
	if(__value){\
		const int __EGserrno = errno;\
		fprintf(stderr,"failed with errno %d, %s\n",__EGserrno, strerror(__EGserrno));\
		__EG_PRINTLOC2__;\
		return __value;}}while(0)
/* ========================================================================= */
/** @brief this function, if the input is non zero, print a message of 
 * function, file and line and then goto the second parameter */
#define CHECKRVALG(A,B) ({\
	if(A){\
		__EG_PRINTLOC2__;\
		goto B;}})

/* ========================================================================= */
/** @brief Define the real rand_max value of (random). In linux machines is 
 * as RAND_MAX, but in SUN is 2^31-1 */
#if OS == LINUX || OS == OSXMAC || OS == CYGWIN || OS == AIX
#define EGRAND_MAX RAND_MAX
#endif
#if OS == SUN
#define EGRAND_MAX ( (1LL << 31) - 1 )
#endif
#ifndef EGRAND_MAX
#error You have to specify the architecture, either SUN or LINUX are supported so far
#endif

/* ========================================================================= */
/** @brief retrieve the data of type 'type' in the structure 'data' that is 
 * located in the offset 'osN'. */
#define EGosGetData(data,osN,type) (*((type*)(((char*)data)+osN)))

/* ========================================================================= */
/** @brief set the data of type 'type' in the structure 'data' that is 
 * located in the offset 'osN' to the value 'val'. */
#define EGosSetData(data,osN,type,val) (EGosGetData(data,osN,type)=val)

/* ========================================================================= */
/** @brief Defione copy functions, these functions
 * return copy of objects but with independent storage space, there are two
 * versions, one that require a memory pool from where to look for memory, and
 * another where we don't care about that.... the place from where the memory
 * was asked for depend on the function, se the function definition for 
 * details.
 * Note that if the is no more memory available the function should call
 * exit(1).
 * This is only intended as a readibility help */
typedef void *(*EGcopy_f) (void *p);

/* ========================================================================= */
/** @brief Define a null copy function */
#define nullCopy ((EGcopy_f)0)

/* ========================================================================= */
/** @name Algorithms Return Status
 * Here we define some general status for algorithms, the exact meaning should
 * be sought in the actual algorithm definition, but the definitions here
 * provide a first overview of their meaning. */
/*@{*/
/** @brief the algorithm finish successfully. */
#define EG_ALGSTAT_SUCCESS 0
/** @brief the algorithm could only partially finish */
#define EG_ALGSTAT_PARTIAL 1
/** @brief the algorithm stop because of some numerical problem */
#define EG_ALGSTAT_NUMERROR 2
/** @brief the algorithm stop because of some unforeseen error */
#define EG_ALGSTAT_ERROR 3
/*@}*/

/* ========================================================================= */
/** @name Mathematical Constants 
 * Here we define some mathematical constants needed in some parts of the code
 * that are of general use */
/*@{*/
/** @brief definition of \f$\pi\f$ as a constant, suitable for quad-IEEE
 * operations. */
#define EG_M_PI 3.1415926535897932384626433832795029L
/*@}*/

/* ========================================================================= */
/** @brief hold the version of the libray */
extern const char EGversionString[1024];
/* ========================================================================= */
/** @brief hold the build date of the libray */
extern const char EGdate[1024];

/* ========================================================================= */
/** @} */
/* end of eg_macros.h */
#endif
