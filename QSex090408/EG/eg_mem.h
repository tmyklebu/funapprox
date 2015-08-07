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
/** @defgroup EGmem EGmem
 *
 * Here we define some usefull macros to deal with memory issues, for example,
 * assert that we always return memory when posible, and if no memory is found,
 * then we just exit to the system (because if there is trully no memory....
 * there is no much else to do... unless we start using shrinkable memory
 * pools, like for example @ref EGmemSlab , but that is still a long way off,
 * it will also perform (if debugging enabled) some allocation / freeing
 * checkings and so on.
 *
 * @version 0.0.1
 * @par History:
 * -2005-09-05
 * 					- Add EGrealloc, wich is a wrapper of realloc but that assures us
 * 					to have memory, if there is no  memory, we exit. The idea of these
 * 					functions is that in the future they would interact with the
 * 					memory pools to use any memory still in the pools.
 * -2005-08-20
 * 					- Move memory align definitions here, and set the aligment of
 * 						memory to 8 bytes (i.e. 64 bits). This is to simplify compilation
 * 						in diferent architectures like Sun, opteron 64 and intel 32.
 * -2005-08-01
 * 					- Fix calloc call to the right type (size_t), and some printing
 * 						issues while compiling on 64-bit architectures.
 * -2005-07-30
 * 					- First Implementation
 * */
/** @file
 * @ingroup EGmem */
/** @addtogroup EGmem */
/** @{ */
#ifndef __EG_MEM_H__
#define __EG_MEM_H__
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include "eg_macros.h"

/* ========================================================================= */
/** @brief size of a normal word in this machine (a word is just big enough to
 * store a pointer) */
#define EG_MEM_WORD_SIZE (sizeof(void*))
/* ========================================================================= */
/** @brief memory aligment used by EG alloc functions. */
#define EG_MEM_ALIGNMENT 8U

/* ========================================================================= */
/** @brief \f$log_2(EG_MEM_ALIGNMENT)\f$. */
#define EG_MEM_ALIGNMENT_SHIFT 3U

/* ========================================================================= */
/** @brief Given a pointer, return it's aligned value. */
#define EG_MEM_ALIGN(ptr) \
	((((size_t)ptr)+EG_MEM_ALIGNMENT-1)&(~(EG_MEM_ALIGNMENT-1)))

/* ========================================================================= */
/** @brief type of the free functions that recive only one parameter */
typedef void (*EGfree_f) (void *);

/* ========================================================================= */
/** @brief this is the the data free that does nothing, use it when you don't 
 * want/need to free the internal list data becouse you will do it 
 * elsewere */
#define nullFree ((EGfree_f)0)

/* ========================================================================= */
/** @brief type for constructor functions. Given a pointer to an element of
 * some type, do the internal initialization necesary so that we can work with
 * the lement, such initialization may include allocating some internal memory
 * needed by the structure (not done by the user). This functions must never
 * fail. if some unexpected error does happen inside, then the function should
 * not return. (a call to exit(1) would do the trick). */
typedef void (*EGconstructor_f) (void *);

/* ========================================================================= */
/** @brief Null constructor function (do nothing) */
#define nullConstructor ((EGconstructor_f)0)

/* ========================================================================= */
/** @brief type for destructor functions. Given a pointer to an element of some
 * type, free all internal memory related to the element allocated during the
 * construction phase. (but not the pointer itself). This function must always
 * succed, if an error happen, the function should never return. (a call to
 * exit(1) would do the trick). */
typedef void (*EGdestructor_f) (void *);

/* ========================================================================= */
/** @brief Null destructor function (do nothing) */
#define nullDestructor ((EGdestructor_f)0)

/* ========================================================================= */
/** @brief this function replace malloc, check if the memory is not zero, if 
 * it is, it exit from the program, and display who called it and how much 
 * memory it tryed to alloc.
 * @param A number of bytes to allocate.
 * @return a void* pointer to the newly allocated memory, note that if the 
 * function returns at all, it will return with the amount of memory required, 
 * so no NULL checking is ever necesary after an EGmalloc call. 
 * */
#define EGmalloc(A) ({\
	size_t const _EGmp_sz_ = (A);\
	void * _EGmp_res_ = 0;\
	/*WARNINGL(0,!_EGmp_sz_,"Allocating 0 bytes");*/\
	if(_EGmp_sz_)\
	{\
		_EGmp_res_ = calloc((size_t)1,_EGmp_sz_);\
		EXIT(!_EGmp_res_,"Not enough memory while allocating %zd bytes",_EGmp_sz_);\
	}\
	_EGmp_res_;})

/* ========================================================================= */
/** @brief This function allocate 'count' elements of type 'type' and return 
 * a pointer of type 'type*'. If the memory is not available the program will 
 * exit indicating where it was trying to get memory and how much, it will also
 * check some common errors like allocating zero bytes.
 * @param type type of the element required.
 * @param count number of contiguous elements of the given type required.
 * @return pointer to the beggining of the allocated array of the apropiate
 * type (so no casting is needed). Note that if this function returns at all,
 * then the memory has been allocated and thus no NULL checking return is
 * necesary. */
#define EGsMalloc(type,count) (type*)EGmalloc(sizeof(type)*(count))

/* ========================================================================= */
/** @brief Realloc a given pointer to the new size, and check that we find
 * enough memory to return. If we don't, we exit the execution.
 * @param __ptr pointer to reallocate.
 * @param __sz new number of bytes to reallocate.
 * @return pointer to the new block of memory */
#define EGrealloc(__ptr,__sz) ({\
	const size_t ____sz = (size_t)(__sz);\
	(__ptr) = realloc((__ptr),____sz);\
	EXIT(!(__ptr)&&(____sz),"not enough memory while reallocating %zd",____sz);\
	(__ptr);})

/* ========================================================================= */
/** @brief this is used to enable malloc/free tracking and extra debugging */
#ifndef __EG_MEM_FREE_CHECK__
#define __EG_MEM_FREE_CHECK__  (1 && DEBUG)
#endif

/* ========================================================================= */
/** @brief This function replace free, the idea of this is to HOPEFULLY later 
 * develop a memory leack checker that tell us who and where asked for memory 
 * and didn't free it, in hte meantime they do nothing.
 * @param A pointer to the piece of memory to be freed, if debuging is enabled,
 * the function will test for freing NULL pointers, for suspicios address
 * freing and so on. note that the given pointer will point to NULL after this
 * call, thus reducing the posibility of freeing multiple times the same piece
 * of memory, or of allocating it after freeing it. */
#if __EG_MEM_FREE_CHECK__
#define EGfree(A) ({\
		EXIT(((A) && !(((size_t)(A))>>19)),"Trying to free pointer "#A\
				" with value %zd\nThis is probably an error",(size_t)(A));\
		if(A) free(A);\
		else WARNING(1,"Trying to free "#A", a NULL pointer");\
		(A) = 0;})
#else
#define EGfree(A) ({free(A);(A)=0;})
#endif

/* ========================================================================= */
/* end of eg_mem.h */
/** @} */
#endif
