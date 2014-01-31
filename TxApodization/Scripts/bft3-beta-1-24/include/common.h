/* $Id: common.h,v 1.37 2011-07-26 23:36:20 jmh Exp $ */
#ifndef COMMON_H
#define COMMON_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#else
#define PACKAGE_VERSION "@PACKAGE_VERSION@"
#endif

#if defined(_WIN32) || defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
#undef HAVE_INTTYPES_H
#endif

// TODO: Move to elsewhere
#ifdef _MSC_VER
 #pragma warning( disable : 4996 ) // suppress warning: 'sprintf' was declared deprecated
 #define _CRT_SECURE_NO_WARNINGS
#endif

#if defined(HAVE_STDDEF_H) || defined(_MSC_VER)
# ifndef _SIZE_T_DEFINED
#  include <cstddef>
# endif
#endif

#ifdef HAVE_INTTYPES_H
# define __STDC_FORMAT_MACROS 1
# define __STDC_LIMIT_MACROS 1
# include <inttypes.h> // <cinttype> requires -std=c++0x
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
  // 16 bit systems use long int for 32 bit integer
  typedef   signed long int int32_t;
  typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
  // Microsoft have their own definition
  typedef   signed __int32  int32_t;
  typedef unsigned __int32 uint32_t;
  typedef   signed __int64  int64_t;
  typedef unsigned __int64 uint64_t;
# ifndef SIZE_MAX
#  ifdef  _WIN64
#   define SIZE_MAX             (18446744073709551615UL)
#  else
#   define SIZE_MAX             (4294967295U)
#  endif
# endif
#else
  // This works with most compilers
  typedef signed int          int32_t;
  typedef unsigned int       uint32_t;
  typedef long long           int64_t;
  typedef unsigned long long uint64_t;
#endif

#ifdef HAVE_MQUEUE_H
# include <mqueue.h>
#endif

// Increase if necessary
#define N_MAX_THREADS 24

#ifndef mxPOINTER_CLASS
# define mxPOINTER_CLASS mxUINT64_CLASS
#endif

#ifndef mxIsPointer
# define mxIsPointer mxIsScalarUInt64
#endif

#ifndef mxPOINTER_TYPE
# define mxPOINTER_TYPE UINT64_T
#endif

#ifndef mxIsPointerArray
# define mxIsPointerArray mxIsUint64
#endif

#ifdef M_SINGLE_PRECISION
# define float_type float
# define mxIsRealFloat mxIsRealSingle
# define mxFLOAT_CLASS mxSINGLE_CLASS
# define mxIsScalarFloat mxIsScalarSingle
# ifndef mxIsFloat
#  define mxIsFloat(mx) mxIsSingle(mx)
# endif
#else
# define float_type double
# define mxFLOAT_CLASS mxDOUBLE_CLASS
# define mxIsScalarFloat mxIsScalarDouble
# ifndef mxIsRealFloat
#  define mxIsRealFloat mxIsRealDouble
# endif
# ifndef mxIsFloat
#  define mxIsFloat(mx) mxIsDouble(mx)
# endif
#endif

#define N_MAX_RADIAL_SAMPLE 10000

#if (defined(_MSC_VER) && defined(_WIN32))
# define ALIGN16_BEGIN __declspec(align(16))
# define ALIGN16_END 
#elif defined(__GNUC__)
# define ALIGN16_BEGIN
# define ALIGN16_END __attribute__((aligned(16)))
#endif

#if (defined(_MSC_VER) && defined(_WIN32))
# define STATIC_INLINE_BEGIN static __forceinline 
# define STATIC_INLINE_END
# define INLINE_BEGIN __forceinline 
# define INLINE_END
#elif (defined(__GNUC__))
# define STATIC_INLINE_BEGIN static inline
# define STATIC_INLINE_END __attribute__ ((always_inline))
# define INLINE_BEGIN inline
# define INLINE_END __attribute__ ((always_inline))
#endif

#ifndef SQUARE
#define SQUARE(z) (z) * (z)
#endif

#define CallErr(fun, arg)  { if ((fun arg)<0)          \
   FailErr(#fun) }

#define CallErrExit(fun, arg, ret)  { if ((fun arg)<0) {      \
   FailErr(#fun);                                             \
   return ret; }}

#define FailErr(msg) {                                                 \
  (void)fprintf(stderr, "FAILED: %s %d: %s (errno=%d, strerror=%s)\n", \
                __FILE__, __LINE__, msg, errno, strerror(errno));      \
  (void)fflush(stderr);}


#endif // end COMMON_H
