/*****************************************************************************
*                                                                            *
* Project              PILOT BEAM FORM                                       *
* Module                                                                     *
*                                                                            *
* $Id: mexarg.h,v 1.21 2011-04-19 14:38:18 jmh Exp $              *
*                                                                            *
* $Author: jmh $                                                             *
*                                                                            *
* $Date: 2011-04-19 14:38:18 $                                               *
*                                                                            *
* $State: Exp $                                                              *
*----------------------------------------------------------------------------*
*                                                                            *
* Revision 1.14  2010/03/18 15:07:30  jmh
* Removed a lot of crap
*
*
*****************************************************************************/

/* TODO: C++ version */

#ifdef _MSC_VER
 #pragma warning( disable : 4996 ) // suppress warning: 'sprintf' was declared deprecated
 #define _CRT_SECURE_NO_WARNINGS
#endif

/* include this instead */
#include "mex_macro.h"

#if defined(MATLAB_MEX_FILE)
#include "mex.h"

#ifdef __cplusplus
 extern "C" {
#endif

extern char *mx_string(const mxArray *mx, const char *);
extern bool mx_string_free(char *);
extern bool mexarg(const int, const mxArray *[]);


char *mxu_string(const mxArray *mx, const char *arg);
bool mxu_string_free(char *arg);
bool mxu_fixed_string(char* buffer, size_t length, const mxArray *mx, const char *arg);
bool mxCheckDimensions(const mxArray* mx_array, size_t n_dim,...);

#define mxGetInt(mx) (*((const int*) mxGetData(mx)))
#define mxGetUInt64(mx) (*((const unsigned long long int *) mxGetData(mx)))
#define mxGetInt64(mx) (*((const long long int *) mxGetData(mx)))

   /* If error here const unsigned long int together with use of size_t for trans_elem_no, we get stack corruption */
   /* Possible error on 64-bit platform */
#define mxGetUInt32(mx) (*((const unsigned int *) mxGetData(mx)))

#define mxIsScalar(mx) \
	( (2 == mxGetNumberOfDimensions(mx)) \
		&& (1 == mxGetM(mx)) && (1 == mxGetN(mx)) )

#define mxIsScalarInt64(mx) \
	( mxIsScalar(mx) && mxIsInt64(mx) )
#define mxIsScalarUInt64(mx) \
	( mxIsScalar(mx) && mxIsUint64(mx) )
#define mxIsScalarInt32(mx) \
	( mxIsScalar(mx) && mxIsInt32(mx) )
#define mxIsScalarUInt32(mx) \
	( mxIsScalar(mx) && mxIsUint32(mx) )
#define mxIsComplexSingle(mx) \
	(mxIsSingle(mx) && mxIsComplex(mx))
#define mxIsComplexDouble(mx) \
	(mxIsDouble(mx) && mxIsComplex(mx))
#define mxIsRealSingle(mx) \
	(mxIsSingle(mx) && !mxIsComplex(mx))
#define mxIsRealDouble(mx) \
	(mxIsDouble(mx) && !mxIsComplex(mx))
#define mxIsScalarSingle(mx) \
	( mxIsScalar(mx) && mxIsRealSingle(mx) )
#define mxIsScalarDouble(mx) \
	( mxIsScalar(mx) && mxIsRealDouble(mx) )

/*
#ifndef mxIsRealFloat
# define mxIsRealFloat(mx) mxIsRealDouble(mx)
#endif

#ifndef mxIsFloat
# define mxIsFloat(mx) mxIsDouble(mx)
#endif
*/

#ifdef __cplusplus
 }
#endif

#endif

