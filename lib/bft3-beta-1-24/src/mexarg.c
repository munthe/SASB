/*****************************************************************************
*                                                                            *
* Project              SIMD BEAMFORM MEX                                     *
* Module               Mex routines                                          *
*                                                                            *
* $Id: mexarg.c,v 1.16 2011-04-16 13:19:36 jmh Exp $  *
*                                                                            *
* $Author: jmh $                                                             *
*                                                                            *
* $Date: 2011-04-16 13:19:36 $                                               *
*                                                                            *
* $State: Exp $                                                              *
*----------------------------------------------------------------------------*
*                                                                            *
* $Log: mexarg.c,v $
* Revision 1.16  2011-04-16 13:19:36  jmh
* *** empty log message ***
*
* Revision 1.15  2011/04/08 11:52:19  jmh
* *** empty log message ***
*
* Revision 1.11  2010/03/21 12:20:27  jmh
* Removed crappy fs,f0, and c
*
* Revision 1.9  2010/03/20 20:31:22  jmh
* Dynamic and manual working
*
* Revision 1.8  2010/03/20 13:52:12  jmh
* 2nd attempt dyn emit
*
* Revision 1.7  2010/03/19 19:57:09  jmh
* First attempt dyn xmt apo
*
* Revision 1.6  2010/01/29 11:54:35  jmh
* Something wrong
*
*****************************************************************************/

#include "mexarg.h"
#include <stdarg.h>

#ifdef _WIN32

#else
# include<inttypes.h>
#endif

#if defined(MATLAB_MEX_FILE)

char *mxu_string(const mxArray *mx, const char *arg) {
  char *cbuffer;
  size_t n = mxGetM(mx) * mxGetN(mx) + 1;

  if (!mxIsChar(mx))
    FailM("%s must be char array", arg)

  Call(cbuffer = (char *) mxCalloc, (n, sizeof(char)))

  if (mxGetString(mx, cbuffer, (int)n))
    Warn("bug with mxGetString")
  return cbuffer;
}

bool mxu_string_free(char *arg) {
  mxFree((void*)arg);
  return true;
}

bool mxu_fixed_string(char* buffer, size_t length, const mxArray *mx, const char *arg) {

  size_t n = mxGetM(mx) * mxGetN(mx) + 1;

  if (!mxIsChar(mx))
    FailM("%s must be char array", arg)

  /* PRIu64*/
  if (n > length)
#ifdef _WIN32
    FailM("String too large: length is %lu > %lu", n, length)
#else
    FailM("String too large: length is %zu > %zu", n, length)
#endif
  if (mxGetString(mx, buffer, (int)n))
    Warn("bug with mxGetString")

  return true;
}


/*
 * mx_string()
 * caller must free using mx_string_free()
*/
char *mx_string(const mxArray *mx, const char *arg)
{
  char	*string;
  /*  int n = (int) mxGetM(mx) * mxGetN(mx) + 1; */
  size_t n = mxGetM(mx) * mxGetN(mx) + 1;

  if (!mxIsChar(mx))
    FailM("%s must be char array", arg)
      
#if 1
      Call(string = (char *) mxCalloc, (n, sizeof(char)))
#else
      Mem0(string, n)
#endif
      if (mxGetString(mx, string, n))
	Warn("bug with mxGetString")
	  return string;
}

bool mx_string_free(char *s)
{
  mxFree(s);
  return true;
}

bool mxCheckDimensions(const mxArray* mx_array, size_t n_dim,...) {

  va_list ap;             /* varargs list traverser */
  
  size_t *dims;           /* dimension list */
  size_t i;
  size_t dim;
  bool retval = true;

  va_start(ap,n_dim);
  dims = (size_t *) malloc(n_dim*sizeof(size_t));
  
  for(i=0;i<n_dim;i++) {
    dims[i] = va_arg(ap,size_t);
    dim  = mxGetDimensions(mx_array)[i];
    if (dim != dims[i])
      retval = false;
  }

  va_end(ap);
  free(dims);

  return retval;
}


static bool mex_showdim(const mxArray *mx)
{
  /*	int ndim;*/
  mwSize ndim;
  const mwSize *dims;
  Call(ndim = mxGetNumberOfDimensions, (mx));
  Call(dims = mxGetDimensions, (mx));
  if (ndim == 1)
    NoteM("dims %d", (int)dims[0]);
  if (ndim == 2)
    NoteM("dims %d %d", (int)dims[0], (int)dims[1]);
  if (ndim == 3)
    NoteM("dims %d %d %d", (int)dims[0], (int)dims[1], (int)dims[2]);
  if (ndim == 4)
    NoteM("dims %d %d %d %d", (int)dims[0], (int)dims[1],
	  (int)dims[2], (int)dims[3]);
  if (ndim > 4)
    NoteM("dims %d %d %d ... %d", (int)dims[0], (int)dims[1],
	  (int)dims[2], (int)dims[ndim-1]);
  return EXIT_SUCCESS;
}


/*
* mexarg()
* Show arguments
*/
bool mexarg(const int nmx, const mxArray *pmx[])
{
  int ii;
  
  NoteM("narg=%d", nmx);
    
  for (ii=0; ii < nmx; ++ii) {
    const mxArray *mx = pmx[ii];
    
    if (mxIsChar(mx)) {
      char *arg;
      Call(arg = mx_string, (mx, ""));
      NoteM("arg %d, char, '%s'", ii, arg);
      Call(mx_string_free, (arg));
    }
    
    else if (mxIsInt32(mx)) {
      int val = *((int *) mxGetData(mx));
      NoteM("arg %d, int32, val[0]=%d", ii, val);
    }
      
    else if (mxIsSingle(mx)) {
      float val = *((const float *) mxGetData(mx));
      NoteM("arg %d, single, val[0]=%g", ii, val);
    }
      
    else if (mxIsDouble(mx)) {
      double val = *((const double *) mxGetData(mx));
      NoteM("arg %d, double, val[0]=%g", ii, val);
    }
    
    else
      NoteM("arg %d UNKNOWN!?", ii);
    
    Call(mex_showdim, (mx));
  }
  
  return EXIT_SUCCESS;
}

#endif
 
