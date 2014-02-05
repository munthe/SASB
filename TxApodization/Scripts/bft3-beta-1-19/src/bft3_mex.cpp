/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: bft3_mex.cpp,v 1.3 2011-04-27 20:35:27 jmh Exp $  *
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-27 20:35:27 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *
 ****************************************************************************/

#include <cstdio>

#include "common.h"

#include "mex_utility.h"

#include "mexarg.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <cmath>

#define Usage   "Usage error. see above"

#include "aperture_mex.h"
#include "apodization_mex.h"
#include "line_mex.h"
#include "image_mex.h"
#include "sampled_image_mex.h"

static void bft3_mex_help(void) {
  printf("Usage for bft3_mex:\n");
}

bool bft3_mex(const int nlhs, mxArray *plhs[], const int nrhs,
							const mxArray *prhs[]) {
  
  char type[256];

  /* Check number of arguments */    
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    bft3_mex_help();
    Fail(Usage);
  }

  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strncmp(type, "line,",5)) {
		line_mex(nlhs,plhs,nrhs,prhs);
  }
  else if (!strncmp(type, "image,",6)) {
		image_mex(nlhs,plhs,nrhs,prhs);
  }
  else if (!strncmp(type, "sampled_image,",14)) {
		sampled_image_mex(nlhs,plhs,nrhs,prhs);
  }
  else if (!strncmp(type, "aperture,",9)) {
		aperture_mex(nlhs,plhs,nrhs,prhs);
  }
  else if (!strncmp(type, "apodization,",12)) {
		apodization_mex(nlhs,plhs,nrhs,prhs);
  }
  else {
    bft3_mex_help();
    Fail(Usage);
  }
 
  return true;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (!nlhs && !nrhs) {
    bft3_mex_help();
    return;
  }
  if (!bft3_mex(nlhs, plhs, nrhs, prhs)) {
    mexErrMsgTxt("bft3_mex()");
  }
}

#undef Usage

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
