/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: line_mex.cpp,v 1.37 2011-07-27 21:04:59 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-27 21:04:59 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 ****************************************************************************/

#include <cstdio>

#include "common.h"

#include "mex_utility.h"
#include "mexarg.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <math.h>

#define Usage   "Usage error. see above"

#include "common.h"

#include "aperture.h"
#include "apodization.h"
#include "line.h"
#include "line_mex.h"

static bool enable_line_destructor = false;

void line_mex_help(void) {
  printf("Information\n");
}

bool line_bft_ctor_manual(mxArray*& plhs,
                                 const mxArray* mx_origin,
                                 const mxArray* mx_direction,
                                 const mxArray* mx_dr,
                                 const mxArray* mx_length) {

  const float_type* origin, *direction;
  float_type dr, length;
  mwSize no_lines;
  Line* p_line;

  // Only works for single line
  Call(mxIsRealFloat, (mx_origin));
  Call(mxIsRealFloat, (mx_direction));
  Call(mxIsRealFloat, (mx_dr));
  Call(mxIsRealFloat, (mx_length));

  no_lines = mxGetM(mx_origin);
  
  Call1(mxCheckDimensions, (mx_origin, 2, no_lines, 3),
        "Positions must be 3D");
  origin = ((const float_type *) mxGetData(mx_origin));
      
  Call1(mxCheckDimensions, (mx_direction, 2, no_lines, 3),
        "Number of origins, directions, dr, and length must be equal");
  direction = ((const float_type *) mxGetData(mx_direction));

  Call1(mxCheckDimensions, (mx_dr, 2, no_lines, 1), "Wrong dimension: dr");
  dr = *((float_type *) mxGetData(mx_dr));

  Call1(mxCheckDimensions, (mx_length, 2, no_lines, 1),
        "Wrong dimension: length");
  length = *((float_type *) mxGetData(mx_length));

  // TODO: Verify unit vector
  p_line = new Line(origin, direction, dr, length);

  plhs = create_handle<Line>(p_line);

  enable_line_destructor = true;
  mexLock();
  return true;
}

bool line_bft_dtor(mxArray *plhs[], const mxArray* mx_handle) {

  if (!enable_line_destructor)
    return true;
  
  Call(mxIsPointer, (mx_handle));
  destroy_object<Line>(mx_handle);
  mexUnlock();
  return true;
}

/* Added extra argument, default is NULL */
bool line_bft_get(mxArray *plhs[], const mxArray* mx_handle,
                         const size_t type, const mxArray* mx_data) {

  const int n_odim = 2;
  mwSize o_dims[n_odim];
  size_t i,j;
  Line* p_line;
  float_type* p_flhs = NULL;
  Apodization* p_apodization = NULL;
	Aperture* p_aperture = NULL;
	bool retval = true;
	size_t i_element = 0;

  Call(mxIsPointer, (mx_handle));
  p_line = &(get_object<Line>(mx_handle));

	if (mx_data) {
    Call(mxIsScalarUInt32, (mx_data));
	}

  switch (type) {
  case 0:
    // Construct output
    o_dims[0] = mwSize(p_line->data->m_npos);
    o_dims[1] = mwSize(3);
    
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));

    for (j=0;j<3;j++) {
      for (i=0;i<p_line->data->m_npos;i++) {
        p_flhs[i+j*p_line->data->m_npos] =
          p_line->data->m_origin[j] + 
          p_line->data->m_direction[j]*p_line->data->m_dr*i; 
      }
    }
    break;
  case 1:
    // XMT Apodization
    if (!p_line->data->m_xmt_apodization.data) {
      Call(plhs[0] = mxCreateString, ("unassigned"));
    }
    else {
      // Call(plhs[0] = mxCreateString, ("not supported yet"));
      p_apodization = new Apodization();
      *p_apodization = p_line->data->m_xmt_apodization;
      plhs[0] = create_handle<Apodization>(p_apodization);    
    }
    break;
  case 2:
    // RCV Apodization
    if (!p_line->data->m_rcv_apodization.data) {
      Call(plhs[0] = mxCreateString, ("unassigned"));
    }
    else {
      p_apodization = new Apodization();
      *p_apodization = p_line->data->m_rcv_apodization;
      plhs[0] = create_handle<Apodization>(p_apodization);    
    }
    break;
  case 3:
		// XMT Apodization values
		if (p_line->data->m_xmt_apodization.data) {

			o_dims[0] = mwSize(1);
			o_dims[1] = mwSize(p_line->data->m_npos);

			i_element = (size_t) mxGetUInt32(mx_data);
			i_element--;

			if (!(i_element < p_line->data->m_xmt_apodization.data->m_aperture.data->m_npos)) {
				Fail("Illegal emission index");
			}
				
			Call(plhs[0] = mxCreateNumericArray,
					 (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

			Call(p_flhs = (float_type*) mxGetData,(plhs[0]));
			// Retrieve transmit apodization values
			p_line->apodize(p_flhs,0,i_element);
		}
		else {
			retval = false;
		}
		break;
	case 4:
		// RCV Apodization values
		if (p_line->data->m_rcv_apodization.data) {
			p_aperture = &(p_line->data->m_rcv_apodization.data->m_aperture);

			o_dims[0] = mwSize(p_aperture->data->m_npos);
			o_dims[1] = mwSize(p_line->data->m_npos);
    
			Call(plhs[0] = mxCreateNumericArray,
					 (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

			Call(p_flhs = (float_type*) mxGetData,(plhs[0]));
			// Retrieve receive apodization values
			p_line->apodize(p_flhs,1);
		}
		else {
			retval = false;
		}
		break;
	case 5:
		// Origin
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(3);
    
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    p_flhs = ((float_type*) mxGetData(plhs[0]));
		memcpy(p_flhs,p_line->data->m_origin,3*sizeof(float_type));
		break;
	case 6:
		// Direction
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(3);
    
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    p_flhs = ((float_type*) mxGetData(plhs[0]));
		memcpy(p_flhs,p_line->data->m_direction,3*sizeof(float_type));
		break;
	case 7:
		// Dr
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    p_flhs = ((float_type*) mxGetData(plhs[0]));
		p_flhs[0] = p_line->data->m_dr;
		break;
  default:
    break;
  }
  return retval;
}

bool line_bft_set(mxArray *plhs[], const mxArray* mx_handle,
												 const mxArray* mx_data, const size_t type) {

  bool retval = false;
  const int n_odim = 1;
  mwSize o_dims[n_odim];

  Line* p_line = NULL;

  Call(mxIsPointer, (mx_handle));
  p_line = &(get_object<Line>(mx_handle));

  o_dims[0] = mwSize(1);
  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  if (type < 2)
		Call(mxIsPointer, (mx_handle));

	switch (type) {
	case 0:
		// XMT apodization
		p_line->data->m_xmt_apodization = get_object<Apodization>(mx_data);
		retval = true;
		break;
	case 1:
		// RCV apodization
		p_line->data->m_rcv_apodization = get_object<Apodization>(mx_data);
		retval = true;
		break;
	default:
		retval = false;
	}

  bool* p_blhs = ((bool*) mxGetData(plhs[0]));
  *p_blhs = retval;

  return retval;
}

extern bool line_mex(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
  
  char type[256];
  char* attribute = NULL;

  /* Check number of arguments */
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    line_mex_help();
    Fail(Usage);
  }

  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strcmp(type, "line,ctor,manual")) {
    if ((nrhs != 5) || (nlhs != 1)) {
      line_mex_help();
      Fail(Usage);
    }
    Call(line_bft_ctor_manual, (plhs[0], prhs[1], prhs[2], prhs[3], prhs[4]));
  }
  else if (!strcmp(type, "line,dtor")) {
    if ((nrhs != 2) || (nlhs != 0)) {
      line_mex_help();
      Fail(Usage);
    }
    Call(line_bft_dtor, (plhs, prhs[1]));
  } // != 2
  else if (!strncmp(type, "line,get,",9)) {
    attribute = &type[9];
    if ((nrhs < 2) || (nlhs != 1) || (strlen(type)<10)) {
      line_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute, "pos")) {
      Call(line_bft_get, (plhs, prhs[1],0));
    }
    else if (!strcmp(attribute, "xmt_apodization")) {
      Call(line_bft_get, (plhs, prhs[1],1));
    }
    else if (!strcmp(attribute, "rcv_apodization")) {
      Call(line_bft_get, (plhs, prhs[1],2));
    }
    else if (!strcmp(attribute, "xmt_apodization_values")) {
      Call(line_bft_get, (plhs, prhs[1],3,prhs[2]));
    }
    else if (!strcmp(attribute, "rcv_apodization_values")) {
      Call(line_bft_get, (plhs, prhs[1],4));
    }
    else if (!strcmp(attribute, "origin")) {
      Call(line_bft_get, (plhs, prhs[1],5));
    }
    else if (!strcmp(attribute, "direction")) {
      Call(line_bft_get, (plhs, prhs[1],6));
    }
    else if (!strcmp(attribute, "dr")) {
      Call(line_bft_get, (plhs, prhs[1],7));
    }
    else {
      line_mex_help();
      Fail(Usage);
    }
  }
  else if (!strncmp(type, "line,set,",9)) {
    attribute = &type[9];
    if ((nrhs != 3) || (nlhs != 0) || (strlen(type)<10)) {
      line_mex_help();
      Fail(Usage);
    }
		else if (!strcmp(attribute, "xmt_apodization")) {
			Call(line_bft_set, (plhs, prhs[1],prhs[2],0));
		}
		else if (!strcmp(attribute, "rcv_apodization")) {
			Call(line_bft_set, (plhs, prhs[1],prhs[2],1));
		}
    else {
      line_mex_help();
      Fail(Usage);
    }
  }
  return true;
}

#if defined(Need_mex_gateway)

/** 
 * Gateway routine
 * 
 * @param nlhs 
 * @param plhs 
 * @param nrhs 
 * @param prhs 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (!nlhs && !nrhs) {
    line_mex_help();
    return;
  }
  if (!line_mex(nlhs, plhs, nrhs, prhs))
    mexErrMsgTxt("aperture_mex()");
}
#endif

#undef Usage

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
