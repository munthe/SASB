/******************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: sampled_image_mex.cpp,v 1.24 2011-07-25 18:07:51 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-25 18:07:51 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 ******************************************************************************/

/****************************************************************************
 * TODO:
 * - Check third dimension of data
 ****************************************************************************/

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
#include "sampled_image.h"
#include "sampled_image_mex.h"

static bool enable_sampled_image_destructor = false;

extern void sampled_image_mex_help(void)
{
  printf("Information\n");
}

// Classes only supported in mex-interface in later versions
#if MX_API_VER > 0x07030000

extern bool sampled_image_bft_ctor(mxArray*& plhs,
																	 const mxArray* mx_handle_Th_t,
																	 const mxArray* mx_handle_Th_r,
																	 const mxArray* mx_handle_Ah_t,
																	 const mxArray* mx_handle_Ah_r,
																	 const mxArray* mx_im_geom) {

  sampled_image_t image;
#ifndef _NDEBUG
  image.nx = 1; image.ny = 1; image.nz = 1;
  image.dx = float_type(1.0); image.dy = float_type(1.0);
  image.dz = float_type(1.0);
  image.offset_x = float_type(1.0); image.offset_y = float_type(1.0);
  image.offset_z = float_type(1.0);
#endif

  const mxArray* mx_data = NULL;
  bool retval = false;
                
  // TODO: Convert object to struct
  if (mxIsClass(mx_im_geom,"bft3_im_geom")) {
    retval = true;
    Call(mx_data = mxGetProperty, (mx_im_geom,0,"nx"));
    if (mx_data)
      image.nx = (size_t) (*((float_type*)mxGetData(mx_data)));

    Call(mx_data = mxGetProperty, (mx_im_geom,0,"nz"));
    if (mx_data) {
      if (!mxIsEmpty(mx_data))
				image.nz = (size_t) (*((float_type*)mxGetData(mx_data)));
      else
				Fail("nz can't be empty");
    }
    Call(mx_data = mxGetProperty, (mx_im_geom,0,"ny"));
    if (mx_data)
      if (!mxIsEmpty(mx_data))
				image.ny = (size_t) (*((float_type*)mxGetData(mx_data)));
    
    Call(mx_data = mxGetProperty, (mx_im_geom,0,"dx"));
    if (mx_data)
      image.dx = (*((float_type*)mxGetData(mx_data)));

    Call(mx_data = mxGetProperty, (mx_im_geom,0,"dz"));
    if (mx_data)
      image.dz = (*((float_type*)mxGetData(mx_data)));

    Call(mx_data = mxGetProperty, (mx_im_geom,0,"dy"));
    if (mx_data)
      if (!mxIsEmpty(mx_data))
				image.dy = (*((float_type*)mxGetData(mx_data)));

    Call(mx_data = mxGetProperty, (mx_im_geom,0,"offset_x"));
    if (mx_data)
      image.offset_x = (*((float_type*)mxGetData(mx_data)));

    Call(mx_data = mxGetProperty, (mx_im_geom,0,"offset_z"));
    if (mx_data)
      image.offset_z = (*((float_type*)mxGetData(mx_data)));

    Call(mx_data = mxGetProperty, (mx_im_geom,0,"offset_y"));
    if (mx_data)
      if (!mxIsEmpty(mx_data))
				image.offset_y = (*((float_type*)mxGetData(mx_data)));
  }

  Call(mxIsPointer, (mx_handle_Th_t));
  Call(mxIsPointer, (mx_handle_Th_r));

  Aperture aperture_Th_t = get_object<Aperture>(mx_handle_Th_t);
  Aperture aperture_Th_r = get_object<Aperture>(mx_handle_Th_r);

  Apodization apodization_Ah_t = get_object<Apodization>(mx_handle_Ah_t);
  Apodization apodization_Ah_r = get_object<Apodization>(mx_handle_Ah_r);

  // Construct image
  SampledImage* p_image =
		new SampledImage(aperture_Th_t, aperture_Th_r,
										 apodization_Ah_t,apodization_Ah_r,
										 image);
	
  plhs = create_handle<SampledImage>(p_image);

  enable_sampled_image_destructor = true;
  mexLock();

  return retval;
}

#endif

extern bool sampled_image_bft_dtor(mxArray *plhs[], const mxArray* mx_handle) {

  if (!enable_sampled_image_destructor)
    return true;
  
  Call(mxIsPointer, (mx_handle));
  destroy_object<SampledImage>(mx_handle);
  mexUnlock();
  return true;
}

extern bool sampled_image_bft_get(mxArray *plhs[],
																	const mxArray* mx_handle, const size_t type) {

  SampledImage* p_image;

  Call(mxIsPointer, (mx_handle));
  
  p_image = &(get_object<SampledImage>(mx_handle));
  
  int interp = 0;

  const int n_odims = 1;
  mwSize o_dims[n_odims];
  int32_t* nthreads;

  // Construct output
  switch (type) {
  case 0: // Interpolation
    interp = p_image->interp_type;
    
    switch (interp) {
    case 0:
      Call(plhs[0] = mxCreateString, ("nearest"));
      break;
    case 1:
      Call(plhs[0] = mxCreateString, ("linear"));
      break;
    case 2:
      Call(plhs[0] = mxCreateString, ("cubic"));
      break;
    case 3:
      Call(plhs[0] = mxCreateString, ("spline"));
      break;
    case 4: 
      Call(plhs[0] = mxCreateString, ("fir"));
      break;
    default:
      Call(plhs[0] = mxCreateString, ("linear"));
      break;
    }
    break;
  case 1:
    o_dims[0] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
         (n_odims, (const mwSize*)o_dims, mxINT32_CLASS, mxREAL));
    nthreads = (int32_t*) mxGetData(plhs[0]);
    *nthreads = (int32_t) p_image->nthreads;
  default:
    nthreads = NULL;
    break;
  }
  return true;
}

extern bool sampled_image_bft_set(mxArray *plhs[], const mxArray* mx_handle,
																	const mxArray* mx_data, const size_t type) {

  char stype[256];
  const int n_odim = 1;
  mwSize o_dims[n_odim];

  if (type==0) {
    Call(mxIsChar, (mx_data));
    Call(mxu_fixed_string, (stype, 256, mx_data, "1st argument"));
  }
  else if (type==1) {
    Call(mxIsScalarInt32, (mx_data));
  }

  SampledImage* p_image;

  bool retval = true;
  int32_t nthreads;

  Call(mxIsPointer, (mx_handle));

  p_image = &(get_object<SampledImage>(mx_handle));

  o_dims[0] = mwSize(1);
  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  switch (type) {
  case 0: {
    if (!strncmp(stype,"nearest",7))
      p_image->interp_type = (SampledImage::nearestNeighbour);
    else if (!strncmp(stype,"linear",6))
      p_image->interp_type = (SampledImage::linear);
    else if (!strncmp(stype,"cubic",5))
      p_image->interp_type = (SampledImage::cubic);
    else if (!strncmp(stype,"spline",6))
      p_image->interp_type = (SampledImage::spline);
    else if (!strncmp(stype,"fir",3))
      p_image->interp_type = (SampledImage::fir);
    break;
  }
  case 1:
    Call(mxIsScalarInt32, (mx_data));
    nthreads = mxGetInt(mx_data);
    if (nthreads > N_MAX_THREADS) {
      Fail(Usage);
    }
    else {
      p_image->nthreads = nthreads;
    }
  default:
    break;
  }
  
  bool* p_blhs = ((bool*) mxGetData(plhs[0]));
  *p_blhs = retval;

  return true;
}

bool sampled_image_compound_slow(mxArray *plhs[], const mxArray* mx_handle, 
																 const mxArray* mx_rf_data, const mxArray* mx_delay, 
																 const mxArray* mx_trans_elem_no,
																 const mxArray* mx_angles) {

  size_t i;
  size_t n_rf_samples, n_recv_channels;

  const float_type* angles; // radians

  const uint32_t* trans_elem_no;

  const float_type *rf_data = NULL;
  const float_type *rf_imag_data = NULL;

  const float_type *delay;

  float_type *image = NULL;
  float_type *imag_image = NULL;

  SampledImage* p_image = NULL;

  const int n_odim = 2;
  mwSize o_dims[n_odim];

  size_t n_angles = 1;

  Call(mxIsPointer, (mx_handle));
  p_image = &(get_object<SampledImage>(mx_handle));

  n_recv_channels = p_image->m_rcv_aperture.n_elements();

  Call(mxIsFloat, (mx_rf_data));
  
  if (n_recv_channels != mxGetDimensions(mx_rf_data)[1]) {
    Fail("Second dimension must match number of receive elements");
  }

  n_rf_samples = mxGetM(mx_rf_data);

	if (n_rf_samples < 4)
		Fail("Interpolation requires at least 4 RF samples");

  rf_data = ((const float_type *) mxGetData(mx_rf_data));

  if (::mxIsComplex(mx_rf_data)) {
    rf_imag_data = ((const float_type*) ::mxGetImagData(mx_rf_data));
  }

  Call(mxIsRealFloat, (mx_delay));
  Call1(mxCheckDimensions, (mx_delay, 2, p_image->m_xmt_aperture.data->m_nemissions, 1),
        "Delays must match number of emission");
  delay = ((const float_type*) mxGetData(mx_delay));

  Call1(mxCheckDimensions, (mx_trans_elem_no, 2, p_image->m_xmt_aperture.data->m_nemissions, 1),
        "Transmit indices must match number of emissions");

  Call(mxIsUint32, (mx_trans_elem_no));
  trans_elem_no = ((const unsigned int *) mxGetData(mx_trans_elem_no));

  for (i=0;i<p_image->m_xmt_aperture.data->m_nemissions;i++) {
    if ((trans_elem_no[i] < 1) ||
				(trans_elem_no[i] > p_image->m_xmt_aperture.n_elements()))
      Fail("Invalid transmit element index");
  }

  Call(mxIsRealFloat, (mx_angles));

  if (2 == mxGetNumberOfDimensions(mx_angles)) {
    n_angles = mxGetDimensions(mx_angles)[1];
  }

  angles = (const float_type*) mxGetData(mx_angles);

#if 0
  for (i=0;i<n_recv_channels;i++) {
    mexPrintf("Delay %lu: %g\n",i,delays[i]);
  }
#endif

  // Construct output, TODO: 3D with nz innermost and nx outermost
  o_dims[0] = mwSize(p_image->im_geom.nz);
  o_dims[1] = mwSize(p_image->im_geom.nx);

  // Allocate Image
  if (rf_imag_data) {
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxCOMPLEX));
  }
  else {
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
  }

  image = ((float_type*) mxGetData(plhs[0]));
  memset(image,0,sizeof(float_type)*o_dims[0]*o_dims[1]);

  if (rf_imag_data) {
    imag_image =  ((float_type*) ::mxGetImagData(plhs[0]));
    memset(imag_image,0,sizeof(float_type)*o_dims[0]*o_dims[1]);
  }

#if 0
  if (1) {
    Note1("Depth samples=%lu", o_dims[0]);
    Note1("First offset=%g", delay[0]);
    Note1("trans_elem_no=%u", trans_elem_no[0]);
    Note1("#RF samples=%lu", n_rf_samples);
  }
#endif

  // n_active_elements, angle
  bool retval = p_image->Beamform(image, rf_data, delay, n_rf_samples,
																	trans_elem_no, n_angles, angles);

  if (rf_imag_data) {
    retval = retval && p_image->Beamform(imag_image, rf_imag_data, delay,
																				 n_rf_samples,
																				 trans_elem_no, n_angles, angles);
  }
  return retval;
}



bool sampled_image_mex(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
  
  char type[256];

  /* Check number of arguments */
        
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    sampled_image_mex_help();
    Fail(Usage);
  }

  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strcmp(type, "sampled_image,ctor")) {
    if ((nrhs != 6) || (nlhs != 1)) {
      sampled_image_mex_help();
      Fail(Usage);
    }
    Call(sampled_image_bft_ctor, (plhs[0], prhs[1], prhs[2], prhs[3], prhs[4], prhs[5]));
  }
  else if (!strcmp(type, "sampled_image,dtor")) {
    if ((nrhs != 2) || (nlhs != 0)) {
      sampled_image_mex_help();
      Fail(Usage);
    }
    Call(sampled_image_bft_dtor, (plhs, prhs[1]));
  }
  else if (!strncmp(type, "sampled_image,get,",18)) {
    if ((nrhs != 2) || (nlhs != 1) || (strlen(type)<24)) {
      sampled_image_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(&type[18],"interp")) {
      Call(sampled_image_bft_get, (plhs, prhs[1],0));
    }
    else if (!strcmp(&type[18],"nthreads")) {
      Call(sampled_image_bft_get, (plhs, prhs[1],1));
    }
  }
  else if (!strncmp(type, "sampled_image,set,",18)) {
    if ((nrhs != 3) || (nlhs != 0) || (strlen(type)<24)) {
      sampled_image_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(&type[18],"interp")) {
      Call(sampled_image_bft_set, (plhs, prhs[1], prhs[2], 0));
    }
    else if (!strcmp(&type[18],"nthreads")) {
      Call(sampled_image_bft_set, (plhs, prhs[1], prhs[2], 1));
    }
  }
  else if (!strcmp(type, "sampled_image,beamform,slow")) {
    if ((nrhs != 6) || (nlhs != 1)) {
      sampled_image_mex_help();
      Fail(Usage);
    }
    Call(sampled_image_compound_slow, (plhs, prhs[1], prhs[2], prhs[3], prhs[4],
																			 prhs[5]));
  }
  else {
    sampled_image_mex_help();
    Fail(Usage);
  }
 
  return true;
}

#if defined(Need_mex_gateway)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (!nlhs && !nrhs) {
    sampled_image_mex_help();
    return;
  }
  if (!sampled_image_mex(nlhs, plhs, nrhs, prhs))
    mexErrMsgTxt("sampled_image_mex()");
}
#endif

#undef Usage

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
