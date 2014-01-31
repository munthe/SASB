/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: image_mex.cpp,v 1.67 2011-07-27 21:57:42 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-27 21:57:42 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Use common headers
 *  - Assign apodization in line ctor
 *  - Check third dimension of data
 ****************************************************************************/

#include "mex_utility.h"
#include "mexarg.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <cmath>

#include <vector>
using std::vector;

#include <algorithm> // generate_n

#define Usage   "Usage error. see above"

#include "common.h"
#include "aperture.h"
#include "image.h"
#include "line.h"

#include "image_mex.h"

static bool enable_image_destructor = false;

extern void image_mex_help(void)
{
  printf("Information\n");
}

extern bool image_bft_ctor_manual(mxArray*& plhs,
                                  const mxArray* mx_handle_xmt,
                                  const mxArray* mx_handle_rcv,
                                  const mxArray* mx_handle_xmt_apo,
                                  const mxArray* mx_handle_rcv_apo,
                                  const mxArray* mx_handle_lines) {

  mwSize n_lines, n_rcv_apodizations, n_xmt_apodizations;

  vector<Line> lines;
  vector<Apodization> xmt_apodizations;
  vector<Apodization> rcv_apodizations;

  Call(mxIsPointer, (mx_handle_xmt));
  Call(mxIsPointer, (mx_handle_rcv));

  Aperture xmt_aperture = get_object<Aperture>(mx_handle_xmt);
  Aperture rcv_aperture = get_object<Aperture>(mx_handle_rcv);

  Call(mxIsPointerArray, (mx_handle_lines));
  Call(2 == mxGetNumberOfDimensions, (mx_handle_lines)); 
  n_lines  = mxGetDimensions(mx_handle_lines)[1];
  Call1(mxCheckDimensions, (mx_handle_lines, 2, 1, n_lines),\
        "Line handles must be a vector");

  Call(mxIsPointerArray, (mx_handle_xmt_apo));
  Call(2 == mxGetNumberOfDimensions, (mx_handle_xmt_apo));
  n_xmt_apodizations = mxGetDimensions(mx_handle_xmt_apo)[1];

  Call(mxIsPointerArray, (mx_handle_rcv_apo));
  Call(2 == mxGetNumberOfDimensions, (mx_handle_rcv_apo));
  n_rcv_apodizations = mxGetDimensions(mx_handle_rcv_apo)[1];
  Call1(mxCheckDimensions, (mx_handle_lines, 2, 1, n_lines),\
        "Number of rcv_apodizations must equal number of lines");

  if ((n_lines != n_xmt_apodizations) ||		\
      (n_lines != n_rcv_apodizations)) {
    Fail("Number of apodizations must equal number of lines");
  }

  // Get line objects
  generate_n(std::back_inserter(lines), n_lines,
             object_generator<Line>(mx_handle_lines));

  // Get transmit apodization objects
  generate_n(std::back_inserter(xmt_apodizations), n_xmt_apodizations,
             object_generator<Apodization>(mx_handle_xmt_apo));

  // Get receive apodization objects
  generate_n(std::back_inserter(rcv_apodizations), n_rcv_apodizations,
             object_generator<Apodization>(mx_handle_rcv_apo));

  // Construct image
  Image* p_image = new Image(xmt_aperture, rcv_aperture,
                             xmt_apodizations,
                             rcv_apodizations, lines);

  plhs = create_handle<Image>(p_image);

  enable_image_destructor = true;
  mexLock();

  return true;
}


extern bool image_bft_dtor(mxArray *plhs[], const mxArray* mx_handle) {

  if (!enable_image_destructor)
    return true;
  
  Call(mxIsPointer, (mx_handle));
  destroy_object<Image>(mx_handle);
  mexUnlock();
  return true;
}

extern bool image_bft_set(mxArray *plhs[], const mxArray* mx_handle,
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
  Image* p_image = NULL;

  bool retval = true;
  int32_t nthreads;

  Call(mxIsPointer, (mx_handle));
  p_image = &(get_object<Image>(mx_handle));

  o_dims[0] = mwSize(1);
  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  switch (type) {
  case 0: {
    if (!strcmp(stype,"nearest"))
      p_image->data->m_interp_type = (Image::nearestNeighbour);
    else if (!strcmp(stype,"linear"))
      p_image->data->m_interp_type = (Image::linear);
    else if (!strcmp(stype,"cubic"))
      p_image->data->m_interp_type = (Image::cubic);
    else if (!strncmp(stype,"spline",6))
      p_image->data->m_interp_type = (Image::spline);
    else if (!strcmp(stype,"fir"))
      p_image->data->m_interp_type = (Image::fir);
    break;
  }
  case 1:
    Call(mxIsScalarInt32, (mx_data));
    nthreads = mxGetInt(mx_data);
    if (nthreads > N_MAX_THREADS) {
      Fail("Maximum number of threads exceeded");
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

extern bool image_bft_get(mxArray *plhs[],
                          const mxArray* mx_handle, const size_t type) {

  Image* p_image;

  Call(mxIsPointer, (mx_handle));
  
  p_image = &(get_object<Image>(mx_handle));
  
  int interp = 0;

  const int n_odims = 1;
  mwSize o_dims[n_odims];
  int32_t* nthreads = NULL;

  // Construct output
  switch (type) {
  case 0: // Interpolation
    interp = p_image->data->m_interp_type;
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
    break;
  default:
    break;
  }
  return true;
}

bool image_beamform_slow(mxArray *plhs[], const mxArray* mx_handle, 
                         const mxArray* mx_rf_data, const mxArray* mx_delay, 
                         const mxArray* mx_trans_elem_no) {

  size_t i;
  size_t n_rf_samples, n_rcv_channels;
  const uint32_t* trans_elem_no;

  const float_type *rf_data = NULL;
  const float_type *rf_imag_data = NULL;

  const float_type *delay;
  float_type *image = NULL;
  float_type *imag_image = NULL;

  Image* p_image = NULL;

  const int n_odim = 2;
  mwSize o_dims[n_odim];

  Call(mxIsPointer, (mx_handle));
  p_image = &(get_object<Image>(mx_handle));

  n_rcv_channels = p_image->data->m_rcv_aperture.n_elements();

  Call(mxIsFloat, (mx_rf_data));
  
  if (n_rcv_channels != mxGetDimensions(mx_rf_data)[1]) {
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
  Call1(mxCheckDimensions, (mx_delay, 2, p_image->data->m_xmt_aperture.data->m_nemissions, 1),
        "Delays must match number of emissions");
  delay = ((const float_type*) mxGetData(mx_delay));

  Call1(mxCheckDimensions, (mx_trans_elem_no, 2,
														p_image->data->m_xmt_aperture.data->m_nemissions, 1),
        "Transmit indices must match number of emissions");

  Call(mxIsUint32, (mx_trans_elem_no));
  trans_elem_no = ((const unsigned int *) mxGetData(mx_trans_elem_no));

  for (i=0;i<p_image->data->m_xmt_aperture.data->m_nemissions;i++) {
    if ((trans_elem_no[i] < 1) ||
				(trans_elem_no[i] > p_image->data->m_xmt_aperture.n_elements()))
      Fail("Invalid transmit element index");
  }

  // Construct output
  o_dims[0] = mwSize(floor(p_image->data->m_lines[0].data->m_length / p_image->data->m_lines[0].data->m_dr));
  o_dims[1] = mwSize(p_image->data->m_lines.size());

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

  bool retval = p_image->Beamform(image, rf_data, delay, n_rf_samples,
																	trans_elem_no);

  if (rf_imag_data) {
    retval = retval && p_image->Beamform(imag_image, rf_imag_data, delay,
																				 n_rf_samples, trans_elem_no);
  }
  return retval;
}

bool image_mex(const int nlhs, mxArray *plhs[], const int nrhs,
							 const mxArray *prhs[]) {
  
  char type[256];

  char* attribute = NULL;

  /* Check number of arguments */
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    image_mex_help();
    Fail(Usage);
  }

  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strcmp(type, "image,ctor,manual")) {
    if ((nrhs != 6) || (nlhs != 1)) {
      image_mex_help();
      Fail(Usage);
    }
    Call(image_bft_ctor_manual, (plhs[0], prhs[1], prhs[2],
																 prhs[3], prhs[4], prhs[5]));
  }
  else if (!strcmp(type, "image,dtor")) {
    if ((nrhs != 2) || (nlhs != 0)) {
      image_mex_help();
      Fail(Usage);
    }
    Call(image_bft_dtor, (plhs, prhs[1]));
  }
  else if (!strncmp(type, "image,get,",10)) {
    attribute = &type[10];
    if ((nrhs != 2) || (nlhs != 1) || (strlen(type)<16)) {
      image_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute,"interp")) {
      Call(image_bft_get, (plhs, prhs[1],0));
    }
    else if (!strcmp(attribute,"nthreads")) {
      Call(image_bft_get, (plhs, prhs[1],1));
    }
  }
  else if (!strncmp(type, "image,set,",10)) {
    attribute = &type[10];
    if ((nrhs != 3) || (nlhs != 0) || (strlen(type)<16)) {
      image_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute,"interp")) {
      Call(image_bft_set, (plhs, prhs[1], prhs[2], 0));
    }
    else if (!strcmp(attribute,"nthreads")) {
      Call(image_bft_set, (plhs, prhs[1], prhs[2], 1));
    }
  }
  // Functions
  else if (!strcmp(type, "image,beamform,slow")) {
    if ((nrhs != 5) || (nlhs != 1)) {
      image_mex_help();
      Fail(Usage);
    }
    Call(image_beamform_slow, (plhs, prhs[1], prhs[2], prhs[3], prhs[4]));
  }
  else {
    image_mex_help();
    Fail(Usage);
  }
 
  return true;
}

#if defined(Need_mex_gateway)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (!nlhs && !nrhs) {
    image_mex_help();
  }
  else if (!image_mex(nlhs, plhs, nrhs, prhs))
    mexErrMsgTxt("image_mex()");
}
#endif

#undef Usage

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
