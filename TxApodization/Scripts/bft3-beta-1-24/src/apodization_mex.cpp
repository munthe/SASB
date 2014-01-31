/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: apodization_mex.cpp,v 1.69 2011-07-11 14:34:58 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-11 14:34:58 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Add apodizations to an emission and a subaperture, default is low
 *    is first element high is last element
 *  - Verify why GC operates when constructing an aperture here.
 *    (scenario is, object from aperture DLL is referenced here,
 *    reference counter is incremented and when deleted using aperture
 *    DLL it is not freed. It is only freed when apodization DLL is freed
 * 
 ****************************************************************************/

#include "mex_utility.h"
#include "mexarg.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#define Usage   "Usage error. see above"

#include "aperture.h"
#include "apodization.h"
#include "apodization_mex.h"

static bool enable_apodization_destructor = false;

void apodization_mex_help(void) {
  printf("Information\n");
}

bool apodization_bft_ctor_fixed(mxArray*& plhs,
																			 const mxArray* mx_handle,
																			 const mxArray* mx_pos_ref,
																			 const mxArray* mx_distances,
																			 const mxArray* mx_values) {

  const float_type  *pos_ref, *distances, *values;
  mwSize n_elements, n_distances;
  Apodization* p_apodization;
  Aperture*    p_aperture;

  Call(mxIsPointer, (mx_handle));
  Call(mxIsRealFloat, (mx_pos_ref));
  Call(mxIsRealFloat, (mx_distances));
  Call(mxIsRealFloat, (mx_values));

  Call1(mxCheckDimensions, (mx_pos_ref, 2, 1, 3),"Positions must be 3D");

  p_aperture = &(get_object<Aperture>(mx_handle));
  n_elements = p_aperture->data->m_npos;
  n_distances = mxGetDimensions(mx_distances)[1];

  Call1(mxCheckDimensions, (mx_distances, 2, 1, n_distances),\
        "Distances must be a vector");
  Call1(mxCheckDimensions, (mx_values, 2, n_elements, n_distances),\
        "Apodization values must have dimension (n_elements x n_distance)");

  pos_ref   = ((const float_type *) mxGetData(mx_pos_ref));
  distances = ((const float_type *) mxGetData(mx_distances));
  values    = ((const float_type *) mxGetData(mx_values));

  // TODO: Move validation into C++ class
  p_apodization = new Apodization(*p_aperture, pos_ref, distances,
                                  n_distances, values, n_elements);

  // Multiple objects using garbage collector
  plhs = create_handle<Apodization>(p_apodization);

  enable_apodization_destructor = true;
  mexLock();
  return true;
}

bool apodization_bft_dtor(mxArray *plhs[], const mxArray* mx_handle) {

  if (!enable_apodization_destructor)
    return true;
  
  Call(mxIsPointer, (mx_handle));
  
  destroy_object<Apodization>(mx_handle);
  mexUnlock();
  return true;
}

bool apodization_bft_get(mxArray *plhs[],
                                const mxArray* mx_handle, const size_t type) {

  const int n_odim = 2;
  mwSize o_dims[n_odim];
  size_t n_outputsize = 0;
  Apodization* p_apodization = NULL;
  Apodization* p_apodization_clone = NULL;
  float_type* p_data = NULL;
  Aperture* p_aperture = NULL;
  size_t* p_ilhs = NULL;
  uint32_t* p_uint32lhs = NULL;
  bool*   p_blhs = NULL;
  float_type*   p_flhs = NULL;
  int apowindow = 0;
  size_t i;

  Call(mxIsPointer, (mx_handle));

  p_apodization = &(get_object<Apodization>(mx_handle));

  // Construct output
  switch (type) {
  case 0:
    // Aperture
    // Memory leak taken care of by garbage collector. Solved when we
    // use a single library
    
    p_aperture = new Aperture();
    *p_aperture = p_apodization->data->m_aperture;
    plhs[0] = create_handle<Aperture>(p_aperture);    
    break;
  case 1:
    // Reference point
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(3);
    n_outputsize = 3;
    p_data = p_apodization->data->m_ref;
    break;
  case 2:
    // Distances
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(p_apodization->data->m_ndistances);
    n_outputsize = p_apodization->data->m_ndistances;
    p_data = p_apodization->data->m_distances;
    break;
  case 3:
    // Values
    o_dims[0] = mwSize(p_apodization->data->m_nelements);
    o_dims[1] = mwSize(p_apodization->data->m_ndistances);
    n_outputsize = p_apodization->data->m_nelements * \
      p_apodization->data->m_ndistances;
    p_data = p_apodization->data->m_values;
    break;
  case 5:
    // Id
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxPOINTER_CLASS, mxREAL));
    // Use ptr_t
    p_ilhs = ((size_t*) mxGetData(plhs[0]));
    p_ilhs[0] = p_apodization->data->m_id;
    break;
  case 6:
    // Dynamic on/off
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxLOGICAL_CLASS, mxREAL));
    p_blhs = ((bool*) mxGetData(plhs[0]));
    p_blhs[0] = p_apodization->data->m_dynamic;
    break;
  case 7:
    // Apodization window
    apowindow = p_apodization->data->m_window;
    switch (apowindow) {
    case Apodization::HAMMING:
      Call(plhs[0] = mxCreateString, ("Hamming"));
      break;
    case Apodization::HANN:
      Call(plhs[0] = mxCreateString, ("Hann"));
      break;
    case Apodization::BLACKMAN:
      Call(plhs[0] = mxCreateString, ("Blackman"));
      break;
    case Apodization::BARTLETT:
      Call(plhs[0] = mxCreateString, ("Bartlett"));
      break;
    case Apodization::RECTWIN:
      Call(plhs[0] = mxCreateString, ("Rectwin"));
      break;
    case Apodization::TUKEY:
      Call(plhs[0] = mxCreateString, ("Tukey"));
      break;
    case Apodization::GAUSSIAN:
      Call(plhs[0] = mxCreateString, ("Gaussian"));
      break;
    default:
      Call(plhs[0] = mxCreateString, ("Unknown"));
    };
    break;
  case 8:
    // Fixed on/off
    // TODO: Support other apertures than linear arrays
    /*
			if (p_apodization->data->m_aperture.data->m_type != Aperture::linear_array) {
      Fail("Fixed apodization only supported for linear array apertures");
			}
    */
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxLOGICAL_CLASS, mxREAL));
    p_blhs = ((bool*) mxGetData(plhs[0]));
    p_blhs[0] = p_apodization->data->m_fixed;
    break;
  case 9:
    // F number
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    
    p_flhs = ((float_type*) mxGetData(plhs[0]));
    p_flhs[0] = p_apodization->data->m_f;
    break;
  case 10:
    // Clone
    p_apodization_clone = p_apodization->Clone2(); // Old behavior
    // new Apodization(); = .Clone();
    plhs[0] = create_handle<Apodization>(p_apodization_clone);
    enable_apodization_destructor = true;
    mexLock();
    break;
  case 11:
    // Window parameter
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    
    p_flhs = ((float_type*) mxGetData(plhs[0]));
    p_flhs[0] = p_apodization->data->m_window_param;
    break;
  case 12:
    // Manual on/off
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
				 (n_odim, (const mwSize*)o_dims, mxLOGICAL_CLASS, mxREAL));
    p_blhs = ((bool*) mxGetData(plhs[0]));
    p_blhs[0] = p_apodization->data->m_manual;
    break;
  case 13:
    // Active elements
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxUINT32_CLASS, mxREAL));
    p_uint32lhs = ((uint32_t*) mxGetData(plhs[0]));
    p_uint32lhs[0] = (uint32_t) p_apodization->data->m_nactive_elements;
    break;
  case 14:
    // Orientation
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(3);
    n_outputsize = 3;
    p_data = p_apodization->data->m_euler;
    break;
  default:
    n_outputsize = 0;
    o_dims[0] = mwSize(0);
    o_dims[1] = mwSize(0);
  }

  if (((type < 4) && (type > 0)) || (type==14)) {
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));
  
    if (p_data) {
      for (i=0;i<n_outputsize;i++) {
        p_flhs[i] = p_data[i];
      }
    }
  }
  return true;
}

bool apodization_bft_set(mxArray *plhs[], const mxArray* mx_handle,
                                const mxArray* mx_data, const size_t type) {

  const int n_odim = 1;
  float_type* p_data = NULL;
  const bool* p_b = NULL;
  const uint32_t* p_uint32 = NULL;
  mwSize o_dims[n_odim];
  Apodization* p_apodization = NULL;
  Aperture* p_aperture = NULL;
  mwSize n_row, n_col, length;
  size_t i;

  char stype[256];

  Call(mxIsPointer, (mx_handle));
  p_apodization = &(get_object<Apodization>(mx_handle));

  if ((type < 4) && (type > 0)) {
    Call(mxIsRealFloat, (mx_data));
    n_row = mxGetM(mx_data);
    n_col = mxGetDimensions(mx_data)[1];
    length = n_row*n_col;
    p_data = ((float_type *) mxGetData(mx_data));
  }

  // Cheating slow way copy everything everytime
  float_type* pos_ref = p_apodization->data->m_ref;
  float_type* distances = p_apodization->data->m_distances;
  float_type* values = p_apodization->data->m_values;
  
  size_t n_distances = p_apodization->data->m_ndistances;
  size_t n_elements = p_apodization->data->m_nelements;

  // TODO: Validate new inputs
  switch (type) {
  case 0:
    // Warning all other references also change aperture, TODO: clone
    Call(mxIsPointer, (mx_data));
    p_aperture = &(get_object<Aperture>(mx_data));
    if (!p_apodization->setAperture(p_aperture))
      Fail("Incompatible aperture");
    break;
  case 1:
    // Reference point
    length = 3;
    pos_ref = p_data;
    // Aperture, for now return reference
    //    Aperture* p_aperture = new Aperture();
    //    p_aperture = p_apodization->data->m_aperture;
    
    length = 3;
    pos_ref = p_data;
    Call1(mxCheckDimensions, (mx_data, 2, 1, 3),"Positions must be 3D");
    break;
  case 2:
    // Distances
    length = n_distances;
    distances = p_data;
    Call1(mxCheckDimensions,
					(mx_data, 2, 1, n_distances),
					"Dimension of distances must match dimension of apodization");
    break;
  case 3:
    // Values
    length = n_distances*n_elements;
    values = p_data;
    Call1(mxCheckDimensions,
					(mx_data, 2, n_elements, n_distances),
					"Must have dimension (# elements x # distances)");
    break;
  case 5:
    Call(mxIsLogical, (mx_data));
    p_b = ((const bool*)  mxGetData(mx_data));
    p_apodization->data->m_dynamic = *p_b;
    break;
  case 6:
    Call(mxIsChar, (mx_data));
    Call(mxu_fixed_string, (stype, 256, mx_data, "1st argument"));
    if (!strcmp(stype,"Hamming"))
      p_apodization->data->m_window = (Apodization::HAMMING);
    else if (!strcmp(stype,"Hann"))
      p_apodization->data->m_window = (Apodization::HANN);
    else if (!strcmp(stype,"Blackman"))
      p_apodization->data->m_window = (Apodization::BLACKMAN);
    else if (!strcmp(stype,"Bartlett"))
      p_apodization->data->m_window = (Apodization::BARTLETT);
    else if (!strcmp(stype,"Rectwin"))
      p_apodization->data->m_window = (Apodization::RECTWIN);
    else if (!strcmp(stype,"Tukey"))
      p_apodization->data->m_window = (Apodization::TUKEY);
    else if (!strcmp(stype,"Gaussian"))
      p_apodization->data->m_window = (Apodization::GAUSSIAN);
    else {
      Fail("Unknown window function, valid window functions are:\nHamming, Hann, Blackman, Bartlett, Rectwin, Tukey, and Gaussian.");
    }
    break;
  case 7:
    // Fixed
    Call(mxIsLogical, (mx_data));
    p_b = ((const bool*)  mxGetData(mx_data));
		/*
    if ((p_apodization->data->m_aperture.data->m_type == Aperture::custom) &&
				(*p_b)){
      Fail("Fixed apodization only supported for linear and convex array apertures");
    }
		*/
    p_apodization->data->m_fixed = *p_b;
    break;
  case 8:
    // F-number
    Call1(mxIsScalarFloat, (mx_data),"F# must be scalar float");
    p_data = ((float_type*)  mxGetData(mx_data));
    p_apodization->data->m_f = *p_data;
    break;
  case 9:
    // Window parameter
    Call1(mxIsScalarFloat, (mx_data), "Window parameter must be scalar float");
    p_data = ((float_type*)  mxGetData(mx_data));
    p_apodization->data->m_window_param = *p_data;
    break;
  case 10:
    // Manual
    Call(mxIsLogical, (mx_data));
    p_b = ((const bool*)  mxGetData(mx_data));
    p_apodization->data->m_manual = *p_b;
    break;
  case 11:
    // n_active_elements
    Call(mxIsScalarUInt32,(mx_data));
    p_uint32 = ((const uint32_t*) mxGetData(mx_data));
    p_apodization->data->m_nactive_elements =
      std::min((size_t)(*p_uint32),p_apodization->data->m_nelements);
    break;
  case 12:
    // Orientation
    Call(mxIsFloat,(mx_data));
    Call1(mxCheckDimensions, (mx_data, 2, 1, 3),"Orientation must be 3D Euler angles");
    p_data = ((float_type*)  mxGetData(mx_data));
    for (i=0;i<3;i++) {
      p_apodization->data->m_euler[i] = p_data[i];
    }
    /*
      p_apodization->data->m_euler[0] = M_PI/4.0;
      p_apodization->data->m_euler[1] = float_type(0.0);
    */
		/*
    WarnM("Only few orientations are supported for bft3_sampled_image class:\n\tFor this class orientation is assumed to be [pi/2 0.0 %f]\n", p_apodization->data->m_euler[2]);
		*/
    break;
  }

  o_dims[0] = mwSize(1);
  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  if (type < 4) {
    // Update apodization, TODO: Move to Apodization class
    Apodization updated =
      Apodization(p_apodization->data->m_aperture, pos_ref,
									distances, n_distances, values, n_elements);
    updated.data->m_f = p_apodization->data->m_f;
    updated.data->m_fixed = p_apodization->data->m_fixed;
    updated.data->m_dynamic = p_apodization->data->m_dynamic;
    updated.data->m_manual = p_apodization->data->m_manual;
    updated.data->m_window = p_apodization->data->m_window;
    *p_apodization = updated;
  }

  bool* p_blhs = ((bool*) mxGetData(plhs[0]));
  *p_blhs = true;

  return true;

}

bool apodization_mex(const int nlhs, mxArray *plhs[],
										 const int nrhs, const mxArray *prhs[]) {
  
  char type[256];

  char* attribute = NULL;

  /* Check number of arguments */
        
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    apodization_mex_help();
    Fail(Usage);
  }

  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strcmp(type, "apodization,ctor,manual")) {
    if ((nrhs != 5) || (nlhs != 1)) {
      apodization_mex_help();
      Fail(Usage)
        }
    Call(apodization_bft_ctor_fixed, (plhs[0], prhs[1], prhs[2],
																			prhs[3], prhs[4]));
  }
  else if (!strncmp(type, "apodization,get,",16)) {
    attribute = &type[16];
    if ((nrhs != 2) || (nlhs != 1) || (strlen(type)<17)) {
      apodization_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute,"aperture")) {
      Call(apodization_bft_get, (plhs, prhs[1],0));
    }
    else if (!strcmp(attribute,"ref")) {
      Call(apodization_bft_get, (plhs, prhs[1],1));
    }
    else if (!strcmp(attribute,"distances")) {
      Call(apodization_bft_get, (plhs, prhs[1],2));
    }
    else if (!strcmp(attribute,"values")) {
      Call(apodization_bft_get, (plhs, prhs[1],3));
    }
    else if (!strcmp(attribute,"id")) {
      Call(apodization_bft_get, (plhs, prhs[1],5));
    }
    else if (!strcmp(attribute,"dynamic")) {
      Call(apodization_bft_get, (plhs, prhs[1],6));
    }
    else if (!strcmp(attribute,"window_parameter")) {
      Call(apodization_bft_get, (plhs, prhs[1],11));
    }
    else if (!strncmp(attribute,"window",6)) {
      Call(apodization_bft_get, (plhs, prhs[1],7));
    }
    else if (!strncmp(attribute,"fixed",5)) {
      Call(apodization_bft_get, (plhs, prhs[1],8));
    }
    else if (!strcmp(attribute,"f")) {
      Call(apodization_bft_get, (plhs, prhs[1],9));
    }
    else if (!strcmp(attribute,"clone")) {
      Call(apodization_bft_get, (plhs, prhs[1],10));
    }
    else if (!strcmp(attribute,"parametric")) {
      Call(apodization_bft_get, (plhs, prhs[1],12));
    }
    else if (!strncmp(attribute,"n_active_elements",17)) {
      Call(apodization_bft_get, (plhs, prhs[1],13));
    }
    else if (!strncmp(attribute,"orientation",11)) {
      Call(apodization_bft_get, (plhs, prhs[1],14));
    }
    /*
			else if (!strncmp(attribute, "focus_delays",12)) {
      Call(apodization_bft_get, (plhs, prhs[1], 16));
			}
    */
    else {
      apodization_mex_help();
      Fail(Usage);     
    }
  }
  else if (!strncmp(type, "apodization,set,",16)) {
    attribute = &type[16];
    if ((nrhs != 3) || (nlhs != 0) || (strlen(type)<17)) {
      apodization_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute,"aperture")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 0));
    }
    else if (!strcmp(attribute,"ref")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 1));
    }
    else if (!strcmp(attribute,"distances")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 2));
    }
    else if (!strcmp(attribute,"values")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 3));
    }
    else if (!strcmp(attribute,"dynamic")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 5));
    }
    else if (!strcmp(attribute,"window_parameter")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 9));
    }
    else if (!strcmp(attribute,"window")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 6));
    }
    else if (!strncmp(attribute,"fixed",5)) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 7));
    }
    else if (!strcmp(attribute,"f")) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 8));
    }
    else if (!strncmp(attribute,"parametric",6)) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 10));
    }
    else if (!strncmp(attribute,"n_active_elements",17)) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 11));
    }
    else if (!strncmp(attribute,"orientation",11)) {
      Call(apodization_bft_set, (plhs, prhs[1], prhs[2], 12));
    }
    else {
      apodization_mex_help();
      Fail(Usage);
    }
  }
  else if (!strcmp(type, "apodization,dtor")) {
    if ((nrhs != 2) || (nlhs != 0)) {
      apodization_mex_help();
      Fail(Usage);
    }
    Call(apodization_bft_dtor, (plhs, prhs[1]));
  }
  else {
    apodization_mex_help();
    Fail(Usage);
  }
 
  return true;
}

#if defined(Need_mex_gateway)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (!nlhs && !nrhs) {
    apodization_mex_help();
    return;
  }
  if (!apodization_mex(nlhs, plhs, nrhs, prhs))
    mexErrMsgTxt("apodization_mex()");
}
#endif

#undef Usage

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
