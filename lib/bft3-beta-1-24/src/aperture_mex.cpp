/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: aperture_mex.cpp,v 1.98 2011-07-18 19:11:23 jmh Exp $  *
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-18 19:11:23 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *   - mxGetClassID - return ClassID from mxArray*
 *   - Move validation to Aperture class
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

#include "aperture.h"
#include "aperture_mex.h"

// TEST

/*

#include <sstream>
#include <iostream>
#include <algorithm>

// Screws up
class Base {
public:
	Base(const mxArray* mx) : _mx(mx) {};
	virtual std::ostream& operator>>(std::ostream& os) const = 0;
	// virtual ~Base(); // If virtual, type() must be pure virtual
	mxClassID type() const {return mxGetClassID(_mx);}
protected:
	const mxArray*& _mx;
};

template <int ClassID> class mex_wrap2;

template<>
class mex_wrap2<mxDOUBLE_CLASS> : public Base {
public:
	mex_wrap2(const mxArray* mx) : Base(mx) {};
	std::ostream& operator>>(std::ostream& os) const {
		os << *((const double*) mxGetData(_mx));
		return os;
	}
	~mex_wrap2() {};
};

Base *newBase(const mxArray* mx) {
	const mxClassID ClassID = mxGetClassID(mx);

	switch (ClassID) {
	case mxDOUBLE_CLASS:
		return new mex_wrap2<mxDOUBLE_CLASS>(mx);
		break;
	default:
		return new mex_wrap2<mxDOUBLE_CLASS>(mx);
		break;
	}
}
*/

/*
template<int ClassID>
class mex_wrap {
public:
	mex_wrap(const mxArray* mx) : _mx(mx) {};
	std::ostream& operator>>(std::ostream& os);
	const mxArray*& _mx;
};

template<>
class mex_wrap<mxDOUBLE_CLASS> {
public:
	mex_wrap(const mxArray* mx) : _mx(mx) {};
	std::ostream& operator>>(std::ostream& os) {
		os << *((const double*) mxGetData(_mx));
		return os;
	}
	const mxArray*& _mx;
};
*/

/*
void bla(std::ostream &os, const mxArray* mx) {
	size_t n_dim = mxGetNumberOfDimensions(mx);

	mwSize dims[2];
	dims[0] = mxGetM(mx);
	dims[1] = mxGetN(mx);
	if (n_dim > 1) {
		if (dims[0] > 1) {
			os << "[";
		}
		const double* p_double = (const double*) mxGetData(mx);
		for (size_t i=0;i<std::min(mwSize(2),dims[0]);i++) {
			if (i>0) {
				os << ";";
			}
			os << "[";
			for (size_t j=0 ; j < std::min(dims[1],mwSize(3)) ; j++) {
				if (j>0)
					os << " ";
				os << p_double[j*dims[0]+i];
			}
			if (dims[1]>3)
				os << "...";
			os << "]";
		}
	}
}
*/

/*
inline std::ostream& operator<<(std::ostream& os, const mxArray* mx) {
	char *cbuffer = NULL;
	size_t n = mxGetM(mx) * mxGetN(mx) + 1;

	if (mxIsChar(mx)) {
		cbuffer = (char*) mxCalloc(n, sizeof(char));
		mxGetString(mx, cbuffer, (int)n);
		os << "'" << cbuffer << "'";
	}
	else if (mxIsScalar(mx)) {
		if (mxIsScalarSingle(mx)) {
			os << *((const float*) mxGetData(mx));
		}
		else if (mxIsScalarDouble(mx)) {
			os << *((const double*) mxGetData(mx));
		}
		else if (mxIsScalarInt32(mx)) {
			os << *((const int32_t*) mxGetData(mx));
		}
		else if (mxIsScalarUInt32(mx)) {
			os << *((const uint32_t*) mxGetData(mx));
		}
		else if (mxIsScalarInt64(mx)) {
			os << *((const int64_t*) mxGetData(mx));
		}
		else if (mxIsScalarUInt64(mx)) {
			os << *((const uint64_t*) mxGetData(mx));
		}
	}
	else {
		bla(os,mx);
	}
	if (cbuffer)
		mxFree((void*)cbuffer);
	
	return os;
}
*/

static bool enable_aperture_destructor = false;

extern void aperture_mex_help(void) {
  printf("Usage for aperture_mex:\n\
     handle = function('aperture,ctor,manual', pos);\n\
                pos is an array of shape [#elements, 3], containing\n\
                position of the aperture elements.\n\
\n\
     handle = function('aperture,ctor,field', type, n_elements,...\n\
                pitch, f0);\n\
                type is a string, e.g. 'linear_array', n_elements is the\n\
                number of elements, pitch is the pitch and f0 is the\n\
                center frequency\n\
\n");
}

extern bool aperture_bft_ctor_manual(mxArray*& plhs,
                                     const mxArray* mx_pos_vec) {

  const float_type* pos_vec;

  size_t n_dim;
  size_t n_elements;

  size_t length;
  Aperture* p_aperture;

  Call(mxIsRealFloat, (mx_pos_vec));

  n_elements = mxGetM(mx_pos_vec);

  Call(2 == mxGetNumberOfDimensions,(mx_pos_vec));
  n_dim      = (size_t) mxGetDimensions(mx_pos_vec)[1];
  length      = n_dim*n_elements;

  if (n_dim != (int) 3)
    Fail("Positions must be 3D");

  pos_vec = ((const float_type *) mxGetData(mx_pos_vec));

  p_aperture = new Aperture(pos_vec, (size_t)length);

  plhs = create_handle<Aperture>(p_aperture);

  enable_aperture_destructor = true;
  mexLock(); // Fixed static variables

  return true;
}

extern bool aperture_bft_ctor_field(mxArray*& plhs,
                                    const mxArray* mx_type,
                                    const mxArray* mx_n_elements,
                                    const mxArray* mx_pitch,
                                    const mxArray* mx_f0) {

  const float_type *pitch;
  mwSize n_elements;
  char type[256];
  Aperture* p_aperture;

  Call(mxIsChar, (mx_type));
  Call(mxu_fixed_string, (type, 256, mx_type, "1st argument"));

  if (strcmp(type,"linear_array")) {
    mexPrintf(type);
    Fail("Unsupported Field II aperture");
  }

  Call(mxIsScalarFloat, (mx_pitch));
  Call(mxIsScalarInt32, (mx_n_elements));

  pitch = ((const float_type *) mxGetData(mx_pitch));
  n_elements = mxGetInt(mx_n_elements);
  Aperture::_f0 = *((const float_type*) mxGetData(mx_f0));

  p_aperture = new Aperture((size_t) n_elements, *pitch);
  plhs = create_handle<Aperture>(p_aperture);

  enable_aperture_destructor = true;
  mexLock();
  return true;
}

extern bool aperture_bft_ctor_field_convex(mxArray* &plhs,
																					 const mxArray* mx_type,
																					 const mxArray* mx_n_elements,
																					 const mxArray* mx_width,
																					 const mxArray* mx_kerf,
																					 const mxArray* mx_radius,
																					 const mxArray* mx_f0) {

  const float_type *width;
  const float_type *kerf;
	const float_type *radius;
  mwSize n_elements;
  char type[256];
  Aperture* p_aperture;

  Call(mxIsChar, (mx_type));
  Call(mxu_fixed_string, (type, 256, mx_type, "1st argument"));

  if (strcmp(type,"convex_array")) {
    mexPrintf(type);
    Fail("Unsupported Field II aperture");
  }

  Call(mxIsScalarFloat, (mx_width));
  Call(mxIsScalarFloat, (mx_kerf));
  Call(mxIsScalarFloat, (mx_radius));
  Call(mxIsScalarFloat, (mx_f0));
  Call(mxIsScalarInt32, (mx_n_elements));

  width  = ((const float_type *) mxGetData(mx_width));
  kerf   = ((const float_type *) mxGetData(mx_kerf));
  radius = ((const float_type *) mxGetData(mx_radius));


  n_elements = mxGetInt(mx_n_elements);
  Aperture::_f0 = *((const float_type*) mxGetData(mx_f0));

  p_aperture = new Aperture((size_t) n_elements, *width, *kerf, *radius);
  plhs = create_handle<Aperture>(p_aperture);

  enable_aperture_destructor = true;
  mexLock();
  return true;
}

extern bool aperture_bft_dtor(mxArray *plhs[],
                              const mxArray* mx_handle) {
  
  if (!enable_aperture_destructor) {
    return true;
  }
  Call(mxIsPointer, (mx_handle));
  destroy_object<Aperture>(mx_handle);
  mexUnlock();
  return true;
}

extern bool aperture_bft_get(mxArray *plhs[], const mxArray* mx_handle,
                             const size_t type) {

  const size_t n_odim = 2;
  mwSize o_dims[n_odim];

  bool retval = true;
  size_t i;
  Aperture* p_aperture;
  Aperture* p_aperture_clone;

  float_type* p_flhs = NULL;
  bool*   p_blhs     = NULL;
  size_t* p_ilhs     = NULL;

  int aperturetype = 0;

  Call(mxIsPointer, (mx_handle));
  p_aperture = &(get_object<Aperture>(mx_handle));

  switch (type) {
  case 0:
    // Positions
    o_dims[0] = mwSize(p_aperture->n_elements());
    o_dims[1] = mwSize(3);
  
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));

    for (i=0 ; i < 3*p_aperture->n_elements() ; i++ ) {
      p_flhs[i] = p_aperture->data->m_pos->m_data[0][i];
    }
    break;
  case 1:
    // Focus point (VS)
    o_dims[0] = mwSize(0);
    o_dims[1] = mwSize(0);

    if (p_aperture->data->m_focus) {
      o_dims[0] = mwSize(1);
      o_dims[1] = mwSize(3);
    }

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    
    p_flhs = ((float_type*) mxGetData(plhs[0]));

    for (i=0 ; i < 3 ; i++ ) {
      if (p_aperture->data->m_focus)
        p_flhs[i] = p_aperture->data->m_focus[i];
    }
    break;
  case 2:
    // Center focus
    o_dims[0] = mwSize(p_aperture->data->m_nemissions);
    o_dims[1] = mwSize(3);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));

    for (i=0 ; i < 3*p_aperture->data->m_nemissions ; i++ ) {
      p_flhs[i] = p_aperture->data->m_center_focus[i];
    }
    break;
  case 3:
    // fs
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
      
    p_flhs = ((float_type*) mxGetData(plhs[0]));

    p_flhs[0] = Aperture::_fs;
    break;
  case 4:
    // f0
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    
    p_flhs = ((float_type*) mxGetData(plhs[0]));

    p_flhs[0] = Aperture::_f0;
    break;
  case 7:
    // c
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));
    
    p_flhs = ((float_type*) mxGetData(plhs[0]));

    p_flhs[0] = Aperture::_c;
    break;
  case 5:
    // Receive delays
    if (p_aperture->data->m_delays) {
      o_dims[0] = mwSize(p_aperture->n_elements());
      o_dims[1] = mwSize(1);
    }
    else {
      o_dims[0] = mwSize(0);
      o_dims[1] = mwSize(0);
    }

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));

    if (p_aperture->data->m_delays)
      for (i=0 ; i < p_aperture->n_elements() ; i++ )
        p_flhs[i] = p_aperture->data->m_delays[i];
    break;
  case 6:
    // Id
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxPOINTER_CLASS, mxREAL));

    // Use ptr_t
    p_ilhs = ((size_t*) mxGetData(plhs[0]));
    p_ilhs[0] = p_aperture->data->m_id;
    break;
  case 8:
    // Clone
    p_aperture_clone = p_aperture->Clone();
    plhs[0] = create_handle<Aperture>(p_aperture_clone);
    enable_aperture_destructor = true;
    mexLock();
    break;
  case 9:
		// Version
    Call(plhs[0] = mxCreateString, (PACKAGE_VERSION));
    break;
  case 10:
    // Aperture type
    aperturetype = p_aperture->data->m_type;
    switch (aperturetype) {
    case Aperture::custom:
      Call(plhs[0] = mxCreateString, ("custom"));
      break;
    case Aperture::linear_array:
      Call(plhs[0] = mxCreateString, ("linear_array"));
      break;
    case Aperture::convex_array:
      Call(plhs[0] = mxCreateString, ("convex_array"));
      break;
    default:
      Call(plhs[0] = mxCreateString, ("unknown"));
    };
    break;
  case 11:
    /// Focus delays
    o_dims[0] = mwSize(p_aperture->n_elements());
    o_dims[1] = mwSize(1);
  
    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));

    retval = p_aperture->getFocusDelays(p_flhs);
    if (!retval) {
      retval = true;
      // Return zero values
      /*
        Fail("No focus delays available for un-focused array");
      */
    }
    break;
  case 12:
    // PP-wave on/off
    o_dims[0] = mwSize(1);
    o_dims[1] = mwSize(1);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxLOGICAL_CLASS, mxREAL));

    p_blhs = ((bool*) mxGetData(plhs[0]));
    p_blhs[0] = p_aperture->data->m_ppwave;
    break;
	case 13:
		// Orientation
    // Center focus
    o_dims[0] = mwSize(p_aperture->data->m_nemissions);
    o_dims[1] = mwSize(3);

    Call(plhs[0] = mxCreateNumericArray,
         (n_odim, (const mwSize*)o_dims, mxFLOAT_CLASS, mxREAL));

    p_flhs = ((float_type*) mxGetData(plhs[0]));

    for (i=0 ; i < 3*p_aperture->data->m_nemissions ; i++ ) {
      p_flhs[i] = p_aperture->data->m_euler[i];
    }
    break;
  }
  return retval;
}

extern bool aperture_bft_set(mxArray *plhs[], const mxArray* mx_handle,
                             const mxArray* mx_data, const size_t type) {

  const int n_odim = 1;
  const float_type* data_vec;
  const bool* p_b = NULL;

  mwSize o_dims[n_odim];

	//  mwSize n_elements, length;
  size_t n_elements, length;

  Aperture* p_aperture = NULL;
  bool* p_blhs = NULL;
  char stype[256];

  Call(mxIsPointer, (mx_handle));

  if (type < 7)
    Call(mxIsRealFloat, (mx_data));

  p_aperture = &(get_object<Aperture>(mx_handle));

  o_dims[0] = mwSize(1);

  Call(plhs[0] = mxCreateLogicalArray,
       (n_odim, (const mwSize*)o_dims));

  n_elements = (mwSize) mxGetM(mx_data);
  
  data_vec = ((const float_type *) mxGetData(mx_data));

  switch (type) {
  case 0:
    // Positions
    Call1(mxCheckDimensions, (mx_data, 2, n_elements, 3),
          "Positions must be 3D array");
    length = (size_t)3*n_elements;
    Call1(p_aperture->setPositions, (data_vec, length), "Incompatible size");
    p_aperture->data->m_type = Aperture::custom;
    break;
  case 1:
    // Focus
    if ((mxGetM(mx_data)==0) && (mxGetN(mx_data)==0)) {
      data_vec = NULL;
    }
    else {
      Call1(mxCheckDimensions, (mx_data, 2, 1, 3), "Positions must be 3D");
    }
    p_aperture->setFocus(data_vec);
    break;
  case 2:
    // Center focus
    Call1(mxCheckDimensions, (mx_data, 2, n_elements, 3),
          "Positions must be 3D");
    p_aperture->setCenterFocus(data_vec, n_elements);
    break;
  case 3:
    // fs
    Call1(mxIsScalarFloat, (mx_data), "fs must be scalar float");
    Aperture::_fs = *((const float_type*) mxGetData(mx_data));
    break;
  case 6:
    // c
    Call1(mxIsScalarFloat, (mx_data), "c must be scalar float");
    Aperture::_c = *((const float_type*) mxGetData(mx_data));
    break;
  case 4:
    // f0
    Call1(mxIsScalarFloat, (mx_data), "f0 must be scalar float");
    Aperture::_f0 = *((const float_type*) mxGetData(mx_data));
  case 5:
    // Delays
    if ((mxGetM(mx_data)==0) && (mxGetN(mx_data)==0)) {
      data_vec = NULL;
    }
    else {
      Call1(mxCheckDimensions, (mx_data, 2, p_aperture->n_elements(), 1),
            "Receive delays must match number of elements");
    }
    p_aperture->setDelays(data_vec);
    break;
  case 7:
    // Aperture type
    Call(mxIsChar, (mx_data));
    Call(mxu_fixed_string, (stype, 256, mx_data, "1st argument"));
    if (!strcmp(stype,"custom"))
      p_aperture->data->m_type = Aperture::custom;
    else if (!strcmp(stype,"linear_array"))
      p_aperture->data->m_type = Aperture::linear_array;
    else if (!strcmp(stype,"convex_array"))
      p_aperture->data->m_type = Aperture::convex_array;
    else {
      Fail("Unknown type");
    }
    break;
  case 8:
    // PP-wave on/off
    Call(mxIsLogical, (mx_data));
    p_b = ((const bool*)  mxGetData(mx_data));
    p_aperture->data->m_ppwave = *p_b;
    break;
	case 9:
		// Orientation
    Call1(mxCheckDimensions, (mx_data, 2, n_elements, 3),
          "Orientations must be 3D Euler angles");
    p_aperture->setOrientation(data_vec, n_elements);
    break;
  default:
    break;
  }
  
  p_blhs = ((bool*) mxGetData(plhs[0]));
  *p_blhs = true;
  return true;
}

extern bool aperture_mex(const int nlhs, mxArray *plhs[], const int nrhs,
                  const mxArray *prhs[]) {
  
  char type[256];
  char* attribute = NULL;

  /* Check number of arguments */    
  if (nrhs < 1 || nlhs < 0 || !mxIsChar(prhs[0])) {
    aperture_mex_help();
    Fail(Usage);
  }

  Call(mxu_fixed_string, (type, 256, prhs[0], "1st argument"));

  if (!strcmp(type, "aperture,ctor,manual")) {
    if ((nrhs != 2) || (nlhs != 1)) {
      aperture_mex_help();
      Fail(Usage);
    }
    Call(aperture_bft_ctor_manual, (plhs[0], prhs[1]));
  }
  else if (!strcmp(type, "aperture,ctor,field")) {
    if ((nrhs != 5) || (nlhs != 1)) {
      aperture_mex_help();
      Fail(Usage);
    }
    Call(aperture_bft_ctor_field, (plhs[0], prhs[1], prhs[2],
                                   prhs[3], prhs[4]));
  }
  else if (!strcmp(type, "aperture,ctor,field,convex")) {
    if ((nrhs != 7) || (nlhs != 1)) {
      aperture_mex_help();
      Fail(Usage);
    }
    Call(aperture_bft_ctor_field_convex, (plhs[0], prhs[1], prhs[2],
																					prhs[3], prhs[4], prhs[5], prhs[6]));
  }
  else if (!strcmp(type, "aperture,dtor")) {
    if ((nrhs != 2) || (nlhs != 0)) {
      aperture_mex_help();
      Fail(Usage);
    }
    Call(aperture_bft_dtor, (plhs, prhs[1]));
  }
  else if (!strncmp(type, "aperture,get,",13)) {
    attribute = &type[13];
    if ((nrhs != 2) || (nlhs != 1) || (strlen(type)<14)) {
      aperture_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute, "pos")) {
      Call(aperture_bft_get, (plhs, prhs[1], 0));
    }
    else if (!strcmp(attribute, "focus_delays")) {
      Call(aperture_bft_get, (plhs, prhs[1], 11));
    }
    else if (!strcmp(attribute, "focus")) {
      Call(aperture_bft_get, (plhs, prhs[1], 1));
    }
    else if (!strcmp(attribute, "center_focus")) {
      Call(aperture_bft_get, (plhs, prhs[1], 2));
    }
    else if (!strcmp(attribute, "fs")) {
      Call(aperture_bft_get, (plhs, prhs[1], 3));
    }
    else if (!strcmp(attribute, "f0")) {
      Call(aperture_bft_get, (plhs, prhs[1], 4));
    }
    else if (!strcmp(attribute, "delays")) {
      Call(aperture_bft_get, (plhs, prhs[1], 5));
    }
    else if (!strcmp(attribute, "id")) {
      Call(aperture_bft_get, (plhs, prhs[1], 6));
    }
    else if (!strcmp(attribute, "c")) {
      Call(aperture_bft_get, (plhs, prhs[1], 7));
    }
    else if (!strcmp(attribute, "clone")) {
      Call(aperture_bft_get, (plhs, prhs[1], 8));
    }
    else if (!strcmp(attribute, "version")) {
      Call(aperture_bft_get, (plhs, prhs[1], 9));
    }
    else if (!strcmp(attribute, "type")) {
      Call(aperture_bft_get, (plhs, prhs[1], 10));
    }
    else if (!strcmp(attribute, "ppwave")) {
      Call(aperture_bft_get, (plhs, prhs[1], 12));
    }
    else if (!strcmp(attribute, "orientation")) {
      Call(aperture_bft_get, (plhs, prhs[1], 13));
    }
    else {
      aperture_mex_help();
      Fail(Usage);
    }
  }
  else if (!strncmp(type, "aperture,set,",13)) {
    attribute = &type[13];
    if ((nrhs != 3) || (nlhs != 0) || (strlen(type)<14)) {
      aperture_mex_help();
      Fail(Usage);
    }
    else if (!strcmp(attribute, "pos")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],0));
    }
    else if (!strcmp(attribute, "focus")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],1));
    }
    else if (!strcmp(attribute, "center_focus")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],2));
    }
    else if (!strcmp(attribute, "fs")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],3));
    }
    else if (!strcmp(attribute, "f0")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],4));
    }
    else if (!strcmp(attribute, "delays")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2], 5));
    }
    else if (!strcmp(attribute, "c")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],6));
    }
    else if (!strcmp(attribute, "type")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],7));
    }
    else if (!strcmp(attribute, "ppwave")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],8));
    }
    else if (!strcmp(attribute, "orientation")) {
      Call(aperture_bft_set, (plhs, prhs[1], prhs[2],9));
    }
    else {
      aperture_mex_help();
      Fail(Usage);
    }
  }
  else {
    aperture_mex_help();
    Fail(Usage);
  }
 
  return true;
}

#if defined(Need_mex_gateway)

extern void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  if (!nlhs && !nrhs) {
    aperture_mex_help();
    return;
  }
  if (!aperture_mex(nlhs, plhs, nrhs, prhs)) {

    // Experiment to write out how the method was called 
    std::ostringstream m;
    m << "aperture_mex(";
		int i;
    for (i=0; i < (nrhs-1) ; i++) {
			m << prhs[i] << ", ";
    }
    m << prhs[i] << ");";
    std::string entryname = m.str();
    mexErrMsgTxt(entryname.c_str());

/*
    mexErrMsgTxt("aperture_mex()");
*/
  }
}
#endif

#undef Usage

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
