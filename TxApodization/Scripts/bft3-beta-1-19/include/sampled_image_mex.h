/******************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: sampled_image_mex.h,v 1.3 2011-04-27 20:35:27 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-27 20:35:27 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 ******************************************************************************/

/****************************************************************************
 * TODO:
 *
 ****************************************************************************/

#ifndef SAMPLED_IMAGE_MEX_H
#define SAMPLED_IMAGE_MEX_H

//#define Need_sampled_image_mex_gateway 1

extern void sampled_image_mex_help(void);

// Classes only supported in mex-interface in later versions
#if MX_API_VER > 0x07030000

extern bool sampled_image_bft_ctor(mxArray*& plhs,
				   const mxArray* mx_handle_Th_t,
				   const mxArray* mx_handle_Th_r,
				   const mxArray* mx_handle_Ah_t,
				   const mxArray* mx_handle_Ah_r,
				   const mxArray* mx_im_geom);
#endif

extern bool sampled_image_bft_dtor(mxArray *plhs[], const mxArray* mx_handle);

extern bool sampled_image_bft_get(mxArray *plhs[],
				  const mxArray* mx_handle, const size_t type);
extern bool sampled_image_bft_set(mxArray *plhs[], const mxArray* mx_handle,
				  const mxArray* mx_data, const size_t type);

bool sampled_image_compound_slow(mxArray *plhs[], const mxArray* mx_handle, 
				 const mxArray* mx_rf_data, const mxArray* mx_delay, 
				 const mxArray* mx_trans_elem_no,
				 const mxArray* mx_angles);

bool sampled_image_mex(const int nlhs, mxArray *plhs[], const int nrhs,
		       const mxArray *prhs[]);

#endif

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
