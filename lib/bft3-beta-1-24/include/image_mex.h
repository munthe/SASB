/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: image_mex.h,v 1.5 2011-04-05 20:10:27 jmh Exp $                       *
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-05 20:10:27 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * $Log: image_mex.h,v $
 * Revision 1.5  2011-04-05 20:10:27  jmh
 * *** empty log message ***
 *
 * Revision 1.4  2011/03/26 14:46:08  jmh
 * *** empty log message ***
 *
 * Revision 1.3  2010/05/05 09:59:26  jmh
 * Typo bug
 *
 * Revision 1.2  2010/05/05 09:50:19  jmh
 * *** empty log message ***
 *
 * Revision 1.1  2010/05/04 11:23:39  jmh
 * *** empty log message ***
 *
 *****************************************************************************/

/*****************************************************************************
 * TODO:
 *  - Move headers to here from image_mex.cpp
 *
 *****************************************************************************/

#ifndef IMAGE_MEX_H
#define IMAGE_MEX_H

//#define Need_image_mex_gateway 1

/** 
 * Print out help
 * 
 */
extern void image_mex_help(void);

/** 
 * Constructor: Manual
 * 
 * @param plhs 
 * @param mx_handle_Th_t 
 * @param mx_handle_Th_r 
 * @param mx_handle_Ah_t 
 * @param mx_handle_Ah_r 
 * @param mx_handle_lines 
 * 
 * @return 
 */
extern bool image_bft_ctor_manual(mxArray*& plhs,
                                  const mxArray* mx_handle_Th_t,
                                  const mxArray* mx_handle_Th_r,
                                  const mxArray* mx_handle_Ah_t,
                                  const mxArray* mx_handle_Ah_r,
                                  const mxArray* mx_handle_lines);
/** 
 * Destructor
 * 
 * @param plhs 
 * @param mx_handle 
 * 
 * @return 
 */
extern bool image_bft_dtor(mxArray *plhs[], const mxArray* mx_handle);

/** 
 * Set Image property
 * 
 * @param plhs 
 * @param mx_handle 
 * @param mx_data 
 * @param type 
 * 
 * @return 
 */
extern bool image_bft_set(mxArray *plhs[], const mxArray* mx_handle,
                          const mxArray* mx_data, const size_t type);

/** 
 * Get Image property
 * 
 * @param plhs 
 * @param mx_handle 
 * @param type 
 * 
 * @return 
 */
extern bool image_bft_get(mxArray *plhs[],
                          const mxArray* mx_handle, const size_t type);

/** 
 * Slow Beam Formation
 * 
 * @param plhs 
 * @param mx_handle 
 * @param mx_rf_data 
 * @param mx_delay 
 * @param mx_trans_elem_no 
 * 
 * @return 
 */
bool image_beamform_slow(mxArray *plhs[], const mxArray* mx_handle, 
                         const mxArray* mx_rf_data, const mxArray* mx_delay, 
                         const mxArray* mx_trans_elem_no);

/** 
 * Function selector
 * 
 * @param nlhs 
 * @param plhs 
 * @param nrhs 
 * @param prhs 
 * 
 * @return 
 */
bool image_mex(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]);

#if defined(Need_image_mex_gateway)

/** 
 * Gateway routine
 * 
 * @param nlhs 
 * @param plhs 
 * @param nrhs 
 * @param prhs 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#endif

#endif
