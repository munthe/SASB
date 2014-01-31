/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: line_mex.h,v 1.5 2011-04-16 12:22:39 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-16 12:22:39 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 ****************************************************************************/

#ifndef LINE_MEX_H
#define LINE_MEX_H

/** 
 * Print out help
 * 
 */
void line_mex_help(void);

/** 
 * Constructor: Manual
 * 
 * @param plhs 
 * @param mx_origin 
 * @param mx_direction 
 * @param mx_dr 
 * @param mx_length 
 * 
 * @return 
 */
bool line_bft_ctor_manual(mxArray*& plhs,
                                 const mxArray* mx_origin,
                                 const mxArray* mx_direction,
                                 const mxArray* mx_dr,
                                 const mxArray* mx_length);

/** 
 * Destructor
 * 
 * @param plhs 
 * @param mx_handle 
 * 
 * @return 
 */
bool line_bft_dtor(mxArray *plhs[], const mxArray* mx_handle);

/** 
 * Get position, make line_bft_get
 * 
 * @param plhs 
 * @param mx_handle 
 * 
 * @return 
 */
bool line_bft_get(mxArray *plhs[], const mxArray* mx_handle,
                         const size_t type, const mxArray* mx_data = NULL);

/** 
 * Function selector
 * 
 * @param nlhs 
 * @param plhs 
 * @param nrhs 
 * @param prhs 
 * 
 * @return success
 */
extern bool line_mex(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]);

#endif

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
