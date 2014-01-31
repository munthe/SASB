/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: apodization_mex.h,v 1.5 2011-04-16 12:22:39 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-16 12:22:39 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *
 * $Log: apodization_mex.h,v $
 * Revision 1.5  2011-04-16 12:22:39  jmh
 * *** empty log message ***
 *
 * Revision 1.4  2011/04/05 20:10:27  jmh
 * *** empty log message ***
 *
 * Revision 1.3  2011/03/26 14:38:00  jmh
 * *** empty log message ***
 *
 * Revision 1.2  2010/07/17 15:00:03  jmh
 * line apodization property
 *
 * Revision 1.1  2010/07/17 01:59:08  jmh
 * *** empty log message ***
 *
 *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 * 
 ****************************************************************************/

#ifndef APODIZATION_MEX_H
#define APODIZATION_MEX_H

//#define Need_apodization_mex_gateway 1

/** 
 * Print out help
 * 
 */
void apodization_mex_help(void);

/** 
 * Constructor: Fixed
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle
 * @param mx_pos_ref Apodization reference position
 * @param mx_distances Distances
 * @param mx_values Windows, one for each distance
 * 
 * @return 
 */
bool apodization_bft_ctor_fixed(mxArray*& plhs,
                                        const mxArray* mx_handle,
                                        const mxArray* mx_pos_ref,
                                        const mxArray* mx_distances,
				        const mxArray* mx_values);

/** 
 * Destructor
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle
 * 
 * @return 
 */
bool apodization_bft_dtor(mxArray *plhs[], const mxArray* mx_handle);

/** 
 * Get properties of Apodization
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle
 * @param type Property index
 * 
 * @return 
 */
bool apodization_bft_get(mxArray *plhs[],
                                const mxArray* mx_handle, const size_t type);

/** 
 * Set Apodization property
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle
 * @param mx_data 
 * @param type Property selector
 * 
 * @return 
 */
bool apodization_bft_set(mxArray *plhs[], const mxArray* mx_handle,
                                const mxArray* mx_data, const size_t type);

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
bool apodization_mex(const int nlhs, mxArray *plhs[],
			    const int nrhs, const mxArray *prhs[]);

#endif
