/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Gateway routines                                      *
 *                                                                            *
 * $Id: aperture_mex.h,v 1.8 2011-04-06 09:47:54 jmh Exp $  *
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-06 09:47:54 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * $Log: aperture_mex.h,v $
 * Revision 1.8  2011-04-06 09:47:54  jmh
 * *** empty log message ***
 *
 * Revision 1.7  2011-04-05 20:10:27  jmh
 * *** empty log message ***
 *
 * Revision 1.6  2011/03/26 14:38:00  jmh
 * *** empty log message ***
 *
 * Revision 1.5  2011/03/13 14:44:03  jmh
 * *** empty log message ***
 *
 * Revision 1.4  2010/07/17 18:06:23  jmh
 * *** empty log message ***
 *
 * Revision 1.3  2010/07/17 15:00:03  jmh
 * line apodization property
 *
 * Revision 1.2  2010/05/05 09:58:30  jmh
 * *** empty log message ***
 *
 * Revision 1.1  2010/05/04 11:28:15  jmh
 * *** empty log message ***
 *
 *****************************************************************************/

/****************************************************************************
 * TODO: - Move modifier to Aperture Class
 *       - Include headers
 ****************************************************************************/

#ifndef APERTURE_MEX_H
#define APERTURE_MEX_H

#define Need_aperture_mex_gateway 1

/** 
 * Print out help
 * 
 */
extern void aperture_mex_help(void);

/** 
 * Constructor: Manual
 * 
 * @param plhs Matlab mex-handle
 * @param mx_pos_vec position vector
 * 
 * @return 
 */
extern bool aperture_bft_ctor_manual(mxArray*& plhs,
				     const mxArray* mx_pos_vec);

/** 
 * Constructor: Field II apertures
 * 
 * @param plhs 
 * @param mx_type Type, e.g. "linear_array" 
 * @param mx_no_elements 
 * @param mx_width 
 * @param mx_kerf 
 * @param mx_fs
 * @param mx_f0
 * 
 * @return 
 */
extern bool aperture_bft_ctor_field(mxArray*& plhs,
				    const mxArray* mx_type,
				    const mxArray* mx_no_elements,
				    const mxArray* mx_pitch,
				    const mxArray* mx_f0);

/** 
 * Destructor
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle 
 * 
 * @return 
 */
extern bool aperture_bft_dtor(mxArray *plhs[],
			      const mxArray* mx_handle);

/** 
 * Get properties of Aperture
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle 
 * @param type Property index
 * 
 * @return 
 */
extern bool aperture_bft_get(mxArray *plhs[], const mxArray* mx_handle, const size_t type);

/** 
 * Set Aperture property
 * 
 * @param plhs 
 * @param mx_handle Matlab mex-handle 
 * @param mx_data 
 * @param type Property index
 * 
 * @return 
 */
extern bool aperture_bft_set(mxArray *plhs[], const mxArray* mx_handle,
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
extern bool aperture_mex(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]);

#endif
