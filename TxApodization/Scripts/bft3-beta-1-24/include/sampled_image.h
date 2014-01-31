#ifndef SAMPLED_IMAGE_H
#define SAMPLED_IMAGE_H

/******************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Image Class                                           *
 *                                                                            *
 * $Id: sampled_image.h,v 1.11 2011-04-28 17:13:22 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-04-28 17:13:22 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * $Log: sampled_image.h,v $
 * Revision 1.11  2011-04-28 17:13:22  jmh
 * Hacks for IEEE
 *
 * Revision 1.9  2010/11/25 20:14:44  jmh
 * Support multiple angles for sampled image
 *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Use common header
 ****************************************************************************/

#ifdef _WIN32
# include <malloc.h>
#else
# include "mm_malloc.h"
#endif

#ifdef __GNUC__
# ifndef HAVE_THREAD
#  define HAVE_THREAD
# endif
# ifndef HAVE_PTHREAD_H
#  define HAVE_PTHREAD_H
# endif
#endif

#include "apodization.h"
#include "common.h"

// Forward declarations
class Aperture;

#include "mex_thread.h"

// Forward declarations
template <class Domain, class Range> class SampleInterpolate;
class Apodization;
class Line;

typedef struct sampled_image {
  size_t nz;
  size_t nx;
  size_t ny;
  float_type dz;
  float_type dx;
  float_type dy;
  float_type offset_z;
  float_type offset_x;
  float_type offset_y;
} sampled_image_t;

class SampledImage {
 public:
  enum Type {
    nearestNeighbour     = 0x00,
    linear               = 0x01,
    cubic                = 0x02,
    spline               = 0x03,
#ifndef _NOFIR
    fir                  = 0x04,
    n_types              = 0x05
#else
    n_types              = 0x04
#endif
  };
  SampledImage();

  SampledImage(const Aperture& Th_t, const Aperture& Th_r,
	       const Apodization& Ah_t, const Apodization& Ah_r,
	       sampled_image_t im);
  bool Beamform(float_type* image, const float_type* rf_data,
		const float_type* delay, const size_t n_rf_samples,
		const uint32_t* trans_no,
		size_t n_angles, const float_type* angles);

#if defined(HAVE_PTHREAD_H)
  void* beamform_threaded(void *ptarg);
#else
  unsigned __stdcall beamform_threaded(void *ptarg);
#endif

  static size_t nthreads;
  sampled_image_t im_geom;
  Aperture m_xmt_aperture;
  Aperture m_rcv_aperture;
  Apodization m_xmt_apodization;
  Apodization m_rcv_apodization;
  Type interp_type;
};

#endif
