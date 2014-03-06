/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Apodization Class                                     *
 *                                                                            *
 * $Id: apodization.h,v 1.59 2011-07-19 20:17:24 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-19 20:17:24 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Introduce common allocator
 *  - Double check by storing reference to all objects
 *  -
 ****************************************************************************/

#ifndef APODIZATION_H
#define APODIZATION_H

#include "common.h"
#include "bfmath.h"

#include "mex.h"

// Forward declarations
class Aperture;

class Apodization {

 public:

  typedef enum Window {
    HAMMING   = 0x00,
    HANN      = 0x01,
    BLACKMAN  = 0x02,
    BARTLETT  = 0x03,
    RECTWIN   = 0x04,
    TUKEY     = 0x05,
    GAUSSIAN  = 0x06, 
    N_WTYPES  = 0x07 
  } Window_t;

  /** 
   * Constructor
   * 
   */
  Apodization();

  Apodization(const Aperture &aperture, const float_type* ref,
              const float_type* distances, const size_t n_distances,
              const float_type* values, const size_t n_elements);
  /** 
   * Copy-constructor
   * 
   * @param other 
   */
  Apodization(const Apodization &other);

  /** 
   * Assignment
   * 
   * @param other 
   * 
   * @return 
   */
  Apodization& operator=(const Apodization &other);

  /** 
   * Destructor
   * 
   */
  ~Apodization();

  /** 
   * 
   * 
   * @param aperture 
   * 
   * @return 
   */
  bool setAperture(const Aperture* aperture);

  /// Representation class
  class DataRep {
  public:

    static size_t next_ID;

    /** 
     * Constructor
     * 
     */
    DataRep();
    
    /** 
     * Constructor
     * 
     * @param ref 
     * @param distances 
     * @param n_distances 
     * @param values 
     * @param n_elements 
     */
    DataRep(const float_type* ref, const float_type* distances,
            const size_t n_distances, const float_type* values,
            const size_t n_elements);

    /** 
     * Destructor
     * 
     */
    ~DataRep();

    /// Aperture
    Aperture m_aperture;
    
    /// Number of distances for parametric apodization
    size_t m_ndistances;

    /// Number of elements
    size_t m_nelements;

    /// Number of active elements
    size_t m_nactive_elements;

    float_type* m_ref;
    float_type* m_euler;
    float_type* m_distances;
    float_type* m_values;
    float_type m_window_param;

    bool m_dynamic;
    bool m_manual;
    bool m_fixed;
    Apodization::Window m_window;

    float_type m_f;
    
    size_t m_id;

    /// Reference counter
    size_t m_cnt;

    ALIGN16_BEGIN float_type m_point_to_aporef[4] ALIGN16_END;
    float_type m_dyn_size;
    ALIGN16_BEGIN float_type m_apo_dir[4] ALIGN16_END;
  };

  DataRep* data;

 public:
  
  /** 
   * Clone
   * 
   * 
   * @return 
   */
  Apodization Clone() const;

  // TODO: Fix this
  Apodization* Clone2() const;

  /** 
   * Dynamic apodization enabled
   * 
   * 
   * @return 
   */
  inline bool dynamic() const {
    return data->m_dynamic;
  }

  /** 
   * Manual apodization enabled
   * 
   * 
   * @return 
   */
  inline bool manual() const {
    return data->m_manual;
  }

  /** 
   * Fixed apodization enabled
   * 
   * 
   * @return 
   */
  inline bool fixed() const {
    return data->m_fixed;
  }

  /** 
   * Apodization reference
   * 
   * @param dim 
   * 
   * @return 
   */
  inline float_type ref(size_t dim) const {
    return data->m_ref[dim];
  }

  /** 
   * operator()
   * 
   * 
   * @return 
   */
  inline float_type operator()(const size_t element,
                               const size_t distance=0) const {
    return data->m_values[element + data->m_nelements*distance];
  }

  /** 
   * Access distances
   * 
   * @param dist 
   * 
   * @return 
   */
  inline float_type distance(size_t dist) const {
    return data->m_distances[dist];
  }

  /** 
   * Set reference position
   * 
   * @param pixel 
   */
  inline void setPosition(float_type* pixel) {
    // TODO: Test including this in dynApo
    data->m_point_to_aporef[0] = pixel[0]-data->m_ref[0];
    data->m_point_to_aporef[1] = pixel[1]-data->m_ref[1];
    data->m_point_to_aporef[2] = pixel[2]-data->m_ref[2];
  }

  /** 
   * Set reference position
   * 
   * @param pixel 
   */
  inline void setRefPosition(float_type* pixel) {
    // TODO: Test including this in dynApo
    data->m_ref[0] = pixel[0];
    data->m_ref[0] = pixel[0];
    data->m_ref[0] = pixel[0];
  }

  /** 
   * Set F number
   * 
   */
  inline void setFnumber() {
    // TODO: Test including this in dynApo
    data->m_dyn_size =
      sqrt(SQUARE(data->m_point_to_aporef[0]) +
	   SQUARE(data->m_point_to_aporef[1]) +
	   SQUARE(data->m_point_to_aporef[2])) / data->m_f;


#ifdef _WIN32
    if (data->m_dyn_size < (data->m_aperture.data->m_pitch/float_type(2.0)))
      data->m_dyn_size = data->m_aperture.data->m_pitch/2.0;
#else
    data->m_dyn_size = std::max(data->m_dyn_size,data->m_aperture.data->m_pitch/float_type(2.0)); 
#endif

}

  /** 
   * fixApo
   * 
   * @param index 
   * @param pixel 
   * 
   * @return 
   */
  INLINE_BEGIN float_type fixApo(const size_t index,
                                 const float_type* pixel,
				 const float_type ap_size) INLINE_END 
  {
    data->m_dyn_size = ap_size;
    return winApo(index,pixel);
  }

  /** 
   * dynApo
   * 
   * @param index 
   * @param pixel 
   * 
   * @return 
   */
  INLINE_BEGIN float_type dynApo(size_t index,
                                 float_type* pixel) __attribute__((pure)) {
    return winApo(index,pixel);
  }

  /** 
   * winApo
   * 
   * @param index 
   * @param pixel 
   * 
   * @return 
   */
  INLINE_BEGIN float_type winApo(const size_t index,
                                 const float_type* pixel) __attribute__((pure)) {

    float_type val = float_type(0.0);
    
    ALIGN16_BEGIN float_type xyz_r[4] ALIGN16_END;
    ALIGN16_BEGIN float_type x2[4] ALIGN16_END;
    
    // Aperture element
    xyz_r[0] = data->m_aperture(0,index);
    xyz_r[1] = data->m_aperture(1,index);
    xyz_r[2] = data->m_aperture(2,index);

    /* Dynamic focus - orthogonal distance from apodization line to
       xmt or rcv element  */

    /* Use point on line tilded according to orientation */
    x2[0] = pixel[0] + data->m_apo_dir[0];
    x2[1] = pixel[1] + data->m_apo_dir[1];
    x2[2] = pixel[2] + data->m_apo_dir[2];

    //mexPrintf("x2[2]: %f\n",x2[2]);
    float_type dist_elm_to_line_pixel_aporef;

    /* TODO: For sampled images, use below or remember to update m_apo_dir */

    /* Distance to line passing through apodization ref and focus point */
    /*
    dist_elm_to_line_pixel_aporef =		\
       dist_point_to_line(xyz_r, data->m_ref, pixel);
    */
    
    /*
    mexPrintf("pixel: %f %f %f\n",pixel[0],pixel[1],pixel[2]);
    mexPrintf("m_ref_x: %f %f %f\n",data->m_ref[0],data->m_ref[1],data->m_ref[2]);
    mexPrintf("x2: %f %f %f\n",x2[0],x2[1],x2[2]);
    */

    /* Distance to line passing through focus point and tilded according to
       orientation */
    
    
    dist_elm_to_line_pixel_aporef =	\
    dist_point_to_line(xyz_r, x2, pixel);
    

    // We assume symmetric apodization
    float_type index_norm = float_type(0.0);

    // Find a way to avoid FLT_EPSILON for single-element apertures
    
#ifndef M_SINGLE_PRECISION
    index_norm = (dist_elm_to_line_pixel_aporef*float_type(2.0) / std::max(data->m_dyn_size,DBL_EPSILON) );
#else
    index_norm = (dist_elm_to_line_pixel_aporef*float_type(2.0) / std::max(data->m_dyn_size,FLT_EPSILON) );
#endif

    switch (data->m_window) {
    case HAMMING:
      if (index_norm <= float_type(1.0)) {
        val = float_type(0.54) + float_type(0.46)*cos(M_PI*index_norm);
      }
      break;
    case HANN:
      if (index_norm < float_type(1.0)) {
        val = float_type(0.5)*(float_type(1.0)+cos(M_PI*index_norm));
      }
      break;
    case BLACKMAN:
      if (index_norm < float_type(1.0))
        val = float_type(0.42) + float_type(0.5)*cos(M_PI*index_norm) + float_type(0.08)*cos(float_type(2.0)*M_PI*index_norm);
      break;
    case BARTLETT:
      if (index_norm < float_type(1.0))
        val = (float_type(1.0) - index_norm);
      break;
    case RECTWIN:
      val = index_norm < float_type(1.0) ? float_type(1.0) : float_type(0.0);
      break;
    case TUKEY:
      if (index_norm < float_type(1.0)) {
	val = float_type(0.5)*(float_type(1.0)+cos(M_PI*(fabs(index_norm)-data->m_window_param)/(1-data->m_window_param)));
	val = fabs(index_norm) < data->m_window_param ? float_type(1.0) : val;
      }
      break;
    case GAUSSIAN:
      if (index_norm < float_type(1.0))
	val = exp(-float_type(0.5)*(SQUARE(index_norm)*SQUARE(data->m_window_param)));
      break;
    default:
      break;
    }
    return val;
  }
  /*  bool getFocusDelays(float_type* focus_delays); */

 private:
  /** 
   * clone
   * 
   */
  void clone();
};

#endif
