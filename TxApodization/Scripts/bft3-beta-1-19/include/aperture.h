/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Aperture Class                                        *
 *                                                                            *
 * $Id: aperture.h,v 1.45 2011-08-04 18:23:38 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-08-04 18:23:38 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * Revision 1.9  2009/11/04 09:45:00  jmh
 * Corrupted double linked list
 *
 *****************************************************************************/

/****************************************************************************
 * TODO: - Add Clone method for Deep copy (multi-threading)
 *       - Introduce common allocator
 *       - Reference all objects
 *       - Declare fs and c extern if globals
 ****************************************************************************/

#ifndef APERTURE_H
#define APERTURE_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "common.h"

#ifdef _WIN32
# include <malloc.h>
#else
# include "mm_malloc.h"
#endif

class Aperture {
public:

  enum Type {
    custom               = 0x00,
    linear_array         = 0x01,
    convex_array         = 0x02,
    two_dim_array        = 0x03,
    n_types              = 0x04
  };

//@{
  /** Static variables */
  static float_type _fs;
  static float_type _f0;
  static float_type _c;
//@}

  /** 
   * Constructor
   * 
   */
  Aperture();

  /** 
   * Contructor: Manual
   * 
   * @param pos_vec 
   * @param length 
   */
  Aperture(const float_type* pos_vec, size_t length);
  
  /** 
   * Constructor: Field II apertures
   * 
   * @param n_elements 
   * @param pitch 
   */
  Aperture(size_t n_elements, float_type pitch);

  /** 
   * Construct: Field II 2D apertures
   * 
   * @param n_xelements 
   * @param n_yelements 
   * @param xpitch 
   * @param ypitch 
   */
  Aperture(size_t n_xelements, size_t n_yelements,
	   float_type xpitch, float_type ypitch);

  /** 
   * Constructor: Field II convex apertures
   * 
   * @param n_elements 
   * @param width 
   * @param kerf 
   * @param radius 
   */
  Aperture(size_t n_elements, float_type width, float_type kerf,
	   float_type radius);

  /** 
   * Copy-constructor
   * 
   * @param other 
   */
  Aperture(const Aperture &other);

  /** 
   * Assignment
   * 
   * @param other 
   * 
   * @return 
   */
  Aperture& operator=(const Aperture &other);

  /** 
   * Destructor:
   * 
   */
  ~Aperture();

  /** 
   * 
   * 
   * @param pos 
   * @param npos 
   * 
   * @return 
   */
  bool setPositions(const float_type* pos, const size_t npos);

  /** 
   * 
   * 
   * @param focus 
   * 
   * @return 
   */
  bool setFocus(const float_type* focus);

  /** 
   * 
   * 
   * @param delays 
   * 
   * @return 
   */
  bool setDelays(const float_type* delays);

  /** 
   * 
   * 
   * @param center_focus 
   * @param n_emissions 
   * 
   * @return 
   */
  bool setCenterFocus(const float_type* center_focus,
		      size_t n_emissions);


  bool setOrientation(const float_type* orientation,
		      size_t n_emissions);
  /** 
   * 
   * 
   * @param focus_delays 
   * 
   * @return 
   */
  bool getFocusDelays(float_type* focus_delays);

  inline bool ppwave() const {
    return data->m_ppwave;
  }

  Aperture* Clone() const;

  // Matrix class
  template <class T> class PMatrix {
  public:

    T** m_data;
 
    PMatrix(size_t nx, size_t ny) { // ny = 3, nx = np
      
      m_data = (T**) _mm_malloc(4*((ny+3)/4)*sizeof(T*),16);
      m_data[0] = (T*) _mm_malloc(4*((nx*ny+3)/4)*sizeof(T),16);
      
      for (size_t i=1 ; i<ny ; i++) {
	/// Row pointers
        m_data[i] = m_data[0] + nx * i;
      }
    }
    ~PMatrix() {
      if (m_data) {
	if (m_data[0]) {
	  _mm_free(m_data[0]);
	}
	_mm_free(m_data);
	m_data = NULL;
      }
    }
  };

  // Representation class
  class DataRep {
  public:

    static size_t next_ID;

    /** 
     * 
     * 
     */
    DataRep();

    /** 
     * 
     * 
     * @param pos_vec 
     * @param length 
     * @param Type 
     */
    DataRep(const float_type* pos_vec, size_t length, size_t n_emissions, Aperture::Type = Aperture::custom);

    /** 
     * 
     * 
     */
    ~DataRep();
 
    /// Number of positions
    size_t m_npos;

    /// Number of emissions
    size_t m_nemissions;

    /// Actual data
    PMatrix<float_type>* m_pos;

    /// Center focus
    float_type* m_center_focus;

    /// Orientation
    float_type* m_euler;

    /// Virtual source or focus
    float_type* m_focus;
    
    /// Fixed delays added to time-of-flight calculation
    float_type* m_delays;

    /// References to static variables
    float_type* m_fs;
    float_type* m_f0;
    float_type* m_c;

    /// Unique identifier
    size_t m_id;

    Aperture::Type m_type;

    bool m_ppwave;

    /// Reference counter
    size_t m_cnt;

    /// Pitch
    float_type m_pitch;

    /// Radius used for convex arrays
    float_type m_radius;
  };

  DataRep* data;

  /** 
   * 
   * 
   * 
   * @return 
   */
  inline size_t n_elements() const {
    return data->m_npos;
  }

  /** 
   * 
   * 
   * 
   * @return 
   */
  inline float_type operator()(const size_t &dim, const size_t &index) const {
    return data->m_pos->m_data[dim][index];
  }

private:
  /** 
   * 
   * 
   */
  void clone();
};

// Template instantiation
template class Aperture::PMatrix<float_type>;

#endif
